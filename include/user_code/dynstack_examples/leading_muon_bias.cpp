/* leading_muon_bias.cpp
* this file is part of Dynstack/RemoteControl for CORSIKA
*
* Copyright (C) <2018> <Jakob van Santen>
*		All rights reserved.
*
* 	This software may be modified and distributed under the terms
* 	of the LGPL license. See the LICENSE file for details.
*/

#include <cmath>
#include <limits>

#include "particle_deduction.h"
#include "basic/basic.h"
#include "basic/header_manager.h"

#include "leading_muon_bias.h"

typedef DeductedParticleType Particle;

// Use CORSIKA-internal RNG stream for reproducibility
extern "C" {
void rmmard_(double *rvec, int *num, int *seq);
}

namespace dynstack {
namespace leading_muon_bias {

/// The state of the shower-kill decision
enum {
	RESET,     // Start of shower, waiting for first primary
	PENDING,   // Kill decision pending; keep muon ancestors ordered by energy
	COMMITTED, // Shower kept; all particles to secondary stack
	KILLED     // Shower terminated; clear stack and continue 
} state = RESET;

/// Energy/nucleon of the primary nucleus
double primaryEnergy = 0;
double max_x = 0;

/// Number of similar showers this event represents
double weight = 1;

/// Stack statistics
size_t n_pops = 0;
size_t n_showers = 0;
size_t n_killed = 0;
size_t n_pops_before_kill = 0;
size_t n_pops_total = 0;
// We store this in the event header, so keep it always in single precision
float bias_factor = 1e-3;

void set_bias_factor(float bias)
{
    bias_factor = bias;
}

class ElbertYield {
public:
	ElbertYield() : a(14.5), p1(0.757+1), p2(5.25), logN0(0)
	{}
	
	ElbertYield(double primaryEnergy, unsigned A, double cos_theta) :
	    a(14.5), p1(0.757+1), p2(5.25),
	    logN0(std::log(a*A*A/primaryEnergy/effective_costheta(cos_theta)))
	{}
	
	double operator()(double x) const
	{
		return std::exp(logYield(std::log(x)));
	}
	
	// Invert the Elbert formula to find the energy (in units of the primary
	// energy/nucleon) above which N muons are expected per shower.
	double invert(double N) const
	{
		const double xtol(1e-5);
		const double y = std::log(N);
		double x0 = std::min((logN0 - y)/p1, -xtol);
		for (int i=0; i < 20; i++) {
			const double xp = std::min(x0 - (logYield(x0)-y)/dlogYield(x0), -xtol);
			if (std::abs(xp-x0) < xtol)
				return std::exp(xp);
			x0 = xp;
		}
		
		throw std::runtime_error("Root solver did not converge");
		return std::exp(x0);
	}

private:
	double a, p1, p2;
	double logN0;
	
	// Effective local atmospheric density correction from [Chirkin]_.
	// 
	// .. [Chirkin] D. Chirkin. Fluxes of atmospheric leptons at 600-GeV - 60-TeV. 2004. http://arxiv.org/abs/hep-ph/0407078
	static double effective_costheta(double x)
	{
		std::array<double,5> p = {0.102573, -0.068287, 0.958633, 0.0407253, 0.817285};
		return std::sqrt((std::pow(x,2) + std::pow(p[0],2) + p[1]*std::pow(x,p[2]) + p[3]*std::pow(x,p[4]))/(1 + std::pow(p[0],2) + p[1] + p[3]));
	}
	
	double logYield(double logx) const
	{
		return logN0 - p1*logx + p2*std::log(1-std::exp(logx));
	}
	
	double dlogYield(double logx) const
	{
		return -p1 - p2*std::exp(logx)/(1-std::exp(logx));
	}
};

class ElbertBias : private ElbertYield {
public:
	ElbertBias() {}
	// Find x_threshold such that the probability of a shower producing
	// at least 1 muon above the threshold is bias_factor
	ElbertBias(double primaryEnergy, unsigned A, double cos_theta, double bias_factor)
		: ElbertYield(primaryEnergy,A,cos_theta),
		  x_threshold(bias_factor < 1 ? invert(-std::log(1-bias_factor)) : 0)
	{}
	
	// NB: the x of the lead muon will be stored in single precision in the
	// event trailer, so we truncate it here to ensure that we can reproduce
	// the result later.
	double operator()(float x) const
	{
		if (x >= x_threshold) {
			// above threshold, accept all showers
			return 1;
		} else {
			// Below threshold, accept according the the probability for the
			// shower to produce at least one muon above threshold, relative to
			// the probability of producing at least 1 muon above the current
			// energy. At small x this converges to the bias factor. This 
			// happens because nearly all showers produce low energy muons; 
			// once we reach this low x regime, we just have to pick 1 out of 
			// 1/bias_factor showers at random.
			return (1-std::exp(-ElbertYield::operator()(x_threshold)))/(1-std::exp(-ElbertYield::operator()(x)));
		}
	}
	
	double threshold() const { return x_threshold; }
private:
	double x_threshold;
};

ElbertBias kill_bias;

/// Probability of accepting a shower if its highest-energy muon carries a
/// fraction `x` of the primary energy/nucleon.
double bias(double x)
{
	return kill_bias(x);
}

void reset() {
	if (state == COMMITTED || state == PENDING) {
		n_pops_total += n_pops;
	} else if (state == KILLED) {
		n_pops_before_kill += n_pops;
	}
	state = RESET;
	primaryEnergy = 0;
	max_x = 0;
	weight = 1;
	n_pops = 0;
}

void close() {
	std::cerr << "(muon-bias) killed " << n_killed << " of " << n_showers << " showers " << std::endl;
	std::cerr << "(muon-bias) killed showers: propagated " << double(n_pops_before_kill)/n_killed << " particles " << std::endl;
	std::cerr << "(muon-bias) complete showers: propagated " << double(n_pops_total)/(n_showers-n_killed) << " particles " << std::endl;
}

int getType(const Particle &p)
{
	return (int)p[Particle::PARTICLE_TYPE];
}

int getNucleonNumber(const Particle &p)
{
	int ptype = getType(p);
	if (ptype == 14) {
		return 1;
	} else if (ptype > 200) {
		return ptype / 100;
	} else {
		return 0;
	}
}

double getEnergy(const DeductedParticleType &p)
{
	double pama = SBasic().particleRestMass(getType(p));
	return p[Particle::GAMMA]*pama;
}

void particleIn(const Particle *p)
{
	if (state == RESET) {
		// First push after reset is the shower primary
		primaryEnergy = getEnergy(*p)/getNucleonNumber(*p);
		if (bias_factor < 1) {
			kill_bias = ElbertBias(getEnergy(*p), getNucleonNumber(*p), (*p)[Particle::COSTHE], bias_factor);
			if (n_showers == 1)
				std::cerr << "(muon-bias) Ep="<<getEnergy(*p)<<" A="<<getNucleonNumber(*p)<<" P(N_mu(x>"<<kill_bias.threshold()<<")>0)="<<bias_factor<<std::endl;
		}
		// Go straight to secondary stack if no threshold set
		state = (bias_factor < 1) ? PENDING : COMMITTED;
		n_showers++;
	}
}

double uniform()
{
	double v;
	int len(1), seq(2);
	rmmard_(&v, &len, &seq);
	return v;
}

void particleOut(const Particle *p)
{
	n_pops++;
	// if in unbiased LIFO mode, need to search for the highest-energy muon
	// explicitly
	if (bias_factor == 1 && (getType(*p) == 5 || getType(*p) == 6)) {
		double x = getEnergy(*p)/primaryEnergy;
		if (x > max_x)
			max_x = x;
	} else if (state == PENDING) {
		// Fraction of primary energy in this particle
		max_x = getEnergy(*p)/primaryEnergy;
		// Probability of acceptance
		double prob = bias(max_x)*weight;
		// Is this the leading muon?
		bool isMuon = (getType(*p) == 5 || getType(*p) == 6);
		
		/// Since the acceptance probability decreases monotonically with 
		/// decreasing `x`, the probability at current `x` is an upper limit on
		/// the acceptance probability for the highest-energy muon in the
		/// shower. Here we make the final decision if we have the highest-
		/// energy muon. Otherwise, we try kill the shower early and record
		/// the probability if the shower survives.
		if ((isMuon || prob < 0.9)) {
			if (uniform() >= prob) {
				// Stop the shower immediately
				weight = 0;
				state = KILLED;
				n_killed++;
			} else {
				// Account for acceptance probability
				weight /= prob;
				if (isMuon) {
					state = COMMITTED;
				}
			}
		}
	}
}

bool isPending(const Particle *p)
{
	switch (getType(*p)) {
		// Defer EM shower and neutrinos to secondary stack
		case   1: // gamma
		case   2: // eplus
		case   3: // eminus
		case  66: // nu_e
		case  67: // nu_ebar
		case  68: // nu_mu
		case  69: // num_mubar
		case 133: // nu_tau
		case 134: // nu_taubar
			return false;
			break;
		// Keep possible muon ancestors sorted by energy until the leading muon
		// is encountered
		default:
			return (state == PENDING);
	}
}

bool isDead()
{
	return (state == KILLED);
}

bool EnergySort::operator()(const Particle &a, const Particle &b)
{
	return getEnergy(a) > getEnergy(b);
}

void header(lib::data::EventHeader &block)
{
	// Write bias factor to field normally used for energy_interesting
	// in ICECUBE1 (original neutrino kill threshold) mode
	block.write(lib::data::EventHeader::index(219), float(bias_factor));
}

void footer(lib::data::EventEnd &block)
{
	// Write prescale bias to an unused field in event trailer
	// FIXME: is there a better place to put these?
	block.write(lib::data::EventEnd::index(266), float(weight));
	block.write(lib::data::EventEnd::index(267), float(max_x));
}

}}


