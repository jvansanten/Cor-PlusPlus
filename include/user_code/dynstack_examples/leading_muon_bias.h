/* leading_muon_bias.h
* this file is part of Dynstack/RemoteControl for CORSIKA
*
* Copyright (C) <2018> <Jakob van Santen>
*		All rights reserved.
*
* 	This software may be modified and distributed under the terms
* 	of the LGPL license. See the LICENSE file for details.
*/
#pragma once

#include <cmath>
#include <limits>

#include "particle_deduction.h"
#include "basic/basic.h"
#include "basic/header_manager.h"

#include "dynstack/stack/storage/ordered_stack.h"
#include "dynstack/stack/storage/lifo_stack.h"
#include "dynstack/stack/wrapper/filter_stack.h"
#include "dynstack/stack/wrapper/in_callback_stack.h"
#include "dynstack/stack/wrapper/out_callback_stack.h"
#include "dynstack/stack/advanced/if_stack.h"
#include "dynstack/stack/advanced/kill_stack.h"

typedef DeductedParticleType Particle;

// Use CORSIKA-internal RNG stream for reproducibility
extern "C" {
void rmmard_(double *rvec, int *num, int *seq);
}

namespace {

/// The state of the shower-kill decision
enum {
	RESET,     // Start of shower, waiting for first primary
	PENDING,   // Kill decision pending; keep muon ancestors ordered by energy
	COMMITTED, // Shower kept; all particles to secondary stack
	KILLED     // Shower terminated; clear stack and continue 
} state = RESET;

/// Energy/nucleon of the primary nucleus
double primaryEnergy = 0;
double max_x = 1;

/// Number of similar showers this event represents
double weight = 1;

/// Stack statistics
size_t n_pops = 0;
size_t n_showers = 0;
size_t n_killed = 0;
size_t n_pops_before_kill = 0;
size_t n_pops_total = 0;
double target_num_muons = 1e-3;
double x_threshold = 0.1;
double bias_power = 2.;

/// Probability of accepting a shower if its highest-energy muon carries a
/// fraction `x` of the primary energy/nucleon.
double bias(double x)
{
	return std::min(std::pow(x/x_threshold, bias_power), 1.);
}

// Effective local atmospheric density correction from [Chirkin]_.
// 
// .. [Chirkin] D. Chirkin. Fluxes of atmospheric leptons at 600-GeV - 60-TeV. 2004. http://arxiv.org/abs/hep-ph/0407078
double effective_costheta(double x)
{
	std::array<double,5> p = {0.102573, -0.068287, 0.958633, 0.0407253, 0.817285};
	return std::sqrt((std::pow(x,2) + std::pow(p[0],2) + p[1]*std::pow(x,p[2]) + p[3]*std::pow(x,p[4]))/(1 + std::pow(p[0],2) + p[1] + p[3]));
}

// Invert the Elbert formula to find the energy (in units of the primary
// energy/nucleon) above which N muons are expected per shower.
double invert_elbert(double N, double primaryEnergy, unsigned A, double cos_theta)
{
	const double a(14.5), p1(0.757+1), p2(5.25);
	const double logN0 = std::log(a*A*A/primaryEnergy/effective_costheta(cos_theta));
	
	auto logYield = [&](double logx) { return logN0 - p1*logx + p2*std::log(1-std::exp(logx)); };
	auto dlogYield = [&](double logx) { return -p1 - p2*std::exp(logx)/(1-std::exp(logx)); };
	
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

void reset() {
	if (state == COMMITTED || state == PENDING) {
		n_pops_total += n_pops;
	} else if (state == KILLED) {
		n_pops_before_kill += n_pops;
	}
	state = RESET;
	primaryEnergy = 0;
	max_x = 1;
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
		x_threshold = invert_elbert(target_num_muons, getEnergy(*p), getNucleonNumber(*p), (*p)[Particle::COSTHE]);
		if (n_showers == 1)
			std::cerr << "(muon-bias) Ep="<<getEnergy(*p)<<" A="<<getNucleonNumber(*p)<<" N(x>"<<x_threshold<<")="<<target_num_muons<<std::endl;
		// Go straight to secondary stack if no threshold set
		state = (x_threshold > 0) ? PENDING : COMMITTED;
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
	if (state == PENDING) {
		// Fraction of primary energy in this particle
		double x = getEnergy(*p)/primaryEnergy;
		// Probability of acceptance
		double prob = bias(x)*weight;
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

struct EnergySort {
	bool operator()(const Particle &a, const Particle &b)
	{
		return getEnergy(a) > getEnergy(b);
	}
};

void footer(lib::data::EventEnd &eh)
{
	// Write prescale bias to an unused field in event trailer
	eh.write(lib::data::EventEnd::index(266), float(weight));
}

typedef 
	dynstack::advanced::KillStack<
		dynstack::wrapper::OutCallbackStack<
			dynstack::wrapper::InCallbackStack<
				// Propagate particles in descending order of energy, switching
				// to the fast, contiguous stack once the shower is committed
				dynstack::advanced::IfStack<
					dynstack::storage::OrderedStack<Particle, EnergySort>,
					dynstack::storage::LIFO_Stack<Particle>,
					isPending>,
				particleIn>,
			particleOut>,
		isDead> stack_type;

template <typename Iterator, typename Target>
bool read_argument(Iterator &word, Iterator &&end, Target &target, Target lo, Target hi)
{
	auto &arg = *word;
	std::stringstream sstr;
	if (++word == end) {
		std::cerr << arg << " requires an argument" << std::endl;
		return false;
	}
	sstr << *word;
	sstr >> target;
	if (!(target >= lo && target <= hi)) {
		std::cerr << arg << " must be in ["<<lo<<","<<hi<<"]" << std::endl;
		return false;
	}
	return true;
}

auto dynstack_setup(std::vector<long> sizes, std::vector< std::list<std::string> > arguments )
{
	if (sizes.empty()) {
		std::cerr << "No stack size set" << std::endl;
		exit(-1);
	}
	if (arguments.size() != 1) {
		std::cerr << "You must provide a DYNSTACK_P line in the steering card, e.g." << std::endl;
		std::cerr << "DYNSTACK_P muon per_shower 1e-3 bias_power 2" << std::endl;
		std::cerr << "to kill showers with a probability (x/x0)^2, where x0 is " << std::endl;
		std::cerr << "chosen so such that the Elbert formula predicts 1e-3 muons" << std::endl;
		std::cerr << "per shower with x>x0" << std::endl;
		std::cerr << arguments.size() << std::endl;
		exit(-1);
	}
	auto &line = arguments.front();
	auto word = line.begin();
	if (word == line.end() || *word != "muon") {
		std::cerr << "First DYNSTACK argument must be \"muon\" (got \""<<*word<<"\")" << std::endl;
		exit(-1);
	}
	while (++word != line.end()) {
		if (*word == "per_shower") {
			if (!read_argument(word, line.end(), target_num_muons, 0., std::numeric_limits<double>::infinity())) {
				exit(-1);
			}
		} else if (*word == "bias_power") {
			if (!read_argument(word, line.end(), bias_power, 0., 1e2)) {
				exit(-1);
			}
		}
	}
	SHeaderManager().register_evte_callback(footer);
	
	auto stack = std::make_unique<stack_type>(sizes.front());
	return std::move(stack);
}

}




