/* default.cpp
 * this file is part of RemoteControl for CORSIKA
 *
 * Copyright (C) <2016> <Dominik Baack>
 *		All rights reserved.
 *
 * 	This software may be modified and distributed under the terms
 * 	of the LGPL license. See the LICENSE file for details.
 */


#include <functional>
#include <iostream>
#include <mutex>
#include <cmath>
#include <condition_variable>

#include "user_code/remotecontrol_examples/default.h"
#include "user_code/dynstack_examples/leading_muon_bias.h"
#include "remote_control/export.h"

struct {
	double type, energy, theta, phi;
	int bias_target;
	double muon_bias;
	std::array<double,4> elcut;
	bool valid;
} primary_particle = { 0, 0, 0, 0, 0, 1, {NAN,NAN,NAN,NAN}, false};

struct {
	std::mutex mutex;
	std::condition_variable cv;
} primary_particle_lock;

remote_control::communication::Packet
recv_primary(std::vector<uint8_t> msg)
{
	if (msg.size() < 4*sizeof(double)) {
		throw std::runtime_error("primary message has the wrong size");
	}
	std::unique_lock<std::mutex> lock(primary_particle_lock.mutex);
	while (primary_particle.valid) {
		primary_particle_lock.cv.wait(lock);
	}
	
	double *data = reinterpret_cast<double*>(msg.data());
	primary_particle.type = data[0];
	primary_particle.energy = data[1];
	primary_particle.theta = data[2];
	primary_particle.phi = data[3];
	if (msg.size() > 4*sizeof(double)) {
		primary_particle.bias_target = int(data[4]);
	} else {
		primary_particle.bias_target = 0;
	}
	if (msg.size() > 5*sizeof(double)) {
		primary_particle.muon_bias = data[5];
	} else {
		primary_particle.muon_bias = 1.;
	}
	if (msg.size() >= 10*sizeof(double)) {
		std::copy(&data[6], &data[6]+4, primary_particle.elcut.begin());
	} else {
		std::fill(primary_particle.elcut.begin(), primary_particle.elcut.end(), NAN);
	}
	primary_particle.valid = true;
	
	primary_particle_lock.cv.notify_one();
	
	return remote_control::communication::Packet();
}

void extprm_(double *type, double *energy, double *theta, double *phi)
{
	std::unique_lock<std::mutex> lock(primary_particle_lock.mutex);
	while (!primary_particle.valid) {
		primary_particle_lock.cv.wait(lock);
	}
	
	*type = primary_particle.type;
	*energy = primary_particle.energy;
	*theta = primary_particle.theta;
	*phi = primary_particle.phi;
	dynstack::leading_muon_bias::set_bias_target(primary_particle.bias_target);
	dynstack::leading_muon_bias::set_bias_factor(primary_particle.muon_bias);
	for (unsigned int i=0; i < 4; i++) {
		if (std::isfinite(primary_particle.elcut[i])) {
			SBasic().setELCUT(Basic::energy_cut(i), primary_particle.elcut[i]);
		}
	}
	primary_particle.valid = false;
	
	primary_particle_lock.cv.notify_one();
}
