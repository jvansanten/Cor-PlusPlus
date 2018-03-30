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

#include "user_code/remotecontrol_examples/default.h"
#include "user_code/dynstack_examples/leading_muon_bias.h"
#include "remote_control/export.h"

struct {
	double type, energy, theta, phi;
	double muon_bias;
	bool valid;
} primary_particle = { 0, 0, 0, 0, 1, false};

struct {
	std::mutex mutex;
	std::condition_variable cv;
} primary_particle_lock;

remote_control::communication::Packet
recv_primary(std::vector<uint8_t> msg)
{
	if (msg.size() != 5*sizeof(double)) {
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
	primary_particle.muon_bias = data[4];
	primary_particle.valid = true;
	
	primary_particle_lock.cv.notify_one();
	
	return remote_control::communication::Packet();
}

void extprm_(double *type, double *energy, double *theta, double *phi)
{
	// std::cerr << " --- getting a primary --- " << primary_particle.valid << std::endl;
	
	std::unique_lock<std::mutex> lock(primary_particle_lock.mutex);
	while (!primary_particle.valid) {
		// std::cerr << " --- waiting in extprm --- " << primary_particle.valid << std::endl;
		primary_particle_lock.cv.wait(lock);
	}
	// std::cerr << " --- done waiting in extprm --- " << primary_particle.valid << std::endl;
	
	
	*type = primary_particle.type;
	*energy = primary_particle.energy;
	*theta = primary_particle.theta;
	*phi = primary_particle.phi;
	dynstack::leading_muon_bias::set_bias_factor(primary_particle.muon_bias);
	primary_particle.valid = false;
	
	primary_particle_lock.cv.notify_one();
	
	
}