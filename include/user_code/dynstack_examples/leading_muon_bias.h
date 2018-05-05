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

#include <list>

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

namespace dynstack {
namespace leading_muon_bias {

void reset();
void close();

void particleIn(const Particle *p);
void particleOut(const Particle *p);
bool isPending(const Particle *p);
bool isDead();

struct EnergySort {
	bool operator()(const Particle &a, const Particle &b);
};

void header(lib::data::EventHeader &block);
void footer(lib::data::EventEnd &block);

void set_bias_factor(float bias_factor);
void set_bias_target(int target);

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

}}

