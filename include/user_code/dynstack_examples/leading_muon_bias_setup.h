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

#include "leading_muon_bias.h"

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

void reset() { dynstack::leading_muon_bias::reset(); }
void close() { dynstack::leading_muon_bias::close(); }

auto dynstack_setup(std::vector<long> sizes, std::vector< std::list<std::string> > arguments )
{
	if (sizes.empty()) {
		std::cerr << "No stack size set" << std::endl;
		exit(-1);
	}
	if (arguments.size() != 1) {
		std::cerr << "You must provide a DYNSTACK_P line in the steering card, e.g." << std::endl;
		std::cerr << "DYNSTACK_P bias_target mu bias_factor 1e-3" << std::endl;
		std::cerr << "to select approximately 0.1% of showers" << std::endl;
		std::cerr << arguments.size() << std::endl;
		exit(-1);
	}
	auto &line = arguments.front();
	for (auto word = line.begin(); word != line.end(); word++) {
		if (*word == "bias_factor") {
			float bias_factor;
			if (!read_argument(word, line.end(), bias_factor, 0.f, 1.f)) {
				exit(-1);
			}
			dynstack::leading_muon_bias::set_bias_factor(bias_factor);
		} else if (*word == "bias_target") {
			if (++word == line.end()) {
				std::cerr << "no argument provided for bias_target (mu, numu, nue)" << std::endl;
				exit(-1);
			}
			if (*word == "mu")
				dynstack::leading_muon_bias::set_bias_target(0);
			else if (*word == "numu")
				dynstack::leading_muon_bias::set_bias_target(1);
			else if (*word == "nue")
				dynstack::leading_muon_bias::set_bias_target(2);
			else {
				std::cerr << "bias_target got unknown argument '"<<*word<<"'" << std::endl;
				exit(-1);
			}
		} else {
			std::cerr << "unknown option '"<<*word<<"'" << std::endl;
			exit(-1);
		}
	}
	SHeaderManager().register_evth_callback(dynstack::leading_muon_bias::header);
	SHeaderManager().register_evte_callback(dynstack::leading_muon_bias::footer);
	
	auto stack = std::make_unique<dynstack::leading_muon_bias::stack_type>(sizes.front());
	return std::move(stack);
}

