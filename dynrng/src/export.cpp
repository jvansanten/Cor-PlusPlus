/* export.cpp
* this file is part of Cor-PlusPlus for CORSIKA
*
* Copyright (C) <2017> <Dominik Baack>
*		All rights reserved.
*
* 	This software may be modified and distributed under the terms
* 	of the LGPL license. See the LICENSE file for details.
*/
#include "dynrng/export.h"

#include <array>
#include <random>
#include <sstream>

#ifdef __RNG_RANDOM123__
#include "Random123/threefry.h"
#include "Random123/MicroURNG.hpp"
typedef r123::MicroURNG<r123::Threefry4x64> RNGType;
std::array<RNGType, 6> generators = {
	RNGType({{}},{{}}), // gruuuh dirty
	RNGType({{}},{{}}),
	RNGType({{}},{{}}),
	RNGType({{}},{{}}),
	RNGType({{}},{{}}),
	RNGType({{}},{{}})
};
#else
std::array<std::mt19937_64, 6> generators;
#endif

extern "C" void dynrng_seed_(unsigned int *gen, unsigned int* seeds)
{
	if (*gen >= generators.size()) {
		std::ostringstream oss;
		oss << "RNG stream " << *gen << " out of range";
		throw std::out_of_range(oss.str().c_str());
	}
#ifdef __RNG_RANDOM123__
	// in default mode, only use the lower 32 bits for key and lower 64 for counter
	RNGType::key_type key = {{0,0,0,seeds[0]}};
	RNGType::ctr_type counter = {{0,0,0,seeds[1]+1000000000*size_t(seeds[2])}};
	generators[*gen-1].reset(counter, key);
#else
	generators[*gen-1].seed(seeds[0]);
	generators[*gen-1].discard(seeds[1]+1000000000*size_t(seeds[2]));
#endif
}

extern "C" void dynrng_getrandom_(unsigned int* gen, int* N, double* numbers)
{
	std::uniform_real_distribution<double> distribution(0.0,1.0);

	if(*gen > 5)	// todo: change to debug mode only
	{
		return;
	}

	for(int i=0; i < *N; i++)
	{
		numbers[i] = distribution(generators[*gen-1]);
	}
}

