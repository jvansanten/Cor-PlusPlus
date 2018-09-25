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

std::array<std::mt19937_64, 6> generators;


extern "C" void dynrng_seed_(unsigned int *gen, unsigned int* seeds)
{
	if (*gen >= generators.size()) {
		std::ostringstream oss;
		oss << "RNG stream " << *gen << " out of range";
		throw std::out_of_range(oss.str().c_str());
	}
	generators[*gen-1].seed(seeds[0]);
	generators[*gen-1].discard(seeds[1] + 1000000000*size_t(seeds[2]));
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

