// JeLLyFysh/ParaSpheres - Multithreaded event-chain Monte Carlo with local times
// - https://github.com/jellyfysh/paraspheres
// Copyright (C) 2020 The JeLLyFysh organization
// (see the AUTHORS file on jellyfysh/paraspheres for the full list of authors)
//
// This file is part of JeLLyFysh/ParaSpheres.
//
// JeLLyFysh/ParaSpheres is free software: you can redistribute it and/or modify
// it under the terms of the GNU General
// Public License as published by the Free Software Foundation, either > version
// 3 of the License, or (at your option)
// any later version.
//
// JeLLyFysh/ParaSpheres is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even 
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// JeLLyFysh/ParaSpheres in the LICENSE file. If not, see
// <https://www.gnu.org/licenses/>.
//
// If you use JeLLyFysh/ParaSpheres in published work, please cite the following
// reference
// (see [Li2020] in References.bib):
// Botao Li, Synge Todo, A. C. Maggs, Werner Krauth
// Multithreaded event-chain Monte Carlo with local times,
// arXiv e-prints: 2004.11040 (2020), https://arxiv.org/abs/2004.11040

// Parameters
// all parameters will be read in, except for benchmark_chain_length
#ifndef _PARAM
#define _PARAM
#include <string>
#include <cmath>

namespace Param{
    const double benchmark_chain_length = 15000;//sum of chain length of all active spheres during benchmark
    static std::string active_string ("../Data.run/active"); //active sphere list
    static std::string ref_string ("../Data.run/lifting"); //reference lifting list
    const std::string constraint_string ("../Data.run/constraints.dat"); //constraint graph
    const std::string init_string ("../Data.run/configuration.dat"); //configuration
    const std::string info_string ("../Data.run/info.dat"); // information about each run

	extern long int number_spheres; //N
	extern double sigma;
	extern int logkmax; //log_2(maximum number of active spheres)
	extern double validation_chain_length;//sum of chain length of all active spheres during validation
};

#endif
