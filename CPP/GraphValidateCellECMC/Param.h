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
//


#ifndef _Param_h
#define _Param_h

namespace Param{
  	const int factor = 2; //multiply standard chain length by this factor
  	const int n_cell_max = 5;
  	const int outcount = 200000; //frequency of counting of not found constraints
  	const bool ccount = true; //count constraints
  	const double box[2] = {1.0, 1.0};    //Box's size
  	const double eta = 0.708;    //density (volume fraction)
  	const std::string constraint_file("../Data.run/constraints.dat");//file of length number_spheres, containt the constraints
  	const std::string position_file("../Data.run/configuration.dat");//initial positions of the particles
	const std::string info_string ("../Data.run/info.dat");

	extern long int number_spheres;
	extern int number_cell[2];
	extern int total_number_cell;
   	extern double lambda_0;
	extern double sigma;
	extern double cell_size[2];
};


#endif
