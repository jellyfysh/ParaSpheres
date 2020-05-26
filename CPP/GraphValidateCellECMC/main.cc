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

// C++ version of the benchmark
// The algorithm is identical to the fortran code
// check whether every lifting take place on the constraint graph

# include <iostream>
# include <fstream>
# include <cmath>
# include <string>
#include <cstring>
#include <cassert>

#include "Checker.h"
using namespace std;

namespace Param {
	long int number_spheres = 0;
    int number_cell[2] = {0, 0};
    int total_number_cell = number_cell[0]*number_cell[1];
    double lambda_0 = 0;
    double sigma = 0;
    double cell_size[2] = {Param::box[0]/number_cell[0], Param::box[1]/number_cell[1]};
}

int main() {
  cout<<"Open\t"<<Param::info_string<<endl;
    ifstream info  (Param::info_string);
    assert(info.good());
    if (info.is_open()) {
        for (string line; getline(info, line); ) {
            istringstream in(line);
            string label_str;
            string value;
            in >> label_str >> value;
            const char *label = label_str.c_str();
            if (strcmp(label, "N") == 0) Param::number_spheres = stoi(value);
            if (strcmp(label, "sigma") == 0) Param::sigma = stod(value);
        }
    }
    Param::number_cell[0] = (int)(sqrt(Param::number_spheres) * 7.0 / 8.0);
    Param::number_cell[1]=(int)(sqrt(Param::number_spheres) * 7.0 / 8.0);
    Param::total_number_cell = Param::number_cell[0]*Param::number_cell[1];

    Param::lambda_0 = 0.07680/sqrt(Param::number_spheres);
    //Param::sigma = sqrt(Param::eta * Param::box[1]*Param::box[0]/((double)Param::number_spheres * M_PI));
    Param::cell_size[0] = Param::box[0]/Param::number_cell[0];
    Param::cell_size[1] = Param::box[1]/Param::number_cell[1];

	cout << Param::number_spheres << endl;
	cout<<"Checker"<<endl;
    Checker c;

    c.run();

    return 0;
}


