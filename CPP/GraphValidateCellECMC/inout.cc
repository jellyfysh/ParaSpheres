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

#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include "Checker.h"

using namespace std;

void
Checker::getconstraint() {
    int np = 0;
    cout<<"Open\t "<< Param::constraint_file <<endl;
    ifstream constraint (Param::constraint_file);
    assert( constraint.good() );
    if (constraint.is_open()) {
        for (string line; getline(constraint, line); np++) {
            istringstream in(line);
            int ic=0;
            int l;
            while( in >> l ){
	            P[np][ic]=l;
	            M[max(l,np) * Param::number_spheres * Param::number_spheres + min(l, np)]=0;//initialize bag for counting constraints
	            ic++;
            }
        }
        assert(np == Param::number_spheres);
        cout << "Reads " << np << " points for constraint graph.\n";
    }
}

void Checker::countzero(long int l,ofstream & os){
    long int nzero=0;
    for (auto& x: M) {
        if(x.second == 0 ) nzero++;
    }
    os << l << "\t" << nzero << endl;
}
