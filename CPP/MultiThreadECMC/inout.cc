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

#include "Simul.h"
# include <fstream>
# include <sstream>
# include <iostream>
#include <cmath>
#include <cassert>

using namespace std;

void Simul::init() {
    int np = 0;
    cout << init_mode << endl;
	string lactive_string=Param::active_string;
	string lref_string=Param::ref_string;
    std::string number_active_string = to_string(number_actives);

    lactive_string += number_active_string;
    lref_string += number_active_string;
    lactive_string += ".dat";
    lref_string += ".dat";
    cout << lactive_string << endl;
    cout << lref_string << endl;
  
    ifstream init (Param::init_string);
    assert( init.good() );

    if (init.is_open()) {
        for (string line; getline(init, line); np++) {
            istringstream in(line);
            in >> spheres[np].x >> spheres[np].y;
        }
    
        cout << "Reads " << np << " points for configuration.\n";
    } else {
        cerr << "Error: no valid initial configuration.\n"; exit(0);
    }

    np = 0;
    ifstream constraint  (Param::constraint_string);
    ofstream of ("contact.dat");

    assert( constraint.good() );
    if (constraint.is_open()) {
        for (string line; getline(constraint, line); np++) {
            istringstream in(line);
            int ic = 0;
            int l;
            while(in >> l ){
                spheres[np].number_arrows = ic;
	            spheres[np].number_arrows++;
                spheres[np].arrow[ic] = &spheres[l] ;
                spheres[np].b[ic] = contact(np, l);
	            of << spheres[np].b[ic] << "\t";
	            ic++;
            }
            of<<endl;
        }
        cout << "Reads " << np << " points for constraint graph.\n";
    } else {
        cerr << "Error: no valid constraint graph.\n"; exit(0);
    }

	cout<<"active string "<<lactive_string<<endl;
    ifstream act (lactive_string);
    if (init_mode != 1){
        assert(act.good() );
        if( act.is_open() ){
            for (string line; getline(act,line); ) {
	        istringstream in(line);
	        int l;
	            while(in>>l) {
	                spheres[l].tag = Node::tag_active;
	                active_spheres.push_back(spheres + l);
	            }
            }
        } else {
            cerr<<"error reading active particles"<<endl; exit(0);
        }
    } else {
        cout << "Generate initial active particles.\n";
        int distance = (int)number_spheres / number_actives;
        for (int i=0; i<number_actives; i++){
            spheres[distance * i].tag = Node::tag_active;
	        active_spheres.push_back(spheres + distance * i);
        }
    }
    cout << "Number of active particles: " << active_spheres.size() << endl;
  
    if (init_mode == 0) {
        ifstream ben (lref_string);
        int cnt=0;
        assert( ben.good() );
        if( ben.is_open() ){
            for(string line; getline(ben,line); ){
	            istringstream in(line);
	            double f; int i, j;
	            cnt++;
	            in>>f;      in>>i;      in>>j;
	            reference_liftings.emplace_back (make_tuple(f, i, j) );
            }
        } else {
            cerr<<"no reference_liftings file"<<endl; exit(0);
        }
    }
}

double Simul::contact(int particle1, int particle2) {
    double distance = remainder(spheres[particle2].y - spheres[particle1].y, 1.0);
    return sqrt(4.0 * Param::sigma*Param::sigma - distance*distance );
}
