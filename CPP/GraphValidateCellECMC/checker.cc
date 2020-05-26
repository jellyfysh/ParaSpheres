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

#include "Checker.h"
# include <iostream>
# include <fstream>
using namespace std;

Checker::Checker(){
    cout<<"Check construct"<<endl;
    cout.precision(11);
    number_liftings=0;
    cell_cur=0;
    k_cur=0;
    count_visits.resize(Param::number_spheres, 0.);
    P.resize(boost::extents[Param::number_spheres][3]);
    cell_ocp.resize(Param::total_number_cell,0.);

    cell_neighbour.resize (boost::extents[Param::total_number_cell][9]);
    cell.resize(boost::extents[Param::total_number_cell][Param::n_cell_max][2]);
    plist.resize (boost::extents[Param::total_number_cell][Param::n_cell_max]);
    positions.resize (boost::extents[Param::number_spheres][2]);
  
    cout << "Number of particles = " << Param::number_spheres << "\n";
    cout << "Radius = " << Param::sigma << "\n";
    cout << "Box size = (" << Param::box[0] << ", " << Param::box[1] << ")\n";
    cout << "Number of cells = " << Param::number_cell[0] << ", " << Param::number_cell[1] << "\n";

    if(Param::ccount) getconstraint( );
    init_cell();
    cout << "Position of initial active particle = (" << cell[cell_cur][k_cur][0] - Param::box[0]/2 + (cell_cur % Param::number_cell[0] + 0.5)*Param::cell_size[0] <<
    ", " << cell[cell_cur][k_cur][1] - Param::box[1]/2 + (floor(cell_cur / Param::number_cell[1]) + 0.5)*Param::cell_size[1] << ")\n";
  
}

Checker::~Checker(){
    cout.precision(5);
    cout << number_liftings << " Liftings, " << double(time_end - time_start) / (double)CLOCKS_PER_SEC << " seconds\n";
    cout << "Estimated number of events per hour: " << (double)number_liftings * (double)CLOCKS_PER_SEC / double(time_end - time_start) * 3600. << endl;
    cout.precision(11);
    cout << "Position of the final active particle = (" << cell[cell_cur][k_cur][0] - Param::box[0]/2 + (cell_cur % Param::number_cell[0] + 0.5)*Param::cell_size[0] <<
    ", " << cell[cell_cur][k_cur][1] - Param::box[1]/2 + (floor(cell_cur / Param::number_cell[1]) + 0.5)*Param::cell_size[1] << ")\n";
  
    if(Param::ccount == 1){
        ofstream con1("foundConstraints.dat");
        ofstream con2("unfoundConstraints.dat");

        for (auto& x: M) {
	        long int lz= Param::number_spheres * Param::number_spheres;
            if(x.second !=0 ) con1 << x.first/(lz) <<"\t"<<x.first%(lz)<< " \t" << x.second << std::endl;
            if(x.second ==0 ) con2 << x.first/(lz) <<"\t"<<x.first%(lz)<< " \t" << x.second << std::endl;
        }
    }
    ofstream cnt("NumberLiftingsEachSphere.dat");
    for(int i=0;i<Param::number_spheres; i++){
        cnt<<i<<"\t"<<count_visits[i]<<endl;
    }
}
