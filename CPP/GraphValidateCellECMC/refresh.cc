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
using namespace std;

void
Checker::refresh_cell(double dx){
//******************************************************************!
// 
// Compute the new cell of k_cur
//
//******************************************************************!
    int cell_new;
    if (dx > Param::cell_size[0]/2) {                                                //disk k_cur is now in cell_new=cell_neighbour(6,cell_cur)
        dx -= Param::cell_size[0];                                                   //position in the new cell
        cell_new = cell_neighbour[cell_cur][5];                                  //new cell
        cell_ocp[cell_new] = cell_ocp[cell_new] + 1;                          //refresh # of disk in the new cell

        cell[cell_new][cell_ocp[cell_new] - 1][0] = dx;                       //x position in the new cell
        cell[cell_new][cell_ocp[cell_new] - 1][1] = cell[cell_cur][k_cur][1]; //y position in the new cell
        plist[cell_new][  cell_ocp[cell_new] - 1] = plist[cell_cur] [k_cur] ;

        cell[cell_cur][k_cur][0] = cell[cell_cur][cell_ocp[cell_cur] - 1][0]; //swap the last disk in the old cell to take k_cur's vacant place
        cell[cell_cur][k_cur][1] = cell[cell_cur][cell_ocp[cell_cur] - 1][1];
        plist[cell_cur][k_cur] = plist [cell_cur][cell_ocp[cell_cur] - 1] ;
        plist [cell_cur][cell_ocp[cell_cur] - 1] =-1;

        int cache = cell_cur;
        cell_ocp[cell_cur] = cell_ocp[cell_cur] - 1;                          //refresh # of disks in the old cell
        cell_cur = cell_new;                                                  //The cell of the "current" disk is now cell_new....
        if(0)cout<<"Old new "<< plist[cache] [k_cur] <<" ";
        k_cur = cell_ocp[cell_new] - 1;                                       //and its place is the last one
        if(0)cout<< plist[cell_cur] [k_cur]<<endl;

        if(0){
            cout<<"Swap "<<endl;
            cout<<"plist out\t"<<cache<<"\t"<<cell_new<<endl;
            for(int i=0;i<5;i++){
	            cout<<plist[cache][i]<<"\t"<<plist[cell_new][i]<<endl;
            } cout<<endl;
        }
    } else {                                                                  //disk k_cur stays in cell_cur
        cell[cell_cur][k_cur][0] = dx;
    }
}
