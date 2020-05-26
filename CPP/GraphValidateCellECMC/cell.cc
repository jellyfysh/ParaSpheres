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

#include "Checker.h"
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unistd.h>

using namespace std;

void Checker::expl_cell( double &L_min, int &k_min, int &cell_min){
  //******************************************************************!
  // 
  // Explore the disks neighbor of k_cur, 6 cell-neighb to explore (2,3,5,6,8,9)
  // Faster without the loop
  //
  //******************************************************************!
  double dx[2];
  double x_cur[2] = {cell[cell_cur][k_cur][0], cell[cell_cur][k_cur][1]};
  int cell_act;
  double L_coll;
  int  j;
  if(0)	cout<<"\nexpl\n";
  int lcase=-1;
  double ox,oy;
  cell_act = cell_neighbour[cell_cur][5];
  for (j = 0; j < cell_ocp[cell_act]; j++) {
    dx[1] = cell[cell_act][j][1] - x_cur[1];
    if (abs(dx[1]) < 2*Param::sigma) {
      dx[0] = cell[cell_act][j][0] - x_cur[0] + Param::cell_size[0];
      L_coll = dx[0] - sqrt(4 *(Param::sigma*Param::sigma) - (dx[1] * dx[1]));
      if (L_coll < L_min) {
	L_min = L_coll;
	k_min = j;
	cell_min = cell_act;
	if(0) cout<<"case5"<<endl;
        lcase=5;
        ox= x_cur[0] ;
        oy=x_cur[1];
      }
    }
  }
  cell_act = cell_neighbour[cell_cur][1];
  for (j = 0; j < cell_ocp[cell_act]; j++) {
    dx[1] =  x_cur[1] + Param::cell_size[1] - cell[cell_act][j][1];
    if (abs(dx[1]) < 2*Param::sigma) {
      dx[0] = cell[cell_act][j][0] - x_cur[0];
      if (dx[0] > 1.e-16 ) {
	L_coll = dx[0] - sqrt(4 *(Param::sigma*Param::sigma) - (dx[1] * dx[1]));
	if (L_coll < L_min) {
	  L_min = L_coll;
	  k_min = j;
	  cell_min = cell_act;
          if(0) cout<<"case1"<<endl;
	  lcase=1;
          ox= x_cur[0] ;
          oy=x_cur[1];
        }
      }
    }
  }
  cell_act = cell_neighbour[cell_cur][7];
  for (j = 0; j < cell_ocp[cell_act]; j++) {
    dx[1] = cell[cell_act][j][1] - x_cur[1] + Param::cell_size[1];
    if (abs(dx[1]) < 2*Param::sigma) {
      dx[0] = cell[cell_act][j][0] - x_cur[0];
      if (dx[0] > 1.e-16) {
	L_coll = dx[0] - sqrt(4 *(Param::sigma*Param::sigma) - (dx[1] * dx[1]));
	if (L_coll < L_min) {
	  L_min = L_coll;
	  k_min = j;
	  cell_min = cell_act;
          if(0) cout<<"case7"<<endl;
	  lcase=7;
          ox= x_cur[0] ;
          oy=x_cur[1];
        }
      }
    }
  }
  cell_act = cell_neighbour[cell_cur][4];
  for (j = 0; j < cell_ocp[cell_act]; j++) {
    dx[1] = cell[cell_act][j][1] - x_cur[1];
    if (abs(dx[1]) < 2*Param::sigma) {
      dx[0] = cell[cell_act][j][0] - x_cur[0];
      if (dx[0] > 1.e-14) {
	L_coll = dx[0] - sqrt(4 *(Param::sigma*Param::sigma) - (dx[1]*dx[1]));
	if (L_coll < L_min) {
	  L_min = L_coll;
	  k_min = j;
	  cell_min = cell_act;
          if(0) cout<<"case4"<<endl;
	  lcase=4;
          ox= x_cur[0] ;
          oy=x_cur[1];
        }
      }
    }
  }
  cell_act = cell_neighbour[cell_cur][2];
  for (j = 0; j < cell_ocp[cell_act]; j++) {
    dx[1] = x_cur[1] - cell[cell_act][j][1] + Param::cell_size[1];
    if (abs(dx[1]) < 2*Param::sigma) {
      dx[0] = cell[cell_act][j][0] - x_cur[0] + Param::cell_size[0];
      L_coll = dx[0] - sqrt(4 *(Param::sigma*Param::sigma) - (dx[1]*dx[1]));
      if (L_coll < L_min) {
	L_min = L_coll;
	k_min = j;
	cell_min = cell_act;
        if(0) cout<<"case2"<<endl;
	lcase=2;
        ox= x_cur[0] ;
        oy=x_cur[1];
      }
    }
  }
  cell_act = cell_neighbour[cell_cur][8];
  for (j = 0; j < cell_ocp[cell_act]; j++) {
    dx[1] = cell[cell_act][j][1] - x_cur[1] + Param::cell_size[1];
    if (abs(dx[1]) < 2*Param::sigma) {
      dx[0] = cell[cell_act][j][0] - x_cur[0] + Param::cell_size[0];
      L_coll = dx[0] - sqrt(4 *(Param::sigma*Param::sigma) - (dx[1] * dx[1]));
      if (L_coll < L_min) {
	L_min = L_coll;
	k_min = j;
	cell_min = cell_act;
	if(0) cout<<"case8"<<endl;
	lcase=8;
        ox= x_cur[0] ;
        oy=x_cur[1];
      }
    }
  }

  if(Param::ccount){
    int first= plist [cell_cur] [k_cur] ;
    int second = plist [cell_min][k_min];
    bool logic=P[first][0]==second || P[first][1]==second || P[first][2]==second || first==second;
    //lcase -1 corresponds to empty cell
    if (!logic && lcase != -1){
      cout<<"first "<<cell_cur<<" "<<k_cur<<"\t"<< cell_min<<" "<<k_min<<endl;
      cout<<"Lcol "<<L_coll<<endl;
      cout<<"Lmin "<<L_min<<endl;
      cout << "P "<< P[first][0]<<" "<<P[first][1]<<" "<<P[first][2]<<endl;
      cout << "P1 P2 "<<first<<" "<<second<<endl;
      cout<< "dx\t"<<dx[0]<<"\t"<<dx[1]<<endl;
      cout<<"found particle\t"<<cell[cell_min][k_min][0]<<"\t"<<cell[cell_min][k_min][1]<<endl;
      cout <<"case="<<lcase<<endl;
      cout<<"ox,oy\t"<<ox<<"\t"<<oy<<endl;
      cout<<"cell size\t"<<Param::cell_size[0]<<"\t"<<Param::cell_size[1]<<endl;
      cout<<endl;
      cout<<"The program has found a problem with the constrain graph"<<endl;
      assert( logic );
    }

    if( lcase != -1 ) {
      count_visits[second]++;
      if(first != second)  M[(long int)(max(first,second)) * Param::number_spheres * Param::number_spheres + min(first, second)] +=1;
    }
  }
}
