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

#ifndef _CHECKER
#define _CHECKER
#include "boost/multi_array.hpp"
#include <unordered_map>
#include "Param.h" //set all parameters here
using namespace std;

class Checker{
  private:
    long int number_liftings;
    clock_t time_start,time_end;
    int cell_cur;
    int k_cur;

    vector<int> cell_ocp;  //# of disks for each cell
    boost::multi_array<int,2> P; //read-in constraints
    unordered_map<long int,long int> M; //bag for found constraints
    boost::multi_array<double,2> positions; //disk positions
    boost::multi_array<double,3> cell; //
    boost::multi_array<int, 2> plist; //added array to follow particle number during cell swaps
    boost::multi_array<int, 2> cell_neighbour;
    vector<int> count_visits;

    void countzero(long int l,ofstream &);
    void init_cell();
    void init_pos();
    void refresh_cell(double dx);
    void getconstraint(   ) ;
    void expl_cell( double &L_min, int &k_min, int &cell_min);
  
  public:
    Checker();
    ~Checker();
    void run();
};
#endif
