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
# include <iostream>
# include <fstream>
void Checker::run(){
  
    int k_first_event, cell_first_event;
    double displacement_max =
    min(min(Param::box[0], Param::box[1])/2, min(Param::cell_size[0], Param::cell_size[1])) - 2*Param::sigma;
    ofstream oszero ("NumberUnvisitedConstraints.dat");
    double distance_to_go = 50000000*Param::lambda_0 *Param::factor; //standard run length

    time_start = clock();
  
    while (true) {   //loop over collisions
        double displacement_first_event = distance_to_go;
        expl_cell( displacement_first_event, k_first_event, cell_first_event);

        displacement_first_event = max(displacement_first_event, .0);
        double new_position =
                cell[cell_cur][k_cur][0] + min(min(displacement_first_event, distance_to_go), displacement_max);

        refresh_cell( new_position );

        if (displacement_max < min(displacement_first_event,distance_to_go)) {
            distance_to_go -= displacement_max;
            continue; //move limited by cell size, the move continues with the same disk and the same direction
        }
        else if (displacement_first_event < distance_to_go) {
            distance_to_go = distance_to_go - displacement_first_event;
            k_cur = k_first_event;
            cell_cur = cell_first_event;
            number_liftings++;
            if(number_liftings % Param::outcount == 0) countzero(number_liftings, oszero);
            continue; //a regular lifting move
        }
        else {
            break; //end of the chain
        }
    }
    time_end = clock();
}
