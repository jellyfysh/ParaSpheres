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

#include "Simul.h"
#include <iostream>
#include <chrono>

using namespace std;
using namespace chrono;

int Simul::single_thread_run(){
    int number_liftings=0;
    double distance = Param::benchmark_chain_length;
    Node* j = spheres;
    double tau_tilde;
    cout << "Running single thread code" << endl;
  
    while (distance > 0.0) {
        Node *i, *j_tilde;
        number_liftings++;
        i = j;
        double tau = distance;
        for (int j_index = 0; j_index < i->number_arrows; j_index++) {
            j_tilde = i->arrow[j_index];
            tau_tilde = j_tilde->x - i->x - i->b[j_index];
            if (tau_tilde > 1) {
                tau_tilde -= 1;
            } else if (tau_tilde < 0) {
                tau_tilde += 1;
            }
            if (tau_tilde < tau) {
                j = j_tilde;
                tau = tau_tilde;
            }
        }
        i->x += tau;
        if (i->x > 0.5) i->x -= 1;
        distance -= tau;
    }
  cout.precision(11);
  cout << "Position of final active particle = (" << j->x << ", " << j->y << ")\n";
  number_liftings--;
  return number_liftings;
}


void Simul::single_timing(){
    cout.precision(11);
    cout << "Position of initial active particle = (" << spheres->x << ", " << spheres->y << ")\n";
    auto time_point = chrono::high_resolution_clock::now();
    int number_liftings = single_thread_run();

    cout.precision(5);
    auto duration = chrono::high_resolution_clock::now() - time_point;
    auto usec = duration_cast<chrono::microseconds>(duration).count();
    cout << "usec single thread " << usec << endl; //consumed time in microseconds
    cout << "Estimated number of events per hour: " << number_liftings * 1.0e+6 / usec * 3600 << endl;
}
