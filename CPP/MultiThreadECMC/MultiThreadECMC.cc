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
#include <limits>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <chrono>
#include <omp.h>
#include <thread>

// https://en.cppreference.com/w/cpp/atomic/memory_order
//https://preshing.com/20140709/the-purpose-of-memory_order_consume-in-cpp11/

#define AcqRel 1

#ifdef AcqRel
const memory_order load_order=memory_order_acquire;
const memory_order store_order=memory_order_release;
const memory_order cas_order=memory_order_seq_cst;
#endif 


#ifdef Acquire
const memory_order load_order=memory_order_acquire;
const memory_order store_order=memory_order_relaxed;
const memory_order cas_order=memory_order_seq_cst;
#endif


#ifdef Release
const memory_order load_order=memory_order_relaxed;
const memory_order store_order=memory_order_release;
const memory_order cas_order=memory_order_seq_cst;
#endif


#ifdef Relax
const memory_order load_order=memory_order_relaxed;
const memory_order store_order=memory_order_relaxed;
const memory_order cas_order=memory_order_relaxed;
#endif

#ifdef SeqCst
const memory_order load_order=memory_order_seq_cst;
const memory_order store_order=memory_order_seq_cst;
const memory_order cas_order=memory_order_seq_cst;
#endif

inline bool comp(const std::tuple<double, int, int>& a, const std::tuple<double, int,int>& b)
        { return (get<0>(a) < get<0>(b)); }
//static double epsilon = 0.0000000001;
static double epsilon = 0.0000001;

// Algorithm 6
// Correspondence to algorithm 5 is noted in the comments
// Some operations, related solely to local variables, are different from their conterparts in algorithm 5

void Simul::multi_thread_run(double new_breakpoint, bool bMakeList) {
    double earliest_horizon_violation = std::numeric_limits<double>::max();
    double lifting_divergence = std::numeric_limits<double>::max();
    cout << "Number of threads: " << number_threads << endl;
    long int number_liftings[number_threads];
    for(int i = 0; i < number_threads; i++) number_liftings[i] = 0;
    auto wall_time_point = chrono::high_resolution_clock::now();
    #pragma omp parallel for reduction(min : earliest_horizon_violation)
    for (int iota = 0; iota < number_threads; ++iota) {
        long int number_liftings_thread = 0;
        // prepare active_spheres for each thread
        std::vector<Node*> active_spheres_local, active_spheres_local_new;
        for (std::size_t i = iota; i < active_spheres.size(); i += number_threads)
            active_spheres_local.push_back(active_spheres[i]);
        while (!active_spheres_local.empty()) {
            for (auto i : active_spheres_local) {
                i->tag.store(iota, store_order);
                auto j = i; //1\iota
                double distance = new_breakpoint - i->local_time;
                double x = 0; //1\iota
                double tau = distance; //1\iota
                bool bStalled = false;
	            for (int j_index = 0; j_index < i->number_arrows; ++j_index) { //2\iota
                    auto j_tilde = i->arrow[j_index];
                    double x_j_tilde = j_tilde->x; //3\iota
                    double tau_ij_tilde = x_j_tilde - i->x; //4\iota
                    if (tau_ij_tilde < 0) tau_ij_tilde += 1;
                    tau_ij_tilde -= i->b[j_index];
                    if (j_tilde->local_time > i->local_time + tau_ij_tilde) //5\iota
                        earliest_horizon_violation =
                                std::min({earliest_horizon_violation, i->local_time, j_tilde->local_time});  //5\iota
                    if (tau_ij_tilde < tau) { //6\iota
                        j = j_tilde; //7\iota
                        x = x_j_tilde; //8\iota
                        tau = tau_ij_tilde; //9\iota
                    }
                }

	            if (j->tag.load(load_order) == Node::tag_stalled) { // if the target is stall, the active sphere also become stalled, and this is treated as a horizon violation
                        earliest_horizon_violation =
                                std::min({earliest_horizon_violation, i->local_time, j->local_time});
	                i->tag.store(Node::tag_stalled, store_order);
                    bStalled = true;
	            }
	            if (bStalled) continue;
                if (j != i) {
                    int expected = Node::tag_static;
                    j->tag.compare_exchange_strong(expected, iota, cas_order);  //10\iota
                }
                if (j->tag.load(load_order) == iota) {
                    if (tau < distance) {  //11\iota
                        if (x == j->x) {  //12\iota
	                        double lifting_time = i->local_time + tau;
	                        if (bMakeList) all_liftings[iota].emplace_back(make_tuple(lifting_time,
	                                i - spheres, j - spheres));
                            if (not bMakeList) number_liftings_thread++;
                            j->local_time = lifting_time; //14\iota
                            i->local_time = lifting_time; //original 15\iota
                            i->x += tau; //16\iota
                            if (i->x > 0.5) i->x -= 1;
                            //std::this_thread::sleep_for (std::chrono::microseconds(1));//necessary for seeing 'not OK'
	                        //i->local_time = lifting_time;//invert 15\iota 16\iota
                            if (j != i) i->tag.store(Node::tag_static, store_order); //17\iota
	                        j->tag.store(Node::tag_active, store_order);
                            active_spheres_local_new.push_back(j);
                        } else {
                            double lifting_time = i->local_time + tau; //fix for 'necklace'
                            i->local_time = lifting_time;
                            i->x += tau;
                            if (i->x > 0.5) i->x -= 1;
	                        i->tag.store(Node::tag_active, store_order);
                            if (j != i) j->tag.store(Node::tag_static, store_order); //20\iota
                            active_spheres_local_new.push_back(i);
                        }
                    } else {
	                    double temp_local_time = i->local_time;
                        i->local_time = new_breakpoint; //21\iota
                        i->x += new_breakpoint - temp_local_time; //22\iota
                        if (i->x > 0.5) i->x -= 1;
	                    i->tag.store(Node::tag_stalled, store_order); //24\iota
                    }
                } else {
                    double event_time = i->local_time + tau;
                    i->local_time = event_time;
                    i->x += tau;
                    if (i->x > 0.5) i->x -= 1;
	                i->tag.store(Node::tag_active, store_order);
                    active_spheres_local_new.push_back(i);
                }
            }
            swap(active_spheres_local, active_spheres_local_new);
            active_spheres_local_new.resize(0);
        }
        if (not bMakeList) number_liftings[iota] = number_liftings_thread;
    }

    if (bMakeList) {
        for (int iota = 1; iota < number_threads; ++iota)
            all_liftings[0].insert(all_liftings[0].end(), all_liftings[iota].begin(),
                    all_liftings[iota].end());
        std::sort(all_liftings[0].begin(), all_liftings[0].end(), comp);
        int l_index = 0;
        for (auto l = all_liftings[0].begin(); l != all_liftings[0].end(); ++l, ++l_index){
            if ((get<0>(*l) - get<0>(reference_liftings[l_index])) > epsilon ||
                (get<1>(*l) != get<1>(reference_liftings[l_index])) ||
                (get<2>(*l) != get<2>(reference_liftings[l_index]))) {
                cout << "divergence" << " " << get<0>(*l) << " " << get<1>(*l) << ' ' << get<2>(*l) << " " << get<1>(
                        reference_liftings[l_index]) << " " << get<2>(reference_liftings[l_index]) << endl;
                lifting_divergence = min(get<0>(*l), get<0>(reference_liftings[l_index]));
                break;
            }
        }
        if(std::numeric_limits<double>::max() != earliest_horizon_violation){
            cout << "earliest_horizon_violation = " << earliest_horizon_violation << endl;
            cout << "lifting_divergence = " << lifting_divergence << endl;
        }
        if (earliest_horizon_violation <= lifting_divergence + epsilon) {
            cout << "OK\n";
        } else {
	    cerr << earliest_horizon_violation << ' ' << lifting_divergence << endl;
            cout << "Not OK\n" <<endl;
        }
    } else {
        long int total_number_liftings = 0;
        for (int i = 0; i < number_threads; i++) total_number_liftings += number_liftings[i];
        auto duration = chrono::high_resolution_clock::now() - wall_time_point;
        auto usec = chrono::duration_cast<chrono::microseconds>(duration).count();
        cout << "EPH: " << 1.0 * 3600 * 1000000 / usec * (double)total_number_liftings << endl;
        cout << "Duration threaded usec \t" << usec << endl;
        cout << "Number of events: " << total_number_liftings << endl;
    }
}

void Simul::advance(double breakpoint) {
    double time = 0;
    unsigned long number_liftings = 0;
    while (get<0>(reference_liftings[number_liftings]) < breakpoint) {
        double time_of_flight = get<0>(reference_liftings[number_liftings]) - time;
        time = get<0>(reference_liftings[number_liftings]);
        for (auto i : active_spheres) {
            i->x += time_of_flight;
            if (i->x > 0.5) i->x -= 1;
        }
        (spheres + get<1>(reference_liftings[number_liftings]))->tag = Node::tag_static;
        (spheres + get<2>(reference_liftings[number_liftings]))->tag = Node::tag_active;
        active_spheres.erase(remove(active_spheres.begin(), active_spheres.end(),
                spheres + get<1>(reference_liftings[number_liftings])), active_spheres.end());
        active_spheres.push_back(spheres + get<2>(reference_liftings[number_liftings]));
        all_liftings[0].emplace_back(make_tuple(time, get<1>(reference_liftings[number_liftings]), get<2>(
                reference_liftings[number_liftings])));
        number_liftings++;
        if (number_liftings == reference_liftings.size()) break;
    }
    for (auto i : active_spheres) {
        i->x += breakpoint - time;
        if (i->x > 0.5) i->x -= 1;
    }
    for (int i=0; i < number_spheres; i++) (spheres + i)->local_time = breakpoint;
}

void Simul::multi_thread_validate(double new_breakpoint){
    //srand((unsigned)time(0));
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    auto t_int =chrono::duration_cast<chrono::nanoseconds>(t1.time_since_epoch()).count();
    srand((unsigned)t_int);
    double breakpoint = (double)rand() / (double)RAND_MAX * new_breakpoint;
    cout << "Advanced to: " << breakpoint << endl;
    advance(breakpoint);
    multi_thread_run(new_breakpoint, true);
}
