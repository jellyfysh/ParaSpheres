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
