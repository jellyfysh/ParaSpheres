#include "Simul.h"
#include <omp.h>
#include <string>
#include <iostream>

using namespace std;

Simul::Simul(int extern_argu, int extern_num_active){
    number_threads = omp_get_max_threads();
    number_actives = extern_num_active;
    all_liftings.resize(number_threads);
    init_mode = extern_argu;

    number_spheres = Param::number_spheres;
    cout << "Number disks " << number_spheres << endl;
    spheres = new Node[number_spheres];
    init();
    cout << "Number threads\t " << number_threads << endl;
}

Simul::~Simul(){
    delete [] spheres;
}
