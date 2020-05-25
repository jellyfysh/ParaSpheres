#include "Simul.h"
#include <iostream>
#include <cstring>
# include <fstream>
# include <sstream>
#include <cassert>

using namespace std;

namespace Param{
    long int number_spheres= 256 * 256;
    double sigma = 0;
    int logkmax=0;
    double validation_chain_length=0;
}

int main(int argc, char* argv[]){
    ifstream info  (Param::info_string);
    assert(info.good());
    if (info.is_open()) {
	    for (string line; getline(info, line); ) {
            istringstream in(line);
            string label_str;
            string value;
            in >> label_str >> value;
            const char *label = label_str.c_str();
	        if (strcmp(label, "N") == 0) Param::number_spheres = stoi(value);
	        if (strcmp(label, "sigma") == 0) Param::sigma = stod(value);
	        if (strcmp(label, "LogKmax") == 0) Param::logkmax = stoi(value);
	        if (strcmp(label, "ChainLength") == 0) Param::validation_chain_length = stod(value);
        }
    }

    if (argc == 3) {
        if (pow(2, Param::logkmax) < stoi(argv[2])) cerr << "Too many active spheres" << endl;
        if (strcmp(argv[1], "validation") == 0) { //do validation
            Simul s (0, stoi(argv[2]));
            cout << "Validation starting from a random time.\n";
            s.multi_thread_validate(Param::validation_chain_length / stoi(argv[2]));
        } else if (strcmp(argv[1], "timing") == 0) { //do benchmark
            Simul s (1, stoi(argv[2]));
            cout << "Timing\n";
            s.multi_thread_run(Param::benchmark_chain_length / stoi(argv[2]), false);
        } else if (strcmp(argv[1], "comparison") == 0) {
            Simul s2 (2, 1);
            Simul s (2, 1);
            cout << "Comparison between atomic and no atomic, single thread.\n";
            s2.single_timing();
            s.multi_thread_run(Param::benchmark_chain_length, false);
        } else {
            cerr << "Wrong argument.\n";
        }
    } else {
        cerr << "wrong number of arguments.\n";
    }
    return 0;
}
