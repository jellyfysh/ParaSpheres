// C++ version of the benchmark
// The algorithm is identical to the fortran code
// check whether every lifting take place on the constraint graph

# include <iostream>
# include <fstream>
# include <cmath>
# include <string>
#include <cstring>
#include <cassert>

#include "Checker.h"
using namespace std;

namespace Param {
	long int number_spheres = 0;
    int number_cell[2] = {0, 0};
    int total_number_cell = number_cell[0]*number_cell[1];
    double lambda_0 = 0;
    double sigma = 0;
    double cell_size[2] = {Param::box[0]/number_cell[0], Param::box[1]/number_cell[1]};
}

int main() {
  cout<<"Open\t"<<Param::info_string<<endl;
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
        }
    }
    Param::number_cell[0] = (int)(sqrt(Param::number_spheres) * 7.0 / 8.0);
    Param::number_cell[1]=(int)(sqrt(Param::number_spheres) * 7.0 / 8.0);
    Param::total_number_cell = Param::number_cell[0]*Param::number_cell[1];

    Param::lambda_0 = 0.07680/sqrt(Param::number_spheres);
    //Param::sigma = sqrt(Param::eta * Param::box[1]*Param::box[0]/((double)Param::number_spheres * M_PI));
    Param::cell_size[0] = Param::box[0]/Param::number_cell[0];
    Param::cell_size[1] = Param::box[1]/Param::number_cell[1];

	cout << Param::number_spheres << endl;
	cout<<"Checker"<<endl;
    Checker c;

    c.run();

    return 0;
}


