#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include "Checker.h"

using namespace std;

void
Checker::getconstraint() {
    int np = 0;
    cout<<"Open\t "<< Param::constraint_file <<endl;
    ifstream constraint (Param::constraint_file);
    assert( constraint.good() );
    if (constraint.is_open()) {
        for (string line; getline(constraint, line); np++) {
            istringstream in(line);
            int ic=0;
            int l;
            while( in >> l ){
	            P[np][ic]=l;
	            M[max(l,np) * Param::number_spheres * Param::number_spheres + min(l, np)]=0;//initialize bag for counting constraints
	            ic++;
            }
        }
        assert(np == Param::number_spheres);
        cout << "Reads " << np << " points for constraint graph.\n";
    }
}

void Checker::countzero(long int l,ofstream & os){
    long int nzero=0;
    for (auto& x: M) {
        if(x.second == 0 ) nzero++;
    }
    os << l << "\t" << nzero << endl;
}
