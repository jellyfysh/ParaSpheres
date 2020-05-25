#include "Checker.h"
# include <iostream>
# include <fstream>
using namespace std;

Checker::Checker(){
    cout<<"Check construct"<<endl;
    cout.precision(11);
    number_liftings=0;
    cell_cur=0;
    k_cur=0;
    count_visits.resize(Param::number_spheres, 0.);
    P.resize(boost::extents[Param::number_spheres][3]);
    cell_ocp.resize(Param::total_number_cell,0.);

    cell_neighbour.resize (boost::extents[Param::total_number_cell][9]);
    cell.resize(boost::extents[Param::total_number_cell][Param::n_cell_max][2]);
    plist.resize (boost::extents[Param::total_number_cell][Param::n_cell_max]);
    positions.resize (boost::extents[Param::number_spheres][2]);
  
    cout << "Number of particles = " << Param::number_spheres << "\n";
    cout << "Radius = " << Param::sigma << "\n";
    cout << "Box size = (" << Param::box[0] << ", " << Param::box[1] << ")\n";
    cout << "Number of cells = " << Param::number_cell[0] << ", " << Param::number_cell[1] << "\n";

    if(Param::ccount) getconstraint( );
    init_cell();
    cout << "Position of initial active particle = (" << cell[cell_cur][k_cur][0] - Param::box[0]/2 + (cell_cur % Param::number_cell[0] + 0.5)*Param::cell_size[0] <<
    ", " << cell[cell_cur][k_cur][1] - Param::box[1]/2 + (floor(cell_cur / Param::number_cell[1]) + 0.5)*Param::cell_size[1] << ")\n";
  
}

Checker::~Checker(){
    cout.precision(5);
    cout << number_liftings << " Liftings, " << double(time_end - time_start) / (double)CLOCKS_PER_SEC << " seconds\n";
    cout << "Estimated number of events per hour: " << (double)number_liftings * (double)CLOCKS_PER_SEC / double(time_end - time_start) * 3600. << endl;
    cout.precision(11);
    cout << "Position of the final active particle = (" << cell[cell_cur][k_cur][0] - Param::box[0]/2 + (cell_cur % Param::number_cell[0] + 0.5)*Param::cell_size[0] <<
    ", " << cell[cell_cur][k_cur][1] - Param::box[1]/2 + (floor(cell_cur / Param::number_cell[1]) + 0.5)*Param::cell_size[1] << ")\n";
  
    if(Param::ccount == 1){
        ofstream con1("foundConstraints.dat");
        ofstream con2("unfoundConstraints.dat");

        for (auto& x: M) {
	        long int lz= Param::number_spheres * Param::number_spheres;
            if(x.second !=0 ) con1 << x.first/(lz) <<"\t"<<x.first%(lz)<< " \t" << x.second << std::endl;
            if(x.second ==0 ) con2 << x.first/(lz) <<"\t"<<x.first%(lz)<< " \t" << x.second << std::endl;
        }
    }
    ofstream cnt("NumberLiftingsEachSphere.dat");
    for(int i=0;i<Param::number_spheres; i++){
        cnt<<i<<"\t"<<count_visits[i]<<endl;
    }
}
