#ifndef _Param_h
#define _Param_h

namespace Param{
  	const int factor = 2; //multiply standard chain length by this factor
  	const int n_cell_max = 5;
  	const int outcount = 200000; //frequency of counting of not found constraints
  	const bool ccount = true; //count constraints
  	const double box[2] = {1.0, 1.0};    //Box's size
  	const double eta = 0.708;    //density (volume fraction)
  	const std::string constraint_file("../Data.run/constraints.dat");//file of length number_spheres, containt the constraints
  	const std::string position_file("../Data.run/configuration.dat");//initial positions of the particles
	const std::string info_string ("../Data.run/info.dat");

	extern long int number_spheres;
	extern int number_cell[2];
	extern int total_number_cell;
   	extern double lambda_0;
	extern double sigma;
	extern double cell_size[2];
};


#endif
