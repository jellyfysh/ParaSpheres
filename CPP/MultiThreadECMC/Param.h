// Parameters
// all parameters will be read in, except for benchmark_chain_length
#ifndef _PARAM
#define _PARAM
#include <string>
#include <cmath>

namespace Param{
    const double benchmark_chain_length = 15000;//sum of chain length of all active spheres during benchmark
    static std::string active_string ("../Data.run/active"); //active sphere list
    static std::string ref_string ("../Data.run/lifting"); //reference lifting list
    const std::string constraint_string ("../Data.run/constraints.dat"); //constraint graph
    const std::string init_string ("../Data.run/configuration.dat"); //configuration
    const std::string info_string ("../Data.run/info.dat"); // information about each run

	extern long int number_spheres; //N
	extern double sigma;
	extern int logkmax; //log_2(maximum number of active spheres)
	extern double validation_chain_length;//sum of chain length of all active spheres during validation
};

#endif
