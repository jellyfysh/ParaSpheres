#ifndef _CHECKER
#define _CHECKER
#include "boost/multi_array.hpp"
#include <unordered_map>
#include "Param.h" //set all parameters here
using namespace std;

class Checker{
  private:
    long int number_liftings;
    clock_t time_start,time_end;
    int cell_cur;
    int k_cur;

    vector<int> cell_ocp;  //# of disks for each cell
    boost::multi_array<int,2> P; //read-in constraints
    unordered_map<long int,long int> M; //bag for found constraints
    boost::multi_array<double,2> positions; //disk positions
    boost::multi_array<double,3> cell; //
    boost::multi_array<int, 2> plist; //added array to follow particle number during cell swaps
    boost::multi_array<int, 2> cell_neighbour;
    vector<int> count_visits;

    void countzero(long int l,ofstream &);
    void init_cell();
    void init_pos();
    void refresh_cell(double dx);
    void getconstraint(   ) ;
    void expl_cell( double &L_min, int &k_min, int &cell_min);
  
  public:
    Checker();
    ~Checker();
    void run();
};
#endif
