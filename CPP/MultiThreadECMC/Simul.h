#ifndef _Simul
#define _Simul

#include "Param.h"
#include <vector>
#include <string>
# include <ctime>
#include <stack>
#include <unordered_set>
#include <tuple>        // std::tuple, std::get, std::tie, std::ignore
#include <atomic>

using namespace std;

class Node {
  public:
    static const int tag_static = -3;
    static const int tag_stalled = -2;
    static const int tag_active = -1;
  //  Indeed, we have four types in the case where the number of active particles is larger than the number of threads...
  //moving active particles (tag = thread id)
  //waiting active particles (tag = _tag_active)
  //stalled particles (tag = tag_stalled)
  //static particles (tag = tag_static)

    Node* arrow[3];
    double b[3], x, y;
    double local_time;
    int number_arrows;
    std::atomic<int> tag;
  // -3 : static
  // -2 : stalled (already arrived at horizon or hitting such sphere)
  // -1 : active but not currently handled by any thread
  // >0 : active and handled by some thread (value represents thread id)
    Node() {local_time = 0; tag = tag_static;}
};

class Simul {
  private:
    int number_threads;
    long int number_spheres; //N
    int number_actives;

    Node *spheres;
    void init();
    double contact(int, int); //b_{ij}
    int init_mode;
    std::vector<Node*> active_spheres;
    vector<vector<tuple<double,int,int>>> all_liftings; //lift of liftings for each thread
    vector<tuple<double,int,int>> reference_liftings;
  public:
    void multi_thread_run(double, bool);
    int single_thread_run();
    void single_timing();
    void advance(double); // advance the configuration to a certain time using the reference lifrings
    void multi_thread_validate(double);
    Simul(int, int);
    ~Simul();
};

#endif
