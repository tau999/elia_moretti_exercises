#ifndef __GENETIC__
#define __GENETIC__

#include "random.h"
#include <vector>

using namespace std;

int seed[4];
Random rnd;




// PARAMETRI 
const int n_cities=32;
int n_iter=5000000;
double T=1;
int type_cities = 0;  // 0 -> CIRCLE    &      1 -> SQUARE

double attempted = 0;
double accepted = 0;

vector<vector<double>> cities;


// STRUCT
struct Population {
    int route[n_cities]={0};
    double fitness=0;
};

Population current_route;
Population new_route;




#endif
