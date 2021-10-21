#ifndef __GENETIC__
#define __GENETIC__

#include "random.h"
#include <vector>

using namespace std;

int seed[4];
Random rnd;




// PARAMETRI 
const int n_cities=30;
int n_iter=200;
int n_pop=1000;
double r_cross = 0.5;
double r_mut = 0.1;	
vector<vector<double>> cities;

// STRUCT
struct Population {
    int route[n_cities]={0};
    double fitness=0;
		int to_exchange=0;
};

Population pop;
vector<Population> current_pop;
vector<Population> new_pop(n_pop);




#endif
