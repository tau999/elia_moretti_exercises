#include "Fitness.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;


// COSTRUTTORE

Fitness::Fitness(vector<int>* route_){
	
	route = route_; 
	distance = 0.;
	fitness_value = 0.;
	
}

// DISTRUTTORE

Fitness::~Fitness(){
}


