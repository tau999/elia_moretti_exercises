#ifndef __Fitness_h__
#define __Fitness_h__

using namespace std;


class Fitness {

public:

	// Costruttori
	Fitness(vector<int>* route_);
	~Fitness(); 




private: 
	vector<int>* route;
	double distance, fitness_value;

};	

	
#endif
	
