#include "City.h"
#include <cmath>



// COSTRUTTORE

City::City(double x_, double y_){
	
	x = x_; 
	y = y_;
	
}

// DISTRUTTORE

City::~City(){
}

// METODI

double City::Distance(double x_1, double y_1){ 
	return pow(pow(x-x_1,2)+pow(y-y_1,2),0.5);
}

double City::getX() const{
	return x;
}

double City::getY() const{
	return y;
}

