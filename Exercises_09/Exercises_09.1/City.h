#ifndef __City_h__
#define __City_h__


class City {

public:

	// Costruttori
	City(double x_, double y_);
	~City(); 

	// METODI
	double Distance(double x_1, double y_1);
	double getX() const;       
  double getY() const;


private: 
	double x, y; 

};	

	
#endif
	
