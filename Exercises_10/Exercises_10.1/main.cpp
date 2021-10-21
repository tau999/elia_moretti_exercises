#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#include "main.h"

using namespace std;

void Input(int,char**);
void InitialRoute(void);
double FitnessRoute(Population);
int Pbc(int);
void NewRoute(void);

Population InversionMutation (Population);
Population PermutationMutation (Population);
Population ShiftMutation (Population);

void PrintRoute(int gen);
void PrintFitness(void);



int main (int argc, char ** argv){
	
	Input(argc, argv);

	InitialRoute();
	current_route.fitness = FitnessRoute(current_route);

	for (int i=0; i<n_iter; i++){

		T =  1 - (double) (i+1)/n_iter;		
		NewRoute();
		

		if (i%100 == 0 )
			PrintFitness();
		if (i==n_iter-1 || i == 0)
			PrintRoute(i+1);

	}


return 0;
}


// ---------------------
// ------> Input <------
// ---------------------

void Input(int argc, char ** argv){

	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){	
			Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			 input >> property;
			 if( property == "RANDOMSEED" ){
			    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			    rnd.SetRandom(seed,p1,p2);
			 }
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	// SET TYPE OF CITIES (CIRCLE OR SQUARE)

	if(argc!=2) 
		cout<<" 	ERROR: 0 -> CIRCLE    &      1 -> SQUARE"<<endl;
	else 
		type_cities = atoi(argv[1]);

	//LOAD CITIES FROM ./Exercises_09.2
	ifstream Cities;
	if(type_cities == 0 )
		Cities.open("/home/elia/Scrivania/SimulazioneNumerica/Exercises_09/Exercises_09.1/cities.dat");
	else 
		Cities.open("/home/elia/Scrivania/SimulazioneNumerica/Exercises_09/Exercises_09.2/cities.dat");

	double x;
	double y;
	for (int i = 0; i < n_cities; i++) { 
      vector<double> a1; 
     	
			Cities >> x >> y;
			a1.push_back(x); 
			a1.push_back(y); 

      cities.push_back(a1); 
    } 
	Cities.close();
	
	// SAVE CITIES TO FILE

	ofstream Cities_s;
	Cities_s.open("cities.dat");
	
	for (unsigned int i = 0; i < cities.size(); i++) { 
        for (unsigned int j = 0; j < cities[i].size(); j++) 
             Cities_s << cities[i][j] << " "; 
        Cities_s << endl; 
    } 

	Cities_s.close();

	// REMOVE OLD FILE
	remove( "fitness.dat.dat" );

}

// ----------------------------
// ------> InitialRoute <------
// ----------------------------

void InitialRoute(void){
	
	
	int i = 1; 
	while(i < n_cities){
		int r = rnd.Rannyu(1,n_cities);
		int diff=0;
		for (int j=0; j<i; j++){
			if(r == current_route.route[j])	diff++;		 
		}

		if(diff == 0 ){
			current_route.route[i]=r; 
			i++;
		}
	}
}


// ---------------------------
// ------> FitnessRoute <-----
// ---------------------------

double FitnessRoute(Population pop_individual){
	
	double sum = 0;

	for (int i = 0; i<n_cities; i ++ )
  	sum += pow(pow(cities[pop_individual.route[Pbc(i)]][0]-cities[pop_individual.route[Pbc(i+1)]][0],2)+pow(cities[pop_individual.route[Pbc(i)]][1]-cities[pop_individual.route[Pbc(i+1)]][1],2),0.5);
       
	return sum;
}



// ------------------------------
// ------> NewRoute <-------
// ------------------------------

void NewRoute(void){
	double select_mutation = rnd.Rannyu();
	if (select_mutation<0.33)
		new_route = InversionMutation(current_route);

	else if(0.33<select_mutation && select_mutation<0.66)
		new_route = PermutationMutation(current_route);

	else
		new_route = ShiftMutation(current_route);

	new_route.fitness = FitnessRoute(new_route);
	attempted++;

	
	if( (new_route.fitness - current_route.fitness) < 0 ){
	  current_route=new_route;
	  accepted++;
	    
	}
	else{
	  double r = rnd.Rannyu();
	  double A = exp(-(1./T) * (new_route.fitness - current_route.fitness));

	  
	  if(r < min(1.,A))
	  {
	    current_route=new_route;
	  	accepted++;
	  }
	}


}




// --------------------------------
// ------> InversionMutation <-----
// --------------------------------

Population InversionMutation (Population child){

	int City_to_swap_1 = int(rnd.Rannyu(1,n_cities));
	int City_to_swap_2 = int(rnd.Rannyu(1,n_cities));

  int city1 = child.route[City_to_swap_1];
	int city2 = child.route[City_to_swap_2];
            
	child.route[City_to_swap_1] = city2;
	child.route[City_to_swap_2] = city1;

	return child;
}


// ----------------------------------
// ------> PermutationMutation <-----
// ----------------------------------

Population PermutationMutation (Population child){

	int City_to_swap_1 = int(rnd.Rannyu(1,n_cities-1));

  int city1 = child.route[City_to_swap_1];
	int city2 = child.route[City_to_swap_1+1];
            
	child.route[City_to_swap_1] = city2;
	child.route[City_to_swap_1+1] = city1;

	return child;
}


// -----------------------------
// ------> ShiftMutation <------
// -----------------------------

Population ShiftMutation (Population child){

	int City_to_swap =  int(rnd.Rannyu(1,n_cities));
	vector<int> appoggio;

	for (int i = 0; i < n_cities; i++)
		appoggio.push_back(child.route[i]);

	rotate(appoggio.begin()+1,appoggio.begin()+City_to_swap,appoggio.end());

	for (int i = 0; i < n_cities; i++)
		child.route[i]=appoggio[i];

	return child;
}



// -------------------------------------------
// ------> Periodic Boundary Conditions <-----
// -------------------------------------------

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= n_cities) i = i - n_cities;
    else if(i < 0) i = i + n_cities;
    return i;
}


// ------------------------------------
// ------> Print Route <---------------
// ------------------------------------

void PrintRoute(int gen){
	ofstream route;
	route.open("route_gen_"+to_string(gen)+".dat");

	for(int i = 0; i<n_cities; i++)
		route<<current_route.route[i]<<endl;
	route<<0<<endl;
}


// ------------------------------------
// ------> Print Fitness <---------------
// ------------------------------------

void PrintFitness(void){
	ofstream fitness;

	fitness.open("fitness.dat",ios::app);

	fitness << current_route.fitness << endl;
  fitness.close();
}

// -------------------------------------------
// ------> Print struct <------------------
// -------------------------------------------
/* 
cout<<endl<<endl;
for (const auto &arr1 : current_pop) {
	cout << "Route: ";
		for (int i =0; i<n_cities; i++){
		cout<< arr1.route[i] << "  ";
		}
	cout<<endl;
	cout<< "FitnessRoute: " << arr1.fitness << endl;
	cout<< "Possible parents: " << arr1.possible_parents << endl;
}
*/



