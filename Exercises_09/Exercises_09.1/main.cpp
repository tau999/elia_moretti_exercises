#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "main.h"

using namespace std;

void Input(void);
Population InitialRoute(void);
void InizialPopulation(void);
double FitnessRoute(Population);
int Pbc(int);
void FitnessPopulation(void);
bool compareByFitness(const Population &a, const Population &b);
vector<Population> RouteCrossover(int parent_A, int parent_B);
Population InversionMutation (Population);
void NewPopulation(void);
void PrintFitnessBestValue(void);
void PrintFitnessMean(void);
void PrintRoute(int gen);
Population PermutationMutation (Population);
Population ShiftMutation (Population);


int main (){
	
	Input();

	InizialPopulation();
	FitnessPopulation();
	PrintFitnessBestValue();
	PrintFitnessMean();

	for (int i=0; i<=n_iter; i++){
		
		NewPopulation();
		FitnessPopulation();
		PrintFitnessMean();
		PrintFitnessBestValue();
		if (i%100 == 0 )
			PrintRoute(i);
	}

return 0;
}


// ---------------------
// ------> Input <------
// ---------------------

void Input(void){

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

	//CREATE CITIES on a CIRCLE
	for (int i = 0; i < n_cities; i++) { 
      vector<double> a1; 
     	
			double r = rnd.Rannyu()*M_PI*2.;
			double x = cos(r);
			double y = sin(r);
			a1.push_back(x); 
			a1.push_back(y); 

      cities.push_back(a1); 
    } 

	// SAVE CITIES TO FILE

	ofstream Cities;
	Cities.open("cities.dat");
	
	for (unsigned int i = 0; i < cities.size(); i++) { 
        for (unsigned int j = 0; j < cities[i].size(); j++) 
             Cities << cities[i][j] << " "; 
        Cities << endl; 
    } 

	Cities.close();	

	// REMOVE OLD FILE
	remove( "fitness_mean.dat" );
	remove( "fitness_best.dat" );

}

// ----------------------------
// ------> InitialRoute <------
// ----------------------------

Population InitialRoute(void){
	
	Population one_route;
	
	int i = 1; 
	while(i < n_cities){
		int r = rnd.Rannyu(1,n_cities);
		int diff=0;
		for (int j=0; j<i; j++){
			if(r == one_route.route[j])	diff++;		 
		}

		if(diff == 0 ){
			one_route.route[i]=r; 
			i++;
		}
	}
	return	one_route;
}

// --------------------------------
// ------> InizialPopulation <-----
// --------------------------------

void InizialPopulation(void){
		
	Population one_route;

	int i =0;
	while (i < n_pop){
		one_route = InitialRoute();
		int equal_route=0;
		for (int j = 0; j < i; j++){
		int equal_city= 0;
			for(int k = 0; k<n_cities;k++){
				if (current_pop[j].route[k] == one_route.route[k]){
					equal_city ++;
				}
			}
			if(equal_city == n_cities){
				equal_route ++;
				break;
			}
		}
		if (equal_route == 0){
				current_pop.push_back(one_route); 
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


// --------------------------------
// ------> FitnessPopulation <-----
// --------------------------------

void FitnessPopulation(void){
	
	// Create a sort population
	for (int i = 0; i < n_pop; i++) 
		current_pop[i].fitness = FitnessRoute(current_pop[i]);

	sort(current_pop.begin(), current_pop.end(), compareByFitness);

}

// ------------------------------
// ------> NewPopulation <-------
// ------------------------------

void NewPopulation(void){

	int n=0;
	vector<Population> child;

	while(n<n_pop-1){
		int parent_A = (int) n_pop*pow(rnd.Rannyu(), 5);
		int parent_B = (int) n_pop*pow(rnd.Rannyu(), 5);

		child.push_back(current_pop[parent_A]);
		child.push_back(current_pop[parent_B]);

		if(rnd.Rannyu()<r_cross)
			child = RouteCrossover(parent_A,parent_B);
		
		if(rnd.Rannyu()<r_mut)
			child[0] = InversionMutation(child[0]);
		if(rnd.Rannyu()<r_mut)
			child[1] = InversionMutation(child[1]);
		
		if(rnd.Rannyu()<r_mut)
			child[0] = PermutationMutation(child[0]);
		if(rnd.Rannyu()<r_mut)
			child[1] = PermutationMutation(child[1]);

		if(rnd.Rannyu()<r_mut)
			child[0] = ShiftMutation(child[0]);
		if(rnd.Rannyu()<r_mut)
			child[1] = ShiftMutation(child[1]);


		new_pop[n] = child[0];
		new_pop[n+1] = child[1];

		n+=2;
	}

	current_pop = new_pop;

}




// ------------------------------
// ------> Route Crossover <-----
// ------------------------------

vector<Population> RouteCrossover(int parent_A, int parent_B){
		
	Population child_A;
	Population child_B;

	int StartCity = int(rnd.Rannyu(1,n_cities));
  
  for (int i=0; i < StartCity; i++)
      child_A.route[i] = current_pop[parent_A].route[i];

	for (int i=0; i < StartCity; i++)
      child_B.route[i] = current_pop[parent_B].route[i];


	int k = StartCity;
	for(int i =0; i < n_cities; i ++){
		int diff=0;
		for (int j = 0; j < n_cities; j++){
				if (current_pop[parent_A].route[i] == child_B.route[j] ){
					diff+=1;
				}
		}
		if (diff==0){
			child_B.route[k]=current_pop[parent_A].route[i];	
			k++;
		}
	}

	k = StartCity;
	for(int i =0; i < n_cities; i ++){
		
		int diff=0;
		for (int j = 0; j < n_cities; j++){
				if (current_pop[parent_B].route[i] == child_A.route[j] ){
					diff+=1;
				}
		}
		if (diff==0){
			child_A.route[k]=current_pop[parent_B].route[i];	
			k++;
		}
	}

	child_A.fitness = 0;
	child_B.fitness = 0;

	vector<Population> child_A_B;
	child_A_B.push_back(child_A);
	child_A_B.push_back(child_B);

	return child_A_B;
}


// ------------------------------
// ------> InversionMutation <-----
// ------------------------------

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

// -------------------------------------------
// ------> compareByFitness <------------------
// -------------------------------------------

bool compareByFitness(const Population &a, const Population &b)
{
    return a.fitness < b.fitness;
}

// -------------------------------------------
// ------> Print Fitness Best Value <---------
// -------------------------------------------

void PrintFitnessBestValue(void){
	ofstream fitness;

	fitness.open("fitness_best.dat",ios::app);

	fitness << current_pop[0].fitness << endl;
  fitness.close();
}

// -------------------------------------------
// ------> Print Fitness Mean <---------------
// -------------------------------------------

void PrintFitnessMean(void){
	ofstream fitness_mean;
	fitness_mean.open("fitness_mean.dat",ios::app);

	int	half_best = n_pop/2;
	double sum = 0;
	for(int i = 0; i<half_best; i++)
		sum += current_pop[i].fitness;
	sum = sum/half_best;

	fitness_mean << sum << endl;
  fitness_mean.close();
}

// ------------------------------------
// ------> Print Route <---------------
// ------------------------------------

void PrintRoute(int gen){
	ofstream route;
	route.open("route_gen_"+to_string(gen)+".dat");

	for(int i = 0; i<n_cities; i++)
		route<<current_pop[0].route[i]<<endl;
	route<<0<<endl;
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



