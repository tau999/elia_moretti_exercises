#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "mpi.h"
#include "main.h"

using namespace std;

void Input(int);
Population InitialRoute(void);
void InizialPopulation(void);
double FitnessRoute(Population);
int Pbc(int);
void FitnessPopulation(void);
bool compareByFitness(const Population &a, const Population &b);
bool compareByToExchange(const Population &a, const Population &b);

vector<Population> RouteCrossover(int parent_A, int parent_B);
void NewPopulation(void);

Population InversionMutation (Population);
Population PermutationMutation (Population);
Population ShiftMutation (Population);

void PrintFitnessBestValue(int);
void PrintFitnessMean(int);
void PrintRoute(int, int );

void ExchangeBest(void);




int main (int argc, char* argv[]){
	
	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	Input(rank);

	InizialPopulation();
	FitnessPopulation();

	PrintFitnessBestValue(rank);
	PrintFitnessMean(rank);


	for (int i=0; i<=n_iter; i++){

		NewPopulation();
		FitnessPopulation();
		PrintFitnessMean(rank);
		PrintFitnessBestValue(rank);
		

		if ((i%100==0))
				PrintRoute(i,rank);

		if(i%50==0)
			ExchangeBest();
				
	}


	MPI_Finalize();
	return 0;
}



// ---------------------
// ------> Input <------
// ---------------------

void Input(int rank){

	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){	
			if (rank==0)
				Primes >> p1 >> p2;
			else if(rank==1)
				Primes >> p1 >> p2 >> p1 >> p2;
			else if(rank==3)
				Primes >> p1 >> p2 >> p1 >> p2 >> p1 >> p2;
			else
				Primes >> p1 >> p2 >> p1 >> p2 >> p1 >> p2 >> p1 >> p2;
			
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

	//LOAD CITIES FROM ./Exercises_09.2
	ifstream Cities("cities.dat");
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
	
	// REMOVE OLD FILE
	remove( "fitness_mean0.dat" );
	remove( "fitness_best0.dat" );
	remove( "fitness_mean1.dat" );
	remove( "fitness_best1.dat" );
	remove( "fitness_mean2.dat" );
	remove( "fitness_best2.dat" );
	remove( "fitness_mean3.dat" );
	remove( "fitness_best3.dat" );

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
		int equal_city = 0;
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
	
	for (int i = 0; i < n_pop; i++) 
		current_pop[i].fitness = FitnessRoute(current_pop[i]);

	// Create a sort population
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

	int StartCity = int(rnd.Rannyu() *(n_cities-1))+1;
  
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

	int City_to_swap_1 = int(rnd.Rannyu() *(n_cities-1))+1;
	int City_to_swap_2 = int(rnd.Rannyu() *(n_cities-1))+1;

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
// ------> compareByToExchange <------------------
// -------------------------------------------

bool compareByToExchange(const Population &a, const Population &b)
{
    return a.to_exchange > b.to_exchange;
}


// -------------------------------------------
// ------> Print Fitness Best Value <---------
// -------------------------------------------

void PrintFitnessBestValue(int rank){
	ofstream fitness;

	fitness.open("fitness_best"+to_string(rank)+".dat",ios::app);

	fitness << current_pop[0].fitness << endl;
  fitness.close();
}

// -------------------------------------------
// ------> Print Fitness Mean <---------------
// -------------------------------------------

void PrintFitnessMean(int rank){
	ofstream fitness_mean;
	fitness_mean.open("fitness_mean"+to_string(rank)+".dat",ios::app);

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

void PrintRoute(int gen, int rank){
	ofstream route;
	route.open("route_gen_"+to_string(gen)+"_"+to_string(rank)+".dat");

	for(int i = 0; i<n_cities; i++)
		route<<current_pop[0].route[i]<<endl;
	route<<0<<endl;
}

// ------------------------------------
// ------> ExchangeBest <--------------
// ------------------------------------

void ExchangeBest(void){

	for (int i =0 ; i<(int)n_pop/4; i++){
		//int best_route = (int) n_pop*pow(rnd.Rannyu(), 5);
		current_pop[i].to_exchange = 1;
	}

	// Create a sort population
	sort(current_pop.begin(), current_pop.end(), compareByToExchange);

	vector<Population> pop_to_cut = current_pop;
	pop_to_cut.erase (pop_to_cut.begin()+(int)n_pop/4, pop_to_cut.end());
	

	// DEFINE NEW MPI_TYPE
	const int nitems=3;
	int          blocklengths[3] = {n_cities,1,1};
	MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_INT};
	MPI_Datatype mpi_pop_type;
	MPI_Aint     offsets[3];

	offsets[0] = offsetof(Population, route);
	offsets[1] = offsetof(Population, fitness);
	offsets[2] = offsetof(Population, to_exchange);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_pop_type);
	MPI_Type_commit(&mpi_pop_type);
	
	MPI_Gather(pop_to_cut.data(), pop_to_cut.size(), mpi_pop_type, current_pop.data(), pop_to_cut.size(), mpi_pop_type, 0, MPI_COMM_WORLD);

	MPI_Bcast(current_pop.data(), current_pop.size(), mpi_pop_type, 0, MPI_COMM_WORLD);

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
	cout<< "Possible to_exchange: " << arr1.to_exchange << endl;
}
*/



