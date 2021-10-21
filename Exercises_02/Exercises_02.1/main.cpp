#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include "random.h"
#include <random>

using namespace std;
 
int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
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

	
	int M = 10000;
	int N = 100;
	int L = M/N;

	double int_unif[N]={0};
	double int_unif2[N]={0};
	double rand;

	
	//Calcolo ora l'integrale attraverso il metodo del campionamento uniforme.	
	for(int i=0;i<N; i+=1){
		for (int j=0; j<L;j+=1){
			rand = rnd.Rannyu();
			int_unif[i] += (M_PI/2)*cos(M_PI*rand/2);
		}
		int_unif[i]=int_unif[i]/L;
		int_unif2[i]=pow(int_unif[i],2);
		
	}

	// faccio ora le medie progressive per avere una stima dell'errore che commetto. 
	double sumprog_mean2[N]={0};
	double sumprog_mean[N]={0};
	double std_mean[N];

	
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=int_unif2[j];
			sumprog_mean[i]+=int_unif[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}


	// Trascrivo ora i risultati su file
	ofstream file_int_unif;
	file_int_unif.open ("file_int_unif.txt");

	for(int i=0;i<N;i+=1){
		file_int_unif << sumprog_mean[i]<<","<<std_mean[i]<<endl;
	}
	
	file_int_unif.close();



	//Calcolo ora l'integrale attraverso il metodo dell'Importance sampling.
	double int_importance[N]={0};
	double int_importance2[N]={0};
	double std_int_importance[N]={0};

	float x;

	ofstream verify;
	verify.open("verify.txt");

	for(int i=0;i<N; i+=1){
		for(int j=0;j<L;j++){
			rand = rnd.Rannyu();
			x =1 - sqrt(1-rand);
			verify<<x<<endl;
			int_importance[i] += (M_PI/2.)*cos(x*M_PI/2)/ ( 2 * ( 1 - x));
		}

		int_importance[i]=int_importance[i]/L;
		int_importance2[i]=pow(int_importance[i],2);

		sumprog_mean2[i]=0;
		sumprog_mean[i]=0;
			
	}
	
	verify.close();


	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=int_importance2[j];
			sumprog_mean[i]+=int_importance[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_int_importance[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}
		
	// Trascrivo ora i risultati su file
	ofstream file_int_importance;
	file_int_importance.open ("file_int_importance.txt");

	for(int i=0;i<N;i+=1){
		file_int_importance << sumprog_mean[i]<<","<<std_int_importance[i]<<endl;
	}
	
	file_int_importance.close();


	rnd.SaveSeed();
	return 0;
}












