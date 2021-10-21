#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include "random.h"

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


	double bar_length = 0.75;
	double line_distance = 1;
	int N_hit = 0;
	int N_throws = 1000;
	int repet = 10000;
	int N_block=100;
	int L = repet/N_block;


	double bar_center_x;
	double bar_extreme_x;
	double bar_direction;
	double pi[repet]={0};

	for(int j =0; j<repet; j+=1){

		N_hit=0;

		for(int i =0; i<N_throws; i+=1){

			bar_center_x = rnd.Rannyu();
			bar_direction = rnd.Rannyu(0,2*M_PI);

			bar_extreme_x = bar_center_x + cos(bar_direction)*bar_length;

			if (bar_extreme_x<0 || bar_extreme_x>1){
				N_hit+=1;
			}
		}
		
		pi[j]= 2*bar_length*N_throws/(N_hit*line_distance);

	}



	// Calcolo la media dei singoli blocchi 
	double mean[N_block]={0};
	double mean2[N_block]={0};

	int count = 0;

	for( int i =0; i<N_block; i+=1){
		double sum_mean = 0;
		for (int j=0; j<L; j+=1){
			sum_mean+=pi[j+count];
		}
		count+=L;
		mean[i]=sum_mean/L;
		mean2[i]=pow(mean[i],2);
	}

	// Calcolo ora la media per un numero crescente di blocchi, partendo da 1 fino a 100. 
	double sumprog_mean2[N_block];
	double sumprog_mean[N_block];
	double std_mean[N_block];

	
	for(int i =0; i<N_block; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=mean2[j];
			sumprog_mean[i]+=mean[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}


	// Trascrivo ora i risultati su file
	ofstream file_mean;
	file_mean.open ("file_mean.txt");

	for(int i=0;i<N_block;i+=1){
		file_mean << sumprog_mean[i]<<","<<std_mean[i]<<endl;
	}
	
	file_mean.close();


	rnd.SaveSeed();
	return 0;
}








