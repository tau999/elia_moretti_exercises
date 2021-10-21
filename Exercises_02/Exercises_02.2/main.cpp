#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include "random.h"
#include <random>

using namespace std;

double Error(double , double , int );
 
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


	//Discrete random walk

	int N_rep = 10000;
	int N_step = 100;

	
	int position[N_rep][3] = {0}; 	//Le tre componenti del vettore descrivono le tre coordinate
	int direction = 0;				// Lo scelgo int poiché può essere x,y o z ossia la componente 0,1 o 2 del vettore position
	double versus = 0; 				
	double MSD[N_step] = {0};		//MSD=mean square displacement
	double std_MSD[N_step] = {0};


	for(int i=1;i<=N_step; i+=1){
		for (int j=0; j<N_rep;j+=1){
			direction = rnd.Rannyu()*3; 			
			versus = rnd.Rannyu()*2-1;

			if(versus<0){	versus = -1;}
			else {	versus = 1;}

			position[j][direction]+=versus;

			MSD[i] += pow(position[j][0],2)+ pow(position[j][1],2) + pow(position[j][2],2);
			std_MSD[i] += pow(position[j][0],2)+ pow(position[j][1],2) + pow(position[j][2],4);
			
		}

		std_MSD[i] = Error(pow(MSD[i],0.5)/(double) N_rep ,pow(std_MSD[i],0.5)/(double) N_rep,i+1);

		MSD[i] = pow(MSD[i]/N_rep,0.5);

		
	}


	// Trascrivo ora i risultati su file
	ofstream file_MSD_discreto;
	file_MSD_discreto.open ("file_MSD_discreto.txt");

	for(int i=0;i<N_step;i+=1){
		file_MSD_discreto << MSD[i]<<","<<std_MSD[i]<<endl;
	}
	
	file_MSD_discreto.close();




	// Continuos random walk

	double direction_theta = 0;
	double direction_phi = 0;

	double position_continuous[N_rep][3] = {0};

	for(int i=1;i<=N_step; i+=1){
		for (int j=0; j<N_rep;j+=1){
			direction_theta = acos(1-2*rnd.Rannyu()); 
			direction_phi = rnd.Rannyu()*2*M_PI;

			position_continuous[j][0]+=sin(direction_theta)*cos(direction_phi);
			position_continuous[j][1]+=sin(direction_theta)*sin(direction_phi);
			position_continuous[j][2]+=cos(direction_theta);

			MSD[i] += pow(position_continuous[j][0],2)+ pow(position_continuous[j][1],2) + pow(position_continuous[j][2],2);
			std_MSD[i] += pow(position_continuous[j][0],2)+ pow(position_continuous[j][1],2) + pow(position_continuous[j][2],4);
		}

		std_MSD[i] = Error(pow(MSD[i],0.5)/(double) N_rep ,pow(std_MSD[i],0.5)/(double) N_rep,i+1);
		MSD[i] = pow(MSD[i]/N_rep,0.5);
		
	}

	// Trascrivo ora i risultati su file
	ofstream file_MSD_continuous;
	file_MSD_continuous.open ("file_MSD_continuous.txt");

	for(int i=0;i<N_step;i+=1){
		file_MSD_continuous << MSD[i]<<","<<std_MSD[i]<<endl;
	}
	
	file_MSD_continuous.close();



	rnd.SaveSeed();
	return 0;
}

double Error(double sum, double sum_sq, int n) {
  if (n == 1)
    return 0;
  else
    return sqrt((sum_sq / double(n) - pow(sum / double(n), 2)) / double(n - 1));
}


