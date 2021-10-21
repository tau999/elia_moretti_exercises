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

	double C_direct[N]={0};
	double C_step[N]={0};
	double P_direct[N]={0};
	double P_step[N]={0};
	
	double C_direct2[N]={0};
	double C_step2[N]={0};
	double P_direct2[N]={0};
	double P_step2[N]={0};

	double rand;
	double price;

	double S_0 = 100;
	double K = 100;
	double T = 1;
	double r = 0.1;
	double sigma = 0.25;
	
	// Calcolo ora il prezzo delle azioni attraverso il ---> metodo DIRETTO <---  
	
	for(int i=0;i<N; i+=1){
		for (int j=0; j<L;j+=1){
			rand = rnd.Gauss(0,1);
			price = S_0*exp((r-sigma*sigma/2)*T+sigma*rand*pow(T,0.5));
			C_direct[i] += exp(-r*T)*max(0.,price-K);
			P_direct[i] += exp(-r*T)*max(0.,K-price);
		}
		C_direct[i]=C_direct[i]/L;
		C_direct2[i]=pow(C_direct[i],2);
		P_direct[i]=P_direct[i]/L;
		P_direct2[i]=pow(P_direct[i],2);
	}
	


	double C_sumprog_mean2[N]={0};
	double C_sumprog_mean[N]={0};
	double C_std_mean[N];
	double P_sumprog_mean2[N]={0};
	double P_sumprog_mean[N]={0};
	double P_std_mean[N];

	
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			C_sumprog_mean2[i]+=C_direct2[j];
			C_sumprog_mean[i]+=C_direct[j];
			P_sumprog_mean2[i]+= P_direct2[j];
			P_sumprog_mean[i]+=P_direct[j];
			
		}
		C_sumprog_mean2[i]=C_sumprog_mean2[i]/(i+1);
		C_sumprog_mean[i]=C_sumprog_mean[i]/(i+1);
		P_sumprog_mean2[i]=P_sumprog_mean2[i]/(i+1);
		P_sumprog_mean[i]=P_sumprog_mean[i]/(i+1);

		C_std_mean[i]=pow((C_sumprog_mean2[i]-pow(C_sumprog_mean[i],2))/(i+1),0.5);
		P_std_mean[i]=pow((P_sumprog_mean2[i]-pow(P_sumprog_mean[i],2))/(i+1),0.5);
	}

	// Trascrivo ora i risultati su file
	ofstream file_direct;
	file_direct.open ("file_direct.txt");

	for(int i=0;i<N;i+=1){

		file_direct << C_sumprog_mean[i]<<","<<C_std_mean[i]<<","<<P_sumprog_mean[i]<<","<<P_std_mean[i]<<endl;
	}
	
	file_direct.close();


	// Calcolo ora il prezzo delle azioni attraverso il ---> metodo DISCRETO <---  
	
	double N_time = 100;
	
	for(int i=0;i<N; i+=1){
		for (int j=0; j<L;j+=1){
				price = S_0;
				for(int k=0;k<N_time;k++){
					rand = rnd.Gauss(0,1);
					price = price*exp((r-sigma*sigma/2)*T/(double) N_time+sigma*rand*pow(T/(double) N_time,0.5));
				}
			C_step[i] += exp(-r*T)*max(0.,price-K);
			P_step[i] += exp(-r*T)*max(0.,K-price);
		}
		C_step[i]=C_step[i]/L;
		C_step2[i]=pow(C_step[i],2);
		P_step[i]=P_step[i]/L;
		P_step2[i]=pow(P_step[i],2);
	}

	for (int i =0; i<N;i++){
		C_sumprog_mean2[i] = 0;
		C_sumprog_mean[i] = 0;
		C_std_mean[i] = 0;
		P_sumprog_mean2[i] = 0;
		P_sumprog_mean[i] = 0;
		P_std_mean[i] = 0;
	}
	
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			C_sumprog_mean2[i]+=C_step2[j];
			C_sumprog_mean[i]+=C_step[j];
			P_sumprog_mean2[i]+= P_step2[j];
			P_sumprog_mean[i]+=P_step[j];
			
		}
		C_sumprog_mean2[i]=C_sumprog_mean2[i]/(i+1);
		C_sumprog_mean[i]=C_sumprog_mean[i]/(i+1);
		P_sumprog_mean2[i]=P_sumprog_mean2[i]/(i+1);
		P_sumprog_mean[i]=P_sumprog_mean[i]/(i+1);

		C_std_mean[i]=pow((C_sumprog_mean2[i]-pow(C_sumprog_mean[i],2))/(i+1),0.5);
		P_std_mean[i]=pow((P_sumprog_mean2[i]-pow(P_sumprog_mean[i],2))/(i+1),0.5);
	}

	// Trascrivo ora i risultati su file
	ofstream file_step;
	file_step.open ("file_step.txt");

	for(int i=0;i<N;i+=1){

		file_step << C_sumprog_mean[i]<<","<<C_std_mean[i]<<","<<P_sumprog_mean[i]<<","<<P_std_mean[i]<<endl;
	}
	
	file_step.close();

	
	rnd.SaveSeed();
	return 0;
}



