#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
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


	// Considero 100000 estrazioni da raggruppare in 100 blocchi da 100 estrazioni
	int M = 100000;
	int N = 100;
	int L = M/N;
	double rand;

	// Definisco le variabili per la media
	double mean[N];
	double mean2[N];
	double std_mean[N];

	// Definisco le variabili per la varianza
	double sigma[N];
	double sigma2[N];
	double std_sigma[N];

	// Apro i canali per l'output
	ofstream file_mean;
	ofstream file_sigma;
	ofstream file_chi_quadro;
	

	// Calcolo la media e la varianza dei singoli blocchi 
	for( int i =0; i<N; i+=1){
		double sum_mean = 0;
		double sum_sigma =0;
		for (int j=0; j<L; j+=1){
			rand = rnd.Rannyu();
			sum_mean+=rand;
			sum_sigma+= pow(rand-0.5,2);
		}
		mean[i]=sum_mean/L;
		mean2[i]=pow(mean[i],2);
		sigma[i]=sum_sigma/L;
		sigma2[i]=pow(sigma[i],2);
	}


	// Calcolo ora la media e la varianza per un numero crescente di blocchi, partendo da 1 fino a 100. 
	double sumprog_mean2[N];
	double sumprog_mean[N];
	double sumprog_sigma2[N];
	double sumprog_sigma[N];
	
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=mean2[j];
			sumprog_mean[i]+=mean[j];
			sumprog_sigma2[i]+=sigma2[j];
			sumprog_sigma[i]+=sigma[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		sumprog_sigma2[i]=sumprog_sigma2[i]/(i+1);
		sumprog_sigma[i]=sumprog_sigma[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
		std_sigma[i]=pow((sumprog_sigma2[i]-pow(sumprog_sigma[i],2))/(i+1),0.5);
	}



	// Trascrivo ora i risultati su file
	file_mean.open ("file_mean.txt");
	file_sigma.open ("file_sigma.txt");

	for(int i=0;i<N;i+=1){
		file_mean << sumprog_mean[i]<<","<<std_mean[i]<<endl;
		file_sigma << sumprog_sigma[i]<<","<<std_sigma[i]<<endl;	
	}
	
	file_mean.close();
	file_sigma.close();

	
	// Verifichiamo ora l'affidabilità del nostro generatore attraverso un test del chi-quadro. Per farlo dividiamo il nostro intervallo in M=100 sottointervalli e, generando n=10000 numeri casuali, vediamo come questi si distribuiscono nei 100 bin. In particolare, ci aspettiamo di ottenere un chi-quadro di 100. Ripetiamo poi questa procedura N=100 volte.

	M = 100;
	int n = 10000;
	N =100;
	double chi_quadro[N]={0};

	

	for(int i = 0; i<N;i+=1){
	
		int occupation_number[M]={0}; 
	
		for(int j =0; j<n; j+=1){
			rand = rnd.Rannyu();
			for(int k =0; k<M ;k+=1){
				if(k/(double)M<rand && rand<(k+1)/(double)M){
					occupation_number[k]+=1;		// Mi dice il numero di estreazioni finite nel M-esimo bin	

				}
			}
		}
	
		// Una volta ottenuto il numero di occupazione per ogni bin, possiamo calcolare il valore di chi-quadro
		for(int p=0; p<M;p+=1){
			chi_quadro[i]+=pow(occupation_number[p]-n/M,2)/(n/M);
		}
	}

	file_chi_quadro.open ("file_chi_quadro.txt");

	for(int i=0;i<N;i+=1){
		file_chi_quadro << chi_quadro[i]<<endl;
	}
	
	file_chi_quadro.close();
	

	
	rnd.SaveSeed();
	return 0;
}













