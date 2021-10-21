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


	// Usiamo un dato con distribuzione uniforme. Quello che facciamo è calcolare 10000 volte la media di N={1,2,10,100} lanci. Scriviamo poi i risultati su un file. 

	double dice;
	int N[4]={1,2,10,100};
	int n_throws=10000;
	double sumprog=0;
	ofstream file_uniform;

	file_uniform.open("file_uniform.txt");


	for(int i=0; i<n_throws; i+=1){			

		for(int j=0;j<4;j+=1){

			sumprog=0;
			
			for(int k=0; k<N[j];k+=1){
				
				dice = rnd.Rannyu();
				sumprog+=dice;
			}
			sumprog=sumprog/N[j];
			file_uniform << sumprog<< " ";


		}
		
		file_uniform<<endl;
		
	}

	file_uniform.close();

	// Usiamo ora lo stesso procedimento per un dado con distribuzione esponenziale. Per estrarre i risultati ho aggiunto un membro alla classe Random che lei ci ha fornito.

	ofstream file_exp;

	file_exp.open("file_exp.txt");


	for(int i=0; i<n_throws; i+=1){			

		for(int j=0;j<4;j+=1){

			sumprog=0;
			
			for(int k=0; k<N[j];k+=1){
				
				dice = rnd.Exp(1);
				sumprog+=dice;
			}
			sumprog=sumprog/N[j];
			file_exp << sumprog<< " ";
		}
		
		file_exp<<endl;
		
	}

	file_exp.close();


	// Usiamo ora lo stesso procedimento per un dado con distribuzione Lorentziana. Per estrarre i risultati ho aggiunto un membro alla classe Random che lei ci ha fornito.

	ofstream file_lorentz;

	file_lorentz.open("file_lorentz.txt");


	for(int i=0; i<n_throws; i+=1){			

		for(int j=0;j<4;j+=1){

			sumprog=0;
			
			for(int k=0; k<N[j];k+=1){
				
				dice = rnd.Lorentzian(0,1);
				sumprog+=dice;
			}
			sumprog=sumprog/N[j];
			file_lorentz << sumprog<< " ";
		}
		
		file_lorentz<<endl;
		
	}

	file_lorentz.close();
	

	return 0;
}






