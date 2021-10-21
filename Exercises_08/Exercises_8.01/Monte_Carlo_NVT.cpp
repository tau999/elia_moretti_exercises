#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include "random.h"
#include <random>
#include <iomanip>

using namespace std;

double psi(double x, double mu, double sigma);
double integral (double x, double mu, double sigma);
double Error(double sum, double sum2, int iblk);

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

// READ PROBLEM PARAMETER
	double  mu, sigma, delta_x, x_n;
	int nstep, nblock;
	
	ifstream ReadInput;
	ReadInput.open("input.dat");
	
	ReadInput >> nstep;
	ReadInput >> nblock;

	ReadInput >> mu;  // 0.76;
	ReadInput >> sigma; // 0.6;
	ReadInput >> delta_x; // = 2;
	ReadInput >> x_n;

	ReadInput.close();

// HISTOGRAM PARAMETERS

	int nbins = 100;
	int bin_integral[nbins] = {0};
	double bin_min = -3;
	double bin_max = 3;
	double bin_size = (bin_max-bin_min)/nbins;
	double glob_av[nbins] = {0};
	double glob_av2[nbins] = {0};
	int tot = 0;
	
	

// INTEGRAL CALCULUS

	double x;
	double integral_n;
	int attempted = 0;
	int accepted = 0;

	integral_n = integral(x_n, mu, sigma);


	double mean_integral[nblock]={0};
	double mean_integral2[nblock]={0};


	for (int i = 0; i < nblock; i+=1){

		for(int k=0; k<nbins; ++k)
		{
			bin_integral[k]=0.0;
			tot = 0;
		}
		attempted = 0;
		accepted = 0;


		for (int j = 0; j<nstep; j+=1){

			x = rnd.Rannyu(x_n+delta_x,x_n-delta_x);
			double A = min(1., psi(x,mu, sigma)*psi(x,mu, sigma)/(psi(x_n,mu, sigma)*psi(x_n,mu, sigma)));

			double rand = rnd.Rannyu();

			if( rand <= A){ 

				integral_n = integral(x, mu, sigma);
				x_n = x;	
				
				accepted ++;
			}

			attempted ++;

			for (int k =0; k<nbins; k+=1){	
				double inf = bin_min+bin_size*k;
				double max = bin_min+bin_size*(k+1); 																				
				if( (inf<x_n) && (x_n<max)){				
					 bin_integral[k]+=1;																																
				}																																													
			}

			mean_integral[i] += integral_n;

		}		
 
		for(int k=0; k<nbins; ++k){ 
			tot += bin_integral[k];
		}

		ofstream histogram;
		histogram.open("output.histo.0");

    for(int k=0; k<nbins; k++)
    {

        double x_position = k * bin_size + bin_min+bin_size*0.5;
        double estimate_histo = bin_integral[k]/(bin_size * tot);
        glob_av[k] += estimate_histo;
        glob_av2[k] += estimate_histo*estimate_histo;

        double err_histo = Error(glob_av[k], glob_av2[k], i+1);

				if(i== nblock-1) {  
					histogram << " " << x_position << " 	" << glob_av[k]/(double)(i+1.0) << " 	" << err_histo << endl;
				}
    }
    
	
	}

	cout<<"Rate di accettazione "<< accepted / (double) attempted<<endl;
	

	for(int j = 0; j<nblock; j++){
		mean_integral[j] = mean_integral[j]/nstep;
		mean_integral2[j] = pow(mean_integral[j],2);

	}

	// faccio ora le medie progressive per avere una stima dell'errore che commetto. 
	double sumprog_mean2[nblock]={0};
	double sumprog_mean[nblock]={0};
	double std_mean[nblock];

	
	for(int i =0; i<nblock; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=mean_integral2[j];
			sumprog_mean[i]+=mean_integral[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}



	// Trascrivo ora i risultati su file
	ofstream file_mean_integral;
	file_mean_integral.open ("file_mean_integral.txt");	

	for(int i=0;i<nblock;i+=1){
		file_mean_integral << i*nstep <<","<< sumprog_mean[i] <<","<< std_mean[i] <<endl;
	}


	
	rnd.SaveSeed();
	return 0;
}




double psi(double x, double mu, double sigma) {

	double psi = exp( - (x - mu)*(x - mu) / (2. * sigma*sigma) ) + exp( - (x + mu)*(x + mu) / (2. * sigma*sigma) );
	return psi;
}


double integral (double x, double mu, double sigma) {
	
	double kinetic = exp( - (x - mu)*(x - mu) / (2. * sigma*sigma) ) * ( (x-mu)*(x-mu) / (sigma*sigma*sigma*sigma) - 1./(sigma*sigma)) + exp( - (x + mu)*(x + mu) / (2. * sigma*sigma) ) * ( (x+mu)*(x+mu) / (sigma*sigma*sigma*sigma) - 1./(sigma*sigma));

	double potential = pow(x,4.)-2.5*pow(x,2.);

	double func = (-0.5 * kinetic  + potential* psi(x, mu, sigma)) / (psi(x, mu, sigma)) ;

	return func; 
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


