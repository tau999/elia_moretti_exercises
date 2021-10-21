#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include "random.h"
#include <random>

using namespace std;

double pdf_1(double r, double r_n);
double pdf_2(double r, double r_n, double z, double z_n);
 
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

	
	double x, y, z;
	double x_n, y_n, z_n;

	

	//Calcolo ora il valor medio di r come integrale attraverso il sampling della pdf.

	

	

	// Come prima cosa valuto i "cammini" del raggio per entrambi gli orbitali. Eseguo 2000 passi. 

	// Considero prima il caso in cui parto dall'origine

	double r;
	double r_n;
	int N = 10000;
	double delta_x=1.2;

	x_n = y_n = z_n = 0;
	r_n = pow(x_n*x_n+y_n*y_n+z_n*z_n,0.5);

	ofstream file_path_1;
	file_path_1.open("file_path_1.txt");


	for (int i = 0; i < N; i+=1){
			x = rnd.Rannyu(x_n+delta_x,x_n-delta_x);
			y = rnd.Rannyu(y_n+delta_x,y_n-delta_x);
			z = rnd.Rannyu(z_n+delta_x,z_n-delta_x);
			
			r = pow(x*x+y*y+z*z,0.5);

			double A = pdf_1(r,r_n);
			double rand = rnd.Rannyu();

			if( rand <= A){ 
			
				r_n = r;
				x_n = x;	
				y_n = y;
				z_n = z;	

			}

			file_path_1 << r_n <<endl;
	}

	file_path_1.close();

	x_n = y_n = z_n = 0;

	ofstream file_path_2;
	file_path_2.open("file_path_2.txt");

	for (int i = 0; i < N; i+=1){
			x = rnd.Rannyu(x_n+delta_x,x_n-delta_x);
			y = rnd.Rannyu(y_n+delta_x,y_n-delta_x);
			z = rnd.Rannyu(z_n+delta_x,z_n-delta_x);
			
			r = pow(x*x+y*y+z*z,0.5);

			double A = pdf_2(r,r_n,z,z_n);
			double rand = rnd.Rannyu();

			if( rand <= A){ 
			
				r_n = r;
				x_n = x;	
				y_n = y;
				z_n = z;	
				
			}

			file_path_2 << r_n <<endl;

	}
	
	file_path_2.close();


	// Considero ora il caso in cui parto lontano dall'origine

	x_n = y_n = z_n = 30;
	r_n = pow(x_n*x_n+y_n*y_n+z_n*z_n,0.5);
	

	file_path_1.open("file_path_1_no_origin.txt");

	for (int i = 0; i < N; i+=1){
			x = rnd.Rannyu(x_n+delta_x,x_n-delta_x);
			y = rnd.Rannyu(y_n+delta_x,y_n-delta_x);
			z = rnd.Rannyu(z_n+delta_x,z_n-delta_x);
			
			r = pow(x*x+y*y+z*z,0.5);

			double A = pdf_1(r,r_n);
			double rand = rnd.Rannyu();

			if( rand <= A){ 
			
				r_n = r;
				x_n = x;	
				y_n = y;
				z_n = z;	
				
			}

			file_path_1 << r_n <<endl;

	}
	
	file_path_1.close();

	x_n = y_n = z_n = 30;
	r_n = pow(x_n*x_n+y_n*y_n+z_n*z_n,0.5);

	file_path_2.open("file_path_2_no_origin.txt");

	for (int i = 0; i < N; i+=1){
			x = rnd.Rannyu(x_n+delta_x,x_n-delta_x);
			y = rnd.Rannyu(y_n+delta_x,y_n-delta_x);
			z = rnd.Rannyu(z_n+delta_x,z_n-delta_x);
			
			r = pow(x*x+y*y+z*z,0.5);

			double A = pdf_2(r,r_n,z,z_n);
			double rand = rnd.Rannyu();

			if( rand <= A){ 
			
				r_n = r;
				x_n = x;	
				y_n = y;
				z_n = z;	
				
			}

			file_path_2 << r_n <<endl;

	}
	
	file_path_2.close();


	// Calcolo ora la lunghezza ideale dei blocchi attraverso il TLC. Ossia vario la lunghezza del blocco in [1,10,100,1000] unità e valuto quando la distribuzione è guassiana. Considero solo il caso in cui si parta da r=0.

	int n_rep = 100000;
	int N_block[4] = {10,100,200,500};
	

	ofstream file_N_block;
	

	for(int k =0; k<4;k+=1){

		file_N_block.open("file_N_block_"+to_string(N_block[k])+".txt");
		r_n = 0;
		x_n = y_n = z_n = 0;

		for (int i = 0; i < n_rep; i+=1){

			double sum_r = 0;

			for (int j = 0; j<N_block[k]; j+=1){

				x = rnd.Rannyu(x_n+delta_x,x_n-delta_x);
				y = rnd.Rannyu(y_n+delta_x,y_n-delta_x);
				z = rnd.Rannyu(z_n+delta_x,z_n-delta_x);
				
				r = pow(x*x+y*y+z*z,0.5);

				double A = pdf_1(r,r_n);
				double rand = rnd.Rannyu();

				if( rand <= A){ 
				
					r_n = r;
					x_n = x;	
					y_n = y;
					z_n = z;	
					
				}

				sum_r += r_n;

			}		

			file_N_block << sum_r/ (double)N_block[k] <<endl;
		}

		file_N_block.close();
	}


	// Calcolo ora invece il valor medio del raggio per entrambi i livelli atomici. Per la stima dell'errore devo considerare dei blocchi con 1000 realizzazioni. 
	
	// 1s

	r_n = 0;
	x_n = y_n = z_n = 0;
	
	int M = 100000;
	int L = 500;
	N = M/L;

	double mean_r[N]={0};
	double mean_r2[N]={0};

	ofstream file_points_1;										// Salvo alcune configurazioni così da avere poi una rappresentazione 3D 
	file_points_1.open("file_points_1.txt");

	for (int i = 0; i < N; i+=1){

		for (int j = 0; j<L; j+=1){

			x = rnd.Rannyu(x_n+delta_x,x_n-delta_x);
			y = rnd.Rannyu(y_n+delta_x,y_n-delta_x);
			z = rnd.Rannyu(z_n+delta_x,z_n-delta_x);
			
			r = pow(x*x+y*y+z*z,0.5);

			double A = pdf_1(r,r_n);
			double rand = rnd.Rannyu();

			if( rand <= A){ 
			
				r_n = r;
				x_n = x;	
				y_n = y;
				z_n = z;	
				
			}

			mean_r[i] += r_n;
	
			if(i%10 == 0){
				file_points_1 << x <<","<<y<<","<<z<<endl;
			}


		}		
	}

	
	file_points_1.close();

	for(int j = 0; j<N; j++){
		mean_r[j] = mean_r[j]/L;
		mean_r2[j] = pow(mean_r[j],2);

	}

	// faccio ora le medie progressive per avere una stima dell'errore che commetto. 
	double sumprog_mean2[N]={0};
	double sumprog_mean[N]={0};
	double std_mean[N];

	
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=mean_r2[j];
			sumprog_mean[i]+=mean_r[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}



	// Trascrivo ora i risultati su file
	ofstream file_mean_r;
	file_mean_r.open ("file_mean_r_1.txt");	

	for(int i=0;i<N;i+=1){
		file_mean_r << sumprog_mean[i]<<","<<std_mean[i]<<endl;
	}


	// 2s

	r_n = 0;
	x_n = y_n = z_n = 0;
	
	ofstream file_points_2;										// Salvo alcune configurazioni così da avere poi una rappresentazione 3D 
	file_points_2.open("file_points_2.txt");


	for (int i = 0; i < N; i+=1){

		mean_r[i] = 0;

		for (int j = 0; j<L; j+=1){

			x = rnd.Rannyu(x_n+delta_x,x_n-delta_x);
			y = rnd.Rannyu(y_n+delta_x,y_n-delta_x);
			z = rnd.Rannyu(z_n+delta_x,z_n-delta_x);
			
			r = pow(x*x+y*y+z*z,0.5);

			double A = pdf_2(r,r_n, z,z_n);
			double rand = rnd.Rannyu();

			if( rand <= A){ 
			
				r_n = r;
				x_n = x;	
				y_n = y;
				z_n = z;	
				
			}

			mean_r[i] += r_n;
	
			if(i%10 == 0){
				file_points_2 << x <<","<<y<<","<<z<<endl;
			}

		}		
	}

	file_points_2.close();

	for(int j = 0; j<N; j++){
		mean_r[j] = mean_r[j]/L;
		mean_r2[j] = pow(mean_r[j],2);
		
		sumprog_mean2[j]=0;
		sumprog_mean[j]=0;

	}

	// faccio ora le medie progressive per avere una stima dell'errore che commetto. 
	

	
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=mean_r2[j];
			sumprog_mean[i]+=mean_r[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}



	// Trascrivo ora i risultati su file
	ofstream file_mean_r_2;
	file_mean_r_2.open ("file_mean_r_2.txt");	

	for(int i=0;i<N;i+=1){
		file_mean_r_2 << sumprog_mean[i]<<","<<std_mean[i]<<endl;
	}


	// Valutiamo come cambia il sistema se utilizziamo una distribuzione normale piuttosto che uniforme

	// 1s

	r_n = 0;
	x_n = y_n = z_n = 0;

	delta_x = 0.7;

	for (int i = 0; i < N; i+=1){
		
		 mean_r[i]=0;

		for (int j = 0; j<L; j+=1){

			x = rnd.Gauss(x_n,delta_x);
			y = rnd.Gauss(y_n,delta_x);
			z = rnd.Gauss(z_n,delta_x);
			
			r = pow(x*x+y*y+z*z,0.5);

			double A = pdf_1(r,r_n);
			double rand = rnd.Rannyu();

			if( rand <= A){ 
			
				r_n = r;
				x_n = x;	
				y_n = y;
				z_n = z;	
			}

			mean_r[i] += r_n;

		}		
	}

	for(int j = 0; j<N; j++){
		mean_r[j] = mean_r[j]/L;
		mean_r2[j] = pow(mean_r[j],2);

		sumprog_mean2[j]=0;
		sumprog_mean[j]=0;

	}

	// faccio ora le medie progressive per avere una stima dell'errore che commetto. 
		
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=mean_r2[j];
			sumprog_mean[i]+=mean_r[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}



	// Trascrivo ora i risultati su file
	ofstream file_mean_norm_r_1;
	file_mean_norm_r_1.open ("file_mean_norm_r_1.txt");	

	for(int i=0;i<N;i+=1){
		file_mean_norm_r_1 << sumprog_mean[i]<<","<<std_mean[i]<<endl;
	}


	// 2s

	r_n = 100;
	x_n = y_n = z_n = 0;

	for (int i = 0; i < N; i+=1){

		mean_r[i] = 0;

		for (int j = 0; j<L; j+=1){

			x = rnd.Gauss(x_n,delta_x);
			y = rnd.Gauss(y_n,delta_x);
			z = rnd.Gauss(z_n,delta_x);
			
			r = pow(x*x+y*y+z*z,0.5);

			double A = pdf_2(r,r_n, z,z_n);
			double rand = rnd.Rannyu();

			if( rand <= A){ 
			
				r_n = r;
				x_n = x;	
				y_n = y;
				z_n = z;	
				
			}

			mean_r[i] += r_n;

		}		
	}

	for(int j = 0; j<N; j++){
		mean_r[j] = mean_r[j]/L;
		mean_r2[j] = pow(mean_r[j],2);
		
		sumprog_mean2[j]=0;
		sumprog_mean[j]=0;

	}

	// faccio ora le medie progressive per avere una stima dell'errore che commetto. 
	

	
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=mean_r2[j];
			sumprog_mean[i]+=mean_r[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}



	// Trascrivo ora i risultati su file
	ofstream file_mean_norm_r_2;
	file_mean_norm_r_2.open ("file_mean_norm_r_2.txt");	

	for(int i=0;i<N;i+=1){
		file_mean_norm_r_2 << sumprog_mean[i]<<","<<std_mean[i]<<endl;
	}


	rnd.SaveSeed();
	return 0;
}


double pdf_1(double r, double r_n) {
   
		double pdf = min(1.,exp(-2.*r)/exp(-2.*r_n));

return pdf; 
}

double pdf_2(double r, double r_n, double z, double z_n){
    double alpha;
    double move;
    
    move = ( z * z * exp(- r) ) / ( z_n * z_n * exp(-r_n));
    
    alpha = min(1., move);
    
    return alpha;
}




