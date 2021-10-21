/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include <string>
#include <iomanip>

using namespace std;

int main(int argc, char ** argv){ 
  Input(argc, argv); 
	Equilibration();            
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {

    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages

    }
         
		
		Save_average(iblk); //Print results for current block
  }
  //ConfFinal(); //Write final configuration

  return 0;
}




void Input(int argc, char ** argv){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  if(argc!=2) 
		cout<<" 	ERROR: solid or liquid or gE"<<endl;
	else{
		input = argv[1];
		ReadInput.open("input."+input); //Read input
	}

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
	ReadInput >> nblk;
  ReadInput >> iprint;
	ReadInput >> start_old;
	ReadInput >> rescale_position;
	ReadInput >> equilibration_time;
	
	cout << "Number of blocks = " << nblk << endl;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
	if(start_old==1){
		cout << "Warning! The program will restart with an old spatial configuration" << endl << endl;
		ReadInput.close();
	}
	
 
  bin_size = (box/2.0)/(double)nbins;
	

//Read initial configuration
	if(start_old==0){
		cout << "Read initial configuration from file config.0 " << endl << endl;
		ReadConf.open("config.0");
		for (int i=0; i<npart; ++i){
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = x[i] * box;
			y[i] = y[i] * box;
			z[i] = z[i] * box;
		}
		ReadConf.close();
	
//Prepare initial velocities
	 cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	 double sumv[3] = {0.0, 0.0, 0.0};
	 for (int i=0; i<npart; ++i){
	   vx[i] = rand()/double(RAND_MAX) - 0.5;
	   vy[i] = rand()/double(RAND_MAX) - 0.5;
	   vz[i] = rand()/double(RAND_MAX) - 0.5;

	   sumv[0] += vx[i];
	   sumv[1] += vy[i];
	   sumv[2] += vz[i];
	 }
	 for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	 double sumv2 = 0.0, fs;
	 for (int i=0; i<npart; ++i){
	   vx[i] = vx[i] - sumv[0];
	   vy[i] = vy[i] - sumv[1];
	   vz[i] = vz[i] - sumv[2];

	   sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	 }
	 sumv2 /= (double)npart;

	 fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	 for (int i=0; i<npart; ++i){
	   vx[i] *= fs;
	   vy[i] *= fs;
	   vz[i] *= fs;

	   xold[i] = Pbc(x[i] - vx[i] * delta);
	   yold[i] = Pbc(y[i] - vy[i] * delta);
	   zold[i] = Pbc(z[i] - vz[i] * delta);
	 }
	}

	// -------> START OLD <--------
	else{
			ReadConf.open("config."+input+".final");

			for (int i=0; i<npart; ++i){
				ReadConf >> x[i] >> y[i] >> z[i];
				x[i] = x[i] * box;
				y[i] = y[i] * box;
				z[i] = z[i] * box;
			}
			ReadConf.close();

			ReadConf.open("old."+input+".final");
			
			for (int i=0; i<npart; ++i){
				ReadConf >> xold[i] >> yold[i] >> zold[i];
				xold[i] = xold[i] * box;
				yold[i] = yold[i] * box;
				zold[i] = zold[i] * box;
			}
			ReadConf.close();
			
			for (int i=0; i<npart; ++i){
			 vx[i] = Pbc(x[i] - xold[i])/(2.0 * delta);
			 vy[i] = Pbc(y[i] - yold[i])/(2.0 * delta);
			 vz[i] = Pbc(z[i] - zold[i])/(2.0 * delta);
			}
	}

	return;
	}


void Equilibration (void){
		cout << "----------------------------" << endl << endl;
	 	cout << "Equilibration time:"<< equilibration_time <<endl<<endl;
	 	cout << "----------------------------" << endl << endl;

		for (int i=0; i<=equilibration_time; i++) {
			Move();	}
	
		equilibration_time=0;
	
}

void Move(){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

	if((rescale_position == 1) & (equilibration_time !=0)){
			double t = 0;
			double correction_factor;

			for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
			stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
			correction_factor = temp/stima_temp;
			for (int i=0; i<npart; ++i){
				vx[i]=pow(correction_factor,0.5)*vx[i];
				vy[i]=pow(correction_factor,0.5)*vy[i];
				vz[i]=pow(correction_factor,0.5)*vz[i];
				xold[i] = Pbc(x[i] - vx[i] * delta);
			 	yold[i] = Pbc(y[i] - vy[i] * delta);
			 	zold[i] = Pbc(z[i] - vz[i] * delta);
			}
	}

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }

  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
 
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pres.open("output_pres.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
	p = 0.0;
	
//reset the hystogram of g(r)
  for (int k=0; k<nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

		//potential energy
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
			 v += vij;
			  
     }

		// Pressure
    pij = 48 * (pow(1 / dr, 12) - 0.5 * pow(1 / dr, 6));
    p += pij;

		//g(r)		
	  for (int k =0; k<nbins; k+=1){
				if( (bin_size*k<dr) && (dr<bin_size*(k+1)) ){
					 walker[k]+=2;

				}
	  }

       
    }          
  }

	
//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
		stima_pres = (rho * t + 1. / (3.*vol) * p) / (double)npart;
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
		Pres << stima_pres << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
		Pres.close();

    return;
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<nbins; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];

   }
   blk_norm = blk_norm + 1.0;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << endl << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

	cout << "Print old configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();


  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<nbins; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<nbins; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Block_average(const char *file_in, const char *file_out, int nblk){
	
	int L = nstep;

	int N = nblk;

	ifstream file_read;
	file_read.open(file_in);
	
	double ave[N] = {0};
	double ave2[N] = {0};
	double x;

	int i = 0;
	int j = 0;

	// Faccio ora le medie a blocchi

	while (file_read >> x) {

      ave[i] += x;
			j++;
			if(j == L){
				ave[i] = ave[i]/(double)L;
				ave2[i] = pow(ave[i],2);
				j = 0;
				i += 1;
			}
    }
	// faccio ora le medie progressive per avere una stima dell'errore che commetto. 



	double sumprog_mean2[N]={0};
	double sumprog_mean[N]={0};
	double std_mean[N];

	
	for(int i =0; i<N; i+=1){
		for(int j =0; j<=i;j+=1){
			sumprog_mean2[i]+=ave2[j];
			sumprog_mean[i]+=ave[j];
		}
		sumprog_mean2[i]=sumprog_mean2[i]/(i+1);
		sumprog_mean[i]=sumprog_mean[i]/(i+1);
		std_mean[i]=pow((sumprog_mean2[i]-pow(sumprog_mean[i],2))/(i+1),0.5);
	}

	// Trascrivo ora i risultati su file
	ofstream file_write;
	file_write.open (file_out);

	for(int i=0;i<N;i+=1){
		file_write << sumprog_mean[i]<<","<<std_mean[i]<<endl;
	}
	
	file_write.close();

	
	remove( file_in );


}

void Save_average(int iblk){

cout<< "Current block   "<< iblk<<endl;

	if (iblk == nblk){
		Block_average("output_epot.dat", "ave_epot.out", nblk);
		Block_average("output_ekin.dat", "ave_ekin.out",nblk);
		Block_average("output_etot.dat", "ave_etot.out",nblk);
		Block_average("output_temp.dat", "ave_temp.out",nblk);
		Block_average("output_pres.dat", "ave_pres.out",nblk);
	}
    	
   ofstream Gave;
   const int wd=12;
    
    Gave.open("output.gave.0");
  
		for (int i = 0; i<nbins; i+=1){

      double r = bin_size*i;
			stima_gofr = blk_av[i]/(blk_norm*(rho * npart * (4. * M_PI / 3.) * (pow(r + bin_size, 3) - pow(r, 3))));
			glob_av[i] += stima_gofr;
    	glob_av2[i] += stima_gofr*stima_gofr;
			

			
			if (iblk==nblk){
				err_gave=Error(glob_av[i],glob_av2[i],iblk);
				Gave << setw(wd) << r <<  setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_gave << endl;
			}
		}


    cout << "----------------------------" << endl << endl;

    
}



double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
