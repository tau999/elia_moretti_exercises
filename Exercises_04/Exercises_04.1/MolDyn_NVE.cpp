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

using namespace std;

int main(int argc, char ** argv){ 
  Input(argc, argv);             //Inizialization
	Equilibration();			//Equilibration
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%1 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
		 Save_average();
  }
  ConfFinal();         //Write final configuration to restart

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
		cout<<" 	ERROR: 0 -> CIRCLE    &      1 -> SQUARE"<<endl;
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
	ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> iprint;
	ReadInput >> rescale_position;
	ReadInput >> start_old;
	ReadInput >> equilibration_time;
	

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

	if(start_old == 0){
		//Read initial configuration
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
	return;
}

void Move(void){ //Move particles with Verlet algorithm
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

  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
	ofstream WriteConf_old;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config."+input+".final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
 WriteConf.close();

	cout << "Print old configuration to file old.final " << endl << endl;
  WriteConf_old.open("old."+input+".final");

  for (int i=0; i<npart; ++i){
    WriteConf_old << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf_old.close();

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

void Save_average(void){
		Block_average("output_epot.dat", "ave_epot.out",nblk);
		Block_average("output_ekin.dat", "ave_ekin.out",nblk);
		Block_average("output_etot.dat", "ave_etot.out",nblk);
		Block_average("output_temp.dat", "ave_temp.out",nblk);
		Block_average("output_pres.dat", "ave_pres.out",nblk);
}

void Block_average(const char *file_in, const char *file_out, int nblk){
	
	int L = nstep/nblk;
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
		file_write << i <<","<< sumprog_mean[i]<<","<<std_mean[i]<<endl;
	}
	
	file_write.close();

	


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
