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

//Constants to go from LJ units to SI units
double sigma=0.34; //nanometers
double boltzmann =1.38064852E-23; //Joule/Kelvin
double epsilon=boltzmann*120; //Joule
double m=39.948; //atomic mass units

int main(int argc, char** argv){ 
  //Choice to restart from the initial config or an old one
  if (argc != 3){
    cerr << "Usage: <" << argv[0] << "> <1 to run from old config, 0 start anew> <1 to rescale velocities>" << endl;
    return -1;
  }
  if (atoi(argv[1]) == 1)
    readoldconf = true;
  else 
    readoldconf = false;

  if (atoi(argv[1]) == 1)
    rescale = true;
  else 
    rescale = false;


  Input();             //Inizialization

  //blocking method objects
  Blocking Epotave(nstep/10, nblocks);
  Blocking Ekinave(nstep/10, nblocks);
  Blocking Etotave(nstep/10, nblocks);
  Blocking Tempave(nstep/10, nblocks);
  

  for(int iblk=1; iblk <= nblocks; ++iblk)
  {
    Reset(iblk);
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0) Measure(Epotave, Ekinave, Etotave, Tempave);
      if(istep%10 == 0) Accumulate();

    }
    Averages(iblk);
  }

  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf, ReadOldConf, ReadOlderConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

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
  ReadInput >> nblocks;
  ReadInput >> iprint;
  ReadInput >> state;


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
   

//measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;



  if (readoldconf == false){
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
  else {
    //Reading the configurations at t and t-dt from old config file
    double sumv2 = 0.0, fs;
    cout << "Read initial configuration from file config.final and old config from old.final " << endl << endl;
    ReadOldConf.open("config.final");
    for (int i=0; i<npart; ++i){
      ReadOldConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;

    }
    ReadConf.close();

    ReadOlderConf.open("old.final");
    for (int i=0; i<npart; ++i){
      ReadOlderConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;

    }
    Move();
    for (int i=0; i<npart; ++i){
      vx[i] = (x[i] - xold[i])/(delta);
      vy[i] = (y[i] - yold[i])/(delta);
      vz[i] = (z[i] - zold[i])/(delta);
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }

    sumv2 /= (double)npart;

    if (rescale == true){
      fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
      for (int i=0; i<npart; ++i){
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs; 
         
        xold[i] = Pbc(x[i] - delta * vx[i]);
        yold[i] = Pbc(y[i] - delta * vy[i]);
        zold[i] = Pbc(z[i] - delta * vz[i]);
      }
    }

  }

  //reset the hystogram of g(r)
  for (int k=5; k<igofr+nbins; ++k) walker[k]=0.0;

   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
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

void Measure(Blocking & Epotave, Blocking & Ekinave, Blocking & Etotave, Blocking & Tempave){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Epotblock, Ekinblock, Etotblock, Tempblock;



  Epot.open("output_epot_" + state + ".dat",ios::app);
  Ekin.open("output_ekin_" + state + ".dat",ios::app);
  Temp.open("output_temp_" + state + ".dat",ios::app);
  Etot.open("output_etot_" + state + ".dat",ios::app);

  Epotblock.open("ave_epot_" + state + ".dat",ios::app);
  Ekinblock.open("ave_ekin_" + state + ".dat",ios::app);
  Tempblock.open("ave_temp_" + state + ".dat",ios::app);
  Etotblock.open("ave_etot_" + state + ".dat",ios::app);

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

//update of the histogram of g(r)
      for(int k=0; k<nbins; k++){
        if(dr<(k+1)*bin_size){
          walker[igofr+k]+=2;
          break;
        }
      }

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

    Epotave.Add(stima_pot*epsilon, Epotblock);
    Ekinave.Add(stima_kin*epsilon, Ekinblock);
    Etotave.Add(stima_etot*epsilon, Etotblock);
    Tempave.Add(stima_temp*epsilon/boltzmann, Tempblock);

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Epotblock.close();
    Ekinblock.close();
    Tempblock.close();
    Etotblock.close();


    return;
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   double r, gdir, delta_vol;
   double gave[nbins];
   ofstream Gofr, Gave;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output." + state + ".gave.0",ios::app);
    
//g(r)
//calcolare i valori medi di blocco per la g(r) per ogni bin e le rispettive fluttuazioni statistiche mano a mano che la simulazione procede
    
    for (int i=0; i<nbins; i++){
      r = i*bin_size;
      Gofr << "Blocco " << i+1 << " di " << nbins << endl;
      delta_vol = 4.*M_PI*(pow(r+bin_size,3)-pow(r,3))/3.;
      gave[i] = blk_av[igofr+i]/(blk_norm*rho*npart*delta_vol);
      glob_av[igofr+i] += gave[i];
      glob_av2[igofr+i] += gave[i]*gave[i];
      Gofr << gave[i] << setw(wd);
      if(iblk==nblocks){
        err_gdir=Error(glob_av[igofr+i],glob_av2[igofr+i],iblk);
        Gave << r << setw(wd) << glob_av[igofr+i]/(double)iblk << setw(wd) << err_gdir << endl;
      }
    }
    Gofr << endl;
    cout << "----------------------------" << endl << endl;
    Gofr.close();
    Gave.close();
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf, WriteOld;

  cout << "Print final configuration to file config.final " << endl << endl;
  cout << "Print final configuration to file old.final " << endl << endl;

  WriteConf.open("config.final");
  WriteOld.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  WriteOld.close();
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

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
