/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(int argc, char** argv)
{ 
  //Choice to restart from the initial config or an old one
  if (argc != 3){
    cerr << "Usage: <" << argv[0] << "> <0 to run from old config, any other number to start anew> <iteration>" << endl;
    return -1;
  }
  if (atoi(argv[1]) == 0)
    oldconf = true;
  else 
    oldconf = false;
 
  temp = atoi(argv[2]) * 0.1 + 0.5;
  Input(); //Inizialization
  if (oldconf == false){
    InitialConfig();
  }
  else InitialConfig("config.final");
  
  //equilibration
  for (int i=0; i<nequilibrium; i++){
    Move(metro);
    Measure();
    if (temp == 0.5 || temp == 1.5){
      PrintInstant(directory + "/Equilibration_" + to_string(temp)+".out");
    }
  }

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  
  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  
  if (metro == 1)
    directory = "Metro";
  else 
    directory = "Gibbs";

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) {
    cout << "The program perform Metropolis moves" << endl;
    directory = "Metro";
  }
  else {
    cout << "The program perform Gibbs moves" << endl;
    directory = "Gibbs";
  }
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  if (oldconf == false){
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
    PrintConf("InitialConf.dat");
    cout << "Sto printando" << endl;
  }
  else {
    ReadConf.open("config.final");
    for (int i=0; i<nspin; ++i)
    {
      ReadConf >> s[i];
    }
    ReadConf.close();
  }
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    sm = s[o];

    if(metro==1) //Metropolis
    {
      double alpha = fmin ( 1, exp (beta *(Boltzmann (s[o], o) - Boltzmann (-s[o], o)) )) ;
      double r = rnd.Rannyu() ;
      attempted ++ ;
      if (r < alpha)
      {
          s[o] *= -1 ;
          accepted ++ ;
      }
    }
    else //Gibbs sampling
    {
      double p = 1. /( 1. + exp ( beta * (Boltzmann(1, o) - Boltzmann(-1, o) ) ) );
      double r = rnd.Rannyu() ;
      if (r < p)
          s[o] = 1 ;
      else
          s[o] = -1 ;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}


void Measure()
{
  double u = 0.0, u2=0.0, spin=0.0, spin2=0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     spin += s[i];     
  }
  spin2 = spin*spin;
  u2 = u*u;
  walker[iu] = u;
  walker[ic] = u2;
  walker[ix] = spin2;
  walker[im] = spin;
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
   attempted = 0;
   accepted = 0;
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
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd = 12;
   string ext_field;
   if (h==0){
     ext_field = "0.00";
   }
   else if (h==0.02){
     ext_field = "0.02";
   }
   else {
     ext_field = "something is wrong";
   }
    
    cout << "Block number " << iblk << endl;
    if (metro == 1)
        cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    //Internal energy
    Ene.open(directory + "/output.ene_" + ext_field,ios::app);
    stima_u = blk_av[iu] / blk_norm; 
    glob_av[iu]  += stima_u; 
    glob_av2[iu] += stima_u * stima_u; 
    err_u = Error (glob_av[iu], glob_av2[iu], iblk) ; 
    
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u/(double)nspin << setw(wd) << glob_av[iu]/(double)iblk/(double)nspin << setw(wd) << err_u/(double)nspin << endl;
    Ene.close();
    
    //Heat capacity
    Heat.open(directory + "/output.heat_" + ext_field,ios::app);
    stima_c = blk_av[ic] / blk_norm ; 
    glob_av[ic] += stima_c ; 
    glob_av2[ic] += stima_c * stima_c ; 
    err_c = Error (glob_av[ic], glob_av2[ic], iblk) ; 
    
    Heat << setw(wd) << iblk << setw(wd) << beta*beta * (stima_c - pow (stima_u, 2)) << setw(wd) << beta*beta * (glob_av[ic]/(double)iblk - pow (glob_av[iu]/(double)iblk, 2))/(double)nspin << setw(wd) << beta * beta * sqrt ( pow ( err_c/(double)nspin, 2 ) + pow ( 2 * fabs(glob_av[iu])/(double)iblk/(double)nspin * err_u /(double) nspin, 2 ) ) << endl ;
    Heat.close() ;
    
    //Magnetisation
    Mag.open(directory + "/output.mag_" + ext_field,ios::app);
    stima_m = blk_av[im] / blk_norm ; 
    glob_av[im] += stima_m ; 
    glob_av2[im] += stima_m * stima_m ; 
    err_m = Error (glob_av[im], glob_av2[im], iblk) ; 
    
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk/(double)nspin << setw(wd) << err_m/(double)nspin << endl;
    Mag.close();
    
    //Magnetic susceptibility
    Chi.open(directory + "/output.chi_" + ext_field,ios::app);
    stima_x = blk_av[ix] / blk_norm ; 
    glob_av[ix] += stima_x ; 
    glob_av2[ix] += stima_x*stima_x ; 
    err_x = Error (glob_av[ix], glob_av2[ix], iblk) ; 
    
    if (h==0)
      Chi << setw(wd) << iblk << setw(wd) << beta * (stima_x) << setw(wd) << beta * glob_av[ix]/(double)iblk/(double)nspin  << setw(wd) << beta * err_x /(double)nspin << endl ;
    else
      Chi << setw(wd) << iblk << setw(wd) << beta * (stima_x) << setw(wd) << beta * (glob_av[ix]/(double)iblk - pow(glob_av[im]/(double)iblk, 2))/(double)nspin << setw(wd) <<  beta * sqrt ( pow(err_x, 2) + pow(2 * fabs(glob_av[im])/(double)iblk * err_m, 2) )/(double)nspin << endl ;

    Chi.close() ;
   
    if (iblk==nblk)
    {
        ofstream EneTot, HeatTot, MagTot, ChiTot;
        
        if (h==0){
          EneTot.open(directory + "/output.ene_tot", ios::app) ;
          EneTot << temp << setw(wd) << glob_av[iu]/(double)iblk/(double)nspin << setw(wd) << err_u/(double)nspin << endl; 
          EneTot.close();

          HeatTot.open(directory + "/output.heat_tot", ios::app) ;
          HeatTot << temp << setw(wd) << beta*beta * (glob_av[ic]/(double)iblk - pow (glob_av[iu]/(double)iblk, 2))/(double)nspin << setw(wd) << beta*beta * sqrt ( pow ( err_c / double (nspin), 2 ) + pow ( 2 * fabs(glob_av[iu])/(double)iblk/double(nspin) * err_u, 2 ) ) << endl; 
          HeatTot.close();

          ChiTot.open(directory + "/output.chi_tot", ios::app) ;
          ChiTot << temp << setw(wd) << beta * glob_av[ix]/(double)iblk/(double)nspin << setw(wd) << beta * err_x /(double)nspin << endl; 
          ChiTot.close();
        }
        else {
          MagTot.open(directory + "/output.mag_tot", ios::app) ;
          MagTot << temp << setw(wd) << glob_av[im]/(double)iblk/(double)nspin << setw(wd) << err_m/(double)nspin << endl; 
          MagTot.close();
        }
    }

    cout << "----------------------------" << endl << endl;
}


void ConfFinal()
{
  ofstream WriteConf;

  cout << "Print final configuration to file " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void PrintConf(string filename)
{
  ofstream WriteConf;

  cout << "Print final configuration to file " << filename << endl << endl;
  WriteConf.open(filename);
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();
}

void PrintInstant(string filename)
{
  ofstream Write;
  Write.open(filename, ios::app);
  Write << walker[iu]/(double)nspin << endl ;
  Write.close();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


void InitialConfig(void)
{
  //initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
}

void InitialConfig(std::string filename)
{
  ifstream ReadConf;
  
  //initial configuration
  ReadConf.open(filename);
  for (int i=0; i<nspin; ++i)
  {
    ReadConf >> s[i];
  }
  ReadConf.close();
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
