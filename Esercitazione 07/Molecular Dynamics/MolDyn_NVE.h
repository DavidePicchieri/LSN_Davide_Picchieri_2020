/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "functions.h"
#include "blocking.h"

//parameters, observables
const int m_props=1000;
int n_props, igofr, nbins;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, bin_size;
double walker[m_props], err_gdir;

// averages
double accepted,attempted;
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;
string state;

// simulation
int nstep, nblocks, iprint, seed;
double delta;

//read from old file
bool readoldconf, rescale;


//functions
void Input(void);
void Move(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(Blocking &, Blocking &, Blocking &, Blocking &);
double Force(int, int);
double Pbc(double);
void GofR(void);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
