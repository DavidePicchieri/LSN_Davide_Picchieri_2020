#ifndef __metropolis__
#define __metropolis__

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "random.h"
#include "PDF.h"


using namespace std;

class Metropolis {

private:
  double m_x, m_y, m_z, m_alpha, m_alpha_sum;
  Random * m_rnd;
  PDF * m_pdf;


protected:

public:
  // constructors
  Metropolis();
  Metropolis(double, double, double, Random *, PDF *);
  Metropolis(double, Random *, PDF *);
  // destructor
  ~Metropolis();
  // methods
  void SetPosizione(double, double, double);
  void Reset();
  double GetX() {return m_x;}
  double GetY() {return m_y;}
  double GetZ() {return m_z;}
  double GetR() {return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);}
  double GetAlpha() {return m_alpha;}
  double GetAlphaSum() {return m_alpha_sum;}
  double GetEloc();
  void PassoCartUniforme(double a);
  void PassoCartGauss(double sigma);
  void Passo1DUniforme(double a);
};

#endif 
