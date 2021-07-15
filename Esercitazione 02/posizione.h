#ifndef __posizione__
#define __posizione__

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

class Posizione {

private:
  double m_x, m_y, m_z;

protected:

public:
  // constructors
  Posizione();
  Posizione(double, double, double);
  // destructor
  ~Posizione();
  // methods
  void SetPosizione(double, double, double);
  double GetX() {return m_x;}
  double GetY() {return m_y;}
  double GetZ() {return m_z;}
  void PassoCart(double, double, double);
  void PassoSfer(double, double, double);
  double DistanzaDaPunto(double, double, double);
  double DistanzaQuadDaPunto(double, double, double);

};

#endif 
