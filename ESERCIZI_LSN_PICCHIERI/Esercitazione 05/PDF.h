#ifndef _PDF_
#define _PDF_

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

class PDF{
private:
  double m_x, m_y, m_z;

protected:

public:
  // constructors
  PDF();
  // destructor
  virtual ~PDF();
  // methods
  virtual double Probability(double, double, double) = 0;

};

class Hydro100 : public PDF{
    public:
        // constructors
        Hydro100();
        // destructor
        ~Hydro100();
        // methods
        virtual double Probability(double, double, double);


};

class Hydro210 : public PDF{
    public:
        // constructors
        Hydro210();
        // destructor
        ~Hydro210();
        // methods
        virtual double Probability(double, double, double);


};

#endif 
