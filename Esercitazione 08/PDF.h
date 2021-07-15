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
  double m_sigma, m_mu;

protected:

public:
  // constructors
  PDF();
  // destructor
  virtual ~PDF();
  // methods
  virtual double Probability(double, double, double) = 0;
  virtual double Eloc(double) = 0;
  void SetSigma(double sigma) {m_sigma=sigma;};
  void SetMu(double mu) {m_mu=mu;};

};

class Hydro100 : public PDF{
    public:
        // constructors
        Hydro100();
        // destructor
        ~Hydro100();
        // methods
        virtual double Probability(double, double, double);
        virtual double Eloc(double) {return 0;};
};

class Hydro210 : public PDF{
    public:
        // constructors
        Hydro210();
        // destructor
        ~Hydro210();
        // methods
        virtual double Probability(double, double, double);
        virtual double Eloc(double) {return 0;};
};

class WaveFunction : public PDF{
    public:
        // constructors
        WaveFunction(double, double);
        // destructor
        ~WaveFunction();
        // methods
        double GetSigma() {return m_sigma;};
        double GetMu() {return m_mu;};
        virtual double Probability(double, double, double);
        virtual double Eloc(double);
};

#endif 
