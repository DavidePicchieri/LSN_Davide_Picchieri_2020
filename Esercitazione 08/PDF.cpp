#include "PDF.h"
#include <cmath>

using namespace std;

PDF :: PDF(){
    m_x = 0;
    m_y = 0;
    m_z = 0;
}

PDF :: ~PDF(){}

Hydro100 :: Hydro100(){}

Hydro100 :: ~Hydro100(){}

double Hydro100 :: Probability(double x, double y, double z){
    return pow(exp(-sqrt(x*x + y*y + z*z))/sqrt(M_PI),2);
}

Hydro210 :: Hydro210(){}

Hydro210 :: ~Hydro210(){}

double Hydro210 :: Probability(double x, double y, double z){
    return pow(z*exp(-sqrt(x*x + y*y + z*z)/2)*sqrt(2./M_PI)/8.,2);
}

WaveFunction :: WaveFunction(double sigma, double mu){
    m_sigma = sigma;
    m_mu = mu;
}

WaveFunction :: ~WaveFunction(){}

double WaveFunction :: Probability(double x, double y, double z){
    return pow((exp(-(x-m_mu)*(x-m_mu)/(2*m_sigma*m_sigma))+exp(-(x+m_mu)*(x+m_mu)/(2*m_sigma*m_sigma))),2);
}

double WaveFunction :: Eloc(double x){
    return -1./(2.*pow(m_sigma,4))*(x*x+m_mu*m_mu-m_sigma*m_sigma-2*x*m_mu*tanh(m_mu*x/(m_sigma*m_sigma)))+pow(x,4)-5./2.*x*x;
}