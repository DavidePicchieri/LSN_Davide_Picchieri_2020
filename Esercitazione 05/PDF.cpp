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