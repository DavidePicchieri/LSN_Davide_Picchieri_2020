#include "posizione.h"

using namespace std;

Posizione :: Posizione(){}

Posizione :: Posizione(double a, double b, double c){
    m_x = a;
    m_y = b;
    m_z = c;
}

Posizione :: ~Posizione(){}

void Posizione :: SetPosizione(double a, double b, double c){
    m_x = a;
    m_y = b;
    m_z = c;
}

void Posizione :: PassoCart(double a, double b, double c){
    m_x += a;
    m_y += b;
    m_z += c;
}
  
void Posizione :: PassoSfer(double r, double phi, double theta){
    m_x += r*sin(theta)*sin(phi);
    m_y += r*sin(theta)*cos(phi);
    m_z += r*cos(theta);
}

double Posizione :: DistanzaDaPunto(double a, double b, double c){
    return sqrt((m_x-a)*(m_x-a) + (m_y-b)*(m_y-b) + (m_z-c)*(m_z-c));
}
double Posizione :: DistanzaQuadDaPunto(double a, double b, double c){
    return (m_x-a)*(m_x-a) + (m_y-b)*(m_y-b) + (m_z-c)*(m_z-c);
}
