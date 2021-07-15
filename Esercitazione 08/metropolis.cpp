#include "metropolis.h"

using namespace std;

Metropolis :: Metropolis(){
    m_x = 0;
    m_y = 0;
    m_z = 0;
    m_alpha = 0;
    m_alpha_sum = 0;
    m_rnd = NULL;
    m_pdf = NULL;
}

Metropolis :: Metropolis(double a, double b, double c, Random * rnd, PDF * pdf){
    m_x = a;
    m_y = b;
    m_z = c;
    m_alpha = 0;
    m_alpha_sum = 0;
    m_rnd = rnd;
    m_pdf = pdf;
}

Metropolis :: Metropolis(double a, Random * rnd, PDF * pdf){
    m_x = a;
    m_alpha = 0;
    m_alpha_sum = 0;
    m_rnd = rnd;
    m_pdf = pdf;
}

Metropolis :: ~Metropolis(){}

void Metropolis :: SetPosizione(double a, double b, double c){
    m_x = a;
    m_y = b;
    m_z = c;
}

void Metropolis :: Reset(){
    m_alpha = 0;
    m_alpha_sum = 0;
}


void Metropolis :: PassoCartUniforme(double a){
    double x = m_x + m_rnd->Rannyu(-a,a);
    double y = m_y + m_rnd->Rannyu(-a,a);
    double z = m_z + m_rnd->Rannyu(-a,a);

    //cout << m_x << "\t" << x << endl;

    m_alpha = fmin(m_pdf->Probability(x, y, z) / m_pdf->Probability(m_x, m_y, m_z), 1);
    m_alpha_sum += m_alpha;
    if (m_rnd->Rannyu() < m_alpha){
        m_x = x;
        m_y = y;
        m_z = z;
    }
}

void Metropolis :: PassoCartGauss(double sigma){
    double x = m_x + m_rnd->Gauss(0, sigma);
    double y = m_y + m_rnd->Gauss(0, sigma);
    double z = m_z + m_rnd->Gauss(0, sigma);

    m_alpha = fmin(m_pdf->Probability(x, y, z) / m_pdf->Probability(m_x, m_y, m_z), 1);
    m_alpha_sum += m_alpha;

    if (m_rnd->Rannyu() < m_alpha){
        m_x = x;
        m_y = y;
        m_z = z;
    }
}

void Metropolis :: Passo1DUniforme(double a){
    double x = m_x + m_rnd->Rannyu(-a,a);

    m_alpha = fmin(m_pdf->Probability(x, 0, 0) / m_pdf->Probability(m_x, 0, 0), 1);
    m_alpha_sum += m_alpha;
    if (m_rnd->Rannyu() < m_alpha){
        m_x = x;
    }
}

double Metropolis :: GetEloc(){
    return m_pdf->Eloc(m_x);
}