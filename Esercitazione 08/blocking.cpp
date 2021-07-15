#include "blocking.h"

using namespace std;

Blocking :: Blocking(){
    m_nblocks = 0;
    m_length = 0;
    m_sum = 0;
    m_sum2 = 0;
    m_ave = 0;
    m_count = 0;
    m_appo = 0;
    m_error = 0;
}

Blocking :: Blocking(int neval, int nblocks){
    if (neval % nblocks != 0){
        cout << "Number of Evaluations MUST be a multiple of the Number of Blocks: It wasn't so the object is initialized as 0,0" << endl;
        m_nblocks = 0;
        m_neval = 0;
        m_length = 0;
    }
    else {
        m_nblocks = nblocks;
        m_neval = neval;
        m_length = neval / nblocks;
    }
    m_sum = 0;
    m_sum2 = 0;
    m_ave = 0;
    m_count = 0;
    m_appo = 0;
    m_error = 0;
}

Blocking :: ~Blocking(){}

void Blocking :: Add(double eval,ofstream & write){
    m_appo += eval;
    m_count ++;
    if (m_count % m_length == 0){
        m_sum += m_appo / m_length;
        m_sum2 += pow (m_appo / m_length, 2);
        m_error = Error(m_sum, m_sum2, (m_count / m_length - 1));
        m_appo = 0;
        m_ave = m_sum / (m_count / m_length);
        write <<  m_count / m_length  << " \t " << m_ave << " \t " << m_error << endl;
    }
}
void Blocking :: Add(double eval){
    m_appo += eval;
    m_count ++;
    if (m_count % m_length == 0){
        m_sum += m_appo / m_length;
        m_sum2 += pow (m_appo / m_length, 2);
        m_error = Error(m_sum, m_sum2, (m_count / m_length - 1));
        m_appo = 0;
        m_ave = m_sum / (m_count / m_length);
    }
}
void Blocking :: SetNblocks(int nblocks){
    m_nblocks = nblocks;
}

void Blocking :: SetLength(int length){
    m_length = length;
}

void Blocking :: SetNeval(int neval){
    m_neval = neval;
}

void Blocking :: SetCount(int count){
    m_count = count;
}

void Blocking :: Reset(){
    m_count = 0;
    m_ave = 0;
    m_sum = 0;
    m_sum2 = 0;
    m_error = 0;
    m_appo = 0;
}