#ifndef __blocking__
#define __blocking__

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "functions.h"

using namespace std;

class Blocking {

private:
  int m_nblocks, m_neval, m_length, m_count;
  double m_sum, m_sum2, m_ave, m_error, m_appo;
protected:

public:
  // constructors
  Blocking();
  Blocking(int neval, int nblocks);
  // destructor
  ~Blocking();
  // methods
  void Reset();
  void Add(double add, ofstream & write);
  void Add(double add);
  double GetAverage() {return m_ave;}
  double GetSum() {return m_sum;}
  double GetSum2() {return m_sum2;}
  double GetNblocks() {return m_nblocks;}
  double GetNeval() {return m_neval;}
  double GetLength() {return m_length;}
  double GetCount() {return m_count;}
  double GetError() {return m_error;}
  void SetNblocks(int);
  void SetNeval(int);
  void SetLength(int);
  void SetCount(int);
};














#endif
