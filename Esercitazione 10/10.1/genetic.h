#ifndef __genetic_h_
#define __genetic_h_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "random.h"

using namespace std;

class Specimen{

private:
  vector<int> m_cities; //Vector containing the cities
  int         m_gen; //Number of generation
  Random *    m_rnd;
  double      m_distance;
  int         Pbc(int i);

public:
  //costructor
  Specimen(int, Random *);
  // destructor
  ~Specimen();
  
  Specimen& operator=(Specimen&);

  void Check();
  void Print(string);
  void SetDistance(double dist){m_distance=dist;}
  
  vector<int> GetCities(){return m_cities;}
  void SetCity(int index, int city){m_cities[index]=city;}
  double GetDistance(){return m_distance;}
  
  //Mutations
  void Swap(int);
  void Shift(int , int , int);
  void SwapGroup(int, int);
  void Inversion(int, int);

};


class Species{

private:

  int                    m_popsize, m_ngen, m_ncities, m_bestgen;
  double                 m_mindist, m_temp, m_ann_coeff;
  Random *               m_rnd;
  vector<Specimen>       m_pop, m_popnew;
  vector<vector<double>> m_pos;
  vector<int>            m_besttrip;
  
public:
  //costructor
  Species(int n, Random * rnd, int dim);
  // destructor
  ~Species();
  
  void Print(int, string);
  void PrintBest(string);
  void PrintBestSA(string);
  void PrintDistance(string);
  void PrintDistance(int, string);
  void PrintAve(int, string);
  void Copy(int i, int j);  
  void SetPosCircle();
  void SetPosSquare();
  void SetCoeff(double);
  void SetTemp(double);
  double GetTemp(){return m_temp;};
  double GetAlpha();
  void Measure();
  void Sort();
  void Selection(double);
  void Crossover(int, double);
  void Mutation(double);
  void Annealing();
  void Accept();
  void Rearrange();
  
  
  
  
};


#endif