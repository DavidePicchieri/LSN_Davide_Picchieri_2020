#ifndef __functions_h_
#define __functions_h_

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "random.h"

using namespace std;

double Slope(double x1, double y1, double x2, double y2);
double Y_Intersection(double m, bool halfplane, double r);
double Price(double t, double S0, double mu, double sigma, double Wt);
double Error(double ave, double ave2, double block_number);
void SetGenerator(Random *);
#endif