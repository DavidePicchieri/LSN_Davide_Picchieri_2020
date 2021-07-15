#include "functions.h"

using namespace std;

double Slope(double x1, double y1, double x2, double y2){
    return (y2-y1)/(x2-x1);
}

double Y_Intersection(double m, bool halfplane, double r){

    if (halfplane == true){return abs(r*m/(sqrt(m*m+1)));}
    else {return -abs(r*m/(sqrt(m*m+1)));} 

}