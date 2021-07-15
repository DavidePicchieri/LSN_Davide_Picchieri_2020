#include "functions.h"

using namespace std;

double Slope(double x1, double y1, double x2, double y2){
    return (y2-y1)/(x2-x1);
}

double Y_Intersection(double m, bool halfplane, double r){

    if (halfplane == true){return abs(r*m/(sqrt(m*m+1)));}
    else {return -abs(r*m/(sqrt(m*m+1)));} 

}

double Price(double t, double S0, double mu, double sigma, double Wt){
    return S0 * exp ((mu - sigma * sigma / 2.) * t + sigma * Wt);
}

double Error(double ave, double ave2, double block_number){
    if (block_number == 0)
        return 0;
    else
        return sqrt( ((ave2 / (block_number + 1)) - pow((ave / (block_number + 1)),2)) / (block_number));
}
