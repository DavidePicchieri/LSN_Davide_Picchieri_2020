#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "genetic.h"

using namespace std;


int main (int argc, char *argv[]){
    bool circle;
    if (argc != 2){
        cerr << "Usage: <" << argv[0] << "> <0 to set the cities on a cirle, any other number to set them randomly in a square>" << endl;
        return -1;
    }
    if (atoi(argv[1]) == 0)
        circle = true;
    else 
        circle = false;


    Random rnd("Primes","seed.in");
    Random * p_rnd = & rnd;
  
    int ncities = 32;
    int populationsize = 2; //The "old" one and the "new" one which might take its place
    double t_init = 10;  //Initial temperature for the annealing
    double t_fin = 0.001; //Final temperature
    double ann_coeff = 0.99999; //Coeffincient which regulates the cooling: T(i+1) = T(i) * ann_coeff
    Species pop(populationsize, p_rnd, ncities);
    if (circle == true){
        pop.SetPosCircle();
    }
    else{
        pop.SetPosSquare();
    }
    pop.Measure();  //Measure the trip length
    pop.SetTemp(t_init); //Set the temperature
    pop.SetCoeff(ann_coeff); //Set the coefficient

    for (int i=0; i<populationsize; i++){ // Print initial configuration
        if (circle == true){
            pop.Print(i, "CircleInitialConfiguration.out");
        }
        else{
            pop.Print(i, "SquareInitialConfiguration.out");
        }
    }
    
    do{
        pop.Annealing(); // Multiply the old temperature fot the coefficient
        pop.Rearrange(); // Performs a swap between cities and an inversion
        pop.Measure(); //Measure the effect of rearrange 
        
        if (circle == true){
            pop.PrintDistance(0, "BestTripGenerationCircle_SA.out");
        }
        else{
            pop.PrintDistance(0, "BestTripGenerationSquare_SA.out");
        }
        
        pop.Accept(); //Decides if the new trip is to be accepted
        pop.Measure(); //Measure again the distance of the 'surviving' trip
        cout << pop.GetTemp() << endl; //To check on the progression of the code
    }
    while (pop.GetTemp() > t_fin);

    if (circle == true){
        pop.PrintBestSA("CircleBestTrip_SA.out"); //Print best configuration
    }
    else{
        pop.PrintBestSA("SquareBestTrip_SA.out");
    }

    return 0;

}
