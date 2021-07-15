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
    int populationsize = 100;
    int generations = 5000;
    Species pop(populationsize, p_rnd, ncities);
    if (circle == true){      //Initialize the cities
        pop.SetPosCircle(); 
    }
    else{
        pop.SetPosSquare();
    }
    pop.Measure();
    pop.Sort();

    for (int i=0; i<populationsize; i++){
        if (circle == true){
            pop.Print(i, "CircleInitialConfiguration.out");
        }
        else{
            pop.Print(i, "SquareInitialConfiguration.out");
        }
    }
    
    for(int i=0; i<generations; i++){
        pop.Selection(2); //The argument of the function is the powerlaw of the selection
        pop.Crossover(16, 0.50); // The first argument is at what point of the trip the crossover should happen, and the second is the chance of the crossover happening at all
        pop.Mutation(0.05); //The argument is the chance of each of the four mutations to happen
        pop.Measure(); //Measures the length of the trips
        pop.Sort();
        if (circle == true){
            pop.PrintDistance(0, "BestTripGenerationCircle.out"); //First argument is which trip of thepopulation print
            pop.PrintAve(50, "AveTripGenerationCircle.out"); //First argument is over how many of the first n trip to average the distance
        }
        else{
            pop.PrintDistance(0, "BestTripGenerationSquare.out");
            pop.PrintAve(50, "AveTripGenerationSquare.out");
        }
    }
    if (circle == true){
        pop.PrintBestSA("CircleBestTrip.out"); //Prints the configuration of the best trip of the population
    }
    else{
        pop.PrintBestSA("SquareBestTrip.out");
    }

    return 0;

}
