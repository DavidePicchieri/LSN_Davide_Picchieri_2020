#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "functions.h"
#include "random.h"
#include "genetic.h"
#include "mpi.h"


using namespace std;


int main (int argc, char *argv[]){

    int ncities = 32;
    int populationsize = 100;
    int generations = 10000;
    int migrations = 600;


    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ParallelRandom Par; //setting random generator so that each core uses a different seed
    Random * rnd;
    if(rank==0)
        rnd=Par.SetGen("Primes", "seed.in", 0);
    if(rank==1)
        rnd=Par.SetGen("Primes", "seed.in", 1);
    if(rank==2)
        rnd=Par.SetGen("Primes", "seed.in", 2);
    if(rank==3)
        rnd=Par.SetGen("Primes", "seed.in", 3);


    Species pop(populationsize, rnd, ncities);
    
    pop.SetPosSquare(); //initialize the cities

    double x[ncities];
    double y[ncities];
    
    if(rank==0){         //Save the initial trip in two vector which are then shared with the other cores to have the same starting point
        for(int i=0; i<ncities; i++){
        x[i]=pop.GetPos(0, i);
        y[i]=pop.GetPos(1, i);
        }
    }
    
    MPI_Bcast(x, ncities, MPI_REAL8, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, ncities, MPI_REAL8, 0, MPI_COMM_WORLD);
    
    if(rank!=0){ 
        for(int i=0; i<ncities; i++){
        pop.SetPos(0, i, x[i]);
        pop.SetPos(1, i, y[i]);
        }
    }


    pop.Measure(); 
    pop.Sort();


    if(rank==0)pop.PrintDistance("Square_Initial_Distance.r0.out"); //Print initial distance
    if(rank==1)pop.PrintDistance("Square_Initial_Distance.r1.out");
    if(rank==2)pop.PrintDistance("Square_Initial_Distance.r2.out");
    if(rank==3)pop.PrintDistance("Square_Initial_Distance.r3.out");
    
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0; i<generations; i++){ 
        pop.Selection(2); //Using the same algorithm and parameters as in exercise 9
        pop.Crossover(16, 0.50);
        pop.Mutation(0.05);
        pop.Measure();
        pop.Sort();

        if(rank==0){
            pop.PrintDistance(0, "Square_Best_Dist.r0.out"); //print the best of each rank and the average of the best half
            pop.PrintAve (populationsize/2, "Square_Ave_Dist.r0.out");
        }
        if(rank==1){
            pop.PrintDistance(0, "Square_Best_Dist.r1.out");
            pop.PrintAve (populationsize/2, "Square_Ave_Dist.r1.out");
        }
        if(rank==2){
            pop.PrintDistance(0, "Square_Best_Dist.r2.out");
            pop.PrintAve (populationsize/2, "Square_Ave_Dist.r2.out");
        }
        if(rank==3){
            pop.PrintDistance(0, "Square_Best_Dist.r3.out");
            pop.PrintAve (populationsize/2, "Square_Ave_Dist.r3.out");
        }

        if( (i+1) % migrations == 0){ //Every 'migrations' generations perform the migration of the best trip
            int vec[2];
            int target;
            if(rank==0) 
                pop.MigrationDirection(vec, target); //deciding which rank communicates with which one and then communicating the decision to the others

            MPI_Bcast(&target, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); 
            MPI_Bcast(vec, 2, MPI_INTEGER, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            pop.Migration(rank, target, vec);  //Performs the migration
        }
    }

    if(rank==0) //print the best path configuration for each rank
        pop.PrintBestRank("Square_BestPath.r0.out");
    if(rank==1)
        pop.PrintBestRank("Square_BestPath.r1.out");
    if(rank==2)
        pop.PrintBestRank("Square_BestPath.r2.out");
    if(rank==3)
        pop.PrintBestRank("Square_BestPath.r3.out");

    pop.Best(rank); //find which one is the best path among the four ranks
    if(rank==0) //Print the best path altogether
        pop.PrintBestRank("Square_BestPath.out");

    pop.Sort();
    pop.PrintDistance(0, "Distance.r" + to_string(rank) + ".out"); //Print the best distance for each rank

    MPI_Finalize();

    return 0;

}
