#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"
#include "functions.h"
#include "metropolis.h"
#include "PDF.h"
#include "blocking.h"

using namespace std;

int main(){
    int N_eval = 1E6; //Numero di samples per valutare il raggio medio
    int N_block = 1E2; //Numero di blocchi in cui dividere i dati
    int N_equilibrium = 1E3; //Numero di passi per arrivare all'equilibrio prima di iniziare il sampling

    Random rnd("Primes","seed.in");
    Random * p_rnd = & rnd;
    PDF * p_pdf1 = new Hydro100();
    PDF * p_pdf2 = new Hydro210();


    Metropolis met100U(50, 50, 50, p_rnd, p_pdf1); //Creo le 4 classi metropoli 2 funzioni d'onda x 2 sampling
    Metropolis met100G(50, 50, 50, p_rnd, p_pdf1);
    Metropolis met210U(50, 50, 50, p_rnd, p_pdf2);
    Metropolis met210G(50, 50, 50, p_rnd, p_pdf2);

    ofstream R_100_Unif, R_100_Gauss, R_210_Unif, R_210_Gauss, Coords100, Coords210;
    
    R_100_Unif.open("R_100_Unif.dat");
	if (!R_100_Unif.is_open()){
		cerr << "PROBLEM: Can't open R_100_Unif.dat" << endl;
		return 1;
	}
    R_100_Gauss.open("R_100_Gauss.dat");
	if (!R_100_Gauss.is_open()){
		cerr << "PROBLEM: Can't open R_100_Gauss.dat" << endl;
		return 1;
	}
    R_210_Unif.open("R_210_Unif.dat");
	if (!R_210_Unif.is_open()){
		cerr << "PROBLEM: Can't open R_210_Unif.dat" << endl;
		return 1;
	}
    R_210_Gauss.open("R_210_Gauss.dat");
	if (!R_210_Gauss.is_open()){
		cerr << "PROBLEM: Can't open R_210_Gauss.dat" << endl;
		return 1;
	}

    Coords100.open("Coords100.dat");
	if (!Coords100.is_open()){
		cerr << "PROBLEM: Can't open Coords100.dat" << endl;
		return 1;
	}

    Coords210.open("Coords210.dat");
	if (!Coords210.is_open()){
		cerr << "PROBLEM: Can't open Coords210.dat" << endl;
		return 1;
	}


    Blocking Unif100(N_eval, N_block); //Class Blocking to perform the blocking average
    Blocking Gauss100(N_eval, N_block);
    Blocking Unif210(N_eval, N_block);
    Blocking Gauss210(N_eval, N_block);

    for (int i=0; i<N_equilibrium; i++){ //Equilibration
        met100U.PassoCartUniforme(1.2);
        met100G.PassoCartGauss(0.75);
        met210U.PassoCartUniforme(3);
        met210G.PassoCartGauss(1.8);    

        //Saving the coords of the points on file for the two uniform distributions in order to plot them

        Coords100 << met100U.GetX() << "\t" << met100U.GetY() << "\t" << met100U.GetZ() << endl;
        Coords210 << met210U.GetX() << "\t" << met210U.GetY() << "\t" << met210U.GetZ() << endl;

    }

    for (int i=0; i<N_eval; i++){  //Actual simulation
        met100U.PassoCartUniforme(1.2);
        met100G.PassoCartGauss(0.75);
        met210U.PassoCartUniforme(3);
        met210G.PassoCartGauss(1.8);    

        Unif100.Add(met100U.GetR(), R_100_Unif);
        Gauss100.Add(met100G.GetR(), R_100_Gauss);
        Unif210.Add(met210U.GetR(), R_210_Unif);
        Gauss210.Add(met210G.GetR(), R_210_Gauss);

        //Saving the coords of the points on file for the two uniform distributions in order to plot them

        if (i % 100 == 0){
            Coords100 << met100U.GetX() << "\t" << met100U.GetY() << "\t" << met100U.GetZ() << endl;
            Coords210 << met210U.GetX() << "\t" << met210U.GetY() << "\t" << met210U.GetZ() << endl;
        }
    }

    cout << met100U.GetAlphaSum() / double (N_eval) << "\t" << met100G.GetAlphaSum() / double (N_eval) << "\t"<< met210U.GetAlphaSum() / double (N_eval) << "\t" << met210G.GetAlphaSum() / double (N_eval) << endl; 

    R_100_Unif.close();
    R_100_Gauss.close();
    R_210_Unif.close();
    R_210_Gauss.close();
    Coords100.close();
    Coords210.close();


    return 0;

}