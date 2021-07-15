#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"
#include "functions.h"

using namespace std;

int main(){

    Random rnd("Primes", "seed.in");
    int Neval = 1E6; //Numero di volte in cui valuto l'integrale
    int Nblocks = 1E2; // Numero di blocchi
    int L = Neval / Nblocks; //Numero di valutazioni per blocco
    double xmin = 0; //Valore minimo di x nell'intervallo
    double xmax =1; //Valore massimo di x nell'intervallo
    double x = 0; //Valore di x in cui valuterò la funzione
    double f = 0; //Dichiaro la mia funzione
    double Integral = 0; //Stima progressiva dell'integrale
    double Integral2 = 0; //Somma delle stime progressive dell'integrale al quadrato
    double error = 0; //Errore trovato tramite metodo dei blocchi
    double y = 0; //Valore di y in cui valuterò la funzione usata per l'importance sampling
    double g = 0; //Dichiaro la mia nuova funzione dell'importance sampling 
    double p = 0; //Nuova distribuzione per l'importance sampling
    double Integral_imp = 0; //Stima progressiva dell'integrale con importance sampling
    double Integral2_imp = 0; //Somma delle stime progressive dell'integrale al quadrato con importance sampling
    double error_imp = 0; //Errore trovato tramite metodo dei blocchi con importance sampling

    ofstream Int_unif, Int_imp;
    
    Int_unif.open("Int_unif.dat");
	if (!Int_unif.is_open()){
		cerr << "PROBLEM: Can't open Int_unif.dat" << endl;
		return 1;
	}
    Int_imp.open("Int_imp.dat");
	if (!Int_imp.is_open()){
		cerr << "PROBLEM: Can't open Int_imp.dat" << endl;
		return 1;
	}

    for (int i=0; i<Nblocks; i++){
        for (int j=0; j<L; j++){
            x = rnd.Rannyu(xmin, xmax);
            f += (M_PI / 2.) * cos( M_PI * x / 2.);

            y = rnd.Rannyu();
            p = 1-sqrt(y);
            g += (M_PI / 2.) * cos( M_PI * p / 2.) / (2 *(1 - p));

        }

        Integral += f / L;
        Integral2 += pow(f / L, 2); 
		error = sqrt( ((Integral2 / (i+1)) - pow((Integral / (i+1)),2)) / (i+1));

        Int_unif << i << " \t " << Integral / (i+1) << " \t " << error << endl;

        Integral_imp += g / L;
        Integral2_imp += pow(g / L, 2); 
		error_imp = sqrt( ((Integral2_imp / (i+1)) - pow((Integral_imp / (i+1)),2)) / (i+1));

        Int_imp << i << " \t " << Integral_imp / (i+1) << " \t " << error_imp << endl;

        f = 0;
        g = 0;
    }


    Int_unif.close();
    Int_imp.close();

    return 0;
}