#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){

    int M = 1E4; //Numero realizzazioni totali
    int N[4] = {1,2,10,100}; //Numero di eventi su cui mediare per ogni realizzazione
    double sumunif = 0; //Risultato dado uniforme
    double sumexp = 0; //Risultato dado esponenziale
    double sumlor = 0; //Risultato dado lorentziano
    double lambda = 1; //Lambda per esponenziale
    double mu = 0; //Mu per lorentziana
    double gamma = 1; //Gamma per lorentziana

    ofstream Risultati;

    //Generatore di numeri casuali

    Random rnd("Primes", "seed.in");

	Risultati.open("Risultati.dat");
	if (!Risultati.is_open()){
		cerr << "PROBLEM: Can't open Risultati.dat" << endl;
		return 1;
	}

    for (int i=0; i<4; i++){
        
        for (int j=0; j<M; j++){
            for (int k=0; k<N[i]; k++){
                sumunif += rnd.Rannyu(0.0,6.0);
                sumexp += rnd.Exponential(lambda);
                sumlor += rnd.Lorentz(mu, gamma);
            }

            Risultati << " \t " << sumunif / N[i] << " \t " << sumexp / N[i] << " \t " << sumlor / N[i] << endl;

            sumunif = 0;
            sumexp = 0;
            sumlor = 0;

        }


    }

    Risultati.close();



    return 0;

}