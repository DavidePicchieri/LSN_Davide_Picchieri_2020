#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "functions.h"
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){

    int Nthrows = 1E5; //Numero lanci
    int Nhits = 0; //Numero di sbarrette che intersecano la linea
    int Nblocks = 1E2; //Numero di blocchi
    int Nthr_blk = Nthrows / Nblocks; //Numero di lanci per blocco
    double D = 1; //Distanza tra le linee
    double L = 0.7; //Lunghezza delle sbarrette
    double ymin = 0; //Valore minimo generazione y casuale
    double ymax = 1E4; //Valore massimo generazione y casuale
    double ysbarretta1 = 0; //Coord y del punto 'iniziale' della sbarretta
    double x_ran = 0; //Valore x secondo punto generato uniforme in circonferenza di raggio L
    double y_ran = 0; //Valore y secondo punto generato uniforme in circonferenza di raggio L
    double ysbarretta2 = 0; //Coord y del punto 'finale' della sbarretta
    bool halfplane = true; // Check su quale metà del piano si trova il secondo punto
    double m = 0; //Coefficiente angolare retta congiungente tra punto 1 e punto 2 generato uniformemente
    int parte_intera_y1 = 0; 
    int parte_intera_y2 = 0; 
    double av_pi = 0; //Stima progressiva di pi
    double av2_pi = 0; //Stima al quadrato di pi progressiva
    double pi_block = 0; //Stima di pi nel singolo blocco
    double error = 0; //Errore calcolato con metodo dei blocchi

    Random rnd("Primes", "seed.in");

    ofstream pi;
    pi.open("pi.dat");

	if (!pi.is_open()){
		cerr << "PROBLEM: Can't open pi.dat" << endl;
		return 1;
	}

    for (int i=0; i<Nblocks; i++){

        for (int j=0; j<Nthr_blk; j++){

            ysbarretta1 = rnd.Rannyu(ymin, ymax);

            do{
                x_ran = rnd.Rannyu(-0.7, 0.7);
                y_ran = rnd.Rannyu(-0.7, 0.7);
            }while(x_ran*x_ran + y_ran*y_ran > L*L);

            m = Slope(0., ysbarretta1, x_ran, ysbarretta1 + y_ran);

            if (y_ran < 0) {halfplane = false;}
            ysbarretta2 = ysbarretta1 + Y_Intersection(m, halfplane, L);

            if (i == 0){
                cout << y_ran << " " << m << " " << Y_Intersection(m, halfplane, L) << endl;
            }

            parte_intera_y1 = trunc(ysbarretta1);
            parte_intera_y2 = trunc(ysbarretta2);
            if (parte_intera_y1 != parte_intera_y2){Nhits++;}
        }

        pi_block = 2. * L * Nthr_blk / ( Nhits * D);
        //cout << pi_block << " " << Nhits << endl;
        av_pi += pi_block;
        av2_pi += pow(pi_block,2);
        error = sqrt( ((av2_pi / (i+1)) - pow((av_pi / (i+1)),2)) / (i+1));

        pi << i << " \t " << av_pi / (i+1)<< " \t " << error << endl;

        Nhits = 0;
        halfplane = true;
    }

    pi.close();



    return 0;
}