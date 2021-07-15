#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "random.h"
#include "functions.h"
#include "posizione.h"

using namespace std;

int main(){
    int Nwalks = 1E5; //Numero di walk diverse
    int Nblocks = 1E2; //Numero di blocchi
    int L = Nwalks / Nblocks; //Numero di walks per blocco
    int Nsteps = 1E2; //Numero di step della singola walk
    double a = 1; //Lunghezza passo del walk
    double direzione = 0; //Per decidere in che direzione si muove nel lattice
    double phi = 0; //Per decidere in che direzione si muove nel continuo
    double theta = 0; //Per decidere in che direzione si muove nel continuo
    vector <double> appo_lat (100,0); //Vettore da usare per push_back lattice
    vector <double> appo_con (100,0); //Vettore da usare per push_back continuo
    vector <double> appo2_lat (100,0); //Vettore da usare per push_back lattice con quad di r quadro
    vector <double> appo2_con (100,0); //Vettore da usare per push_back continuo con quad di r quadro
    vector <double> error_lat (100,0); // Errore nel calcolo del r quadro lattice
    vector <double> error_con (100,0); // Errore nel calcolo del r quadro lattice
    vector <double> r2_lat (100,0); //Salvo i risultati di r2 dei blocchi per il lattice
    vector <double> r2_con (100,0); //Salvo i risultati di r2 dei blocchi per il continuo
    vector <double> r2_lat_2 (100,0); //Salvo i risultati di r2 al quadrato dei blocchi per il lattice
    vector <double> r2_con_2 (100,0); //Salvo i risultati di r2 al quadrato dei blocchi per il continuo



    Posizione p_lat(0., 0., 0.), p_con(0., 0., 0.); //Classe in cui salvo la posizione
    Random rnd("Primes", "seed.in");


    ofstream Lattice, Continuum;

    Lattice.open("Lattice.dat");
	if (!Lattice.is_open()){
		cerr << "PROBLEM: Can't open Lattice.dat" << endl;
		return 1;
	}
    Continuum.open("Continuum.dat");
	if (!Lattice.is_open()){
		cerr << "PROBLEM: Can't open Continuum.dat" << endl;
		return 1;
	}
    for (int i=0; i< Nblocks; i++){
        for (int j=0; j<L; j++){
            for (int k=0; k<Nsteps; k++){
                direzione = rnd.Rannyu(0.,6.);
                phi = 2 * M_PI * rnd.Rannyu();
                theta = acos(1- 2 * rnd.Rannyu());

                if (direzione < 1)
                    p_lat.PassoCart(a,0.,0.);

                else if (direzione < 2)
                    p_lat.PassoCart(-a,0.,0.);

                else if (direzione < 3)
                    p_lat.PassoCart(0.,a,0.);

                else if (direzione < 4)
                    p_lat.PassoCart(0.,-a,0.);

                else if (direzione < 5)
                    p_lat.PassoCart(0.,0.,a);

                else
                    p_lat.PassoCart(0.,0.,-a);

                p_con.PassoSfer(a, phi, theta);

                appo_lat[k] += p_lat.DistanzaQuadDaPunto(0.,0.,0.);
                appo_con[k] += p_con.DistanzaQuadDaPunto(0.,0.,0.);

                appo2_lat[k] += pow(p_lat.DistanzaQuadDaPunto(0.,0.,0.), 2);
                appo2_con[k] += pow(p_con.DistanzaQuadDaPunto(0.,0.,0.), 2);

            }

            p_lat.SetPosizione(0., 0., 0.);
            p_con.SetPosizione(0., 0., 0.);



        }                     

        for (int k=0; k<Nsteps; k++){

            r2_lat[k] += appo_lat[k] / L;
            r2_con[k] += appo_con[k] / L;

            r2_lat_2[k] += appo2_lat[k] / L;
            r2_con_2[k] += appo2_con[k] / L;


            error_lat[k] = sqrt( ((r2_lat_2[k] / (i + 1)) - pow((r2_lat[k] / (i + 1)),2)) / (i + 1));
		    error_con[k] = sqrt( ((r2_con_2[k] / (i + 1)) - pow((r2_lat[k] / (i + 1)),2)) / (i + 1));


            appo_lat[k] = 0;
            appo_con[k] = 0;

            appo2_lat[k] = 0;
            appo2_con[k] = 0;


        }



    }

    for (int i=0; i<Nsteps; i++){



        Lattice << i << " \t " << sqrt (r2_lat[i] / Nblocks) <<  " \t " << error_lat[i] * ( 1./ (2 * sqrt (r2_lat[i] / Nblocks))) << endl;
        Continuum << i << " \t " << sqrt (r2_con[i] / Nblocks) <<  " \t " << error_con[i] * ( 1./ (2 * sqrt (r2_con[i] / Nblocks))) << endl;


    }

    Lattice.close();
    Continuum.close();

    return 0;
}
