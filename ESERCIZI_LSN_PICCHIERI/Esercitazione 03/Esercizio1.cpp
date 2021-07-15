#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"
#include "functions.h"

using namespace std;

int main(){
    int N_eval = 1E6; //Numero di samples per valutare il prezzo medio
    int N_block = 1E2; //Numero di blocchi in cui dividere i dati
    int L = N_eval / N_block; //Numero di samples per blocco
    double S0 = 100; //Prezzo al tempo t = 0
    double T = 1; //Tempo della consegna
    double N_int = 1E2; //Numero di intervalli di tempo 
    double dt = T / N_int; //Lasso di tempo tra due momenti in cui viene valutato St
    double St = 0; //Prezzo al tempo t
    double K = 100; //Strike price
    double r = 0.1; //Tasso d'interesse a rischio zero
    double sigma = 0.25; //Volatilità
    double Wt = 0; // Valore Wiener process al tempo t
    double Zi = 0; //Variabile che verrà generata casualmente tra 0 e 1

    //Variabili usate per costo call option

    double C_dir = 0; //Stima costo della call-option con sample diretto di S(T)
    double C_disc = 0; //Stima costo della call-option con sample discreto dei vari S(t_i)
    double C_dir_ave = 0; //Media delle stime costo della call-option con sample diretto di S(T)
    double C_dir_ave2 = 0; //Media dei quadrati delle stime costo della call-option con sample diretto di S(T)
    double C_disc_ave = 0; //Media delle stime costo della call-option con sample discreto dei vari S(t_i)
    double C_disc_ave2 = 0; //Media dei quadrati delle stime costo della call-option con sample discreto dei vari S(t_i)
    double error_C_dir = 0; //Stima incertezza sample diretto
    double error_C_disc = 0; //Stima incertezza sample discreto

    //Variabili usate per costo put option

    double P_dir = 0; //Stima costo della put-option con sample diretto di S(T)
    double P_disc = 0; //Stima costo della put-option con sample discreto dei vari S(t_i)
    double P_dir_ave = 0; //Media delle stime costo della put-option con sample diretto di S(T)
    double P_dir_ave2 = 0; //Media dei quadrati delle stime costo della put-option con sample diretto di S(T)
    double P_disc_ave = 0; //Media delle stime costo della put-option con sample discreto dei vari S(t_i)
    double P_disc_ave2 = 0; //Media dei quadrati delle stime costo della put-option con sample discreto dei vari S(t_i)
    double error_P_dir = 0; //Stima incertezza sample diretto put-option
    double error_P_disc = 0; //Stima incertezza sample discreto put-option



    Random rnd("Primes","seed.in");

    ofstream Call_dir, Call_disc, Put_dir, Put_disc;
    
    Call_dir.open("Call_dir.dat");
	if (!Call_dir.is_open()){
		cerr << "PROBLEM: Can't open Call_dir.dat" << endl;
		return 1;
	}
    Call_disc.open("Call_disc.dat");
	if (!Call_disc.is_open()){
		cerr << "PROBLEM: Can't open Call_disc.dat" << endl;
		return 1;
	}
    Put_dir.open("Put_dir.dat");
	if (!Put_dir.is_open()){
		cerr << "PROBLEM: Can't open Put_dir.dat" << endl;
		return 1;
	}
    Put_disc.open("Put_disc.dat");
	if (!Put_disc.is_open()){
		cerr << "PROBLEM: Can't open Put_disc.dat" << endl;
		return 1;
	}


    for (int i=0; i<N_block; i++){
        for (int j=0; j<L; j++){

            //Calcolo direttamente S(T)
     
            Wt = rnd.Gauss(0.0,T);
            St = Price(T, S0, r, sigma, Wt); //Prezzo calcolato direttamente al momento T

            C_dir += exp(-r*T) * max(0.,St-K);
            P_dir += exp(-r*T) * max(0.,K-St);
            
            //Calcolo S(t_i) per 100 passi

            St = S0; //Ri-inizializzo il prezzo

            for (int k=0; k<N_int; k++){
                Zi = rnd.Gauss(0.0,1);
                St = Price(dt, St, r, sigma, Zi*sqrt(dt));
            }            

            C_disc += exp(-r*T) * max(0.,St-K);
            P_disc += exp(-r*T) * max(0.,K-St);
        }

        //calcolo valori medi dopo ogni blocco per usare metodo dei blocchi e salvo su file

        C_dir /= L;
        C_disc /= L;
        P_dir /= L;
        P_disc /= L;

        C_dir_ave += C_dir;
        C_dir_ave2 += pow(C_dir,2);
        P_dir_ave += P_dir;
        P_dir_ave2 += pow(P_dir,2);

        error_C_dir = Error(C_dir_ave, C_dir_ave2, i); 
        error_P_dir = Error(P_dir_ave, P_dir_ave2, i); 

        Call_dir << i << " \t " << C_dir_ave / (i + 1) << " \t " << error_C_dir << endl;
        Put_dir << i << " \t " << P_dir_ave / (i + 1) << " \t " << error_P_dir << endl;

        C_disc_ave += C_disc;
        C_disc_ave2 += pow(C_disc,2);
        P_disc_ave += P_disc;
        P_disc_ave2 += pow(P_disc,2);

        error_C_disc = Error(C_disc_ave, C_disc_ave2, i); 
        error_P_disc = Error(P_disc_ave, P_disc_ave2, i); 

        Call_disc << i << " \t " << C_disc_ave / (i + 1) << " \t " << error_C_disc << endl;
        Put_disc << i << " \t " << P_disc_ave / (i + 1) << " \t " << error_P_disc << endl;

        //ri-inizializzo le variabili

        C_dir = 0;
        C_disc = 0;
        P_dir = 0;
        P_disc = 0;
    }

    Call_dir.close();
    Call_disc.close();
    Put_dir.close();
    Put_disc.close();

    return 0;
}