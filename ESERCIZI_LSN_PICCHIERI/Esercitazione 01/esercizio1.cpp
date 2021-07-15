#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

int main() {


       
	int M = 1E4; //Numero di lanci
	int N = 1E2; //Numero di blocchi
	int L = M/N;  //Numero di lanci per blocco 
	double sum = 0; //Somma del singolo blocco
	double sum2 = 0; //Somma al quadrato del singolo blocco
	double sumvar = 0; //Somma  del singolo blocco per la varianza
	double sumvar2 = 0; //Somma al quadrato del singolo blocco per la varianza
	double media = 0; //Media cumulativa
	double media2 = 0; //Somma cumulativa dei quadrati della media 
	double variance = 0 ; //Varianza cumulativa
	double variance2 = 0; //Somma cumulativa dei quadrati della varianza 
	double errore = 0; //Errore statistico
	double random = 0; //Numero random tra 0 e 1
	int bin = 0; //Individua bin per test chi quadro
	int Z = 1E2; //Numero di test chi quadro
	int Y = 1E4; //Numero numeri random per test chi quadro
	const int X = 1E2; //Numero di bin tra 0 e 1 per test chi quadro
	int chidist[X] = {0}; //Vettore per test chi quadro
	double chi2bin = 0; //Chi quadro del singolo bin
	double chi2 = 0; //Chi quadro 

   
	Random rnd;
   	int seed[4];
   	int p1, p2;
   	ifstream Primes("Primes");
   	if (Primes.is_open()){
      		Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
   	Primes.close();

   	ifstream input("seed.in");
   	string property;
   	if (input.is_open()){
      		while ( !input.eof() ){
         		input >> property;
         		if( property == "RANDOMSEED" ){
            			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            			rnd.SetRandom(seed,p1,p2);
         		}
      		}
      	input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	ofstream Media, Varianza, Chi2;

	Media.open("Media.dat");
	if (!Media.is_open()){
		cerr << "PROBLEM: Can't open Media.dat" << endl;
		return 1;
	}

	Varianza.open("Varianza.dat");
	if (!Varianza.is_open()){
		cerr << "PROBLEM: Can't open Varianza.dat" << endl;
		return 1;
	}

	Chi2.open("Chi2.dat");
	if (!Chi2.is_open()){
		cerr << "PROBLEM: Can't open Chi2.dat" << endl;
		return 1;
	}

	for (int i=0; i<N; i++) {

        	for(int j=0; j<L; j++){
                	random = rnd.Rannyu();
					sum +=random;						//Per calcolo media
					sumvar += pow(random-0.5,2);		//Per calcolo varianza

        	}
		
		//Calcolo media e errore

		sum /= L;		        
		sum2 = pow(sum,2);
		


		media += sum;
		media2 += sum2;
		errore = sqrt( ((media2 / (i+1)) - pow((media / (i+1)),2)) / (i+1));
		Media << i << " \t " << media / (i+1) << " \t " << errore << endl;
		
		//Calcolo varianza e errore

		sumvar /= L;		
		sumvar2 = pow(sumvar,2);
		variance += sumvar;
		variance2 += sumvar2;
		
		errore = sqrt( ((variance2 / (i+1)) - pow((variance / (i+1)),2)) / (i+1));

		Varianza << i << " \t " << variance / (i+1) << " \t " << errore << endl;

		//Calcolo chi quadro

		
		

		random = 0;
		sum = 0;
		sum2 = 0;
		sumvar = 0;
		sumvar2 = 0;
	}

	//Calcolo chi quadro

	for (int i=0; i<Z; i++){
		for (int j=0; j<Y; j++){

			random = rnd.Rannyu();
			random *= 100;						
			bin = trunc(random);
			chidist[bin]++;
		}
		for (int j=0; j<X; j++){

			chi2bin = (pow((chidist[j]-Y/Z),2)/(Y/X));
			chi2 += chi2bin;
			chidist[j] = 0;   //Rimetto a zero il vettore per il ciclo successivo
		}

		Chi2 << i << " \t " << chi2 << endl;
		chi2 = 0;
	}

	Media.close();
	Varianza.close();
	Chi2.close();

	return 0;
}
