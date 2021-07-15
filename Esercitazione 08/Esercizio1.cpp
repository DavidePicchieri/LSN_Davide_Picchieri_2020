#include <cmath>
#include <iostream>
#include <fstream>
#include "random.h"
#include "metropolis.h"
#include "PDF.h"
#include "blocking.h"

using namespace std;

int main(){
    double sigma = 0.8; //Sigma for wavefunction
    double mu = 0.6; //Mu for wavefunction
    int N_eval = 1E6; //Number of samples
    int N_block = 1E3; //Number of blocks
    int N_equilibration = 1E2; //Number of steps to equilibrate the system

    Random rnd("Primes","seed.in");
    Random * p_rnd = & rnd;
    WaveFunction pdf(sigma, mu);
    PDF *p_pdf = &pdf;

    Metropolis WFSampling(20, p_rnd, p_pdf);
 
    Blocking WFBlocks(N_eval, N_block);

    ofstream WFEvalue, WFPosition, WFParametersOptimization, WFEvalueBest, WFhistostream;
    
    WFEvalue.open("ExpectationValue.dat");
	if (!WFEvalue.is_open()){
		cerr << "PROBLEM: Can't open ExpectationValue.dat" << endl;
		return 1;
	}
    WFPosition.open("Xcoord.dat");
	if (!WFPosition.is_open()){
		cerr << "PROBLEM: Can't open Xcoord.dat" << endl;
		return 1;
	}
    WFParametersOptimization.open("ValuesParameteres.dat");
	if (!WFParametersOptimization.is_open()){
		cerr << "PROBLEM: Can't open ValuesParameteres.dat" << endl;
		return 1;
	}

    //Creating a method to evaluate <H> after the equilibration [8.01] 

    for (int i=0; i<N_equilibration; i++){
        if (i < 100){
            WFPosition << WFSampling.GetX() << endl;
        }
        WFSampling.Passo1DUniforme(3.3);
    }



    for (int i=0; i<N_eval; i++){
        WFSampling.Passo1DUniforme(3.3);
        WFBlocks.Add(WFSampling.GetEloc(), WFEvalue); 
    }

    //Now I try to minimize the parameters mu and sigma [8.02]

    double step_mu = 0.05, range_mu = 2, step_sigma = 0.05, range_sigma = 1.5, mu_best = mu, sigma_best = sigma;
    double LowestEnergy = WFBlocks.GetAverage();
    int sign = 1;
    mu = 0;
    sigma = 0.05;
 
    WFBlocks.Reset();

    //Scanning on the sigma-mu plane to find the minimum
    while (mu < range_mu){
        p_pdf->SetSigma(sigma);
        p_pdf->SetMu(mu);      
        while (sigma <= range_sigma && sigma >= step_sigma){   
            for (int i=0; i<N_eval; i++){
                WFSampling.Passo1DUniforme(3.3);
                WFBlocks.Add(WFSampling.GetEloc());
            }
            WFParametersOptimization << sigma << "\t" << mu << "\t" << WFBlocks.GetAverage() << endl;
            if (WFBlocks.GetAverage() < LowestEnergy){
                sigma_best = sigma;
                mu_best = mu;
                LowestEnergy = WFBlocks.GetAverage();
            }
            WFBlocks.Reset();
            //WFSampling.Reset();
            sigma += step_sigma*sign;
            cout << "sigma = " << sigma << "\t mu = " << mu << endl;
            p_pdf->SetSigma(sigma);
        }
        sign *= -1;
        mu += step_mu;
        sigma += step_sigma*sign;
        p_pdf->SetMu(mu);      
        p_pdf->SetSigma(sigma);
        for (int i=0; i<N_eval; i++){
            WFSampling.Passo1DUniforme(3.3);
            WFBlocks.Add(WFSampling.GetEloc());
        }
        WFParametersOptimization << sigma << "\t" << mu << "\t" << WFBlocks.GetAverage() << endl;
        if (WFBlocks.GetAverage() < LowestEnergy){
            sigma_best = sigma;
            mu_best = mu;
            LowestEnergy = WFBlocks.GetAverage();
        }
        WFBlocks.Reset();
        //WFSampling.Reset();
        cout << "sigma = " << sigma << "\t mu = " << mu << endl;
    }

    mu = mu_best - step_mu;
    double mu_max = mu_best + step_mu;
    double sigma_min = sigma_best - step_sigma;
    sigma = sigma_min;
    double sigma_max = sigma_best + step_sigma;
    step_mu /= 10;
    step_sigma /= 10;
    double xmin = -2, xmax = 2, deltax = 0.1;
    int nbins = (xmax - xmin) / deltax;
    int WFhisto[nbins];

    for (int i=0; i<nbins; i++){
        WFhisto[i] = 0;
    }

    //scanning again but on a smaller "area" and with smaller steps

    while (mu < mu_max){
        p_pdf->SetSigma(sigma);
        p_pdf->SetMu(mu);      
        while (sigma <= sigma_max && sigma >= sigma_min){   
            for (int i=0; i<N_eval; i++){
                WFSampling.Passo1DUniforme(3.3);
                WFBlocks.Add(WFSampling.GetEloc());
            }
            WFParametersOptimization << sigma << "\t" << mu << "\t" << WFBlocks.GetAverage() << endl;
            if (WFBlocks.GetAverage() < LowestEnergy){
                sigma_best = sigma;
                mu_best = mu;
                LowestEnergy = WFBlocks.GetAverage();
            }
            WFBlocks.Reset();
            //WFSampling.Reset();
            sigma += step_sigma*sign;
            cout << "sigma = " << sigma << "\t mu = " << mu << endl;
            p_pdf->SetSigma(sigma);
        }
        sign *= -1;
        mu += step_mu;
        sigma += step_sigma*sign;
        p_pdf->SetMu(mu);      
        p_pdf->SetSigma(sigma);
        for (int i=0; i<N_eval; i++){
            WFSampling.Passo1DUniforme(3.3);
            WFBlocks.Add(WFSampling.GetEloc());
        }
        WFParametersOptimization << sigma << "\t" << mu << "\t" << WFBlocks.GetAverage() << endl;
        if (WFBlocks.GetAverage() < LowestEnergy){
            sigma_best = sigma;
            mu_best = mu;
            LowestEnergy = WFBlocks.GetAverage();
        }
        WFBlocks.Reset();
        //WFSampling.Reset();
        cout << "sigma = " << sigma << "\t mu = " << mu << endl;
    }
    
    cout << WFSampling.GetAlphaSum() / double(N_eval) << endl;
    cout << "sigma_best = " << sigma_best << "\t mu_best = " << mu_best << "\t Lowest Energy = "<< LowestEnergy <<  endl;

    p_pdf->SetMu(mu_best);      
    p_pdf->SetSigma(sigma_best);
    
    WFEvalueBest.open("BestExpectationValue.dat");
	if (!WFEvalueBest.is_open()){
		cerr << "PROBLEM: Can't open BestExpectationValue.dat" << endl;
		return 1;
	}

    WFhistostream.open("Histogram.dat");
	if (!WFhistostream.is_open()){
		cerr << "PROBLEM: Can't open Histogram.dat" << endl;
		return 1;
	}


    for (int i=0; i<N_equilibration; i++){
        WFSampling.Passo1DUniforme(3.3);
    }

    for (int i=0; i<N_eval; i++){
        WFSampling.Passo1DUniforme(3.3);
        WFBlocks.Add(WFSampling.GetEloc(), WFEvalueBest);
        WFhistostream << WFSampling.GetX() << endl;
    }

    WFEvalue.close();
    WFPosition.close();
    WFParametersOptimization.close();
    WFEvalueBest.close();
    WFhistostream.close();

    return 0;
}