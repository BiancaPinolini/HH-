#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <THStack.h>
#include <TApplication.h>
#include <fstream>
#include "delta.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <algorithm>
#include <TRandom.h>

using namespace std;

void DrawHisto	(TTree *tree){	
	
	TTreeReader reader (tree);

	TTreeReaderValue<TLorentzVector> b1_I(reader, "b1_I"); 
	TTreeReaderValue<TLorentzVector> b1_F(reader, "b1_F");

	int nbin = 100;
	double min = -0.5, max = 0.5;
	double prima = 0, dopo = 0, delta = 0;
	string title = "delta";

	TCanvas *can = new TCanvas("can", "Confronto tra b1_Pt");
 	TH1D *h = new TH1D( title.c_str(), "Confronto tra b1_Pt", nbin, min, max); 
	
	//SCRITTURA ISTOGRAMMA DISTRIBUZIONE
	while ( reader.Next() ) {
		prima = b1_I -> Pt();
		dopo = b1_F -> Pt();
		delta = (dopo - prima) / prima;
		h->Fill( delta );
	}
	
	gStyle->SetOptStat(0);
	h -> Draw();
	TF1 *fit = new TF1("gaus", "gaus(0)", -2., 12.);
	fit -> FixParameter(1, 0.);
	fit -> FixParameter(2, 0.1);

	h -> Fit(fit);

	double mean, de_mean; 
	mean = fit -> GetParameter(1);
	de_mean = fit -> GetParError(1);

	cout << " mean = " << mean << " pm " << de_mean << endl;

    double sigma, de_sigma, width, de_width;
	sigma = fit -> GetParameter(2);
	de_sigma = fit -> GetParError(2);

    width = 2.355 * (sigma);
	de_width = 2.355 * de_sigma;
	cout << "width = " << width << " pm " << de_width << endl;

	return;
}

