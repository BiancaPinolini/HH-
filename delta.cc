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
	
	TTreeReader reader( tree );

	TTreeReaderValue<TLorentzVector> lep_I(reader, "lep_I"); 
	TTreeReaderValue<TLorentzVector> lep_F(reader, "lep_F");

	int nbin = 100;
	double min = -0.05, max = 0.05;
	double prima = 0, dopo = 0, delta = 0;
	string title = "delta";

	TCanvas *can = new TCanvas("can", "Confronto tra lep_Pt");
 	TH1D *h = new TH1D( title.c_str(), "Confronto tra lep_Pt", nbin, min, max); 
	
	//SCRITTURA ISTOGRAMMA DISTRIBUZIONE
	while ( reader.Next() ) {
		prima = lep_I -> Pt();
		dopo = lep_F -> Pt();
		delta = (dopo - prima) / prima;
		h->Fill( delta );
	}
	
	h -> Draw();
	TF1 *fit = new TF1("gaus", "gaus(0)");

	h -> Fit(fit);

	double mean, de_mean; 
	mean = fit -> GetParameter(1);
	de_mean = fit -> GetParError(1);

	if ((mean - de_mean) < 0 && (mean + de_mean) > 0) {
		cout << "La media è accettabile" << endl;
	} else {
		cout << "La media non è accettabile" << endl;
	}
	cout << " mean = " << mean << " pm " << de_mean << endl;

    double sigma, de_sigma, width, de_width;
	sigma = fit -> GetParameter(2);
	de_sigma = fit -> GetParError(2);

    width = 2.355 * (sigma);
	de_width = 2.355 * de_sigma;

	if ((width - de_width) < 0.01 && (width + de_width) > 0.01) {
		cout << "La width è accettabile" << endl;
	} else {
		cout << "La width non è accettabile" << endl;
	}
	cout << "width = " << width << " pm " << de_width << endl;

	return;
}

