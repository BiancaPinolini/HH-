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
#include "TLorentzVector.h"
#include "TString.h"
#include <algorithm>
#include <TRandom.h>

using namespace std;

void norm1	(TTree *signal, TTree *background, int nbin, float min, float max) {	

    TTreeReader reader_S (signal);
    TTreeReader reader_B (background);

    TTreeReaderValue<TLorentzVector> WW_S(reader_S, "WW_F");
    TTreeReaderValue<TLorentzVector> WW_B(reader_B, "WW_F");

    string title1 = "signal";
	string title2 = "background";

    TCanvas *can = new TCanvas("can", "Massa Invariante");
    THStack * dist = new THStack("dist", "Distribution histo");
	TH1D * hs = new TH1D( title1.c_str(), "Segnale", nbin, min, max);
 	TH1D * hb = new TH1D( title2.c_str(), "Background", nbin, min, max); 

	//SCRITTURA ISTOGRAMMA DISTRIBUZIONE
	double mWW_S = 0, mWW_B = 0;

	while (reader_S.Next()) {
        mWW_S = WW_S -> M();
		hs -> Fill(mWW_S);
	}
	
    while (reader_B.Next()) {
        mWW_B = WW_B -> M();
        hb -> Fill(mWW_B) ;
    }

    //NORMALIZZAZIONE ISTOGRAMMI A 1
	double int_s = hs -> Integral();
	hs -> Scale(1/int_s);
	cout << "int_s = " << hs -> Integral() << endl;
	double int_b = hb -> Integral();
	hb -> Scale(1/int_b);
	cout << "int_b = " << hb -> Integral() << endl;
	
    //ROBA DI GRAFICA
	hb->SetFillStyle(3003);
	hs->SetFillStyle(3004);
    hs->SetFillColor(2);
	hb->SetFillColor(4);
	hs->SetLineColor(2);
	hb->SetLineColor(4);

	//SOVRAPPONGO GRAFICI SEGNALE E FONDO
	dist->Add(hs);		
	dist->Add(hb);
    dist->Draw("histo nostack");

	gStyle->SetOptStat(0);	
	gPad->SetGrid(1,1);

	TLegend *legend1 = new TLegend(0.8,0.2,0.98,0.38);
    legend1->AddEntry(hs,"Signal", "f");
    legend1->AddEntry(hb,"Background", "f");
	legend1->Draw("SAME");

	dist -> SetTitle("Confronto tra mWW S e B; Mass[GeV/c^{2}]; Frequency"); 
	return;
}

