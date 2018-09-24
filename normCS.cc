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

//CROSS SECTION
#define cr_s 0.01997
#define cr_b1 137.5
#define cr_b2 0.003464
#define cr_b3 0.002649


void normCS	(TTree *signal, TTree *bg1, TTree *bg2, TTree *bg3, int nbin, double min, double max) {

    TTreeReader reader_S (signal);
    TTreeReader reader_B1 (bg1);
	TTreeReader reader_B2 (bg2);
	TTreeReader reader_B3 (bg3);

    TTreeReaderValue<TLorentzVector> WW_S(reader_S, "WW_F");
    TTreeReaderValue<TLorentzVector> WW_B1(reader_B1, "WW_F");
    TTreeReaderValue<TLorentzVector> WW_B2(reader_B2, "WW_F");
    TTreeReaderValue<TLorentzVector> WW_B3(reader_B3, "WW_F");

    string title = "signal";
	string title1 = "background_1";
	string title2 = "background_2";
	string title3 = "background_3";

	TCanvas *can = new TCanvas("can", "Massa Invariante");
    THStack * dist = new THStack("dist", "Distribution histo");
	TH1D * hs = new TH1D (title.c_str(), "Segnale", nbin, min, max);
 	TH1D * hb1 = new TH1D (title1.c_str(), "Background_1", nbin, min, max); 
	TH1D * hb2 = new TH1D (title2.c_str(), "Background_2", nbin, min, max);
	TH1D * hb3 = new TH1D (title3.c_str(), "Background_3", nbin, min, max);

	//SCRITTURA ISTOGRAMMA DISTRIBUZIONE
	double mWW_S = 0, mWW_B1 = 0, mWW_B2 = 0, mWW_B3 = 0;

	while (reader_S.Next()) {
        mWW_S = WW_S -> M();
		hs -> Fill(mWW_S);
	}
	
    while (reader_B1.Next()) {
        mWW_B1 = WW_B1 -> M();
        hb1 -> Fill (mWW_B1);
    }

	while (reader_B2.Next()) {
		mWW_B2 = WW_B2 -> M();
		hb2 -> Fill(mWW_B2);
	}

	while (reader_B3.Next()) {
		mWW_B3 = WW_B3 -> M();
		hb3 -> Fill(mWW_B3);
	}

    //NORMALIZZAZIONE ISTOGRAMMI ALLA SEZIONE D'URTO
	double int_s = hs -> Integral();
	hs -> Scale(cr_s/int_s);
	cout << "int_s = " << hs -> Integral() << endl;

	double int_b1 = hb1 -> Integral();
	hb1 -> Scale(cr_b1/int_b1);
	cout << "int_b1 = " << hb1 -> Integral() << endl;

	double int_b2 = hb2 -> Integral();
	hb2 -> Scale(cr_b2/int_b2);
	cout << "int_b2 = " << hb2 -> Integral() << endl;

	double int_b3 = hb3 -> Integral();
	hb3 -> Scale(cr_b3/int_b3);
	cout << "int_b3 = " << hb3 -> Integral() << endl;

	//ROBA DI GRAFICA
	hs->SetFillStyle(3004);
	hb1->SetFillStyle(3003);
	hb2->SetFillStyle(3003);
	hb3->SetFillStyle(3003);
    hs->SetFillColor(2);
	hb1->SetFillColor(4);
	hb2->SetFillColor(4);
	hb3->SetFillColor(4);
	hs->SetLineColor(2);
	hb1->SetLineColor(4);
	hb2->SetLineColor(4);
	hb3->SetLineColor(4);

	dist->Add(hs);
	//dist->Add(hb1);
	dist->Add(hb2);
	dist->Add(hb3);

    dist->Draw("histo nostack");

	gStyle->SetOptStat(0);	
	gPad->SetGrid(1,1);

	TLegend *legend1 = new TLegend(0.8,0.2,0.98,0.38);
    legend1->AddEntry(hs,"Signal", "f");
    legend1->AddEntry(hb2,"Background", "f");
	legend1->Draw("same");

	dist -> SetTitle("Confronto tra mWW S e B; Mass[GeV/c^{2}]; Frequency");

	return;
}

