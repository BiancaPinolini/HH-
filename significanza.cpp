#include <iostream>
#include <TH1F.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TStyle.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TApplication.h>
#include <fstream>
#include "TLorentzVector.h"

using namespace std;

//CROSS SECTION
#define cr_s 0.01997
#define cr_b1 137.5
#define cr_b2 0.003464
#define cr_b3 0.002649
//NUMERO INIZIALE DI EVENTI
#define N_s 7709
#define N_b1 337452
#define N_b2 299199
#define N_b3 130574
//EFFICIENZE
#define e_s 0.08
#define e_b1 0.34
#define e_b2 0.30
#define e_b3 0.13
//LUMINOSITA'
#define lumi 1500

int main(){
	TApplication * Grafica = new TApplication("App", 0, 0);

	//LETTURA DEI TTree
    TFile * sinput = TFile::Open ("HH_cut.root");
    TFile * b1input = TFile::Open ("ttbar_cut.root");
    TFile * b2input = TFile::Open ("WWZ_cut.root");
    TFile * b3input = TFile::Open ("WZZ_cut.root");
    
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * bg1 = (TTree*)b1input->Get("tree");
    TTree * bg2 = (TTree*)b2input->Get("tree");
    TTree * bg3 = (TTree*)b3input->Get("tree");
  
	int nbin = 50;
	//INSERISCO MINIMO E MASSIMO DELL'ISTOGRAMMA
	double min = 20, max = 300;

	TTreeReader reader_S (signal);
	TTreeReader reader_B1 (bg1);
	TTreeReader reader_B2 (bg2);
	TTreeReader reader_B3 (bg3);

	TTreeReaderValue<TLorentzVector> var_S(reader_S, "MET_F");
	TTreeReaderValue<TLorentzVector> var_B1(reader_B1, "MET_F");
	TTreeReaderValue<TLorentzVector> var_B2(reader_B2, "MET_F");
	TTreeReaderValue<TLorentzVector> var_B3(reader_B3, "MET_F");

	// TTreeReaderValue<float> var_S(reader_S, "bb_Phi");
	// TTreeReaderValue<float> var_B1(reader_B1, "bb_Phi");
	// TTreeReaderValue<float> var_B2(reader_B2, "bb_Phi");
	// TTreeReaderValue<float> var_B3(reader_B3, "bb_Phi");

	string title = "signal";
	string title1 = "ttbar";
	string title2 = "WWZ";
	string title3 = "WZZ";

	TCanvas *can = new TCanvas("can", "Massa Invariante");
	THStack * TH_bg = new THStack("TH_bg", "Background histo");
	THStack * TH_s = new THStack("TH_s", "Signal histo");
	TH1F * hs = new TH1F (title.c_str(), "Segnale", nbin, min, max);
	TH1F * hb1 = new TH1F (title1.c_str(), "ttbar", nbin, min, max); 
	TH1F * hb2 = new TH1F (title2.c_str(), "WWZ", nbin, min, max);
	TH1F * hb3 = new TH1F (title3.c_str(), "WZZ", nbin, min, max);

	//SCRITTURA ISTOGRAMMA DI DIRIBUZIONE
	float var = 0;
	float s = 0, b1 = 0, b2 = 0, b3 = 0;

	while (reader_S.Next()) {
		var =  var_S->Pt();
		// var = *var_S;
		hs -> Fill(var);
		s++;
	}

	while (reader_B1.Next()){
		var =  var_B1->Pt();
		// var = *var_B1;
		hb1 -> Fill (var);
		b1++;
	}

	while (reader_B2.Next()){
		var =  var_B2->Pt();
		// var = *var_B2;
		hb2 -> Fill(var);
		b2++;
	}

	while (reader_B3.Next()){
		var =  var_B3->Pt();
		// var = *var_B3;
		hb3 -> Fill(var);
		b3++;
	}

	//CALCOLO IL NUMERO DI EVENTI SOPRAVVISSUTI AL TAGLIO
	double int_s = hs -> Integral();
	double int_b1 = hb1 -> Integral();
	double int_b2 = hb2 -> Integral();
	double int_b3 = hb3 -> Integral();

	cout << "Efficienza di segnale = " << int_s / N_s << endl;
	cout << "Efficienza di fondo = " << (int_b1 + int_b2 + int_b3) / (N_b3 + N_b2 + N_b1)  << endl << endl;
/*
	//CALCOLO S e B COME EFFICIENZA * EFFICIENZA DEL PRIMO TAGLIO * CS * LUMINOSITA'
	float S = int_s * cr_s * lumi * e_s / N_s;
	float B = (int_b1 * cr_b1 * e_b1 / N_b1 + int_b2 * cr_b2 * e_b2 / N_b2 + int_b3 * cr_b3 * e_b3 / N_b3) * lumi;

	//CALCOLO SIGNIFICATIVITA'
	cout << "Q = " << S / sqrt(B) << endl;

	//NORMALIZZAZIONE ALLA CS
	//Moltiplico il segnale per un fattore 1000 altrimenti non si vede
	double norm_s = cr_s * lumi * e_s / N_s; 
	hs -> Scale (norm_s*1000);
	double norm_b1 = cr_b1 * lumi * e_b1 / N_b1;
	hb1 -> Scale (norm_b1);
	double norm_b2 = cr_b2 * lumi * e_b2 / N_b2;
	hb2 -> Scale (norm_b2);
	double norm_b3 = cr_b3 * lumi * e_b3 / N_b3;
	hb3 -> Scale (norm_b3);
*/
	//NORMALIZZAZIONE ISTOGRAMMI A 1
	hs -> Scale(1/s);
	hb1 -> Scale(1/b1);
	hb2 -> Scale(1/b2);
	hb3 -> Scale(1/b3);
	cout << "hs = " << hs->Integral() << endl;
	cout << "hb1 = " << hb1->Integral() << endl;

	//ROBA DI GRAFICA
	hs->SetLineColor(kAzure-6);	
	hs->SetLineWidth(3);

	hb1->SetFillStyle(3002);
	hb2->SetFillStyle(3002);
	hb3->SetFillStyle(3002);

	hb1->SetLineWidth(3);
	hb2->SetLineWidth(3);
	hb3->SetLineWidth(3);

	hb1->SetFillColor(kOrange-2);
	hb2->SetFillColor(kOrange+7);
	hb3->SetFillColor(kRed-3);

	hb1->SetLineColor(kOrange-2);
	hb2->SetLineColor(kOrange+7);
	hb3->SetLineColor(kRed-3);

	//STACK
	// TH_bg->Add(hb2);
	// TH_bg->Add(hb3);
	// TH_bg->Add(hb1);
	// TH_bg->Add(hs);

    // TH_bg->Draw("histo");

	//NoSTACK
	TH_bg->Add(hb2);
	TH_bg->Add(hb3);
	TH_bg->Add(hb1);
	TH_bg->Add(hs);

    TH_bg->Draw("histo nostack");

	//FONDO STACK ma SEGNALE NoSTACK
	// TH_s->Add(hs);
	// TH_bg->Add(hb2);
	// TH_bg->Add(hb3);
	// TH_bg->Add(hb1);

    // TH_bg->Draw("histo");
	// TH_s->Draw("histo same");

	gStyle->SetOptStat(0);	
	//gPad->SetLogy();

	TLegend *legend1 = new TLegend(0.8,0.2,0.98,0.38);
	legend1->AddEntry(hs,"Signal", "f");
	legend1->AddEntry(hb1,"ttbar", "f");
	legend1->AddEntry(hb2, "WWZ", "f");
	legend1->AddEntry(hb3, "WZZ", "f");
	legend1->Draw("same");

	TH_bg -> SetTitle("MET_{Pt}; Mass [GeV]; Frequency");
/*
	//TAGLIO E SIGNIFICANZA
	double s_sx = 0, s_dx = 0, b_sx = 0, b_dx = 0;
	double eff_s = 0, eff_b = 0;
	int bin = 0;
	double Ql = 0, Qr = 0, Ql_max = 0, Qr_max = 0, cut_l = 0, cut_r = 0;

	for(int j = 1; j <= nbin; j++){
		s_sx = hs->Integral(1,j);
		b_sx = hb1->Integral(1,j) + hb2->Integral(1,j) + hb3->Integral(1,j);

		s_dx = hs->Integral(j,nbin);
		b_dx = hb1->Integral(j,nbin) + hb2->Integral(j,nbin) + hb3->Integral(j,nbin);

		Ql = s_sx / sqrt(b_sx);
		Qr = s_dx / sqrt(b_dx);

		if (Ql > Ql_max) {
			Ql_max = Ql;
			cut_l = j;
		}

		if (Qr > Qr_max) {
			Qr_max = Qr;
			cut_r = j;
		}
	}
	
	double r = (cut_r * max) / 100, l = (cut_l * max) / 100;
	cout << "Taglio: " << r << " < var < " << l << endl << endl;

	double n_s = hs->Integral(cut_r,cut_l) * N_s;
	double n_b1 = hb1->Integral(cut_r, cut_l) * N_b1;
	double n_b2 = hb2->Integral(cut_r, cut_l) * N_b2;
	double n_b3 = hb3->Integral(cut_r, cut_l) * N_b3;

	cout << "Numero di eventi sopravvissuti al taglio: "  << endl;
	cout << "segnale: " << n_s << " su " << N_s << endl;
	cout << "fondo: " << n_b1 + n_b2 + n_b3 << " su " << N_b1 + N_b2 + N_b3 <<endl << endl; 

	cout << "Efficienza di segnale = " << n_s / N_s << endl;
	cout << "Efficienza di fondo = " << (n_b1 + n_b2 + n_b3) / (N_b3 + N_b2 + N_b1)  << endl << endl;

	float S = n_s * cr_s * lumi * e_s / N_s;
	float B = (n_b1 * cr_b1 * e_b1 / N_b1 + n_b2 * cr_b2 * e_b2 / N_b2 + n_b3 * cr_b3 *e_b3 / N_b3) * lumi;

	cout << "Q = " << S / sqrt(B) << endl;
*/
	Grafica->Run();
	return 0;
}
