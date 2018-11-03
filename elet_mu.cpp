//  c++ -o elet_mu elet_mu.cpp `root-config --cflags --glibs`
//  ./elet_mu HH.lhe

#include "LHEF.h"
#include<iostream>
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include <fstream>
#include <algorithm>
#include "TLorentzVector.h"
#include <TF1.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLegend.h>
#include <THStack.h>
#include <TApplication.h>
#include "TLorentzVector.h"

using namespace std;

int main(int argc, char** argv){

	TApplication * Grafica = new TApplication("App", 0, 0);

    TLorentzVector elet;
    TLorentzVector mu;

    int iEv = 0, n = 1;

    int ID = 0;

    ifstream ifs(argv[1]);
    LHEF::Reader reader(ifs);

    int nbin = 250;
	double min = -0.002, max = 0.002;
    double var = 0;

	TCanvas *c1 = new TCanvas("c1","multipads");
	gStyle->SetOptStat(0);
   	c1->Divide(1,2,0.001,0.001);
    TH1F *h_e = new TH1F("histo mass electrons", "e_{M}", nbin, min, max); 
    TH1F *h_m = new TH1F("histo mass muons", "#mu_{M}", nbin, min, max); 
    
    //Loop over all the events 
    while (reader.readEvent()){
        iEv++;

        vector<int> electrons;
        vector<int> muons;

        //per controllare che vada tutto bene stampo un #evento ogni 1000
        if (iEv%10000 == 0){
            cout << "Event " << iEv << endl;
        }

        //Working on info of final state particles only
        for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size(); iPart++){

            if (reader.hepeup.ISTUP.at (iPart) == 1) {
                ID = reader.hepeup.IDUP.at(iPart);

                //ELETTRONI
                if (abs(ID) == 11) { 
                    electrons.push_back(iPart);

                    elet.SetPxPyPzE(
                        reader.hepeup.PUP.at (iPart).at (0), //px
                        reader.hepeup.PUP.at (iPart).at (1), //py
                        reader.hepeup.PUP.at (iPart).at (2), //pz
                        reader.hepeup.PUP.at (iPart).at (3)  //E
                    );

                    var = elet.M();
                    h_e -> Fill(var);
                }
                
                //MUONI
                if (abs(ID) == 13 ){ 
                    muons.push_back(iPart);

                    mu.SetPxPyPzE(
                        reader.hepeup.PUP.at (iPart).at (0), //px
                        reader.hepeup.PUP.at (iPart).at (1), //py
                        reader.hepeup.PUP.at (iPart).at (2), //pz
                        reader.hepeup.PUP.at (iPart).at (3)  //E
                    );

                    var = mu.M();
                    h_m -> Fill(var);
                }
            }
        }        
    }
    double int_e = h_e->Integral();
	double int_m = h_m->Integral();

    //NORMALIZZAZIONE ISTOGRAMMI A 1
	h_e -> Scale(1/int_e);
	h_m -> Scale(1/int_m);

    //ROBA DI GRAFICA
	h_e->SetLineColor(kAzure-6);	
	h_e->SetLineWidth(3);

	h_m->SetLineColor(kAzure-6);	
	h_m->SetLineWidth(3);

    c1->cd(1);
    h_m->Draw("histo");

    c1->cd(2);
    h_e->Draw("histo");

    ifs.close(); 

    Grafica->Run();
    return 0;
}