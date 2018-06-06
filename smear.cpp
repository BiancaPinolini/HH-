#include "LHEF.h"
#include<iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include <fstream>
#include <algorithm>
#include "TLorentzVector.h"
#include <TF1.h>
#include <TRandom.h>
#include "smear.h"

using namespace std;

int main(int argc, char** argv){
    if (argc < 3){
        cout << "Usage: " << argv[0] << " output.root file1.lhe.." << endl;
        return 1;
    }

    int max = -1;

    if (argc > 3) max = atoi(argv[3]);

    char* rootfile = argv[1];
    TFile output(rootfile, "RECREATE");
    TTree *tree = new TTree ("tree", "Smearing on signal events");

    double P;
    double lep_PxI, lep_PyI, b1_PxI, b1_PyI, b2_PxI, b2_PyI, j1_PxI, j1_PyI, j2_PxI, j2_PyI, nu_PxI, nu_PyI;
    double somma_PxI = 0, somma_PyI = 0;

    tree->Branch("lep_PxI", &lep_PxI);
    tree->Branch("b1_PxI", &b1_PxI);
    tree->Branch("b2_PxI", &b2_PxI);
    tree->Branch("j1_PxI", &j1_PxI);
    tree->Branch("j2_PxI", &j2_PxI);
    tree->Branch("nu_PxI", &nu_PxI);

    tree->Branch("lep_PyI", &lep_PyI);
    tree->Branch("b1_PyI", &b1_PyI);
    tree->Branch("b2_PyI", &b2_PyI);
    tree->Branch("j1_PyI", &j1_PyI);
    tree->Branch("j2_PyI", &j2_PyI);
    tree->Branch("nu_PyI", &nu_PyI);

    double lep_PxF, lep_PyF, b1_PxF, b1_PyF, b2_PxF, b2_PyF, j1_PxF, j1_PyF, j2_PxF, j2_PyF, nu_PxF, nu_PyF;
    double somma_PxF = 0, somma_PyF = 0;
    double lep_PzF, lep_E, b1_PzF, b1_E, b2_PzF, b2_E, j1_PzF, j1_E, j2_PzF, j2_E;

    tree->Branch("lep_PxF", &lep_PxF);
    tree->Branch("b1_PxF", &b1_PxF);
    tree->Branch("b2_PxF", &b2_PxF);
    tree->Branch("j1_PxF", &j1_PxF);
    tree->Branch("j2_PxF", &j2_PxF);
    tree->Branch("nu_PxF", &nu_PxF);

    tree->Branch("lep_PyF", &lep_PyF);
    tree->Branch("b1_PyF", &b1_PyF);
    tree->Branch("b2_PyF", &b2_PyF);
    tree->Branch("j1_PyF", &j1_PyF);
    tree->Branch("j2_PyF", &j2_PyF);
    tree->Branch("nu_PyF", &nu_PyF);

    int iEv = 0, n = 1;
    int tau = 0;

    int ID = 0;

    TLorentzVector jet_j1_mom, jet_j2_mom;

    //Save quadrimomentum of all particles mapped with position in the event
    map<int, TLorentzVector> momenta;

    ifstream ifs(argv[2]);
    LHEF::Reader reader(ifs);
    
    //Loop over all the events 
    while ( reader.readEvent() ){
        iEv++;
        if (max > 0) {
            if (iEv > max) break ; 
        }
        vector<int> charged_leptons;
        vector<int> leptons;
        vector<int> quarks_b;
        vector<int> photons;
        vector<int> quarks_jet;
        vector<int> jets;

        //per controllare che vada tutto bene stampo un #evento ogni 1000
        //if (iEv%1000 == 0){
            cout << "Event " << iEv << endl;
        //}

        //Working on info of final state particles only
        for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size(); iPart++){

            if (reader.hepeup.ISTUP.at (iPart) == 1) {
                ID = reader.hepeup.IDUP.at(iPart);

                TLorentzVector momentum(
                    reader.hepeup.PUP.at (iPart).at (0), //PG px
                    reader.hepeup.PUP.at (iPart).at (1), //PG py
                    reader.hepeup.PUP.at (iPart).at (2), //PG pz
                    reader.hepeup.PUP.at (iPart).at (3) //PG E
                );

                momenta[iPart] =  momentum;

                //ELETTRONI E MUONI
                if( abs(ID) == 11 || abs(ID) == 13 ){ 
                    charged_leptons.push_back(iPart);
                    leptons.push_back(iPart);

                    P = momentum.P();

                    lep_PxI = momentum.Px();
                    lep_PyI = momentum.Py();

                    somma_PxI = somma_PxI + lep_PxI;
                    somma_PyI = somma_PyI + lep_PyI;

                    P = smear(P, abs(ID));
                    lep_PxF = P * TMath::Sin(momentum.Theta()) * TMath::Cos(momentum.Phi());
                    lep_PyF = P * TMath::Sin(momentum.Theta()) * TMath::Sin(momentum.Phi());

                    somma_PxF = somma_PxF + lep_PxF;
                    somma_PyF = somma_PyF + lep_PyF;

                    //CALCOLO ALTRE COMPONENTI DEL QUADRIMOMENTO
                    lep_PzF = P * TMath::Cos(momentum.Theta());
                    lep_E = sqrt(pow(P,2)+pow(momentum.M(),2));
                    
                    momentum.SetPxPyPzE(lep_PxF, lep_PyF, lep_PzF, lep_E);
                    //da questo posso ricavare tutte le variabili che voglio dopo lo smearing mantenendo le correlazioni
                }
                
                //QUARK BEAUTY
                if ( ID == 5){  
                    quarks_b.push_back(iPart);

                    P = momentum.P();
                    b1_PxI = momentum.Px();
                    b1_PyI = momentum.Py();

                    somma_PxI = somma_PxI + b1_PxI;
                    somma_PyI = somma_PyI + b1_PyI;

                    P = smear(P, abs(ID));
                    b1_PxF = P * TMath::Sin(momentum.Theta()) * TMath::Cos(momentum.Phi());
                    b1_PyF = P * TMath::Sin(momentum.Theta()) * TMath::Sin(momentum.Phi());

                    somma_PxF = somma_PxF + b1_PxF;
                    somma_PyF = somma_PyF + b1_PyF;

                    b1_PzF = P * TMath::Cos(momentum.Theta());
                    b1_E = sqrt(pow(P,2)+pow(momentum.M(),2));
                    
                    momentum.SetPxPyPzE(b1_PxF, b1_PyF, b1_PzF, b1_E);
                }

                //QUARK ANTIBEAUTY
                if ( ID == -5 ){  
                    quarks_b.push_back(iPart);

                    P = momentum.P();
                    b2_PxI = momentum.Px();
                    b2_PyI = momentum.Py();

                    somma_PxI = somma_PxI + b2_PxI;
                    somma_PyI = somma_PyI + b2_PyI;

                    P = smear(P, abs(ID));
                    b2_PxF = P * TMath::Sin(momentum.Theta()) * TMath::Cos(momentum.Phi());
                    b2_PyF = P * TMath::Sin(momentum.Theta()) * TMath::Sin(momentum.Phi());

                    somma_PxF = somma_PxF + b2_PxF;
                    somma_PyF = somma_PyF + b2_PyF;

                    b2_PzF = P * TMath::Cos(momentum.Theta());
                    b2_E = sqrt(pow(P,2)+pow(momentum.M(),2));
                    
                    momentum.SetPxPyPzE(b2_PxF, b2_PyF, b2_PzF, b2_E);
                }

                // Other quarks in final state are from the jet 
                if ( abs(ID) < 7 && abs(ID) != 5 ){  
                    quarks_jet.push_back(iPart);
                    jets.push_back(iPart);
                } 

                if (abs(ID) == 15)
                    tau = 1;
                
            }
        }

        if(tau == 0) {
            jet_j1_mom = momenta[jets.at(0)];  
            jet_j2_mom = momenta[jets.at(1)];

            //PRIMO QUARK
            P = jet_j1_mom.P();
            j1_PxI = jet_j1_mom.Px();
            j1_PyI = jet_j1_mom.Py();

            somma_PxI = somma_PxI + j1_PxI;
            somma_PyI = somma_PyI + j1_PyI;
        
            P = smear(P, 1);
            j1_PxF = P * TMath::Sin(jet_j1_mom.Theta()) * TMath::Cos(jet_j1_mom.Phi());;
            j1_PyF = P * TMath::Sin(jet_j1_mom.Theta()) * TMath::Sin(jet_j1_mom.Phi());

            somma_PxF = somma_PxF + j1_PxF;
            somma_PyF = somma_PyF + j1_PyF;

            j1_PzF = P * TMath::Cos(jet_j1_mom.Theta());
            j1_E = sqrt(pow(P,2)+pow(jet_j1_mom.M(),2));
            
            jet_j1_mom.SetPxPyPzE(j1_PxF, j1_PyF, j1_PzF, j1_E);

            //SECONDO QUARK
            P = jet_j2_mom.P();
            j2_PxI = jet_j2_mom.Px();
            j2_PyI = jet_j2_mom.Py();

            somma_PxI = somma_PxI + j2_PxI;
            somma_PyI = somma_PyI + j2_PyI;
        
            P = smear(P, 1);
            j2_PxF = P * TMath::Sin(jet_j2_mom.Theta()) * TMath::Cos(jet_j2_mom.Phi());;
            j2_PyF = P * TMath::Sin(jet_j2_mom.Theta()) * TMath::Sin(jet_j2_mom.Phi());

            somma_PxF = somma_PxF + j2_PxF;
            somma_PyF = somma_PyF + j2_PyF;

            j2_PzF = P * TMath::Cos(jet_j2_mom.Theta());
            j2_E = sqrt(pow(P,2)+pow(jet_j2_mom.M(),2));
            
            jet_j2_mom.SetPxPyPzE(j2_PxF, j2_PyF, j2_PzF, j2_E);

            //MISSING ENERGY
            nu_PxI = - somma_PxI;
            nu_PyI = - somma_PyI;

            nu_PxF = - somma_PxF;
            nu_PyF = - somma_PyF;

            //RIEMPIO IL TREE
            tree->Fill();

            somma_PxI = 0;
            somma_PyI = 0;

            somma_PxF = 0;
            somma_PyF = 0;

        } else tau = 0;
    }
    ifs.close();
    output.Write();
    output.Close();    
    return 0;
}
