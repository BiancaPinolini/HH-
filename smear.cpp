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
    //in questo punto pongo "signal" or "background" a seconda del file da analizzare
    TTree *tree = new TTree ("tree", "Smearing on signal events");

    double P;
    double lep_PxI, lep_PyI, lep_PxF, lep_PyF, lep_PzF, lep_E; 
    double b1_PxI, b1_PyI, b1_PxF, b1_PyF, b1_PzF, b1_E;
    double b2_PxI, b2_PyI, b2_PxF, b2_PyF, b2_PzF, b2_E;
    double j1_PxI, j1_PyI, j1_PxF, j1_PyF, j1_PzF, j1_E;
    double j2_PxI, j2_PyI, j2_PxF, j2_PyF, j2_PzF, j2_E;
    double nu_PxI, nu_PyI, nu_PxF, nu_PyF;
    
    double somma_PxI = 0, somma_PyI = 0, somma_PxF = 0, somma_PyF = 0;

    TLorentzVector lep_I, lep_F;
    TLorentzVector b1_I, b1_F;
    TLorentzVector b2_I, b2_F;
    TLorentzVector j1_I, j2_I, j1_F, j2_F;

    //branches con i TLorentzVector associati a ogni particella
    tree->Branch("lep_I", &lep_I);
    tree->Branch("lep_F", &lep_F);
    tree->Branch("b1_I", &b1_I);
    tree->Branch("b1_F", &b1_F);
    tree->Branch("b2_I", &b2_I);
    tree->Branch("b2_F", &b2_F);
    tree->Branch("j1_I", &j1_I);
    tree->Branch("j1_F", &j1_F);
    tree->Branch("j2_I", &j2_I);
    tree->Branch("j2_F", &j2_F);

    //branches con Px e Py pre e post smearing del neutrino
    tree->Branch("nu_PxI", &nu_PxI);
    tree->Branch("nu_PxF", &nu_PxF);
    tree->Branch("nu_PyI", &nu_PyI);
    tree->Branch("nu_PyF", &nu_PyF);

    int iEv = 0, n = 1;
    int tau = 0, l = 0, q = 0, b = 0, b2 = 0;

    int ID = 0;

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
        if (iEv%1000 == 0){
            cout << "Event " << iEv << endl;
        }

        //Working on info of final state particles only
        for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size(); iPart++){

            if (reader.hepeup.ISTUP.at (iPart) == 1) {
                ID = reader.hepeup.IDUP.at(iPart);

                //ELETTRONI E MUONI
                if (abs(ID) == 11 || abs(ID) == 13 ){ 
                    l ++;
                    charged_leptons.push_back(iPart);
                    leptons.push_back(iPart);

                    lep_I.SetPxPyPzE(
                        reader.hepeup.PUP.at (iPart).at (0), //px
                        reader.hepeup.PUP.at (iPart).at (1), //py
                        reader.hepeup.PUP.at (iPart).at (2), //pz
                        reader.hepeup.PUP.at (iPart).at (3)  //E
                    );

                    P = lep_I.P();
                    lep_PxI = lep_I.Px();
                    lep_PyI = lep_I.Py();

                    somma_PxI = somma_PxI + lep_PxI;
                    somma_PyI = somma_PyI + lep_PyI;

                    P = smear(P, abs(ID));
                    lep_PxF = P * TMath::Sin(lep_I.Theta()) * TMath::Cos(lep_I.Phi());
                    lep_PyF = P * TMath::Sin(lep_I.Theta()) * TMath::Sin(lep_I.Phi());

                    somma_PxF = somma_PxF + lep_PxF;
                    somma_PyF = somma_PyF + lep_PyF;

                    //CALCOLO ALTRE COMPONENTI DEL QUADRIMOMENTO
                    lep_PzF = P * TMath::Cos(lep_I.Theta());
                    lep_E = sqrt(pow(P,2)+pow(lep_I.M(),2));
                    
                    lep_F.SetPxPyPzE(lep_PxF, lep_PyF, lep_PzF, lep_E);

                    if (l == 2) tau = 1;
                }

                //QUARK BEAUTY
                if (ID == 5){
                    b++;
                    quarks_b.push_back(iPart);

                    b1_I.SetPxPyPzE(
                        reader.hepeup.PUP.at (iPart).at (0), //px
                        reader.hepeup.PUP.at (iPart).at (1), //py
                        reader.hepeup.PUP.at (iPart).at (2), //pz
                        reader.hepeup.PUP.at (iPart).at (3)  //E
                    );

                    P = b1_I.P();
                    b1_PxI = b1_I.Px();
                    b1_PyI = b1_I.Py();

                    somma_PxI = somma_PxI + b1_PxI;
                    somma_PyI = somma_PyI + b1_PyI;

                    P = smear(P, abs(ID));
                    b1_PxF = P * TMath::Sin(b1_I.Theta()) * TMath::Cos(b1_I.Phi());
                    b1_PyF = P * TMath::Sin(b1_I.Theta()) * TMath::Sin(b1_I.Phi());

                    somma_PxF = somma_PxF + b1_PxF;
                    somma_PyF = somma_PyF + b1_PyF;

                    b1_PzF = P * TMath::Cos(b1_I.Theta());
                    b1_E = sqrt(pow(P,2)+pow(b1_I.M(),2));
                        
                    b1_F.SetPxPyPzE(b1_PxF, b1_PyF, b1_PzF, b1_E);

                    if (b == 2) tau = 1;
                }

                //QUARK ANTIBEAUTY
                if ( ID == -5 ){  
                    b2 ++;
                    quarks_b.push_back(iPart);
                    
                    b2_I.SetPxPyPzE(
                        reader.hepeup.PUP.at (iPart).at (0), //px
                        reader.hepeup.PUP.at (iPart).at (1), //py
                        reader.hepeup.PUP.at (iPart).at (2), //pz
                        reader.hepeup.PUP.at (iPart).at (3)  //E
                    );

                    P = b2_I.P();
                    b2_PxI = b2_I.Px();
                    b2_PyI = b2_I.Py();

                    somma_PxI = somma_PxI + b2_PxI;
                    somma_PyI = somma_PyI + b2_PyI;

                    P = smear(P, abs(ID));
                    b2_PxF = P * TMath::Sin(b2_I.Theta()) * TMath::Cos(b2_I.Phi());
                    b2_PyF = P * TMath::Sin(b2_I.Theta()) * TMath::Sin(b2_I.Phi());

                    somma_PxF = somma_PxF + b2_PxF;
                    somma_PyF = somma_PyF + b2_PyF;

                    b2_PzF = P * TMath::Cos(b2_I.Theta());
                    b2_E = sqrt(pow(P,2)+pow(b2_I.M(),2));
                        
                    b2_F.SetPxPyPzE(b2_PxF, b2_PyF, b2_PzF, b2_E);

                    if (b2 == 2) tau = 1;
                }

                // Other quarks in final state are from the jet 
                if ( abs(ID) < 7 && abs(ID) != 5 ){  
                    q++;

                    quarks_jet.push_back(iPart);
                    jets.push_back(iPart);

                    //PRIMO QUARK
                    if (q == 1) { 
                        j1_I.SetPxPyPzE(
                            reader.hepeup.PUP.at (iPart).at (0), //px
                            reader.hepeup.PUP.at (iPart).at (1), //py
                            reader.hepeup.PUP.at (iPart).at (2), //pz
                            reader.hepeup.PUP.at (iPart).at (3)  //E
                        );

                        P = j1_I.P();
                        j1_PxI = j1_I.Px();
                        j1_PyI = j1_I.Py();

                        somma_PxI = somma_PxI + j1_PxI;
                        somma_PyI = somma_PyI + j1_PyI;
                    
                        P = smear(P, 1);
                        j1_PxF = P * TMath::Sin(j1_I.Theta()) * TMath::Cos(j1_I.Phi());;
                        j1_PyF = P * TMath::Sin(j1_I.Theta()) * TMath::Sin(j1_I.Phi());

                        somma_PxF = somma_PxF + j1_PxF;
                        somma_PyF = somma_PyF + j1_PyF;

                        j1_PzF = P * TMath::Cos(j1_I.Theta());
                        j1_E = sqrt(pow(P,2)+pow(j1_I.M(),2));
                        
                        j1_F.SetPxPyPzE(j1_PxF, j1_PyF, j1_PzF, j1_E);
                    }
                    
                    //SECONDO QUARK
                    if (q == 2) { 
                        j2_I.SetPxPyPzE(
                            reader.hepeup.PUP.at (iPart).at (0), //px
                            reader.hepeup.PUP.at (iPart).at (1), //py
                            reader.hepeup.PUP.at (iPart).at (2), //pz
                            reader.hepeup.PUP.at (iPart).at (3)  //E
                        );

                        P = j2_I.P();
                        j1_PxI = j2_I.Px();
                        j1_PyI = j2_I.Py();

                        somma_PxI = somma_PxI + j1_PxI;
                        somma_PyI = somma_PyI + j1_PyI;
                    
                        P = smear(P, 1);
                        j1_PxF = P * TMath::Sin(j2_I.Theta()) * TMath::Cos(j2_I.Phi());;
                        j1_PyF = P * TMath::Sin(j2_I.Theta()) * TMath::Sin(j2_I.Phi());

                        somma_PxF = somma_PxF + j1_PxF;
                        somma_PyF = somma_PyF + j1_PyF;

                        j1_PzF = P * TMath::Cos(j2_I.Theta());
                        j1_E = sqrt(pow(P,2)+pow(j2_I.M(),2));
                        
                        j1_F.SetPxPyPzE(j1_PxF, j1_PyF, j1_PzF, j1_E);
                    }

                    if (q == 3) tau = 1;
                } 

                if (abs(ID) == 15)
                    tau = 1;
                
            }
        }

        if(tau == 0) {
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

        } else {
            tau = 0;
            q = 0;
            b = 0;
            b2 = 0;
            l = 0;
        }
    }
    ifs.close();
    output.Write();
    output.Close();    
    return 0;
}
