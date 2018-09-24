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
#include "normCS.h"

using namespace std;

int main(){
	TApplication * Grafica = new TApplication("App", 0, 0);

	//LETTURA DEI TTree
    TFile * sinput = TFile::Open ("HH.root");
    TFile * b1input = TFile::Open ("ttbar.root");
    TFile * b2input = TFile::Open ("WWZ.root");
    TFile * b3input = TFile::Open ("WZZ.root");
    
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background1 = (TTree*)b1input->Get("tree");
    TTree * background2 = (TTree*)b2input->Get("tree");
    TTree * background3 = (TTree*)b3input->Get("tree");
  
	int nbin = 100;
	double min = 0, max = 150;
	normCS (signal, background1, background2, background3, nbin, min, max);
	
    Grafica->Run();
	return 0;
}
