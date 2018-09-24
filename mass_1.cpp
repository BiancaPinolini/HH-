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
#include "norm1.h"

using namespace std;

int main(int argc, char** argv){
    if (argc < 3){
        cout << "Usage: " << argv[0] << "HH.root" << endl;
        return 1;
    } 

	TApplication * Grafica = new TApplication("App", 0, 0);
	//LETTURA DEI TTree
    TFile * sinput = TFile::Open( argv[1] );
    TFile * binput = TFile::Open( argv[2] );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
  
	int nbin = 100;
	double min = 0, max = 120;
	norm1 (signal, background, nbin, min, max);
	
    Grafica->Run();
	return 0;
}
