#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include <TVirtualFitter.h>

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "boundaries.h"
#include "plot.h"

using namespace std;

const int digi=3;

void plot_Validatin_DatavsMC(int radius = 4,
			     std::string coll = "PbPb",
			     std::string algo = "Pu",
			     std::string jetType = "PF",
			     std::string DataFile = "Data.root",
			     std::string MCFile = "MC.root")
{

  // get the data and MC histograms
  // these are simple histograms
  TFile * fData = TFile::Open(DataFile.c_str());
  TH1F * pt2overpt1_Data = (TH1F*)fData->Get("pt2overpt1");
  TH1F * hJetEta_Data = (TH1F*)fData->Get("hJetEta");
  TH1F * hJetPhi_Data = (TH1F*)fData->Get("hJetPhi");
  TH1F * hJetpT_Data = (TH1F*)fData->Get("hJetpT");
  TH1F * hAj_Data = (TH1F*)fData->Get("hAj");
  
  TFile * fMC = TFile::Open(MCFile.c_str());
  TH1F * pt2overpt1_MC = (TH1F*)fMC->Get("pt2overpt1");
  TH1F * hJetEta_MC = (TH1F*)fMC->Get("hJetEta");
  TH1F * hJetPhi_MC = (TH1F*)fMC->Get("hJetPhi");
  TH1F * hJetpT_MC = (TH1F*)fMC->Get("hJetpT");
  TH1F * hAj_MC = (TH1F*)fMC->Get("hAj");
  
  // plot them on top of each other. 
  TCanvas * cAj = new TCanvas("Aj","",800,600);
  
    

}
