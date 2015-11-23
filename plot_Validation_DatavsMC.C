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
			     std::string coll = "PP",
			     std::string algo = "",
			     std::string jetType = "PF",
			     std::string DataFile = "PromptForestPP_DataJet80_ak4PF.root",
			     std::string MCFile = "PromptForestPP_MC_ak4PF.root")
{

  // get the data and MC histograms
  // these are simple histograms
  TFile * fData = TFile::Open(DataFile.c_str());
  TH1F * pt2overpt1_Data = (TH1F*)fData->Get("pt2overpt1");
  TH1F * hJetEta_Data = (TH1F*)fData->Get("hJetEta");
  TH1F * hJetPhi_Data = (TH1F*)fData->Get("hJetPhi");
  TH1F * hJetpT_Data = (TH1F*)fData->Get("hJetpT");
  TH1F * hAj_Data = (TH1F*)fData->Get("hAj");
  TH1F * hDeltaPhi_Data = (TH11F*)fData->Get("hDeltaPhi");
  
  TFile * fMC = TFile::Open(MCFile.c_str());
  TH1F * pt2overpt1_MC = (TH1F*)fMC->Get("pt2overpt1");
  TH1F * hJetEta_MC = (TH1F*)fMC->Get("hJetEta");
  TH1F * hJetPhi_MC = (TH1F*)fMC->Get("hJetPhi");
  TH1F * hJetpT_MC = (TH1F*)fMC->Get("hJetpT");
  TH1F * hAj_MC = (TH1F*)fMC->Get("hAj");
  TH1F * hDeltaPhi_MC = (TH11F*)fMC->Get("hDeltaPhi");
  
  // plot them on top of each other. 
  TCanvas * cAj = new TCanvas("Aj","",800,600);
  hAj_Data->SetMarkerStyle(20);
  hAj_Data->SetMarkerColor(kBlack);
  hAj_Data->SetAxisRange(0.0, 1.0, "X");
  hAj_Data->SetYTitle("Event Fraction");
  hAj_Data->SetXTitle("A_{j}");
  hAj_Data->DrawNormalized();
  
  hAj_MC->SetMarkerStyle(25);
  hAj_MC->SetMarkerColor(kRed);
  hAj_MC->DrawNormalized("same");

  putCMSPrel();
  putPPLumi();
  
  TLegend * lAj = myLegend(0.4,0.6,0.7,0.9);
  lAj->AddEntry("","p_{T}^{lead} > 90 GeV","");
  lAj->AddEntry("","p_{T}^{sublead}>20GeV","");
  lAj->AddEntry(hAj_Data,"Data (Jet 80 trigger)","pl");
  lAj->AddEntry(hAj_MC,"MC (pthat 80)","pl");
  lAj->Draw();

  cAj->SaveAs(Form("Aj_datavsMC_%s_ak%s%d%s.pdf",coll.c_str(), algo.c_str(), radius, jetType.c_str()),"RECREATE");

  // plot them on top of each other. 
  TCanvas * cDeltaPhi = new TCanvas("DeltaPhi","",800,600);
  hDeltaPhi_Data->SetMarkerStyle(20);
  hDeltaPhi_Data->SetMarkerColor(kBlack);
  hDeltaPhi_Data->SetAxisRange(0.0, 2.5, "X");
  hDeltaPhi_Data->SetYTitle("Event Fraction");
  hDeltaPhi_Data->SetXTitle("DeltaPhi");
  hDeltaPhi_Data->DrawNormalized();
  
  hDeltaPhi_MC->SetMarkerStyle(25);
  hDeltaPhi_MC->SetMarkerColor(kRed);
  hDeltaPhi_MC->DrawNormalized("same");

  putCMSPrel();
  putPPLumi();
  
  TLegend * lDeltaPhi = myLegend(0.4,0.6,0.7,0.9);
  lDeltaPhi->AddEntry("","p_{T}^{lead} > 90 GeV","");
  lDeltaPhi->AddEntry("","p_{T}^{sublead}>20GeV","");
  lDeltaPhi->AddEntry(hDeltaPhi_Data,"Data (Jet 80 trigger)","pl");
  lDeltaPhi->AddEntry(hDeltaPhi_MC,"MC (pthat 80)","pl");
  lDeltaPhi->Draw();

  cDeltaPhi->SaveAs(Form("DeltaPhi_datavsMC_%s_ak%s%d%s.pdf",coll.c_str(), algo.c_str(), radius, jetType.c_str()),"RECREATE");

    // plot them on top of each other. 
  TCanvas * cpt2overpt1 = new TCanvas("pt2overpt1","",800,600);
  hpt2overpt1_Data->SetMarkerStyle(20);
  hpt2overpt1_Data->SetMarkerColor(kBlack);
  hpt2overpt1_Data->SetAxisRange(0.0, 1, "X");
  hpt2overpt1_Data->SetYTitle("Event Fraction");
  hpt2overpt1_Data->SetXTitle("p_{T}^{subLead}/p_{T}^{Lead}");
  hpt2overpt1_Data->DrawNormalized();
  
  hpt2overpt1_MC->SetMarkerStyle(25);
  hpt2overpt1_MC->SetMarkerColor(kRed);
  //hpt2overpt1_MC->SetAxisRange(0.0, 1.0, "X");
  hpt2overpt1_MC->DrawNormalized("same");

  putCMSPrel();
  putPPLumi();
  
  TLegend * lpt2overpt1 = myLegend(0.4,0.6,0.7,0.9);
  lpt2overpt1->AddEntry("","p_{T}^{lead} > 90 GeV","");
  lpt2overpt1->AddEntry("","p_{T}^{sublead}>20GeV","");
  lpt2overpt1->AddEntry(hpt2overpt1_Data,"Data (Jet 80 trigger)","pl");
  lpt2overpt1->AddEntry(hpt2overpt1_MC,"MC (pthat 80)","pl");
  lpt2overpt1->Draw();

  cpt2overpt1->SaveAs(Form("pt2overpt1_datavsMC_%s_ak%s%d%s.pdf",coll.c_str(), algo.c_str(), radius, jetType.c_str()),"RECREATE");


  
}
