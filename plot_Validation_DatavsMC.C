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
  hJetpT_Data->Print("base");
  TH1F * hAj_Data = (TH1F*)fData->Get("hAj");
  TH1F * hDeltaPhi_Data = (TH1F*)fData->Get("hDeltaPhi");

  TH1F * hevent = (TH1F*)fData->Get("hRunN_Vs_NJets");
  hevent->Print("base");
  TH1F * hJet80 = (TH1F*)fData->Get("hJet80");
  hJet80->Print("base");
  
  TFile * fMC = TFile::Open(MCFile.c_str());
  TH1F * pt2overpt1_MC = (TH1F*)fMC->Get("pt2overpt1");
  TH1F * hJetEta_MC = (TH1F*)fMC->Get("hJetEta");
  TH1F * hJetPhi_MC = (TH1F*)fMC->Get("hJetPhi");
  TH1F * hJetpT_MC = (TH1F*)fMC->Get("hJetpT");
  TH1F * hAj_MC = (TH1F*)fMC->Get("hAj");
  TH1F * hDeltaPhi_MC = (TH1F*)fMC->Get("hDeltaPhi");

  // get marta's histograms
  TFile * fMarta = TFile::Open("AnaResultsPFvsCaloJets4M.root");
  TList *lst = (TList*)fMarta->Get("anaPFvsCaloJet");
  TH3F * hpT3 = (TH3F*)lst->FindObject("fh3PtTrueEtaDeltaPt");
  int etabinlow = hpT3->GetYaxis()->FindBin(-2.00001);
  int etabinhigh = hpT3->GetYaxis()->FindBin(2.00001);
  TH1F * hPFCaloMatchpT = (TH1F*)hpT3->ProjectionX("hPFCaloMatchpT",etabinlow, etabinhigh);
  
  hPFCaloMatchpT->Print("base");
  
  // plot them on top of each other. 
  TCanvas * cAj = new TCanvas("Aj","",800,600);
  hAj_Data->SetMarkerStyle(20);
  hAj_Data->SetMarkerColor(kBlack);
  hAj_Data->SetAxisRange(0.0, 1.0, "X");
  hAj_Data->SetYTitle("Event Fraction");
  hAj_Data->SetXTitle("A_{j}");
  hAj_Data->SetTitle(" ");
  hAj_Data->DrawNormalized();
  
  hAj_MC->SetMarkerStyle(25);
  hAj_MC->SetMarkerColor(kRed);
  hAj_MC->DrawNormalized("same");

  putCMSPrel();
  //putPPLumi();
  
  TLegend * lAj = myLegend(0.4,0.6,0.7,0.9);
  lAj->AddEntry("","p_{T}^{lead}   > 100 GeV/c","");
  lAj->AddEntry("","p_{T}^{sublead}> 40  GeV/c","");
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
  hDeltaPhi_Data->SetTitle(" ");
  hDeltaPhi_Data->DrawNormalized();
  
  hDeltaPhi_MC->SetMarkerStyle(25);
  hDeltaPhi_MC->SetMarkerColor(kRed);
  hDeltaPhi_MC->DrawNormalized("same");

  putCMSPrel();
  //putPPLumi();
  
  TLegend * lDeltaPhi = myLegend(0.4,0.6,0.7,0.9);
  lDeltaPhi->AddEntry("","p_{T}^{lead}   > 100 GeV/c","");
  lDeltaPhi->AddEntry("","p_{T}^{sublead}> 40  GeV/c","");
  lDeltaPhi->AddEntry(hDeltaPhi_Data,"Data (Jet 80 trigger)","pl");
  lDeltaPhi->AddEntry(hDeltaPhi_MC,"MC (pthat 80)","pl");
  lDeltaPhi->Draw();

  cDeltaPhi->SaveAs(Form("DeltaPhi_datavsMC_%s_ak%s%d%s.pdf",coll.c_str(), algo.c_str(), radius, jetType.c_str()),"RECREATE");

    // plot them on top of each other. 
  TCanvas * cpt2overpt1 = new TCanvas("pt2overpt1","",800,600);
  pt2overpt1_Data->SetMarkerStyle(20);
  pt2overpt1_Data->SetMarkerColor(kBlack);
  pt2overpt1_Data->SetAxisRange(0.0, 1, "X");
  pt2overpt1_Data->SetYTitle("Event Fraction");
  pt2overpt1_Data->SetXTitle("p_{T}^{subLead}/p_{T}^{Lead}");
  pt2overpt1_Data->DrawNormalized();
 
  pt2overpt1_MC->SetMarkerStyle(25);
  pt2overpt1_MC->SetMarkerColor(kRed);
  //pt2overpt1_MC->SetAxisRange(0.0, 1.0, "X");
  pt2overpt1_MC->DrawNormalized("same");

  putCMSPrel();
  //putPPLumi();
  
  TLegend * lpt2overpt1 = myLegend(0.4,0.6,0.7,0.9);
  lpt2overpt1->AddEntry("","p_{T}^{lead}   > 100 GeV/c","");
  lpt2overpt1->AddEntry("","p_{T}^{sublead}> 40  GeV/c","");
  lpt2overpt1->AddEntry(pt2overpt1_Data,"Data (Jet 80 trigger)","pl");
  lpt2overpt1->AddEntry(pt2overpt1_MC,"MC (pthat 80)","pl");
  lpt2overpt1->Draw();

  cpt2overpt1->SaveAs(Form("pt2overpt1_datavsMC_%s_ak%s%d%s.pdf",coll.c_str(), algo.c_str(), radius, jetType.c_str()),"RECREATE");


  // plot them on top of each other. 
  TCanvas * chJetpT = new TCanvas("hJetpT","",800,600);
  chJetpT->SetLogy();
  
  hJetpT_Data->SetMarkerStyle(20);
  hJetpT_Data->SetMarkerColor(kBlack);
  hJetpT_Data->Rebin(5);
  //divideBinWidth(hJetpT_Data);
  hJetpT_Data->SetAxisRange(90.0, 500, "X");
  hJetpT_Data->GetXaxis()->SetTitleOffset(1.3);
  hJetpT_Data->SetYTitle("counts");
  hJetpT_Data->SetXTitle("Jet p_{T} (GeV/c)");
  hJetpT_Data->SetLineColor(kBlack);
  //hJetpT_Data->Scale(1./2426425);
  hJetpT_Data->Draw();

  hPFCaloMatchpT->SetMarkerStyle(33);
  hPFCaloMatchpT->SetMarkerColor(kBlue);
  divideBinWidth(hPFCaloMatchpT);
  //hPFCaloMatchpT->Scale(1./);
  //hPFCaloMatchpT->Draw("same");
  
  // hJetpT_MC->SetMarkerStyle(25);
  // hJetpT_MC->SetMarkerColor(kRed);
  // //hJetpT_MC->SetAxisRange(0.0, 1.0, "X");
  // hJetpT_MC->DrawNormalized("same");

  putCMSPrel(0.15,0.93,0.04);
  //putPPLumi();
  
  TLegend * lhJetpT = myLegend(0.30,0.6,0.7,0.9);
  // lhJetpT->AddEntry("","p_{T}^{lead} > 90 GeV","");
  lhJetpT->AddEntry("","Express data, pp #sqrt{s}=5.02 TeV","");
  lhJetpT->AddEntry("",Form("anti-k_{T} R=0.%d Particle Flow Jets",radius),"");
  lhJetpT->AddEntry("","|#eta|<2.0","");
  //lhJetpT->AddEntry("","","");
  lhJetpT->AddEntry(hJetpT_Data,"Data (HLT Jet80 trigger)","pl");
  //lhJetpT->AddEntry(hPFCaloMatchpT,"PF-Calo Matching","pl");
  lhJetpT->SetTextSize(0.04);
  lhJetpT->Draw();

  chJetpT->SaveAs(Form("hJetpT_%s_ak%s%d%s.pdf",coll.c_str(), algo.c_str(), radius, jetType.c_str()),"RECREATE");

  
}
