#ifndef __boundaries_h_
#define __boundaries_h_

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

#define NOBJECT_MAX 16384

const bool printDebug=false;
const bool doBjets = false;

// boundaries of the pt bins, cent bins and eta bins for the runForest and plot macros.

static const double pthat[12] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540, 2000};
static const double xsecs[12] = {5.269E-01, 3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05 , 1.223E-05, 3.031E-06 , 7.746E-07, 1.410E-07 , 3.216E-08, 1.001E-08 , 0.0};

static const double weight_xsec[9] = { 7.20357e-07, 4.51655e-08, 2.6964e-09, 2.77274e-10, 3.1878e-11, 3.87126e-12, 1.62138e-12, 1.09471e-12, 4.40012e-13};
static const int nentries_file[9] = { 0, 333206, 250567, 395126, 368126, 366982, 392206, 181018, 50455};

const int ptbins[] = {50, 80, 120, 170, 220, 300};
const double ptbins_bound[] = {50, 80, 120, 170, 220, 300};
const int nbins_pt = sizeof(ptbins)/sizeof(int) -1;


const double etabins[] = {-5.191, -2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522, 0.522, 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 5.191};
const int nbins_eta = sizeof(etabins)/sizeof(double) -1;

/* const int ncen=10; */
/* const int centbins[ncen+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}; */
/* const char *cdir[ncen]  = {"010","1020","2030","3040","4050","5060","6070","7080","8090","90100"}; */
/* const char *ccent[ncen] = {"0-10%","10-20%","20-30%","30-50%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"}; */

const int knj = 1;
std::string srad[knj]={"4"};

double xmin=ptbins[0];
double xmax=ptbins[nbins_pt];

int PbPbcentbins[5] = {0, 20, 60, 100, 200};
string PbPbcdir[5]  = {"010","1030","3050","50100", "PP"};
string PbPbccent[5] = {"0-10%","10-30%","30-50%","50-100%", "PP"};
	

// Adding PF candidates plots
// 1) 2D histograms for eta vs pT for the candidate types ( for PbPb this needs to be before and after Vs subtraction)
// 2) 1D histograms for candidate pT/ jet pT for candidates inside the jet radius, for each candidate type
// Particle::pdgId_ PFCandidate::particleId_
// PFCandidate::ParticleType Particle
// 0           0  X          unknown, or dummy 
// +211, -211  1  h          charged hadron 
// +11, -11    2  e          electron 
// +13, -13    3  mu         muon 
// 22          4  gamma      photon 
// 130         5  h0         neutral hadron 
// 130         6  h_HF       hadronic energy in an HF tower 
// 22          7  egamma_HF  electromagnetic energy in an HF tower
  
const int PFType = 8;
std::string PFCandType[PFType] = {"unknown",
				  "chargedHadron",
				  "electron",
				  "muon",
				  "photon",
				  "neutralHadron",
				  "HadEnergyinHF",
				  "EMEnergyinHF"};


/*
int findBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins               
  if(bin<20)ibin=0; //! 0-10%
  else if(bin>=20  && bin<40 )ibin=1; //! 10-20%
  else if(bin>=40  && bin<60 )ibin=2; //! 20-30%
  else if(bin>=60  && bin<80 )ibin=3; //! 30-40%
  else if(bin>=80  && bin<100 )ibin=4; //! 40-50%
  else if(bin>=100  && bin<120 )ibin=5; //! 50-60%
  else if(bin>=120  && bin<140 )ibin=6; //! 60-70%
  else if(bin>=140  && bin<160 )ibin=7; //! 70-80%
  else if(bin>=160  && bin<180 )ibin=8; //! 80-90%
  else if(bin>=180  && bin<200 )ibin=9; //! 90-100%
  return ibin;

}*/

int findBin(int bin)
{
	int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins    
  if(bin<20)ibin=0; //! 0-10%
  else if(bin>=20  && bin<60 )ibin=1; //! 10-30%
  else if(bin>=60  && bin<100 )ibin=2; //! 30-50%
  else if(bin>=100  && bin<200 )ibin=3; //! 50-100%
  return ibin;
}



#define pi 3.14159265

float deltaphi(float phi1, float phi2)
{
  //float pi=TMath::Pi();
  
  float dphi = TMath::Abs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;

  return TMath::Abs(dphi);
}
//static const int nbins_pt = 14;
//static const double boundaries_pt[nbins_pt+1] = {
//  49, 56, 64, 74, 84, 97, 114, 133,
//  153, 174, 196, 220, 245, 272, 300
//};

//these are the only radii we are interested for the RAA analysis: 2,3,4,5
//static const int no_radius = 3; 
//static const int list_radius[no_radius] = {2,3,4};


float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  float dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}

// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

struct Jet{
  int id;
  float pt;
};
bool compare_pt(Jet jet1, Jet jet2);
bool compare_pt(Jet jet1, Jet jet2){
  return jet1.pt > jet2.pt;
}


float LjCut = 80.0;
float SbjCut = 30.0;

#endif
