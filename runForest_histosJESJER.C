// Raghav Kunnawalkam Elayavalli
// Nov 11th 2015
// Rutgers
// for questions or comments: raghav.k.e at CERN dot CH

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
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

#include "boundaries.h"

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
}

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

using namespace std;

void runForest_histosJESJER(int startfile = 0,
			    int endfile = 1,
			    int radius = 3,
			    std::string algo = "Pu",
			    std::string jetType= "PF",
			    std::string kFoname="test_output.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  infile_Forest = "mergedfile.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 5; //6

  TChain * jetpbpb[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%s%d%sJetAnalyzer",algo.c_str(), radius, jetType.c_str());
  dir[3] = "akPu3CaloJetAnalyzer";
  dir[4] = "hiEvtAnalyzer";
  // dir[4] = "hltobject";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "t",
    "HiTree"
    // , "jetObjTree"
  };

  for(int t = 0;t<N;t++){
    jetpbpb[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends
  
  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    jetpbpb[0]->Add(filename_Forest.c_str());
    jetpbpb[1]->Add(filename_Forest.c_str());
    jetpbpb[2]->Add(filename_Forest.c_str());
    jetpbpb[3]->Add(filename_Forest.c_str());
    jetpbpb[4]->Add(filename_Forest.c_str());

    cout<<"filename: "<<filename_Forest<<endl;
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[3]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[4]->GetEntries() << endl;

    cout<<"Total number of events loaded in HiForest = "<<jetpbpb[2]->GetEntries()<<endl;

  }
  
  jetpbpb[2]->AddFriend(jetpbpb[0]);
  jetpbpb[2]->AddFriend(jetpbpb[1]);
  jetpbpb[2]->AddFriend(jetpbpb[4]);
  jetpbpb[3]->AddFriend(jetpbpb[0]);
  jetpbpb[3]->AddFriend(jetpbpb[1]);
  jetpbpb[3]->AddFriend(jetpbpb[4]);
  
  // Forest files 
  int nref_F;
  float pt_F[1000];
  float refpt_F[1000];
  float rawpt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float chMax_F[1000];
  float trkMax_F[1000];
  float chSum_F[1000];
  float phSum_F[1000];
  float neSum_F[1000];
  float trkSum_F[1000];
  float phMax_F[1000];
  float neMax_F[1000];
  float eMax_F[1000];
  float muMax_F[1000];
  float eSum_F[1000];
  float muSum_F[1000];
  float jtpu_F[1000];
  int   subid_F[1000];
  float refdrjt_F[1000];
  int refparton_F[1000];
  float pthat_F;
  int jet55_F;
  int jet65_F;
  int jet80_F;
  int L1_sj36_F;
  int L1_sj52_F;
  int L1_sj36_p_F;
  int L1_sj52_p_F;
  int jet55_p_F;
  int jet65_p_F;
  int jet80_p_F;
  float vz_F;
  int evt_F;
  int run_F;
  int lumi_F;
  int hiNpix_F;
  int hiNpixelTracks_F;
  int hiBin_F;
  float hiHF_F;
  int hiNtracks_F;
  int hiNtracksPtCut_F;
  int hiNtracksEtaCut_F;
  int hiNtracksEtaPtCut_F;
  int pcollisionEventSelection_F;
  int pHBHENoiseFilter_F;
  int pprimaryvertexFilter_F;
  int pVertexFilterCutGplus_F;

  float calopt_F[1000];
  jetpbpb[3]->SetBranchAddress("jtpt",&calopt_F);
  
  jetpbpb[4]->SetBranchAddress("evt",&evt_F);
  jetpbpb[4]->SetBranchAddress("run",&run_F);
  jetpbpb[4]->SetBranchAddress("lumi",&lumi_F);
  jetpbpb[4]->SetBranchAddress("hiBin",&hiBin_F);
  jetpbpb[4]->SetBranchAddress("hiHF", &hiHF_F);
  jetpbpb[4]->SetBranchAddress("hiNpix",&hiNpix_F);
  jetpbpb[4]->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks_F);
  jetpbpb[4]->SetBranchAddress("hiNtracks",&hiNtracks_F);
  jetpbpb[4]->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut_F);
  jetpbpb[4]->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut_F);
  jetpbpb[4]->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut_F);
  jetpbpb[4]->SetBranchAddress("vz",&vz_F);
  jetpbpb[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
  jetpbpb[0]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  jetpbpb[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  jetpbpb[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  jetpbpb[2]->SetBranchAddress("pthat",&pthat_F);
  jetpbpb[2]->SetBranchAddress("nref",&nref_F);
  // jetpbpb[2]->SetBranchAddress("subid",subid_F);
  jetpbpb[2]->SetBranchAddress("refdrjt",refdrjt_F);
  jetpbpb[2]->SetBranchAddress("refparton_flavor",refparton_F);
  jetpbpb[2]->SetBranchAddress("refpt",refpt_F);
  jetpbpb[2]->SetBranchAddress("jtpt",pt_F);
  jetpbpb[2]->SetBranchAddress("jteta",eta_F);
  jetpbpb[2]->SetBranchAddress("jtphi",phi_F);
  jetpbpb[2]->SetBranchAddress("rawpt",rawpt_F);
  jetpbpb[2]->SetBranchAddress("jtpu",jtpu_F);
  jetpbpb[2]->SetBranchAddress("chargedMax",chMax_F);
  jetpbpb[2]->SetBranchAddress("chargedSum",chSum_F);
  jetpbpb[2]->SetBranchAddress("trackMax",trkMax_F);
  jetpbpb[2]->SetBranchAddress("trackSum",trkSum_F);
  jetpbpb[2]->SetBranchAddress("photonMax",phMax_F);
  jetpbpb[2]->SetBranchAddress("photonSum",phSum_F);
  jetpbpb[2]->SetBranchAddress("neutralMax",neMax_F);
  jetpbpb[2]->SetBranchAddress("neutralSum",neSum_F);
  jetpbpb[2]->SetBranchAddress("eSum",eSum_F);
  jetpbpb[2]->SetBranchAddress("eMax",eMax_F);
  jetpbpb[2]->SetBranchAddress("muSum",muSum_F);
  jetpbpb[2]->SetBranchAddress("muMax",muMax_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v7",&jet55_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet55_v7_Prescl",&jet55_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v7",&jet65_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet65_v7_Prescl",&jet65_p_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v7",&jet80_F);
  jetpbpb[0]->SetBranchAddress("HLT_HIJet80_v7_Prescl",&jet80_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);

  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();

  TH1F * hJER_pt[nbins_pt][ncen];
  TH1F * hJER_eta[nbins_eta][ncen];
  TH1F * hpthat[ncen];
  TH1F * hpT[ncen];
  TH2F * hresponse_matrix[ncen];

  for(int i = 0;i<ncen;++i){

    hpthat[i] = new TH1F(Form("hpthat_cent%d",i),"",1000, 0, 1000);
    hpT[i] = new TH1F(Form("hpT_cent%d",i),"",1000, 0, 1000);
    hresponse_matrix[i] = new TH2F(Form("hresponse_matrix_cent%d",i),"",1000, 0, 1000, 1000, 0, 1000);
    
    for(int bin = 0; bin<nbins_pt; ++bin){
      hJER_pt[bin][i] = new TH1F(Form("hJER_ptbin_%d_cent%d",bin, i),Form("%d < pt < %d, %d-%d centrality", ptbins[bin], ptbins[bin+1], centbins[i], centbins[i+1]),150, 0, 3);
    }
    for(int bin = 0; bin<nbins_eta; ++bin){
      hJER_eta[bin][i] = new TH1F(Form("hJER_etabin_%d_cent%d", bin, i),Form("%2.2f < eta < %2.2f, %d-%d centrality", etabins[bin], etabins[bin+1], centbins[i], centbins[i+1]),150, 0, 3);
    }
    
  }

  
  // now start the event loop for each file. 
  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetpbpb[0]->GetEntries();
  if(printDebug) nentries = 10;
  TRandom rnd;

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    jetpbpb[0]->GetEntry(nEvt);
    jetpbpb[1]->GetEntry(nEvt);
    jetpbpb[2]->GetEntry(nEvt);
    jetpbpb[4]->GetEntry(nEvt);
    jetpbpb[3]->GetEntry(nEvt);
    
    if(pcollisionEventSelection_F==0) continue; 
    if(fabs(vz_F)>15) continue;
    
    int cBin = findBin(hiBin_F);//tells us the centrality of the event. 
  
    for( int jet = 0; jet<nref_F; jet++ ){

      if( fabs(eta_F[jet]) > 2.0 ) continue;
      //if( subid_F[jet] != 0) continue;
      if( refdrjt_F[jet] > 0.3 ) continue;
      if( pt_F[jet] > 3.0 * pthat_F ) continue;

      Float_t genpt = refpt_F[jet];
      Float_t recpt = pt_F[jet];
      Float_t rawpt = rawpt_F[jet];

      hresponse_matrix[cBin]->Fill(genpt, recpt);
      
      int ptbin = 0;
      for(int bin = 0; bin<nbins_pt; ++bin){
	if(genpt > ptbins[bin]) ptbin = bin;
      }
      hJER_pt[ptbin][cBin]->Fill((float)recpt/genpt);
      
      int etabin = 0;
      for(int bin = 0; bin<nbins_eta; ++bin){
	if(eta_F[jet] > etabins[bin]) etabin = bin;
      }
      hJER_eta[etabin][cBin]->Fill((float)recpt/genpt);

    }// jet loop
    
  }// event loop

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
