// macro to read in the latest hiForests and make plots to check the Jets.
// Author: Raghav Kunnawalkam Elayavalli
//         Rutgers,@ CERN for Run2
//         Nov 19th 2015

#include "boundaries.h"

using namespace std;

void Validate_Data_JetsPP(int startfile = 0,
			  int endfile = 1,
			  int radius = 4,
			  std::string algo = "Pu",
			  std::string jetType= "Calo",
			  std::string kFoname="test_output.root")
{

  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  infile_Forest = "pthat80_globalfit_marta.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 4; //6

  TChain * jetpbpb[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%s%d%sJetAnalyzer",algo.c_str(), radius, jetType.c_str());
  //dir[3] = "akPu3CaloJetAnalyzer";
  dir[3] = "hiEvtAnalyzer";
  // dir[4] = "hltobject";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    // "t",
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
    //jetpbpb[4]->Add(filename_Forest.c_str());

    cout<<"filename: "<<filename_Forest<<endl;
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jetpbpb[3]->GetEntries() << endl;
    //if(printDebug)cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
    //if(printDebug)cout << "Entries : " << jetpbpb[4]->GetEntries() << endl;

    cout<<"Total number of events loaded in HiForest = "<<jetpbpb[2]->GetEntries()<<endl;

  }
  
  jetpbpb[2]->AddFriend(jetpbpb[0]);
  jetpbpb[2]->AddFriend(jetpbpb[1]);
  jetpbpb[2]->AddFriend(jetpbpb[3]);
  //jetpbpb[3]->AddFriend(jetpbpb[0]);
  //jetpbpb[3]->AddFriend(jetpbpb[1]);
  //jetpbpb[3]->AddFriend(jetpbpb[4]);
  
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
  int jet40_F;
  int jet60_F;
  int jet80_F;
  int jet100_F;
  int L1_sj36_F;
  int L1_sj52_F;
  int L1_sj36_p_F;
  int L1_sj52_p_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  int jet100_p_F;
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

  //float calopt_F[1000];
  //jetpbpb[3]->SetBranchAddress("jtpt",&calopt_F);
  
  jetpbpb[3]->SetBranchAddress("evt",&evt_F);
  jetpbpb[3]->SetBranchAddress("run",&run_F);
  jetpbpb[3]->SetBranchAddress("lumi",&lumi_F);
  jetpbpb[3]->SetBranchAddress("hiBin",&hiBin_F);
  jetpbpb[3]->SetBranchAddress("hiHF", &hiHF_F);
  jetpbpb[3]->SetBranchAddress("hiNpix",&hiNpix_F);
  jetpbpb[3]->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks_F);
  jetpbpb[3]->SetBranchAddress("hiNtracks",&hiNtracks_F);
  jetpbpb[3]->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut_F);
  jetpbpb[3]->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut_F);
  jetpbpb[3]->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut_F);
  jetpbpb[3]->SetBranchAddress("vz",&vz_F);
  jetpbpb[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
  // jetpbpb[0]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  // jetpbpb[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  // jetpbpb[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
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
  jetpbpb[0]->SetBranchAddress("",&jet40_F);
  jetpbpb[0]->SetBranchAddress("",&jet40_p_F);
  jetpbpb[0]->SetBranchAddress("",&jet60_F);
  jetpbpb[0]->SetBranchAddress("",&jet60_p_F);
  jetpbpb[0]->SetBranchAddress("",&jet80_F);
  jetpbpb[0]->SetBranchAddress("",&jet80_p_F);
  jetpbpb[0]->SetBranchAddress("",&jet100_F);
  jetpbpb[0]->SetBranchAddress("",&jet100_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  jetpbpb[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);

  TFile *fout = new TFile(kFoname.c_str(),"RECREATE");
  fout->cd();




  
  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  

}
