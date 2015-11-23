// macro to read in the latest hiForests and make plots to check the Jets.
// Author: Raghav Kunnawalkam Elayavalli, Kurt Jung 
//         @ CERN for Run2
//         Nov 19th 2015

#include "boundaries.h"

using namespace std;

void Validate_Jets(int startfile = 0,
		   int endfile = 1,
		   int radius = 4,
		   std::string coll= "PP",
		   std::string run= "Data",
		   std::string jetType= "Calo",
		   std::string algo= "",
		   std::string kFoname="PromptForest")
{

  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool doBjets = false;
  bool skipPho50 = true;
  bool printDebug = false;
  bool doDijetImbalance = true;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  infile_Forest = Form("%s_%s_forests.txt", coll.c_str(), run.c_str());
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 4; //6

  TChain * jtTree[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = "test";
  if(coll == "PbPb") dir[2] = Form("ak%s%d%sJetAnalyzer", algo.c_str(), radius, jetType.c_str());
  if(coll == "PP") dir[2] = Form("ak%d%sJetAnalyzer", radius, jetType.c_str());
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
    jtTree[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends
  
  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    jtTree[0]->Add(filename_Forest.c_str());
    jtTree[1]->Add(filename_Forest.c_str());
    jtTree[2]->Add(filename_Forest.c_str());
    jtTree[3]->Add(filename_Forest.c_str());
    //jtTree[4]->Add(filename_Forest.c_str());

    cout<<"filename: "<<filename_Forest<<endl;
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[3]->GetEntries() << endl;
    //if(printDebug)cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
    //if(printDebug)cout << "Entries : " << jtTree[4]->GetEntries() << endl;

    cout<<"Total number of events loaded in HiForest = "<<jtTree[2]->GetEntries()<<endl;

  }
  
  jtTree[2]->AddFriend(jtTree[0]);
  jtTree[2]->AddFriend(jtTree[1]);
  jtTree[2]->AddFriend(jtTree[3]);
  //jtTree[3]->AddFriend(jtTree[0]);
  //jtTree[3]->AddFriend(jtTree[1]);
  //jtTree[3]->AddFriend(jtTree[4]);
  
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
  float discr_ssvHighEff_F[1000];
  float discr_ssvHighPur_F[1000];
  float discr_csv_F[1000];
  int jet40_F;
  int jet60_F;
  int jet80_F;
  int jet100_F;
  int jetMB_F;
  int L1_sj36_F;
  int L1_sj52_F;
  int L1_sj36_p_F;
  int L1_sj52_p_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  int jet100_p_F;
  int jetMB_p_F;
  int photon50_F;
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
  //jtTree[3]->SetBranchAddress("jtpt",&calopt_F);
  
  jtTree[3]->SetBranchAddress("evt",&evt_F);
  jtTree[3]->SetBranchAddress("run",&run_F);
  jtTree[3]->SetBranchAddress("lumi",&lumi_F);
  jtTree[3]->SetBranchAddress("hiBin",&hiBin_F);
  jtTree[3]->SetBranchAddress("hiHF", &hiHF_F);
  jtTree[3]->SetBranchAddress("hiNpix",&hiNpix_F);
  jtTree[3]->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks_F);
  jtTree[3]->SetBranchAddress("hiNtracks",&hiNtracks_F);
  jtTree[3]->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut_F);
  jtTree[3]->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut_F);
  jtTree[3]->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut_F);
  jtTree[3]->SetBranchAddress("vz",&vz_F);
  jtTree[1]->SetBranchAddress("PAcollisionEventSelection",&pcollisionEventSelection_F);
  // jtTree[0]->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter_F);
  // jtTree[0]->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter_F);
  // jtTree[0]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  if(run == "MC") jtTree[2]->SetBranchAddress("pthat",&pthat_F);
  jtTree[2]->SetBranchAddress("nref",&nref_F);
  // jtTree[2]->SetBranchAddress("subid",subid_F);
  if(run == "MC") jtTree[2]->SetBranchAddress("refdrjt",refdrjt_F);
  if(run == "MC") jtTree[2]->SetBranchAddress("refparton_flavor",refparton_F);
  if(run == "MC") jtTree[2]->SetBranchAddress("refpt",refpt_F);
  if(doBjets){
    jtTree[2]->SetBranchAddress("discr_ssvHighEff",discr_ssvHighEff_F);
    jtTree[2]->SetBranchAddress("discr_ssvHighPur",discr_ssvHighPur_F);
    jtTree[2]->SetBranchAddress("discr_csvSimple",discr_csv_F);
  }
  jtTree[2]->SetBranchAddress("jtpt",pt_F);
  jtTree[2]->SetBranchAddress("jteta",eta_F);
  jtTree[2]->SetBranchAddress("jtphi",phi_F);
  jtTree[2]->SetBranchAddress("rawpt",rawpt_F);
  jtTree[2]->SetBranchAddress("jtpu",jtpu_F);
  jtTree[2]->SetBranchAddress("chargedMax",chMax_F);
  jtTree[2]->SetBranchAddress("chargedSum",chSum_F);
  jtTree[2]->SetBranchAddress("trackMax",trkMax_F);
  jtTree[2]->SetBranchAddress("trackSum",trkSum_F);
  jtTree[2]->SetBranchAddress("photonMax",phMax_F);
  jtTree[2]->SetBranchAddress("photonSum",phSum_F);
  jtTree[2]->SetBranchAddress("neutralMax",neMax_F);
  jtTree[2]->SetBranchAddress("neutralSum",neSum_F);
  jtTree[2]->SetBranchAddress("eSum",eSum_F);
  jtTree[2]->SetBranchAddress("eMax",eMax_F);
  jtTree[2]->SetBranchAddress("muSum",muSum_F);
  jtTree[2]->SetBranchAddress("muMax",muMax_F);
  jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part1_v1",&jetMB_F);
  //jtTree[0]->SetBranchAddress("",&jetMB_p_F);  
  jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet40_Eta5p1_v1", jetType.c_str()),&jet40_F);
  //jtTree[0]->SetBranchAddress("",&jet40_p_F);
  jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet60_Eta5p1_v1", jetType.c_str()),&jet60_F);
  //jtTree[0]->SetBranchAddress("",&jet60_p_F);
  jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet80_Eta5p1_v1", jetType.c_str()),&jet80_F);
  //jtTree[0]->SetBranchAddress("",&jet80_p_F);
  jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet100_Eta5p1_v1", jetType.c_str()),&jet100_F);
  //jtTree[0]->SetBranchAddress("",&jet100_p_F);
  jtTree[0]->SetBranchAddress("HLT_HISinglePhoton50_Eta1p5_v1",&photon50_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);

  std::string rad = Form("%d",radius);
  std::string end = Form("%d",endfile);
  
  TFile *fout = new TFile((kFoname+coll+"_"+run+"_ak"+algo+rad+jetType+"_"+end+".root").c_str(),"RECREATE");
  fout->cd();

  // Add the histograms necessary for the validation,
  // Aj (eta bins), Relative Response (eta bins), PF cands, rechits,
  // dijets are selected by delta phi > 2 pi / 3
  // there is also an alpha cut off which is pt third / pt avg < 0.2

  // TH1F * hRelResponse[nbins_pt][nbins_eta];
  // TH1F * hAj[nbins_pt][nbins_eta];

  // for(int npt = 0; npt<nbins_pt; ++npt){
  //   for(int neta = 0; neta<nbins_eta; ++neta){
  //     hRelResponse[npt][neta] = new TH1F(Form("hRelResponse_ptbin%d_etabin%d",npt, neta),Form("Relative Response in %d < p_{T}^{avg} < %d, %2.2f < #eta_{LeadJet} < %2.2f", ptbins[npt], ptbins[npt+1], etabins[neta], etabins[neta+1]),300, -3, 3);
  //     hAj[npt][neta] = new TH1F(Form("hAj_ptbin%d_etabin%d",npt, neta),Form("Aj in %d < p_{T}^{avg} < %d, %2.2f < #eta_{LeadJet} < %2.2f", ptbins[npt], ptbins[npt+1], etabins[neta], etabins[neta+1]),300, -3, 3);
  //   }
  // }

  // Add the Jet pT spectra histograms
  //TH1F * hJtpt[nbins_cent+1] = new TH1F("hJtpt","Jet Spectra",1000, 0, 1000);

  //TH1F * hRelResponse_40_pt_80 = new TH1F("hRelResponse_40_pt_80","Relative Response, 40 < p_{T}^{avg} < 80, |#eta|<2;#frac{2+B}{2-B}, B = #frac{p_{T}^{Lead} - p_{T}^{subLead}}{p_{T}^{avg}};counts",300, -3, 3);
  //TH1F * hRelResponse_80_pt_120 = new TH1F("hRelResponse_80_pt_120","Relative Response, 80 < p_{T}^{avg} < 120, |#eta|<2;#frac{2+B}{2-B}, B = #frac{p_{T}^{Lead} - p_{T}^{subLead}}{p_{T}^{avg}};counts",300, -3, 3);
  //TH1F * hRelResponse_120_pt_300 = new TH1F("hRelResponse_120_pt_300","Relative Response, 120 < p_{T}^{avg} < 300, |#eta|<2;#frac{2+B}{2-B}, B = #frac{p_{T}^{Lead} - p_{T}^{subLead}}{p_{T}^{avg}};counts",300, -3, 3);
  
  TH1F *hRelResponse[nbins_pt];
  TH1F *hRelResponse_outereta[nbins_pt];
  TH1F *hRelResponse_innereta[nbins_pt];

  TH2F * hRunN_vs_NJets = new TH2F("hRunN_Vs_NJets","",719, 261445, 262164, 50, 0, 50);
  
  for(int npt = 0; npt<nbins_pt; ++npt){
    hRelResponse[npt] = new TH1F(Form("hRelResponse_ptbin%d",npt),Form("Relative Response in %d < p_{T}^{avg} < %d", ptbins[npt], ptbins[npt+1]), 1000, -3, 3);
    hRelResponse_outereta[npt] = new TH1F(Form("hRelResponse_outereta_ptbin%d",npt),Form("Relative Response in %d < p_{T}^{avg} < %d", ptbins[npt], ptbins[npt+1]), 1000, -3, 3);
    hRelResponse_innereta[npt] = new TH1F(Form("hRelResponse_innereta_ptbin%d",npt),Form("Relative Response in %d < p_{T}^{avg} < %d", ptbins[npt], ptbins[npt+1]), 1000, -3, 3);
  }

  vector<Float_t> j_array_pt;
  vector<Float_t> j_array_phi;
  vector<Float_t> j_array_eta;

  TH1F * hJEC_eta0_05 = new TH1F("hJEC_eta0_05","JEC applied in the forest, |#eta|<0.5",100, 0, 5);
  TH2F * hJEC_vs_rawpT_eta0_05 = new TH2F("hJEC_vs_rawpT_eta0_05","JEC applied in the forest vs raw pT, |#eta|<0.5",200, 0, 200, 300, 0.5, 2);

  TH1F * hJEC_eta05_10 = new TH1F("hJEC_eta05_10","JEC applied in the forest, 0.5<|#eta|<1.0, ",100, 0, 5);
  TH2F * hJEC_vs_rawpT_eta05_10 = new TH2F("hJEC_vs_rawpT_eta05_10","JEC applied in the forest vs raw pT, 0.5<|#eta|<1.0",200, 0, 200, 300, 0.5, 2);

  TH1F * hJEC_eta10_15 = new TH1F("hJEC_eta10_15","JEC applied in the forest, 1.0<|#eta|<1.5, ",100, 0, 5);
  TH2F * hJEC_vs_rawpT_eta10_15 = new TH2F("hJEC_vs_rawpT_eta10_15","JEC applied in the forest vs raw pT, 1.0<|#eta|<1.5",200, 0, 200, 300, 0.5, 2);

  TH1F * hJEC_eta15_20 = new TH1F("hJEC_eta15_20","JEC applied in the forest, 1.5<|#eta|<2.0, ",100, 0, 5);
  TH2F * hJEC_vs_rawpT_eta15_20 = new TH2F("hJEC_vs_rawpT_eta15_20","JEC applied in the forest vs raw pT, 1.5<|#eta|<2.0",200, 0, 200, 100, 0.5, 2);

  TH1F * hJEC_eta20_25 = new TH1F("hJEC_eta20_25","JEC applied in the forest, 2.0<|#eta|<2.5, ",100, 0, 5);
  TH2F * hJEC_vs_rawpT_eta20_25 = new TH2F("hJEC_vs_rawpT_eta20_25","JEC applied in the forest vs raw pT, 2.0<|#eta|<2.5",200, 0, 200, 100, 0.5, 2);

  TH1F * hJEC_eta25_30 = new TH1F("hJEC_eta25_30","JEC applied in the forest, 2.5<|#eta|<3.0, ",100, 0, 5);
  TH2F * hJEC_vs_rawpT_eta25_30 = new TH2F("hJEC_vs_rawpT_eta25_30","JEC applied in the forest vs raw pT, 2.5<|#eta|<3.0",200, 0, 200, 100, 0.5, 2);

  TH1F * hJEC_eta30_35 = new TH1F("hJEC_eta30_35","JEC applied in the forest, 3.0<|#eta|<3.5, ",100, 0, 5);
  TH2F * hJEC_vs_rawpT_eta30_35 = new TH2F("hJEC_vs_rawpT_eta30_35","JEC applied in the forest vs raw pT, 3.0<|#eta|<3.5",200, 0, 200, 100, 0.5, 2);

  TH1F * hJEC_eta35_40 = new TH1F("hJEC_eta35_40","JEC applied in the forest, 3.5<|#eta|<4.0, ",100, 0, 5);
  TH2F * hJEC_vs_rawpT_eta35_40 = new TH2F("hJEC_vs_rawpT_eta35_40","JEC applied in the forest vs raw pT, 3.5<|#eta|<4.0",200, 0, 200, 100, 0.5, 2);
  
  Float_t dphi = 0;
  Float_t eta_cut_min = -3;
  Float_t eta_cut_max = 3;
  
  TH1F * hAj = new TH1F("hAj","A_{j};A_{j} = #frac{p_{T}^{Lead} - p_{T}^{subLead}}{p_{T}^{Lead} + p_{T}^{subLead}};counts",300, 0, 2);
  TH1F * hVz = new TH1F("hVz","Primary Vertex Z position;v_{z};counts",200, -20, 20);
  
  TH1F * hMBSpectra = new TH1F("hMBSpectra","MB Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
  TH1F * hJet40andMB = new TH1F("hJet40andMB","MB and Jet 40 Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
  TH1F * hJet60andMB = new TH1F("hJet60andMB","MB and Jet 60 Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
  TH1F * hJet80andMB = new TH1F("hJet80andMB","MB and Jet 80 Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
  TH1F * hJet100andMB = new TH1F("hJet100andMB","MB and Jet 100 Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);

  TH1F * hJet40 = new TH1F("hJet40","Jet 40 Leading Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
  TH1F * hJet60 = new TH1F("hJet60","Jet 60 Leading Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
  TH1F * hJet80 = new TH1F("hJet80","Jet 80 Leading Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
  TH1F * hJet100 = new TH1F("hJet100","Jet 100 Leading Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);  

  TH1F *hJet40All = new TH1F("hJet40All","",50,0,300);
  TH1F *hJet60All = new TH1F("hJet60All","",50,0,300);
  TH1F *hJet80All = new TH1F("hJet80All","",50,0,300);
  TH1F *hJet100All = new TH1F("hJet100All","",50,0,300);

  TH1F *discrSSVHE = new TH1F("discrSSVHE","",30,0,6);
  TH1F *discrSSVHP = new TH1F("discrSSVHP","",30,0,6);
  TH1F *discrCSV = new TH1F("discrCSV","",30,0,1);
  
  TH1F * pt2overpt1 = new TH1F("pt2overpt1","pt2/pt1",100, 0, 2);
  TH1F * hJetEta = new TH1F("hJetEta","",60, -5, +5);
  TH1F * hJetPhi = new TH1F("hJetPhi","",60, -5, +5);
  TH1F * hJetpT = new TH1F("hJetpT","",400, 0, 600);
  TH1F * hDeltaPhi = new TH1F("hDeltaPhi","delta phi leading and subleading jets",350, 0, +3);
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jtTree[0]->GetEntries();
  if(printDebug) nentries = 500;
  TRandom rnd;

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    jtTree[0]->GetEntry(nEvt);
    jtTree[1]->GetEntry(nEvt);
    jtTree[2]->GetEntry(nEvt);
    //jtTree[4]->GetEntry(nEvt);
    jtTree[3]->GetEntry(nEvt);

    if(run == "Data"){
      if(skipPho50 && photon50_F) continue;
      if(pcollisionEventSelection_F==0) continue;
      //if(pHBHENoiseFilter_F == 0) continue;
      if(fabs(vz_F)>15) continue;
    }
    hRunN_vs_NJets->Fill(run_F, nref_F);
    
    if(jetMB_F) hMBSpectra->Fill(pt_F[0]);
    if(jetMB_F && jet40_F) hJet40andMB->Fill(pt_F[0]);
    if(jetMB_F && jet60_F) hJet60andMB->Fill(pt_F[0]);
    if(jetMB_F && jet80_F) hJet80andMB->Fill(pt_F[0]);
    if(jetMB_F && jet100_F) hJet100andMB->Fill(pt_F[0]);

    if(jet40_F) hJet40->Fill(pt_F[0]);
    if(jet60_F) hJet60->Fill(pt_F[0]);
    if(jet80_F) hJet80->Fill(pt_F[0]);

    if(jet100_F) hJet100->Fill(pt_F[0]);

    if(nref_F >=3 && doDijetImbalance) {

      j_array_pt.clear();
      j_array_phi.clear();
      j_array_eta.clear();

      Float_t dijetbalanceparameter = 0;
      Float_t referencept = 0;
      Float_t probept = 0;
      Float_t probephi = 0;
      Float_t referencephi = 0;
      Float_t refeta = 0;
      Float_t probeeta = 0;
      
      for (int  n = 0; n < nref_F; n++) { //APPLYING CUTS AND FILLS NEW ARRAY 
	if (eta_F[n] > eta_cut_min && eta_F[n] < eta_cut_max) {
	  j_array_pt.push_back(pt_F[n]);
	  j_array_phi.push_back(phi_F[n]);
	  j_array_eta.push_back(eta_F[n]);
	}
      }
            
      for (int ii = 0; ii < j_array_pt.size(); ii++)   {
	if (j_array_eta[ii] > -1.3 && j_array_eta[ii] < 1.3)  {
	  if (j_array_pt[ii] > referencept) {
	    referencept = j_array_pt[ii];
	    referencephi = j_array_phi[ii];
	    refeta = j_array_eta[ii];
	  }
	  if (j_array_pt[0]!=referencept) {
	    probept = j_array_pt[0];
	    probephi = j_array_phi[0];
	    probeeta = j_array_eta[0];
	  }
	  else {
	    probept = j_array_pt[1];
	    probephi = j_array_phi[1];
	    probeeta = j_array_eta[1];
	  }
	}
      }// pt array 

      if(printDebug) cout<<"Probe pT       = "<<probept<<endl;
      if(printDebug) cout<<"reference pT   = "<<referencept<<endl;
      if(printDebug) cout<<"leading jet pT = "<<j_array_pt[0]<<endl;
      if(printDebug) cout<<"SubLead jet pT = "<<j_array_pt[1]<<endl;
      
      Float_t alpha = 2*j_array_pt[2]/(j_array_pt[0] + j_array_pt[1]);
      if(printDebug) cout<<"alpha = "<<alpha<<endl;
      if (alpha > 0.2) continue;
      if(referencept == probept && probept != 0)  cout<< "There is a problem!    "<<probept<<endl;
      if (referencept == 0 || probept == 0) continue;     
      Float_t averagept = (float)(probept + referencept)/2;
      float delPhi = deltaphi(probephi, referencephi);
      if (delPhi < 2.5) continue;
      dijetbalanceparameter = 2*(probept - referencept)/(probept + referencept);

      if(printDebug) cout<<"dijet imbalance parameter = "<<dijetbalanceparameter<<endl;
      
      int avgpTbin = -1;
      for(int npt = 0; npt<nbins_pt; ++npt){
	if(averagept > ptbins[npt]) avgpTbin = npt;
      }
      if(avgpTbin !=-1){

	hRelResponse[avgpTbin]->Fill(dijetbalanceparameter);
    
	//FILL OUTERETA HISTOGRAMS
	if (probeeta < -1.3 || probeeta > 1.3){
	  hRelResponse_outereta[avgpTbin]->Fill(dijetbalanceparameter);
	}

	if (refeta > -1.3 && refeta < 1.3 && probeeta > -1.3 && probeeta < 1.3) {
	  Int_t v1 = rand();
	  Int_t v2 = v1%2;
	  if (v2 == 1) {
	    dijetbalanceparameter = 2*(probept - referencept)/(probept + referencept);
	  }
	  else  {
	    dijetbalanceparameter = 2*(referencept - probept)/(probept + referencept);
	  }
	  //FILL INNERETA HISTOGRAMS
	  hRelResponse_innereta[avgpTbin]->Fill(dijetbalanceparameter);
	}//refeta, probeeta if statement
      
      }// avg pT bin

      j_array_pt.clear();
      j_array_phi.clear();
      j_array_eta.clear();

    }//dijet imbalance
    
    for(int ijet=0; ijet<nref_F; ijet++){
      hJet40All->Fill(pt_F[ijet]);
      hJet60All->Fill(pt_F[ijet]);
      hJet80All->Fill(pt_F[ijet]);
      hJet100All->Fill(pt_F[ijet]);

      if(doBjets){
	discrSSVHE->Fill(discr_ssvHighEff_F[ijet]);
	discrSSVHP->Fill(discr_ssvHighPur_F[ijet]);
	discrCSV->Fill(discr_csv_F[ijet]);
      }

      if(nref_F >= 2) hDeltaPhi->Fill(deltaphi(phi_F[0], phi_F[1]));
      
      if(jet80_F && pt_F[0] >90.0 && pt_F[1] >20.0) { 
	hJetEta->Fill(eta_F[ijet]);
	hJetPhi->Fill(phi_F[ijet]);
	hJetpT->Fill(pt_F[ijet]);
      }
      
      if(fabs(eta_F[ijet]) < 0.5){
	hJEC_eta0_05->Fill((float)pt_F[ijet]/rawpt_F[ijet]);
	hJEC_vs_rawpT_eta0_05->Fill(rawpt_F[ijet],(float)pt_F[ijet]/rawpt_F[ijet]);
      }
      if(fabs(eta_F[ijet]) > 0.5 && fabs(eta_F[ijet])<1.0){
	hJEC_eta05_10->Fill((float)pt_F[ijet]/rawpt_F[ijet]);
	hJEC_vs_rawpT_eta05_10->Fill(rawpt_F[ijet],(float)pt_F[ijet]/rawpt_F[ijet]);
      }
      if(fabs(eta_F[ijet]) > 1.0 && fabs(eta_F[ijet])<1.5){
	hJEC_eta10_15->Fill((float)pt_F[ijet]/rawpt_F[ijet]);
	hJEC_vs_rawpT_eta10_15->Fill(rawpt_F[ijet],(float)pt_F[ijet]/rawpt_F[ijet]);
      }
      if(fabs(eta_F[ijet]) > 1.5 && fabs(eta_F[ijet])<2){
	hJEC_eta15_20->Fill((float)pt_F[ijet]/rawpt_F[ijet]);
	hJEC_vs_rawpT_eta15_20->Fill(rawpt_F[ijet],(float)pt_F[ijet]/rawpt_F[ijet]);
      }
      if(fabs(eta_F[ijet]) > 2 && fabs(eta_F[ijet])<2.5){
	hJEC_eta20_25->Fill((float)pt_F[ijet]/rawpt_F[ijet]);
	hJEC_vs_rawpT_eta20_25->Fill(rawpt_F[ijet],(float)pt_F[ijet]/rawpt_F[ijet]);
      }
      if(fabs(eta_F[ijet]) > 2.5 && fabs(eta_F[ijet])<3){
	hJEC_eta25_30->Fill((float)pt_F[ijet]/rawpt_F[ijet]);
	hJEC_vs_rawpT_eta25_30->Fill(rawpt_F[ijet],(float)pt_F[ijet]/rawpt_F[ijet]);
      }
      if(fabs(eta_F[ijet]) > 3 && fabs(eta_F[ijet])<3.5){
	hJEC_eta30_35->Fill((float)pt_F[ijet]/rawpt_F[ijet]);
	hJEC_vs_rawpT_eta30_35->Fill(rawpt_F[ijet],(float)pt_F[ijet]/rawpt_F[ijet]);
      }
      if(fabs(eta_F[ijet]) > 3.5 && fabs(eta_F[ijet])<4){
	hJEC_eta35_40->Fill((float)pt_F[ijet]/rawpt_F[ijet]);
	hJEC_vs_rawpT_eta35_40->Fill(rawpt_F[ijet],(float)pt_F[ijet]/rawpt_F[ijet]);
      }
    }// njets

    if(run == "MC"){
      if(pt_F[0]>90.0 && pt_F[1] >20.0){
	float Aj = (float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]);
	pt2overpt1->Fill((float)pt_F[1]/pt_F[0]);
	hAj->Fill(Aj);
      }
    }//mc

    if(run == "Data"){
      if(jet80_F && pt_F[0]>90.0 && pt_F[1] >20.0){
	float Aj = (float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]);
	pt2overpt1->Fill((float)pt_F[1]/pt_F[0]);
	hAj->Fill(Aj);
      }
      // int binpt = -1, bineta = -1;
      // for(int npt = 0; npt<nbins_pt; ++npt){
      //   if((pt_F[0]+pt_F[1]) > ptbins[npt]) binpt = npt;
      // }
      // if(binpt == -1) continue;
      // for(int neta = 0; neta<nbins_eta; ++neta){
      //   if(eta_F[0] > etabins[neta]) bineta = neta;      
      // }
      // if(bineta == -1) continue;

      // if((float)(pt_F[2]/ptAvg) < 0.2 && deltaphi(phi_F[0], phi_F[1]) > (float)2*pi/3)
      //   hRelResponse[binpt][bineta]->Fill(DijetRel);
      // hAj[binpt][bineta]->Fill(Aj);

    } //data

  }// nevents
  // TH1F * hJet40Turnon = (TH1F*)hJet40andMB->Clone("hJet40Turnon");
  // hJet40Turnon->Divide(hMBSpectra);
  // TH1F * hJet60Turnon = (TH1F*)hJet60andMB->Clone("hJet60Turnon");
  // hJet60Turnon->Divide(hMBSpectra);
  // TH1F * hJet80Turnon = (TH1F*)hJet80andMB->Clone("hJet80Turnon");
  // hJet80Turnon->Divide(hMBSpectra);
  // TH1F * hJet100Turnon = (TH1F*)hJet100andMB->Clone("hJet100Turnon");
  // hJet100Turnon->Divide(hMBSpectra);
  
  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;  
  
}
