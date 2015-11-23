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
		   std::string jetType= "PF",
		   std::string algo= "",
		   std::string kFoname="PromptForest")
{

  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool doBjets = true;
  std::string bJetString = "";
  if(doBjets) bJetString = "_withBjets";
  bool skipPho30 = true;
  bool printDebug = true;
  bool doDijetImbalance = true;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  infile_Forest = Form("%s_%s_ExpressForest.txt", coll.c_str(), run.c_str());
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  
  const int N = 6; // data does not have the runAnalyzer

  TChain * jtTree[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = "test";
  if(coll == "PbPb") dir[2] = Form("ak%s%d%sJetAnalyzer", algo.c_str(), radius, jetType.c_str());
  if(coll == "PP") dir[2] = Form("ak%d%sJetAnalyzer", radius, jetType.c_str());
  //dir[3] = "akPu3CaloJetAnalyzer";
  dir[3] = "hiEvtAnalyzer";
  if(jetType == "Calo") dir[4] = "rechitanalyzer" ;
  if(jetType == "PF") dir[4] = "pfcandAnalyzer" ;
  if(run == "MC") dir[5] = "runAnalyzer";
  // dir[4] = "hltobject";

  string trees[N];
  trees[0] = "HltTree";
  trees[1] = "HltTree";
  trees[2] = "t";
  // trees[3] = "t";
  trees[3] = "HiTree";
    // , "jetObjTree"
  if(jetType == "Calo") trees[4] = "tower" ;
  if(jetType == "PF") trees[4] = "pfTree" ;
  if(run == "MC") trees[5] = "run";

  int NLoop = 5;
  if(run == "MC") NLoop = 6;
  
  for(int t = 0;t<NLoop;t++){
    jtTree[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends
  
  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    jtTree[0]->Add(filename_Forest.c_str());
    jtTree[1]->Add(filename_Forest.c_str());
    jtTree[2]->Add(filename_Forest.c_str());
    jtTree[3]->Add(filename_Forest.c_str());
    jtTree[4]->Add(filename_Forest.c_str());
    if(run == "MC") jtTree[5]->Add(filename_Forest.c_str());
    
    cout<<"filename: "<<filename_Forest<<endl;
    
    if(printDebug)cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[0]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[1]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[2]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[3]->GetEntries() << endl;
    if(printDebug)cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
    if(printDebug)cout << "Entries : " << jtTree[4]->GetEntries() << endl;

    if(run == "MC") {
      if(printDebug)cout << "Tree loaded  " << string(dir[5]+"/"+trees[5]).data() << endl;
      if(printDebug)cout << "Entries : " << jtTree[5]->GetEntries() << endl;
    }
    
    cout<<"Total number of events loaded in HiForest = "<<jtTree[2]->GetEntries()<<endl;
  }
  
  jtTree[2]->AddFriend(jtTree[0]);
  jtTree[2]->AddFriend(jtTree[1]);
  jtTree[2]->AddFriend(jtTree[3]);
  jtTree[2]->AddFriend(jtTree[4]);
  if(run == "MC") jtTree[2]->AddFriend(jtTree[5]);
  
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
  int photon30_F;
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
  float crossSection_F;
  int ssvTrg60_F;
  int ssvTrg80_F;
  int csvTrg60_F;
  int csvTrg80_F;

  Int_t nPFpart_F;
  Int_t pfId_F[NOBJECT_MAX];
  Float_t pfPt_F[NOBJECT_MAX];
  Float_t pfVsPtInitial_F[NOBJECT_MAX];
  Float_t pfEta_F[NOBJECT_MAX];
  Float_t pfPhi_F[NOBJECT_MAX];


  if(jetType == "Calo") {
    jtTree[4]->SetBranchAddress("n", &nPFpart_F);
    jtTree[4]->SetBranchAddress("et", pfPt_F);
    jtTree[4]->SetBranchAddress("eta", pfEta_F);
    jtTree[4]->SetBranchAddress("phi", pfPhi_F);
  }
  if(jetType == "PF") {
    jtTree[4]->SetBranchAddress("nPFpart", &nPFpart_F);
    jtTree[4]->SetBranchAddress("pfId", pfId_F);
    jtTree[4]->SetBranchAddress("pfPt", pfPt_F);
    jtTree[4]->SetBranchAddress("pfVsPtInitial", pfVsPtInitial_F);
    jtTree[4]->SetBranchAddress("pfEta", pfEta_F);
    jtTree[4]->SetBranchAddress("pfPhi", pfPhi_F);
  }

  if(run == "MC"){
    jtTree[5]->SetBranchAddress("xsec",&crossSection_F);
  }

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
  jtTree[1]->SetBranchAddress("pHBHENoiseFilterResultProducer",&pHBHENoiseFilter_F);
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
  jtTree[0]->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5_v1",&photon30_F);
  jtTree[0]->SetBranchAddress("HLT_AK4PFBJetBSSV60_Eta2p1_v1",&ssvTrg60_F);
  jtTree[0]->SetBranchAddress("HLT_AK4PFBJetBSSV80_Eta2p1_v1",&ssvTrg80_F);
  jtTree[0]->SetBranchAddress("HLT_AK4PFBJetBCSV60_Eta2p1_v1",&csvTrg60_F);
  jtTree[0]->SetBranchAddress("HLT_AK4PFBJetBCSV80_Eta2p1_v1",&csvTrg80_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);

  std::string rad = Form("%d",radius);
  std::string end = Form("%d",endfile);
  
  TFile *fout = new TFile((kFoname+coll+"_"+run+"_ak"+algo+rad+jetType+bJetString+"_"+end+".root").c_str(),"RECREATE");
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

  vector<Float_t> j_array_pt;
  vector<Float_t> j_array_phi;
  vector<Float_t> j_array_eta;

  
  Float_t dphi = 0;
  Float_t eta_cut_min = -3;
  Float_t eta_cut_max = 3;

  TH2F * hRunN_vs_NJets = new TH2F("hRunN_Vs_NJets","",719, 261445, 262164, 50, 0, 50);
  TH1F * hVz = new TH1F("hVz","Primary Vertex Z position;v_{z};counts",200, -20, 20);

  TH1F *hRelResponse[nbins_pt][ncen+1];
  TH1F *hRelResponse_outereta[nbins_pt][ncen+1];
  TH1F *hRelResponse_innereta[nbins_pt][ncen+1];

  TH3F * hJEC[ncen+1];

  TH1F * hAj[ncen+1];
  TH1F * pt2overpt1[ncen+1];
  TH1F * hDeltaPhi[ncen+1];

  TH1F * hMBSpectra[ncen+1];
  TH1F * hJet40andMB[ncen+1];
  TH1F * hJet60andMB[ncen+1];
  TH1F * hJet80andMB[ncen+1];
  TH1F * hJet100andMB[ncen+1];
  TH1F * hJet40[ncen+1];
  TH1F * hJet60[ncen+1];
  TH1F * hJet80[ncen+1];
  TH1F * hJet100[ncen+1];
  TH1F * hJet40All[ncen+1];
  TH1F * hJet60All[ncen+1];
  TH1F * hJet80All[ncen+1];
  TH1F * hJet100All[ncen+1];
  TH1F * hJetEta[ncen+1];
  TH1F * hJetPhi[ncen+1];
  TH1F * hJetpT[ncen+1];
    
  TH2F * hPFCand_eta_vs_pT[PFType][ncen+1];
  TH1F * hPFCand_pTscale_insideJet[PFType][ncen+1];

  TH1F * discrSSVHE[ncen+1];
  TH1F * discrSSVHP[ncen+1];
  TH1F * discrCSV[ncen+1];
  TH1F * csvTrg80[ncen+1];
  TH1F * ssvTrg80[ncen+1];
  TH1F * csvTrg60[ncen+1];
  TH1F * ssvTrg60[ncen+1];
  TH1F * csvTrg80withCSV[ncen+1];
  TH1F * csvTrg60withCSV[ncen+1];
  TH1F * ssvTrg80withSSVHP[ncen+1];
  TH1F * ssvTrg60withSSVHP[ncen+1];

  TH1F * csvDistr[ncen+1];
  TH1F * ssvDistr[ncen+1];

  
  for(int icen = 0; icen<=ncen; ++icen){
    for(int npt = 0; npt<nbins_pt; ++npt){
      hRelResponse[npt][icen] = new TH1F(Form("hRelResponse_ptbin%d_%s",npt, cdir[icen]),Form("Relative Response in %d < p_{T}^{avg} < %d", ptbins[npt], ptbins[npt+1]), 1000, -3, 3);
      hRelResponse_outereta[npt][icen] = new TH1F(Form("hRelResponse_outereta_ptbin%d_%s",npt, cdir[icen]),Form("Relative Response in %d < p_{T}^{avg} < %d", ptbins[npt], ptbins[npt+1]), 1000, -3, 3);
      hRelResponse_innereta[npt][icen] = new TH1F(Form("hRelResponse_innereta_ptbin%d_%s",npt, cdir[icen]),Form("Relative Response in %d < p_{T}^{avg} < %d", ptbins[npt], ptbins[npt+1]), 1000, -3, 3);
    }

    hJEC[icen] = new TH3F(Form("hJEC_%s",cdir[icen]),";raw p_{T};#eta;JEC",500, 0, 500, 200, -5, +5, 300, 0, 5);

    hAj[icen] = new TH1F(Form("hAj_%s", cdir[icen]),"A_{j};A_{j} = #frac{p_{T}^{Lead} - p_{T}^{subLead}}{p_{T}^{Lead} + p_{T}^{subLead}};counts",300, 0, 2);
  
    hMBSpectra[icen] = new TH1F(Form("hMBSpectra_%s", cdir[icen]),"MB Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
    hJet40andMB[icen] = new TH1F(Form("hJet40andMB_%s", cdir[icen]),"MB and Jet 40 Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
    hJet60andMB[icen] = new TH1F(Form("hJet60andMB_%s", cdir[icen]),"MB and Jet 60 Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
    hJet80andMB[icen] = new TH1F(Form("hJet80andMB_%s", cdir[icen]),"MB and Jet 80 Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
    hJet100andMB[icen] = new TH1F(Form("hJet100andMB_%s", cdir[icen]),"MB and Jet 100 Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);

    hJet40[icen] = new TH1F(Form("hJet40_%s", cdir[icen]),"Jet 40 Leading Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
    hJet60[icen] = new TH1F(Form("hJet60_%s", cdir[icen]),"Jet 60 Leading Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
    hJet80[icen] = new TH1F(Form("hJet80_%s", cdir[icen]),"Jet 80 Leading Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);
    hJet100[icen] = new TH1F(Form("hJet100_%s", cdir[icen]),"Jet 100 Leading Spectra;Jet p_{T} GeV/c;counts",100, 0, 500);  

    hJet40All[icen] = new TH1F(Form("hJet40All_%s", cdir[icen]),"",50,0,300);
    hJet60All[icen] = new TH1F(Form("hJet60All_%s", cdir[icen]),"",50,0,300);
    hJet80All[icen] = new TH1F(Form("hJet80All_%s", cdir[icen]),"",50,0,300);
    hJet100All[icen] = new TH1F(Form("hJet100All_%s", cdir[icen]),"",50,0,300);

    discrSSVHE[icen] = new TH1F(Form("discrSSVHE_%s", cdir[icen]),"",30,0,6);
    discrSSVHP[icen] = new TH1F(Form("discrSSVHP_%s", cdir[icen]),"",30,0,6);
    discrCSV[icen] = new TH1F(Form("discrCSV_%s", cdir[icen]),"",30,0,1);
  
    csvTrg80[icen] = new TH1F(Form("csvTrg80_%s", cdir[icen]),"",60,0,300);
    ssvTrg80[icen] = new TH1F(Form("ssvTrg80_%s", cdir[icen]),"",60,0,300);
    csvTrg60[icen] = new TH1F(Form("csvTrg60_%s", cdir[icen]),"",60,0,300);
    ssvTrg60[icen] = new TH1F(Form("ssvTrg60_%s", cdir[icen]),"",60,0,300);

    csvTrg80withCSV[icen] = new TH1F(Form("csvTrg80withCSV_%s", cdir[icen]),"",60,0,300);
    csvTrg60withCSV[icen] = new TH1F(Form("csvTrg60withCSV_%s", cdir[icen]),"",60,0,300);
    ssvTrg80withSSVHP[icen] = new TH1F(Form("ssvTrg80withSSVHP_%s", cdir[icen]),"",60,0,300);
    ssvTrg60withSSVHP[icen] = new TH1F(Form("ssvTrg60withSSVHP_%s", cdir[icen]),"",60,0,300);

    csvDistr[icen] = new TH1F(Form("csvDistr_%s", cdir[icen]),"",60,0,300);
    ssvDistr[icen] = new TH1F(Form("ssvDistr_%s", cdir[icen]),"",60,0,300);

    pt2overpt1[icen] = new TH1F(Form("pt2overpt1_%s", cdir[icen]),"pt2/pt1",100, 0, 2);
    hJetEta[icen] = new TH1F(Form("hJetEta_%s", cdir[icen]),"",60, -5, +5);
    hJetPhi[icen] = new TH1F(Form("hJetPhi_%s", cdir[icen]),"",60, -5, +5);
    hJetpT[icen] = new TH1F(Form("hJetpT_%s", cdir[icen]),"",400, 0, 600);
    hDeltaPhi[icen] = new TH1F(Form("hDeltaPhi_%s", cdir[icen]),"delta phi leading and subleading jets",350, 0, +3);
    
    if(jetType == "PF"){
      for(int ipf = 0; ipf<PFType; ++ipf){
	hPFCand_eta_vs_pT[ipf][icen] = new TH2F(Form("hPFCand_eta_vs_pT_%s_%s",PFCandType[ipf].c_str(), cdir[icen]),Form("%s;#eta;candidate p_{T}",PFCandType[ipf].c_str()),120, -5, +5, 300, 0, 300);
	hPFCand_pTscale_insideJet[ipf][icen] = new TH1F(Form("hPFCand_pTscale_insideJet_%s_%s",PFCandType[ipf].c_str(), cdir[icen]),Form("%s;candidate p_{T}/jet p_{T};counts",PFCandType[ipf].c_str()),100, 0, 1);
      }
    }
    if(jetType == "Calo"){
      hPFCand_eta_vs_pT[0][icen] = new TH2F(Form("hPFCand_eta_vs_pT_CaloCand_%s",cdir[icen]),"Calo Cand ;#eta;candidate p_{T}",120, -5, +5, 300, 0, 300);
      hPFCand_pTscale_insideJet[0][icen] = new TH1F(Form("hPFCand_pTscale_insideJet_CaloCand_%s", cdir[icen]),"Calo Cand ;candidate p_{T}/jet p_{T};counts",100, 0, 1);
    }
    
    
  }// cent loop
    
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
    jtTree[3]->GetEntry(nEvt);
    jtTree[4]->GetEntry(nEvt);
    if(run == "MC") jtTree[5]->GetEntry(nEvt);
    
    if(run == "Data"){
      if(skipPho30 && photon30_F) continue;
      if(pcollisionEventSelection_F==0) continue;
      if(pHBHENoiseFilter_F == 0) continue;
      if(fabs(vz_F)>15) continue;
    }
    hRunN_vs_NJets->Fill(run_F, nref_F);
    hVz->Fill(vz_F);

    // find the centrality bin.

    int centbin = -1;
    if(coll == "PP") centbin = ncen;
    if(coll == "PbPb") centbin = findBin(hiBin_F);
    if(centbin == -1) continue;
    
    for(int npf = 0; npf<nPFpart_F; ++npf){

      if(jetType == "Calo"){
	hPFCand_eta_vs_pT[0][centbin]->Fill(pfEta_F[npf], pfPt_F[npf]);
	for(int njet = 0; njet<nref_F; ++njet){
	  float delR = deltaR(pfEta_F[npf], pfPhi_F[npf], eta_F[njet], phi_F[njet]);
	  if(delR < (float)radius/10)
	    hPFCand_pTscale_insideJet[0][centbin]->Fill((float)(pfPt_F[npf]/pt_F[njet]));	
	}
      }

      if(jetType == "PF"){
	hPFCand_eta_vs_pT[pfId_F[npf]][centbin]->Fill(pfEta_F[npf], pfPt_F[npf]);
	for(int njet = 0; njet<nref_F; ++njet){
	  float delR = deltaR(pfEta_F[npf], pfPhi_F[npf], eta_F[njet], phi_F[njet]);
	  if(delR < (float)radius/10)
	    hPFCand_pTscale_insideJet[pfId_F[npf]][centbin]->Fill((float)(pfPt_F[npf]/pt_F[njet])); 
	}
      }
      
    }// pf cands loop
    
    
    if(jetMB_F) hMBSpectra[centbin]->Fill(pt_F[0]);
    if(jetMB_F && jet40_F) hJet40andMB[centbin]->Fill(pt_F[0]);
    if(jetMB_F && jet60_F) hJet60andMB[centbin]->Fill(pt_F[0]);
    if(jetMB_F && jet80_F) hJet80andMB[centbin]->Fill(pt_F[0]);
    if(jetMB_F && jet100_F) hJet100andMB[centbin]->Fill(pt_F[0]);

    if(jet40_F) hJet40[centbin]->Fill(pt_F[0]);
    if(jet60_F) hJet60[centbin]->Fill(pt_F[0]);
    if(jet80_F) hJet80[centbin]->Fill(pt_F[0]);
    if(jet100_F) hJet100[centbin]->Fill(pt_F[0]);
    
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
            
      for (unsigned ii = 0; ii < j_array_pt.size(); ii++)   {
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
      if (alpha < 0.2){
	if(referencept == probept && probept != 0)  cout<< "There is a problem!    "<<probept<<endl;
	if (referencept != 0 && probept != 0) {     
	  Float_t averagept = (float)(probept + referencept)/2;
	  float delPhi = deltaphi(probephi, referencephi);
	  if (delPhi > (float)(2*pi/3)) {
	    dijetbalanceparameter = 2*(probept - referencept)/(probept + referencept);
	    
	    if(printDebug) cout<<"dijet imbalance parameter = "<<dijetbalanceparameter<<endl;
      
	    int avgpTbin = -1;
	    for(int npt = 0; npt<nbins_pt; ++npt){
	      if(averagept > ptbins[npt]) avgpTbin = npt;
	    }
	    if(avgpTbin !=-1){
	      
	      hRelResponse[avgpTbin][centbin]->Fill(dijetbalanceparameter);
	      
	      //FILL OUTERETA HISTOGRAMS
	      if (probeeta < -1.3 || probeeta > 1.3){
		hRelResponse_outereta[avgpTbin][centbin]->Fill(dijetbalanceparameter);
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
		hRelResponse_innereta[avgpTbin][centbin]->Fill(dijetbalanceparameter);
	      }//refeta, probeeta if statement
	      
	    }// avg pT bin

	  }// delta phi selection

	}// refpt!=-, probpt!=0

      }// alpha selection, suppression of third jet
      j_array_pt.clear();
      j_array_phi.clear();
      j_array_eta.clear();
      
    }
    
    for(int ijet=0; ijet<nref_F; ijet++){

      if(jet40_F) hJet40All[centbin]->Fill(pt_F[ijet]);
      if(jet60_F) hJet60All[centbin]->Fill(pt_F[ijet]);
      if(jet80_F) hJet80All[centbin]->Fill(pt_F[ijet]);
      if(jet100_F) hJet100All[centbin]->Fill(pt_F[ijet]);

      if(doBjets){
	discrSSVHE[centbin]->Fill(discr_ssvHighEff_F[ijet]);
	discrSSVHP[centbin]->Fill(discr_ssvHighPur_F[ijet]);
	discrCSV[centbin]->Fill(discr_csv_F[ijet]);

	if(csvTrg60_F){
		csvTrg60[centbin]->Fill(pt_F[ijet]);
		if(discr_csv_F[ijet]>0.9) csvTrg60withCSV[centbin]->Fill(pt_F[ijet]);
	}
	if(csvTrg80_F){
		csvTrg80[centbin]->Fill(pt_F[ijet]);
		if(discr_csv_F[ijet]>0.9) csvTrg80withCSV[centbin]->Fill(pt_F[ijet]);
	}
	if(ssvTrg60_F){
		ssvTrg60[centbin]->Fill(pt_F[ijet]);
		if(discr_ssvHighPur_F[ijet]>1.2) ssvTrg60withSSVHP[centbin]->Fill(pt_F[ijet]);
	}
	if(ssvTrg80_F){
		ssvTrg80[centbin]->Fill(pt_F[ijet]);
		if(discr_ssvHighPur_F[ijet]>1.2) ssvTrg80withSSVHP[centbin]->Fill(pt_F[ijet]);
	}

	if(discr_csv_F[ijet]>0.9) csvDistr[centbin]->Fill(pt_F[ijet]);
	if(discr_ssvHighPur_F[ijet]>1.2) ssvDistr[centbin]->Fill(pt_F[ijet]);

      }

      if(nref_F >= 2) hDeltaPhi[centbin]->Fill(deltaphi(phi_F[0], phi_F[1]));
      
      if(jet80_F && pt_F[0] >LjCut && pt_F[1] >SbjCut) { 
	hJetEta[centbin]->Fill(eta_F[ijet]);
	hJetPhi[centbin]->Fill(phi_F[ijet]);
	hDeltaPhi[centbin]->Fill(deltaphi(phi_F[0], phi_F[1]));
      }

      hJetpT[centbin]->Fill(pt_F[ijet]);

      hJEC[centbin]->Fill(rawpt_F[ijet], eta_F[ijet], (float)(pt_F[ijet]/rawpt_F[ijet]));
      
    }// njets

    if(run == "MC"){
      if(pt_F[0]>LjCut && pt_F[1] >SbjCut){
	float Aj = (float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]);
	pt2overpt1[centbin]->Fill((float)pt_F[1]/pt_F[0]);
	hAj[centbin]->Fill(Aj);
      }
    }//mc

    if(run == "Data"){
      if(jet80_F && pt_F[0]>LjCut && pt_F[1] >SbjCut){
	float Aj = (float)(pt_F[0]-pt_F[1])/(pt_F[0]+pt_F[1]);
	pt2overpt1[centbin]->Fill((float)pt_F[1]/pt_F[0]);
	hAj[centbin]->Fill(Aj);
      }

    } //data

  }// nevents
  
  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;  
  
}
