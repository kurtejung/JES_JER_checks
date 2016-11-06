
#ifndef __FORESTCONTAINER__
#define __FORESTCONTAINER__

#include "boundaries.h"

using namespace std;

class ForestContainer{

public:

	ForestContainer(int startfile, int endfile, int radius, string coll, string run, string jetType, string algo){

		std::string infile_Forest;
		infile_Forest = Form("%s_%s_forests.txt", coll.c_str(), run.c_str());
		std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
		std::string filename_Forest;

		if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
        if(!instr_Forest){ cout << "Cannot find " << infile_Forest << "!" << endl; exit(0); }

		for(int ifile = 0;ifile<startfile;ifile++){
			instr_Forest>>filename_Forest;
		}

        string dir[N];
        dir[0] = "hltanalysis";
        dir[1] = "skimanalysis";
        dir[2] = "";
        if(coll == "PbPb") dir[2] = Form("ak%s%d%sJetAnalyzer", algo.c_str(), radius, jetType.c_str());
        if(coll == "pPb") dir[2] = Form("ak%s%d%sJetAnalyzer", algo.c_str(), radius, jetType.c_str());
        if(coll == "PP") dir[2] = Form("ak%d%sJetAnalyzer", radius, jetType.c_str());
        dir[3] = "hiEvtAnalyzer";
        if(jetType == "Calo") dir[4] = "rechitanalyzer" ;
        if(jetType == "PF") dir[4] = "pfcandAnalyzer" ;
        if(run == "MC") dir[5] = "runAnalyzer";
        if(coll != "PP") dir[6] = "hiFJRhoAnalyzer";

        string trees[N];
        trees[0] = "HltTree";
        trees[1] = "HltTree";
        trees[2] = "t";
  // trees[3] = "t";
        trees[3] = "HiTree";
        if(jetType == "Calo") trees[4] = "tower" ;
        if(jetType == "PF") trees[4] = "pfTree" ;
  //trees[5] = "jetObjTree";
        if(run == "MC") trees[5] = "run";
        if(coll != "PP") trees[6] = "t";
        

        int NLoop = N;

        for(int t = 0;t<NLoop;t++){
            if(t==5 && run!= "MC") continue;
        	jtTree[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
        }//tree loop ends

        for(int ifile = startfile; ifile<endfile; ++ifile){

        	instr_Forest>>filename_Forest;

        	jtTree[0]->Add(filename_Forest.c_str());
        	jtTree[1]->Add(filename_Forest.c_str());
        	jtTree[2]->Add(filename_Forest.c_str());
        	jtTree[3]->Add(filename_Forest.c_str());
        	jtTree[4]->Add(filename_Forest.c_str());
    //jtTree[5]->Add(filename_Forest.c_str());
        	if(run == "MC") jtTree[5]->Add(filename_Forest.c_str());

            jtTree[6]->Add(filename_Forest.c_str());

        	cout<<"filename: "<<filename_Forest<<endl;

        	if(printDebug){
                cout << "Tree loaded  " << string(dir[0]+"/"+trees[0]).data() << endl;
                cout << "Entries : " << jtTree[0]->GetEntries() << endl;
                cout << "Tree loaded  " << string(dir[1]+"/"+trees[1]).data() << endl;
                cout << "Entries : " << jtTree[1]->GetEntries() << endl;
                cout << "Tree loaded  " << string(dir[2]+"/"+trees[2]).data() << endl;
                cout << "Entries : " << jtTree[2]->GetEntries() << endl;
                cout << "Tree loaded  " << string(dir[3]+"/"+trees[3]).data() << endl;
                cout << "Entries : " << jtTree[3]->GetEntries() << endl;
                cout << "Tree loaded  " << string(dir[4]+"/"+trees[4]).data() << endl;
                cout << "Entries : " << jtTree[4]->GetEntries() << endl;
    //if(printDebug)cout << "Tree loaded  " << string(dir[5]+"/"+trees[5]).data() << endl;
    //if(printDebug)cout << "Entries : " << jtTree[5]->GetEntries() << endl;
                
                if(run == "MC") {
                    cout << "Tree loaded  " << string(dir[5]+"/"+trees[5]).data() << endl;
                    cout << "Entries : " << jtTree[5]->GetEntries() << endl;
                }
                cout << "Tree loaded  " << string(dir[6]+"/"+trees[6]).data() << endl;
                cout << "Entries : " << jtTree[6]->GetEntries() << endl;
            }
              
            cout<<"Total number of events loaded in HiForest = "<<jtTree[2]->GetEntries()<<endl;
        }

        jtTree[2]->AddFriend(jtTree[0]);
        jtTree[2]->AddFriend(jtTree[1]);
        jtTree[2]->AddFriend(jtTree[3]);
        jtTree[2]->AddFriend(jtTree[4]);
  //jtTree[2]->AddFriend(jtTree[5]);
        if(run == "MC") jtTree[2]->AddFriend(jtTree[5]);

  //jtTree[3]->AddFriend(jtTree[0]);
  //jtTree[3]->AddFriend(jtTree[1]);
  //jtTree[3]->AddFriend(jtTree[4]);

        if(jetType == "Calo") {
        	jtTree[4]->SetBranchAddress("n", &nPFpart_F);
        	jtTree[4]->SetBranchAddress("et", pfPt_F);
        	jtTree[4]->SetBranchAddress("eta", pfEta_F);
        	jtTree[4]->SetBranchAddress("phi", pfPhi_F);
        }
        if(jetType == "PF") {
        	jtTree[4]->SetBranchAddress("nPFpart", &nPFpart_F);
        	jtTree[4]->SetBranchAddress("pfId", &pfId_F);
        	jtTree[4]->SetBranchAddress("pfPt", &pfPt_F);
        	jtTree[4]->SetBranchAddress("pfVsPtInitial", &pfVsPtInitial_F);
        	jtTree[4]->SetBranchAddress("pfEta", &pfEta_F);
        	jtTree[4]->SetBranchAddress("pfPhi", &pfPhi_F);

        }

        if(run == "MC"){
        	jtTree[5]->SetBranchAddress("xsec",&crossSection_F);
        }
        if(coll != "PP"){
            jtTree[6]->SetBranchAddress("rho",&rho);
            jtTree[6]->SetBranchAddress("rhom",&rhoM);
            jtTree[6]->SetBranchAddress("etaMax",&etaMax);
        }

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
        if(coll == "PP")jtTree[1]->SetBranchAddress("PAcollisionEventSelection",&pcollisionEventSelection_F);
        
        ///WARNING THIS IS NOT RIGHT!!!!
        if(coll=="pPb")jtTree[1]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pcollisionEventSelection_F);
        
        
        if(coll == "PbPb")jtTree[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
        jtTree[1]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);

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

        if(coll == "PP"){

        	jtTree[0]->SetBranchAddress("L1_MinimumBiasHF1_OR",&L1_MB_F);

        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part0_v1",&jetMB_p0_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part1_v1",&jetMB_p1_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part2_v1",&jetMB_p2_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part3_v1",&jetMB_p3_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part4_v1",&jetMB_p4_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1",&jetMB_p5_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part6_v1",&jetMB_p6_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part7_v1",&jetMB_p7_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part8_v1",&jetMB_p8_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part9_v1",&jetMB_p9_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part10_v1",&jetMB_p10_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part11_v1",&jetMB_p11_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part12_v1",&jetMB_p12_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part13_v1",&jetMB_p13_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part14_v1",&jetMB_p14_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part15_v1",&jetMB_p15_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part16_v1",&jetMB_p16_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part17_v1",&jetMB_p17_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part18_v1",&jetMB_p18_F);  
        	jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part19_v1",&jetMB_p19_F);  

        	jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet40_Eta5p1_v1", jetType.c_str()),&jet40_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet40_Eta5p1_v1_Prescl",jetType.c_str()),&jet40_p_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet60_Eta5p1_v1", jetType.c_str()),&jet60_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet60_Eta5p1_v1_Prescl",jetType.c_str()),&jet60_p_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet80_Eta5p1_v1", jetType.c_str()),&jet80_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet80_Eta5p1_v1_Prescl",jetType.c_str()),&jet80_p_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet100_Eta5p1_v1", jetType.c_str()),&jet100_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_AK4%sJet100_Eta5p1_v1_Prescl",jetType.c_str()),&jet100_p_F);
        	jtTree[0]->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5_v1",&photon30_F);
        	jtTree[0]->SetBranchAddress("HLT_AK4PFBJetBSSV60_Eta2p1_v1",&ssvTrg60_F);
        	jtTree[0]->SetBranchAddress("HLT_AK4PFBJetBSSV80_Eta2p1_v1",&ssvTrg80_F);
        	jtTree[0]->SetBranchAddress("HLT_AK4PFBJetBCSV60_Eta2p1_v1",&csvTrg60_F);
        	jtTree[0]->SetBranchAddress("HLT_AK4PFBJetBCSV80_Eta2p1_v1",&csvTrg80_F);
            jtTree[0]->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5_v1",&photon30_F);
                    
        }

        if(coll == "pPb"){
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1",&jetMB_p1_F);  
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1",&jetMB_p2_F);  
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1",&jetMB_p3_F);  
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1",&jetMB_p4_F);  
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1",&jetMB_p5_F);  
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1",&jetMB_p6_F);  
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1",&jetMB_p7_F);  
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1",&jetMB_p8_F);  
            jtTree[0]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_v1",&jetMB_p0_F);
            
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloJet40_Eta5p1_v2",&jet40_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloJet40_Eta5p1_v2_Prescl",&jet40_p_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v2",&jet60_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v2_Prescl",&jet60_p_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v2",&jet80_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v2_Prescl",&jet80_p_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v2",&jet100_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v2_Prescl",&jet100_p_F);
            jtTree[0]->SetBranchAddress("HLT_PASinglePhoton30_Eta3p1_v2",&photon30_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloBJetCSV60_Eta2p1_v1",&csvTrg60_F);
            jtTree[0]->SetBranchAddress("HLT_PAAK4CaloBJetCSV80_Eta2p1_v1",&csvTrg80_F);
        }
        if(coll == "PbPb"){
        	jtTree[0]->SetBranchAddress("",&jetMB_F);
    //jtTree[0]->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part1_v1_Prescl",&jetMB_p_F);  
        	jtTree[0]->SetBranchAddress(Form("HLT_HIPuAK4%sJet40_Eta5p1_v1", jetType.c_str()),&jet40_F);
    //jtTree[0]->SetBranchAddress(Form("HLT_HIPuAK4%sJet40_Eta5p1_v1_Prescl"jetType.c_str()),&jet40_p_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_HIPuAK4%sJet60_Eta5p1_v1", jetType.c_str()),&jet60_F);
    //jtTree[0]->SetBranchAddress(Form("HLT_HIPuAK4%sJet60_Eta5p1_v1_Prescl"jetType.c_str()),&jet60_p_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_HIPuAK4%sJet80_Eta5p1_v1", jetType.c_str()),&jet80_F);
    //jtTree[0]->SetBranchAddress(Form("HLT_HIPuAK4%sJet80_Eta5p1_v1_Prescl"jetType.c_str()),&jet80_p_F);
        	jtTree[0]->SetBranchAddress(Form("HLT_HIPuAK4%sJet100_Eta5p1_v1", jetType.c_str()),&jet100_F);
    //jtTree[0]->SetBranchAddress(Form("HLT_HIPuAK4%sJet100_Eta5p1_v1_Prescl"jetType.c_str()),&jet100_p_F);
        	jtTree[0]->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5_v1",&photon30_F);
        	jtTree[0]->SetBranchAddress("HLT_HIPuAK4PFBJetBSSV60_Eta2p1_v1",&ssvTrg60_F);
        	jtTree[0]->SetBranchAddress("HLT_HIPuAK4PFBJetBSSV80_Eta2p1_v1",&ssvTrg80_F);
        	jtTree[0]->SetBranchAddress("HLT_HIPuAK4PFBJetBCSV60_Eta2p1_v1",&csvTrg60_F);
        	jtTree[0]->SetBranchAddress("HLT_HIPuAK4PFBJetBCSV80_Eta2p1_v1",&csvTrg80_F);
            jtTree[0]->SetBranchAddress("HLT_HISinglePhoton30_Eta1p5_v1",&photon30_F);
        }

  // jtTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_F);
  // jtTree[0]->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_F);

    };
    void LoadEntries(int nEvt, string run){
    	jtTree[0]->GetEntry(nEvt);
    	jtTree[1]->GetEntry(nEvt);
    	jtTree[2]->GetEntry(nEvt);
    	jtTree[3]->GetEntry(nEvt);
    	jtTree[4]->GetEntry(nEvt);
    	if(run == "MC") jtTree[5]->GetEntry(nEvt);
        jtTree[6]->GetEntry(nEvt);
    };
public:
	const static int N = 7; // data does not have the runAnalyzer
	TChain * jtTree[N];
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

  // for pp prompt reco
	int jetMB_p0_F;
	int jetMB_p1_F;
	int jetMB_p2_F;
	int jetMB_p3_F;
	int jetMB_p4_F;
	int jetMB_p5_F;
	int jetMB_p6_F;
	int jetMB_p7_F;
	int jetMB_p8_F;
	int jetMB_p9_F;
	int jetMB_p10_F;
	int jetMB_p11_F;
	int jetMB_p12_F;
	int jetMB_p13_F;
	int jetMB_p14_F;
	int jetMB_p15_F;
	int jetMB_p16_F;
	int jetMB_p17_F;
	int jetMB_p18_F;
	int jetMB_p19_F;

	int L1_MB_F;

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
	ULong64_t evt_F;
	unsigned int run_F;
	unsigned int lumi_F;
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
	vector<int> *pfId_F=0;
	vector<float> *pfPt_F=0;
	vector<float> *pfVsPtInitial_F=0;
	vector<float> *pfEta_F=0;
	vector<float> *pfPhi_F=0;
    
    vector<double> *rho=0;
    vector<double> *rhoM=0;
    vector<double> *etaMax=0;

};

#endif

