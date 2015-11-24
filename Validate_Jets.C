// macro to read in the latest hiForests and make plots to check the Jets.
// Author: Raghav Kunnawalkam Elayavalli, Kurt Jung 
//         @ CERN for Run2
//         Nov 19th 2015

#include "boundaries.h"
#include "ForestContainer.h"

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

  //printDebug and doBjets moved to boundaries.h

  std::string bJetString = "";
  if(doBjets) bJetString = "_withBjets";

  bool skipPho30 = false;
  bool doDijetImbalance = true;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  ForestContainer jt(startfile,endfile,radius,coll,run,jetType,algo);

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

  
  //Float_t dphi = 0;
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
  //Long64_t nentries = jtTree[0]->GetEntries();
  Long64_t nentries = 1000000;
  if(printDebug) nentries = 500;
  TRandom rnd;
  

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;

    jt.LoadEntries(nEvt, run);
    
    if(run == "Data"){
      if(skipPho30 && jt.photon30_F) continue;
      if(jt.pcollisionEventSelection_F==0) continue;
      if(jt.pHBHENoiseFilter_F == 0) continue;
      if(fabs(jt.vz_F)>15) continue;
      if(jt.jet80_F) continue;
    }
    hRunN_vs_NJets->Fill(jt.run_F, jt.nref_F);
    hVz->Fill(jt.vz_F);

    // find the centrality bin.

    int centbin = -1;
    if(coll == "PP") centbin = ncen;
    if(coll == "PbPb") centbin = findBin(jt.hiBin_F);
    if(centbin == -1) continue;
    
    for(int npf = 0; npf<jt.nPFpart_F; ++npf){

      if(jetType == "Calo"){
	hPFCand_eta_vs_pT[0][centbin]->Fill(jt.pfEta_F[npf], jt.pfPt_F[npf]);
	for(int njet = 0; njet<jt.nref_F; ++njet){
	  float delR = deltaR(jt.pfEta_F[npf], jt.pfPhi_F[npf], jt.eta_F[njet], jt.phi_F[njet]);
	  if(delR < (float)radius/10)
	    hPFCand_pTscale_insideJet[0][centbin]->Fill((float)(jt.pfPt_F[npf]/jt.pt_F[njet]));	
	}
      }

      if(jetType == "PF"){
	hPFCand_eta_vs_pT[jt.pfId_F[npf]][centbin]->Fill(jt.pfEta_F[npf], jt.pfPt_F[npf]);
	for(int njet = 0; njet<jt.nref_F; ++njet){
	  float delR = deltaR(jt.pfEta_F[npf], jt.pfPhi_F[npf], jt.eta_F[njet], jt.phi_F[njet]);
	  if(delR < (float)radius/10)
	    hPFCand_pTscale_insideJet[jt.pfId_F[npf]][centbin]->Fill((float)(jt.pfPt_F[npf]/jt.pt_F[njet])); 
	}
      }
      
    }// pf cands loop
    
    // selecting on all the MB bits to get one MB bit. 
    if(jt.L1_MB_F && (jt.jetMB_p0_F || jt.jetMB_p1_F || jt.jetMB_p2_F || jt.jetMB_p3_F || jt.jetMB_p4_F || jt.jetMB_p5_F || jt.jetMB_p6_F || jt.jetMB_p7_F || jt.jetMB_p8_F || jt.jetMB_p9_F || jt.jetMB_p10_F || jt.jetMB_p11_F || jt.jetMB_p12_F || jt.jetMB_p13_F || jt.jetMB_p14_F || jt.jetMB_p15_F || jt.jetMB_p16_F || jt.jetMB_p17_F || jt.jetMB_p18_F || jt.jetMB_p19_F))
      jt.jetMB_F = 1;

    if(printDebug) cout<<"MinBias bit = "<<jt.jetMB_F<<endl;

    if(jt.jetMB_F)
      hMBSpectra[centbin]->Fill(jt.pt_F[0]);
    if(jt.jetMB_F && jt.jet40_F) hJet40andMB[centbin]->Fill(jt.pt_F[0], jt.jet40_p_F);
    if(jt.jetMB_F && jt.jet60_F) hJet60andMB[centbin]->Fill(jt.pt_F[0], jt.jet60_p_F);
    if(jt.jetMB_F && jt.jet80_F) hJet80andMB[centbin]->Fill(jt.pt_F[0], jt.jet80_p_F);
    if(jt.jetMB_F && jt.jet100_F) hJet100andMB[centbin]->Fill(jt.pt_F[0], jt.jet100_p_F);

    if(jt.jet40_F) hJet40[centbin]->Fill(jt.pt_F[0], jt.jet40_p_F);
    if(jt.jet60_F) hJet60[centbin]->Fill(jt.pt_F[0], jt.jet60_p_F);
    if(jt.jet80_F) hJet80[centbin]->Fill(jt.pt_F[0], jt.jet80_p_F);
    if(jt.jet100_F) hJet100[centbin]->Fill(jt.pt_F[0], jt.jet100_p_F);
    
    if(jt.nref_F >=3 && doDijetImbalance) {

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
      
      for (int  n = 0; n < jt.nref_F; n++) { //APPLYING CUTS AND FILLS NEW ARRAY 
	if (jt.eta_F[n] > eta_cut_min && jt.eta_F[n] < eta_cut_max) {
	  j_array_pt.push_back(jt.pt_F[n]);
	  j_array_phi.push_back(jt.phi_F[n]);
	  j_array_eta.push_back(jt.eta_F[n]);
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
	      
	      if(jt.jet80_F) hRelResponse[avgpTbin][centbin]->Fill(dijetbalanceparameter);
	      
	      //FILL OUTERETA HISTOGRAMS
	      if (probeeta < -1.3 || probeeta > 1.3){
		if(jt.jet80_F) hRelResponse_outereta[avgpTbin][centbin]->Fill(dijetbalanceparameter);
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
		if(jt.jet80_F) hRelResponse_innereta[avgpTbin][centbin]->Fill(dijetbalanceparameter);
	      }//refeta, probeeta if statement
	      
	    }// avg pT bin

	  }// delta phi selection

	}// refpt!=-, probpt!=0

      }// alpha selection, suppression of third jet
      j_array_pt.clear();
      j_array_phi.clear();
      j_array_eta.clear();
      
    }
    
    for(int ijet=0; ijet<jt.nref_F; ijet++){

      if(jt.jet40_F) hJet40All[centbin]->Fill(jt.pt_F[ijet], jt.jet40_p_F);
      if(jt.jet60_F) hJet60All[centbin]->Fill(jt.pt_F[ijet], jt.jet60_p_F);
      if(jt.jet80_F) hJet80All[centbin]->Fill(jt.pt_F[ijet], jt.jet80_p_F);
      if(jt.jet100_F) hJet100All[centbin]->Fill(jt.pt_F[ijet], jt.jet100_p_F);

      if(doBjets){
	discrSSVHE[centbin]->Fill(jt.discr_ssvHighEff_F[ijet]);
	discrSSVHP[centbin]->Fill(jt.discr_ssvHighPur_F[ijet]);
	discrCSV[centbin]->Fill(jt.discr_csv_F[ijet]);

	if(jt.csvTrg60_F){
		csvTrg60[centbin]->Fill(jt.pt_F[ijet]);
		if(jt.discr_csv_F[ijet]>0.9) csvTrg60withCSV[centbin]->Fill(jt.pt_F[ijet]);
	}
	if(jt.csvTrg80_F){
		csvTrg80[centbin]->Fill(jt.pt_F[ijet]);
		if(jt.discr_csv_F[ijet]>0.9) csvTrg80withCSV[centbin]->Fill(jt.pt_F[ijet]);
	}
	if(jt.ssvTrg60_F){
		ssvTrg60[centbin]->Fill(jt.pt_F[ijet]);
		if(jt.discr_ssvHighPur_F[ijet]>1.2) ssvTrg60withSSVHP[centbin]->Fill(jt.pt_F[ijet]);
	}
	if(jt.ssvTrg80_F){
		ssvTrg80[centbin]->Fill(jt.pt_F[ijet]);
		if(jt.discr_ssvHighPur_F[ijet]>1.2) ssvTrg80withSSVHP[centbin]->Fill(jt.pt_F[ijet]);
	}

	if(jt.discr_csv_F[ijet]>0.9) csvDistr[centbin]->Fill(jt.pt_F[ijet]);
	if(jt.discr_ssvHighPur_F[ijet]>1.2) ssvDistr[centbin]->Fill(jt.pt_F[ijet]);

      }

      if(jt.nref_F >= 2) hDeltaPhi[centbin]->Fill(deltaphi(jt.phi_F[0], jt.phi_F[1]));
      
      if(jt.jet80_F && jt.pt_F[0] >LjCut && jt.pt_F[1] >SbjCut) { 
	hJetEta[centbin]->Fill(jt.eta_F[ijet]);
	hJetPhi[centbin]->Fill(jt.phi_F[ijet]);
	hDeltaPhi[centbin]->Fill(deltaphi(jt.phi_F[0], jt.phi_F[1]));
      }

      hJetpT[centbin]->Fill(jt.pt_F[ijet]);

      hJEC[centbin]->Fill(jt.rawpt_F[ijet], jt.eta_F[ijet], (float)(jt.pt_F[ijet]/jt.rawpt_F[ijet]));
      
    }// njets

    if(run == "MC"){
      if(jt.pt_F[0]>LjCut && jt.pt_F[1] >SbjCut){
	float Aj = (float)(jt.pt_F[0]-jt.pt_F[1])/(jt.pt_F[0]+jt.pt_F[1]);
	pt2overpt1[centbin]->Fill((float)jt.pt_F[1]/jt.pt_F[0]);
	hAj[centbin]->Fill(Aj);
      }
    }//mc

    if(run == "Data"){
      if(jt.jet80_F && jt.pt_F[0]>LjCut && jt.pt_F[1] >SbjCut){
	float Aj = (float)(jt.pt_F[0]-jt.pt_F[1])/(jt.pt_F[0]+jt.pt_F[1]);
	pt2overpt1[centbin]->Fill((float)jt.pt_F[1]/jt.pt_F[0]);
	hAj[centbin]->Fill(Aj);
      }

    } //data

  }// nevents

  TH1F * hJet40Turnon[ncen+1];
  TH1F * hJet60Turnon[ncen+1];
  TH1F * hJet80Turnon[ncen+1];
  TH1F * hJet100Turnon[ncen+1];
  
  if(run == "Data"){

    if(coll == "PbPb"){
      for(int icen = 0; icen<ncen; ++icen){
	hJet40Turnon[icen] = (TH1F*)hJet40andMB[icen]->Clone(Form("hJet40Turnon_%s",cdir[icen]));
	hJet60Turnon[icen] = (TH1F*)hJet60andMB[icen]->Clone(Form("hJet60Turnon_%s",cdir[icen]));
	hJet80Turnon[icen] = (TH1F*)hJet80andMB[icen]->Clone(Form("hJet80Turnon_%s",cdir[icen]));
	hJet100Turnon[icen] = (TH1F*)hJet100andMB[icen]->Clone(Form("hJet100Turnon_%s",cdir[icen])); 
      hJet40Turnon[icen]->Divide(hMBSpectra[icen]);
      hJet60Turnon[icen]->Divide(hMBSpectra[icen]);
      hJet80Turnon[icen]->Divide(hMBSpectra[icen]);
      hJet100Turnon[icen]->Divide(hMBSpectra[icen]);
      }
    }

    if(coll == "PP"){
      hJet40Turnon[ncen] = (TH1F*)hJet40andMB[ncen]->Clone("hJet40Turnon_PP");
      hJet60Turnon[ncen] = (TH1F*)hJet60andMB[ncen]->Clone("hJet60Turnon_PP");
      hJet80Turnon[ncen] = (TH1F*)hJet80andMB[ncen]->Clone("hJet80Turnon_PP");
      hJet100Turnon[ncen] = (TH1F*)hJet100andMB[ncen]->Clone("hJet100Turnon_PP"); 
      hJet40Turnon[ncen]->Divide(hMBSpectra[ncen]);
      hJet60Turnon[ncen]->Divide(hMBSpectra[ncen]);
      hJet80Turnon[ncen]->Divide(hMBSpectra[ncen]);
      hJet100Turnon[ncen]->Divide(hMBSpectra[ncen]);
    }
  }
    
  
  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;  
  
}
