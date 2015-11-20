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

  
  

  // plot them on top of each other. 


}
