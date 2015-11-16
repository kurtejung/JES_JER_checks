// boundaries of the pt bins, cent bins and eta bins for the runForest and plot macros.

static const double pthat[12] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540, 2000};
static const double xsecs[12] = {5.269E-01, 3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05 , 1.223E-05, 3.031E-06 , 7.746E-07, 1.410E-07 , 3.216E-08, 1.001E-08 , 0.0};

static const double weight_xsec[9] = { 7.20357e-07, 4.51655e-08, 2.6964e-09, 2.77274e-10, 3.1878e-11, 3.87126e-12, 1.62138e-12, 1.09471e-12, 4.40012e-13};
static const int nentries_file[9] = { 0, 333206, 250567, 395126, 368126, 366982, 392206, 181018, 50455};

const int ptbins[] = {50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};
const int nbins_pt = sizeof(ptbins)/sizeof(int) -1;

const double etabins[] = {-5.191, -2.650, -2.043, -1.740, -1.479, -1.131, -0.783, -0.522, 0.522, 0.783, 1.131, 1.479, 1.740, 2.043, 2.650, 5.191};
const int nbins_eta = sizeof(etabins)/sizeof(double) -1;

const int ncen=10;
const int centbins[ncen] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
const char *cdir[ncen]  = {"010","1020","2030","3040","4050","5060","6070","7080","8090","90100"};
const char *ccent[ncen] = {"0-10%","10-20%","20-30%","30-50%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};

const int knj = 1;
std::string srad[knj]={"4"};

double xmin=ptbins[0];
double xmax=ptbins[nbins_pt];
