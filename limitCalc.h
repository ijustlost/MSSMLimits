#include "TString.h"

class TH2F;
class TFile;
class TString;
class TGraphAsymmErrors;

void makeLimitPlot(TH2F *hin, TString tag, 
		   TGraphAsymmErrors *g_lim_obs,
		   TGraphAsymmErrors *g_lim_exp_1s,  TGraphAsymmErrors *g_lim_exp_2s,
		   TFile *foutput, TString folderName, TH2F *hfrac=0, bool useFrac=false
		   );

void makeLimitPlotMass(TH2F *hin_cx, TH2F *hin_br, TH2F *hmass, TString tag, 
		       TGraphAsymmErrors *g_lim_obs, 
		       TGraphAsymmErrors *g_lim_exp_1s,  TGraphAsymmErrors *g_lim_exp_2s,
		       TFile *foutput, TString folderName, double upper_mass=1000., 
		       TH2F *hfrac=0, bool useFrac=false, bool use_bin_centre=true
		       ) ;

void get_lims(TString tfilename="hMSSM_8TeV_1000_May18.root", TString ofilename="hmssm.root", float upper_mass=400 );
