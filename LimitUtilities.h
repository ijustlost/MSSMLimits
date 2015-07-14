#ifndef LimitUtilities_h
#define LimitUtilities_h

#include<vector>
#include "TString.h"

class TGraph;
class TGraph2D;
class TH2F;

struct TLimit_Info {
  double mA;
  double tanb;
  int upcross;
  int extrap;
};

void setupTLimit_Info(double mA, double tanb, int upcross, int extrap, TLimit_Info &test);

int interpolate(TGraph *g, int i, double target, double &result, 
		bool m_debug=false);

int crossing(TGraph *g, int i, double target1, double target2, double &result, 
		bool m_debug=false);

bool Interpolator(TGraph *g, double target,  TString tag, 
		  std::vector< TLimit_Info > &v_liminfo, 
		  double mass, bool quiet=false, bool m_debug=false);

bool Interpolator(TGraph *g, TGraph *g_mass, TGraph *g_target,  TString tag, 
		  std::vector< TLimit_Info > &v_liminfo, 
		  double mass, bool quiet=false, bool m_debug=false, bool useM11=false, TGraph2D *g_out=0);

bool Interpolator_TauTau(TGraph *g, TGraph *gfrac, int level, TString tag, 
			 std::vector< TLimit_Info > &v_liminfo, 
			 double mass, bool quiet=false, bool m_debug=false);

void add_th2f(TH2F *h1, TH2F *h2, TH2F *h_res, double c1 = 1., double c2=1.);

void multiply_th2f(TH2F *h1, TH2F *h2, TH2F *h_res, double c1 = 1.);

void addmult_th2f(TH2F *h1, TH2F *h2, TH2F *h3, TH2F *h_res, double c1 = 1.);

void getfraction_th2f(TH2F *h1, TH2F *h2, TH2F *h_res);

void getMAplot(TH2F *h_mH,  TH2F *h_res);
void getSantanderMatched(TH2F *h_4f, TH2F *h_5f, TH2F *h_mass, TH2F *h_res);
bool getTGraph(TH2F *hin, double mA, TGraph &g);
bool getTGraph(TH2F *hin_cx, TH2F* hin_br, TH2F *hmass, double mA, TGraph &g, TGraph &g_mass, double tb_low=0.7, bool use_bin_centre=true);
void add_lower_points( std::vector< TLimit_Info > &v_liminfo );
void remove_lowest_point( std::vector< TLimit_Info > &v_liminfo );
#endif
