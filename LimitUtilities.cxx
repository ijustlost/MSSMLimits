
//
// header file to include some limit utilities
//


#include <map>
#include <utility>
#include <iostream>
#include <cmath>

#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TCutG.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"

#include "LimitUtilities.h"
#include "LimitUtilitiesTauTau.h"

using namespace std;

void setupTLimit_Info(double mA, double tanb, int upcross, int extrap, TLimit_Info &test) {
  test.mA = mA;
  test.tanb = tanb;
  test.upcross = upcross;
  test.extrap = extrap;
}

// interpolate a TGraph to find target
// code return :   -999  failure
//                   -1  down crossing, interpolation
//                    1  up crossing, interpolation
//                  -10 (-11)  down crossing, extrapolation with target < y1 (>y2)
//                   10  (11)  up crossing, extrapolation with target < y1 (>y2)
// target is the limit cx
// g is the model cx graph
int interpolate(TGraph *g, int i, double target, double &result, 
		bool m_debug) {
  int N = g->GetN();
  if (i <0 || i+1 >= N) {
    return -999;
  }

  double x1 = -999., x2 = -999., y1 = -999., y2 = -999.;
  g->GetPoint(i,   x1, y1);
  g->GetPoint(i+1, x2, y2);

  // up or down crossing
  int sign = -1;
  if (y1 < y2) { sign = 1; }

  result =    x1 + (x2-x1)*(target - y1) / (y2 - y1);

  int extrap = 0;
  if (result < std::min(x1,x2) ) { extrap = 9;}
  else  if (result > std::max(x1,x2) ) { extrap = 10;}
  if (m_debug) {
  cout << "debug: (" << x1 << ", " << y1 << "), (" << x2 << ", " << y2 << ") : " 
       << target << "  -> " << result << "  " << (sign * (1 + extrap))
       << endl;
  }
  return (sign * (1 + extrap));
  

}

// find crossing between points i and i+1 of the TGraph with
// the line segment defined by [(x[i], target1), (x[i+1], target2)]
// code return :   -999  failure
//                   -1  down crossing, interpolation
//                    1  up crossing, interpolation
//                  -10 (-11)  down crossing, extrapolation with target < y1 (>y2)
//                   10  (11)  up crossing, extrapolation with target < y1 (>y2)
// target is the limit cx
// g is the model cx graph
int crossing(TGraph *g, int i, double target1, double target2, double &result, 
		bool m_debug) {
  int N = g->GetN();
  if (i <0 || i+1 >= N) {
    return -999;
  }

  double x1 = -999., x2 = -999., y1 = -999., y2 = -999.;
  g->GetPoint(i,   x1, y1);
  g->GetPoint(i+1, x2, y2);

  // up or down crossing
  int sign = -1;
  if (target1 > y1) { sign = 1; }

  float dx = x2-x1;
  float dy = y2-y1;
  float dt = target2 - target1;
  float xc = x1 + (y1-target1)*dx/(dt-dy);
  //float yc = y1 + (y1-target1)/(dt-dy) * dy / dx;

  //result =    x1 + (x2-x1)*(target - y1) / (y2 - y1);
  result = xc;

  int extrap = 0;
  if (result < std::min(x1,x2) ) { extrap = 9;}
  else  if (result > std::max(x1,x2) ) { extrap = 10;}
  if (m_debug) {
  cout << "debug: (" << x1 << ", " << y1 << "), (" << x2 << ", " << y2 << ") : " 
       << "(" << target1 <<", " << target2 <<") " << "  -> " << result << "  " << (sign * (1 + extrap))
       << endl;
  }
  return (sign * (1 + extrap));
  

}

//
// Scan the TGraph for possible points
bool Interpolator(TGraph *g, double target,  TString tag, 
		  std::vector< TLimit_Info > &v_liminfo, 
		  double mass, bool quiet, bool m_debug) {
  int N = g->GetN();
  TString s_result("");
  v_liminfo.clear();
  for (int i =0; i<N-1; ++i) {
    double result=0;
    int int_res = interpolate(g, i, target, result, (m_debug && !quiet));
    bool willPrint = (  (i == 0 &&abs(int_res) == 10) || 
			(i == N-2 && abs(int_res) == 11) || 
			(abs(int_res) < 5)  );
    if (willPrint) {
      int upcrossing = (int_res>0)?1:0;
      int extrap     = (abs(int_res)<5)?1:0;
      //cout << "DEBUG:  " << TString::Format("%6.3f %i%i     ", result, upcrossing,extrap) << endl;
      s_result += TString::Format("%6.3f %i%i     ", result, upcrossing,extrap);
      TLimit_Info info;
      double mA, yx;
      g->GetPoint(i, mA, yx);
      info.mA = mass;
      info.tanb = result;
      info.upcross = upcrossing;
      info.extrap = extrap;
      v_liminfo.push_back(info);
    }
  }
  
  if (s_result != "") {
    s_result = tag +" : "+ s_result;
    if (not quiet) {
      cout << s_result << endl;
    }
    return true;
  }
  if (not quiet) {
    //cout << tag << " : no result" << endl;
  }
  return false;
  
}

// another interpolator: use the mass from a different source
bool Interpolator(TGraph *g, TGraph *g_mass, TGraph *g_target,  TString tag, 
		  std::vector< TLimit_Info > &v_liminfo, 
		  double mass, bool quiet, bool m_debug, bool useM11, TGraph2D *g_out) {
  cout << "Interpolator " << tag << " quiet " << quiet << " m_debug " << m_debug << endl;
  int N = g->GetN();
  TString s_result("");
  v_liminfo.clear();
  for (int i =0; i<N-1; ++i) {
    double result=0;
    double current_mass_1, current_mass_2, current_mass, tb_1, tb_2;
    g_mass->GetPoint(i, tb_1, current_mass_1);
    g_mass->GetPoint(i+1, tb_2, current_mass_2);
    current_mass = 0.5* (current_mass_1 + current_mass_2);
    double target = g_target->Eval(current_mass);
    double target1 = g_target->Eval(current_mass_1);
    double target2 = g_target->Eval(current_mass_2);
    if (m_debug) {
      double d_mA, d_yx; g->GetPoint(i, d_mA, d_yx);
      printf("%2i %3.0f %4.2f :  %6.2f  %12.6f vs %12.6f \n", i, mass, tb_1, current_mass, d_yx, target);
    }
    if (g_out) {
      double tb = -999., cx = -999;
      g->GetPoint(i, tb, cx);
      if (tb<10) {
        int iOut = g_out->GetN();
        std::cout  << "TGraph2D mass " << mass << " tb " << tb << " " << cx - target << std::endl;
        g_out->SetPoint(iOut, mass, tb, cx - target);
      }
    }

    // result is limit (tanB). target are cx limits
    int int_res = crossing(g, i, target1, target2, result, (m_debug && !quiet));
    //int int_res = interpolate(g, i, target, result, (m_debug && !quiet));
    /*
    bool willPrint = (  (i == 0 &&abs(int_res) == 10) || 
			((i == N-2) && abs(int_res) == 11) || (useM11 && int_res==-11) ||
			(abs(int_res) < 5)  );
      */
    bool willPrint =( (abs(int_res) < 5)  );
    if (willPrint) {
      // First check we're not doing a crazy extrapolation
      double x1 = -999., x2 = -999., y1 = -999., y2 = -999.;
      g->GetPoint(i,   x1, y1);
      g->GetPoint(i+1, x2, y2);
      //cout << " result " << result << " y2 " << y2 << endl;
      if ( fabs(result-x2)<((x2-x1)*10)) {
        int upcrossing = (int_res>0)?1:0;
        int extrap     = (abs(int_res)<5)?1:0;
        //cout << "DEBUG:  " << TString::Format("%6.3f %i%i     ", result, upcrossing,extrap) << endl;
        s_result += TString::Format("%6.3f %i%i     ", result, upcrossing,extrap);
        TLimit_Info info;
        double mA, yx;
        g->GetPoint(i, mA, yx);
        info.mA = mass;
        info.tanb = result;
        info.upcross = upcrossing;
        info.extrap = extrap;
        v_liminfo.push_back(info);
      }
    }
  }
  if (s_result != "") {
    s_result = tag +" : "+ s_result;
    if (not quiet) {
      cout << s_result << endl;
    }
    return true;
  } else {
    //cout << tag << " : " << " no crossings found! " << endl;
  }
  /* if ((not quiet) && (not useM11)) { */
  /*   cout << tag << " : no result, will useM11" << endl; */
  /*   Interpolator(g, g_mass, g_target, tag, v_liminfo, mass, quiet, m_debug, true); */
  /* } */
  if (not quiet) {
    cout << tag << " : no result" << endl;
  }
  return false;
  
}


//
// TauTau interpolator:
// TGraph *g :    model cross section X BR
// TGraph *gfrac: gluon fusion fraction
// int level:  0 expected; +/- n: number of sigmas, up to 2 ; 99 observed
bool Interpolator_TauTau(TGraph *g, TGraph *gfrac, int level, TString tag, 
			 std::vector< TLimit_Info > &v_liminfo, 
			 double mass, bool quiet, bool m_debug) {
  if (m_debug) {
    cout << "DEBUG: inside " << tag << " mass: " << mass << " : " 
	 << level << endl;
  }
  int N = g->GetN();
  TString s_result("");
  v_liminfo.clear();
  for (int i =0; i<N-1; ++i) {
    double result=0;

    // calculate the target approximate from the average fraction
    // do it with Eval to make sure that the order is the same
    double x, tb1, tb2;
    g->GetPoint(i, tb1, x);
    g->GetPoint(i, tb2, x);
    double fraction = gfrac->Eval(0.5 * (tb1 + tb2));

    double target = getLimit( mass, fraction,  level);

    int int_res = interpolate(g, i, target, result, (m_debug&&quiet));
    bool willPrint = (  (i == 0 &&abs(int_res) == 10) || 
			(i == N-2 && abs(int_res) == 11) || 
			(abs(int_res) < 5)  );
    if (willPrint) {
      int upcrossing = (int_res>0)?1:0;
      int extrap     = (abs(int_res)<5)?1:0;
      //cout << "DEBUG:  " << TString::Format("%6.3f %i%i     ", result, upcrossing,extrap) << endl;
      s_result += TString::Format("%6.3f %i%i     ", result, upcrossing,extrap);
      TLimit_Info info;
      double mA, yx;
      g->GetPoint(i, mA, yx);
      info.mA = mass;
      info.tanb = result;
      info.upcross = upcrossing;
      info.extrap = extrap;
      v_liminfo.push_back(info);
    }
  }
  
  if (s_result != "") {
    s_result = tag +" : "+ s_result;
    if (not quiet) {
      cout << s_result << endl;
    }
    if (m_debug) {
      cout << "DEBUG: end " << tag << " mass: " << mass << " : " 
	   << level << " --> returning true" <<  endl;
    }

    return true;
  }
  //cout << tag << " no result" << endl;

    if (m_debug) {
      cout << "DEBUG: end " << tag << " mass: " << mass << " : " 
	   << level << " --> returning false" <<  endl;
    }

  return false;
  
}


///////
///////  Real helper functions here
///////

void add_th2f(TH2F *h1, TH2F *h2, TH2F *h_res, double c1, double c2) {
  int nx = h1->GetNbinsX();
  int ny = h1->GetNbinsY();

  for (int i=1; i<= nx; ++i) {
    for (int j=1; j<= ny; ++j) {
      double y1 = h1->GetBinContent(i,j);
      double y2 = h2->GetBinContent(i,j);

      h_res->SetBinContent(i,j,  c1*y1 + c2*y2);

    }
  }

}

void multiply_th2f(TH2F *h1, TH2F *h2, TH2F *h_res, double c1 ) {
  cout << "Inside multiply: " << (!h1) << (!h2) << (!h_res) << endl;
  if (h1 && h2 && h_res) {
    cout << "all histos checked and they are ok" << endl;
  } else {
    cout << "histos are faulty" << endl;
    return;
  }


  int nx = h1->GetNbinsX();
  int ny = h1->GetNbinsY();

  cout << "Trying to multiply histo1(" << nx << ","<<ny << ") with hist2("
       << h2->GetNbinsX() << "," << h2->GetNbinsY() << ")  "  << h2->GetName()  << endl;

  for (int i=1; i<= nx; ++i) {
    for (int j=1; j<= ny; ++j) {
      double y1 = h1->GetBinContent(i,j);
      double y2 = h2->GetBinContent(i,j);

      h_res->SetBinContent(i,j,  c1*y1*y2);

    }
  }

}



void addmult_th2f(TH2F *h1, TH2F *h2, TH2F *h3, TH2F *h_res, double c1) {
  int nx = h1->GetNbinsX();
  int ny = h1->GetNbinsY();

  for (int i=1; i<= nx; ++i) {
    for (int j=1; j<= ny; ++j) {
      double y1 = h1->GetBinContent(i,j);
      double y2 = h2->GetBinContent(i,j);
      double y3 = h3->GetBinContent(i,j);

      h_res->SetBinContent(i,j,  c1*(y1+y2)*y3);

    }
  }

}


void getfraction_th2f(TH2F *h1, TH2F *h2, TH2F *h_res) {
  int nx = h1->GetNbinsX();
  int ny = h1->GetNbinsY();

  for (int i=1; i<= nx; ++i) {
    for (int j=1; j<= ny; ++j) {
      double y1 = h1->GetBinContent(i,j);
      double y2 = h2->GetBinContent(i,j);

      h_res->SetBinContent(i,j,  y1/(y1+y2));

    }
  }

}

void getMAplot(TH2F *h_mH,  TH2F *h_res) {

  int nx = h_mH->GetNbinsX();
  int ny = h_mH->GetNbinsY();
  
  for (int i=1; i<= nx; ++i) {
    double mass = h_mH->GetXaxis()->GetBinCenter(i);
    for (int j=1; j<= ny; ++j) {
      h_res->SetBinContent(i,j, mass );
    }
  }
  
  

}


void getSantanderMatched(TH2F *h_4f, TH2F *h_5f, TH2F *h_mass, TH2F *h_res) {

  if (!h_4f) {
    cout << "ERROR!  Invalid pointer 4FS" << endl;
    return;
  }
  if (!h_5f) {
    cout << "ERROR!  Invalid pointer 4FS" << endl;
    return;
  }


  int nx = h_4f->GetNbinsX();
  int ny = h_4f->GetNbinsY();
  cout << "Inside Santander matched for " << h_4f->GetName() << endl;

  if (!h_mass) {
    cout << "Invalid pointer for h_mass: will continue with mA" << endl;
    //return;
  }

  if ((nx != h_5f->GetNbinsX()) || (ny != h_5f->GetNbinsY())) {
    cout << "ERROR: dimension mismatch in  getSantanderMatched (5FS)" << endl;
    return;
  } 
  else if (h_mass) {
    if ((nx != h_mass->GetNbinsX()) || (ny != h_mass->GetNbinsY())) {
      cout<<"ERROR: dimension mismatch in  getSantanderMatched (MASS)" << endl;
    return;    
    }
  }
  
  cout << "getSantanderMatched::INFO will work with " 
       << nx << " X " << ny << endl;
  
  //return;
  for (int i=1; i<= nx; ++i) {
    for (int j=1; j<= ny; ++j) {
      double y1 = h_4f->GetBinContent(i,j);
      double y2 = h_5f->GetBinContent(i,j);
      //double tanbeta = h_4f->GetYaxis()->GetBinLowEdge(j);
      double mass;
      /* if (j==1) */
      /* cout << "DEBUG " << i << ", " << j << " : " << y1 << " , " << y2 ; */
      if (!h_mass) {
	//cout << "inside loop for i = " << i << endl;
	//int globalbin = h_4f->GetBin(i,j);
	mass = h_4f->GetXaxis()->GetBinCenter(i);
      } else {
	mass = h_mass->GetBinContent(i,j);
	//mass = 100;
      }
      /* if (j==1) */
      /* cout <<  "  " << mass << endl; */
      if (mass < 0) {
	cout << "ERROR in getSantanderMatched! mass is " << mass << endl;
	return;
      }
      //cout << " " << mass << " , " << log(mass/4.75)-2. << endl;
      double t = log(mass/4.75)-2.;
      double s_santander= (1/(1+t)) * (y1 + t*y2);
      h_res->SetBinContent(i,j,  s_santander);
      // TEST:
      //      if (s_santander < 0 || 
      //	  (fabs(s_santander-y1)/y1 > 0.25 && fabs(s_santander-y2)/y2 > 0.25)) {
      //	cout << "getSantanderMatched::WARNING " << h_4f->GetName() << " "
      //	     << mass << "," << tanbeta << " t=" << t  << " : " << s_santander 
      //	     << " vs 4F: " << y1 << " 5F: " << y2  << endl;
      //      }

    }
  } 
  
  
}


bool getTGraph(TH2F *hin, double mA, TGraph &g) {
  // get the tanbeta list
  int N = hin->GetYaxis()->GetNbins();
  if (N <=0) return false;

  int counter=0;
  for (int i =1; i<= N; ++i) {
    double tb = hin->GetYaxis()->GetBinLowEdge(i);
    int ibin = hin->FindBin(mA, tb);
    double y = hin->GetBinContent(ibin);
    if ( y == y ) { // use only numbers that aren't NaN
      if ( y > 0 && tb >= 0.7) { // use non zero values
	g.SetPoint(counter, tb, y);
	++counter;
      }
    }
    
  }

  return true;
}

bool getTGraph(TH2F *hin_cx, TH2F* hin_br, TH2F *hmass, double mA, TGraph &g, TGraph &g_mass, double tb_low, bool use_bin_centre) {
  // get the tanbeta list
  int N = hin_br->GetYaxis()->GetNbins();
  if (N <=0) return false;

  int counter=0;
  for (int i =1; i<= N; ++i) {
    double tb;
    if (use_bin_centre) {
      tb = hin_br->GetYaxis()->GetBinCenter(i);
    } else {
      tb = hin_br->GetYaxis()->GetBinCenter(i);
    }
    //double tb = hin_br->GetYaxis()->GetBinLowEdge(i);
    int ibin = hin_br->FindBin(mA, tb);
    double cx = hin_cx->GetBinContent(ibin);
    double br = hin_br->GetBinContent(ibin);
    double mass = hmass->GetBinContent(  ibin   );
    cout <<  hin_br->GetYaxis()->GetBinLowEdge(i) << " y " << cx*br <<endl;
    if ( cx == cx and br==br and cx<1000.) { // use only numbers that aren't NaN or crazy
      if ( br>0 and cx > 0 && tb >= tb_low) { // use non zero values
        cout << " Adding point " << tb << " " << cx*br << endl;
	g.SetPoint(counter, tb, cx*br);
	g_mass.SetPoint(counter, tb, mass);
	++counter;
      }
    }
    
  }

  return true;
}



void add_lower_points( std::vector< TLimit_Info > &v_liminfo ) {
  if (v_liminfo.size()!=1) return;
  if (v_liminfo[0].upcross==0) {
      v_liminfo.push_back(v_liminfo[0]);
      v_liminfo[0].tanb = -1;
      v_liminfo[0].upcross = 1;
      v_liminfo[0].extrap = 0;
  }
}

void remove_lowest_point( std::vector< TLimit_Info > &v_liminfo ) {
  if (v_liminfo.size()<1) return;
  v_liminfo.erase( v_liminfo.begin() );
}

