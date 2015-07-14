// root macro for simple limit calculations
//
//
#include <map>
#include <utility>
#include <iostream>

#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TCutG.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "BlobUtils.h"

#include "LimitUtilities.h"
#include "BlobUtils.h"
#include "LimitCalc.h"

// this is to define the cross section
// limits TGraphAsymmErrors if you don't have them in
// a single file
//#include "LimitUtilitiesSetupLimits.h"

using namespace std;

//const bool quiet=true;
//const bool dodebug=false;

/*
void makeLimitPlot(TH2F *hin, TString tag, 
		   TGraphAsymmErrors *g_lim_obs,
		   TGraphAsymmErrors *g_lim_exp_1s,  TGraphAsymmErrors *g_lim_exp_2s,
		   TFile *foutput, TString folderName, TH2F *hfrac, bool useFrac
		   );
       */
void makeLimitPlotMass(TH2F *hin, TH2F *hmass, TString tag, 
		       TGraphAsymmErrors *g_lim_obs, 
		       TGraphAsymmErrors *g_lim_exp_1s,  TGraphAsymmErrors *g_lim_exp_2s,
		       TFile *foutput, TString folderName, double upper_mass, 
		       TH2F *hfrac, bool useFrac
		       );
void save( std::vector<Blob> v_blobs, TString name);

void save( std::vector<Blob> v_blobs, TString name) {
  std::cout << "Found "<< v_blobs.size() << " blobs for " << name  << std::endl;
  int iBlob=0;
  for ( std::vector<Blob>::iterator it = v_blobs.begin() ;
      it != v_blobs.end(); it++ ) {
    (*it).Print();
    std::vector<TCutG*> tcuts = (*it).TCutGs(Form("%s_%i",name.Data(), iBlob), 0.49);
    for (auto tc : tcuts) {
      tc->Write();
    }
    iBlob++;
  }
}


// TH2F which is mA,tanbeta versus xsec X BR
//
// Idea of the function: give in a TH2F the sigma X BR and compare at each point 
// with the limit, which is given in TGraphAsymmErrors format. Write output in a
// given root file.
/*
void makeLimitPlot(TH2F *hin, TString tag, 
		   TGraphAsymmErrors *g_lim_obs,
		   TGraphAsymmErrors *g_lim_exp_1s,  TGraphAsymmErrors *g_lim_exp_2s,
		   TFile *foutput, TString folderName, TH2F *hfrac, bool useFrac
		   ) {
  // output that you want: 
  TGraphAsymmErrors *g_high   = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_high_0 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_high_1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_high_2 = new TGraphAsymmErrors();
  TCutG *gc_high   ;//= new TCutG();
  TCutG *gc_high_0 ;//= new TCutG();
  g_high->SetName("Observed");
  g_high_0->SetName("Expected");
  g_high_1->SetName("Expected1sigma");
  g_high_2->SetName("Expected2sigma");

  TGraphAsymmErrors *g_low   = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_low_0 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_low_1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_low_2 = new TGraphAsymmErrors();
  TCutG *gc_low   = new TCutG();
  TCutG *gc_low_0 = new TCutG();
  g_low->SetName("Observed");
  g_low_0->SetName("Expected");
  g_low_1->SetName("Expected1sigma");
  g_low_2->SetName("Expected2sigma");

  int counter_high_obs = 0;
  int counter_high_exp = 0;
  int counter_low_obs = 0;
  int counter_low_exp = 0;

  // for each point in the TGraphAsymmErrors do the interpolatio
  int N = g_lim_obs->GetN();
  //cout << "Will do for " << N << endl;
  for (int i=0; i < N; ++i) {
    // cout << "Looping: " << i  << endl;

    std::vector< TLimit_Info > v_info_obs;
    std::vector< TLimit_Info > v_info_exp;
    std::vector< TLimit_Info > v_info_p2s;
    std::vector< TLimit_Info > v_info_p1s;
    std::vector< TLimit_Info > v_info_m1s;
    std::vector< TLimit_Info > v_info_m2s;


    double mass, lim_obs, lim_exp, lim_p1s, lim_m1s, lim_p2s, lim_m2s;
    g_lim_obs->GetPoint(i, mass, lim_obs);

    TGraph g;
    getTGraph(hin, mass, g);
    TString name = tag + " " + TString::Format("%4.0f", mass);
    // do not use the actual input limits, but just do the tautau lims
    if (useFrac) {
      //cout << "Will run interpolator for " << name << endl;
      TGraph gfrac, g1;
      getTGraph(hin, hfrac,  mass, g1, gfrac);
      //cout << "ready to run" << endl;
      Interpolator_TauTau(&g1, &gfrac, 99, name + " obs", v_info_obs, mass, false,dodebug);
      Interpolator_TauTau(&g1, &gfrac,  0, name + " exp", v_info_exp, mass, quiet,dodebug);
      Interpolator_TauTau(&g1, &gfrac,  2, name + " +2s", v_info_p2s, mass, quiet,dodebug);
      Interpolator_TauTau(&g1, &gfrac,  1, name + " +1s", v_info_p1s, mass, quiet,dodebug);
      Interpolator_TauTau(&g1, &gfrac, -1, name + " -1s", v_info_m1s, mass, quiet,dodebug);
      Interpolator_TauTau(&g1, &gfrac, -2, name + " -2s", v_info_m2s, mass, quiet,dodebug);
      //cout << "Interpolator run successfully" << endl;
    }
    else { // default method: use the actual input limits
      g_lim_exp_1s->GetPoint(i, mass, lim_exp);
      lim_p1s = lim_exp + g_lim_exp_1s->GetErrorYhigh(i);
      lim_m1s = lim_exp + g_lim_exp_1s->GetErrorYlow(i);
      lim_p2s = lim_exp + g_lim_exp_2s->GetErrorYhigh(i);
      lim_m2s = lim_exp + g_lim_exp_2s->GetErrorYlow(i);
      //dodebug = true;
      Interpolator(&g, lim_obs, name + " obs", v_info_obs, mass, false,dodebug);
      Interpolator(&g, lim_exp, name + " exp", v_info_exp, mass, quiet,dodebug);
      Interpolator(&g, lim_p2s, name + " +2s", v_info_p2s, mass, quiet,dodebug);
      Interpolator(&g, lim_p1s, name + " +1s", v_info_p1s, mass, quiet,dodebug);
      Interpolator(&g, lim_m1s, name + " -1s", v_info_m1s, mass, quiet,dodebug);
      Interpolator(&g, lim_m2s, name + " -2s", v_info_m2s, mass, quiet,dodebug);
    }


    int ncur=0;
    // fill in the TGraphAsymmErrors
    //cout << "Fill in the observed" << endl;
    /////////////////////////////////////////////////////////////
    ncur = v_info_obs.size();
    double ma_prev = 0; // this is in order to keep the previous mass value; if you see a mass twice, keep the second value
    for (int i=0; i < ncur; ++i) {
      TLimit_Info info = v_info_obs[i];
      //cout << "DEBUG: " << i << ", " << info.upcross << info.extrap << " mA=" << info.mA << " tb=" << info.tanb << endl;
      // see the last that is 11 or 10 and fill the high histo
      if ((i == ncur - 1)  && (info.upcross == 1) && (info.extrap == 1 || info.extrap == 0)) {
        g_high->SetPoint(counter_high_obs, info.mA, info.tanb);
        ++counter_high_obs;
      }
      if (( (info.upcross == 0) && (info.extrap == 1) ) ||  (  (info.upcross == 0) && (info.extrap == 0) ) ) {
        if (ma_prev == info.mA) {
          g_low->SetPoint(counter_low_obs-1, info.mA, info.tanb);
        } else {
          g_low->SetPoint(counter_low_obs, info.mA, info.tanb);
          ++counter_low_obs;	
        }
        ma_prev = info.mA;
      }
    }
    /////////////////////////////////////////////////////////////
    ncur = v_info_exp.size(); 

    // here you assume that all the bands exist!
    bool found_discrepancy = false;
    if ( ncur != int(v_info_p2s.size())  || ncur != int(v_info_p1s.size()) 
        || ncur != int(v_info_m1s.size()) || ncur != int(v_info_m2s.size()) ) {
      // seems to be happening in the +2sigma at mH=230??
      cout << "WARNING!!! expected has: " << ncur << " entries and the bands: " << endl;
      cout << "+2 sigma:  " << v_info_p2s.size();
      cout << "  +1 sigma:  " << v_info_p1s.size();
      cout << "  -1 sigma:  " << v_info_m1s.size();
      cout << "  -2 sigma:  " << v_info_m2s.size() << endl;
      found_discrepancy = true;
    }
    ma_prev = 0;
    for (int i=0; i < ncur; ++i) {
      TLimit_Info info = v_info_exp[i];
      // if there is a expected limit, then there is also the band
      TLimit_Info info_p2s; //= v_info_p2s[i];
      TLimit_Info info_p1s; //= v_info_p1s[i];
      TLimit_Info info_m1s; //= v_info_m1s[i];
      TLimit_Info info_m2s; //= v_info_m2s[i];
      if (!found_discrepancy) {
        info_p2s = v_info_p2s[i];
        info_p1s = v_info_p1s[i];
        info_m1s = v_info_m1s[i];
        info_m2s = v_info_m2s[i];
      }
      //cout << "debug: " << i << " : " << info.upcross << " " << info.extrap << " :  " << counter_high_exp << endl;
      // see the last that is 11 or 10 and fill the high histo

      if ((i == ncur - 1)  && (info.upcross == 1) && (info.extrap == 1 || info.extrap == 0)) {
        g_high_0->SetPoint(counter_high_exp, info.mA, info.tanb);
        g_high_1->SetPoint(counter_high_exp, info.mA, info.tanb);
        g_high_2->SetPoint(counter_high_exp, info.mA, info.tanb);
        if (!found_discrepancy) {
          g_high_1->SetPointError(counter_high_exp, 0.,0., 
              info.tanb - info_m1s.tanb,info_p1s.tanb -info.tanb);
          g_high_2->SetPointError(counter_high_exp, 0.,0., 
              info.tanb - info_m2s.tanb,info_p2s.tanb -info.tanb);
        }
        ++counter_high_exp;
      }
      if ( ( (info.upcross == 0) && (info.extrap == 1) )  || ( (info.upcross == 0) && (info.extrap == 0) )  ) {
        if (ma_prev == info.mA) {
          g_low_0->SetPoint(counter_low_exp-1, info.mA, info.tanb);
          g_low_1->SetPoint(counter_low_exp-1, info.mA, info.tanb);
          g_low_2->SetPoint(counter_low_exp-1, info.mA, info.tanb);
          if (!found_discrepancy) {
            g_low_1->SetPointError(counter_low_exp-1, 0.,0., 
                info.tanb - info_m1s.tanb,info_p1s.tanb -info.tanb);
            g_low_2->SetPointError(counter_low_exp-1, 0.,0., 
                info.tanb - info_m2s.tanb,info_p2s.tanb -info.tanb);
          }

        } else {
          g_low_0->SetPoint(counter_low_exp, info.mA, info.tanb);
          g_low_1->SetPoint(counter_low_exp, info.mA, info.tanb);
          g_low_2->SetPoint(counter_low_exp, info.mA, info.tanb);
          if (!found_discrepancy) {
            g_low_1->SetPointError(counter_low_exp, 0.,0., 
                info.tanb - info_m1s.tanb,info_p1s.tanb -info.tanb);
            g_low_2->SetPointError(counter_low_exp, 0.,0., 
                info.tanb - info_m2s.tanb,info_p2s.tanb -info.tanb);
          }
          ++counter_low_exp;	
        }
        ma_prev = info.mA;
      }
    }
    //cout << "finished with filling in limits" << endl;
  }
  //cout << "DEBUG: obs points " << counter_high_obs << ", " << counter_low_obs << endl;
  //cout << "DEBUG: exp points " << counter_high_exp << ", " << counter_low_exp << endl;

  // now make the TGCut plots
  gc_high  = new TCutG("gc_high", counter_high_obs+3);
  gc_high_0  = new TCutG("gc_high_0", counter_high_exp+3);

  gc_low  = new TCutG("gc_low", counter_low_obs+3);
  gc_low_0  = new TCutG("gc_low_0", counter_low_exp+3);


  gc_high->SetName("Observed_CutG");
  gc_high_0->SetName("Expected_CutG");
  gc_low->SetName("Observed_CutG");
  gc_low_0->SetName("Expected_CutG");

  double m_init, dump, m_init_0;
  //////////////////////////////////////////////////////
  g_high->GetPoint(0, m_init, dump);
  gc_high->SetPoint(0, m_init, 100.);  
  m_init_0 = m_init;
  for (int i=0; i < g_high->GetN(); ++i) {
    g_high->GetPoint(i, m_init, dump);
    gc_high->SetPoint(i+1, m_init, dump);
  }
  gc_high->SetPoint(counter_high_obs+1, m_init, 100.);
  gc_high->SetPoint(counter_high_obs+2, m_init_0, 100.);
  //  --> this is not correct in the edges, need to do one more interpolation
  g_low->GetPoint(0, m_init, dump);
  gc_high->SetPoint(0, m_init, 0.);  
  m_init_0 = m_init;
  for (int i=0; i < g_low->GetN(); ++i) {
    g_low->GetPoint(i, m_init, dump);
    gc_low->SetPoint(i+1, m_init, dump);
  }
  gc_low->SetPoint(counter_low_obs+1, m_init, 0.);
  gc_low->SetPoint(counter_low_obs+2, m_init_0, 0.);
  //////////////////////////////////////////////////////
  g_high_0->GetPoint(0, m_init, dump);
  gc_high_0->SetPoint(0, m_init, 100.); 
  m_init_0 = m_init;
  for (int i=0; i < g_high_0->GetN(); ++i) {
    g_high_0->GetPoint(i, m_init, dump);
    gc_high_0->SetPoint(i+1, m_init, dump);
  }
  gc_high_0->SetPoint(counter_high_exp+1, m_init, 100.);
  gc_high_0->SetPoint(counter_high_exp+2, m_init_0, 100.);
  // ................................................
  g_low_0->GetPoint(0, m_init, dump);
  gc_high_0->SetPoint(0, m_init, 0.);  
  m_init_0 = m_init;
  for (int i=0; i < g_low_0->GetN(); ++i) {
    g_low_0->GetPoint(i, m_init, dump);
    gc_low_0->SetPoint(i+1, m_init, dump);
  }
  gc_low_0->SetPoint(counter_low_obs+1, m_init, 0.);
  gc_low_0->SetPoint(counter_low_obs+2, m_init_0, 0.);
  //////////////////////////////////////////////////////

  gc_low_0->SetFillColor(1);
  gc_low_0->SetFillStyle(3004);

  gc_high_0->SetFillColor(1);
  gc_high_0->SetFillStyle(3004);



  // second part: write things in the output file
  foutput->mkdir(folderName+"_high");
  foutput->mkdir(folderName+"_low");
  foutput->cd(folderName+"_high");
  g_high->Write("Observed");
  g_high_0->Write("Expected");
  g_high_1->Write("Expected1sigma"); 
  g_high_2->Write("Expected2sigma"); 
  gc_high->Write("Observed_CutG"); 
  gc_high_0->Write("Expected_CutG");

  foutput->cd("../");
  foutput->cd(folderName+"_low");
  g_low->Write("Observed");
  g_low_0->Write("Expected");
  g_low_1->Write("Expected1sigma"); 
  g_low_2->Write("Expected2sigma"); 
  gc_low->Write("Observed_CutG"); 
  gc_low_0->Write("Expected_CutG");


  foutput->cd("../");


} // end function: makeLimitPlot( ... )
*/


//
// this is to address the fact that you need an mA-tanbeta limit, but you
// have to constrain a different mass, e.g. mH;  
// The method also scans all the available mass points in the hMSSM ntuple grid
// and compares with the linearly interpolated limit from the TGraph.
// This may avoid the artificial behaviour of the limit around thresholds
// hin   : the sigmaxBR for mA--tanbeta
// hmass : the mass for mA--tanbeta 
void makeLimitPlotMass(TH2F *hin_cx, TH2F *hin_br, TH2F *hmass, TString tag, 
		       TGraphAsymmErrors *g_lim_obs, 
		       TGraphAsymmErrors *g_lim_exp_1s,  TGraphAsymmErrors *g_lim_exp_2s,
		       TFile *foutput, TString folderName, double upper_mass, 
		       TH2F *hfrac, bool useFrac, bool use_bin_centre
		       ) {
  // some preparation
  TGraphAsymmErrors *g_lim_exp = new TGraphAsymmErrors(*g_lim_exp_1s);

  TGraphAsymmErrors *g_lim_exp_p1s = new TGraphAsymmErrors(*g_lim_exp_1s);
  TGraphAsymmErrors *g_lim_exp_m1s = new TGraphAsymmErrors(*g_lim_exp_1s);
  for (int i=0; i<g_lim_exp_1s->GetN(); ++i) {
    double x, y, ymin, ymax;
    g_lim_exp_1s->GetPoint(i, x, y);
    ymin = g_lim_exp_1s->GetErrorYlow(i);
    ymax = g_lim_exp_1s->GetErrorYhigh(i);
    g_lim_exp_p1s->SetPoint(i, x, y+ymax);
    g_lim_exp_m1s->SetPoint(i, x, y-ymin);
  }

  TGraphAsymmErrors *g_lim_exp_p2s = new TGraphAsymmErrors(*g_lim_exp_2s);
  TGraphAsymmErrors *g_lim_exp_m2s = new TGraphAsymmErrors(*g_lim_exp_2s);
  for (int i=0; i<g_lim_exp_2s->GetN(); ++i) {
    double x, y, ymin, ymax;
    g_lim_exp_2s->GetPoint(i, x, y);
    ymin = g_lim_exp_2s->GetErrorYlow(i);
    ymax = g_lim_exp_2s->GetErrorYhigh(i);
    g_lim_exp_p2s->SetPoint(i, x, y+ymax);
    g_lim_exp_m2s->SetPoint(i, x, y-ymin);
  }

  // output that you want: 
  TGraphAsymmErrors *g_high   = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_high_0 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_high_1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_high_2 = new TGraphAsymmErrors();
  //TCutG *gc_high   ;//= new TCutG();
  //TCutG *gc_high_0 ;//= new TCutG();
  g_high->SetName("Observed");
  g_high_0->SetName("Expected");
  g_high_1->SetName("Expected1sigma");
  g_high_2->SetName("Expected2sigma");

  //TGraphAsymmErrors *g_low   = new TGraphAsymmErrors();
  //TGraphAsymmErrors *g_low_0 = new TGraphAsymmErrors();
  //TGraphAsymmErrors *g_low_1 = new TGraphAsymmErrors();
  //TGraphAsymmErrors *g_low_2 = new TGraphAsymmErrors();
  //TCutG *gc_low   = new TCutG();
  //TCutG *gc_low_0 = new TCutG();
  //TCutG *gc_low_m1s = new TCutG();
  //TCutG *gc_low_m2s = new TCutG();
  //TCutG *gc_low_p1s = new TCutG();
  //TCutG *gc_low_p2s = new TCutG();
  //g_low->SetName("Observed");
  //g_low_0->SetName("Expected");
  //g_low_1->SetName("Expected1sigma");
  //g_low_2->SetName("Expected2sigma");

  //int counter_high_obs = 0;
  //int counter_high_exp = 0;
  //int counter_low_obs = 0;
  //int counter_low_exp = 0;

  // for each point calculate the limit
  int N = hin_cx->GetXaxis()->GetNbins();

  int n_limit = g_lim_obs->GetN();
  double mass_low, mass_high, lim; 
  g_lim_obs->GetPoint(0, mass_low, lim);
  g_lim_obs->GetPoint(n_limit-1, mass_high, lim);

  double tb_low = 0.99;
  // Assume equal spacing in mH
  float deltaM = hin_cx->GetXaxis()->GetBinWidth(1);

  std::vector< Blob > v_blobs_obs;
  std::vector< Blob > v_blobs_exp;
  std::vector< Blob > v_blobs_p2s;
  std::vector< Blob > v_blobs_p1s;
  std::vector< Blob > v_blobs_m1s;
  std::vector< Blob > v_blobs_m2s;
  
  TGraph2D g_out_p2s;
  g_out_p2s.SetName("g2d_p2s");
  TGraph2D g_out_m2s;
  g_out_m2s.SetName("g2d_m2s");
  TGraph2D g_out_p1s;
  g_out_p1s.SetName("g2d_p1s");
  TGraph2D g_out_m1s;
  g_out_m1s.SetName("g2d_m1s");
  TGraph2D g_out_exp;
  g_out_exp.SetName("g2d_exp");
  TGraph2D g_out_obs;
  g_out_obs.SetName("g2d_obs");

  for (int i=1; i<=N; ++i) {
    double mass = hin_cx->GetXaxis()->GetBinCenter(i);
    if (mass < 139 || mass > upper_mass) continue;
    //if (tag.Contains("A->Zh") && mass < 220.) continue;
    //if (mass < mass_low || mass > mass_high) continue;

    if (tag.Contains("H->WW") && mass < 270.) continue;

    TGraph g_XsecBR;
    TGraph g_mass;

    //if (tag.Contains("H->ZZ") && mass==240) tb_low = 1.0;
    getTGraph(hin_cx, hin_br, hmass,  mass, g_XsecBR, g_mass, tb_low, use_bin_centre);
    g_XsecBR.SetName( Form("cxbr_%.0f", mass));
    g_mass.SetName( Form("mass_%.0f", mass));
    foutput->cd();
    g_XsecBR.Write();
    g_mass.Write();
    g_lim_obs->Write();
    
    std::vector< TLimit_Info > v_info_obs;
    std::vector< TLimit_Info > v_info_exp;
    std::vector< TLimit_Info > v_info_p2s;
    std::vector< TLimit_Info > v_info_p1s;
    std::vector< TLimit_Info > v_info_m1s;
    std::vector< TLimit_Info > v_info_m2s;

    //double lim_mass, lim_obs, lim_exp, lim_p1s, lim_m1s, lim_p2s, lim_m2s;

    TString name = tag + " " + TString::Format("%4.0f", mass);
    bool quiet=false;
    bool dodebug = true;
    if (useFrac) {
      // here you will assume mH = mA
      TGraph gfrac, g_XsecBR1;
      getTGraph(hin_cx, hin_br, hfrac,  mass, g_XsecBR1, gfrac, use_bin_centre);
      Interpolator_TauTau(&g_XsecBR1, &gfrac, 99, name + " obs", v_info_obs, mass, false,dodebug);
      Interpolator_TauTau(&g_XsecBR1, &gfrac,  0, name + " exp", v_info_exp, mass, quiet,dodebug);
      Interpolator_TauTau(&g_XsecBR1, &gfrac,  2, name + " +2s", v_info_p2s, mass, quiet,dodebug);
      Interpolator_TauTau(&g_XsecBR1, &gfrac,  1, name + " +1s", v_info_p1s, mass, quiet,dodebug);
      Interpolator_TauTau(&g_XsecBR1, &gfrac, -1, name + " -1s", v_info_m1s, mass, quiet,dodebug);
      Interpolator_TauTau(&g_XsecBR1, &gfrac, -2, name + " -2s", v_info_m2s, mass, quiet,dodebug);
    } else {
      // HZZ interpolator
      bool forceM11 = false;
      // if (tag.Contains("H->ZZ") && (mass == 350)) {forceM11 = true;}
      // if (tag.Contains("H->WW") && (mass == 280)) {forceM11 = true;}
      Interpolator(&g_XsecBR, &g_mass, g_lim_obs,name + " obs", v_info_obs, mass, false,dodebug, forceM11, &g_out_obs);

      // forceM11 = false;
      // if (tag.Contains("H->ZZ") && (mass == 240)) {forceM11 = true;}
      // if (tag.Contains("H->WW") && (mass == 280)) {forceM11 = true;}
      Interpolator(&g_XsecBR, &g_mass, g_lim_exp,name + " exp", v_info_exp, mass, false,dodebug, forceM11, &g_out_exp);
      // dodebug = false;
      // forceM11 = false;
      Interpolator(&g_XsecBR, &g_mass, g_lim_exp_p2s, name + " +2s", v_info_p2s, mass, quiet,dodebug, forceM11, &g_out_p2s);
      Interpolator(&g_XsecBR, &g_mass, g_lim_exp_p1s, name + " +1s", v_info_p1s, mass, quiet,dodebug, forceM11, &g_out_p1s);
      Interpolator(&g_XsecBR, &g_mass, g_lim_exp_m1s, name + " -1s", v_info_m1s, mass, quiet,dodebug, forceM11, &g_out_m1s);
      Interpolator(&g_XsecBR, &g_mass, g_lim_exp_m2s, name + " -2s", v_info_m2s, mass, quiet,dodebug, forceM11, &g_out_m2s);


      if (tag.Contains("A->Zh") && mass < 220. && v_info_obs.size() > 0) {      
        TLimit_Info info_obs;
        setupTLimit_Info(220., v_info_obs[0].tanb, v_info_obs[0].upcross, v_info_obs[0].extrap, info_obs);
        v_info_obs.clear();
        v_info_obs.push_back(info_obs);
      }
      if (tag.Contains("A->Zh") && mass < 220. && v_info_exp.size() > 0) {      
        TLimit_Info info_exp;
        setupTLimit_Info(220., v_info_exp[0].tanb, v_info_exp[0].upcross, v_info_exp[0].extrap, info_exp);
        v_info_exp.clear();
        v_info_exp.push_back(info_exp);
      }
       

    }

    if (dodebug and not quiet) cout << "Done interp" << endl;

    // Form segments
    update_blobs( v_info_obs, v_blobs_obs, deltaM, tb_low );
    update_blobs( v_info_exp, v_blobs_exp, deltaM, tb_low );
    update_blobs( v_info_p2s, v_blobs_p2s, deltaM, tb_low );
    update_blobs( v_info_m2s, v_blobs_m2s, deltaM, tb_low );
    update_blobs( v_info_p1s, v_blobs_p1s, deltaM, tb_low );
    update_blobs( v_info_m1s, v_blobs_m1s, deltaM, tb_low );

    if (dodebug and not quiet) cout << "Updated blobs" << endl;

  } // End loop over masses

  g_out_p2s.Write();
  g_out_p1s.Write();
  g_out_m2s.Write();
  g_out_m1s.Write();
  g_out_obs.Write();
  g_out_exp.Write();

  // Save blobs
  save(v_blobs_obs, "Obs");
  save(v_blobs_p2s, "Plus2sigma");
  save(v_blobs_p1s, "Plus1sigma");
  save(v_blobs_exp, "Exp");
  save(v_blobs_m1s, "Minus1sigma");
  save(v_blobs_m2s, "Minus2sigma");


} // end function: makeLimitPlotMass( ... )



void get_lims(TString tfilename, TString ofilename, float upper_mass ) {

  //TString tfilename = "hMSSM_8TeV_1000_May18.root";

  TFile *f_xsecs = new TFile(tfilename,"read");

  TH2F *h_mh = (TH2F*) f_xsecs->Get("h_mh");
  TH2F *h_mH = (TH2F*) f_xsecs->Get("h_mH");
  TH2F *h_mA = (TH2F*) f_xsecs->Get("h_mh");
  getMAplot(h_mh, h_mA);

  // branching ratios from the file

  TH2F *h_brZZ_H = (TH2F*) f_xsecs->Get("h_brZZ_H");
  //TH2F *h_brbb_h = (TH2F*) f_xsecs->Get("h_brbb_h");

  // production cross sections
  //TH2F *h_ggF_xsec_A = (TH2F*) f_xsecs->Get("h_ggF_xsec_A");
  TH2F *h_ggF_xsec_H = (TH2F*) f_xsecs->Get("h_ggF_xsec_H");

  // gluon-fusion H->ZZ
  TH2F h_ggHZZ  = TH2F(*h_ggF_xsec_H);
  multiply_th2f(h_ggF_xsec_H,  h_brZZ_H, &h_ggHZZ);


  //setupLimits_ZZ();
  TFile rf_lims("graphs_combined_final_ggF.root");
  auto g_zz = (TGraphAsymmErrors*) rf_lims.Get("obs");
  auto g_zz_1 = (TGraphAsymmErrors*) rf_lims.Get("exp1s");
  auto g_zz_2 = (TGraphAsymmErrors*) rf_lims.Get("exp2s");

  TFile *foutput1 = new TFile(ofilename,"recreate");
  
  cout << "ZZ limits test ---- MASS method:" << endl;
  // hin, hmass, tag, obs, exp1s, exp2s, fout, folder, upper_mass
  bool use_bin_centre=true;
  if (tfilename.Contains("hMSSM")) use_bin_centre=false;
  makeLimitPlotMass(h_ggF_xsec_H, h_brZZ_H, h_mH, "gg->H->ZZ", g_zz, g_zz_1, g_zz_2, foutput1, "HZZ", upper_mass, 0, 0, use_bin_centre);
  
  
  foutput1->Close();
  
}
