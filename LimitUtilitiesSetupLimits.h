#include "TGraphAsymmErrors.h"


TGraphAsymmErrors *g_tata  ;//= new TGraphAsymmErrors();
TGraphAsymmErrors *g_tata_0;//= new TGraphAsymmErrors();
TGraphAsymmErrors *g_tata_1;//= new TGraphAsymmErrors();
TGraphAsymmErrors *g_tata_2;//= new TGraphAsymmErrors();

TGraphAsymmErrors *g_zz  ;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zz_0;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zz_1;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zz_2;// = new TGraphAsymmErrors();

TGraphAsymmErrors *g_zh  ;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_0;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_1;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_2;// = new TGraphAsymmErrors();

TGraphAsymmErrors *g_zh_ta  ;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_ta_0;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_ta_1;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_ta_2;// = new TGraphAsymmErrors();

TGraphAsymmErrors *g_zh_bbta  ;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_bbta_0;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_bbta_1;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_zh_bbta_2;// = new TGraphAsymmErrors();

TGraphAsymmErrors *g_ww_gf  ;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_ww_gf_0;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_ww_gf_1;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_ww_gf_2;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_ww_vbf  ;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_ww_vbf_0;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_ww_vbf_1;// = new TGraphAsymmErrors();
TGraphAsymmErrors *g_ww_vbf_2;// = new TGraphAsymmErrors();

void setupLimits_WW() {
  g_ww_gf   = new TGraphAsymmErrors();
  g_ww_gf_0 = new TGraphAsymmErrors();
  g_ww_gf_1 = new TGraphAsymmErrors();
  g_ww_gf_2 = new TGraphAsymmErrors();

  // mass obs exp, gluon fusion
  //  300  0.906925   1.07369 
  //  400  0.350416   0.449319 
  //  500  0.193205   0.21124 
  //  600  0.119645   0.127064 
  //  700  0.047786   0.0826859 
  //  800  0.029902   0.0578101 
  //  900  0.03264    0.0439012 
  // 1000 0.0195125  0.0298902
  // 1100 0.0183116  0.0244631
  // 1200 0.0262699  0.0230397
  // 1300 0.0253825  0.0185742
  // 1400 0.0229941  0.0170418
  // 1500 0.0207188  0.016298 

  double vm[5] = { 299, 300, 400, 500, 600 };
  double v[5]  = {10000.,0.906925, 0.350416, 0.193205, 0.119645};
  double v0[5] = {10000.,1.07369,  0.449319, 0.21124, 0.127064};
   for (int i=0; i<5; ++i) {
    g_ww_gf->SetPoint(i, vm[i], v[i]);  
    g_ww_gf_0->SetPoint(i, vm[i], v0[i]); 
    g_ww_gf_1->SetPoint(i, vm[i], v0[i]); 
    g_ww_gf_2->SetPoint(i, vm[i], v0[i]);
    g_ww_gf_1->SetPointError(i, 0,0,   v0[i]*0.1, v0[i]*0.10);
    g_ww_gf_2->SetPointError(i, 0,0,   v0[i]*0.1, v0[i]*0.10);
  }

  // VBF
  // 
  // 300  0.230469   0.213375 
  // 400  0.14667    0.112445 
  // 500  0.111352   0.0686241 
  // 600  0.0806416  0.0426829 
  // 700  0.0340431  0.0281964 
  // 800  0.0236066  0.0203046 
  // 900  0.0187602  0.0155217 
  // 1000 0.0136101  0.0107314 
  // 1100 0.00994469 0.0108211 
  // 1200 0.00765705 0.00856942
  // 1300 0.00661648 0.00746784
  // 1400 0.00655076 0.00720692
  // 1500 0.00597646 0.00700654

  g_ww_vbf   = new TGraphAsymmErrors();
  g_ww_vbf_0 = new TGraphAsymmErrors();
  g_ww_vbf_1 = new TGraphAsymmErrors();
  g_ww_vbf_2 = new TGraphAsymmErrors();

  double v1m[4] = { 300, 400, 500, 600 };
  double v1[4]  = {0.230469,0.14667, 0.111352, 0.0806416};
  double v10[4] = {0.213375,0.112445 ,0.0686241,0.0426829};

   for (int i=0; i<4; ++i) {
    g_ww_vbf->SetPoint(i, v1m[i], v1[i]);  
    g_ww_vbf_0->SetPoint(i, v1m[i], v10[i]); 
    g_ww_vbf_1->SetPoint(i, v1m[i], v10[i]); 
    g_ww_vbf_2->SetPoint(i, v1m[i], v10[i]);
    g_ww_vbf_1->SetPointError(i, 0,0,   v10[i]*0.1, v10[i]*0.10);
    g_ww_vbf_2->SetPointError(i, 0,0,   v10[i]*0.1, v10[i]*0.10);
  }

}



void setupLimits_AZh() {

  g_zh   = new TGraphAsymmErrors();
  g_zh_0 = new TGraphAsymmErrors();
  g_zh_1 = new TGraphAsymmErrors();
  g_zh_2 = new TGraphAsymmErrors();


  double vm[11] = {      219, 220,  240,   260,   300 , 340 , 350, 400 ,  450 ,  500 };
  double v[11]  = {10000000.,0.34 ,0.35 , 0.31 , 0.23 , 0.16,0.14, 0.092, 0.065 ,0.046};
  double v0[11] = {10000000.,0.57 ,0.29 , 0.38 , 0.21 , 0.18,0.17, 0.15 , 0.043 ,0.034};
 

  for (int i=0; i<11; ++i) {
    g_zh->SetPoint(i, vm[i], v[i]);  
    g_zh_0->SetPoint(i, vm[i], v0[i]); g_zh_1->SetPoint(i, vm[i], v0[i]); g_zh_2->SetPoint(i, vm[i], v0[i]);
    g_zh_1->SetPointError(i, 0,0,   v0[i]*0.1, v0[i]*0.10);
    g_zh_2->SetPointError(i, 0,0,   v0[i]*0.1, v0[i]*0.10);
  }

  // tau+bb limits assuming SM BRs
  g_zh_ta   = new TGraphAsymmErrors();
  g_zh_ta_0 = new TGraphAsymmErrors();
  g_zh_ta_1 = new TGraphAsymmErrors();
  g_zh_ta_2 = new TGraphAsymmErrors();
  
  double vt1m[9] = {  219, 220,  240,   260,    300 , 340 , 350,  400 ,  500 };
  double vt1[9]  = {10000000.,   0.098 , 0.079 , 0.12, 0.066 , 0.039 , 0.038 , 0.042 , 0.018};
  double vt10[9] = {10000000.,   0.11, 0.094 , 0.083 , 0.070 , 0.052 , 0.049 , 0.047 , 0.031};

  for (int i=0; i<9; ++i) {
    g_zh_ta->SetPoint(i, vt1m[i], vt1[i]);  
    g_zh_ta_0->SetPoint(i, vt1m[i], vt10[i]); 
    g_zh_ta_1->SetPoint(i, vt1m[i], vt10[i]); g_zh_ta_2->SetPoint(i, vt1m[i], vt10[i]);
    g_zh_ta_1->SetPointError(i, 0,0,   vt10[i]*0.1, vt10[i]*0.10);
    g_zh_ta_2->SetPointError(i, 0,0,   vt10[i]*0.1, vt10[i]*0.10);
  }


  // tau+bb limits assuming SM BRs
  g_zh_bbta   = new TGraphAsymmErrors();
  g_zh_bbta_0 = new TGraphAsymmErrors();
  g_zh_bbta_1 = new TGraphAsymmErrors();
  g_zh_bbta_2 = new TGraphAsymmErrors();
  
  double v1m[9] = {  219, 220,  240,   260,    300 , 340 , 350,  400 ,  500 };
  double v1[9]  = {10000000., 0.89 , 0.45 , 0.64 , 0.30 , 0.26 , 0.23 , 0.23, 0.049};
  double v10[9] = {10000000.,  0.54, 0.52 , 0.47 , 0.35 , 0.24 , 0.22 , 0.15, 0.076};

  for (int i=0; i<9; ++i) {
    g_zh_bbta->SetPoint(i, v1m[i], v1[i]);  
    g_zh_bbta_0->SetPoint(i, v1m[i], v10[i]); 
    g_zh_bbta_1->SetPoint(i, v1m[i], v10[i]); g_zh_bbta_2->SetPoint(i, v1m[i], v10[i]);
    g_zh_bbta_1->SetPointError(i, 0,0,   v10[i]*0.1, v10[i]*0.10);
    g_zh_bbta_2->SetPointError(i, 0,0,   v10[i]*0.1, v10[i]*0.10);
  }



}


// cross section limits for ZZ
//  cross section limits, gluon fusion, high mass
//
void setupLimits_ZZ() {
  



  g_zz   = new TGraphAsymmErrors();
  g_zz_0 = new TGraphAsymmErrors();
  g_zz_1 = new TGraphAsymmErrors();
  g_zz_2 = new TGraphAsymmErrors();
  

  g_zz->SetPoint(0, 200, 0.32428);  
g_zz_0->SetPoint(0, 200, 0.35979); g_zz_1->SetPoint(0, 200, 0.35979); g_zz_2->SetPoint(0, 200, 0.35979);
g_zz_1->SetPointError(0, 0,0,   0.35979 - 0.23366, 0.46334 -  0.35979);
g_zz_2->SetPointError(0, 0,0,   0.35979 - 0.17405, 0.64839 -  0.35979);
g_zz->SetPoint(1, 220, 0.2964);  
g_zz_0->SetPoint(1, 220, 0.36228); g_zz_1->SetPoint(1, 220, 0.36228); g_zz_2->SetPoint(1, 220, 0.36228);
g_zz_1->SetPointError(1, 0,0,   0.36228 - 0.21358, 0.42263 -  0.36228);
g_zz_2->SetPointError(1, 0,0,   0.36228 - 0.15909, 0.59386 -  0.36228);
g_zz->SetPoint(2, 240, 0.2643);  
g_zz_0->SetPoint(2, 240, 0.24713); g_zz_1->SetPoint(2, 240, 0.24713); g_zz_2->SetPoint(2, 240, 0.24713);
g_zz_1->SetPointError(2, 0,0,   0.24713 - 0.19044, 0.37564 -  0.24713);
g_zz_2->SetPointError(2, 0,0,   0.24713 - 0.14186, 0.52888 -  0.24713);
g_zz->SetPoint(3, 260, 0.22584);  
g_zz_0->SetPoint(3, 260, 0.34451); g_zz_1->SetPoint(3, 260, 0.34451); g_zz_2->SetPoint(3, 260, 0.34451);
g_zz_1->SetPointError(3, 0,0,   0.34451 - 0.16273, 0.3222 -  0.34451);
g_zz_2->SetPointError(3, 0,0,   0.34451 - 0.12121, 0.45553 -  0.34451);
g_zz->SetPoint(4, 280, 0.19552);  
g_zz_0->SetPoint(4, 280, 0.22516); g_zz_1->SetPoint(4, 280, 0.22516); g_zz_2->SetPoint(4, 280, 0.22516);
g_zz_1->SetPointError(4, 0,0,   0.22516 - 0.14088, 0.27889 -  0.22516);
g_zz_2->SetPointError(4, 0,0,   0.22516 - 0.10494, 0.39304 -  0.22516);
g_zz->SetPoint(5, 300, 0.17403);  
g_zz_0->SetPoint(5, 300, 0.25519); g_zz_1->SetPoint(5, 300, 0.25519); g_zz_2->SetPoint(5, 300, 0.25519);
g_zz_1->SetPointError(5, 0,0,   0.25519 - 0.1254, 0.24831 -  0.25519);
g_zz_2->SetPointError(5, 0,0,   0.25519 - 0.09341, 0.35006 -  0.25519);
g_zz->SetPoint(6, 320, 0.1518);  
g_zz_0->SetPoint(6, 320, 0.09262); g_zz_1->SetPoint(6, 320, 0.09262); g_zz_2->SetPoint(6, 320, 0.09262);
g_zz_1->SetPointError(6, 0,0,   0.09262 - 0.10938, 0.21654 -  0.09262);
g_zz_2->SetPointError(6, 0,0,   0.09262 - 0.08148, 0.30498 -  0.09262);
g_zz->SetPoint(7, 340, 0.13109);  
g_zz_0->SetPoint(7, 340, 0.10741); g_zz_1->SetPoint(7, 340, 0.10741); g_zz_2->SetPoint(7, 340, 0.10741);
g_zz_1->SetPointError(7, 0,0,   0.10741 - 0.09446, 0.18667 -  0.10741);
g_zz_2->SetPointError(7, 0,0,   0.10741 - 0.07036, 0.2622 -  0.10741);
g_zz->SetPoint(8, 360, 0.11659);  
g_zz_0->SetPoint(8, 360, 0.0987); g_zz_1->SetPoint(8, 360, 0.0987); g_zz_2->SetPoint(8, 360, 0.0987);
g_zz_1->SetPointError(8, 0,0,   0.0987 - 0.08401, 0.16644 -  0.0987);
g_zz_2->SetPointError(8, 0,0,   0.0987 - 0.06258, 0.23356 -  0.0987);
g_zz->SetPoint(9, 380, 0.1001);  
g_zz_0->SetPoint(9, 380, 0.12462); g_zz_1->SetPoint(9, 380, 0.12462); g_zz_2->SetPoint(9, 380, 0.12462);
g_zz_1->SetPointError(9, 0,0,   0.12462 - 0.07213, 0.14281 -  0.12462);
g_zz_2->SetPointError(9, 0,0,   0.12462 - 0.05373, 0.20025 -  0.12462);
g_zz->SetPoint(10, 400, 0.08458);  
g_zz_0->SetPoint(10, 400, 0.06498); g_zz_1->SetPoint(10, 400, 0.06498); g_zz_2->SetPoint(10, 400, 0.06498);
g_zz_1->SetPointError(10, 0,0,   0.06498 - 0.06094, 0.12056 -  0.06498);
g_zz_2->SetPointError(10, 0,0,   0.06498 - 0.0454, 0.16869 -  0.06498);
g_zz->SetPoint(11, 420, 0.07216);  
g_zz_0->SetPoint(11, 420, 0.06054); g_zz_1->SetPoint(11, 420, 0.06054); g_zz_2->SetPoint(11, 420, 0.06054);
g_zz_1->SetPointError(11, 0,0,   0.06054 - 0.05199, 0.10269 -  0.06054);
g_zz_2->SetPointError(11, 0,0,   0.06054 - 0.03873, 0.14384 -  0.06054);
g_zz->SetPoint(12, 440, 0.06505);  
g_zz_0->SetPoint(12, 440, 0.07115); g_zz_1->SetPoint(12, 440, 0.07115); g_zz_2->SetPoint(12, 440, 0.07115);
g_zz_1->SetPointError(12, 0,0,   0.07115 - 0.04687, 0.0926 -  0.07115);
g_zz_2->SetPointError(12, 0,0,   0.07115 - 0.03491, 0.12934 -  0.07115);
g_zz->SetPoint(13, 460, 0.05789);  
g_zz_0->SetPoint(13, 460, 0.04842); g_zz_1->SetPoint(13, 460, 0.04842); g_zz_2->SetPoint(13, 460, 0.04842);
g_zz_1->SetPointError(13, 0,0,   0.04842 - 0.04171, 0.08238 -  0.04842);
g_zz_2->SetPointError(13, 0,0,   0.04842 - 0.03107, 0.11507 -  0.04842);
g_zz->SetPoint(14, 480, 0.05269);  
g_zz_0->SetPoint(14, 480, 0.04178); g_zz_1->SetPoint(14, 480, 0.04178); g_zz_2->SetPoint(14, 480, 0.04178);
g_zz_1->SetPointError(14, 0,0,   0.04178 - 0.03797, 0.07503 -  0.04178);
g_zz_2->SetPointError(14, 0,0,   0.04178 - 0.02828, 0.10502 -  0.04178);
g_zz->SetPoint(15, 500, 0.04708);  
g_zz_0->SetPoint(15, 500, 0.04432); g_zz_1->SetPoint(15, 500, 0.04432); g_zz_2->SetPoint(15, 500, 0.04432);
g_zz_1->SetPointError(15, 0,0,   0.04432 - 0.03393, 0.0671 -  0.04432);
g_zz_2->SetPointError(15, 0,0,   0.04432 - 0.02527, 0.09404 -  0.04432);
g_zz->SetPoint(16, 520, 0.04272);  
g_zz_0->SetPoint(16, 520, 0.03905); g_zz_1->SetPoint(16, 520, 0.03905); g_zz_2->SetPoint(16, 520, 0.03905);
g_zz_1->SetPointError(16, 0,0,   0.03905 - 0.03078, 0.06095 -  0.03905);
g_zz_2->SetPointError(16, 0,0,   0.03905 - 0.02293, 0.08565 -  0.03905);
g_zz->SetPoint(17, 540, 0.03904);  
g_zz_0->SetPoint(17, 540, 0.03223); g_zz_1->SetPoint(17, 540, 0.03223); g_zz_2->SetPoint(17, 540, 0.03223);
g_zz_1->SetPointError(17, 0,0,   0.03223 - 0.02813, 0.05582 -  0.03223);
g_zz_2->SetPointError(17, 0,0,   0.03223 - 0.02095, 0.07863 -  0.03223);
g_zz->SetPoint(18, 560, 0.03665);  
g_zz_0->SetPoint(18, 560, 0.02969); g_zz_1->SetPoint(18, 560, 0.02969); g_zz_2->SetPoint(18, 560, 0.02969);
g_zz_1->SetPointError(18, 0,0,   0.02969 - 0.02641, 0.05236 -  0.02969);
g_zz_2->SetPointError(18, 0,0,   0.02969 - 0.01967, 0.07377 -  0.02969);
g_zz->SetPoint(19, 580, 0.03395);  
g_zz_0->SetPoint(19, 580, 0.02531); g_zz_1->SetPoint(19, 580, 0.02531); g_zz_2->SetPoint(19, 580, 0.02531);
g_zz_1->SetPointError(19, 0,0,   0.02531 - 0.02446, 0.04855 -  0.02531);
g_zz_2->SetPointError(19, 0,0,   0.02531 - 0.01822, 0.06856 -  0.02531);
g_zz->SetPoint(20, 600, 0.0321);  
g_zz_0->SetPoint(20, 600, 0.02165); g_zz_1->SetPoint(20, 600, 0.02165); g_zz_2->SetPoint(20, 600, 0.02165);
g_zz_1->SetPointError(20, 0,0,   0.02165 - 0.02313, 0.04595 -  0.02165);
g_zz_2->SetPointError(20, 0,0,   0.02165 - 0.01723, 0.06501 -  0.02165);
g_zz->SetPoint(21, 650, 0.02698);  
g_zz_0->SetPoint(21, 650, 0.02647); g_zz_1->SetPoint(21, 650, 0.02647); g_zz_2->SetPoint(21, 650, 0.02647);
g_zz_1->SetPointError(21, 0,0,   0.02647 - 0.01944, 0.03867 -  0.02647);
g_zz_2->SetPointError(21, 0,0,   0.02647 - 0.01448, 0.05481 -  0.02647);
g_zz->SetPoint(22, 700, 0.02251);  
g_zz_0->SetPoint(22, 700, 0.0199); g_zz_1->SetPoint(22, 700, 0.0199); g_zz_2->SetPoint(22, 700, 0.0199);
g_zz_1->SetPointError(22, 0,0,   0.0199 - 0.01622, 0.03245 -  0.0199);
g_zz_2->SetPointError(22, 0,0,   0.0199 - 0.01208, 0.04636 -  0.0199);
g_zz->SetPoint(23, 750, 0.01948);  
g_zz_0->SetPoint(23, 750, 0.01169); g_zz_1->SetPoint(23, 750, 0.01169); g_zz_2->SetPoint(23, 750, 0.01169);
g_zz_1->SetPointError(23, 0,0,   0.01169 - 0.01404, 0.02817 -  0.01169);
g_zz_2->SetPointError(23, 0,0,   0.01169 - 0.01046, 0.04041 -  0.01169);
g_zz->SetPoint(24, 800, 0.01716);  
g_zz_0->SetPoint(24, 800, 0.01198); g_zz_1->SetPoint(24, 800, 0.01198); g_zz_2->SetPoint(24, 800, 0.01198);
g_zz_1->SetPointError(24, 0,0,   0.01198 - 0.01236, 0.02481 -  0.01198);
g_zz_2->SetPointError(24, 0,0,   0.01198 - 0.00921, 0.03561 -  0.01198);
g_zz->SetPoint(25, 850, 0.01536);  
g_zz_0->SetPoint(25, 850, 0.01131); g_zz_1->SetPoint(25, 850, 0.01131); g_zz_2->SetPoint(25, 850, 0.01131);
g_zz_1->SetPointError(25, 0,0,   0.01131 - 0.01107, 0.02225 -  0.01131);
g_zz_2->SetPointError(25, 0,0,   0.01131 - 0.00824, 0.03211 -  0.01131);
g_zz->SetPoint(26, 900, 0.01367);  
g_zz_0->SetPoint(26, 900, 0.01388); g_zz_1->SetPoint(26, 900, 0.01388); g_zz_2->SetPoint(26, 900, 0.01388);
g_zz_1->SetPointError(26, 0,0,   0.01388 - 0.00985, 0.01985 -  0.01388);
g_zz_2->SetPointError(26, 0,0,   0.01388 - 0.00733, 0.02882 -  0.01388);
g_zz->SetPoint(27, 950, 0.01264);  
g_zz_0->SetPoint(27, 950, 0.00762); g_zz_1->SetPoint(27, 950, 0.00762); g_zz_2->SetPoint(27, 950, 0.00762);
g_zz_1->SetPointError(27, 0,0,   0.00762 - 0.00911, 0.01842 -  0.00762);
g_zz_2->SetPointError(27, 0,0,   0.00762 - 0.00679, 0.02697 -  0.00762);
g_zz->SetPoint(28, 1000, 0.01166);  
g_zz_0->SetPoint(28, 1000, 0.01088); g_zz_1->SetPoint(28, 1000, 0.01088); g_zz_2->SetPoint(28, 1000, 0.01088);
g_zz_1->SetPointError(28, 0,0,   0.01088 - 0.0084, 0.01704 -  0.01088);
g_zz_2->SetPointError(28, 0,0,   0.01088 - 0.00626, 0.02507 -  0.01088);


}

//  cross section limits, gluon fusion, high mass
//
void setupLimits_tata() {

  g_tata  = new TGraphAsymmErrors();
  g_tata_0= new TGraphAsymmErrors();
  g_tata_1= new TGraphAsymmErrors();
  g_tata_2= new TGraphAsymmErrors();
  

 // Summary: 200 0.8823 0.7435 1.4126 1.0397 0.5358 0.3991
 g_tata->SetPoint(0, 200,  0.8823);
 g_tata_0->SetPoint(0, 200, 0.7435); g_tata_1->SetPoint(0, 200, 0.7435);g_tata_2->SetPoint(0, 200, 0.7435);
 g_tata_1->SetPointError(0, 0,0, 0.7435-0.5358, 1.0397-0.7435);
 g_tata_2->SetPointError(0, 0,0, 0.7435-0.3991, 1.4126-0.7435);
 // Summary: 250 0.2798 0.3516 0.6702 0.4925 0.2534 0.1887
 g_tata->SetPoint(1, 250,  0.2798);
 g_tata_0->SetPoint(1, 250, 0.3516); g_tata_1->SetPoint(1, 250, 0.3516);g_tata_2->SetPoint(1, 250, 0.3516);
 g_tata_1->SetPointError(1, 0,0, 0.3516-0.2534, 0.4925-0.3516);
 g_tata_2->SetPointError(1, 0,0, 0.3516-0.1887, 0.6702-0.3516);
 // Summary: 300 0.1358 0.1980 0.3789 0.2781 0.1427 0.1063
 g_tata->SetPoint(2, 300,  0.1358);
 g_tata_0->SetPoint(2, 300, 0.1980); g_tata_1->SetPoint(2, 300, 0.1980);g_tata_2->SetPoint(2, 300, 0.1980);
 g_tata_1->SetPointError(2, 0,0, 0.1980-0.1427, 0.2781-0.1980);
 g_tata_2->SetPointError(2, 0,0, 0.1980-0.1063, 0.3789-0.1980);
 // Summary: 350 0.1250 0.1156 0.2276 0.1646 0.0833 0.0620
 g_tata->SetPoint(3, 350,  0.1250);
 g_tata_0->SetPoint(3, 350, 0.1156); g_tata_1->SetPoint(3, 350, 0.1156);g_tata_2->SetPoint(3, 350, 0.1156);
 g_tata_1->SetPointError(3, 0,0, 0.1156-0.0833, 0.1646-0.1156);
 g_tata_2->SetPointError(3, 0,0, 0.1156-0.0620, 0.2276-0.1156);
 // Summary: 400 0.0863 0.0766 0.1509 0.1091 0.0552 0.0411
 g_tata->SetPoint(4, 400,  0.0863);
 g_tata_0->SetPoint(4, 400, 0.0766); g_tata_1->SetPoint(4, 400, 0.0766);g_tata_2->SetPoint(4, 400, 0.0766);
 g_tata_1->SetPointError(4, 0,0, 0.0766-0.0552, 0.1091-0.0766);
 g_tata_2->SetPointError(4, 0,0, 0.0766-0.0411, 0.1509-0.0766);
 // Summary: 450 0.0447 0.0542 0.1075 0.0774 0.0390 0.0291
 g_tata->SetPoint(5, 450,  0.0447);
 g_tata_0->SetPoint(5, 450, 0.0542); g_tata_1->SetPoint(5, 450, 0.0542);g_tata_2->SetPoint(5, 450, 0.0542);
 g_tata_1->SetPointError(5, 0,0, 0.0542-0.0390, 0.0774-0.0542);
 g_tata_2->SetPointError(5, 0,0, 0.0542-0.0291, 0.1075-0.0542);
 // Summary: 500 0.0262 0.0395 0.0789 0.0565 0.0284 0.0212
 g_tata->SetPoint(6, 500,  0.0262);
 g_tata_0->SetPoint(6, 500, 0.0395); g_tata_1->SetPoint(6, 500, 0.0395);g_tata_2->SetPoint(6, 500, 0.0395);
 g_tata_1->SetPointError(6, 0,0, 0.0395-0.0284, 0.0565-0.0395);
 g_tata_2->SetPointError(6, 0,0, 0.0395-0.0212, 0.0789-0.0395);
 // Summary: 600 0.0194 0.0259 0.0518 0.0371 0.0187 0.0139
 g_tata->SetPoint(7, 600,  0.0194);
 g_tata_0->SetPoint(7, 600, 0.0259); g_tata_1->SetPoint(7, 600, 0.0259);g_tata_2->SetPoint(7, 600, 0.0259);
 g_tata_1->SetPointError(7, 0,0, 0.0259-0.0187, 0.0371-0.0259);
 g_tata_2->SetPointError(7, 0,0, 0.0259-0.0139, 0.0518-0.0259);
 // Summary: 700 0.0137 0.0173 0.0350 0.0249 0.0125 0.0093
 g_tata->SetPoint(8, 700,  0.0137);
 g_tata_0->SetPoint(8, 700, 0.0173); g_tata_1->SetPoint(8, 700, 0.0173);g_tata_2->SetPoint(8, 700, 0.0173);
 g_tata_1->SetPointError(8, 0,0, 0.0173-0.0125, 0.0249-0.0173);
 g_tata_2->SetPointError(8, 0,0, 0.0173-0.0093, 0.0350-0.0173);
 // Summary: 800 0.0104 0.0129 0.0263 0.0186 0.0093 0.0069
 g_tata->SetPoint(9, 800,  0.0104);
 g_tata_0->SetPoint(9, 800, 0.0129); g_tata_1->SetPoint(9, 800, 0.0129);g_tata_2->SetPoint(9, 800, 0.0129);
 g_tata_1->SetPointError(9, 0,0, 0.0129-0.0093, 0.0186-0.0129);
 g_tata_2->SetPointError(9, 0,0, 0.0129-0.0069, 0.0263-0.0129);
 // Summary: 900 0.0088 0.0108 0.0219 0.0155 0.0078 0.0058
 g_tata->SetPoint(10, 900,  0.0088);
 g_tata_0->SetPoint(10, 900, 0.0108); g_tata_1->SetPoint(10, 900, 0.0108);g_tata_2->SetPoint(10, 900, 0.0108);
 g_tata_1->SetPointError(10, 0,0, 0.0108-0.0078, 0.0155-0.0108);
 g_tata_2->SetPointError(10, 0,0, 0.0108-0.0058, 0.0219-0.0108);
 // Summary: 1000 0.0076 0.0089 0.0178 0.0127 0.0064 0.0048
 g_tata->SetPoint(11, 1000,  0.0076);
 g_tata_0->SetPoint(11, 1000, 0.0089); g_tata_1->SetPoint(11, 1000, 0.0089);g_tata_2->SetPoint(11, 1000, 0.0089);
 g_tata_1->SetPointError(11, 0,0, 0.0089-0.0064, 0.0127-0.0089);
 g_tata_2->SetPointError(11, 0,0, 0.0089-0.0048, 0.0178-0.0089);



}
