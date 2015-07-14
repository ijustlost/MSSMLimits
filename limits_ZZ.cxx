#include "limitCalc.h"

int main() {
    get_lims("out.mhmodm-8TeV-tanbHigh-nnlo.root", "root-files/mhmodm.root", 400 );
    get_lims("out.mhmodp-8TeV-tanbHigh-nnlo.root", "root-files/mhmodp.root" , 400);
    get_lims("out.lightstopmod-8TeV-tanbHigh-nnlo.root", "root-files/lightstop.root", 400 );
    get_lims("out.lightstau1-8TeV-tanbHigh-nnlo.root", "root-files/lightstau.root", 400 );
    get_lims("hMSSM_8TeV_1000_May18.root", "root-files/hmssm.root", 400 );
    get_lims("tauphobic_8TeV_tanbHigh_ataueqat_mu2000.root", "root-files/tauphobic_at.root", 400 );
    get_lims("tauphobic_8TeV_tanbHigh_ataueq0_mu2000.root", "root-files/tauphobic0.root", 400 );
}
