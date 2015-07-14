import ROOT as R
import sys

#rf = R.TFile("tauphobic_8TeV_tanbHigh_ataueqat_mu2000.root")
rf = R.TFile("out.lightstopmod-8TeV-tanbHigh-nnlo.root")
#rf = R.TFile("out.mhmodp-8TeV-tanbHigh-nnlo.root")
#rf = R.TFile("out.lightstau1-8TeV-tanbHigh-nnlo.root")

h_cx = rf.Get("h_ggF_xsec_H")
h_br = rf.Get("h_brZZ_H")

def print_cx_br(mass):
    iM = h_cx.GetXaxis().FindBin(mass)
    print mass
    for iTb in range(1, h_cx.GetYaxis().GetNbins()+1):
        tb_low = h_cx.GetYaxis().GetBinLowEdge(iTb)
        if tb_low>30: continue
        tb_high = tb_low+h_cx.GetYaxis().GetBinWidth(iTb)
        cx = h_cx.GetBinContent(iM, iTb)
        br = h_br.GetBinContent(iM, iTb)
        #print "{}~{}: {:.2g} ".format(tb_low, tb_high, cx, )
        print "{}~{}: {:.2g} {:.2g} {:.2g}".format(tb_low, tb_high, cx, br, cx*br)

print_cx_br(float(sys.argv[1]))
