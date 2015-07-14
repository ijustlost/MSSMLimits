import AtlasStyle
import ROOT as R
import PyROOTUtils
import AtlasUtil

R.gROOT.SetBatch()

#R.gStyle.SetHatchesLineWidth(10)

heap = []

def draw_graph(rf, color, style=):
    cut_obs = rf.Get("Obs_{}_excl".format(number))
    if not cut_obs:
        return False
    cut_obs.Draw("L")
    # format observed
    cut_obs.SetFillStyle(3636)
    #cut_obs.SetFillColor(R.kRed-7)
    cut_obs.SetFillColor(R.kRed-7)
    cut_obs.SetLineColor(R.kRed-7)
    cut_obs.SetLineWidth(3)
    #cut_obs.SetFillStyle(3004)
    cut_obs.Draw("F")
    heap.append(cut_obs)
    return True

def draw_blob(rf, number):
    # get objects
    cut_exp = rf.Get("Exp_{}_excl".format(number))
    cut_m1s = rf.Get("Minus1sigma_{}_excl".format(number))
    cut_m2s = rf.Get("Minus2sigma_{}_excl".format(number))
    cut_p1s = rf.Get("Plus1sigma_{}_excl".format(number))
    cut_p2s = rf.Get("Plus2sigma_{}_excl".format(number))

    if not cut_m2s:
        return False

    # format expected
    cut_exp.SetLineWidth(3)
    cut_exp.SetLineColor(R.kBlue)
    cut_exp.SetLineStyle(R.kDashed)

    # format bands
    cut_m2s.SetFillColor(R.kYellow)
    cut_m1s.SetFillColor(R.kGreen)
    cut_p1s.SetFillColor(R.kYellow)
    cut_p2s.SetFillColor(R.kWhite)

    cut_p2s.SetFillStyle(10001)
    cut_p1s.SetFillStyle(10001)
    cut_m1s.SetFillStyle(10001)
    cut_m2s.SetFillStyle(10001)

    # draw elements
    cut_m2s.Draw("F")
    cut_m1s.Draw("F")
    cut_p1s.Draw("F")
    cut_p2s.Draw("F")
    cut_exp.Draw("L")

    heap.append(cut_m2s)
    heap.append(cut_m1s)
    heap.append(cut_p1s)
    heap.append(cut_p2s)
    heap.append(cut_exp)

    return True



def draw( rfname, title, model, xmin, xmax, ymin, ymax, blob_numbers=[0,], blob_numbers_obs=[0,]):

    rf = R.TFile(rfname)

    # canvas setup
    #c = R.TCanvas("", "", 1100, 1300)
    c = R.TCanvas("", "", 600,600)
    c.SetLogy()

    # draw frame
    hframe = R.TH2D("frame", ";m_{A} [GeV];tan #beta", 1, xmin, xmax, 1, ymin, ymax)
    hframe.GetYaxis().SetRangeUser(ymin, ymax)
    hframe.Draw()

    #for blob in blob_numbers: draw_blob(rf, blob)
    #for blob in blob_numbers_obs: draw_blob_obs(rf, blob)
    i=0
    while True:
        if not draw_blob(rf, i): break
        i+=1
    i=0
    while True:
        if not draw_blob_obs(rf, i): break
        i+=1

    hframe.Draw("axis,same")

    # Decorations
    col_l=0.19
    top=0.85
    AtlasUtil.AtlasLabel(col_l, top, size=0.055)
    #AtlasUtil.AtlasLabel(col_l, top, size=0.06)
    AtlasUtil.DrawLuminosityFbEcmFirst(col_l, top-0.06, 20.3, size=0.045, sqrts=8, twoLines=False)
    AtlasUtil.DrawText(col_l, top-0.115, "#font[62]{H#rightarrowZZ}", size=0.045, font=62)
    #AtlasUtil.DrawText(col_l, top-0.11, "#font[12]{H#rightarrowZZ#rightarrowllll+ll#nu#nu+llqq+#nu#nuqq}", size=0.04)
    #AtlasUtil.DrawText(col_l, top-0.16, model, size=0.045, font=62)
    AtlasUtil.DrawText(col_l, top-0.17, model, size=0.045, font=62)

    # legend
    col_r = 0.55

    # pseudo-hist for obs
    g_0s = R.TH1D("", "", 1,0,1);
    g_0s.SetFillColor(R.kRed-7);
    g_0s.SetLineColor(R.kRed-7);
    g_0s.SetFillStyle(3636);
    # pseudo-hist for m1s
    g_exp = R.TH1D("", "", 1,0,1);
    #g_exp.SetFillColor(R.kGreen);
    g_exp.SetLineColor(R.kBlue);
    g_exp.SetLineStyle(R.kDashed);
    # pseudo-hist for m1s
    g_m1s = R.TH1D("", "", 1,0,1);
    g_m1s.SetFillColor(R.kGreen);
    g_m1s.SetLineColor(R.kBlack);
    g_m1s.SetFillStyle(1001);
    # pseudo-hist for m1s
    g_m2s = R.TH1D("", "", 1,0,1);
    g_m2s.SetFillColor(R.kYellow);
    g_m2s.SetLineColor(R.kBlack);
    g_m2s.SetFillStyle(1001);

    # make legend
    leg = R.TLegend(0.57, 0.62, 0.92, 0.9)
    #leg = R.TLegend(0.47, 0.7, 0.92, 0.9)
    #leg.SetNColumns(2)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetTextFont(62)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.AddEntry(g_0s, "Obs 95% CL limit", "l")
    leg.AddEntry(g_exp, "Exp 95% CL limit", "l")
    leg.AddEntry(g_0s, "Excluded", "f")
    leg.AddEntry(g_m1s, "#pm1#sigma band", "f")
    leg.AddEntry(g_m2s, "#pm2#sigma band", "f")
    leg.Draw()

    c.SaveAs(title+".pdf")
    c.SaveAs(title+".png")
    c.SaveAs(title+".eps")

#draw("root-files/mhmodm.root", "mhmodm", "mhmod-", 140, 400, 1, 20, [0,1,2,3], [0,1,2,3,4,5,6])
#draw("root-files/lightstau.root", "lightstau", "lightstau", 140, 400, 1, 20)
#draw("root-files/lightstop.root", "lightstop", "lightstop", 140, 400, 1, 20)
draw("root-files/hmssm.root", "hMSSM", "hMSSM", 141, 400, 1, 20)
#draw("root-files/tauphobic.root", "tauphobic", "tauphobic", 140, 400, 1, 20)
#draw("root-files/mhmodp.root", "mhmodp", "mhmod+", 140, 400, 1, 20)


