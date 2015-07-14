import AtlasStyle
import ROOT as R
import PyROOTUtils
import AtlasUtil
import os

R.gROOT.SetBatch()

#R.gStyle.SetHatchesLineWidth(10)

heap = []

# merge two TCutG. All of the points in cg2 are inserted to cg1 after
# the first point after mass_low
def merge_blobs(cg1, cg2, mass_low):
    cgnew = R.TCutG()
    cgnew.SetName(cg1.GetName()+"_"+cg2.GetName())
    m, tb = R.Double(), R.Double()
    # Add points from cg1 up to mass
    for iP in range(cg1.GetN()):
        cg1.GetPoint(iP, m, tb)
        print cgnew.GetN(), m, tb
        cgnew.SetPoint(cgnew.GetN(), m, tb)
        if m==mass_low:
            break
    # Add points from cg2 - but not last one
    for iP2 in range(cg2.GetN()-1):
        cg2.GetPoint(iP2, m, tb)
        print cgnew.GetN(), m, tb
        cgnew.SetPoint(cgnew.GetN(), m, tb)
    # Add points remaining cg1 points  (skip one because this is the extra one added at turnover)
    for iP in range(iP+2, cg1.GetN()):
        cg1.GetPoint(iP, m, tb)
        print cgnew.GetN(), m, tb
        cgnew.SetPoint(cgnew.GetN(), m, tb)
    return cgnew


def draw_blob_obs(rf, number, model):
    cut_obs = rf.Get("Obs_{}_excl".format(number))
    if not cut_obs:
        return False
    if model=="lightstau" and number==1:
        return
    if model=="lightstau" and number==0:
        cut_obs = merge_blobs(cut_obs, rf.Get("Obs_{}_excl".format(1)), 314)
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

def draw_graph2d(rf, canv, name, colour, contour_numbers, extra_points=None, line=None, fill=None):
    gr2d = rf.Get(name)
    #gr2d.SetNpx(260)
    #gr2d.SetNpy(80)
    ct = R.TCanvas()
    gr2d.Draw("ACOLZ")
    ct.Update()
    c = gr2d.GetContourList(0)
    for iC in contour_numbers:
        print "contour", iC
        gr = c.At(iC)
        if extra_points:
            if extra_points[iC]:
                gr.SetPoint(gr.GetN(), extra_points[iC][0], extra_points[iC][1] )
        canv.cd()
        if fill:
            gr.SetFillStyle(fill)
            gr.SetFillColor(colour)
            gr.Draw("F")
        if line:
            gr.Draw("L")
            gr.SetLineStyle(line)
            gr.SetLineWidth(3)
            gr.SetLineColor(colour)
        heap.append(gr)

def draw_blob(rf, number, model):
    # get objects
    cut_exp = rf.Get("Exp_{}_excl".format(number))
    cut_m1s = rf.Get("Minus1sigma_{}_excl".format(number))
    cut_m2s = rf.Get("Minus2sigma_{}_excl".format(number))
    cut_p1s = rf.Get("Plus1sigma_{}_excl".format(number))
    cut_p2s = rf.Get("Plus2sigma_{}_excl".format(number))

    if not cut_m2s and not cut_m1s and not cut_m2s and not cut_p1s and not cut_p2s:
        return False

    if model=="mhmodp" and number==0:
        cut_m2s = merge_blobs(cut_m2s, rf.Get("Minus2sigma_{}_excl".format(1)), 315)
    if model=="lightstau" and number==0:
        cut_m2s = merge_blobs(cut_m2s, rf.Get("Minus2sigma_{}_excl".format(1)), 314)
        cut_m1s = merge_blobs(cut_m1s, rf.Get("Minus1sigma_{}_excl".format(1)), 314)
        cut_exp = merge_blobs(cut_exp, rf.Get("Exp_{}_excl".format(1)), 314)
        #cut_p2s = merge_blobs(cut_p2s, rf.Get("Plus2sigma_{}_excl".format(1)), 314)
        #cut_p1s = merge_blobs(cut_p1s, rf.Get("Plus1sigma_{}_excl".format(1)), 314)


    # format bands
    if cut_m2s:
        if (model=="mhmodp" or model=="lightstau") and number==1:
            pass
        else:
            cut_m2s.SetFillStyle(10001)
            cut_m2s.SetFillColor(R.kYellow)
            cut_m2s.Draw("F")
            heap.append(cut_m2s)

    if cut_m1s:
        if (model=="lightstau") and number==1:
            pass
        else:
            cut_m1s.SetFillColor(R.kGreen)
            cut_m1s.SetFillStyle(10001)
            cut_m1s.Draw("F")
            heap.append(cut_m1s)

    if cut_p1s:
        cut_p1s.SetFillColor(R.kYellow)
        if (model=="lightstau") and number==1:
            pass
        else:
            cut_p1s.SetFillStyle(10001)
            cut_p1s.Draw("F")
            heap.append(cut_p1s)

    if cut_p2s:
        cut_p2s.SetFillColor(R.kWhite)
        print model
        if model=="tauphobic_at" and number==1: 
            cut_p2s.SetFillColor(R.kGreen)
        if (model=="lightstau") and number==1:
            pass
        else:
            cut_p2s.SetFillStyle(10001)
            cut_p2s.Draw("F")
            heap.append(cut_p2s)

    # format expected and draw
    if cut_exp:
        if (model=="lightstau") and number==1:
            pass
        else:
            cut_exp.SetLineWidth(3)
            cut_exp.SetLineColor(R.kBlue)
            cut_exp.SetLineStyle(R.kDashed)
            cut_exp.Draw("L")
            heap.append(cut_exp)


    return True



def draw( rfname, title, model, xmin, xmax, ymin, ymax, blob_numbers=[0,], blob_numbers_obs=[0,], mode=0):

    rf = R.TFile(rfname)

    # canvas setup
    #c = R.TCanvas("", "", 1100, 1300)
    c = R.TCanvas("", "", 600,600)
    c.SetLogy()

    # draw frame
    hframe = R.TH2D("frame", ";m_{A} [GeV];tan #beta", 1, xmin, xmax, 1, ymin, ymax)
    hframe.GetYaxis().SetRangeUser(ymin, ymax)
    hframe.GetXaxis().SetNdivisions(9,5,2)
    hframe.Draw()

    #for blob in blob_numbers: draw_blob(rf, blob)
    #for blob in blob_numbers_obs: draw_blob_obs(rf, blob)
    if mode==0:
        i=0
        while True:
            if not draw_blob(rf, i, title): break
            i+=1
        i=0
        while True:
            if not draw_blob_obs(rf, i, title): break
            i+=1
    else:
        draw_graph2d(rf, c, "g2d_m2s", R.kYellow, [0,], [(140,ymin)], fill=10001)
        draw_graph2d(rf, c, "g2d_m1s", R.kGreen, [0,], [(140,ymin)], fill=10001)
        draw_graph2d(rf, c, "g2d_p1s", R.kYellow, [0,], [(140,ymin)], fill=10001)
        draw_graph2d(rf, c, "g2d_p2s", R.kWhite, [0,], [(140,ymin)], fill=10001)
        draw_graph2d(rf, c, "g2d_exp", R.kBlue, [0,], line=R.kDashed)
        draw_graph2d(rf, c, "g2d_obs", R.kRed-7, [0,], [(140,ymin)], fill=3663, line=1)

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

    c.SaveAs("plots/"+tag+"/"+title+".pdf")
    c.SaveAs("plots/"+tag+"/"+title+".png")
    c.SaveAs("plots/"+tag+"/"+title+".eps")


tag = "cutg_0p5_cut1000_full"
mintb = 1
maxtb=15
mmax=330
os.system("mkdir -p plots/"+tag)

draw("root-files/mhmodm.root", "mhmodm", "mhmod-", 140, mmax, mintb, maxtb, [0,1,2,3], [0,1,2,3,4,5,6])
draw("root-files/lightstau.root", "lightstau", "Light Stau", 140, mmax, mintb, maxtb)
draw("root-files/lightstop.root", "lightstop", "Light Stop", 140, mmax, mintb, maxtb)
draw("root-files/hmssm.root", "hMSSM", "hMSSM", 160, mmax, mmax, maxtb)
draw("root-files/tauphobic_0.root", "tauphobic_0", "#tau-phobic A_{#tau}=0", 140, mmax, 0.99, maxtb)
draw("root-files/tauphobic_at.root", "tauphobic_at", "#tau-phobic A_{#tau}=A_{t}", 140, mmax, 0.99, maxtb)
draw("root-files/mhmodp.root", "mhmodp", "mhmod+", 140, mmax, mintb, maxtb)


