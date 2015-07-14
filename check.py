import ROOT as R

rf = R.TFile("root-files/hmssm.root")
g = rf.Get("g2d_p2s")

x, y, z = [g.GetX(), g.GetY(), g.GetZ()]
for i in range(1, g.GetN()+1):
    print x[i], y[i], z[i]
