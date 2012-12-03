#!/usr/bin/env python

from ROOT import *
from TopAnalysis.TTbarDilepton.PlotTool import *

gROOT.ProcessLine(".x rootlogon.C")

mcProdName = "Summer12"
realDataName = "Run2012"
mode = "mm"
catOS = "m_OS"
catSS = "m_SS"
histName = "hMuon2_relIso"
#histName = "hMuon1_relIso"
#histName = "hElectron1_mva"

## Real data
p = PlotBuilder(realDataName, mcProdName, mode, "hEvents")
hRDOS = p.buildRDHist(catOS, histName)
hRDSS = p.buildRDHist(catSS, histName)
hMCOS, hMCOSs = p.buildMCHists(catOS, histName)
hMCSS, hMCSSs = p.buildMCHists(catSS, histName)

hQCDOS = hRDOS.Clone("hQCDOS")
hQCDSS = hRDSS.Clone("hQCDSS")

for name, h in hMCOSs:
    hQCDOS.Add(h, -1.0)
for name, h in hMCSSs:
    hQCDSS.Add(h, -1.0)

c = TCanvas("c", "c", 800, 800)
c.Divide(2,2)
c.cd(1)
hRDOS.SetMinimum(0)
hRDOS.Draw()
hMCOS.Draw("same")
hRDOS.Draw("sameaxis")

c.cd(2)
hRDSS.SetMinimum(0)
hRDSS.Draw()
hMCSS.Draw("same")
hRDSS.Draw("sameaxis")

c.cd(3)
hQCDOS.SetMinimum(0)
hQCDOS.Draw()

c.cd(4)
hQCDSS.SetMinimum(0)
hQCDSS.Draw()
