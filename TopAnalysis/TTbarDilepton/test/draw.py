#!/usr/bin/env python

import sys, os

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")

from TopAnalysis.TTbarDilepton.PlotTool import *

#plotTool = PlotTool("Run2012", "Summer12", "hEventCounter")
plotTool = PlotTool("Run2012", "Summer12", "EventSummary")
plotTool.draw("image")
plotTool.printCutFlow()
