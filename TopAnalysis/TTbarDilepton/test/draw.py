#!/usr/bin/env python

gMCProd = "Summer12"
gRealData = ["Run2012A", "Run2012B", "Run2012C"]
gLumi = 30.0*1e3

import sys, os
from xml.dom.minidom import parse
from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")

class Channel:
    def __init__(self):
        self.name = ""
        self.color = kBlack
        self.samples = []

class Sample:
    def __init__(self):
        self.name = ""
        self.xsec = 0.0
        self.file = None

class PlotTool:
    def __init__(self):
        self.modes = ["mm", "ee", "me", "ll"]

        self.channels = []
        channels = self.channels

        self.histNames = set()
        self.steps = set()

        self.cachedObjs = []
        self.cachedCanvases = {}

        dsxml = parse("%s/src/TopAnalysis/TTbarDilepton/data/dataset.xml" % os.environ["CMSSW_BASE"])
        dataset = dsxml.getElementById(gMCProd)

        if not dataset:
            print "Dataset not found"
            return

        for channel in dataset.getElementsByTagName("channel"):
            channelName = channel.getAttribute("name")
            color = channel.getAttribute("color")
            if not color or len(color) == 0:
                print "Color is not set for channel, ", channelName
                print "Skipping this channel"
                continue
            color = gROOT.ProcessLine(color+";")

            samples = channel.getElementsByTagName("sample")
            if len(samples) == 0:
                print "Warning : no sample assigned to channel. Ignore this channel"
                continue

            ch = Channel()
            ch.name = channelName
            ch.color = color
            for sample in samples:
                sampleName = sample.getAttribute("name")
                xsec = sample.getAttribute("xsec")
                if not xsec or len(xsec) == 0:
                    print "Cross-section is not set for sample, ", sampleName
                    print "Skipping this sample"
                    continue
                xsec = 1.0*gROOT.ProcessLine(xsec+";")

                sample = Sample()
                sample.name = sampleName
                sample.xsec = xsec
                sample.file_mm = TFile("hist/%s-%s_mm.root" % (gMCProd, sampleName))
                sample.file_ee = TFile("hist/%s-%s_ee.root" % (gMCProd, sampleName))
                sample.file_me = TFile("hist/%s-%s_me.root" % (gMCProd, sampleName))
                sample.nEvent = sample.file_mm.Get("Step 0/hNEvent").GetBinContent(1) ## User number of event in the first mode

                ## Make list of cut steps and histograms
                stepKeys = sample.file_mm.GetListOfKeys()
                for stepKey in stepKeys:
                    subDirName = stepKey.GetName()
                    subDir = sample.file_mm.Get(subDirName)
                    if not subDir.IsA().InheritsFrom("TDirectory"): continue
                    self.steps.add(subDirName)

                    histKeys = subDir.GetListOfKeys()
                    for histKey in histKeys:
                        histName = histKey.GetName()
                        hist = subDir.Get(histName)
                        if not hist.IsA().InheritsFrom("TH1"): continue

                        self.histNames.add(histName)

                ch.samples.append(sample)

            channels.append(ch)

    def makeChannelMergedHistogram(self, histName, channelName, step, options):
        channel = None
        for ch in self.channels:
            if ch.name == channelName:
                channel = ch
                break

        hists = {}
        for mode in self.modes[:-1]:
            hist = None

            for sample in channel.samples:
                hSrc = getattr(sample, "file_%s" % mode).Get("%s/%s" % (step, histName))

                if not hist:
                    hist = hSrc.Clone("hch_%s_%s_%s_%s" % (histName, channelName, mode, step))
                    hist.SetFillColor(channel.color)
                    hist.Reset()
                    self.cachedObjs.append(hist)

                hist.Add(hist, hSrc, 1, sample.xsec/sample.nEvent*gLumi)

            nbin = hist.GetNbinsX()
            hist.AddBinContent(nbin, hist.GetBinContent(nbin+1))
            hists[mode] = hist

        hist = hists["mm"].Clone("hch_%s_%s_%s_%s" % (histName, channelName, self.modes[-1], step))
        hist.Reset()
        self.cachedObjs.append(hist)
        for mode in hists:
            hist.Add(hist, hists[mode], 1, 1)
        hists[self.modes[-1]] = hist

        return hists

    def makeStackWithLegend(self, histName, step, options):
        plotElements = {}
        for mode in self.modes:
            hStack = THStack("hs_%s_%s_%s" % (histName, mode, step), "%s %s %s" % (histName, mode, step))
            legend = TLegend(0.70, 0.70, 0.93, 0.93)
            legend.SetFillColor(kWhite)
            plotElements[mode] = (hStack, legend)

        hList = {}
        for channel in self.channels:
            hists = self.makeChannelMergedHistogram(histName, channel.name, step, options)
            for mode in self.modes:
                hStack, legend = plotElements[mode]
                if "ymin" in options:
                    hStack.SetMinimum(float(options["ymin"]))
                if "ymax" in options:
                    hStack.SetMaximum(float(options["ymax"]))

                hStack.Add(hists[mode])

                if mode not in hList:
                    hList[mode] = []
                hList[mode].append((channel.name, hists[mode]))

        #for mode in hList:
        #    hList[mode].reverse()

        for mode in hList:
            for label, h in hList[mode]:
                hStack, legend = plotElements[mode]
                legend.AddEntry(h, label, "f")

        return plotElements

    def draw(self):
        ## Build plot list from plot.xml
        plotList = {}
        psxml = parse("%s/src/TopAnalysis/TTbarDilepton/data/plots.xml" % os.environ["CMSSW_BASE"])
        plotCutSteps = psxml.getElementsByTagName("step")
        for plotCutStep in plotCutSteps:
            stepName = plotCutStep.getAttribute("name")

            ## Step name should be in the cut step list
            if stepName not in self.steps: continue
            plotList[stepName] = []

            ymin = float(plotCutStep.getAttribute("ymin"))
            ymax = float(plotCutStep.getAttribute("ymax"))
            logy = plotCutStep.getAttribute("logy").lower() == "true"

            for plot in plotCutStep.getElementsByTagName("plot"):
                plotName = plot.getAttribute("name")
                plotNames = []

                if '*' != plotName[-1]:
                    if plotName not in self.histNames: continue
                    plotNames.append(plotName)
                else:
                    plotName = plotName[:-1]
                    for histName in self.histNames:
                        if histName[:len(plotName)] == plotName: plotNames.append(histName)

                plotOpt = {"ymin":ymin, "ymax":ymax, "logy":logy}

                if plot.hasAttribute("ymin"): plotOpt["ymin"] = float(plot.getAttribute("ymin"))
                if plot.hasAttribute("ymax"): plotOpt["ymax"] = float(plot.getAttribute("ymax"))
                if plot.hasAttribute("logy"): plotOpt["logy"] = plot.getAttribute("logy").lower() == "true"
                if plot.hasAttribute("mode"): plotOpt["mode"] = plot.getAttribute("mode").split(",")

                for plotName in plotNames:
                    plotList[stepName].append( (plotName, plotOpt) )

        for step in plotList:
            for name, opt in plotList[step]:
                plotElement = self.makeStackWithLegend(name, step, opt)
                plotName = plotElement[self.modes[-1]][0].GetName()
                plotName.replace(" ", "_")
                if plotName in self.cachedCanvases: plotName += "_"
                c = TCanvas("c_%s" % plotName, "%s" % plotName, 800, 800)
                c.Divide(2,2)
                self.cachedCanvases[plotName] = c
                for i, mode in enumerate(self.modes):
                    pad = c.cd(i+1)

                    if opt["logy"] : pad.SetLogy()
                    hist, legend = plotElement[mode]

                    hist.Draw()
                    legend.Draw()

                    self.cachedObjs.extend([hist, legend])

    def printCutFlow(self):
        outBorder = "="*(22+10*len(self.steps)+2)
        seperator = " "+"-"*(22+10*len(self.steps))+" "

        for mode in self.modes[:-1]:
            print outBorder
            header = " %20s |" % mode
            header += "".join(["%10s" % step for step in self.steps])

            print header
            print seperator

            ## Get cut steps and build contents
            for channel in self.channels:
                for sample in channel.samples:
                    file = getattr(sample, "file_%s" % mode)

                    line = " %20s |" % sample.name
                    for step in self.steps:
                        nPassing = file.Get("%s/hNEvent" % step).GetBinContent(3)
                        line += "%10.2f" % nPassing
                    line += "\n"

                    print line,

        print outBorder

plotTool = PlotTool()
plotTool.draw()
plotTool.printCutFlow()
