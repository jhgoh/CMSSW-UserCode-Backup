#!/usr/bin/env python

import sys, os
from ROOT import *

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

class DatasetXML:
    def __init__(self, fileName, realDataName, mcProdName):
        from xml.dom.minidom import parse
        dsxml = parse(fileName)
        realDataset = dsxml.getElementById(realDataName)
        simDataset  = dsxml.getElementById(mcProdName  )

        if not realDataset:
            print "Real dataset not found"
            return
        if realDataset.getAttribute("type") != 'Data':
            print "Wrong real data input"
            return

        if not simDataset:
            print "Dataset not found"
            return
        if simDataset.getAttribute("type") != 'MC':
            print "Wrong MC sample input"
            return

        self.realDataName = realDataName
        self.mcProdName = mcProdName 

        self.lumi = {}
        self.dsNames = {}
        self.dataFiles = {}

        self.channels = []

        for dataset in realDataset.getElementsByTagName("dataset"):
            dsName = dataset.getAttribute("name")
            lumi = float(dataset.getAttribute("lumi"))
            mode = dataset.getAttribute("mode")

            if mode not in self.lumi:
                self.lumi[mode] = 0.
                self.dsNames[mode] = []
                self.dataFiles[mode] = []
            self.lumi[mode] += lumi
            self.dsNames[mode].append(dsName)
            
        for channel in simDataset.getElementsByTagName("channel"):
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
                ch.samples.append(sample)

            self.channels.append(ch)

class PlotTool:
    def __init__(self, realDataName, mcProdName, hNEventName = "Step 0/hNEvent"):
        self.realDataName = realDataName

        self.modes = ["mm", "ee", "me", "ll"]

        self.histNames = set()
        self.steps = set()

        self.cachedObjs = []
        self.cachedCanvases = {}

        xmlPath = "%s/src/TopAnalysis/TTbarDilepton/data/dataset.xml" % os.environ["CMSSW_BASE"]
        self.dset = DatasetXML(xmlPath, realDataName, mcProdName)
        from xml.dom.minidom import parse
        psxml = parse("%s/src/TopAnalysis/TTbarDilepton/data/plots.xml" % os.environ["CMSSW_BASE"])
        plotSteps = [s.getAttribute("name") for s in psxml.getElementsByTagName("step")]

        for mode in self.dset.dataFiles:
            for dsName in self.dset.dsNames[mode]:
                self.dset.dataFiles[mode].append(TFile("hist/%s_%s.root" % (dsName, mode)))

        for channel in self.dset.channels:
            for sample in channel.samples:
                sampleName = sample.name
                sample.file_mm = TFile("hist/%s-%s_mm.root" % (mcProdName, sampleName))
                sample.file_ee = TFile("hist/%s-%s_ee.root" % (mcProdName, sampleName))
                sample.file_me = TFile("hist/%s-%s_me.root" % (mcProdName, sampleName))
                sample.nEvent = sample.file_mm.Get(hNEventName).GetBinContent(1) ## User number of event in the first mode

                ## Make list of cut steps and histograms
                stepKeys = sample.file_mm.GetListOfKeys()
                for stepKey in stepKeys:
                    subDirName = stepKey.GetName()
                    if subDirName not in plotSteps: continue
                    subDir = sample.file_mm.Get(subDirName)
                    if not subDir.IsA().InheritsFrom("TDirectory"): continue
                    self.steps.add(subDirName)

                    histKeys = subDir.GetListOfKeys()
                    for histKey in histKeys:
                        histName = histKey.GetName()
                        hist = subDir.Get(histName)
                        if not hist.IsA().InheritsFrom("TH1"): continue

                        self.histNames.add(histName)

    def makeDataMergedHistogram(self, histName, step, options):
        hists = {}
        for mode in self.modes[:-1]:
            hist = None

            for f in self.dset.dataFiles[mode]:
                hSrc = f.Get("%s/%s" % (step, histName))
                if not hSrc: continue
                hSrc.Sumw2()

                if not hist:
                    hist = hSrc.Clone("hdata_%s_%s_%s" % (histName, mode, step))
                    hist.Reset()

                    self.cachedObjs.append(hist)

                hist.Add(hSrc) 

            if not hist: continue
            nbin = hist.GetNbinsX()
            hist.AddBinContent(nbin, hist.GetBinContent(nbin+1))
            hists[mode] = hist

        hist = hists[self.modes[0]].Clone("hdata_%s_%s_%s" % (histName, self.modes[-1], step))
        hist.Reset()
        self.cachedObjs.append(hist)
        for mode in hists:
            hist.Add(hists[mode])
        hists[self.modes[-1]] = hist

        for mode in hists:
            hist = hists[mode]
            if "ymin" in options:
                hist.SetMinimum(float(options["ymin"]))
            if "ymax" in options:
                hist.SetMaximum(float(options["ymax"]))

        return hists

    def makeChannelMergedHistogram(self, histName, channelName, step, options):
        channel = None
        for ch in self.dset.channels:
            if ch.name == channelName:
                channel = ch
                break

        hists = {}
        for mode in self.modes[:-1]:
            hist = None

            for sample in channel.samples:
                hSrc = getattr(sample, "file_%s" % mode).Get("%s/%s" % (step, histName))
                if not hSrc: continue

                if not hist:
                    hist = hSrc.Clone("hch_%s_%s_%s_%s" % (histName, channelName, mode, step))
                    hist.SetFillColor(channel.color)
                    hist.Reset()
                    self.cachedObjs.append(hist)

                hist.Add(hist, hSrc, 1, sample.xsec/sample.nEvent*self.dset.lumi[mode])

            if not hist: continue
            nbin = hist.GetNbinsX()
            hist.AddBinContent(nbin, hist.GetBinContent(nbin+1))
            hists[mode] = hist

        hist = hists[self.modes[0]].Clone("hch_%s_%s_%s_%s" % (histName, channelName, self.modes[-1], step))
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
            legend.SetLineColor(kWhite)
            plotElements[mode] = (hStack, legend)

        hList = {}
        for channel in self.dset.channels:
            hists = self.makeChannelMergedHistogram(histName, channel.name, step, options)
            for mode in self.modes:
                hStack, legend = plotElements[mode]
                if "ymin" in options:
                    hStack.SetMinimum(float(options["ymin"]))
                if "ymax" in options:
                    hStack.SetMaximum(float(options["ymax"]))

                if mode not in hists: continue
                hStack.Add(hists[mode])

                if mode not in hList:
                    hList[mode] = []
                hList[mode].append((channel.name, hists[mode]))

        for mode in hList:
            hList[mode].reverse()

        for mode in hList:
            for label, h in hList[mode]:
                hStack, legend = plotElements[mode]
                legend.AddEntry(h, label, "f")

        return plotElements

    def draw(self, outDir):
        ## Build plot list from plot.xml
        plotList = {}

        from xml.dom.minidom import parse
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
                dataHists = self.makeDataMergedHistogram(name, step, opt)
                plotElement = self.makeStackWithLegend(name, step, opt)
                plotName = plotElement[self.modes[-1]][0].GetName()
                plotName = plotName.replace(" ", "_")
                if plotName in self.cachedCanvases: plotName += "_"
                c = TCanvas("c_%s" % plotName, "%s" % plotName, 800, 800)
                c.Divide(2,2)
                self.cachedCanvases[plotName] = c
                for i, mode in enumerate(self.modes):
                    pad = c.cd(i+1)

                    if opt["logy"] : pad.SetLogy()

                    if mode not in dataHists: continue
                    hData = dataHists[mode]
                    hist, legend = plotElement[mode]
                    legend.AddEntry(hData, self.realDataName, "lp")

                    hData.Draw()
                    legend.Draw()
                    hist.Draw("same")
                    hData.Draw("same")
                    hData.Draw("sameaxis")

                    self.cachedObjs.extend([hData, hist, legend])

                c.Print("image/"+c.GetName()+".pdf")
                c.Print("image/"+c.GetName()+".png")

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
            for channel in self.dset.channels:
                for sample in channel.samples:
                    file = getattr(sample, "file_%s" % mode)

                    line = " %20s |" % sample.name
                    for step in self.steps:
                        nPassing = file.Get("%s/hNEvent" % step).GetBinContent(3)
                        line += "%10.2f" % (nPassing*sample.xsec/sample.nEvent*self.dset.lumi[mode])
                    line += "\n"

                    print line,

        print outBorder

class PlotBuilder:
    def __init__(self, realDataName, mcProdName, mode, hNEventName="hEvents", outFile=gROOT):
        self.mode = mode
        self.outFile = outFile

        xmlPath = "%s/src/TopAnalysis/TTbarDilepton/data/dataset.xml" % os.environ["CMSSW_BASE"]
        self.dset = DatasetXML(xmlPath, realDataName, mcProdName)

        self.lumi = self.dset.lumi[mode]
        for dsName in self.dset.dsNames[mode]:
            self.dset.dataFiles[mode].append(TFile("hist/%s_%s.root" % (dsName, mode)))

        for channel in self.dset.channels:
            for sample in channel.samples:
                sample.file = TFile("hist/%s-%s_%s.root" % (mcProdName, sample.name, mode))
                sample.nEvent = sample.file.Get(hNEventName).GetBinContent(1)

    def buildRDHist(self, category, histName):
        hRD = None
        for file in self.dset.dataFiles[self.mode]:
            hRD_src = file.Get("%s/%s" % (category, histName))
            if not hRD:
                self.outFile.cd()
                hRD = hRD_src.Clone("%s_%s_%s" % (self.mode, category, histName))
                hRD.Reset()
                hRD.Sumw2()
            hRD.Add(hRD_src)

        return hRD

    def buildMCHists(self, category, histName):
        hStack = THStack("hs_%s_%s_%s" % (self.mode, category, histName), "hStack %s %s %s" % (self.mode, category, histName))
        hMCs = []
        for channel in self.dset.channels:
            hMC = None
            for sample in channel.samples:
                hMC_src = sample.file.Get("%s/%s" % (category, histName))
                if not hMC:
                    self.outFile.cd()
                    hMC = hMC_src.Clone("%s_%s_%s" % (self.mode, category, histName))
                    hMC.Reset()
                    hMC.SetFillColor(channel.color)
                hMC.Add(hMC_src, self.lumi*sample.xsec/sample.nEvent)
            if not hMC: continue

            hStack.Add(hMC)
            hMCs.append((channel.name, hMC))

        return hStack, hMCs

class NtupleAnalyzerLite:
    def __init__(self, realDataName, mcProdName, mode, dsName, srcDir, hNEventName="hEvents", treeName="tree"):
        self.categories = []
        self.hists = {}
        self.inputFile = None
        self.inputTree = None
        self.outFile = None

        xmlPath = "%s/src/TopAnalysis/TTbarDilepton/data/dataset.xml" % os.environ["CMSSW_BASE"]
        self.dset = DatasetXML(xmlPath, realDataName, mcProdName)

        if dsName in self.dset.dsNames[mode]:
            name = "%s_%s" % (dsName, mode)
            self.inputFile = TFile("%s/%s.root" % (srcDir, dsName))
            #self.inputTree = self.inputFile.Get("%s/%s" % (mode, treeName))
            self.inputTree = self.inputFile.Get("%s" % treeName)
            self.outFile = TFile("hist/%s.root" % name, "RECREATE")
            #hNEvent = self.inputFile.Get("%s/%s" % (mode, hNEventName)).Clone(hNEventName.replace("/", "_"))
            hNEvent = self.inputFile.Get("%s" % hNEventName).Clone(os.path.basename(hNEventName))
            self.nEvent = hNEvent.GetBinContent(1)
            hNEvent.Write()
        else:
            samples = {}
            for channel in self.dset.channels:
                for sample in channel.samples:
                    samples[sample.name] = sample

            if dsName in samples.keys():
                sample = samples[dsName]
                name = "%s-%s_%s" % (mcProdName, sample.name, mode)
                self.inputFile = TFile("%s/%s-%s.root" % (srcDir, mcProdName, dsName))
                #self.inputTree = self.inputFile.Get("%s/%s" % (mode, treeName))
                self.inputTree = self.inputFile.Get("%s" % treeName)
                self.outFile = TFile("hist/%s.root" % name, "RECREATE")
                #hNEvent = self.inputFile.Get("%s/%s" % (mode, hNEventName)).Clone(hNEventName.replace("/", "_"))
                hNEvent = self.inputFile.Get("%s" % hNEventName).Clone(os.path.basename(hNEventName))
                self.nEvent = hNEvent.GetBinContent(1)
                hNEvent.Write()
            else:
                print "Cannot find sample", dsName
                

    def addCategory(self, name, cut):
        self.categories.append((name, cut))

    def addHistogram(self, varexp, name, title, nbin, xmin, xmax):
        self.hists[name] = (varexp, title, nbin, xmin, xmax)

    def scanCutFlow(self):
        for categoryName, cut in self.categories:
            cutStepDir = self.outFile.mkdir(categoryName)
            cutStepDir.cd()

            h = TH1F("hNEvent", "Number of events", 3, 1, 4)
            h.Fill(1, self.nEvent)
            self.inputTree.Draw("2>>+hNEvent", "(%s)" % cut, "goff")
            self.inputTree.Draw("3>>+hNEvent", "(puweight)*(%s)" % cut, "goff")
            h.Write()

            for histName in self.hists:
                varexp, title, nbin, xmin, xmax = self.hists[histName]

                h = TH1F(histName, title, nbin, xmin, xmax)

                self.inputTree.Draw("%s>>%s" % (varexp, histName), "(puweight)*(%s)" % cut, "goff")

                h.Write()
