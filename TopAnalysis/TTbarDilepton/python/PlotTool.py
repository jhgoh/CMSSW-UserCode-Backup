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

class PlotsetXML:
    def __init__(self, fileName):
        self.steps = []
        self.plots = {}

        from xml.dom.minidom import parse
        plxml = parse(fileName)
        for step in plxml.getElementsByTagName("step"):
            stepName = step.getAttribute("name")
            ymin = float(step.getAttribute("ymin"))
            ymax = float(step.getAttribute("ymax"))
            logy = step.getAttribute("logy").lower() == "true"
            
            self.steps.append(stepName)
            self.plots[stepName] = []
            for plot in step.getElementsByTagName("plot"):
                plotName = plot.getAttribute("name")
                plotOpt = {"ymin":ymin, "ymax":ymax, "logy":logy}

                if plot.hasAttribute("ymin"): plotOpt["ymin"] = float(plot.getAttribute("ymin"))
                if plot.hasAttribute("ymax"): plotOpt["ymax"] = float(plot.getAttribute("ymax"))
                if plot.hasAttribute("logy"): plotOpt["logy"] = plot.getAttribute("logy").lower() == "true"
                if plot.hasAttribute("mode"): plotOpt["mode"] = plot.getAttribute("mode").split(",")

                self.plots[stepName].append( (plotName, plotOpt) )

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
    def __init__(self, realDataName, mcProdName, hNEventName = "Step 0/hNEvent", histDir="hist"):
        self.realDataName = realDataName
        self.cachedCanvases = {}
        self.cachedObjects = []

        self.modes = ["mm", "ee", "me", "ll"]

        self.steps = []
        self.plots = {}
        ## Build plotting interfaces and check consistency
        self.plotters = {}
        for mode in self.modes[:-1]:
            self.plotters[mode] = PlotBuilder(realDataName, mcProdName, mode, hNEventName, histDir=histDir)
            steps = self.plotters[mode].plotset.steps
            plots = self.plotters[mode].plotset.plots
            if len(self.steps) == 0:
                for step in steps:
                    if step not in self.plotters[mode].steps: continue
                    self.steps.append(step)

                    self.plots[step] = []
                    for plotName, plotOpt in plots[step]:
                        if plotName not in self.plotters[mode].plots: continue
                        self.plots[step].append((plotName, plotOpt))

    def makeMergedPlot(self, category, histName, hRD_srcs, hMCs_srcs):
        hRD = None
        labels = []
        hMCs = []

        for mode in hRD_srcs:
            hRD_src = hRD_srcs[mode]
            if not hRD:
                hRD = hRD_src.Clone("hrd_ll_%s_%s" % (category, histName))
                hRD.Reset()
            hRD.Add(hRD_src)

            for i, (label, hMC_src) in enumerate(hMCs_srcs[mode]):
                if i >= len(hMCs):
                    labels.append(label)
                    hMC = hMC_src.Clone("hch%d_ll_%s_%s" % (i, category, histName))
                    hMCs.append((label, hMC))
                else:
                    label, hMC = hMCs[labels.index(label)]
                    hMC.Add(hMC_src)

        return (hRD, hMCs)

    def buildStack(self, mode, category, histName, hMCs):
        hsName = "hs_%s_%s_%s" % (mode, category, histName)
        hstack = THStack(hsName, hsName)
        for label, hMC in hMCs:
            hstack.Add(hMC)

        return hstack
        
    def buildLegend(self, labelAndHists, reverse=True):
        lh = labelAndHists[:]
        if reverse: lh.reverse()
        legend = TLegend(0.70, 0.70, 0.93, 0.93, "", "NDC")
        legend.SetFillStyle(0)
        #legend.SetLineStyle(0)
        legend.SetLineColor(kWhite)
        for label, hist in lh:
            legend.AddEntry(hist, label, "f")
        return legend

    def draw(self, outDir):
        for step in self.steps:
            for plotName, plotOpt in self.plots[step]:
                hRD = {}
                hsMC = {}
                hMCs = {}
                legend = {}
                for mode in self.modes[:-1]:
                    plotter = self.plotters[mode]
                    hRD[mode] = plotter.buildRDHist(step, plotName)
                    hMCs[mode] = plotter.buildMCHists(step, plotName)
                    hsMC[mode] = self.buildStack(mode, step, plotName, hMCs[mode])

                    self.cachedObjects.extend([hRD[mode], hsMC[mode], hMCs[mode]])
                hRD["ll"], hMCs["ll"] = self.makeMergedPlot(step, plotName, hRD, hMCs)
                hsMC["ll"] = self.buildStack("ll", step, plotName, hMCs["ll"])
                self.cachedObjects.extend([hRD["ll"], hsMC["ll"], hMCs["ll"]])

                for mode in self.modes:
                    legend[mode] = self.buildLegend(hMCs[mode])
                    legend[mode].AddEntry(hRD[mode], self.realDataName, "lp")

                    self.cachedObjects.extend([legend[mode]])

                canvasName = "c_%s_%s" % (step.replace(" ", "_"), plotName.replace(" ", "_"))
                while canvasName in self.cachedCanvases:
                    canvasName += "_"
                c = TCanvas(canvasName, canvasName, 600, 600)
                c.Divide(2,2)
                self.cachedCanvases[canvasName] = c

                for i, mode in enumerate(self.modes):
                    pad = c.cd(i+1)

                    if plotOpt["logy"] : pad.SetLogy()

                    if mode not in hRD: continue
                    hRD[mode].SetMinimum(plotOpt["ymin"])
                    hRD[mode].SetMaximum(plotOpt["ymax"])
                    hRD[mode].Draw()
                    hsMC[mode].Draw("same")
                    hRD[mode].Draw("same")
                    hRD[mode].Draw("sameaxis")
                    legend[mode].Draw()

                if outDir is not None:
                    c.Print(os.path.join(outDir, c.GetName()+".pdf"))
                    c.Print(os.path.join(outDir, c.GetName()+".png"))

    def printCutFlow(self):
        outBorder = "="*(22+10*len(self.steps)+2)
        separator = " "+"-"*(22+10*len(self.steps))+" "

        for mode in self.modes[:-1]:
            print outBorder
            header = " %20s |" % mode
            header += "".join(["%10s" % step for step in self.steps])

            print header
            print separator

            plotter = self.plotters[mode]

            sumMCs = [0.]*len(self.steps)
            ## Get cut steps and build contents
            for channel in plotter.dataset.channels:
                for sample in channel.samples:
                    line = " %20s |" % sample.name
                    for step in self.steps:
                        nPass = sample.file.Get("%s/hNEvent" % step).GetBinContent(4)
                        #nPassNorm = nPass*sample.xsec/sample.nEvent*plotter.dataset.lumi[mode]
                        nPassNorm = nPass*sample.xsec*plotter.dataset.lumi[mode]
                        sumMCs[self.steps.index(step)] += nPassNorm
                        line += "%10.2f" % nPassNorm
                    line += "\n"

                    print line,

            print separator
            line = " %20s |" % "MC total"
            for sumMC in sumMCs: line += "%10.2f" % sumMC
            print line
            print separator
            line = " %20s |" % self.realDataName
            dataFiles = plotter.dataset.dataFiles[mode]
            for step in self.steps:
                nPass = sum([f.Get("%s/hNEvent" % step).GetBinContent(4) for f in dataFiles])
                line += "%10.2f" % nPass
            print line

        print outBorder

class PlotBuilder:
    def __init__(self, realDataName, mcProdName, mode, hNEventName="hNEvent", histDir="hist", outFile=gROOT):
        self.mode = mode
        self.outFile = outFile

        datasetXMLPath = "%s/src/TopAnalysis/TTbarDilepton/data/dataset.xml" % os.environ["CMSSW_BASE"]
        self.dataset = DatasetXML(datasetXMLPath, realDataName, mcProdName)

        plotsetXMLPath = "%s/src/TopAnalysis/TTbarDilepton/data/plots.xml" % os.environ["CMSSW_BASE"]
        self.plotset = PlotsetXML(plotsetXMLPath)

        self.lumi = self.dataset.lumi[mode]
        for dsName in self.dataset.dsNames[mode]:
            self.dataset.dataFiles[mode].append(TFile("%s/%s-%s.root" % (histDir, dsName, mode)))

        for channel in self.dataset.channels:
            for sample in channel.samples:
                sample.file = TFile("%s/%s-%s-%s.root" % (histDir, mcProdName, sample.name, mode))
                sample.nEvent = sample.file.Get(hNEventName).GetBinContent(1)

        ## Get list of steps in the file using data file and xml interface
        self.steps = []
        self.plots = []
        dataFile = self.dataset.dataFiles[mode][0]
        for stepName in [d.GetName() for d in dataFile.GetListOfKeys()]:
            stepDir = dataFile.Get(stepName)
            if not stepDir.IsA().InheritsFrom("TDirectory"): continue

            self.steps.append(stepName)

            for plotName in [k.GetName() for k in stepDir.GetListOfKeys()]:
                plot = stepDir.Get(plotName)
                if not plot.IsA().InheritsFrom("TH1"): continue

                self.plots.append(plotName)

    def buildRDHist(self, category, histName):
        hRD = None
        for file in self.dataset.dataFiles[self.mode]:
            hRD_src = file.Get("%s/%s" % (category, histName))
            if not hRD_src:
                print "Cannot find", category, histName
            if not hRD:
                self.outFile.cd()
                hRD = hRD_src.Clone("hrd_%s_%s_%s" % (self.mode, category, histName))
                hRD.Reset()
                hRD.Sumw2()
            hRD.Add(hRD_src)

        return hRD

    def buildMCHists(self, category, histName):
        hMCs = []
        for i, channel in enumerate(self.dataset.channels):
            hMC = None
            for sample in channel.samples:
                hMC_src = sample.file.Get("%s/%s" % (category, histName))
                if not hMC:
                    self.outFile.cd()
                    hMC = hMC_src.Clone("hch%d_%s_%s_%s" % (i, self.mode, category, histName))
                    hMC.Reset()
                    hMC.SetFillColor(channel.color)
                hMC.Add(hMC_src, self.lumi*sample.xsec)
            if not hMC: continue

            hMCs.append((channel.name, hMC))

        return hMCs

class NtupleAnalyzerLite:
    def __init__(self, inputTreePath, outputFilePath, hNEventPath, normFactor, weightVar="puweight"):
        self.weightVar = weightVar
        self.hists = {}
        self.categories = []

        inputFileName, inputTreeName = inputTreePath.split(":")
        self.inputFile = TFile(inputFileName)
        self.inputTree = self.inputFile.Get(inputTreeName)

        self.outputFile = TFile(outputFilePath, "RECREATE")
        self.hNEvent = self.inputFile.Get(hNEventPath).Clone("hNEvent")
        self.hNEvent.Write()

        ## Get overall normalization factor
        if type(normFactor) in (type(0.0), type(0)):
            self.nEvent = 0
            self.scale = float(normFactor)
        elif type(normFactor) == type("") and ':' in normFactor:
            ## Try to get event counter histogram
            ## The form should be "Some/directory/histogram:binNumber"
            normHistPath, binNumber = normFactor.split(":")
            self.nEvent = self.inputFile.Get(normHistPath).GetBinContent(int(binNumber))
            self.scale = 1.0/self.nEvent
        else:
            print "@@ Invalid normFactor,", normFactor, inputTreePath
            return

    def addCategory(self, name, cut, histNames, options):
        if type(options) == type(''): options = [options]
        self.categories.append( (name, cut, histNames, options) )

    def addHistogram(self, varexp, name, title, nbin, xmin, xmax):
        self.hists[name] = (varexp, title, nbin, xmin, xmax)

    def scanCutFlow(self):
        for categoryName, cut, histNames, options in self.categories:
            categoryDir = self.outputFile.mkdir(categoryName)
            categoryDir.cd()

            tree = self.inputTree

            if "renewCut" in options: tree.SetEventList(0)

            eventList = TEventList("eventList", "")
            tree.Draw(">>eventList", cut)
            tree.SetEventList(eventList)
            cut = "1"

            h = TH1F("hNEvent", "Number of events", 4, 1, 5)
            h.Fill(1, self.hNEvent.GetBinContent(1))
            tree.Draw("2>>+hNEvent", "%s" % cut, "goff")
            tree.Draw("3>>+hNEvent", "(%s)*(%s)" % (self.weightVar, cut), "goff")
            tree.Draw("4>>+hNEvent", "(%s*%g)*(%s)" % (self.weightVar, self.scale, cut), "goff")

            h.GetXaxis().SetBinLabel(1, "Generated event")
            h.GetXaxis().SetBinLabel(2, "Event after cut")
            h.GetXaxis().SetBinLabel(3, "Event after cut with weight")
            h.GetXaxis().SetBinLabel(4, "Normalized event after cut with weight")
            h.Write()
        
            for histName in histNames:
                if histName not in self.hists: continue

                varexp = self.hists[histName][0]
                options = self.hists[histName][1:]
                h = TH1F(histName, *options)

                tree.Draw("min(%s,%f)>>%s" % (varexp, h.GetXaxis().GetXmax()-1e-9, histName), "(%s)*(%s)" % (self.weightVar, cut), "goff")
                h.Scale(self.scale)

                h.Write()

