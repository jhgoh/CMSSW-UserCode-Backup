#!/usr/bin/python

import os, tempfile
import cgi

def drawDefault():
	print "Content-type: image/gif\n\n",
	print open("cmsLogo.gif").read()

# Process CGI variables
form = cgi.FieldStorage()

try:
	rootFileName = form["ROOTFile"].value
	objPath = form["ObjectPath"].value
except:
	drawDefault();
	exit(0);

# ROOT stuffs
from ROOT import gROOT, gStyle
from ROOT import TCanvas, TFile, TH1, TH2
from ROOT import TF1

gROOT.Reset()
gStyle.SetPalette(1)

canvas = TCanvas("canvas", "ROOT canvas", 200, 10, 700, 500)

try:
	rootFile = TFile(rootFileName)
	drawObject = rootFile.Get(objPath)
	drawObject.Draw()
except:
	drawDefault()
	exit(0)

imgFileName = tempfile.mktemp('.gif')

# Start to make it visible
print "Content-type: image/gif\n\n",

canvas.Print(imgFileName)

img = open(imgFileName, "rb")
print img.read()

os.remove(imgFileName)
