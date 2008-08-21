#!/usr/bin/python

#import os
from os import remove
from tempfile import mktemp
from sys import path
from mod_python import util
from mod_python import apache

dataDir = '/users/jhgoh/public_html/CMS/Validation/data';
path.append('/usr/local/root/root_v5.19.04/lib/root')
# ROOT stuffs
from ROOT import gROOT, gStyle
from ROOT import TCanvas, TFile, TH1F, TH2F

def setValue(form, key, outdict, func):
	if key in form.keys() and key in outdict.keys():
		outdict[key] = func(form[key].value)

def index(req):
	# Process variables
	firstRun = True;
	objHolder = [];
	lineColor = 2;

	defaults = {'normalize':False, 
				'height':500, 'width':600,
				'xmin':None, 'xmax':None, 'ymin':None, 'ymax':None,
				'opt':''}

	if 'h' not in req.form.keys():
		return apache.HTTP_NOT_FOUND

	if 'norm' in req.form.keys() and req.form['norm'].value == 'yes':
		defaults['normalize'] = True

	setValue(req.form, 'opt', defaults, lambda x:string(x))
	setValue(req.form, 'width', defaults, lambda x:int(x))
	setValue(req.form, 'height', defaults, lambda x:int(x))
	setValue(req.form, 'xmin', defaults, lambda x:float(x))
	setValue(req.form, 'xmax', defaults, lambda x:float(x))
	setValue(req.form, 'ymin', defaults, lambda x:float(x))
	setValue(req.form, 'ymax', defaults, lambda x:float(x))

	ext = 'gif'
	formats = {'gif':'image/gif','png':'image/png','jpg':'image/jpg',
			   'ps':'application/postscript','pdf':'application/pdf'}
	if 'ext' in req.form.keys() and req.form['ext'] in formats.keys():
		ext = req.form['ext']
	req.content_type = formats[ext]
	
	gROOT.Reset()
	gStyle.SetPalette(1)

	canvas = TCanvas('canvas', 'ROOT canvas', 200, 10, defaults['width'], defaults['height'])

	if 'logx' in req.form.keys() and req.form['logx'].value == 'yes':
		canvas.SetLogx();
	if 'logy' in req.form.keys() and req.form['logy'].value == 'yes':
		canvas.SetLogy();

	for fullPath in req.form['h'].value.split(','):
		if fullPath == None or fullPath == '':
			continue;

		filePath = dataDir+'/'+fullPath.split(':')[0];
		histPath = fullPath.split(':')[1];

		try:
			rootFile = TFile(filePath);
			obj = rootFile.Get(histPath);
			opt = defaults['opt']

			if ( obj.IsA().InheritsFrom('TH2F') ):
				rootHist = TH2F(obj);
				opt += 'colz'
			elif ( obj.IsA().InheritsFrom('TH1F') ):
				rootHist = TH1F(obj);
				# SetMaximum/Minimum not working. why?
				if defaults['ymax']:rootHist.SetMaximum(defaults['ymax'])
				if defaults['ymin']:rootHist.SetMinimum(defaults['ymin'])
			else:
				return apache.HTTP_NOT_FOUND

			rootHist.SetLineColor(lineColor);
			lineColor = lineColor+1;

			if not defaults['normalize']:
				rootHist.Draw(opt);
			else:
				rootHist.DrawNormalized(opt);

			objHolder.append( (rootFile,rootHist,opt) );

		except Exception, what:
			return apache.HTTP_NOT_FOUND

		if firstRun : 
			defaults['opt'] = defaults['opt']+'same';
			firstRun = False;

	imgFileName = mktemp('.'+ext);

	# Start to make it visible
	canvas.Print(imgFileName);

	img = open(imgFileName, 'rb');
	imgContent = img.read();

	remove(imgFileName);

	return imgContent
	#return apache.OK
