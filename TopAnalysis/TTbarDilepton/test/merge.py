#!/usr/bin/env python

import sys, os

datasetDir = "%s/src/TopAnalysis/TTbarDilepton/data" % os.environ["CMSSW_BASE"]
logDir = "log"
ntupleDir = "ntuple/unmerged"
outputDir = "ntuple/merged"

## Get list of files
datasets = {}
for fileName in os.listdir(datasetDir):
    if len(fileName) < 4: continue
    ext = fileName.split('.')[-1]
    if ext != 'txt': continue

    dataset = os.path.basename(fileName).replace('dataset-', '')
    dataset = dataset.replace('.txt', '')

    datasets[dataset] = []
    for line in open(datasetDir+"/"+fileName).readlines():
        line = line.strip()
        if '.root' not in line: continue
        if '#' == line[0]: continue

        line = line.replace('root://eoscms', '')
        line = line.replace('/eos/cms', '')
        line = line.replace('//', '/')

        datasets[dataset].append(line)

## Analyze log files to get file processing history
logs = {}
for fileName in os.listdir(logDir):
    if len(fileName) < 5: continue
    ext = fileName.split('.')[-1]

    if ext not in ('log', 'txt'): continue

    sample = '.'.join(fileName.split('.')[:-1])
    #index = sample.split('_')[-1].lstrip('0')
    index = int(sample.split('_')[-1])
    sample = '%s_%d' % ('_'.join(sample.split('_')[:-1]), index)

    if sample not in logs:
        logs[sample] = []

    logs[sample].append(fileName)

## List up ntuples
ntuples = {}
for fileName in os.listdir(ntupleDir):
    if len(fileName) < 6: continue
    ext = fileName.split('.')[-1]

    if ext != 'root': continue

    sample = '.'.join(fileName.split('.')[:-1])
    #index = int(sample.split('_')[-1].lstrip('0'))
    index = int(sample.split('_')[-1])
    sample = '%s_%d' % ('_'.join(sample.split('_')[:-1]), index)

    if sample not in ntuples:
        ntuples[sample] = []

    ntuples[sample].append(fileName)

errors = {}
successFiles = []

samples = list(set(logs.keys()) | set(ntuples.keys()))
for sample in samples:
    errors[sample] = []

    if sample not in logs.keys():
        print "!!! Sample %s is not in log file list" % sample
        continue

    if sample not in ntuples.keys():
        print "!!! Sample %s is not in ntuple list" % sample
        continue

    openedFiles = {}
    for log in logs[sample]:
        for line in open(logDir+"/"+log).readlines():
            if '.root' not in line: continue
            line = line.strip()

            fileName = line.split()[-1].replace('root://eoscms', '')
            fileName = fileName.replace('/eos/cms', '')
            fileName = fileName.replace('//', '/')
            state = ' '.join(line.split()[-3:-2])

            if fileName not in openedFiles:
                openedFiles[fileName] = []

            openedFiles[fileName].append(state)

    for fileName in openedFiles:
        state = openedFiles[fileName]

        if 'open' not in state:
            errors[sample].append(fileName)
            print "Init error", fileName
            continue
        if 'opened' not in state:
            print "Open error", fileName
            errors[sample].append(fileName)
            continue
        if 'Closed' not in state:
            print "Close error", fileName
            errors[sample].append(fileName)
            continue

        successFiles.append(fileName)

isValid = True
for dataset in datasets:
    for file in datasets[dataset]:
        if file not in successFiles:
            print "!!! File %s missing in log fies" % file
            isValid = False

for sample in samples:
    if len(errors[sample]) != 0: isValid = False
    for file in errors[sample]:
        print sample
        print file

if not isValid:
    print "File processing error"
    sys.exit(1)
else:
    print "All files are successfully processed"

## Start merging
for dataset in datasets:
    print "Merging %s" % dataset
    os.system("hadd %s/%s.root %s/%s_*.root" % (outputDir, dataset, ntupleDir, dataset))

