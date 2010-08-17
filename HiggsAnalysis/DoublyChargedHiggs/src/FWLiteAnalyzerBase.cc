#include "HiggsAnalysis/DoublyChargedHiggs/interface/FWLiteAnalyzerBase.h"

#include "TPRegexp.h"
#include "TSystem.h"
#include "TString.h"
#include "TH1F.h"

#include <iostream>

using namespace std;

FWLiteAnalyzerBase::FWLiteAnalyzerBase(const string outFileName, const bool verbose):
  verbose_(verbose)
{
  outFile_ = TFile::Open(outFileName.c_str(), "RECREATE");
}

void FWLiteAnalyzerBase::AddSignal(const string name, const string inputFile, const double xsec)
{
  signalXSecTable_[name] = xsec;
  AddFile(inputFile, signalFiles_[name]);
}

void FWLiteAnalyzerBase::AddBackground(const string name, const string inputFile, const double xsec)
{
  backgroundXSecTable_[name] = xsec;
  AddFile(inputFile, backgroundFiles_[name]);
}

void FWLiteAnalyzerBase::AddFile(const string inputName, vector<string>& inputFiles)
{
  TPRegexp rootFilePattern("\\.root$");
  TPRegexp plainFileExtPattern("\\.[a-zA-Z0-9]\\+$");

  const bool isRfio = TPRegexp("^rfio:").Match(inputName);

  if ( rootFilePattern.Match(inputName) )
  {
    inputFiles.push_back(inputName);
  }
  else if ( plainFileExtPattern.Match(inputName) )
  {
    if ( verbose_ ) cout << "@@ Only root files are supported\n";
    return;
  }
  else // Directories
  {
    if ( verbose_ ) cout << "@@ Reading root files in the directory " << inputName << endl;

    string dirCmd;
    string baseDir;

    if ( isRfio )
    {
      dirCmd = "nsls ";
      baseDir = inputName.substr(5, inputName.size());
    }
    else
    {
      dirCmd = "ls ";
      baseDir = inputName;
    }

    // Get list of files and collect root files
    if ( verbose_ ) cout << "@@@@ Executing directory listing : " << dirCmd << baseDir << endl;
    FILE* pipe = gSystem->OpenPipe((dirCmd+baseDir).c_str(), "r");
    const int nBufferSize = 200;
    char buffer[nBufferSize];
    while ( !feof(pipe) and fgets(buffer, nBufferSize, pipe) != NULL )
    {
      const string line(buffer);
      const string foundFileName = line.substr(0, line.size()-1);

      if ( verbose_ ) cout << "@@@@ Testing fileName " << foundFileName << endl;
      if ( !rootFilePattern.Match(foundFileName) ) continue;

      string fileFullPath;
      if ( isRfio ) fileFullPath = "rfio:";
      fileFullPath += baseDir;
      if ( baseDir[baseDir.size()-1] != '/' ) fileFullPath += "/";
      fileFullPath += foundFileName;

      inputFiles.push_back(fileFullPath);
    }
  }
}

void FWLiteAnalyzerBase::ListDataFiles()
{
  cout << "@@ List of signal samples\n";
  for ( FileMap::const_iterator dataMapIter = signalFiles_.begin();
        dataMapIter != signalFiles_.end(); ++dataMapIter )
  {
    const string name = dataMapIter->first;
    const vector<string> files = dataMapIter->second;

    cout << "@@@ Dataset = " << name << endl;
    for ( vector<string>::const_iterator filesIter = files.begin();
          filesIter != files.end(); ++filesIter )
    {
      cout << *filesIter << endl;
    }
  }

  cout << "@@ List of background samples\n";
  for ( FileMap::const_iterator dataMapIter = backgroundFiles_.begin();
        dataMapIter != backgroundFiles_.end(); ++dataMapIter )
  {
    string name = dataMapIter->first;
    vector<string> files = dataMapIter->second;

    cout << "@@@ Dataset = " << name << endl;
    for ( vector<string>::const_iterator filesIter = files.begin();
          filesIter != files.end(); ++filesIter )
    {
      cout << *filesIter << endl;
    }
  }
}

void FWLiteAnalyzerBase::ProcessEvent()
{
  if ( verbose_ ) cout << "@@ Processing signal events\n";
  for ( FileMap::const_iterator fileItem = signalFiles_.begin();
        fileItem != signalFiles_.end(); ++fileItem )
  {
    const string name = fileItem->first;
    const vector<string> files = fileItem->second;

    if ( verbose_ ) cout << "@@@ Processing dataset " << name << endl;
    Analyze(name, files);
  }

  if ( verbose_ ) cout << "@@ Processing background events\n";
  for ( FileMap::const_iterator fileItem = backgroundFiles_.begin();
        fileItem != backgroundFiles_.end(); ++fileItem )
  {
    const string name = fileItem->first;
    const vector<string> files = fileItem->second;

    if ( verbose_ ) cout << "@@@ Processing dataset " << name << endl;
    Analyze(name, files);
  }
}

