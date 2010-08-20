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
  lumi_ = 1.0;
}

void FWLiteAnalyzerBase::SetLumi(const double lumi)
{
  lumi_ = lumi;
}

void FWLiteAnalyzerBase::AddMCSample(const string name, const string inputName, const double xsec, const long nGenEvent)
{
  vector<string>& inputFiles = mcSamples_[name];
  mcScaleFactors_[name] = xsec/nGenEvent;

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
  cout << "@@ List of samples\n";
  for ( FileMap::const_iterator dataMapIter = mcSamples_.begin();
        dataMapIter != mcSamples_.end(); ++dataMapIter )
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
}

void FWLiteAnalyzerBase::ProcessEvent()
{
  if ( verbose_ ) cout << "@@ Processing signal events\n";
  for ( FileMap::const_iterator fileItem = mcSamples_.begin();
        fileItem != mcSamples_.end(); ++fileItem )
  {
    const string name = fileItem->first;
    const vector<string> files = fileItem->second;

    if ( verbose_ ) cout << "@@@ Processing dataset " << name << endl;
    Analyze(name, files);
  }
}

TDirectory* FWLiteAnalyzerBase::MakeDirectory(const std::string& path)
{
  TDirectory* dir = outFile_;
  string::size_type pos1 = 0, pos2 = 0;

  while ( pos1 != string::npos )
  {
    if ( !dir ) 
    {
      cout << "@@@@ MakeDirectory failed. directory is Null\n";
      return 0;
    }

    pos1 = path.find_first_not_of('/', pos2);
    pos2 = path.find_first_of('/', pos1);

    if ( pos1 == pos2 ) break;

    string subDir = path.substr(pos1, pos2-pos1);
    if ( verbose_ ) cout << "@@@@ mkdir " << subDir << endl;

    TObject* obj = dir->Get(subDir.c_str());
    if ( !obj ) 
    {
      dir = dir->mkdir(subDir.c_str(), subDir.c_str());
    }
    else
    {
      dir = dynamic_cast<TDirectory*>(obj);
    }

    if ( !dir )
    {
      cout << "@@@ MakekDirectory failed. Object already exists and not a directory\n";
      return 0;
    }
  }

  return dir;
}

