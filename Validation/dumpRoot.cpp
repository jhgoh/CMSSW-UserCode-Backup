/**
 *
 *  Class : dumpRoot
 *
 *  Description:
 *  Dump root histograms to gif images.
 *  All objects will be stored in the same
 *  directory paths of the source root file.
 *
 *  You can specify sub-directory of interest
 *  appending full path to TDirectory with deliminator ':'
 *
 *  dumpRoot myfile.root:/my/directory/ ./out
 *
 *  $Date: 2008/07/19 14:00:38 $
 *  $Revision: 1.1 $
 *
 *  Authors:
 *  \authro Junghwan Goh - SungKyunKwan University
 *
 **/
#include <iostream>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using namespace std;

vector<string> exts;

void process(TObject* obj, TString outPath);

int main(int argc, char** argv)
{
  using namespace boost;

  program_options::options_description desc("Options for dumpRoot");
  desc.add_options()("infile,i", program_options::value<string>(), "Input file name")
                    ("outdir,o", program_options::value<string>(), "output directory (current directory by default)")
                    ("help,h", "Print this message")
                    ("ext,x", program_options::value<string>(), "output extensions (gif by default)");

  stringstream ss_usage;
  ss_usage << "Usage : " << argv[0] << " -i source.root -o destination_directory\n"
           << "        " << argv[0] << " -i source.root:/my/sub/directory\n"
           << "        " << argv[0] << " -i srouce.root:/my/dir -x gif,pdf\n"
           << "Options : \n"
           << "   --ext=gif,png,pdf,... : List up output extensions (gif turned on by default)\n";
  const string usage = ss_usage.str();

  program_options::positional_options_description pos;
  program_options::variables_map vmap;

  try {
    program_options::store(program_options::command_line_parser(argc, argv).options(desc).positional(pos).run(), vmap);
  }
  catch ( program_options::error const& x ) {
    cerr << "Unable to parse options:\n" << x.what() << "\n\n"
         << desc << usage << endl;
    return 1;
  }

  program_options::notify(vmap);
  if ( vmap.count("help") ) {
    cout << desc << usage << endl;
    return 0;
  }

  string srcName;
  if ( vmap.count("infile") ) {
    srcName = vmap["infile"].as<string>();
  }
  else {
    cout << desc << usage << endl;
    return 0;
  }

  string outDirName;
  if ( vmap.count("outdir") ) {
    outDirName = vmap["outdir"].as<string>();
  }
  else {
    outDirName = ".";
  }

  if ( !vmap.count("ext") ) {
    exts.push_back("gif");
  }
  else {
    string expopts = vmap["ext"].as<string>();
    algorithm::split(exts, expopts, algorithm::is_any_of(","));
  }

  // Set source root file and TDirectory
  TFile* srcFile = 0;
  TDirectory* srcTDir = 0;

  int delimPos = srcName.rfind(".root:");
  if ( delimPos > 0 ) {
    string srcFileName = srcName.substr(0, delimPos+5);
    string srcTDirName = srcName.substr(delimPos+6, srcName.size());

    srcFile = new TFile(srcFileName.c_str(), "r");
    if ( !srcFile || srcFile->IsZombie() ) {
      cout << "Cannot open ROOT file : " << srcFileName << endl;
      return 2;
    }
    srcTDir = srcFile->GetDirectory(srcTDirName.c_str());
  }
  else {
    srcFile = new TFile(srcName.c_str(), "r");
    srcTDir = srcFile;
  }

  if ( !srcTDir || srcTDir->IsZombie() ) {
    cout << "Invalid source TDirectory : " << srcName << endl;
    return 1;
  }

  process(srcTDir, outDirName);

  if ( srcFile ) {
    srcFile->Close();
    delete srcFile;
  }

  return 0;
}

void process(TObject* obj, TString outPath)
{
  if ( !obj ) return;

  TClass* objInfo = obj->IsA();

  if ( objInfo->InheritsFrom("TH1") ||
       objInfo->InheritsFrom("TH2") ) {
    // Draw histograms to an image file
    TCanvas canvas(TString(obj->GetName())+"_canvas", obj->GetTitle());
    obj->Draw();
    for(int i=0; i<exts.size(); i++) {
      canvas.Print(outPath+"/"+obj->GetName()+"."+exts[i].c_str());
    }
  }
  else if ( objInfo->InheritsFrom("TCanvas") ) {
    // If stored object is TCanvas, just draw it
    TCanvas* canvas = dynamic_cast<TCanvas*>(obj);
    for(int i=0; i<exts.size(); i++) {
      canvas->Print(outPath+"/"+obj->GetName()+"."+exts[i].c_str());
    }
  }
  else if ( objInfo->InheritsFrom("TDirectory") ) {
    // Make new directory and loop over all objects under this directory

    TDirectory* dir = dynamic_cast<TDirectory*>(obj);

    // Make new directory
    TString newPath = outPath+"/";
    if ( ! TString(dir->GetName()).EndsWith(".root") ) {
      // Do not create directory with input file's name
      newPath+=dir->GetName();
    }
    gSystem->mkdir(newPath, true);

    // Loop over all objects in this directory
    TIter nextkey(dir->GetListOfKeys());
    TKey* key;
    while ( (key = dynamic_cast<TKey*>(nextkey())) ) {
      // Process recursive
      process(key->ReadObj(), newPath);
    }
  }
  else if ( objInfo->InheritsFrom("TTree") ) {
    // We cannot draw TTree or TNtuple cannot be simply
    cout << "I don't draw TTree" << endl;
    return;
  }
  else {
    cout << "Unknown object :"
         << " type="  << objInfo->GetName()
         << " name=" << obj->GetName()
         << " title=" << obj->GetTitle() << endl;
    return;
  }
}
