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
 *  $Date: 2008/07/19 19:26:00 $
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

using namespace std;

void process(TObject* obj, TString outPath);

int main(int argc, char** argv)
{
  if ( argc != 3 ) {
    cout << "Usage : " << argv[0] << " source.root destination_directory\n";
    cout << "        " << argv[0] << " source.root:/my/sub/directory ./my/output/directory\n";
    return 1;
  }

  string srcName(argv[1]);
  string outDirName(argv[2]);

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
    canvas.Print(outPath+"/"+obj->GetName()+".gif");
  }
  else if ( objInfo->InheritsFrom("TCanvas") ) {
    // If stored object is TCanvas, just draw it
    TCanvas* canvas = dynamic_cast<TCanvas*>(obj);
    canvas->Print(outPath+"/"+canvas->GetName()+".gif");
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
