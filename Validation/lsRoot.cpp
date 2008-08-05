/**
 *
 *  Class : lsRoot
 *
 *  Description:
 *  List up content of root file
 *
 *  lsRoot myfile.root
 *
 *  $Date: 2008/08/05 13:02:41 $
 *  $Revision: 1.1 $
 *
 *  Authors:
 *  \authro Junghwan Goh - SungKyunKwan University
 *
 **/
#include <iostream>

#include "TROOT.h"
#include "TClass.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TKey.h"
#include "Riostream.h"

using namespace std;

void process(TObject* obj, TString path);

int main(int argc, char** argv)
{
  if ( argc < 2 ) {
    cout << "Usage : " << argv[0] << " source.root\n"
         << " Options :\n" 
         << "           --format=(txt|json)    default is txt\n";
    return 1;
  }

  string srcName(argv[1]);

  // Set source root file and TDirectory
  TFile* srcFile = 0;

  srcFile = new TFile(srcName.c_str(), "r");

  if ( !srcFile || srcFile->IsZombie() ) {
    cout << "Invalid source root file: " << srcName << endl;
    return 1;
  }

  process(srcFile, "/");

  if ( srcFile ) {
    srcFile->Close();
    delete srcFile;
  }

  return 0;
}

void process(TObject* obj, TString path)
{
  if ( !obj ) return;

  TClass* objInfo = obj->IsA();

  cout << objInfo->GetName() << "\t" << path << obj->GetName() << endl;

  if ( objInfo->InheritsFrom("TDirectory") ) {
    // Make new directory and loop over all objects under this directory

    TDirectory* dir = dynamic_cast<TDirectory*>(obj);

    TString newPath = path;
    if ( ! TString(dir->GetName()).EndsWith(".root") ) {
      newPath=path+dir->GetName()+"/";
    }

    // Loop over all objects in this directory
    TIter nextkey(dir->GetListOfKeys());
    TKey* key;
    while ( (key = dynamic_cast<TKey*>(nextkey())) ) {
      // Process recursive
      process(key->ReadObj(), newPath);
    }
  }
}
