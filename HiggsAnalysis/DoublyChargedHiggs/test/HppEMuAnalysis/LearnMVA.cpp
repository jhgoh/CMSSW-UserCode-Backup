#ifndef __CINT__

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"

//#include "TMVA/macros/TMVAGui.C"

#endif

using namespace std;

class LearnMVA
{
public:
  LearnMVA(const TString outFileName)
  {
    outFile_ = TFile::Open(outFileName, "RECREATE");

    factory_ = new TMVA::Factory("LearnMVA", outFile_,
                                 "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G;D");
    nSig_ = 0;
    nBkg_ = 0;
  }

  ~LearnMVA() 
  {
    if ( outFile_ && outFile_->IsOpen() ) outFile_->Close();
  };

  void SetSignal(TChain* eventTree, const double weight);
  void SetBackground(TChain* eventTree, const double weight);

  void Run();
private:
  int nSig_, nBkg_;

  TFile* outFile_;
  TMVA::Factory* factory_;

};

void LearnMVA::SetSignal(TChain* eventTree, const double weight)
{
  factory_->AddSignalTree(eventTree, weight);
  ++nSig_;
}

void LearnMVA::SetBackground(TChain* eventTree, const double weight)
{
  factory_->AddBackgroundTree(eventTree, weight);
  ++nBkg_;
}

void LearnMVA::Run()
{
  if ( nSig_ + nBkg_ < 1 )
  {
    cerr << "No signal or background set up" << endl;
    return;
  }

  // Define input variables
  factory_->AddVariable("ntp.MaxElectronPt", 'F');
  factory_->AddVariable("ntp.MinElectronPt", 'F');
  factory_->AddVariable("ntp.MaxElectronIso", 'F');
  factory_->AddVariable("ntp.MinElectronIso", 'F');

  factory_->AddVariable("ntp.MaxMuonPt", 'F');
  factory_->AddVariable("ntp.MinMuonPt", 'F');
  factory_->AddVariable("ntp.MaxMuonIso", 'F');
  factory_->AddVariable("ntp.MinMuonIso", 'F');

  factory_->AddVariable("ntp.EEMass", 'F');
  factory_->AddVariable("ntp.MuMuMass", 'F');
  factory_->AddVariable("ntp.MaxHiggsPt", 'F');
  factory_->AddVariable("ntp.MinHiggsPt", 'F');

  //factory_->AddSpectator("ntp.HmmMass", "H^{--} Mass", "GeV/c^{2}", 'F');
  //factory_->AddSpectator("ntp.HppMass", "H^{++} Mass", "GeV/c^{2}", 'F');

  TCut cut = "";
      
  factory_->PrepareTrainingAndTestTree(cut, cut, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

/*
  factory_->BookMethod(TMVA::Types::kCuts, "Cuts", 
                       "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
  factory_->BookMethod(TMVA::Types::kCuts, "CutsD", 
                       "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );
  factory_->BookMethod(TMVA::Types::kCuts, "CutsPCA", 
                       "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );
  factory_->BookMethod(TMVA::Types::kCuts, "CutsGA",
                       "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );
  factory_->BookMethod(TMVA::Types::kCuts, "CutsSA",
                       "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
*/

  factory_->BookMethod( TMVA::Types::kLikelihood, "Likelihood", 
                       "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" ); 
  factory_->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD", 
                       "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" ); 
  factory_->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA", 
                       "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 
  //factory_->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE", 
  //                     "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 
  //factory_->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX", 
  //                    "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 
  factory_->BookMethod( TMVA::Types::kKNN, "KNN", 
                       "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
  factory_->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" ); 
  factory_->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=60:NsmoothMVAPdf=10" );
  factory_->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N*2:TestRate=5" );
  factory_->BookMethod( TMVA::Types::kBDT, "BDT", 
                       "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

  factory_->TrainAllMethods();
  factory_->TestAllMethods();
  factory_->EvaluateAllMethods();

}

