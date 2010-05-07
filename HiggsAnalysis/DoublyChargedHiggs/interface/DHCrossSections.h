#ifndef HiggsAnalysis_DoublyChargedHiggs_DHCrossSections_H
#define HiggsAnalysis_DoublyChargedHiggs_DHCrossSections_H

#include "TGraph.h"
#include "TCanvas.h"
//#include "Math/Interpolator.h"
//#include "Math/Polynomial.h"

#include <cmath>

class DHCrossSections
{
public:
  enum BeamEnergyType
  {
    _7TeV, _10TeV, _14TeV
  };

public:
  DHCrossSections(const BeamEnergyType beamEnergyType = _10TeV)
  {
    switch ( beamEnergyType )
    {
      case _7TeV : n_ = n7_ ; x_ = x7_ ; y_ = y7_ ; break;
      case _10TeV: n_ = n10_; x_ = x10_; y_ = y10_; break;
      case _14TeV: n_ = n14_; x_ = x14_; y_ = y14_; break;
      default: n_ = 0; x_ = 0; y_ = 0;
    }
  };

  TGraph* GetTGraph(const unsigned int nSampling=1);
  double Eval(const double x);

private:
  //ROOT::Math::Interpolator logInterp_;
  
  int n_;
  const double* x_, * y_;
  
  const static int n7_ = 20;
  const static double x7_[n7_];
  const static double y7_[n7_];

  const static int n10_ = 20;
  const static double x10_[n10_];
  const static double y10_[n10_];

  const static int n14_ = 20;
  const static double x14_[n14_];
  const static double y14_[n14_];
};

// Cross-section table to be filled up
const double DHCrossSections::x14_[] = { 50, 100, 150, 200, 250,
                                        300, 350, 400, 450, 500, 
                                        550, 600, 650, 700, 750, 
                                        800, 850, 900, 950, 1000};
const double DHCrossSections::y14_[] = {};

// PYTHIA calculation by Mario Kadastik
const double DHCrossSections::x7_[] = { 60,  70,  80,  90, 100,
                                       110, 120, 130, 140, 150,
                                       160};

const double DHCrossSections::y7_[] = {2342, 1265.93, 768.03, 499.53, 335.1,
                                       236.87, 168.83, 124.64, 93.84, 71.73, 
                                       55.52};

// Calculated by M. Spira
const double DHCrossSections::x10_[] = { 50, 100, 150, 200, 250,
                                        300, 350, 400, 450, 500, 
                                        550, 600, 650, 700, 750, 
                                        800, 850, 900, 950, 1000};

const double DHCrossSections::y10_[] = {11183, 706.10, 158.00, 52.039, 21.039, 
                                        9.6677, 4.8495, 2.5927, 1.4551, 0.84814, 
                                        0.50981, 0.31425, 0.19782, 0.12677, 0.82489E-01, 
                                        0.54386E-01, 0.36275E-01, 0.24439E-01, 0.16612E-01, 0.11380E-01};

double DHCrossSections::Eval(const double x0)
{
  if ( n_ < 2 ) return -1;

  // Find x bin position
  int bin1 = 0, bin2 = n_-1;
  if ( x0 < x_[bin1] )
  {
    bin2 = 1;
  }
  else if ( x0 > x_[bin2] )
  {
    bin1 = n_-2;
    bin2 = n_-1;
  }
  else
  {
    // Start binary search
    while ( bin2-bin1 > 1 )
    {
      const int testBin = (bin2+bin1)/2;
      const int xAtTestBin = x_[testBin];
      if ( x0 < xAtTestBin ) bin2 = testBin;
      else if ( x0 == xAtTestBin ) return y_[testBin];
      else bin1 = testBin;
    }
  }

  // Now do the log-linear interpolation
  const double x1 = x_[bin1];
  const double x2 = x_[bin2];
  const double logy1 = std::log(y_[bin1]);
  const double logy2 = std::log(y_[bin2]);

  const double logy = logy1+(x0-x1)/(x2-x1)*(logy2-logy1);

  return std::exp(logy);
}

TGraph* DHCrossSections::GetTGraph(const unsigned int nSampling)
{
  const unsigned int n = n_*nSampling;
  TGraph* grp = 0;

  if ( nSampling < 2 )
  {
    grp = new TGraph(n, x_, y_);
  }
  else
  {
    grp = new TGraph(n);

    const double xMin = x_[0];
    const double xMax = x_[n_-1];
    const double dx = (xMax-xMin)/n;

    for ( unsigned int i=0; i<n; ++i )
    {
      const double x = xMin + dx*i;
      const double y = Eval(x);

      grp->SetPoint(i, x, y);
    }
  }

  return grp;
}

#endif
