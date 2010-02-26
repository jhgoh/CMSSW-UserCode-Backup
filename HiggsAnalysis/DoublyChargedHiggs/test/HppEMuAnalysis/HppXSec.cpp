#include "TGraph.h"
#include "TCanvas.h"
#include "Math/Interpolator.h"
#include "Math/Polynomial.h"

#include <vector>

class HppCurve
{
public:
  HppCurve();
  TGraph* GetTGraph(const unsigned int nSampling=1);
  inline double Eval(const double x);
private:
  ROOT::Math::Interpolator logInterp_;

  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> yLog_;
};

HppCurve::HppCurve()
{
  // Drell-Yan H++ cross section at 10TeV
  // Calculated by M. Spira
  x_.push_back(50.000); y_.push_back(11183.     );
  x_.push_back(100.00); y_.push_back(706.10     );
  x_.push_back(150.00); y_.push_back(158.00     );
  x_.push_back(200.00); y_.push_back(52.039     );
  x_.push_back(250.00); y_.push_back(21.039     );
  x_.push_back(300.00); y_.push_back(9.6677     );
  x_.push_back(350.00); y_.push_back(4.8495     );
  x_.push_back(400.00); y_.push_back(2.5927     );
  x_.push_back(450.00); y_.push_back(1.4551     );
  x_.push_back(500.00); y_.push_back(0.84814    );
  x_.push_back(550.00); y_.push_back(0.50981    );
  x_.push_back(600.00); y_.push_back(0.31425    );
  x_.push_back(650.00); y_.push_back(0.19782    );
  x_.push_back(700.00); y_.push_back(0.12677    );
  x_.push_back(750.00); y_.push_back(0.82489E-01);
  x_.push_back(800.00); y_.push_back(0.54386E-01);
  x_.push_back(850.00); y_.push_back(0.36275E-01);
  x_.push_back(900.00); y_.push_back(0.24439E-01);
  x_.push_back(950.00); y_.push_back(0.16612E-01);
  x_.push_back(1000.0); y_.push_back(0.11380E-01);

  const unsigned int n = y_.size()+1;
  yLog_.reserve(n);
  for(unsigned int i=0; i<n; ++i)
  {
    yLog_.push_back(log(y_[i]));
  }

  logInterp_.SetData(x_, yLog_);
}

TGraph* HppCurve::GetTGraph(const unsigned int nSampling)
{
  const unsigned int n = x_.size()*nSampling;
  TGraph* grp;

  if ( nSampling < 2 )
  {
    // if Number of sampling < 2 then
    // return TGraph with the original data points
    grp = new TGraph(n, &x_[0], &y_[0]);
  }
  else
  {
    // Do the interpolation
    grp = new TGraph(n);

    const double xMin = x_[0];
    const double xMax = x_[x_.size()-1];
    const double dx = (xMax-xMin)/n;

    for(unsigned int i=0; i<n; ++i)
    {
      const double x = xMin + dx*i;
      const double y = Eval(x);

      grp->SetPoint(i, x, y);
    }
  }

  return grp;
}

double HppCurve::Eval(const double x)
{
  return exp(logInterp_.Eval(x));
}

///

TGraph* HppXSec()
{
  HppCurve hppXSecNLO;

  TCanvas* c = new TCanvas;
  c->SetLogy();

  TGraph* grpXSec = hppXSecNLO.GetTGraph(5);
  grpXSec->SetLineColor(kBlue);
  grpXSec->SetMarkerColor(kBlue);
  grpXSec->SetMarkerStyle(7);
  grpXSec->SetMarkerSize(2);

  TGraph* grpXSecOrig = hppXSecNLO.GetTGraph(1);
  grpXSecOrig->SetLineColor(kRed);
  grpXSecOrig->SetMarkerColor(kRed);
  grpXSecOrig->SetMarkerStyle(2);
  grpXSecOrig->SetMarkerSize(2);

  grpXSec->Draw("AP");
  grpXSecOrig->Draw("P");

  return grpXSec;
}

double getHppXSec(const double hppMass)
{
  static HppCurve hppCurve;
  return hppCurve.Eval(hppMass);
};

