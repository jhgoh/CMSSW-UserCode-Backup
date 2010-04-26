#include <string>

class LeptonTypes
{
public:
  enum
  {
    None = 1, Electron = 2, Muon = 4, Tau = 8
  };
  
  static int getType(const std::string& leptonType)
  {
    const char leptonInitial = leptonType[0];
    if ( leptonInitial == 'e' ) return Electron;
    else if ( leptonInitial == 'm' ) return Muon;
    else if ( leptonInitial == 't' ) return Tau;
    else return None;
  };

  static int getType(const int pdgId)
  {
    const int absPdgId = abs(pdgId);
    if ( absPdgId == 11 ) return Electron;
    else if ( absPdgId == 13 ) return Muon;
    else if ( absPdgId == 15 ) return Tau;
    else return None;
  };
};

