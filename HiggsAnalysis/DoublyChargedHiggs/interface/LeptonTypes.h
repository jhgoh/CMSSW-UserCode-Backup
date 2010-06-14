#ifndef HiggsAnalysis_DoublyChargedHiggs_LeptonTypes_H
#define HiggsAnalysis_DoublyChargedHiggs_LeptonTypes_H

#include <string>

class LeptonTypes
{
public:
  enum
  {
    None = 1<<0, Electron = 1<<1, Muon = 1<<2, Tau = 1<<3
  };

  static int getType(const std::string& leptonName)
  {
    const char leptonInitial = std::toupper(leptonName[0]);

    switch ( leptonInitial )
    {
      case 'E': return Electron; break;
      case 'M': return Muon    ; break;
      case 'T': return Tau     ; break;
      default : return None;
    }
  };

  static int getType(const int pdgId)
  {
    const int absPdgId = abs(pdgId);

    switch ( absPdgId )
    {
      case 11: return Electron; break;
      case 13: return Muon    ; break;
      case 15: return Tau     ; break;
      default: return None;
    }
  };
};

#endif
