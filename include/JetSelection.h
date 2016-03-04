#ifndef JetSelection_h
#define JetSelection_h

#include <iostream>
using namespace std;

#include "TLorentzVector.h"
#include <vector>
#include "Jet.h"
#include "Lepton.h"
#include "OtherFunctions.h"

class JJ {

  TLorentzVector vJet;
  Bool_t jetIsOK;
  Double_t pt_cut_min, pt_cut_max, eta_cut, bdisc_cut;
  
 public:
  JJ();
  ~JJ();
 
  void JetSelection(std::vector<Bool_t> pfJet, std::vector<Double_t> Pt, std::vector<Double_t> Eta, std::vector<Double_t> Phi, std::vector<Double_t> E, std::vector<Double_t> BTag, std::vector<Jet>& jetColl);
  
  void JetSelectionLeptonVeto(std::vector<Bool_t> pfJet, std::vector<Double_t> Pt, std::vector<Double_t> Eta, std::vector<Double_t> Phi, std::vector<Double_t> E, std::vector<Double_t> BTag, std::vector<Lepton>& leptonColl1, std::vector<Lepton>& leptonColl2, std::vector<Jet>& jetColl);

  void SetPt(Double_t minPt, Double_t maxPt);
  void SetPt(Double_t minPt);
  void SetEta(Double_t Eta);
  void SetBdisc(Double_t Bdisc);
 
};

#endif
