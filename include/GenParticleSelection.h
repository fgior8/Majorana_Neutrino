#ifndef GenParticleSelection_h
#define GenParticleSelection_h

#include "OtherFunctions.h"
#include "TLorentzVector.h"
#include "Jet.h"
#include "LeptonSelection.h"
#include <iostream>
using namespace std;
//#include <vector>

class GenSelection : public Lep {

  TLorentzVector vPart;
  Bool_t partOK, notinclude;
  Double_t mindex, dummy, eta;
  Int_t charge;
  Lepton::FakeType fakeType;
  Lepton::LooseTight looseTight;
  Lepton::LeptonType leptonType;
  vector<Lepton> leptonColl_tmp;

 public:
  GenSelection();
  ~GenSelection();

  void GenLepSelection();
  void GenLepSelection(std::vector<float> Eta, std::vector<float> Phi, std::vector<float> Pt, std::vector<float> Energy, std::vector<int> PdgId, std::vector<int> Status, std::vector<int> MotherId, std::vector<bool> Ishardprocess, std::vector<bool> Fromhardprocess, std::vector<Lepton>& leptonColl);

};

#endif
