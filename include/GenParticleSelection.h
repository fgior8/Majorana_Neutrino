#ifndef GenParticleSelection_h
#define GenParticleSelection_h

#include <iostream>
using namespace std;
#include "OtherFunctions.h"
#include "TLorentzVector.h"
#include "Jet.h"
#include "LeptonSelection.h"
//#include <vector>

class GenSelection : public Lep {

  TLorentzVector vPart;
  Bool_t partOK;
  Double_t btag, eta;
  Int_t charge;
  Lepton::FakeType fakeType;
  Lepton::LooseTight looseTight;
  Lepton::LeptonType leptonType;

 public:
  GenSelection();
  ~GenSelection();

/*
  void GenPartSelection(Int_t nGen, Int_t *Gen_pdgId, Int_t *Gen_status, Int_t * Gen_motherId, Int_t * Gen_grammaId, Int_t * Gen_sourceId, Double_t *Gen_charge, Double_t *Gen_pt, Double_t *Gen_eta, Double_t *Gen_phi, Double_t *Gen_mass, std::vector<Lepton>& electronColl, std::vector<Lepton>& electronNuColl, std::vector<Lepton>& muonColl, std::vector<Lepton>& muonNuColl, std::vector<Jet>& bquarkColl, std::vector<TLorentzVector>& lightquarkColl);
  
  void GenJetSelection(Double_t *Gen_px, Double_t *Gen_py, Double_t *Gen_pz, Double_t *Gen_energy, Int_t *Gen_Flavour, vector<Jet>& bquarkColl, vector<Jet>& RecojetColl, vector<Jet>& GenjetColl);
*/

  void GenLepSelection(std::vector<float> Eta, std::vector<float> Phi, std::vector<float> Pt, std::vector<float> Energy, std::vector<int> PdgId, std::vector<int> Status, std::vector<int> MotherId, std::vector<Lepton>& leptonColl);

};

#endif
