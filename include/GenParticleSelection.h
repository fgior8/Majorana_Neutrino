#ifndef GenParticleSelection_h
#define GenParticleSelection_h

#include <iostream>
using namespace std;
#include "OtherFunctions.h"
#include "TLorentzVector.h"
#include "Jet.h"
#include "Lepton.h"
#include <vector>

class GenSelection {

  Double_t pt_cut_min, pt_cut_max, eta_cut;
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

  void GenPartSelection(Int_t nGen, Int_t *Gen_pdgId, Int_t *Gen_status, Int_t * Gen_motherId, Int_t * Gen_grammaId, Int_t * Gen_sourceId, Double_t *Gen_charge, Double_t *Gen_pt, Double_t *Gen_eta, Double_t *Gen_phi, Double_t *Gen_mass, std::vector<Lepton>& electronColl, std::vector<Lepton>& electronNuColl, std::vector<Lepton>& muonColl, std::vector<Lepton>& muonNuColl, std::vector<Jet>& bquarkColl, std::vector<TLorentzVector>& lightquarkColl);
  
  void GenJetSelection(Double_t *Gen_px, Double_t *Gen_py, Double_t *Gen_pz, Double_t *Gen_energy, Int_t *Gen_Flavour, vector<Jet>& bquarkColl, vector<Jet>& RecojetColl, vector<Jet>& GenjetColl);

};

#endif
