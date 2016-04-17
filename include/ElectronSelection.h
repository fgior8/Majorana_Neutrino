#ifndef ElectronSelection_h
#define ElectronSelection_h

#include "LeptonSelection.h"

class ElectronSel : public Lep {
  
  Bool_t ElectronID;
  Double_t ElTkIso, ElEcalIso, ElHcalIso;

 public:
  ElectronSel();
  ~ElectronSel();

  void ElectronSelection(std::vector<Double_t> scEta, std::vector<Double_t> Pt, std::vector<Double_t> Eta, std::vector<Double_t> Phi, std::vector<Double_t> E, std::vector<Double_t> PhotonIso, std::vector<Double_t> NeutralIso, std::vector<Double_t> ChargeIso, std::vector<Double_t> PUpt, std::vector<Int_t> Charge, std::vector<Bool_t> passConversionVeto, std::vector<Bool_t> snuID, std::vector<Double_t> Dxy, std::vector<Double_t> Dz, std::vector<Lepton>& leptonColl); 


};

#endif
