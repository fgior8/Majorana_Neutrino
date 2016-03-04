#ifndef MuonSelection_h
#define MuonSelection_h

#include "LeptonSelection.h"

class MuonSel : public Lep {
  Int_t numVer, leptoni;
  Double_t ECalDeposit_max, HCalDeposit_max, ECalDeposit_min, HCalDeposit_min;

 public:
  MuonSel();
  ~MuonSel();

  void MuonSelection(std::vector<Bool_t> IsPF, std::vector<Bool_t> IsGlobal, std::vector<Double_t> Pt, std::vector<Double_t> Eta,  std::vector<Double_t> Phi, std::vector<Double_t> E, std::vector<Double_t> RelIso, std::vector<Int_t> Charge, std::vector<Int_t> ValidHits, std::vector<Int_t> PixelValidHits, std::vector<Int_t> ValidStations, std::vector<Int_t> LayersWithMeasurement, std::vector<Double_t> GlobalChi2, std::vector<Double_t> Dxy, std::vector<Double_t> Dz, std::vector<Lepton>& leptonColl);

  void SetDeposits(Double_t ECalDeposit, Double_t HCalDeposit);
  void SetDeposits(Double_t ECalDeposit1 , Double_t HCalDeposit1, Double_t ECalDeposit2 , Double_t HCalDeposit2);

};

#endif
