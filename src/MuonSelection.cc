#include "MuonSelection.h"

MuonSel::MuonSel() {};


MuonSel::~MuonSel() {};

void MuonSel::MuonSelection(std::vector<Bool_t> IsPF, std::vector<Bool_t> IsGlobal, std::vector<Double_t> Pt, std::vector<Double_t> Eta,  std::vector<Double_t> Phi, std::vector<Double_t> E, std::vector<Double_t> RelIso, std::vector<Int_t> Charge, std::vector<Int_t> ValidHits, std::vector<Int_t> PixelValidHits, std::vector<Int_t> ValidStations, std::vector<Int_t> LayersWithMeasurement, std::vector<Double_t> GlobalChi2, std::vector<Double_t> Dxy, std::vector<Double_t> Dz, std::vector<Lepton>& leptonColl) {

  for (UInt_t ilep=0; ilep<Pt.size(); ilep++) {
    LeptonchiNdof = GlobalChi2[ilep]; 
    dz = Dxy[ilep];
    dxy = Dz[ilep];   
    vLepton.SetPtEtaPhiE(Pt[ilep], Eta[ilep], Phi[ilep], E[ilep]);

    fakeType = Lepton::unknown;
    looseTight = Lepton::Other;
    leptonType = Lepton::Muon;
 
    if (Pt[ilep] > 0.01)
      LeptonRelIso = RelIso[ilep];
    else LeptonRelIso = 9999.;
    
    (IsPF[ilep] && IsGlobal[ilep] && ValidHits[ilep]>0 && PixelValidHits[ilep]>0 && ValidStations[ilep]>1 && LayersWithMeasurement[ilep]>5) ? individual = true :individual = false;

    (fabs(Eta[ilep]) < eta_cut && Pt[ilep] >= pt_cut_min && Pt[ilep] < pt_cut_max) ? etaPt=true : etaPt =false;

    (LeptonchiNdof<chiNdof_cut && LeptonRelIso < relIso_cut && fabs(dz)<dz_cut && fabs(dxy)<dxy_cut && ( LeptonRelIso >= relIsoMIN_cut || LeptonchiNdof>=chiNdofMIN_cut || fabs(dxy)>=dxyMIN_cut) ) ? RelIsod0Chi2=true : RelIsod0Chi2=false;
    int charge = Charge[ilep];
    if (etaPt && RelIsod0Chi2 && individual) 
      leptonColl.push_back( Lepton(leptonType, ilep, vLepton, Eta[ilep], LeptonchiNdof, dxy, dz, charge, fakeType, looseTight, LeptonRelIso) );
   
  }
  
  std::sort( leptonColl.begin(), leptonColl.end(), LeptonPTSorter );
  
}

void MuonSel::SetDeposits(Double_t ECalDeposit , Double_t HCalDeposit) {
    ECalDeposit ? ECalDeposit_max = ECalDeposit : ECalDeposit_max=400.0;
    HCalDeposit ? HCalDeposit_max = HCalDeposit : HCalDeposit_max=600.0;
    ECalDeposit_min = 0.0;
    HCalDeposit_min = 0.0;
}

void MuonSel::SetDeposits(Double_t ECalDeposit1 , Double_t HCalDeposit1, Double_t ECalDeposit2 , Double_t HCalDeposit2) {
    ECalDeposit1 ? ECalDeposit_min = ECalDeposit1 : ECalDeposit_min=0.0;
    HCalDeposit1 ? HCalDeposit_min = HCalDeposit1 : HCalDeposit_min=0.0;
    ECalDeposit2 ? ECalDeposit_max = ECalDeposit2 : ECalDeposit_max=4.0;
    HCalDeposit2 ? HCalDeposit_max = HCalDeposit2 : HCalDeposit_max=6.0;
}

