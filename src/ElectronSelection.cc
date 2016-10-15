#include "ElectronSelection.h"

ElectronSel::ElectronSel() {};


ElectronSel::~ElectronSel() {};


void ElectronSel::ElectronSelection(std::vector<Double_t> scEta, std::vector<Double_t> Pt, std::vector<Double_t> Eta, std::vector<Double_t> Phi, std::vector<Double_t> E, std::vector<Double_t> PhotonIso, std::vector<Double_t> NeutralIso, std::vector<Double_t> ChargeIso, std::vector<Double_t> PUpt, std::vector<Int_t> Charge, std::vector<Bool_t> passConversionVeto, std::vector<Bool_t> snuID, std::vector<Double_t> Dxy, std::vector<Double_t> Dz, std::vector<Lepton>& leptonColl) {

  for (UInt_t ilep=0; ilep<Pt.size(); ilep++) {
    LeptonchiNdof = -99.;
    if ( fabs(scEta[ilep])>1.4442 && fabs(scEta[ilep])<1.566 ) continue;

    vLepton.SetPtEtaPhiE(Pt[ilep], Eta[ilep], Phi[ilep], E[ilep]);

    dz = Dz[ilep];
    dxy = Dxy[ilep];
    
    fakeType = Lepton::unknown;
    looseTight = Lepton::Other;
    leptonType = Lepton::Electron;

    
    if (Pt[ilep] > 0.01)
      LeptonRelIso = ( Charge[ilep] + max( NeutralIso[ilep] + PhotonIso[ilep] - 0.5 * PUpt[ilep], 0.) ) / Pt[ilep];
    else LeptonRelIso = 9999.;
    
    //// electron ID medium WP /// 
    
    snuID[ilep] > 0 ? ElectronID = true : ElectronID = false;
 
    (fabs(Eta[ilep]) < eta_cut && Pt[ilep] >= pt_cut_min && Pt[ilep] < pt_cut_max) ? etaPt=true : etaPt =false;
    
    (LeptonRelIso < relIso_cut && fabs(dz)<dz_cut && fabs(dxy)<dxy_cut && ( LeptonRelIso >= relIsoMIN_cut || fabs(dxy)>=dxyMIN_cut)) ? RelIsod0=true : RelIsod0=false;
    int charge = Charge[ilep];
    if (ElectronID && etaPt && RelIsod0) {
      leptonColl.push_back( Lepton(leptonType, ilep, vLepton, Eta[ilep], LeptonchiNdof, dxy, dz, charge, fakeType, looseTight, LeptonRelIso) );
    }
  }
  
  std::sort( leptonColl.begin(), leptonColl.end(), LeptonPTSorter );
  
}
