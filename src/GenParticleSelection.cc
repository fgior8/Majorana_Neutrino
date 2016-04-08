#include "GenParticleSelection.h"

GenSelection::GenSelection() { }
GenSelection::~GenSelection() { }

void GenSelection::GenPartSelection(std::vector<float> Eta, std::vector<float> Phi, std::vector<float> Pt, std::vector<float> Energy, std::vector<int> PdgId, std::vector<int> Status, std::vector<int> MotherId, double Mass, std::vector<Lepton>& leptonColl) {

  fakeType = Lepton::unknown;
  looseTight = Lepton::Other;
  btag = 9999.;
	
  for (UInt_t ipart = 0; ipart < Pt.size(); ++ipart) {
    //  if ((nthdigit( abs(Gen_pdgId[ipart]), 0)==5 || nthdigit( abs(Gen_pdgId[ipart]), 1)==5 || nthdigit( abs(Gen_pdgId[ipart]), 2)==5) )
    //cout << "Gen_pdgId[ipart] " << Gen_pdgId[ipart] << " Gen_mass[ipart] " << Gen_mass[ipart] << " Gen_motherId[ipart] " << Gen_motherId[ipart] <<endl;
    pt_cut_min = 10.;
    pt_cut_max = 1000000000.;
    eta_cut = 3.0;
    charge = PdgId[ipart]/abs(PdgId[ipart]);
    eta = Eta[ipart];
    if (Pt[ipart] >= pt_cut_min && Pt[ipart] < pt_cut_max && fabs(Eta[ipart]) < eta_cut) { 
      vPart.SetPtEtaPhiM(Pt[ipart], Eta[ipart], Phi[ipart], Mass);
      if (fabs(PdgId[ipart])==13 && Status[ipart]==1) {
	leptonType = Lepton::Muon;
	leptonColl.push_back( Lepton(leptonType, ipart, vPart, eta, btag, btag, btag, charge, fakeType, looseTight, btag) );
      }
    }
  }
  
  // std::sort( jetColl.begin(), jetColl.end(), JetPTSorter );
}


/*
void GenSelection::GenJetSelection(Double_t *Gen_px, Double_t *Gen_py, Double_t *Gen_pz, Double_t *Gen_energy, Int_t *Gen_Flavour, vector<Jet>& bquarkColl, vector<Jet>& RecojetColl, vector<Jet>& GenjetColl) {

  for (UInt_t i = 0; i < RecojetColl.size(); i++) {
    partOK = false;
    UInt_t idx = RecojetColl[i].ijet();
    if (Gen_energy[idx]<=0.) continue;
    vPart.SetPxPyPzE(Gen_px[idx], Gen_py[idx], Gen_pz[idx], Gen_energy[idx]);
    eta = vPart.Eta();
    btag = Gen_Flavour[idx];
    for (UInt_t j=0; j<bquarkColl.size(); j++)
      if (vPart.DeltaR( bquarkColl[j].lorentzVec() ) < 0.1 )
	partOK = true;
    if (fabs(eta)<6.0 && partOK)
      GenjetColl.push_back( Jet(vPart, eta, btag, idx) );
  }
    
}
*/
