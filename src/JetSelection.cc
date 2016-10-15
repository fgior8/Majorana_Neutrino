#include "JetSelection.h"

JJ::JJ() { }
JJ::~JJ() { }

void JJ::JetSelection (std::vector<Bool_t> pfJet, std::vector<Double_t> Pt, std::vector<Double_t> Eta, std::vector<Double_t> Phi, std::vector<Double_t> E, std::vector<Double_t> BTag, std::vector<Int_t> Flavour, std::vector<Jet>& jetColl) {

  for (UInt_t ijet = 0; ijet < Pt.size(); ++ijet) {
    if (Pt[ijet] >= pt_cut_min && Pt[ijet] < pt_cut_max	&& fabs(Eta[ijet]) < eta_cut && pfJet[ijet]>0) {
      vJet.SetPtEtaPhiE(Pt[ijet], Eta[ijet], Phi[ijet], E[ijet]);;
      jetColl.push_back( Jet(vJet, Flavour[ijet], BTag[ijet], ijet) );
    }
  }
  std::sort( jetColl.begin(), jetColl.end(), JetPTSorter );
}

void JJ::JetSelectionLeptonVeto(std::vector<Bool_t> pfJet, std::vector<Double_t> Pt, std::vector<Double_t> Eta, std::vector<Double_t> Phi, std::vector<Double_t> E, std::vector<Double_t> BTag, std::vector<Int_t> Flavour, std::vector<Lepton>& leptonColl1, std::vector<Lepton>& leptonColl2, std::vector<Jet>& jetColl) {
  std::vector<Jet> pre_jetColl;

  for (UInt_t ijet = 0; ijet < Pt.size(); ijet++) {
    if (Pt[ijet] >= pt_cut_min && Pt[ijet] < pt_cut_max	&& fabs(Eta[ijet]) < eta_cut && pfJet[ijet]>0) {
      vJet.SetPtEtaPhiE(Pt[ijet], Eta[ijet], Phi[ijet], E[ijet]);;
      pre_jetColl.push_back( Jet(vJet, Flavour[ijet], BTag[ijet], ijet) );
    }
  }

  for (UInt_t ijet = 0; ijet < pre_jetColl.size(); ijet++) {
    jetIsOK = true;
    for (UInt_t ilep = 0; ilep < leptonColl1.size(); ilep++) {
      if (pre_jetColl[ijet].lorentzVec().DeltaR( leptonColl1[ilep].lorentzVec() ) < 0.4) {
	jetIsOK = false;
	break;
      }
    }
    for (UInt_t ilep = 0; ilep < leptonColl2.size(); ilep++) {
      if (pre_jetColl[ijet].lorentzVec().DeltaR( leptonColl2[ilep].lorentzVec() ) < 0.4 ) {
	jetIsOK = false;
	break;
      }
    }
    
    if (jetIsOK)
      jetColl.push_back( pre_jetColl[ijet] );
  }
  std::sort( jetColl.begin(), jetColl.end(), JetPTSorter );
}

void JJ::SetPt(Double_t minPt) {
  minPt ? pt_cut_min=minPt : pt_cut_min=0.0;
  pt_cut_max=10000.0;
}

void JJ::SetPt(Double_t minPt, Double_t maxPt) {
  minPt ? pt_cut_min=minPt : pt_cut_min=0.0;
  maxPt ? pt_cut_max=maxPt : pt_cut_max=10000.0;
}

void JJ::SetEta(Double_t Eta) {
  Eta ? eta_cut=Eta : eta_cut=3.0;
}

void JJ::SetBdisc(Double_t Bdisc) {
  Bdisc ? bdisc_cut=Bdisc : bdisc_cut=100000.0;
}
