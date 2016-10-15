#include "GenParticleSelection.h"

GenSelection::GenSelection() { }
GenSelection::~GenSelection() { }


void GenSelection::GenLepSelection(std::vector<float> Eta, std::vector<float> Phi, std::vector<float> Pt, std::vector<float> Energy, std::vector<int> PdgId, std::vector<int> Status, std::vector<int> MotherId, std::vector<bool> Ishardprocess, std::vector<bool> Fromhardprocess, std::vector<Lepton>& leptonColl) {
  fakeType = Lepton::unknown;
  looseTight = Lepton::Other;
  dummy = 9999.;
  leptonColl_tmp.clear();
  for (UInt_t ipart=0; ipart<Pt.size(); ipart++) {
    //simple checks to avoid crashes   
    if(MotherId[ipart] <= 0) continue;
    if(MotherId[ipart] >= Pt.size()) continue;
    if(Pt[ipart] < 0.1) continue;
    //we want status 1 leptons
    if(Status[ipart]!=1) continue;
    //only muons and electrons
    if(fabs(PdgId[ipart])==13) leptonType = Lepton::Muon;
    //else if (fabs(PdgId[ipart])==11) leptonType = Lepton::Electron;
    else continue;
    ///looking at ancestor only muon 
    mindex = ipart;
    while ( fabs(PdgId[mindex]) == 13) {
      if(MotherId[mindex] >= Pt.size()) break;
      mindex = MotherId[mindex];
    }
    //if () nDauther++;
    if (Status[mindex]==2 && fabs(PdgId[mindex])>50) continue;
    if (!(Status[mindex]==2 && fabs(PdgId[mindex])==15) && (Fromhardprocess[ipart]==0 && Fromhardprocess[mindex]==0)) continue;
    charge = PdgId[ipart]/abs(PdgId[ipart]);
    eta = Eta[ipart];
    double ID = PdgId[ipart];
    if (Pt[ipart] >= pt_cut_min && Pt[ipart] < pt_cut_max && fabs(Eta[ipart]) < eta_cut) { 
      vPart.SetPtEtaPhiE(Pt[ipart], Eta[ipart], Phi[ipart], Energy[ipart]);
      leptonColl_tmp.push_back( Lepton(leptonType, ipart, vPart, eta, mindex, dummy, dummy, charge, fakeType, looseTight, ID ) );
    }
  }
  if (leptonColl_tmp.size()>1) {
    for (UInt_t i=0; i<leptonColl_tmp.size()-1; i++) {
      notinclude = false;
      for (UInt_t j=i+1; j<leptonColl_tmp.size(); j++) {
        if (leptonColl_tmp[i].chiNdof()==leptonColl_tmp[j].chiNdof()) {
          if (Fromhardprocess[leptonColl_tmp[i].ilepton()] == 0 )
            notinclude = false; 
        }
      }
      if (!notinclude) leptonColl.push_back(leptonColl_tmp[i]);
    }
  }
  else if (leptonColl_tmp.size()==1)
    leptonColl.push_back(leptonColl_tmp[0]);
  else
    cout << " no gen " << endl;
  // std::sort( jetColl.begin(), jetColl.end(), JetPTSorter );
}

