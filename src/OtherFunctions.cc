#include "OtherFunctions.h"

float getMT2bb(std::vector<Jet>& jets, std::vector<Lepton>& leptons, float MET, float MET_phi) {
  float METx = MET*TMath::Cos(MET_phi);
  float METy = MET*TMath::Sin(MET_phi);
  TLorentzVector Lep0 = leptons[0].lorentzVec();
  TLorentzVector Lep1 = leptons[1].lorentzVec();
  TLorentzVector BtagJet0 = jets[0].lorentzVec();
  TLorentzVector BtagJet1 = jets[1].lorentzVec();
  return getMT2_80( BtagJet0, BtagJet1, sqrt(pow(METx+(Lep0+Lep1).Px(),2)+pow(METy+(Lep0+Lep1).Py(),2)), TMath::ATan2(METy+(Lep0+Lep1).Py(),METx+(Lep0+Lep1).Px()) );
}

float getMT2lblb(std::vector<Jet>& jets, std::vector<Lepton>& leptons, float MET, float MET_phi) {
  float MT2llbb00, MT2llbb01, MT2llbb;
  float METx = MET*TMath::Cos(MET_phi);
  float METy = MET*TMath::Sin(MET_phi);
  TLorentzVector Lep0 = leptons[0].lorentzVec();
  TLorentzVector Lep1 = leptons[1].lorentzVec();
  TLorentzVector BtagJet0 = jets[0].lorentzVec();
  TLorentzVector BtagJet1 = jets[1].lorentzVec();
  TLorentzVector LepPlusBtagJet00 = Lep0 + BtagJet0;
  TLorentzVector LepPlusBtagJet10 = Lep1 + BtagJet0;
  TLorentzVector LepPlusBtagJet11 = Lep1 + BtagJet1;
  TLorentzVector LepPlusBtagJet01 = Lep0 + BtagJet1;
  if (LepPlusBtagJet11.M()<173 && LepPlusBtagJet00.M()<173 && (LepPlusBtagJet10.M()>173 || LepPlusBtagJet01.M()>173))
    MT2llbb=getMT2(LepPlusBtagJet00, LepPlusBtagJet11, MET, MET_phi);
  else if ((LepPlusBtagJet11.M()>173 || LepPlusBtagJet00.M()>173) && LepPlusBtagJet10.M()<173 && LepPlusBtagJet01.M()<173)
    MT2llbb=getMT2(LepPlusBtagJet01, LepPlusBtagJet10, MET, MET_phi);
  else if (LepPlusBtagJet11.M()<173 && LepPlusBtagJet00.M()<173 && LepPlusBtagJet10.M()<173 && LepPlusBtagJet01.M()<173) {
    if ( fabs(LepPlusBtagJet11.M()-LepPlusBtagJet00.M()) < fabs(LepPlusBtagJet10.M()-LepPlusBtagJet01.M()) )
      MT2llbb=getMT2(LepPlusBtagJet00, LepPlusBtagJet11, MET, MET_phi);
    else
      MT2llbb=getMT2(LepPlusBtagJet01, LepPlusBtagJet10, MET, MET_phi);
  }
  else
    MT2llbb=0;
  return MT2llbb;
}

float getMT2(TLorentzVector lept1, TLorentzVector lept2, float theMET, float theMETphi) {
  // Calculate MT2 variable for two leptons and missing energy, assuming zero testmass
  double pa[3];
  double pb[3];
  double pmiss[3];

  TLorentzVector pmet;
  pmet.SetPtEtaPhiM(theMET, 0., theMETphi, 0.);
  pmiss[0] = 0.; // irrelevant
  pmiss[1] = pmet.Px();
  pmiss[2] = pmet.Py();

  pa[0] = lept1.M(); // setting the mass, just in case
  pa[1] = lept1.Px();
  pa[2] = lept1.Py();

  pb[0] = lept2.M(); // setting the mass, just in case
  pb[1] = lept2.Px();
  pb[2] = lept2.Py();

  mt2bisect* MT2bisect = new mt2bisect();
  MT2bisect->set_verbose(0);
  MT2bisect->set_momenta(pa, pb, pmiss);
  MT2bisect->set_mn(0.); // testmass
  double MT2 = MT2bisect->get_mt2();
  delete MT2bisect;
  return MT2;
}

float getMT2_80(TLorentzVector lept1, TLorentzVector lept2, float theMET, float theMETphi) {
  // Calculate MT2 variable for two leptons and missing energy, assuming zero testmass
  double pa[3];
  double pb[3];
  double pmiss[3];

  TLorentzVector pmet;
  pmet.SetPtEtaPhiM(theMET, 0., theMETphi, 0.);
  pmiss[0] = 0.; // irrelevant
  pmiss[1] = pmet.Px();
  pmiss[2] = pmet.Py();

  pa[0] = lept1.M(); // setting the mass, just in case
  pa[1] = lept1.Px();
  pa[2] = lept1.Py();

  pb[0] = lept2.M(); // setting the mass, just in case
  pb[1] = lept2.Px();
  pb[2] = lept2.Py();

  mt2bisect* MT2bisect = new mt2bisect();
  MT2bisect->set_verbose(0);
  MT2bisect->set_momenta(pa, pb, pmiss);
  MT2bisect->set_mn(80.398); // testmass
  double MT2 = MT2bisect->get_mt2();
  delete MT2bisect;
  return MT2;
}

bool LeptonPTSorter(Lepton lep1, Lepton lep2) {
  return lep1.lorentzVec().Pt() > lep2.lorentzVec().Pt();
}

bool LeptonIsoSorter(Lepton lep1, Lepton lep2) {
  return lep1.relIso() < lep2.relIso();
}

bool JetPTSorter(Jet jet1, Jet jet2) {
  return jet1.lorentzVec().Pt() > jet2.lorentzVec().Pt();
}
//bool LeptonIsolationSorter(Lepton lep1, Lepton lep2) {
//////  return (lep1.relIso() < 0.1) > (lep2.relIso() < 0.1);
////// true if lep1 earlier-sorted lepton is isolated
//////}


void getFakerate(TH2F* h1, TH2F* h2, TH2F* out, int nbinX, int nbinY) {
   double frate,ferr;
   for (int i=1; i<nbinX+1; i++)
      for (int j=1; j<nbinY+1; j++){
         double a = h1->GetBinContent(i,j);
         double b = h2->GetBinContent(i,j);
         if (b){
            frate = ((double) a)/((double) b);
            ferr = sqrt( frate*(1-frate)/(double) b);
            out->SetBinContent(i,j,frate);
            out->SetBinError(i,j,ferr);
         }
         else {
            out->SetBinContent(i,j,0);
         }
   }
}

void getFakerate(TH1F* h1, TH1F* h2, TH1F* out, int nbinX) {
   double frate,ferr;
   for (int i=1; i<nbinX+1; i++) {
     double a = h1->GetBinContent(i);
     double b = h2->GetBinContent(i);
       if (b){
         frate = ((double) a)/((double) b);
         ferr = sqrt( frate*(1-frate)/(double) b);
         out->SetBinContent(i,frate);
         out->SetBinError(i,ferr);
       }
       else {
         out->SetBinContent(i,0);
       }
   }
}

double DoubleTOSinglebkg(TH2F* fakerate, std::vector<Lepton>& leptonColl, UInt_t &ilep, UInt_t &jlep) {
  double bkg=0;

  int i0=1; int j0=1;
  int i1=1; int j1=1;
  double eta0= fabs(leptonColl[ilep].eta());
  double Pt0 = leptonColl[ilep].lorentzVec().Pt();
  double eta1= fabs(leptonColl[jlep].eta());
  double Pt1 = leptonColl[jlep].lorentzVec().Pt();
  if (eta0<0.0 || eta0>=2.5) {cout<<"CACCHIO eta!!!! "<<eta0<<endl; eta0<0.0 ? eta0=0.0 : eta0=2.49;}
  if (Pt0<10.0) {cout<<"CACCHIO Pt!!!! "<<Pt0<<endl; Pt0=10.0;}
  if (eta1<0.0 || eta1>=2.5) {cout<<"CACCHIO eta!!!! "<<eta1<<endl; eta1<0.0 ? eta1=0.0 : eta1=2.49;}
  if (Pt1<10.0) {cout<<"CACCHIO Pt!!!! "<<Pt1<<endl; Pt1=10.0;}

  while(1) {
    if( arrayeta[(i0-1)%(ninteta+1)]<=eta0 && eta0<arrayeta[i0%(ninteta+1)] )
      break;
    i0++;
  }
  if (Pt0>=arraypT[nintpT]) j0=nintpT;
  else {
    while(1) {
      if( arraypT[(j0-1)%(nintpT+1)]<=Pt0 && Pt0<arraypT[j0%(nintpT+1)] )
        break;
      j0++;
    }
  }

  while(1) {
    if( arrayeta[(i1-1)%(ninteta+1)]<=eta1 && eta1<arrayeta[i1%(ninteta+1)] )
      break;
    i1++;
  }
  if (Pt1>=arraypT[nintpT]) j1=nintpT;
  else {
    while(1) {
      if( arraypT[(j1-1)%(nintpT+1)]<=Pt1 && Pt1<arraypT[j1%(nintpT+1)] )
        break;
      j1++;
    }
  }

  double A = fakerate->GetBinContent(i0,j0);
  double B = fakerate->GetBinContent(i1,j1);

  bkg = (-A/(1-A))* (A*(1-B)+B*(1-A)) / (( 1-A )*( 1-B ));

  return bkg;
}

double SinglebackGround(TH2F* fakerate, std::vector<Lepton>& leptonColl, UInt_t &ilep, Double_t ***fakeN, UInt_t &type, Double_t weight) {
  double bkg=0;

  int i=1; int j=1;
  double eta= fabs(leptonColl[ilep].eta());
  double Pt = leptonColl[ilep].lorentzVec().Pt();
  if (eta<0.0 || eta>=2.5) {cout<<"CACCHIO eta!!!! "<<eta<<endl; eta<0.0 ? eta=0.0 : eta=2.49;}
  if (Pt<10.0) {cout<<"CACCHIO Pt!!!! "<<Pt<<endl; Pt=10.0;}

  while(1) {
    if( arrayeta[(i-1)%(ninteta+1)]<=eta && eta<arrayeta[i%(ninteta+1)] )
      break;
    i++;
  }
  if (Pt>=arraypT[nintpT]) j=nintpT;
  else {
    while(1) {
      if( arraypT[(j-1)%(nintpT+1)]<=Pt && Pt<arraypT[j%(nintpT+1)] )
        break;
      j++;
    }
  }
  fakeN[0][i-1][j-1]+=weight;
  fakeN[type][i-1][j-1]+=weight;

  bkg = fakerate->GetBinContent(i,j) /( 1-fakerate->GetBinContent(i,j) );
  //errbkg = fakerate->GetBinError(i,j) / pow( 1-fakerate->GetBinContent(i,j),2 );

  return bkg;
}

double DoublebackGround(TH2F* fakerate, std::vector<Lepton>& leptonColl, UInt_t &ilep, UInt_t &jlep, Double_t *****fakeN, UInt_t &type, Double_t weight) {
  double bkg=0;

  int i0=1; int j0=1;
  int i1=1; int j1=1;
  double eta0= fabs(leptonColl[ilep].eta());
  double Pt0 = leptonColl[ilep].lorentzVec().Pt();
  double eta1= fabs(leptonColl[jlep].eta());
  double Pt1 = leptonColl[jlep].lorentzVec().Pt();
  if (eta0<0.0 || eta0>=2.5) {cout<<"CACCHIO eta!!!! "<<eta0<<endl; eta0<0.0 ? eta0=0.0 : eta0=2.49;}
  if (Pt0<10.0) {cout<<"CACCHIO Pt!!!! "<<Pt0<<endl; Pt0=10.0;}
  if (eta1<0.0 || eta1>=2.5) {cout<<"CACCHIO eta!!!! "<<eta1<<endl; eta1<0.0 ? eta1=0.0 : eta1=2.49;}
  if (Pt1<10.0) {cout<<"CACCHIO Pt!!!! "<<Pt1<<endl; Pt1=10.0;}

  while(1) {
    if( arrayeta[(i0-1)%(ninteta+1)]<=eta0 && eta0<arrayeta[i0%(ninteta+1)] )
      break;
    i0++;
  }
  if (Pt0>=arraypT[nintpT]) j0=nintpT;
  else {
    while(1) {
      if( arraypT[(j0-1)%(nintpT+1)]<=Pt0 && Pt0<arraypT[j0%(nintpT+1)] )
        break;
      j0++;
    }
  }

  while(1) {
    if( arrayeta[(i1-1)%(ninteta+1)]<=eta1 && eta1<arrayeta[i1%(ninteta+1)] )
      break;
    i1++;
  }
  if (Pt1>=arraypT[nintpT]) j1=nintpT;
  else {
    while(1) {
      if( arraypT[(j1-1)%(nintpT+1)]<=Pt1 && Pt1<arraypT[j1%(nintpT+1)] )
        break;
      j1++;
    }
  }

  double A = fakerate->GetBinContent(i0,j0);
  double B = fakerate->GetBinContent(i1,j1);
  //double deltaA = fakerate->GetBinError(i0,j0);
  //double deltaB = fakerate->GetBinError(i1,j1);

  fakeN[0][i0-1][j0-1][i1-1][j1-1]+=weight;
  fakeN[type][i0-1][j0-1][i1-1][j1-1]+=weight;

  bkg = A*B / (( 1-A )*( 1-B ));
  //errbkg = sqrt ( ( B*(1-B)*deltaA+A*(1-A)*deltaB )/( pow(1-A,2)*pow(1-B,2) ) );

  return bkg;
}

void DoubleANDSinglebkg(std::vector<Lepton>& leptonColli, UInt_t &ilep, std::vector<Lepton>& leptonCollj, UInt_t &jlep, Double_t *****fakeN, UInt_t &type) {

  int i0=1; int j0=1;
  int i1=1; int j1=1;
  double eta0= fabs(leptonColli[ilep].eta());
  double Pt0 = leptonColli[ilep].lorentzVec().Pt();
  double eta1= fabs(leptonCollj[jlep].eta());
  double Pt1 = leptonCollj[jlep].lorentzVec().Pt();
  if (eta0<0.0 || eta0>=2.5) {cout<<"CACCHIO eta!!!! "<<eta0<<endl; eta0<0.0 ? eta0=0.0 : eta0=2.49;}
  if (Pt0<10.0) {cout<<"CACCHIO Pt!!!! "<<Pt0<<endl; Pt0=10.0;}
  if (eta1<0.0 || eta1>=2.5) {cout<<"CACCHIO eta!!!! "<<eta1<<endl; eta1<0.0 ? eta1=0.0 : eta1=2.49;}
  if (Pt1<10.0) {cout<<"CACCHIO Pt!!!! "<<Pt1<<endl; Pt1=10.0;}

  while(1) {
    if( arrayeta[(i0-1)%(ninteta+1)]<=eta0 && arrayeta[i0%(ninteta+1)]>eta0 )
      break;
    i0++;
  }
  if (Pt0>=arraypT[nintpT]) j0=nintpT;
  else {
    while(1) {
      if( arraypT[(j0-1)%(nintpT+1)]<=Pt0 && arraypT[j0%(nintpT+1)]>Pt0 )
        break;
      j0++;
    }
  }

  while(1) {
    if( arrayeta[(i1-1)%(ninteta+1)]<=eta1 && arrayeta[i1%(ninteta+1)]>eta1 )
      break;
    i1++;
  }
  if (Pt1>=arraypT[nintpT]) j1=nintpT;
  else {
    while(1) {
      if( arraypT[(j1-1)%(nintpT+1)]<=Pt1 && arraypT[j1%(nintpT+1)]>Pt1 )
        break;
      j1++;
    }
  }
  fakeN[0][i1-1][j1-1][i0-1][j0-1]+=1;
  fakeN[type][i1-1][j1-1][i0-1][j0-1]+=1;
}

