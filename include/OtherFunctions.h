#ifndef OtherFunctions_h
#define OtherFunctions_h

#include <iostream>
using namespace std;
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "Data.h"
#include "Lepton.h"
#include "Jet.h"
#include "mt2bisect.h"

float getMT2(TLorentzVector lept1, TLorentzVector lept2, float theMET, float theMETphi);

float getMT2_80(TLorentzVector lept1, TLorentzVector lept2, float theMET, float theMETphi);

float getMT2bb(std::vector<Jet>& jets, std::vector<Lepton>& leptons, float theMET, float theMETphi);

float getMT2lblb(std::vector<Jet>& jets, std::vector<Lepton>& leptons, float theMET, float theMETphi);

bool LeptonPTSorter(Lepton lep1, Lepton lep2);

bool LeptonIsoSorter(Lepton lep1, Lepton lep2);

bool JetPTSorter(Jet jet1, Jet jet2);

void getFakerate(TH2F* h1, TH2F* h2, TH2F* out, int nbinX, int nbinY);

void getFakerate(TH1F* h1, TH1F* h2, TH1F* out, int nbinX);

double DoubleTOSinglebkg(TH2F* fakerate, std::vector<Lepton>& leptonColl, UInt_t &ilep, UInt_t &jlep);

double SinglebackGround(TH2F* fakerate, std::vector<Lepton>& leptonColl, UInt_t &ilep, Double_t ***fakeN, UInt_t &type, Double_t weight);

double DoublebackGround(TH2F* fakerate, std::vector<Lepton>& leptonColl, UInt_t &ilep, UInt_t &jlep, Double_t *****fakeN, UInt_t &type, Double_t weight);

void DoubleANDSinglebkg(std::vector<Lepton>& leptonColli, UInt_t &ilep, std::vector<Lepton>& leptonCollj, UInt_t &jlep, Double_t *****fakeN, UInt_t &type);


static const double arrayeta[] = {0.0,0.8,1.479,2.0,2.5};
static const double arraypT [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};
static const Int_t nintpT=7;
static const Int_t ninteta=4;

#endif
