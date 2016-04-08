#ifndef Analyzer_h
#define Analyzer_h

#include <set>
#include "TTree.h"
#include "Data.h"
#include "TH2F.h"
#include "Reweight.h"
#include "BTagSFUtil.h"
#include "SignalPlots.h"
#include "MuonPlots.h"
#include "ElectronPlots.h"
#include "JetPlots.h"
#include "XSlist.h"

#include "Lepton.h"
#include "OtherFunctions.h"
#include "SelectionFunctions.h"
#include "MuonSelection.h"
#include "ElectronSelection.h"
#include "JetSelection.h"
#include "GenParticleSelection.h"

#include <iostream>
#include <cmath>

using namespace std;

class Analyzer : public Data {

  const Bool_t debug = false; 
  //const Double_t integratedlumi = 199.149; // for Fakes
  const Double_t integratedlumi = 2318.267; // Signal
  const Bool_t GenMatch = true;
  const Double_t Mass_Z = 91.1876;
  const Double_t Mass_W = 80.398;
  const Double_t Mass_Mu = 0.105658; // GeV

  //SF parametrization
  TFile *MuSF_trig, *ElSF_trig, *MuElSF_trig;
  TFile *MuSF_ID, *MuSF_ISO, *ElSF_IDISO;
  TH2F *hmuIDSF, *hmuISOSF, *hmumuTriggerSF;
  TH2F *heIDSF, *heeTriggerSF;
  TH2F *hmueTriggerSF;
  //SF

  Bool_t isData, MCatNLO;
  TString completename, treename;
  TFile *outfile;
  Long64_t entrieslimit;
  ReweightPU *reweightPU;
  BTagSFUtil *fBTagSF;
  
  TDirectory *Dir;
  Int_t cut, channel, index;
  Bool_t *goodVerticies;
  Double_t HT, MET, MET_phi;
  UInt_t numberVertices, VertexN;
  Bool_t Zveto, triggerOK;
  TString trigger;
  Double_t MCweight, weight;
  
  MuonSel Muon;
  GenSelection Gen;
  ElectronSel Electron;
  JJ Jets;

  MuonPlots ***h_muons, ***h_muons_singlefakes, ***h_muons_doublefakes, ***h_muons_totalfakes;
  ElectronPlots ***h_electrons;
  JetPlots ***h_jets, ***h_bjets;
  SignalPlots ***h_signal, ***h_singlefakes, ***h_doublefakes, ***h_totalfakes;

  Double_t *****doubleFake; Double_t ***singleFake; Double_t *****doubleANDsingleFake;
  Double_t *finalbkg1, *finalbkgerror1, *finalbkg2, *finalbkgerror2, *realsingle, *realsingleerror, *realdouble, *realtotal, *doubletosingle, *errdoubletosingle;

  UInt_t ncuts, nchannels;
  vector<Int_t> selectionStep;
  vector<Int_t> selectChannel;
  std::vector<Lepton> muonLooseNotTightColl;
  std::vector<Lepton> muonLooseColl;
  std::vector<Lepton> muonColl;
  std::vector<Jet> jetColl;
  std::vector<Lepton> muonSelected;
  std::vector<Lepton> electronColl;
  std::vector<Lepton> genColl;

  static const Int_t nintpT=9;
  Double_t *arraypT;
  static const Int_t ninteta=4;
  TH1F *h_dileptMass, *h_MT, *h_MET, *h_HT, *h_METPhi, *h_dr, *h_dPhi, *h_ptMuOverptJet;
  TH2F *h_nEvents, *h_nEventsFO, *h_FOrate;
  MuonPlots *h_muonsFR, *h_muonsLooseFR, *h_TLnum, *h_TLden;
  ElectronPlots *h_electronsFR;
  JetPlots *h_jetsFR;
    
  TH1F *h_prova;
  TH1F *h_VertexNoReweight, *h_VertexPostReweight;

  ///Tree for optimization//////
  TFile *outfileTree;
  // Trees
  TTree *AnalysisTree;
  // Branches
  int TNPV;
  int TNJets;     
  int TNJetsBtag;
  int TNMuon;
  int TNElec;
  float TWeight;
  float TMET;
  float TMET_Phi;
  float THT;
  float TMT2ll;
  float TMT2bb;
  float TMT2lblb;
  TLorentzVector TMuon_Lorentz[20];
  TLorentzVector TElec_Lorentz[20];
  TLorentzVector TJet_Lorentz[50];
  Float_t TJet_discriminant[50];
  TLorentzVector TBJet_Lorentz[30];
  Float_t TMuon_Charge[20];
  Float_t TElec_Charge[20];

  //done optimization tree variables

 public:
  static const Bool_t MC_pu = false; 


  Analyzer();
  ~Analyzer();
  void Loop();
  void LoopFR();
  void LoopQFlip();
  void SetWeight(TString name);
  void SetName(TString name, Int_t version);
  void SetEvtN(Long64_t events);

};
#endif
