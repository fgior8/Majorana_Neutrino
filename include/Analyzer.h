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

  //SF parametrization
  TFile *MuSF_trig, *ElSF_trig, *MuElSF_trig;
  TFile *MuSF_ID, *MuSF_ISO, *ElSF_IDISO;
  TH2F *hmuIDSF, *hmuISOSF, *hmumuTriggerSF;
  TH2F *heIDSF, *heeTriggerSF;
  TH2F *hmueTriggerSF;
  //SF

  Bool_t isData, MCatNLO;
  Bool_t isBtag, isBtagVeto;
  TString completename, treename;
  TFile *outfile;
  Long64_t entrieslimit;
  ReweightPU *reweightPU;
  BTagSFUtil *lBTagSF, *hBTagSF, *slBTagSF, *shBTagSF;
  
  TDirectory *Dir;
  Int_t cut, channel, index;
  Bool_t *goodVerticies;
  Double_t HT, ST, MET, MET_phi;
  Double_t temp_pt;
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
  std::vector<Jet> jetColl; std::vector<Jet> jetCollTop;
  std::vector<Lepton> muonSelected;
  std::vector<Lepton> electronColl;
  std::vector<Lepton> genColl;
  std::vector<Lepton> muonCollMatch;

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

 public:

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
