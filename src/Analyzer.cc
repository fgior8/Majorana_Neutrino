#include "Analyzer.h"

Analyzer::Analyzer() {

  if (debug) cout<<"inizio"<<endl;
  
  h_prova = new TH1F("h_prova","p_T",100,0,1000);
  h_VertexNoReweight = new TH1F("h_VertexNoReweight","n Vertices no reweighted", 60,0,60);
  h_VertexPostReweight = new TH1F("h_VertexPostReweight","n Vertices reweighted", 60,0,60);
 
  ncuts=9;
  nchannels=5;
  TString cut_name[] = {"particle", "twoLeptons", "twoJets", "ss_Leptons", "ss_oneJet", "ss_twoJets", "ss_Ctrl_MET", "ss_Ctrl_B", "Signal"};
  TString channel_name[] = {"all", "sameF","ee", "mumu", "emu"};
  // 0-all, 1-SF, 2 ee, 3 mumu, 4 emu, (5 mue if needed)
  h_muons = new MuonPlots**[ncuts];
  h_electrons = new ElectronPlots**[ncuts];
  h_jets = new JetPlots**[ncuts];
  //h_bjets = new JetPlots**[ncuts];
  h_signal = new SignalPlots**[ncuts];
  h_singlefakes = new SignalPlots**[ncuts];
  h_doublefakes = new SignalPlots**[ncuts];
  h_totalfakes = new SignalPlots**[ncuts];
  
  for (UInt_t i=0;i<ncuts;i++) {
    h_muons[i] = new MuonPlots*[nchannels];
    h_electrons[i] = new ElectronPlots*[nchannels];
    h_jets[i] = new JetPlots*[nchannels];
    //h_bjets[i] = new JetPlots*[nchannels];
    h_signal[i] = new SignalPlots*[nchannels];
    h_singlefakes[i] = new SignalPlots*[nchannels];
    h_doublefakes[i] = new SignalPlots*[nchannels];
    h_totalfakes[i] = new SignalPlots*[nchannels];
  }
  
  for (UInt_t i=0;i<ncuts;i++)
    for (UInt_t j=0;j<nchannels;j++) {
    h_muons[i][j]     = new MuonPlots    ("muons_"    +cut_name[i]+"_"+channel_name[j]);
    h_electrons[i][j] = new ElectronPlots("electrons_"+cut_name[i]+"_"+channel_name[j]);
    h_jets[i][j]      = new JetPlots     ("jets_"     +cut_name[i]+"_"+channel_name[j]);
    //h_bjets[i][j]     = new JetPlots     ("bjets_"    +cut_name[i]+"_"+channel_name[j]);
    h_signal[i][j]    = new SignalPlots  ("signal_"   +cut_name[i]+"_"+channel_name[j]);
    h_singlefakes[i][j]    = new SignalPlots  ("signal_" +cut_name[i]+"_"+channel_name[j]+"_sf");
    h_doublefakes[i][j]    = new SignalPlots  ("signal_" +cut_name[i]+"_"+channel_name[j]+"_df");
    h_totalfakes[i][j]     = new SignalPlots  ("signal_" +cut_name[i]+"_"+channel_name[j]+"_tf");
      
  }
    
  if (debug) cout<<"fine"<<endl;
}

Analyzer::~Analyzer() { }

void Analyzer::SetName(TString name, Int_t version) {
  completename = name + "_";
  completename += version;
  completename += ".root";
  outfile = new TFile(completename,"RECREATE");
}

void Analyzer::SetWeight(TString name) {

  MCweight = integratedlumi * getXS(name);
  // lumi *  cs(pb) * gen filter efficiency / MCevents
  name.Contains("amcatnlo") ? MCatNLO=true : MCatNLO=false; 
  name.Contains("Data") ? isData=true : isData=false;
  if(isData) MCweight = 1;
  cout<<"MCweight = "<<MCweight<<endl;
}

void Analyzer::SetEvtN(Long64_t events) {
  entrieslimit = events;
  cout<<"Processing "<<events<<endl;
}

void Analyzer::Loop() {
  TH2F *FRhisto;
  TFile *infile = new TFile("histoFR/Total_FRcorr40_1.root");
  infile->cd();
  TDirectory *dir=gDirectory;
  dir->GetObject("h_FOrate3",FRhisto);
  UInt_t nbinX=FRhisto->GetNbinsX(); UInt_t nbinY=FRhisto->GetNbinsY(); UInt_t nSplit=4;

  Double_t SingleFake=0; Double_t DoubleFake=0; Double_t Single_Double=0;
  Int_t nSingleFake=0; Int_t nDoubleFake=0;

  if (debug) cout<< "Something wrong reading the FR histo" <<endl;

  singleFake=new Double_t**[nSplit];
  doubleFake=new Double_t****[nSplit];
  doubleANDsingleFake=new Double_t ****[nSplit];
  finalbkg1=new Double_t[nSplit];
  finalbkgerror1=new Double_t[nSplit];
  finalbkg2=new Double_t[nSplit];
  finalbkgerror2=new Double_t[nSplit];
  realsingle=new Double_t[nSplit];
  realsingleerror=new Double_t[nSplit];
  realdouble=new Double_t[nSplit];
  realtotal=new Double_t[nSplit];
  doubletosingle=new Double_t[nSplit];
  errdoubletosingle=new Double_t[nSplit];
  for (UInt_t z=0; z<nSplit; z++) {
    singleFake[z]=new Double_t*[nbinX];
    doubleFake[z]=new Double_t***[nbinX];
    doubleANDsingleFake[z]=new Double_t***[nbinX];
    finalbkg1[z]=0;
    finalbkgerror1[z]=0;
    finalbkg2[z]=0;
    finalbkgerror2[z]=0;
    realsingle[z]=0;
    realsingleerror[z]=0;
    realdouble[z]=0;
    realtotal[z]=0;
    doubletosingle[z]=0;
    errdoubletosingle[z]=0;
  }
  for (UInt_t z=0; z<nSplit; z++)
    for (UInt_t i=0; i<nbinX; i++) {
      singleFake[z][i]=new Double_t[nbinY];
      doubleFake[z][i]=new Double_t**[nbinY];
      doubleANDsingleFake[z][i]=new Double_t**[nbinY];
    }
  for (UInt_t z=0; z<nSplit; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++) {
        singleFake[z][i][j]=0;
        doubleFake[z][i][j]=new Double_t*[nbinX];
        doubleANDsingleFake[z][i][j]=new Double_t*[nbinX];
      }
  for (UInt_t z=0; z<nSplit; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++)
        for (UInt_t m=0; m<nbinX; m++) {
          doubleFake[z][i][j][m]=new Double_t[nbinY];
          doubleANDsingleFake[z][i][j][m]=new Double_t[nbinY];
        }
  for (UInt_t z=0; z<nSplit; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++)
        for (UInt_t m=0; m<nbinX; m++)
          for (UInt_t n=0; n<nbinY; n++) {
            doubleFake[z][i][j][m][n]=0;
            doubleANDsingleFake[z][i][j][m][n]=0;
          }
  
  cout << "total number of entries " <<nentries<<endl;

  if (debug) cout<< "loop begins" <<endl;

  fBTagSF = new BTagSFUtil("CSVM");

  // once we have data we must look a the pileup
  //  reweightPU = new ReweightPU("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/MyDataPileupHistogram_69400.root");

  if (debug) cout<< "PU histos loaded" <<endl;


  if(!MCweight) MCweight=1; 

  if (fChain == 0) 
    cout << "Ciao!" << endl;

  if (entrieslimit != -1)
    nentries=entrieslimit;

  if (debug) cout<< "at the loop" <<endl;
  std::set<int> runs;
  for (Long64_t jentry = 0; jentry < nentries; jentry++ ) {
    //clearing vectors
    selectionStep.clear();    selectChannel.clear();
    muonColl.clear(); muonLooseColl.clear(); muonLooseNotTightColl.clear();
    electronColl.clear(); jetColl.clear();
    
    if (debug) cout<< "Event number " <<jentry<<endl;
    if (debug) cout<<"begin loop"<<endl;
    if (!(jentry % 50000)) 
      cout << jentry << endl;

    if (!fChain) cout<<"problems with the input file"<<endl;
    fChain->GetEntry(jentry);

    triggerOK = false;
    for(UInt_t t=0; t<vtrignames->size(); t++) {
      trigger = vtrignames->at(t);
      Int_t ps = vtrigps->at(t);
      if ( trigger.Contains("Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") && ps>0) {
        triggerOK = true;
        break;
      }
    }
    if(!triggerOK) continue;
    
    // MET filters for now all OFF
    //----------------------------------------------------------------------------
    if (!(HBHENoiseFilter && CSCTightHaloFilter && eeBadScFilter && EcalDeadCellTriggerPrimitiveFilter)) continue;
    
    weight=MCweight;
    //MC@NLO weight
    if(MCatNLO)
      genWeight>=0 ? weight*=1. : weight*=-1.;
    
    // Vertex Select
    if ( !goodVertices ) continue;

    if(debug) cout<< "object selection" <<endl;
    
    Muon.SetPt(15);
    Muon.SetEta(2.4);
    Muon.SetRelIso(0.15);
    Muon.SetChiNdof(10);
    Muon.SetBSdxy(0.20);
    Muon.SetBSdz(0.50);
    Muon.MuonSelection(*muon_isPF, *muon_isGlobal, *muon_pt, *muon_eta, *muon_phi, *muon_energy, *muon_relIso04, *muon_q, *muon_validhits, *muon_validpixhits, *muon_matchedstations, *muon_trackerlayers, *muon_normchi, *muon_dxy, *muon_dz, muonColl);

    Muon.SetPt(15);
    Muon.SetEta(2.4);
    Muon.SetRelIso(0.6);
    Muon.SetChiNdof(50);
    Muon.SetBSdxy(0.20);
    Muon.SetBSdz(0.50);
    Muon.MuonSelection(*muon_isPF, *muon_isGlobal, *muon_pt, *muon_eta, *muon_phi, *muon_energy, *muon_relIso04, *muon_q, *muon_validhits, *muon_validpixhits, *muon_matchedstations, *muon_trackerlayers, *muon_normchi, *muon_dxy, *muon_dz, muonLooseColl);

    Muon.SetPt(15);
    Muon.SetEta(2.4);
    Muon.SetRelIso(0.15,0.6);
    Muon.SetChiNdof(10,50);
    Muon.SetBSdxy(0.20,0.20);
    Muon.SetBSdz(0.50);
    Muon.MuonSelection(*muon_isPF, *muon_isGlobal, *muon_pt, *muon_eta, *muon_phi, *muon_energy, *muon_relIso04, *muon_q, *muon_validhits, *muon_validpixhits, *muon_matchedstations, *muon_trackerlayers, *muon_normchi, *muon_dxy, *muon_dz, muonLooseNotTightColl);
    
    Electron.SetPt(15);
    Electron.SetEta(2.5);
    Electron.SetRelIso(0.15);
    Electron.SetBSdxy(0.02);
    Electron.SetBSdz(0.10);
    Electron.ElectronSelection(*electrons_scEta, *electrons_pt, *electrons_eta, *electrons_phi, *electrons_energy, *electrons_phIso03, *electrons_nhIso03, *electrons_chIso03, *electrons_puChIso03, *electrons_q, *electrons_passConversionVeto, *electrons_electronID_snu, *electrons_dxy, *electrons_dz, electronColl);
    
    Jets.SetPt(30);
    Jets.SetEta(2.4);
    Jets.JetSelectionLeptonVeto(*jets_isTight, *jets_pt, *jets_eta, *jets_phi, *jets_energy, *jets_CSVInclV2, electronColl, muonColl, jetColl);

    if(debug) cout<< "DONE object selection" <<endl;

    if(muonColl.size()>0)
      for (Int_t i=0; i<muonColl.size(); i++) {
	index=muonColl[i].ilepton();
	h_muons[0][0]->Fill(weight, (Int_t) muonColl.size(), muonColl[i].lorentzVec(), muonColl[i].charge(), muonColl[i].relIso(), muonColl[i].chiNdof(), muonColl[i].dxy_BS(), muonColl[i].dz_BS());
      }
    if (electronColl.size() > 0) {
      for (UInt_t i=0; i<electronColl.size(); i++) {
	index=electronColl[i].ilepton();
	h_electrons[0][0]->Fill(weight, (Int_t) electronColl.size(), electronColl[i].lorentzVec(), electronColl[i].charge(), electronColl[i].relIso(), electronColl[i].dxy_BS(), electronColl[i].dz_BS());
      }
    }
    if(jetColl.size()>0)
      for (int i=0; i<jetColl.size(); i++) {
	index=jetColl[i].ijet();
	h_jets[0][0]->Fill(weight, (Int_t) jetColl.size(), jetColl[i].lorentzVec(), jets_CSVInclV2->at(index), jets_vtx3DSig->at(index) );
      }

    MET = metNoHF_pt->at(0);
    MET_phi =  metNoHF_phi->at(0);
    
    if(debug) cout<< "generic plots FILLED" <<endl;

      if (muonColl.size()>=4 && jetColl.size()>=2)
	h_prova->Fill((muonColl[0].lorentzVec()+muonColl[1].lorentzVec()+muonColl[2].lorentzVec()+muonColl[3].lorentzVec()+jetColl[0].lorentzVec()+jetColl[1].lorentzVec()).M(),weight);
      else if (electronColl.size()>=4 && jetColl.size()>=2)
	h_prova->Fill((electronColl[0].lorentzVec()+electronColl[1].lorentzVec()+electronColl[2].lorentzVec()+electronColl[3].lorentzVec()+jetColl[0].lorentzVec()+jetColl[1].lorentzVec()).M(),weight);
      else if (muonColl.size()>=2 && electronColl.size()>=2 && jetColl.size()>=2)
	h_prova->Fill((muonColl[0].lorentzVec()+muonColl[1].lorentzVec()+electronColl[0].lorentzVec()+electronColl[1].lorentzVec()+jetColl[0].lorentzVec()+jetColl[1].lorentzVec()).M(),weight);
      
    
    if (muonColl.size()>=2)
      selectionStep.push_back(1);
    if (muonColl.size()>=2 && jetColl.size()>=2)
      selectionStep.push_back(2);	
    if (electronColl.size()==0 && muonLooseColl.size()==2) {
      if (muonLooseColl[0].charge()==muonLooseColl[1].charge() && (muonLooseColl[0].lorentzVec()+muonLooseColl[1].lorentzVec()).M()>10) {
	selectionStep.push_back(3);
	if (jetColl.size()==1)
	  selectionStep.push_back(4);
	if (jetColl.size()>1) {
	  selectionStep.push_back(5);
	}
      }
    }
    for (UInt_t m=0;m<selectionStep.size();m++) {
      cut = selectionStep[m];
      channel = 0;
      UInt_t dataType = 0;
      UInt_t lep0 = 0;
      UInt_t lep1 = 1;
      if (muonColl.size()==2) {
	h_signal[cut][channel]->Fill(nGoodPV, MET, MET_phi, muonColl, jetColl, weight, channel, cut);
	for (Int_t i=0; i<muonColl.size(); i++) 
	  h_muons[cut][channel]->Fill(weight, (Int_t) muonColl.size(), muonColl[i].lorentzVec(), muonColl[i].charge(), muonColl[i].relIso(), muonColl[i].chiNdof(), muonColl[i].dxy_BS(), muonColl[i].dz_BS());
	for (UInt_t i=0; i<electronColl.size(); i++) 
	  h_electrons[cut][channel]->Fill(weight, (Int_t) electronColl.size(), electronColl[i].lorentzVec(), electronColl[i].charge(), electronColl[i].relIso(), electronColl[i].dxy_BS(), electronColl[i].dz_BS());   
	for (int i=0; i<jetColl.size(); i++) {
	  index=jetColl[i].ijet();
	  h_jets[cut][channel]->Fill(weight, (Int_t) jetColl.size(), jetColl[i].lorentzVec(), jets_CSVInclV2->at(index), jets_vtx3DSig->at(index) );
	}
      }
      else if (muonColl.size()==1) {
	if (muonColl[0].lorentzVec().Pt() > muonLooseNotTightColl[0].lorentzVec().Pt()) {
	  muonSelected.push_back(muonColl[0]);
	  muonSelected.push_back(muonLooseNotTightColl[0]);
	}
	else {
	  muonSelected.push_back(muonLooseNotTightColl[0]);
	  muonSelected.push_back(muonColl[0]);
	}	
	SingleFake=SinglebackGround(FRhisto, muonLooseNotTightColl, lep0, singleFake, dataType, weight);
	DoubleANDSinglebkg(muonColl, lep0, muonLooseNotTightColl, lep1, doubleANDsingleFake, dataType);
	h_singlefakes[cut][channel]->Fill(nGoodPV, MET, MET_phi, muonLooseColl, jetColl, SingleFake*weight, channel, cut);
	h_totalfakes[cut][channel]->Fill(nGoodPV, MET, MET_phi, muonLooseColl, jetColl, SingleFake*weight, channel, cut);
	muonSelected.clear();
      }
      else {
	DoubleFake=DoublebackGround(FRhisto, muonLooseColl, lep0, lep1, doubleFake, dataType, weight);
	Single_Double=DoubleTOSinglebkg(FRhisto, muonLooseColl, lep0, lep1);
	h_doublefakes[cut][channel]->Fill(nGoodPV, MET, MET_phi, muonLooseColl, jetColl, DoubleFake*weight, channel, cut);
	h_totalfakes[cut][channel]->Fill(nGoodPV, MET, MET_phi, muonLooseColl, jetColl, (DoubleFake+Single_Double)*weight, channel, cut);
      }
    }
    ///Filling standard particle plots END

  }
  if(debug) cout<< "out of the loop" <<endl;
  cout << "writing histos" << endl;
  outfile->cd();
  h_prova->Write();
  Dir = outfile->mkdir("Signal");
  Dir = outfile->mkdir("SingleFakes");
  Dir = outfile->mkdir("DoubleFakes");
  Dir = outfile->mkdir("TotalFakes");
  Dir = outfile->mkdir("Muons");
  Dir = outfile->mkdir("Electrons");
  Dir = outfile->mkdir("Jets");

  for(UInt_t i=0;i<ncuts;i++)
    for(UInt_t j=0;j<nchannels;j++){
    outfile->cd( "Signal" );
    h_signal[i][j]->Write();
    outfile->cd( "SingleFakes" );
    h_singlefakes[i][j]->Write();
    outfile->cd( "DoubleFakes" );
    h_doublefakes[i][j]->Write();
    outfile->cd( "TotalFakes" );
    h_totalfakes[i][j]->Write();
    outfile->cd( "Muons" );
    h_muons[i][j]->Write();
    outfile->cd( "Electrons" );
    h_electrons[i][j]->Write();
    outfile->cd( "Jets" );
    h_jets[i][j]->Write();
    }
  outfile->cd();
  outfile->Close();
  cout<<"histo written."<<endl;
  
}










void Analyzer::LoopFR() {

  cout << "total number of entries " <<nentries<<endl;

  if (debug) cout<< "loop begins" <<endl;

  fBTagSF = new BTagSFUtil("CSVM");

  // once we have data we must look a the pileup
  //  reweightPU = new ReweightPU("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/MyDataPileupHistogram_69400.root");

  if (debug) cout<< "PU histos loaded" <<endl;



  Double_t arrayeta [] = {0.0,0.8,1.479,2.0,2.5};
  Double_t arraypT [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  h_MT = new TH1F ("h_MT", "Transverse mass", 50, 0.0, 250.0);
  h_dileptMass = new TH1F ("h_dileptMass", "dilepton mass", 100, 0.0, 500.0);
  h_MET = new TH1F ("h_MET", "Missing energy", 100, 0.0, 500.0);  
  h_METPhi = new TH1F ("h_METPhi", "MET #phi", 100,-3.1415926535,3.1415926535);
  h_HT = new TH1F ("h_HT", "Sum p_{T} of the jets", 100,  0.0, 500.0);
  h_dr = new TH1F ("h_dr", "#Delta R between muon and jet", 100, 0.0, 10.0);
  h_dPhi = new TH1F ("h_dPhi", "#Delta #phi between muon and jet", 100, 0.0, 3.15);;
  h_ptMuOverptJet = new TH1F ("h_ptMuOverptJet", "#mu p_{T} over jet p_{T}", 100, 0.0, 10.0);
  h_nEvents = new TH2F ("h_nEvents", "Number of Events",ninteta,arrayeta,nintpT,arraypT);
  h_nEventsFO = new TH2F ("h_nEventsFO", "Number of Events FO",ninteta,arrayeta,nintpT,arraypT);
  h_FOrate = new TH2F ("h_FOrate", "FO rate",ninteta,arrayeta,nintpT,arraypT);

  h_TLnum = new MuonPlots("TL_numerator");
  h_TLden = new MuonPlots("TL_denominator");
  h_muonsFR = new MuonPlots("tight_muons");
  h_muonsLooseFR = new MuonPlots("loose_muons");
  h_electronsFR = new ElectronPlots("electrons");
  h_jetsFR = new JetPlots("jets");

  if(!MCweight) MCweight=1; 

  if (fChain == 0) 
    cout << "Ciao!" << endl;

  if (entrieslimit != -1)
    nentries=entrieslimit;

  if (debug) cout<< "at the loop" <<endl;
  std::set<int> runs;
  for (Long64_t jentry = 0; jentry < nentries; jentry++ ) {
    //clearing vectors
    selectionStep.clear();    selectChannel.clear();

    
    if (debug) cout<< "Event number " <<jentry<<endl;
    if (debug) cout<<"begin loop"<<endl;
    if (!(jentry % 50000)) 
      cout << jentry << endl;

    if (!fChain) cout<<"problems with the input file"<<endl;
    fChain->GetEntry(jentry);

    // Vertex Select
    if ( !goodVertices ) continue;
    // GOLD JSON
    if (lumiMaskGold<1) continue;
    
    triggerOK = false;
    for(UInt_t t=0; t<vtrignames->size(); t++) {
      trigger = vtrignames->at(t);
      Int_t ps = vtrigps->at(t);
      if ( (trigger.Contains("HLT_Mu8_TrkIsoVVL_v") || trigger.Contains("HLT_Mu17_TrkIsoVVL_v")) && ps>0) {
        triggerOK = true;
        break;
      }
    }
    if(!triggerOK) continue;

    // MET filters for now all OFF
    //----------------------------------------------------------------------------
    if (!(HBHENoiseFilter && CSCTightHaloFilter && eeBadScFilter && EcalDeadCellTriggerPrimitiveFilter)) continue;
    
    weight=MCweight;
    //MC@NLO weight
    if(MCatNLO)
      genWeight>=0 ? weight*=1. : weight*=-1.;
    
    h_VertexNoReweight->Fill(nGoodPV,weight);
    if (!IsData)
      weight*=puWeightGold;
    h_VertexPostReweight->Fill(nGoodPV,weight);

    if(debug) cout<< "object selection" <<endl;
    
    std::vector<Lepton> muonColl;
    Muon.SetPt(15);
    Muon.SetEta(2.4);
    Muon.SetRelIso(0.15);
    Muon.SetChiNdof(10);
    Muon.SetBSdxy(0.20);
    Muon.SetBSdz(0.50);
    Muon.MuonSelection(*muon_isPF, *muon_isGlobal, *muon_pt, *muon_eta, *muon_phi, *muon_energy, *muon_relIso04, *muon_q, *muon_validhits, *muon_validpixhits, *muon_matchedstations, *muon_trackerlayers, *muon_normchi, *muon_dxy, *muon_dz, muonColl);

    std::vector<Lepton> muonLooseColl;
    Muon.SetPt(15);
    Muon.SetEta(2.4);
    Muon.SetRelIso(0.6);
    Muon.SetChiNdof(50);
    Muon.SetBSdxy(0.20);
    Muon.SetBSdz(0.50);
    Muon.MuonSelection(*muon_isPF, *muon_isGlobal, *muon_pt, *muon_eta, *muon_phi, *muon_energy, *muon_relIso04, *muon_q, *muon_validhits, *muon_validpixhits, *muon_matchedstations, *muon_trackerlayers, *muon_normchi, *muon_dxy, *muon_dz, muonLooseColl);
    
    std::vector<Lepton> electronColl;
    Electron.SetPt(15);
    Electron.SetEta(2.5);
    Electron.SetRelIso(0.15);
    Electron.SetBSdxy(0.02);
    Electron.SetBSdz(0.10);
    Electron.ElectronSelection(*electrons_scEta, *electrons_pt, *electrons_eta, *electrons_phi, *electrons_energy, *electrons_phIso03, *electrons_nhIso03, *electrons_chIso03, *electrons_puChIso03, *electrons_q, *electrons_passConversionVeto, *electrons_electronID_snu, *electrons_dxy, *electrons_dz, electronColl);
    
    std::vector<Jet> jetColl;
    Jets.SetPt(40);
    Jets.SetEta(2.4);
    Jets.JetSelectionLeptonVeto(*jets_isTight, *jets_pt, *jets_eta, *jets_phi, *jets_energy, *jets_CSVInclV2, electronColl, muonColl, jetColl);

    if(debug) cout<< "DONE object selection" <<endl;

    if(muonColl.size()>0)
      for (Int_t i=0; i<muonColl.size(); i++)
	h_muonsFR->Fill(weight, (Int_t) muonColl.size(), muonColl[i].lorentzVec(), muonColl[i].charge(), muonColl[i].relIso(), muonColl[i].chiNdof(), muonColl[i].dxy_BS(), muonColl[i].dz_BS());
    if(muonLooseColl.size()>0)
      for (Int_t i=0; i<muonLooseColl.size(); i++) 
	h_muonsLooseFR->Fill(weight, (Int_t) muonLooseColl.size(), muonLooseColl[i].lorentzVec(), muonLooseColl[i].charge(), muonLooseColl[i].relIso(), muonLooseColl[i].chiNdof(), muonLooseColl[i].dxy_BS(), muonLooseColl[i].dz_BS());
    if (electronColl.size() > 0) 
      for (UInt_t i=0; i<electronColl.size(); i++) 
	h_electronsFR->Fill(weight, (Int_t) electronColl.size(), electronColl[i].lorentzVec(), electronColl[i].charge(), electronColl[i].relIso(), electronColl[i].dxy_BS(), electronColl[i].dz_BS());
    if(jetColl.size()>0)
      for (int i=0; i<jetColl.size(); i++) {
	index=jetColl[i].ijet();
	h_jetsFR->Fill(weight, (Int_t) jetColl.size(), jetColl[i].lorentzVec(), jets_CSVInclV2->at(index), jets_vtx3DSig->at(index) );
      }
    
    if(debug) cout<< "setting MET" <<endl;

    MET = met_pt->at(0);
    MET_phi =  met_phi->at(0);
    
    if(debug) cout<< "generic plots FILLED" <<endl;

    if (muonColl.size()==2)
      h_dileptMass->Fill( (muonColl[0].lorentzVec()+muonColl[1].lorentzVec()).M(), weight );
    if (muonColl.size()==1) {
      h_MT->Fill( sqrt(2.*muonColl[0].lorentzVec().Pt()*MET* (1 - cos(muonColl[0].lorentzVec().Phi()-MET_phi)) ), weight);
      h_MET->Fill(MET, weight);
      h_METPhi->Fill(MET_phi, weight);
    }

    if ( ZandWveto(muonLooseColl, MET, MET_phi) ) continue;

    if (muonLooseColl.size() != 1) continue;
    if (electronColl.size() > 0) continue;


    ///// FAKE RATES /////

    Double_t dr=-999.9;
    Double_t dPhi=-999.9;
    Double_t ptMuOverptJet=-999.9;

    if (debug) cout<< "denominator" <<endl;

    if (muonLooseColl.size() == 1 && jetColl.size() >= 1) {
      for (UInt_t iii=0; iii<muonLooseColl.size(); iii++) {
        for (UInt_t yyy=0; yyy<jetColl.size(); yyy++) {
          dr = muonLooseColl[iii].lorentzVec().DeltaR( jetColl[yyy].lorentzVec() );
          dPhi = fabs(muonLooseColl[iii].lorentzVec().DeltaPhi(jetColl[yyy].lorentzVec()));
          ptMuOverptJet = muonLooseColl[iii].lorentzVec().Pt()/jetColl[yyy].lorentzVec().Pt();
          //if (dPhi > 2.5 && ptMuOverptJet < 1.0)
	  h_dr->Fill(dr, weight);
          if (dr > 1.0 && ptMuOverptJet < 1.0)
            h_dPhi->Fill(dPhi,weight);
          if (dr > 1.0)
            h_ptMuOverptJet->Fill(ptMuOverptJet, weight);
          if (dr > 1.0 && ptMuOverptJet < 1.0 && dPhi>2.5) {
            //h_dPhi->Fill(dPhi,weight);
            h_nEventsFO->Fill(fabs(muonLooseColl[iii].eta()),muonLooseColl[iii].lorentzVec().Pt(), weight);
            h_TLden->Fill(weight, (Int_t) muonLooseColl.size(), muonLooseColl[iii].lorentzVec(), muonLooseColl[iii].charge(), muonLooseColl[iii].relIso(), muonLooseColl[iii].chiNdof(), muonLooseColl[iii].dxy_BS(), muonLooseColl[iii].dz_BS());
            if (muonColl.size() == 1) {
              h_nEvents->Fill(fabs(muonColl[iii].eta()),muonColl[iii].lorentzVec().Pt(), weight);
              h_TLnum->Fill(weight, (Int_t) muonColl.size(), muonColl[iii].lorentzVec(), muonColl[iii].charge(), muonColl[iii].relIso(), muonColl[iii].chiNdof(), muonColl[iii].dxy_BS(), muonColl[iii].dz_BS());
            }
            goto fineFO;
          }
        }
      }
    }
  fineFO:
    ;

    ///Filling standard particle plots END

  }
  if(debug) cout<< "out of the loop" <<endl;
  getFakerate(h_nEvents,h_nEventsFO,h_FOrate,ninteta,nintpT);
  //  getFakerate(h_nVertex1,h_nVertex0,h_nVertex2,50);

  outfile->cd();
  h_HT->Write();
  h_MT->Write();
  h_MET->Write();
  h_dileptMass->Write();
  h_METPhi->Write();
  h_dr->Write();
  h_dPhi->Write();
  h_ptMuOverptJet->Write();
  h_nEvents->Write();
  h_nEventsFO->Write();
  h_FOrate->Write();
  h_VertexNoReweight->Write();
  h_VertexPostReweight->Write();

  Dir = outfile->mkdir("Muons");
  outfile->cd( Dir->GetName() );
  h_muonsFR->Write();
  h_muonsLooseFR->Write();
  h_TLden->Write();
  h_TLnum->Write();
  outfile->cd();

  Dir = outfile->mkdir("Electrons");
  outfile->cd( Dir->GetName() );
  h_electronsFR->Write();
  outfile->cd();

  Dir = outfile->mkdir("Jets_with_veto");
  outfile->cd( Dir->GetName() );
  h_jetsFR->Write();
  outfile->cd();

  outfile->Close();

  cout<<"histo written."<<endl;
  

}
