#include <vector>

#include "SinglePlot.h"
#include "Single2dPlot.h"

void loadCFO(std::vector<TString>& filename, std::vector<TString>& legendname, std::vector<TString>& plotlabel, std::vector<int>& color, std::vector<int>& linecol, std::vector<std::string>& type, std::vector<SinglePlot>& hist1d, std::vector<Single2dPlot>& hist2d, std::vector<double>& weight, std::vector<bool>& legend) {

  const double luminosity = 1.0;
  //  const double luminosity = 0.004;

  const TString directory = "/uscms/home/byates/CMSSW_7_4_12/src/Majorana_Neutrino/histo/";
  std::vector<TString> classe;
  enum logbool {      nolog,         log };
  enum normbool {     nonorm,        norm };
  enum normfirstbool {nonormToFirst, normToFirst };
  enum stackbool {    nostack,       stack };
  enum overflowbool { nooverflow,    overflow };

  Bool_t ttbar=false; Bool_t ZDY_jets=false; Bool_t W_jets=false; Bool_t QCD=false;
  Bool_t tW = false; Bool_t ttW=false; Bool_t ttZ=false;
  Bool_t WW=true; Bool_t WZ=true; Bool_t ZZ=true;  Bool_t WpWp=false; 
  Bool_t data=false;
  Bool_t signal=false;
  Bool_t FR=false;
  Bool_t QFlip = true;
  TString cut[] = {"two_mu","probe","no_jets","low_MET","SS","OS"};
  TString channel[] = {"mumu"};


  if(false) {
    classe.push_back("muons");

    for(UInt_t iii=0; iii<classe.size(); iii++) {
      cout<<classe[iii]<<endl;
      if (false) {
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_N", 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_pt", 1,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "muon p_{T} (GeV)","Events/10 GeV", "") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_eta", 2,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack,"","", "muon #eta") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_phi", 2,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack,"","", "") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_charge", 1,  nolog, nonorm, nonormToFirst, 2.0, nooverflow, stack,"","", "") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_GlbChi2", 1,  nolog, nonorm, nonormToFirst, 50.0, nooverflow, stack,"","", "") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_PF_RelIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack,"","", "") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_dxy", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack,"","", "") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_dz", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack,"","", "") );
      }
      if (false) {
	hist1d.push_back( SinglePlot("Electrons/h_N_"+classe[iii], 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Electrons/h_"+classe[iii]+"_Detector_RelIso_rho", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Electrons/h_"+classe[iii]+"_pt", 5,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "","muon p_{T} (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Electrons/h_"+classe[iii]+"_eta", 2,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "","muon #eta","Events") );
	hist1d.push_back( SinglePlot("Electrons/h_"+classe[iii]+"_phi", 2,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Electrons/h_"+classe[iii]+"_charge", 1,  nolog, nonorm, nonormToFirst, 2.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Electrons/h_"+classe[iii]+"_dxy", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Electrons/h_"+classe[iii]+"_dz", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","") );
      }
    }
  if (false) {
      hist1d.push_back( SinglePlot("h_VertexNoReweight", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("h_VertexPostReweight", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("h_dr", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("h_dPhi", 2, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("h_ptMuOverptJet", 1, log, nonorm, nonormToFirst, 5.0, nooverflow, stack, "","#mu(p_{T})/jet(p_{T})","") );
      hist1d.push_back( SinglePlot("h_HT", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("h_MET", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","#slash{E}_{T} (GeV)","Events/2 GeV") );
      hist1d.push_back( SinglePlot("h_METPhi", 1, nolog, nonorm, nonormToFirst,1.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("h_MT", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","M_{T} (GeV)","Events/2 GeV") );
      hist1d.push_back( SinglePlot("h_dileptMass", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","m(ll) (GeV)","Events/2 GeV") );
      hist1d.push_back( SinglePlot("Jets_with_veto/h_jets_N", 1,  log, nonorm, nonormToFirst, 30.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("Jets_with_veto/h_jets_pt", 1,  nolog, nonorm, nonormToFirst, 250.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("Jets_with_veto/h_jets_eta", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
      hist1d.push_back( SinglePlot("Jets_with_veto/h_jets_phi", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
    }
  }
  if(QFlip) {
    for (UInt_t i=0; i < sizeof(cut)/sizeof(TString); i++)
      for(UInt_t j=0; j < sizeof(channel)/sizeof(TString); j++) {
        cout << "starting plots" << endl;
/*
  	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_N", 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "N","","") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_pt", 1,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "Pt","muon p_{T} (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_eta", 1,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "eta","muon #eta","Events") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_phi", 1,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","muon #phi") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_charge", 1,  nolog, nonorm, nonormToFirst, 2.0, nooverflow, stack, "","","muon charge") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_GlbChi2", 1,  nolog, nonorm, nonormToFirst, 50.0, nooverflow, stack, "","","muon global #Chi^{2}") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_PF_RelIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","muon relIso") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_dxy", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","muon d_{xy}") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[i]+"_"+channel[j]+"_dz", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","muon d_{z}") );
*/
/*
        hist1d.push_back( SinglePlot("Muons/h_charge_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "charge","","") );
        hist1d.push_back( SinglePlot("Muons/h_mass_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "m(ll)","","") );
        hist1d.push_back( SinglePlot("Muons/h_MET_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "MET","","") );
        hist1d.push_back( SinglePlot("Muons/h_nvtx_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "nvtx","","") );
        hist1d.push_back( SinglePlot("Muons/h_probePt_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Probe P_{T}","","") );
        hist1d.push_back( SinglePlot("Muons/h_tagPt_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Tag P_{T}","","") );
        hist1d.push_back( SinglePlot("Muons/h_invPt_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "1/P_{T}","","") );
*/
        hist1d.push_back( SinglePlot("Muons/h_charge_OS_mumu", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "charge", "charge", "events") );
      }
  }
  
  if(signal) {
    //hist1d.push_back( SinglePlot("h_prova", 1,  log, nonorm, nonormToFirst, 1.0, overflow, stack, "", "m(ll)","Events/10 GeV") );
    for (UInt_t k=0;k<9;k++)
      for (UInt_t w=0;w<1;w++) {
      if (k!=0) {
	hist1d.push_back( SinglePlot("h_MET_signal_"+cut[k]+"_"+channel[w], 2,  log, nonorm, nonormToFirst, 1.0, overflow, stack, "", "Missing Energy (GeV)","Events/10 GeV") );	
	hist1d.push_back( SinglePlot("h_PV_signal_"+cut[k]+"_"+channel[w], 2,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Number of Vertices","Events") );
	hist1d.push_back( SinglePlot("h_MET_phi_signal_"+cut[k]+"_"+channel[w], 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Missing Energy #phi","Events") );
	hist1d.push_back( SinglePlot("h_MT2ll_signal_"+cut[k]+"_"+channel[w], 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "MT_{2}(ll) (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("h_llmass_signal_"+cut[k]+"_"+channel[w], 2,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "di-lepton invariant mass (GeV)","Events/10 GeV") );
	if (k>1 && k!=3 && k!=4) {
	  hist1d.push_back( SinglePlot("h_MT2bb_signal_"+cut[k]+"_"+channel[w], 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "MT_{2}(bb) (GeV)","Events/10 GeV") );
	  hist1d.push_back( SinglePlot("h_MT2lblb_signal_"+cut[k]+"_"+channel[w], 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "MT_{2}(lblb) (GeV)","Events/10 GeV") );
	  hist1d.push_back( SinglePlot("h_l1jjmass_signal_"+cut[k]+"_"+channel[w], 5,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading lepton plus jets invariant mass (GeV)","Events/10 GeV") );
	  hist1d.push_back( SinglePlot("h_l2jjmass_signal_"+cut[k]+"_"+channel[w], 5,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "trailing lepton plus jets invariant mass (GeV)","Events/10 GeV") );
	  hist1d.push_back( SinglePlot("h_dijetsmass_signal_"+cut[k]+"_"+channel[w], 5,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "di-jet invariant mass (GeV)","Events/10 GeV") );
	  hist1d.push_back( SinglePlot("h_lljjmass_signal_"+cut[k]+"_"+channel[w], 5,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "4 particles invariant mass (GeV)","Events/10 GeV") );
	}
      }
      if (((k==0 && w==0) || k>0) && false) {
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[k]+"_"+channel[w]+"_N", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Number of Muons","Events") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[k]+"_"+channel[w]+"_pt", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Muons p_{T}","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[k]+"_"+channel[w]+"_eta", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Muons #eta","Events") );
	hist1d.push_back( SinglePlot("Muons/h_muons_"+cut[k]+"_"+channel[w]+"_phi", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Muons #phi","Events") );
      
	hist1d.push_back( SinglePlot("Electrons/h_electrons_"+cut[k]+"_"+channel[w]+"_N", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Number of Electrons","Events") );
	hist1d.push_back( SinglePlot("Electrons/h_electrons_"+cut[k]+"_"+channel[w]+"_pt", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Electrons p_{T}","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Electrons/h_electrons_"+cut[k]+"_"+channel[w]+"_eta", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Electrons #eta","Events") );
	hist1d.push_back( SinglePlot("Electrons/h_electrons_"+cut[k]+"_"+channel[w]+"_phi", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Electrons #phi","Events") );
      
	hist1d.push_back( SinglePlot("Jets/h_jets_"+cut[k]+"_"+channel[w]+"_N", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Number of Jets","Events") );
	hist1d.push_back( SinglePlot("Jets/h_jets_"+cut[k]+"_"+channel[w]+"_pt", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Jets p_{T}","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Jets/h_jets_"+cut[k]+"_"+channel[w]+"_eta", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Jets #eta","Events") );
	hist1d.push_back( SinglePlot("Jets/h_jets_"+cut[k]+"_"+channel[w]+"_phi", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "Jets #phi","Events") );
      }
    }
    // hist2d.push_back( Single2dPlot("h_leadingJetPtvsMassLog_control", "h_leadingJetPtvsMassLog_control", "COLZ", 1, 1) );
    //hist2d.push_back( Single2dPlot("h_leadingJetPtvsMassLog_smoothing", "h_leadingJetPtvsMassLog_smoothing", "COLZ", 1, 1) );
  }
  
  if (signal && ttZ) {
    filename.push_back(directory+"TTZ_1.root");
    legendname.push_back("t#bar{t}Z");
    plotlabel.push_back("t#bar{t}Z");
    color.push_back(kOrange); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
   // weight.push_back(0.70*1000./249275.);
    weight.push_back(luminosity);
  }
  
  if (signal && ttW) {
    filename.push_back(directory+"TTWJets_1.root");
    legendname.push_back("t#bar{t}W");
    plotlabel.push_back("t#bar{t}W");
    color.push_back(kOrange+7); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.72*1000./246521.);
    weight.push_back(luminosity);
  }

  if (signal && WpWp) {
    filename.push_back(directory+"WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8_1.root");
    legendname.push_back("W^{#pm}W^{#pm}");
    plotlabel.push_back("W^{#pm}W^{#pm}");
    color.push_back(kGreen+3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(12.5*1000./899900.);
    weight.push_back(luminosity);
  }
    
  if (ZZ) {
    filename.push_back(directory+"ZZ_TuneCUETP8M1_13TeV-pythia8_1.root");
    legendname.push_back("ZZ");
    plotlabel.push_back("ZZ");
    color.push_back(kSpring+10); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(2.2*1000./237484.);
    weight.push_back(luminosity);
  }
  
  if (WZ) {
    filename.push_back(directory+"WZ_TuneCUETP8M1_13TeV-pythia8_1.root");
    //filename.push_back(directory+"WZTo3LNu_14.root");
    legendname.push_back("WZ");
    plotlabel.push_back("WZ");
    color.push_back(kSpring); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(2.2*1000./237484.);
    weight.push_back(luminosity);
  }
  
  if (WW) {
    filename.push_back(directory+"WW_TuneCUETP8M1_13TeV-pythia8_1.root");
    //filename.push_back(directory+"WWTo2L2Nu_Powheg_14.root");
    legendname.push_back("WW");
    plotlabel.push_back("WW");
    color.push_back(kGreen-3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(12.5*1000./899900.);
    weight.push_back(luminosity);
  }
  
  if (signal && tW) {
    filename.push_back(directory+"ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_1.root");
    legendname.push_back("tW");
    plotlabel.push_back("tW");
    color.push_back(kPink-7); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.03610181523);
    weight.push_back(luminosity);
  }

  if (signal && tW) {
    filename.push_back(directory+"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_1.root");
    legendname.push_back("#bar{t}W");
    plotlabel.push_back("#bar{t}W");
    color.push_back(kPink-9); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.03663305207);
    weight.push_back(luminosity);
  }
  
  if (signal && ttbar) {
    filename.push_back(directory+"TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1.root");
    //filename.push_back(directory+"TTJets_MadSpin_14.root ");
    legendname.push_back("t#bar{t}");
    plotlabel.push_back("t#bar{t}");
    color.push_back(kRed); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.03269);
    weight.push_back(luminosity);
  }
  
  if (signal && ZDY_jets) {
    filename.push_back(directory+"ZDY_1.root");
    //filename.push_back(directory+"ZJets_Madgraph_14.root ");
    legendname.push_back("Zjets");
    plotlabel.push_back("Zjets");
    color.push_back(kYellow); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(2.129675056);
    weight.push_back(luminosity);
  }
  
  if (signal && W_jets) {
    filename.push_back(directory+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1.root");
    //filename.push_back(directory+"WJets_Madgraph_14.root");
    legendname.push_back("Wjets");
    plotlabel.push_back("Wjets");
    color.push_back(kGreen); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(6.141944936);
    weight.push_back(luminosity);
  }
  
  if (signal && QCD) {
    filename.push_back(directory+"QCD_1.root");
    legendname.push_back("QCD");
    plotlabel.push_back("QCD");
    color.push_back(kOrange); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.3965);
    weight.push_back(luminosity);
  }

  if (signal && data) {
    filename.push_back(directory+"DoubleMuon_1.root");
    legendname.push_back("Data");
    plotlabel.push_back("Data");
    color.push_back(kBlack); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("data");
    weight.push_back(1.0);
  }

  if (signal && data) {
    filename.push_back("/Users/ferdi/cernbox/HNfiles/vuoto.root");
    legendname.push_back("Misid. muon bkgd.");
    plotlabel.push_back("Misid. muon bkgd.");
    color.push_back(kAzure+10); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    weight.push_back(1.0);
  }


  /////////////////////
  //   FR plots      //
  /////////////////////

  if (FR && ZZ) {
    filename.push_back(directory+"ZZ_TuneCUETP8M1_13TeV-pythia8_1.root");
    legendname.push_back("ZZ");
    plotlabel.push_back("ZZ");
    color.push_back(kSpring+10); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(2.2*1000./237484.);
    weight.push_back(luminosity);
  }
  
  if (FR && WZ) {
    filename.push_back(directory+"WZ_TuneCUETP8M1_13TeV-pythia8_1.root");
    //filename.push_back(directory+"WZTo3LNu_14.root");
    legendname.push_back("WZ");
    plotlabel.push_back("WZ");
    color.push_back(kSpring); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(2.2*1000./237484.);
    weight.push_back(luminosity);
  }
  
  if (FR && WW) {
    filename.push_back(directory+"WW_TuneCUETP8M1_13TeV-pythia8_1.root");
    legendname.push_back("WW");
    plotlabel.push_back("WW");
    color.push_back(kGreen-3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(12.5*1000./899900.);
    weight.push_back(luminosity);
  }
  
  if (FR && ttbar) {
    filename.push_back(directory+"TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1.root");
    //filename.push_back(directory+"TTJets_MadSpin_14.root ");
    legendname.push_back("t#bar{t}");
    plotlabel.push_back("t#bar{t}");
    color.push_back(kRed); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.03269);
    weight.push_back(luminosity);
  }
  
  if (FR && ZDY_jets) {
    filename.push_back(directory+"ZDY_1.root");
    //filename.push_back(directory+"ZJets_Madgraph_14.root ");
    legendname.push_back("Zjets");
    plotlabel.push_back("Zjets");
    color.push_back(kYellow); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(2.129675056);
    weight.push_back(luminosity);
  }
  
  if (FR && W_jets) {
    filename.push_back(directory+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1.root");
    //filename.push_back(directory+"WJets_Madgraph_14.root");
    legendname.push_back("Wjets");
    plotlabel.push_back("Wjets");
    color.push_back(kGreen); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(6.141944936);
    weight.push_back(luminosity);
  }

  if (FR && QCD) {
    filename.push_back(directory+"QCD_1.root");
    legendname.push_back("QCD");
    plotlabel.push_back("QCD");
    color.push_back(kOrange); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.3965);
    weight.push_back(luminosity);
  }
  
  if (FR && data) {
    filename.push_back(directory+"DoubleMuon_1.root");
    legendname.push_back("Data");
    plotlabel.push_back("Data");
    color.push_back(kBlack); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("data");
    weight.push_back(1.0);
  }
  
  /////////////////////
  //   Formatting    //
  /////////////////////

  // SinglePlot(std::string name, unsigned int rebin, bool log, bool normalize, bool normToFirst, double scaleXmax,
  //            bool overflowbin, bool stacked, TString title)

  // Single2dPlot(std::string name, std::string title, std::string drawOption, unsigned int rebinX, unsigned int rebinY)

} 

