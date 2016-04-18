#include <vector>

#include "SinglePlot.h"
#include "Single2dPlot.h"

void loadCFO(std::vector<TString>& filename, std::vector<TString>& legendname, std::vector<TString>& plotlabel, std::vector<int>& color, std::vector<int>& linecol, std::vector<std::string>& type, std::vector<SinglePlot>& hist1d, std::vector<Single2dPlot>& hist2d, std::vector<double>& weight, std::vector<bool>& legend) {

  const double luminosity = 1.0;
  //  const double luminosity = 0.004;
  TString ver = "9";

  const TString directory = "/uscms/home/byates/CMSSW_7_6_4/src/Majorana_Neutrino/test/histo/";
  std::vector<TString> classe;
  enum logbool {      nolog,         log };
  enum normbool {     nonorm,        norm };
  enum normfirstbool {nonormToFirst, normToFirst };
  enum stackbool {    nostack,       stack };
  enum overflowbool { nooverflow,    overflow };

  Bool_t ttbar=true; Bool_t ZDY_jets=false; Bool_t W_jets=true; Bool_t QCD=false;
  Bool_t DY=true; Bool_t tW = false; Bool_t ttW=true; Bool_t ttZ=true;
  Bool_t WW=true; Bool_t WZ=true; Bool_t ZZ=true;  Bool_t WpWp=true; 
  Bool_t data=true;
  Bool_t signal=true;
  Bool_t FR=false;
  Bool_t QFlip = true;
  TString cut[] = {"two_mu","probe","no_jets","low_MET","SS","OS"};
  TString channel[] = {"mumu"};


  if(QFlip) {
    for (UInt_t i=0; i < sizeof(cut)/sizeof(TString); i++)
      for(UInt_t j=0; j < sizeof(channel)/sizeof(TString); j++) {
        cout << "starting plots" << endl;
        cout << cut[i] + " " + channel[j] << endl;
/*
  	hist1d.push_back( SinglePlot("Muons/h_N_muons_"+cut[i]+"_"+channel[j], 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "N","","") );
  	hist1d.push_back( SinglePlot("Muons/h_charge_muons_"+cut[i]+"_"+channel[j], 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Charge","muon charge","") );
	hist1d.push_back( SinglePlot("Muons/h_pt_muons_"+cut[i]+"_"+channel[j], 1,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "Pt","muon p_{T} (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Muons/h_eta_muons_"+cut[i]+"_"+channel[j], 1,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "eta","muon #eta","Events") );
	hist1d.push_back( SinglePlot("Muons/h_phi_muons_"+cut[i]+"_"+channel[j], 1,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","muon #phi") );
	hist1d.push_back( SinglePlot("Muons/h_charge_muons_"+cut[i]+"_"+channel[j], 1,  nolog, nonorm, nonormToFirst, 2.0, nooverflow, stack, "","","muon charge") );
	hist1d.push_back( SinglePlot("Muons/h_GlbChi2_muons_"+cut[i]+"_"+channel[j], 1,  nolog, nonorm, nonormToFirst, 50.0, nooverflow, stack, "","","muon global #Chi^{2}") );
	hist1d.push_back( SinglePlot("Muons/h_PF_RelIso_muons_"+cut[i]+"_"+channel[j], 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","muon relIso") );
	hist1d.push_back( SinglePlot("Muons/h_dxy_muons_"+cut[i]+"_"+channel[j], 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","muon d_{xy}") );
	hist1d.push_back( SinglePlot("Muons/h_dz_muons_"+cut[i]+"_"+channel[j], 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","muon d_{z}") );
*/
        hist1d.push_back( SinglePlot("Muons/h_charge_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "charge","charge","") );
        hist1d.push_back( SinglePlot("Muons/h_mass_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "m(ll)","m(ll) (GeV)","") );
        hist1d.push_back( SinglePlot("Muons/h_MET_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "MET","#slash{E}_{T} (GeV)","") );
        hist1d.push_back( SinglePlot("Muons/h_nvtx_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "nvtx","nvtx","") );
        hist1d.push_back( SinglePlot("Muons/h_probePt_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Probe P_{T}","P_{T} (GeV)","") );
        hist1d.push_back( SinglePlot("Muons/h_tagPt_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Tag P_{T}","P_{T} (GeV)","") );
        hist1d.push_back( SinglePlot("Muons/h_invPt_"+cut[i]+"_"+channel[j], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "1/P_{T}","1/P_{T} (1/Gev)","") );
      }
  }
  
  
  if (signal && ttZ) {
    //filename.push_back(directory+"TTZ_"+ver+".root");
    filename.push_back(directory+"TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_"+ver+".root");
    legendname.push_back("t#bar{t}Z");
    plotlabel.push_back("t#bar{t}Z");
    color.push_back(kOrange); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
   // weight.push_back(0.70*1000./249275.);
    weight.push_back(luminosity);
  }
  
  if (signal && ttW) {
    //filename.push_back(directory+"TTWJets_"+ver+".root");
    filename.push_back(directory+"TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_"+ver+".root");
    legendname.push_back("t#bar{t}W"); plotlabel.push_back("t#bar{t}W");
    color.push_back(kOrange+7); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.72*1000./246521.);
    weight.push_back(luminosity);
  }

  if (signal && WpWp) {
    filename.push_back(directory+"WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8_"+ver+".root");
    legendname.push_back("W^{#pm}W^{#pm}");
    plotlabel.push_back("W^{#pm}W^{#pm}");
    color.push_back(kGreen+3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(12.5*1000./899900.);
    weight.push_back(luminosity);
  }
    
  if (signal && ZZ) {
    filename.push_back(directory+"ZZ_TuneCUETP8M1_13TeV-pythia8_"+ver+".root");
    legendname.push_back("ZZ");
    plotlabel.push_back("ZZ");
    color.push_back(kSpring+10); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(2.2*1000./237484.);
    weight.push_back(luminosity);
  }
  
  if (signal && WZ) {
    filename.push_back(directory+"WZ_TuneCUETP8M1_13TeV-pythia8_"+ver+".root");
    //filename.push_back(directory+"WZTo3LNu_14.root");
    legendname.push_back("WZ");
    plotlabel.push_back("WZ");
    color.push_back(kSpring); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(2.2*1000./237484.);
    weight.push_back(luminosity);
  }
  
  if (signal && WW) {
    filename.push_back(directory+"WW_TuneCUETP8M1_13TeV-pythia8_"+ver+".root");
    //filename.push_back(directory+"WWTo2L2Nu_Powheg_14.root");
    legendname.push_back("WW");
    plotlabel.push_back("WW");
    color.push_back(kGreen-3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(12.5*1000./899900.);
    weight.push_back(luminosity);
  }

  if (signal && DY) {
    filename.push_back(directory+"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_"+ver+".root");
    //filename.push_back(directory+"WWTo2L2Nu_Powheg_14.root");
    legendname.push_back("DY_10-50");
    plotlabel.push_back("DY_10-50");
    color.push_back(kOrange-3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(12.5*1000./899900.);
    weight.push_back(luminosity);
  }

  if (signal && DY) {
    filename.push_back(directory+"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_"+ver+".root");
    //filename.push_back(directory+"WWTo2L2Nu_Powheg_14.root");
    legendname.push_back("DY_50");
    plotlabel.push_back("DY_50");
    color.push_back(kOrange); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(12.5*1000./899900.);
    weight.push_back(luminosity);
  }

  if (signal && tW) {
    filename.push_back(directory+"ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_"+ver+".root");
    legendname.push_back("tW");
    plotlabel.push_back("tW");
    color.push_back(kPink-7); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.03610181523);
    weight.push_back(luminosity);
  }

  if (signal && tW) {
    filename.push_back(directory+"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_"+ver+".root");
    legendname.push_back("#bar{t}W");
    plotlabel.push_back("#bar{t}W");
    color.push_back(kPink-9); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.03663305207);
    weight.push_back(luminosity);
  }
  
  if (signal && ttbar) {
    filename.push_back(directory+"TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_"+ver+".root");
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
    filename.push_back(directory+"ZDY_"+ver+".root");
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
    //filename.push_back(directory+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_"+ver+".root");
    filename.push_back(directory+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_10.root");
    //filename.push_back(directory+"WJets_Madgraph_14.root");
    legendname.push_back("Wjets");
    plotlabel.push_back("Wjets");
    color.push_back(kBlue-3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(6.141944936);
    weight.push_back(luminosity);
  }
  
  if (signal && QCD) {
    filename.push_back(directory+"QCD_"+ver+".root");
    legendname.push_back("QCD");
    plotlabel.push_back("QCD");
    color.push_back(kOrange); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.3965);
    weight.push_back(luminosity);
  }

  if (signal && data) {
    //filename.push_back(directory+"SingleMuon_"+ver+".root");
    filename.push_back(directory+"SingleMuon_9.root");
    legendname.push_back("Data");
    plotlabel.push_back("Data");
    color.push_back(kBlack); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("data");
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

