#include "XSlist.h"

float getXS(TString MCsample) {
  float eventXS=-1.;
  
  if (MCsample == "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8") 
    //eventXS=815.96/11339232.;
    eventXS=831.8/11344206.0;
  if (MCsample == "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1") 
    eventXS=35.6/1000000.;
  if (MCsample == "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1") 
    eventXS=35.6/995600.;
  if (MCsample == "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") 
    //eventXS=1./3210.8;
    eventXS=6024.2/26653508.0;
  if (MCsample == "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") 
    //eventXS=1./1219.8;
    eventXS=18610.0/29639808.0;
  if (MCsample == "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") 
    //eventXS=1./268.863;
    //cs = 61526.7
    eventXS=61526.7/23577660.0;
  if (MCsample == "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") 
    eventXS=0.2043/252908.;
  if (MCsample == "TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") 
    eventXS=0.4062/833964.;
  if (MCsample == "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8") 
    eventXS=0.2529/398000.;
  if (MCsample == "TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8") 
    eventXS=0.5297/749800.;
  
  if (MCsample == "WW_TuneCUETP8M1_13TeV-pythia8") 
    //eventXS=110.8/994416.;
    //cs = 113.8
    eventXS=113.8/993640.0;
  if (MCsample == "WZ_TuneCUETP8M1_13TeV-pythia8") 
    //eventXS=66.1/991232.;
    //cs = 47.1
    eventXS=47.1/978512.0;
  if (MCsample == "ZZ_TuneCUETP8M1_13TeV-pythia8") 
    //eventXS=15.4/996168.;
    //cs = 16.9
    eventXS=16.9/996944.0;
  
  if (MCsample == "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=558528000./10611979.;
  if (MCsample == "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=139803000./9842165.;
  if (MCsample == "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=19222500./5069469.;
  if (MCsample == "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=2758420./2926805.;
  if (MCsample == "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=469797./4026104.;
  if (MCsample == "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=117989./3942640.;
  if (MCsample == "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=7820.25/3910268.;
  if (MCsample == "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=645.528/1928421.;
  if (MCsample == "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=187.109/1983363.;
  if (MCsample == "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=32.3486/1982314.;
  if (MCsample == "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=10.4305/1981954.;
  
  if (eventXS == -1.) {
    cerr<<endl<<"NO CROSS-SECTION FOUND, PLEASE CHECK YOUR FILENAME OR LOOK INTO XSlis.cc"<<endl<<endl;
    return 0;
  }
  else
    return eventXS;
}
