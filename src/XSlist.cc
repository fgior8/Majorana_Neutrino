#include "XSlist.h"

float getXS(TString MCsample) {
  float eventXS=-1.;
  
  if (MCsample == "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8") 
    //eventXS=815.96/11339232.;
    eventXS=831.8/10215131.0;
  if (MCsample == "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1") 
    eventXS=35.6/1000000.;
  if (MCsample == "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1") 
    eventXS=35.6/995600.;
  if (MCsample == "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") 
    //eventXS=1./3210.8;
    eventXS=6024.2/81236728.0;
  if (MCsample == "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") 
    //eventXS=1./1219.8;
    eventXS=18610.0/22606898.0;
  if (MCsample == "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") 
    //eventXS=1./268.863;
    eventXS=61526.7/16521035.0;
  if (MCsample == "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") 
    eventXS=0.2043/129001.0;
  if (MCsample == "TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8") 
    eventXS=0.4062/429599.0;
  if (MCsample == "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8") 
    eventXS=0.2529/184990.0;
  if (MCsample == "TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8") 
    eventXS=0.5297/351398.0;
  
  if (MCsample == "WW_TuneCUETP8M1_13TeV-pythia8") 
    //eventXS=110.8/994416.;
    eventXS=113.8/988418.0;
  if (MCsample == "WZ_TuneCUETP8M1_13TeV-pythia8") 
    //eventXS=66.1/991232.;
    eventXS=47.1/1000000.0;
  if (MCsample == "ZZ_TuneCUETP8M1_13TeV-pythia8") 
    //eventXS=15.4/996168.;
    eventXS=16.9/985600.0;
  if (MCsample == "WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8")
    eventXS=0.02064/145800.0; 
  if (MCsample == "WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8")
    eventXS=0.01538/120000.0;
  if (MCsample == "WW_DoubleScattering_13TeV-pythia8")
    eventXS=1.640/844954.0;
 
  if (MCsample == "QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=2960198.4/31680404.0;
  if (MCsample == "QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=1652471.5/29938364.0;
  if (MCsample == "QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=437504.1/20378392.0;
  if (MCsample == "QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=106033.7/13730591.0;
  if (MCsample == "QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=25190.5/7971018.0;
  if (MCsample == "QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=8654.5/7910182.0;
  if (MCsample == "QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=797.4/7845620.0;
  if (MCsample == "QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=79.0/3841262.0;
  if (MCsample == "QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=25.1/3984898.0;
  if (MCsample == "QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=4.7/3666110.0;
  if (MCsample == "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8") 
    eventXS=1.6/3938782.0;
  
  if (eventXS == -1.) {
    cerr<<endl<<"NO CROSS-SECTION FOUND, PLEASE CHECK YOUR FILENAME OR LOOK INTO XSlis.cc"<<endl<<endl;
    return 0;
  }
  else
    return eventXS;
}
