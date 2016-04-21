#!/bin/bash
dir="../DataSetList/"

#./runAnalyzer.exe -i WZ_TuneCUETP8M1_13TeV-pythia8.txt -o WZ_TuneCUETP8M1_13TeV-pythia8 -d $dir -v $1
#./runAnalyzer.exe -i ZZ_TuneCUETP8M1_13TeV-pythia8.txt -o ZZ_TuneCUETP8M1_13TeV-pythia8 -d $dir -v $1
#./runAnalyzer.exe -i WW_TuneCUETP8M1_13TeV-pythia8.txt -o WW_TuneCUETP8M1_13TeV-pythia8 -d $dir -v $1

./runAnalyzer.exe -i DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt -o DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 -d $dir -v $1
./runAnalyzer.exe -i DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt -o DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 -d $dir -v $1
./runAnalyzer.exe -i WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt -o WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 -d $dir -v $1 

./runAnalyzer.exe -i TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.txt -o TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 -d $dir -v $1


./runAnalyzer.exe -i SingleMuon.txt -o SingleMuon -d $dir -v $1
