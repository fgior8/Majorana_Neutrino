#!/bin/bash
#dir="../DataSetList/"

./runAnalyzer.exe -i WZ_TuneCUETP8M1_13TeV-pythia8.txt -o WZ_TuneCUETP8M1_13TeV-pythia8 -v $1
./runAnalyzer.exe -i ZZ_TuneCUETP8M1_13TeV-pythia8.txt -o ZZ_TuneCUETP8M1_13TeV-pythia8 -v $1
./runAnalyzer.exe -i WW_TuneCUETP8M1_13TeV-pythia8.txt -o WW_TuneCUETP8M1_13TeV-pythia8 -v $1
./runAnalyzer.exe -i WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8.txt -o WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8 -v $1
./runAnalyzer.exe -i WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.txt -o WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8 -v $1
./runAnalyzer.exe -i WW_DoubleScattering_13TeV-pythia8.txt -o WW_DoubleScattering_13TeV-pythia8 -v $1

./runAnalyzer.exe -i TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.txt -o TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8 -v $1
./runAnalyzer.exe -i TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.txt -o TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8 -v $1
#./runAnalyzer.exe -i TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt -o TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8 -v $1
./runAnalyzer.exe -i TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.txt -o TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8 -v $1

#./runAnalyzer.exe -i MajoranaNeutrinoToMM_M-40_TuneZ2star_13TeV-alpgen.txt -o MajoranaNeutrinoToMM_M-40_TuneZ2star_13TeV-alpgen -v $1
#./runAnalyzer.exe -i MajoranaNeutrinoToMM_M-100_TuneZ2star_13TeV-alpgen.txt -o MajoranaNeutrinoToMM_M-100_TuneZ2star_13TeV-alpgen -v $1
#./runAnalyzer.exe -i MajoranaNeutrinoToMM_M-500_TuneZ2star_13TeV-alpgen.txt -o MajoranaNeutrinoToMM_M-500_TuneZ2star_13TeV-alpgen -v $1
#./runAnalyzer.exe -i MajoranaNeutrinoToMM_M-1500_TuneZ2star_13TeV-alpgen.txt -o MajoranaNeutrinoToMM_M-1500_TuneZ2star_13TeV-alpgen -v $1
