#ifndef Reweight_h
#define Reweight_h

#include <stdio.h>
#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"

class ReweightPU {

 public:
  ReweightPU(TString filenameData);
  ~ReweightPU();
  Double_t GetWeight(Int_t nvtx);
  
 private:
  
  TFile* fileData_;
  TFile* fileMC_; 
  TH1D* h_MCmod_;
  TH1D* h_Data_;
      
};


#endif

