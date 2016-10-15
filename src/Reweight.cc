#include "Reweight.h"

ReweightPU::ReweightPU(TString filenameData) {
  fileData_ = 0;
  fileData_ = new TFile(filenameData, "READ");
  fileMC_   = new TFile("MCPileUp.root", "READ"); 
  if (!fileData_) {
    std::cout << "\n\n Data file of the Nvtx reweighting could not be opened!" << std::endl;
  }
  
  h_Data_ = 0;
  h_Data_ = (TH1D*)fileData_->Get("pileup");
  if (!h_Data_) std::cout << "Can't open PU data reweight histo";
  
  // Distribution used for Summer2012 MC.
  
  Double_t Summer2012_S10[60] = {
    2.560E-06,
    5.239E-06,
    1.420E-05,
    5.005E-05,
    1.001E-04,
    2.705E-04,
    1.999E-03,
    6.097E-03,
    1.046E-02,
    1.383E-02,
    1.685E-02,
    2.055E-02,
    2.572E-02,
    3.262E-02,
    4.121E-02,
    4.977E-02,
    5.539E-02,
    5.725E-02,
    5.607E-02,
    5.312E-02,
    5.008E-02,
    4.763E-02,
    4.558E-02,
    4.363E-02,
    4.159E-02,
    3.933E-02,
    3.681E-02,
    3.406E-02,
    3.116E-02,
    2.818E-02,
    2.519E-02,
    2.226E-02,
    1.946E-02,
    1.682E-02,
    1.437E-02,
    1.215E-02,
    1.016E-02,
    8.400E-03,
    6.873E-03,
    5.564E-03,
    4.457E-03,
    3.533E-03,
    2.772E-03,
    2.154E-03,
    1.656E-03,
    1.261E-03,
    9.513E-04,
    7.107E-04,
    5.259E-04,
    3.856E-04,
    2.801E-04,
    2.017E-04,
    1.439E-04,
    1.017E-04,
    7.126E-05,
    4.948E-05,
    3.405E-05,
    2.322E-05,
    1.570E-05,
    5.005E-06};

 Double_t Spring2016[60] = {
		0.000829312873542,
 		0.00124276120498,
 		0.00339329181587,
 		0.00408224735376,
 		0.00383036590008,
		0.00659159288946,
 		0.00816022734493,
 		0.00943640833116,
 		0.0137777376066,
 		0.017059392038,
 		0.0213193035468,
 		0.0247343174676,
 		0.0280848773878,
 		0.0323308476564,
 		0.0370394341409,
 		0.0456917721191,
 		0.0558762890594,
 		0.0576956187107,
 		0.0625325287017,
 		0.0591603758776,
 		0.0656650815128,
 		0.0678329011676,
 		0.0625142146389,
 		0.0548068448797,
 		0.0503893295063,
 		0.040209818868,
 		0.0374446988111,
 		0.0299661572042,
 		0.0272024759921,
 		0.0219328403791,
 		0.0179586571619,
 		0.0142926728247,
 		0.00839941654725,
 		0.00522366397213,
 		0.00224457976761,
 		0.000779274977993,
 		0.000197066585944,
 		7.16031761328e-05,
 		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
 		0.0,
 		0.0,
		0.0};
  
  h_MCmod_ = (TH1D*)h_Data_->Clone("h_MCmod_");
  for (Int_t i = 1; i < 61; i++) {
    h_MCmod_->SetBinContent(i, Spring2016[i-1]);
  }


  //h_MCmod_ = (TH1D*)fileMC_->Get("h_VertexNoReweight");
  double int_MC_ = h_MCmod_->Integral();
  double int_Data_ = h_Data_->Integral();

  h_Data_->Divide(h_MCmod_);
  h_Data_->Scale(int_MC_ / int_Data_);

  
    std::cout << std::endl;
    for (Int_t i = 1; i < 61; i++)
    std::cout << h_Data_->GetBinContent(i) <<std::endl;

}

ReweightPU::~ReweightPU() {
  delete fileData_;
  delete h_MCmod_;
  delete h_Data_;
}

Double_t ReweightPU::GetWeight(Int_t nvtx) {
  return h_Data_->GetBinContent( h_Data_->FindBin(nvtx) );
}
