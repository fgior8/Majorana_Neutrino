{
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  gROOT->ProcessLine(".L CMS_lumi.C+g");
  gStyle->SetPalette(1);
  gROOT->ProcessLine(".L HNmultiplot.C+g");
  multiplot();
  // gROOT->ProcessLine(".L multiplot_yields.C+g");
  // multiplot_yields();
}
