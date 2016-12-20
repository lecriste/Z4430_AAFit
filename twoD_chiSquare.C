TH2F* twoD_chiSquare(const TH2* histo, const TH2* histoErr, const RooAbsPdf* fit, const RooRealVar* x, const RooRealVar* y) {
  
  TH2F* pullHist = (TH2F*)histo->Clone(TString::Format("Pull_%s",fit->GetName())); pullHist->SetTitle(TString::Format("Pull of %s",fit->GetTitle()));
  
  /*
  Float_t xMin = pullHist->GetXaxis()->GetXmin(), xMax = pullHist->GetXaxis()->GetXmax(); Int_t xBins = pullHist->GetNbinsX();
  Float_t yMin = pullHist->GetYaxis()->GetXmin(), yMax = pullHist->GetYaxis()->GetXmax(); Int_t yBins = pullHist->GetNbinsY();
  RooBinning* xRooBinning = new RooBinning(xBins,xMin,xMax);
  RooBinning* yRooBinning = new RooBinning(yBins,yMin,yMax); 
  //
  const TArrayD* xArray = pullHist->GetXaxis()->GetXbins(); const Double_t* xBinning = 0;
  const TArrayD* yArray = pullHist->GetYaxis()->GetXbins(); const Double_t* yBinning = 0;
  if (xArray) {
    xBinning = xArray->GetArray(); xMin = xBinning[0]; xMax = xBinning[xBins];
    if (xBinning) { xRooBinning = new RooBinning(xBins,xBinning); }
  }
  if (yArray) {
    yBinning = yArray->GetArray(); yMin = yBinning[0]; yMax = yBinning[yBins];
    if (yBinning) { yRooBinning = new RooBinning(yBins,yBinning); }
  }

  TH2F* pdfHist = (TH2F*)fit->createHistogram("pdfHist", *x, Binning(*xRooBinning), YVar(*y, Binning(*yRooBinning)) ) ;
  pullHist->Add(pdfHist, -1);
  */

  TF2* tf2 = (TF2*)fit->asTF(RooArgList(*x,*y));
  // chi2 = ((hist-fit)^2) / fit
  // chi2 = ((hist-fit)^2) / hist_err
  pullHist->Add(tf2, -1, "I"); // (hist-fit)
  pullHist->Multiply(pullHist); // ((hist-fit)^2)
  //pullHist->Divide(tf2); // ((hist-fit)^2) / fit
  pullHist->Divide(histoErr); // ((hist-fit)^2) / hist_err

  return pullHist;

}
