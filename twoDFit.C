‰#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooFitResult.h"

#include "sqDalitzToMassesPdf.h"

const Int_t nVars = 2;

Int_t varOrder[nVars] = {0,0};

//Double_t varLimit[] = {0,0}; // {xMin/Max, yMin/Max}
Double_t varLimit[] = {0,0,0,0}; // {xMin/Max,yMin/Max,xMin/Max,yMin/Max}

const Int_t nThrPar = sizeof(varLimit)/sizeof(varLimit[0]) ;


// with RooFit
RooAbsPdf* twoDFit(RooAbsReal& x, RooAbsReal& y, RooDataHist* hist) {

  RooArgList xCoeff("xCoeff");
  RooArgList yCoeff("yCoeff");

  for (Int_t i=0; i<varOrder[0]; ++i) {
    TString coeffName = TString::Format("a_%d",i);
    RooRealVar* a_i = new RooRealVar(coeffName,coeffName,0,-10,10); // 
    xCoeff.add( *a_i );
  } //xCoeff.printMultiline(cout,varOrder[0],kTRUE); xCoeff.printValue(cout);
  for (Int_t j=0; j<varOrder[1]; ++j) {
    TString coeffName = TString::Format("b_%d",j);
    RooRealVar* b_j = new RooRealVar(coeffName,coeffName,0,-10,10); //
    yCoeff.add( *b_j );
  }

  //RooBernstein* xPdf = new RooBernstein("xBernstein", "xBernstein", x, xCoeff) ;
  //RooBernstein* yPdf = new RooBernstein("yBernstein", "yBernstein", y, yCoeff) ;
  RooChebychev* xPdf = new RooChebychev("xChebychev", "xChebychev", x, xCoeff) ;
  RooChebychev* yPdf = new RooChebychev("yChebychev", "yChebychev", y, yCoeff) ;

  RooProdPdf* pdf = new RooProdPdf( TString::Format("%s_prod_%s",xPdf->GetTitle(),yPdf->GetTitle()), TString::Format("%s * %s",xPdf->GetTitle(),yPdf->GetTitle()), RooArgSet(*xPdf,*yPdf) ); 
  //pdf->printMetaArgs(cout); pdf->printMultiline(cout,10,kTRUE);
  //cout <<"pdf->getVal() = " <<pdf->getVal() <<endl;
  ((RooRealVar&)x).setRange("xRange",-1,+1); ((RooRealVar&)y).setRange("yRange",-TMath::Pi(),+TMath::Pi());
  RooFitResult* fitres = pdf->fitTo(*hist, RooFit::SumW2Error(kTRUE), RooFit::Minos(kTRUE), RooFit::Save(kTRUE), RooFit::Verbose(kTRUE), RooFit::PrintLevel(3), RooFit::Range("xRange","yRange")); // solution for bug: RooFit::Optimize(1), RooFit::Minimizer("Minuit2")
  fitres->Print("v");

  return pdf ;
}


// with ROOT
#include "TF2.h"
#include "RooTFnPdfBinding.h"

Double_t polXY(Double_t* var, Double_t* par) {
  Double_t func = par[0];
  for (Int_t i=1; i<=varOrder[0]; ++i) 
    func += par[i]*TMath::Power(var[0],i);
  for (Int_t j=1; j<=varOrder[1]; ++j) { Int_t jOrder = j*(varOrder[0] + 1);
    func += par[jOrder]*TMath::Power(var[1],j);
    for (Int_t i=1; i<=varOrder[0]; ++i) 
      func += par[jOrder + i]*TMath::Power(var[0],i)*TMath::Power(var[1],j);
  }

  return func;
}

Double_t thresholdXY(Double_t* var, Double_t* par) {
  
  Double_t func = 1.;
  
  for (Int_t iLimit=0; iLimit<nThrPar; ++iLimit) { 
    Float_t distance = par[iLimit]*fabs(var[iLimit%nVars]-varLimit[iLimit]) ;
    if (distance < TMath::Pi()/2) 
      func *= TMath::Sin(distance);
  }

  return func ;
}

Double_t polXY_times_thresholdXY(Double_t* var, Double_t* par) {
  
  return thresholdXY(var, &par[0]) * polXY(var, &par[nThrPar]) ;
}

RooAbsPdf* twoDFit(RooAbsReal& x, RooAbsReal& y, TH2F* hist, const Int_t psi_nS, const Int_t xOrder, const Int_t yOrder, Float_t& chi2N) {

  varOrder[0] = xOrder;
  varOrder[nVars-1] = yOrder;

  Double_t MPsi_nS = 0; 
  if (psi_nS == 1) 
    MPsi_nS = MJpsi;
  else if (psi_nS == 2)
    MPsi_nS = MPsi2S;
  else {
    cout <<"psi_nS = " <<psi_nS <<" not allowed at the moment. Aborting" <<endl;
    abort();
  }

  Float_t xLow(0), xHigh(0);
  Float_t yLow(0), yHigh(0);

  TString effName = hist->GetName();
  if (effName.Contains("angle",TString::kIgnoreCase)) {
    if (effName.Contains("helicityAngle_vs_KPiSq",TString::kIgnoreCase)) {
      xLow = TMath::Power(MKaon + MPion,2); xHigh = TMath::Power(MBd - MPsi_nS,2);
      yHigh = +1; yLow = -yHigh;
    } else {
      xHigh = +1; xLow = -xHigh;
      yHigh = +TMath::Pi(); yLow = -yHigh;
    }
  }
  else {
    xLow = MKaon + MPion; xHigh = MBd - MPsi_nS;
    yLow = MPsi_nS + MPion; yHigh = MBd - MKaon;
  }
  cout <<"Fitting in range [" <<xLow <<"," <<xHigh <<"] * [" <<yLow <<"," <<yHigh <<"]" <<endl;

  Float_t xMin = hist->GetXaxis()->GetXmin(); Float_t xMax = hist->GetXaxis()->GetXmax();
  Float_t yMin = hist->GetYaxis()->GetXmin(); Float_t yMax = hist->GetYaxis()->GetXmax();
  
  TF2* finalTF2 = 0;


  // pol(x,y)
  TF2 *funcTF2 = new TF2("funcTF2", polXY, xMin, xMax, yMin, yMax, (varOrder[0]+1)*(varOrder[1]+1));
  funcTF2->SetParName(0,"c_{0}");
  // Must keep this loop synchronized with the one in the polXY function
  for (Int_t i=1; i<=varOrder[0]; ++i)
    funcTF2->SetParName(i,TString::Format("a_{%d}",i));
  for (Int_t j=1; j<=varOrder[1]; ++j) { Int_t jOrder = j*(varOrder[0] + 1);
    funcTF2->SetParName(jOrder,TString::Format("b_{%d}",j));
    for (Int_t i=1; i<=varOrder[0]; ++i)  { Int_t ijOrder = jOrder + i;
      funcTF2->SetParName(ijOrder,TString::Format("c_{%d,%d}",i,j)); }
  }
  
  funcTF2->SetRange(xLow, yLow, xHigh, yHigh);
  hist->Fit("funcTF2","WLR"); //V
  finalTF2 = funcTF2;

  // threshold(x,y)
  if (nThrPar > 0) { varLimit[0] = xHigh;
    if (nThrPar > 1) {  varLimit[1] = yHigh;
      if (nThrPar > 2) { varLimit[2] = xLow;
	if (nThrPar > 3) {  varLimit[3] = yLow;
	} } } }
  
  TF2 *thresholdTF2 = new TF2("thresholdTF2", thresholdXY, xMin, xMax, yMin, yMax, nThrPar);
  
  TString smearingName = "smearing range [#pi/2 fraction]";
  if (nThrPar > 0) { thresholdTF2->SetParName(0,"x upper "+smearingName);
    if (nThrPar > 1) { thresholdTF2->SetParName(1,"y upper smearing range [#pi/2 fraction]");
      if (nThrPar > 2) { thresholdTF2->SetParName(2,"x lower smearing range [#pi/2 fraction]"); 
	if (nThrPar > 3) { thresholdTF2->SetParName(3,"y lower smearing range [#pi/2 fraction]"); 
	} } } }
  Float_t rangeMax = 4; // in [#pi/2 fraction]
  for (Int_t iPar=0; iPar<nThrPar; ++iPar) {
    thresholdTF2->SetParLimits(iPar,rangeMax,9999); thresholdTF2->SetParameter(iPar,100);
  }
  
  thresholdTF2->SetRange(xLow, yLow, xHigh, yHigh);
  //hist->Fit("thresholdTF2","WLVR"); finalTF2 = thresholdTF2;
  
  
  // pol(x,y) * threshold(x,y)
  TF2* polTimesThr = new TF2("polXY_times_thresholdXY", polXY_times_thresholdXY, xMin, xMax, yMin, yMax, funcTF2->GetNumberFreeParameters() + nThrPar);
  for (Int_t iPar=0; iPar<polTimesThr->GetNumberFreeParameters(); ++iPar) {
    if (iPar < nThrPar) {
      polTimesThr->SetParameter(iPar,100);
      polTimesThr->SetParName(iPar,thresholdTF2->GetParName(iPar));
      Double_t parMin(0), parMax(0); thresholdTF2->GetParLimits(iPar, parMin, parMax); 
      polTimesThr->SetParLimits(iPar, parMin, parMax); 
    } else {
      polTimesThr->SetParameter(iPar,1);			     
      polTimesThr->SetParName(iPar,funcTF2->GetParName(iPar-nThrPar));
    }
  }
  //hist->Fit("polXY_times_thresholdXY","WLVR"); finalTF2 = polTimesThr;

  //RooProdPdf* pdf = new RooProdPdf( TString::Format("%s_prod_%s",xPdf->GetTitle(),yPdf->GetTitle()), TString::Format("%s * %s",xPdf->GetTitle(),yPdf->GetTitle()), RooArgSet(*xPdf,*yPdf) ); 
  //pdf->printMetaArgs(cout); pdf->printMultiline(cout,10,kTRUE);
  //cout <<"pdf->getVal() = " <<pdf->getVal() <<endl;
  //RooFitResult* fitres = pdf->fitTo(*hist, RooFit::SumW2Error(kTRUE), RooFit::Minos(kTRUE), RooFit::Save(kTRUE), RooFit::Verbose(kTRUE), RooFit::PrintLevel(3)); // solution for bug: RooFit::Optimize(1), RooFit::Minimizer("Minuit2")
  //fitres->Print("v");

  chi2N = finalTF2->GetChisquare()/finalTF2->GetNDF();
  cout <<"\nNormalized chi2 = " <<chi2N <<endl; 

  return RooFit::bindPdf(finalTF2,x,y) ;
}

