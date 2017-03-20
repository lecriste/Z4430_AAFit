//
//  Analysis.C
//
//  A simple AA analysis code.
//  Created by Ivan Heredia de la Cruz on 4/25/16.
//  Developed by Leonardo Cristella
//
// root -l
// .x myPDF.cxx+
// .x allOtherNeededClasses.cxx+
// .x Analysis.C+
//
// or
// time root -l -b -q run_Analysys.sh > log.txt
//
// or
// root -l -b load_Analysys.sh
// to load libraries and run with custom flags
//


#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "TMath.h"
#include "TH2.h"

#include "RooFitResult.h"
#include <vector>
#include <utility> // std::make_pair
#include "TLegend.h"
#include "RooConstVar.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
//#include "RooNDKeysPdf.h"
#include "RooBinning.h"

#include <sys/time.h> // for timeval
#include <sys/times.h> // for tms
#include "TSystem.h" // to get number of CPUs
#include "TStyle.h" // to use gStyle
#include <TFile.h>
#include <TNtupleD.h>
#include "TPaveText.h"
#include "TPaletteAxis.h"

#include "myPDF.h"
#include "Dalitz_contour.h"
#include "effMasses.h"
#include "Angles_contour.h"
#include "sqDalitz_contour.h"
#include "twoDFit.C"
#include "twoD_chiSquare.C"

#include "initialAmpVal.h"

using namespace RooFit ;

timeval start, stop;
clock_t startCPU, stopCPU;
tms startProc, stopProc;

Float_t TH2_offset = 1.6;


void setBinning(const TH2* hist, Float_t& xMin, Float_t& xMax, Float_t& yMin, Float_t& yMax, RooBinning*& xRooBinning, RooBinning*& yRooBinning) {
  xMin = hist->GetXaxis()->GetXmin(); xMax = hist->GetXaxis()->GetXmax(); Int_t xBins = hist->GetNbinsX();
  yMin = hist->GetYaxis()->GetXmin(); yMax = hist->GetYaxis()->GetXmax(); Int_t yBins = hist->GetNbinsY();

  const TArrayD* xArray = hist->GetXaxis()->GetXbins(); const Double_t* xBinning = 0;
  const TArrayD* yArray = hist->GetYaxis()->GetXbins(); const Double_t* yBinning = 0;
  if (xArray) {
    xBinning = xArray->GetArray();
    if (xBinning) {
      xMin = xBinning[0]; xMax = xBinning[xBins];
    }
  }
  if (yArray) {
    yBinning = yArray->GetArray();
    if (yBinning) {
      yMin = yBinning[0]; yMax = yBinning[yBins];
    }
  }

  xRooBinning = new RooBinning(xBins,xMin,xMax);
  if (xBinning) xRooBinning = new RooBinning(xBins,xBinning);
  yRooBinning = new RooBinning(yBins,yMin,yMax);
  if (yBinning) yRooBinning = new RooBinning(yBins,yBinning);
}


void setXY(const RooArgSet* varsSet, RooRealVar*& x, RooRealVar*& y) {
  TIterator* iter = varsSet->createIterator(); // now iter does not point to x yet
  x = (RooRealVar*)varsSet->first(); iter->Next(); // now iter points to x
  y = (RooRealVar*)iter->Next();
}


void chi2N_hist(TFile* file, const TString errTH2_name, const TH2F* TH2, const RooAbsPdf* pdf, const RooRealVar* x, const RooRealVar* y, const TString method, const TString dir, const TString extension) {

  const TH2F* errTH2 = (TH2F*)file->Get( errTH2_name ) ;
  if (!errTH2) {
    cout <<"\nHistogram \"" <<errTH2_name <<"\" not found in TFile \"" <<file->GetName() <<"\"\nChi2 histogram will not be calculated." <<endl;
  } else {
    TH2F* chi2N_TH2 = twoD_chiSquare(TH2, errTH2, pdf, x, y) ;
    TCanvas* chi2N_C = new TCanvas("chi2N_"+errTH2_name+"_C",TString::Format("chi2N for %s",errTH2->GetTitle()),800,600) ;
    chi2N_TH2->Draw("colz"); //chi2N_hist->SetMaximum(1.5);
    chi2N_C->Update();
    TPaletteAxis *palette = (TPaletteAxis*)chi2N_TH2->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.91); palette->SetX2NDC(0.96);
    chi2N_C->Modified(); chi2N_C->Update();
    chi2N_C->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),chi2N_TH2->GetName(),method.Data(),extension.Data()));
  }
}


void plotting(const RooDataHist* hist, const TString name, const RooRealVar* x, const RooRealVar* y, const RooBinning* xRooBinning, const RooBinning* yRooBinning, const RooAbsPdf* pdf, const TString pdfTitle, const Int_t xOrder, const TString xName, const TString xTitle, const Int_t yOrder, const TString yName, const TString yTitle, const Float_t chi2N, const TString method, const TString dir, const TString extension) {

  TString pdfName = name+"PDF";
  Int_t xBins = xRooBinning->numBins(); Float_t xMin = xRooBinning->binLow(0); Float_t xMax = xRooBinning->binHigh(xBins-1);
  Int_t yBins = yRooBinning->numBins(); Float_t yMin = yRooBinning->binLow(0); Float_t yMax = yRooBinning->binHigh(yBins-1);
  TH2F* TH2_fromHist = (TH2F*)hist->createHistogram(name+"TH2_fromHist", *x, Binning(*xRooBinning), YVar(*y, Binning(*yRooBinning)) ) ;
  TH2_fromHist->SetTitle(hist->GetTitle()); TH2_fromHist->SetTitleOffset(TH2_offset,"XY");
  TH2_fromHist->Draw("LEGO");
  Int_t binDivision = 40;
  cout <<"Using " <<binDivision <<" binDivision to plot " <<pdfTitle <<endl;
  TH2F* pdf_TH2 = (TH2F*)pdf->createHistogram("", *x, Binning(binDivision*xBins,xMin,xMax), YVar(*y, Binning(binDivision*yBins,yMin,yMax)) ) ;
  //TF2* tf2 = (TF2*)pdf->asTF( RooArgSet(*x,*y) ); TH2F* pdf_TH2 = (TH2F*)tf2->GetHistogram(); // Moneta's fix
  pdf_TH2->SetName(pdfName+"_TH2"); pdf_TH2->SetTitle(hist->GetTitle()); pdf_TH2->SetTitleOffset(TH2_offset,"XY");
  cout <<"Plotting " <<pdf_TH2->GetName() <<endl;
  pdf_TH2->SetLineColor(kRed);
  //pdf_TH2->Draw("SURF2"); TH2_fromHist->Draw("LEGO same");
  //pdf_TH2->Draw("LEGO2"); TH2_fromHist->Draw("LEGO same");
  TH2_fromHist->Draw("LEGO"); pdf_TH2->Draw("SURF2 same"); TH2_fromHist->Draw("LEGO same");
  //TH2_fromHist->Draw("LEGO"); pdf_TH2->Draw("SURF3 same"); //TH2_fromHist->SetMaximum( 1.2*TH2_fromHist->GetMaximum() );
  //pdf_TH2->Draw("SURF3"); TH2_fromHist->Draw("LEGO same"); pdf_TH2->SetMaximum( 1.5*pdf_TH2->GetMaximum() );
  if (method.Contains("fit",TString::kIgnoreCase)) {
    TPaveText* fitInfo = new TPaveText(0.78,0.8,0.99,0.9,"NDC");
    fitInfo->AddText("Fit function:");
    fitInfo->AddText(TString::Format("- pol_{%d}[%s] * pol_{%d}[%s]",xOrder,xTitle.Data(),yOrder,yTitle.Data()));
    if (chi2N)
      fitInfo->AddText(TString::Format("- #chi^{2}/n.d.f. = %.3f",chi2N));
    fitInfo->SetTextAlign(12); fitInfo->SetShadowColor(0); fitInfo->SetFillColor(0);
    fitInfo->Draw();
  }
  gPad->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),pdf_TH2->GetName(),method.Data(),extension.Data()));
  //return;
  // Projections
  //TH1D* TH2_proj[] = {relEffTH2->ProjectionX(name+"_"+xName,1,xBins), relEffTH2->ProjectionY(name+"_"+yName,1,xBins)};
  TH1D* TH2_proj[] = {TH2_fromHist->ProjectionX(name+"_"+xName), TH2_fromHist->ProjectionY(name+"_"+yName)};
  //(TH1F*)pdf->createHistogram(name+"_"xName, x, Binning(1*xBins,xMin,xMax) ) ; // not working with message: function value is NAN
  TString projectionName = pdfName+"_";
  TH1D* pdf_proj[] = {pdf_TH2->ProjectionX(projectionName+xName),pdf_TH2->ProjectionY(projectionName+yName)};
  for (Int_t iProj=0; iProj<2; ++iProj) {
    TH2_proj[iProj]->SetTitle(TString::Format("Projection of %s",TH2_fromHist->GetTitle()));
    TH2_proj[iProj]->SetMarkerStyle(20); TH2_proj[iProj]->SetLineColor(kBlack); TH2_proj[iProj]->SetTitleOffset(1);
    //TH2_proj[iProj]->SetMinimum(0);
    TH2_proj[iProj]->Draw();
    pdf_proj[iProj]->Scale((1/pdf_proj[iProj]->Integral())*binDivision*(TH2_proj[iProj]->Integral()));
    pdf_proj[iProj]->SetLineColor(kRed); pdf_proj[iProj]->Draw("same");
    gPad->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),TH2_proj[iProj]->GetName(),method.Data(),extension.Data()));
  }
}


void Analysis(Int_t nEvt = 100000, Bool_t generating = kFALSE, Bool_t bkgFlag = kFALSE, Bool_t effFlag = kFALSE, Bool_t B0BarFlag = kFALSE, Int_t bkgMassOrd = 1, Int_t bkgAngOrd = 1, Int_t effMassOrd = 1, Int_t effAngOrd = 1)
{
  cout <<"With cut-based efficiency the linear interpolation (effOrd=1) of masses does not work" <<endl;

  SysInfo_t* s = new SysInfo_t();
  gSystem->GetSysInfo(s);
  Int_t nCPU = s->fCpus;

  /*
    LHCb strategy:
    As in the Belle analysis, our amplitude model includes all known K*0 -> K+pi- resonances with nominal mass within or slightly above the kinematic limit (1593 MeV) in B0 -> psi' K+ pi- decays: K*_0(800), K*_0(1430) for J = 0; K*(892), K*(1410) and K*(1680) for J = 1; K*_2(1430) for J = 2; and K*_3(1780) for J = 3.
  */


  // Helicity amplitudes values
  // Ivan values
  //const Double_t aIvan = 0.64; const Double_t bIvan = 0.;

  //vector< pair<TString, pair< pair<Double_t, Double_t>, pair<Double_t, Double_t> > > > Kstar_spin;
  vector< pair<TString, pair<const Double_t, const Double_t> > > Kstar_spin;
  map< TString, pair<Double_t, Double_t> > helJ_map;

  cout <<"Adding K*(892)..." <<endl;
  Kstar_spin.push_back( make_pair("892_1", make_pair(M892,G892) ) ) ;
  helJ_map["892_1_0"] = make_pair(K892_1_0_a,K892_1_0_b); helJ_map["892_1_p1"] = make_pair(K892_1_p1_a,K892_1_p1_b); helJ_map["892_1_m1"] = make_pair(K892_1_m1_a,K892_1_m1_b); // from Belle
  //helJ_map["892_1_0"] = make_pair(0.775,0.); helJ_map["892_1_p1"] = make_pair(0.159,1.563); helJ_map["892_1_m1"] = make_pair(0.612,2.712); // from EvtGen

  cout <<"Adding K*(800)..." <<endl;
  Kstar_spin.push_back( make_pair("800_0", make_pair(M800,G800) ) ) ;
  helJ_map["800_0_0"] = make_pair(K800_0_0_a,K800_0_0_b);

  cout <<"Adding K*(1410)..." <<endl;
  Kstar_spin.push_back( make_pair("1410_1", make_pair(M1410,G1410) ) ) ;
  helJ_map["1410_1_0"] = make_pair(K1410_1_0_a,K1410_1_0_b); helJ_map["1410_1_p1"] = make_pair(K1410_1_p1_a,K1410_1_p1_b); helJ_map["1410_1_m1"] = make_pair(K1410_1_m1_a,K1410_1_m1_b);

  cout <<"Adding K*(1430)_0..." <<endl;
  Kstar_spin.push_back( make_pair("1430_0", make_pair(M1430_0,G1430_0) ) ) ;
  helJ_map["1430_0_0"] = make_pair(K1430_0_0_a,K1430_0_0_b);
  //helJ_map["1430_0_0"] = make_pair(1.,0.);

  cout <<"Adding K*(1430)_2..." <<endl;
  Kstar_spin.push_back( make_pair("1430_2", make_pair(M1430_2,G1430_2) ) ) ;
  helJ_map["1430_2_0"] = make_pair(K1430_2_0_a,K1430_2_0_b); helJ_map["1430_2_p1"] = make_pair(K1430_2_p1_a,K1430_2_p1_b); helJ_map["1430_2_m1"] = make_pair(K1430_2_m1_a,K1430_2_m1_b);
  /*
    cout <<"Adding K*(1780)_3..." <<endl;
    Kstar_spin.push_back( make_pair("1780_3", make_pair(M1780_3,G1780_3) ) ) ;
    helJ_map["1780_3_0"] = make_pair(K1780_3_0_a,K1780_3_0_b); helJ_map["1780_3_p1"] = make_pair(K1780_3_p1_a,K1780_3_p1_b); helJ_map["1780_3_m1"] = make_pair(K1780_3_m1_a,K1780_3_m1_b);
  */
  /*
    cout <<"Adding K*(2380)_5..." <<endl;
    Kstar_spin.push_back( make_pair("2380_5", make_pair(M2380_5,G2380_5) ) ) ;
    helJ_map["2380_5_0"] = make_pair(1.,0.); helJ_map["2380_5_p1"] = make_pair(0.,0.); helJ_map["2380_5_m1"] = make_pair(0.,0.);
  */
  TString Hel = ""; //Hel = "_hel0"; //Hel = "_noHel0";
  if (Hel.Contains("_hel0"))
    cout <<"with helicity=0 amplitude only\n" <<endl;
  else if (Hel.Contains("_noHel0"))
    cout <<"without helicity=0 amplitude\n" <<endl;

  vector< RooRealVar* > amplitudeRooRealVar;
  vector< TString > varNames;
  RooArgSet amplitudeVars("amplitudeVars_set");
  Double_t aMax = +9999.;

  // set boundaries:
  Double_t aMin = 0.; //aMin = -aMax;
  Double_t bMin = -9999.; //bMin = -TMath::Pi();

  Double_t bMax = -bMin;
  TString helJ[] = {"m1","0","p1"} ;

  Int_t nKstars = Kstar_spin.size();
  for (Int_t iKstar_S=0; iKstar_S<nKstars; ++iKstar_S)
    for (Int_t iHelJ=0; iHelJ<3; ++iHelJ) {
      if (Kstar_spin[iKstar_S].first.Contains("_0") && !helJ[iHelJ].EqualTo("0")) continue ;
      TString name = Kstar_spin[iKstar_S].first + "_" + helJ[iHelJ] ;
      if (helJ_map.find(name) != helJ_map.end()) {

	if (Hel.Contains("_hel0") && !helJ[iHelJ].EqualTo("0"))
	  helJ_map.find(name)->second.first = 0.; // switching off any amplitude with helicity != 0
	else if (Hel.Contains("_noHel0") && helJ[iHelJ].EqualTo("0"))
	  helJ_map.find(name)->second.first = 0.; // switching off the amplitude with helicity == 0

	pair<Double_t, Double_t> a_b = helJ_map.find(name)->second ;
	TString aName = "a"+name; TString bName = "b"+name;
	//a_b.first = aIvan; a_b.second = bIvan;
	amplitudeRooRealVar.push_back( new RooRealVar(aName,aName, a_b.first, aMin, aMax) ) ; varNames.push_back( aName ) ;
	amplitudeRooRealVar.push_back( new RooRealVar(bName,bName, a_b.second, bMin, bMax) ); varNames.push_back( bName ) ;
      }
      else {
	cout <<"Element \"" <<name <<"\" not found in map helJ_map, please check map filling." <<endl;
	return ;
      }
    }

  Int_t NAmpl = amplitudeRooRealVar.size();
  if (NAmpl > 0)
    for (Int_t iVar=0; iVar<NAmpl; ++iVar) {
      amplitudeVars.add( *amplitudeRooRealVar[iVar] ) ;
    }
  else {
    cout <<"amplitudeRooRealVar is empty! Please check" <<endl;
    return;
  }


  // B^0 -> psi(nS) K* -> mu+ mu- K- pi+
  TString massKPi_name = "massKPi"; TString mass2KPi_name = massKPi_name; mass2KPi_name.ReplaceAll("mass","mass2");
  TString massKPi_title = "m(K^{-}#pi^{+})";
  Float_t massKPi_min = 0.6, massKPi_max = 2.2;
  RooRealVar massKPi(massKPi_name, massKPi_title+" [GeV]", TMath::Sqrt(0.7),massKPi_min,massKPi_max);
  //RooRealVar massKPi(massKPi_name,"m(K^{-}#pi^{+}) [GeV]",TMath::Sqrt(0.7),0.,9.2);
  TString mass2KPi_title = massKPi_title; mass2KPi_title.ReplaceAll("m(","m^{2}(");
  RooFormulaVar mass2KPiFor(mass2KPi_name+"For",mass2KPi_title+" [GeV^{2}]","pow(massKPi,2)",massKPi);
  RooRealVar mass2KPi(mass2KPi_name,mass2KPiFor.getTitle(),TMath::Power(massKPi.getVal(),2),TMath::Power(massKPi.getMin(),2),TMath::Power(massKPi.getMax(),2));

  TString massPsiPi_name = "massMuMuPi"; TString mass2PsiPi_name = massPsiPi_name; mass2PsiPi_name.ReplaceAll("mass","mass2");
  TString massPsiPi_title = "m(#psi#pi^{+})";
  Float_t massPsiPi_min = 3.2, massPsiPi_max = 4.9;
  RooRealVar massPsiPi(massPsiPi_name,massPsiPi_title+" [GeV]",4,massPsiPi_min,massPsiPi_max);
  //RooRealVar massPsiPi(massPsiPi_name,massPsiPi_title+" [GeV]",4,0.,99.9);
  TString mass2PsiPi_title = massPsiPi_title; mass2PsiPi_title.ReplaceAll("m(","m^{2}(");
  RooFormulaVar mass2PsiPiFor(mass2PsiPi_name+"For",mass2PsiPi_title+" [GeV^{2}]","pow(massMuMuPi,2)",massPsiPi);
  RooRealVar mass2PsiPi(mass2PsiPi_name,mass2PsiPiFor.getTitle(),TMath::Power(massPsiPi.getVal(),2),TMath::Power(massPsiPi.getMin(),2),TMath::Power(massPsiPi.getMax(),2));

  TString massVars_name = "massVars"; RooArgSet massVars(massKPi, massPsiPi, massVars_name);
  RooArgSet mass2Fors(mass2KPiFor, mass2PsiPiFor);
  TString mass2Vars_name = "mass2Vars"; RooArgSet mass2Vars(mass2KPi, mass2PsiPi, mass2Vars_name);

  TString cosMuMu_title = "cos(#theta_{#psi})";
  RooRealVar cosMuMu("cosMuMu",cosMuMu_title,0.,-1,1); // cosine of the psi(nS) helicity angle
  TString cosKstar_title = "cos(#theta_{K*})";
  RooRealVar cosKstar("cosKstar",cosKstar_title,0.,-1,1); // cosine of the K* helicity angle
  TString phi_title = "#phi";
  RooRealVar phi("phi",phi_title,0.25,-TMath::Pi(),TMath::Pi());
  //RooRealVar phi("phi","#phi",0.25,-2*TMath::Pi(),2*TMath::Pi());
  TString angleVars_name = "angleVars";
  RooArgSet angleVars(cosMuMu, phi, angleVars_name);

  RooArgSet kinematicVars(massVars, angleVars);
  RooArgSet kinematicVars_m2(mass2Vars, angleVars);
  TString sqDalitz_name = "sqDalitz"; RooArgSet sqDalitz(mass2KPi, cosKstar, sqDalitz_name);
  TString sqDalitz1_name = "sqDalitz1"; RooArgSet sqDalitz1(massKPi, cosKstar, sqDalitz1_name);
  TString sqDalitz_v2_name = "sqDalitz_v2"; RooArgSet sqDalitz_v2(mass2PsiPi, cosKstar, sqDalitz_v2_name);
  TString sqDalitz1_v2_name = "sqDalitz1_v2"; RooArgSet sqDalitz1_v2(massPsiPi, cosKstar, sqDalitz1_v2_name);


  // not accessible on cmssusyX
  TString path = "/lustrehome/cristella/work/Z_analysis/exclusive/clean_14ott/original/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/test/sanjay/selector/";
  //path = "/lustre/home/adrianodif/RootFiles/Z4430/";
  //TString inputFileName = "MC_K892_JPsi_Bd2MuMuKpi_2p0Sig_4p0to6p0SB.root"; Bool_t tmva = kFALSE;
  TString inputFileName = "MC_K892_JPsi_Bd2MuMuKpi_B0massConstraint.root"; Bool_t tmva = kFALSE;
  inputFileName = "Data_JPsi_2p0Sig_6p0to9p0SB.root"; tmva = kFALSE;


  TString bdtCut = "0p00"; //bdtCut = "-0p03";
  TString bdtCut_long = "withBDTCutAt"+bdtCut;
  //inputFileName = "TMVApp_data_"+bdtCut_long+"_JPsi_2p0Sig_6p0to9p0SB.root"; 
  if (inputFileName.Contains("TMVA")) {
    tmva = kTRUE;
    path.Append("TMVA/");
  }
  TString prefix = "", postfix = "";
  if (inputFileName.Contains("MC")) {
    prefix = "Generated "; postfix = "Gen";
  }
  TString fullInputFileName = path+inputFileName ;
  TFile *inputFile = TFile::Open( fullInputFileName ); //inputFile = 0;
  //TFile *inputFile = TFile::Open( inputFileName );
  TString datasetsPath = "datasets";
  RooDataSet *dataToFit = 0;
  RooRealVar B0beauty("B0beauty","B^{0} beauty",0,-1.5,1.5);
  RooArgSet kinematicVars_withBeauty(kinematicVars, B0beauty, TString::Format("%s_with%s",kinematicVars.GetName(),B0beauty.GetName())) ;

  inputFile = 0;
  if (!inputFile) {
    cout <<"Warning: unable to open file \"" <<fullInputFileName <<"\"" <<endl;
  } else {
    //TTree *BkgTree = (TTree*)f->Get("BkgTree");
    TNtupleD* dataNTuple = (TNtupleD*)inputFile->Get("AAVars");
    //dataNTuple->Print() ;
    RooDataSet *data = new RooDataSet("data"+postfix, prefix+"data", dataNTuple, kinematicVars);
    cout <<"\nImported TTree with " <<data->numEntries() <<" entries and the following variables:" <<endl ;
    data->printArgs(cout) ; cout <<"\n" ;
    //cout <<"\n"; kinematicVars.Print("extras") ;
    cout <<"\nImporting dataNTuple..." <<endl;
    RooDataSet *data_B0 = new RooDataSet("data"+postfix+"_B0_B0massConstraint",prefix+"B0 data", dataNTuple, kinematicVars_withBeauty, "B0beauty > 0");
    cout <<"\nImported TTree with " <<data_B0->numEntries() <<" B0" <<endl ;
    data_B0->write(TString::Format("%s/%s.txt",datasetsPath.Data(),data_B0->GetName()));
    RooDataSet *data_B0bar = new RooDataSet("data"+postfix+"_B0bar_B0massConstraint",prefix+"B0bar data", dataNTuple, kinematicVars_withBeauty, "B0beauty < 0");
    cout <<"\nImported TTree with " <<data_B0bar->numEntries() <<" B0bar" <<endl ;
    data_B0bar->write(TString::Format("%s/%s.txt",datasetsPath.Data(),data_B0bar->GetName()));
    /*
      RooRealVar testVar = phi;
      RooPlot* test_frame = testVar.frame(Name(testVar.getTitle()+"_frame"),Title("Projection of "+testVar.getTitle())) ;
      data_B0->plotOn(test_frame); data_B0bar->plotOn(test_frame, LineColor(kRed)); test_frame->Draw() ; return;
    */
    dataToFit = data;
    //dataToFit = data_B0;
    //dataToFit = data_B0bat;
  }
  //return ;

  const Double_t dRadB0 = 5.0; const Double_t dRadKs = 1.5;

  myPDF* sigPDF = 0;
  /*
    sigPDF = new myPDF("signal_pdf","Signal pdf", massKPi, cosMuMu, cosKstar, phi,
    a892m1, b892m1, a892z, b892z, a892p1, b892p1) ;
  */
  TString psi_nS = "1"; //psi_nS = "2";

  TString sigName = "Kstar_", sigTitle = "K*s signal";
  if (nKstars == 1) {
    sigName.Append(Kstar_spin.front().first);
    sigTitle.ReplaceAll("K*s","K*");
  } else {
    sigName.ReplaceAll("Kstar_","Kstars");
    for (Int_t iKstar_S=0; iKstar_S<nKstars; ++iKstar_S)
      sigName.Append("__"+Kstar_spin[iKstar_S].first);
  }
  sigName.Append(Hel);

  pair<TString, TString> sigPDF_varNameTitle[4] = {make_pair(massKPi.GetName(),massKPi_title), make_pair(cosMuMu.GetName(),cosMuMu.GetTitle()), make_pair(massPsiPi.GetName(),massPsiPi_title), make_pair(phi.GetName(),phi.GetTitle())};
  sigPDF = new myPDF(sigName, sigTitle,
		     //massKPi, cosMuMu, massPsiPi, phi,
		     (RooRealVar&)(kinematicVars[sigPDF_varNameTitle[0].first]), (RooRealVar&)(kinematicVars[sigPDF_varNameTitle[1].first]), (RooRealVar&)(kinematicVars[sigPDF_varNameTitle[2].first]), (RooRealVar&)(kinematicVars[sigPDF_varNameTitle[3].first]),
                     B0beauty,
		     Kstar_spin, varNames, amplitudeVars, psi_nS, dRadB0, dRadKs) ;

  if (sigPDF && !nKstars) {
    cout <<"sigPDF set up with no K*! Please check" <<endl;
    return ;
  }
  //cout <<"\nsigPDF->getVal() = " <<sigPDF->getVal() <<endl; return;


  RooConstVar mBd("mBd", "m(B^{0})", MBd) ;
  RooConstVar m2Bd("m2Bd", "m^{2}(B^{0})", MBd2) ;
  RooConstVar mKaon("mKaon", "m(K^{-})", MKaon) ;
  RooConstVar m2Kaon("m2Kaon", "m^{2}(K^{-})", MKaon2) ;
  RooConstVar mPion("mPion", "m(#pi^{+})", MPion) ;
  RooConstVar m2Pion("m2Pion", "m^{2}(#pi^{+})", MPion2) ;

  Double_t massMuMu = 0. ;
  if (psi_nS.EqualTo("1")) {
    massMuMu = MJpsi ;
    //massPsiPi.setMax(4.8);
    massPsiPi_title.ReplaceAll("#psi","J/#psi");
    massPsiPi.SetTitle( massPsiPi.getTitle().ReplaceAll("#psi","J/#psi") );
    cosMuMu_title.ReplaceAll("#psi","J/#psi");
    cosMuMu.SetTitle(cosMuMu_title);
  } else if (psi_nS.EqualTo("2")) {
    massMuMu = MPsi2S ;
    //massPsiPi.setMin(3.7);
    massPsiPi_title.ReplaceAll("#psi","#psi'");
    massPsiPi.SetTitle( massPsiPi.getTitle().ReplaceAll("#psi","#psi'") );
  } else {
    cout <<"psi_nS is neither 1 nor 2, please check it." <<endl;
    return ;
  }
  RooConstVar mPsi("mPsi", "m(#mu^{+}#mu^{-})", massMuMu);
  RooConstVar m2Psi("m2Psi", "m^{2}(#mu^{+}#mu^{-})", TMath::Power(massMuMu,2));
  const Double_t smearing = 0. ;
  RooConstVar smear("smear", "smear", smearing) ;

  // B^{0} -> psi(nS) #pi^{+} K^{-}
  //RooAbsPdf* BdToPsiPiK_PHSP = new RooGenericPdf("BdToPsiPiK_PHSP","3-body PHSP","sqrt( pow(massKPi,4) + pow(mPion,4) + pow(mKaon,4) - 2*pow(massKPi,2)*pow(mPion,2) - 2*pow(massKPi,2)*pow(mKaon,2) - 2*pow(mPion,2)*pow(mKaon,2) ) * sqrt( pow(mBd,4) + pow(massKPi,4) + pow(mPsi,4) - 2*pow(mBd,2)*pow(massKPi,2) - 2*pow(mBd,2)*pow(mPsi,2) - 2*pow(massKPi,2)*pow(mPsi,2) ) / (massKPi)", RooArgSet(massKPi,mPion,mKaon,mBd,mPsi)); // variables name used in the formula must be = name of the RooVariables in the list
  //cout <<"\nBdToPsiPiK_PHSP.getVal() =\n" <<BdToPsiPiK_PHSP->getVal() <<endl; return;
  RooAbsPdf* BdToPsiPiK_PHSP = new RooGenericPdf("BdToPsiPiK_PHSP","3-body PHSP","sqrt( pow(mass2KPiFor,2) + pow(m2Pion,2) + pow(m2Kaon,2) - 2*mass2KPiFor*m2Pion - 2*mass2KPiFor*m2Kaon - 2*m2Pion*m2Kaon ) * sqrt( pow(m2Bd,2) + pow(mass2KPiFor,2) + pow(m2Psi,2) - 2*m2Bd*mass2KPiFor - 2*m2Bd*m2Psi - 2*mass2KPiFor*m2Psi ) / sqrt(mass2KPiFor)", RooArgSet(mass2KPiFor,m2Pion,m2Kaon,m2Bd,m2Psi)); // variables name used in the formula must be = RooVariables name in the RooArgSet
  //cout <<"\nBdToPsiPiK_PHSP.getVal() =\n" <<BdToPsiPiK_PHSP->getVal() <<endl; return;

  RooAbsPdf* bkgPDF = BdToPsiPiK_PHSP; bkgPDF = 0;

  Double_t totEvents = 2000; // Generation time does not scale with number of events up to at least 10k events, from 100k yes
  //totEvents *= 2;
  //totEvents *= 2.5;
  //totEvents *= 5;
  //totEvents *= 10;
  totEvents *= 10; totEvents *= 5;
  //totEvents *= 10; totEvents *= 3;
  //totEvents /= 2; totEvents /= 100;
  RooRealVar nSig("nSig", "n_{SIG}", 0, 0, 1E6);
  //nSig.setVal( 10*nSig.getVal() ); // Does not work on the fly
  Float_t purity = 0.8;
  if (tmva) purity = 0.75;
  RooRealVar nBkg("nBkg", "n_{BKG}", nSig.getVal() * (1-purity), 0, 1E6);
  //RooExtendPdf *extendedBkgPDF = new RooExtendPdf("extendedBkgPDF", "Signal 0 PDF", *bkgPDF, nBkg) ;
  //RooPlot* test_frame = massKPi.frame() ; test_frame->SetTitle( "Projection of "+massKPi.getTitle() ); extendedBkgPDF->plotOn(test_frame) ; test_frame->Draw() ; return;

  TString modelName, modelTitle;
  RooAbsPdf* model = 0;

  if (sigPDF) {
    cout <<"Building " <<sigPDF->GetTitle() <<endl;
    nSig.setVal( totEvents );
    model = (RooAbsPdf*)sigPDF;
    if (bkgPDF) { // in case of efficiency model will be overridden
      cout <<"\nAdding " <<bkgPDF->GetTitle() <<endl;
      nSig.setVal( totEvents/2 ); nBkg.setVal( totEvents/2 );
      model = (RooAbsPdf*) new RooAddPdf("","",RooArgList(*sigPDF,*bkgPDF),RooArgList(nSig,nBkg)) ;
      model->SetName( TString::Format("%s__plus__%s",sigPDF->GetName(),bkgPDF->GetName()) );
      model->SetTitle( TString::Format("%s + %s",sigPDF->GetTitle(),bkgPDF->GetTitle()) );
    }
  } else if (bkgPDF) {
    cout <<"Building " <<bkgPDF->GetTitle() <<endl;
    nBkg.setVal( totEvents );
    //modelTitle = TString::Format("%s*(%s)",nBkg.GetTitle(),bkgPDF->GetTitle());
    model = bkgPDF ;
  } else {
    cout <<"Neither sigPDF nor bkgPDF is != 0! Please check";
    return;
  }
  modelName = model->GetName();

  RooRealVar nEvents("nEvents","nEvents",nSig.getVal() + nBkg.getVal()) ; nEvents.setConstant(kTRUE);
  RooFormulaVar sigFrac("sigFraction",TString::Format("%s fraction",nSig.GetTitle()),"nSig/nEvents",RooArgSet(nSig,nEvents));
  RooFormulaVar bkgFrac("bkgFraction",TString::Format("%s fraction",nBkg.GetTitle()),"nBkg/nEvents",RooArgSet(nBkg,nEvents));


  TString noKinConstr = "_noKinConstr";
  //model->SetName( model->GetName()+noKinConstr );
  // Dalitz boundary
  //Dalitz_contour* kinCheck = new Dalitz_contour("kinCheck","kinematic check", massKPi, massPsiPi, kFALSE, psi_nS) ;
  //RooProdPdf* modelWithKinCheck = new RooProdPdf(modelName,model->GetTitle(),RooArgSet(*kinCheck,*model)) ; model = modelWithKinCheck;


  TString dir = "./plots/";
  if (psi_nS.EqualTo("1"))
    dir.Append("JPsi");
  else if (psi_nS.EqualTo("2"))
    dir.Append("psi2S");
  TString extension = ".png"; //extension.Prepend("_test");

  gStyle->SetOptStat( 10 ) ;


  TString anglesScatt_name = "planesAngle_vs_cos_psi2S_helicityAngle";
  Bool_t DalitzEff = kFALSE;

  RooAbsPdf* null = 0;
  RooProdPdf* sbsModel = 0;
  RooConstVar half = RooConstVar("half", "half", 0.5);
  RooConstVar fixSig = RooConstVar("fixSig", "fixSig", purity);

  if (bkgFlag) { // Model background from sidebands
    //TString bkgFileName = path+"Data_JPsi_2p0Sig_4p0to6p0SB.root";
    //TString bkgFileName = path+"Data_JPsi_2p0Sig_5p0to9p0SB.root";
    TString bkgFileName = path+"Data_JPsi_2p0Sig_6p0to9p0SB.root";
    if (tmva)
      bkgFileName = path+"TMVApp_data_"+bdtCut_long+"_JPsi_2p0Sig_6p0to9p0SB.root";
    cout <<"\nOpening \"" <<bkgFileName <<endl;
    TFile *bkgFile = TFile::Open(bkgFileName);
    const Int_t nVars = 2;
    TString sb_name[] = {"sbs","leftSb","rightSb"};
    pair<RooAbsPdf*, Float_t> sbPdf[nVars][3] = {{make_pair(null,0.),make_pair(null,0.),make_pair(null,0.)},{make_pair(null,0.),make_pair(null,0.),make_pair(null,0.)}};
    /*
    // with fit function
    const Int_t m2KPi_order_bkg = 6, m2PsiPi_order_bkg = 4;
    const Int_t cosMuMu_order_bkg = 6, phi_order_bkg = 6;
    const Int_t mKPi_order_bkg = 5, mPsiPi_order_bkg = 4, cosKstar_order_bkg = 7;
    */
    // with RooHistPdf
    const Int_t m2KPi_order_bkg = bkgMassOrd, m2PsiPi_order_bkg = bkgMassOrd;
    const Int_t mKPi_order_bkg = bkgMassOrd, mPsiPi_order_bkg = bkgMassOrd;
    const Int_t cosMuMu_order_bkg = bkgAngOrd, phi_order_bkg = bkgAngOrd;

    pair< pair<TString, pair<RooArgSet*,pair<Int_t,Int_t> > >, pair<TString,pair< pair<TString,TString>,pair<TString,TString> > > > bkgHisto_names[] = {make_pair( make_pair("psi2SPi_vs_KPi_dalitz",make_pair(&mass2Vars,make_pair(m2KPi_order_bkg,m2PsiPi_order_bkg))), make_pair("bkgDalitz",make_pair(make_pair("m2KPi",mass2KPi_title),make_pair("m2PsiPi",mass2PsiPi_title)))), make_pair(make_pair(anglesScatt_name,make_pair(&angleVars,make_pair(cosMuMu_order_bkg,phi_order_bkg))), make_pair("bkgAngles",make_pair(make_pair("cosMuMu",cosMuMu_title),make_pair("phi",phi_title))))};
    // if you use &mass2Fors you get "ERROR:InputArguments -- RooAbsDataStore::initialize(RelEff_psi2SPi_vs_KPi_B0constr): Data set cannot contain non-fundamental types"
    bkgHisto_names[0] = make_pair( make_pair("psi2SPi_vs_KPi_masses",make_pair(&massVars,make_pair(mKPi_order_bkg,mPsiPi_order_bkg))), make_pair("bkgMasses",make_pair(make_pair("mKPi",massKPi_title),make_pair("mPsiPi",massPsiPi_title)))); // need to apply Dalitz border in TMVAClassificationApplication.C
    //bkgHisto_names[0] = make_pair(make_pair("cos_Kstar_helicityAngle_fromMasses_vs_KPiMassSq",make_pair(&sqDalitz,make_pair(m2KPi_order_bkg,cosKstar_order_bkg))), make_pair("bkgSqDalitz",make_pair(make_pair("m2KPi",mass2KPi_title),make_pair("cosKstar",cosKstar_title)))); DalitzEff = kTRUE;
    //bkgHisto_names[0] = make_pair(make_pair("cos_Kstar_helicityAngle_fromMasses_vs_KPiMass",make_pair(&sqDalitz1,make_pair(mKPi_order_bkg,cosKstar_order_bkg))), make_pair("bkgSqDalitz1",make_pair(make_pair("mKPi",massKPi_title),make_pair("cosKstar",cosKstar_title)))); DalitzEff = kFALSE;
    //
    //bkgHisto_names[0] = make_pair(make_pair("cos_Kstar_helicityAngle_fromMasses_vs_psiPiMassSq",make_pair(&sqDalitz_v2,make_pair(m2PsiPi_order_bkg,cosKstar_order_bkg))), make_pair("bkgSqDalitz_v2",make_pair(make_pair("m2PsiPi",mass2PsiPi_title),make_pair("cosKstar",cosKstar_title)))); DalitzEff = kTRUE;
    //bkgHisto_names[0] = make_pair(make_pair("cos_Kstar_helicityAngle_fromMasses_vs_psiPiMass",make_pair(&sqDalitz1_v2,make_pair(mPsiPi_order_bkg,cosKstar_order_bkg))), make_pair("bkgSqDalitz1_v2",make_pair(make_pair("mPsiPi",massPsiPi_title),make_pair("cosKstar",cosKstar_title)))); DalitzEff = kFALSE;
    
    // bkgFile = 0;
    if (bkgFile)
      for (Int_t iVars = 0; iVars < 2; ++iVars) {
	
	RooArgSet* bkgVars = bkgHisto_names[iVars].first.second.first ;
	// Set x and y vars
	RooRealVar* x = 0, *y = 0;
	setXY(bkgVars, x, y);
	//x->printMultiline(cout,99); y->printMultiline(cout,99); //return;
	
	for (Int_t iSb=0; iSb < 1; ++iSb) {
	  TString histName = bkgHisto_names[iVars].first.first+"_"+sb_name[iSb];
	  if (tmva) histName.Append("_BDT");
	  else {
	    histName.ReplaceAll("_masses","");
	    histName.ReplaceAll(sb_name[iSb],"");
	    histName.Append("hardCuts_1B0_sidebands_B0massC");
	  }
	  
	  const TH2F* sbTH2 = (TH2F*)bkgFile->Get( histName ) ;
	  if (!sbTH2) {
	    cout <<"WARNING! No TH2F \"" <<histName <<"\" found in TFile \"" <<bkgFile->GetName() <<"\".\nSkipping " <<histName <<" evaluation" <<endl;
	    continue; }
	  
	  Float_t xMin = 0, xMax = 0; RooBinning* xRooBinning = 0;
	  Float_t yMin = 0, yMax = 0; RooBinning* yRooBinning = 0;
	  setBinning(sbTH2,xMin,xMax,yMin,yMax,xRooBinning,yRooBinning);
	  
	  cout <<"Setting TH2 range to vars ..." <<endl;
	  x->setRange(xMin,xMax);
	  y->setRange(yMin,yMax);
	  
	  // with RooHistPDF
	  cout <<"Creating RooDataHist from " <<sbTH2->GetName() <<" ..." <<endl;
	  RooDataHist* bkgHist = new RooDataHist(sbTH2->GetName(), sbTH2->GetTitle(), *bkgVars, sbTH2) ;
	  TString bkgName = bkgHisto_names[iVars].second.first+"_"+sb_name[iSb];
	  TString bkgType = bkgName; bkgType.ReplaceAll("bkg","");
	  TString first = TString(bkgType,1);
	  if (!bkgType.EqualTo("Dalitz")) first.ToLower();
	  TString bkgtype = bkgType; bkgtype.Remove(0,1); bkgtype.Prepend(first);
	  
	  TString pdfTitle = "bkg("+bkgtype+") pdf";
	  TString method = "map";
	  Int_t xOrder = bkgHisto_names[iVars].first.second.second.first;
	  Int_t yOrder = bkgHisto_names[iVars].first.second.second.second;
	  if (xOrder) method = "interp"+TString::Itoa(xOrder,10); // with RooHistPdf xOrder = yOrder (see previous comment in this commit)
	  RooHistPdf* bkgHistPdf = new RooHistPdf(bkgName+"PDF_"+method, pdfTitle, *bkgVars, *bkgHist, xOrder) ;
	  //RooHistPdf* bkgHistPdf = new RooHistPdf(bkgName+"PDF", pdfTitle, *bkgVars, *bkgHist, xOrder) ;
	  //If last argument is zero, the weight for the bin enclosing the coordinates contained in 'bin' is returned. For higher values, the result is interpolated in the real dimensions of the dataset with an order of interpolation equal to the value provided (more than ? does not work for Dalitz efficiencies, ? for masses efficiencies, ? for angles)
	  bkgHistPdf->setUnitNorm(kTRUE);
	  
	  sbPdf[iVars][iSb].first = bkgHistPdf;
	  // sbPdf[iVars][iSb].first = twoDFit(*x, *y, sbTH2, psi_nS.Atoi(), xOrder, yOrder, sbPdf[iVars][iSb].second); method = "ROOTfit"; // ROOT fit
	  //sbPdf[iVars][iSb].first = twoDFit(*x, *y, bkgHist, psi_nS.Atoi(), xOrder, yOrder, sbPdf[iVars][iSb].second); method = "RooFit"; // RooFit fit
	  const Float_t chi2N = sbPdf[iVars][iSb].second;
	  
	  RooAbsPdf* kinematicCheck(0), *bkgWithKinCheck(0);
	  if (bkgVars->GetName() == massVars_name  ||  bkgVars->GetName() == mass2Vars_name)
	    kinematicCheck = new Dalitz_contour("Dalitz_kinCheck","kinematic check for Dalitz", *x, *y, DalitzEff, psi_nS) ;
	  else if (bkgVars->GetName() == sqDalitz_name  ||  bkgVars->GetName() == sqDalitz1_name)
	    kinematicCheck = new sqDalitz_contour("sqDalitz_kinCheck","kinematic check for square Dalitz", *x, *y, DalitzEff, psi_nS.Atoi()) ;
	  else if (bkgVars->GetName() == angleVars_name)
	    kinematicCheck = new Angles_contour("angles_kinCheck","kinematic check for angles", *x, *y) ;
	  //
	  if (kinematicCheck) {
	    cout <<"Multiplying " <<bkgHistPdf->GetTitle() <<" by " <<kinematicCheck->GetTitle() <<endl;
	    bkgWithKinCheck = new RooProdPdf(TString::Format("%s_withKinCheck",bkgHistPdf->GetName()),TString::Format("%s with kinematic check",bkgHistPdf->GetTitle()),RooArgSet(*kinematicCheck,*sbPdf[iVars][iSb].first)) ;
	    if (bkgWithKinCheck)
	      sbPdf[iVars][iSb].first = bkgWithKinCheck;
	  }
	  
	  TString sbErrTH2_name = sbTH2->GetName(); sbErrTH2_name.Append("_err");
	  // chi2N_hist(bkgFile, sbErrTH2_name, sbTH2, sbPdf[iVars][iSb].first, x, y, method, dir+"/bkg", extension);

	  plotting(bkgHist, bkgName, x, y, xRooBinning, yRooBinning, sbPdf[iVars][iSb].first, pdfTitle, xOrder, bkgHisto_names[iVars].second.second.first.first, bkgHisto_names[iVars].second.second.first.second, yOrder, bkgHisto_names[iVars].second.second.second.first, bkgHisto_names[iVars].second.second.second.second, chi2N, method, dir+"/bkg", extension);
	  
	} // for (Int_t iSb=0; iSb < 3; ++iSb)
      } // for (Int_t iVars=0; iVars < nVars; ++iVars)
    else
      cout <<"WARNING! TFile \"" <<bkgFileName <<"\" could not be opened.\nSkipping background computation" <<endl;
    //return;

    RooAbsPdf* sbsPdf[] = {sbPdf[0][0].first, sbPdf[1][0].first};
    TString sbsName[] = {"Masses","Angles"};
    /*
      for (Int_t iVar=0; iVar<2; ++iVar) // this to use leftSb + rightSb instead of totalSbs pdf
      if (sbPdf[iVar][1].first && sbPdf[iVar][2].first)
      sbsPdf[iVar] = new RooAddPdf("sbs"+sbsName[iVar]+"Pdf","sidebands "+sbsName[iVar]+" p.d.f.",*sbPdf[iVar][1].first,*sbPdf[iVar][2].first,half) ;
    */
    if (sbsPdf[0] && sbsPdf[1]) {
      sbsModel = new RooProdPdf(TString::Format("%s_X_%s",sbsPdf[0]->GetName(),sbsPdf[1]->GetName()),"sidebands p.d.f.",RooArgSet(*sbsPdf[0],*sbsPdf[1]));
      //cout <<"\nsbsModel->getVal() = " <<sbsModel->getVal() <<endl; return;
    }
    //return;
  } // if (bkgFlag)

  // Masses and angles efficiencies
  //TString effFileName = path+"officialMC_noPtEtaCuts_JPsi_Bd2MuMuKPi_2p0Sig_4p0to6p0SB.root";
  //TString effFileName = path+"officialMC_noPtEtaCuts_JPsi_Bd2MuMuKPi_2p0Sig_5p0to9p0SB.root";
  TString effFileName = path+"officialMC_noPtEtaCuts_JPsi_Bd2MuMuKPi_2p0Sig_6p0to9p0SB.root";
  if (tmva) effFileName = path+"TMVApp_MC_"+bdtCut_long+"_JPsi_2p0Sig_6p0to9p0SB.root";
  cout <<"\nOpening \"" <<effFileName <<endl;
  TFile *effFile = TFile::Open(effFileName);
  RooAbsPdf* pdfToCorrect = sigPDF;
  TString pdfToCorr_name = pdfToCorrect->GetName();
  RooProdPdf* modelWithEff = 0;
  /*
  // with fit function
  const Int_t m2KPi_order_relEff = 5, m2PsiPi_order_relEff = 4, cosMuMu_order_relEff = 5, phi_order_relEff = 5;
  const Int_t mKPi_order_relEff = 4, mPsiPi_order_relEff = 4, cosKstar_order_relEff = 4;
  */
  // with RooHistPdf
  const Int_t m2KPi_order_relEff = effMassOrd, m2PsiPi_order_relEff = effMassOrd;
  const Int_t mKPi_order_relEff = effMassOrd, mPsiPi_order_relEff = effMassOrd;
  const Int_t cosMuMu_order_relEff = effAngOrd, phi_order_relEff = effAngOrd;

  pair< pair<TString, pair<RooArgSet*,pair<Int_t,Int_t> > >, pair<TString,pair< pair<TString,TString>,pair<TString,TString> > > > effHisto_names[] = {make_pair( make_pair("RelEff_psi2SPi_vs_KPi_B0constr_Dalitz",make_pair(&mass2Vars,make_pair(m2KPi_order_relEff,m2PsiPi_order_relEff))), make_pair("relEffDalitz",make_pair(make_pair("m2KPi",mass2KPi_title),make_pair("m2PsiPi",mass2PsiPi_title)))), make_pair(make_pair("RelEff_"+anglesScatt_name,make_pair(&angleVars,make_pair(cosMuMu_order_relEff,phi_order_relEff))), make_pair("relEffAngles",make_pair(make_pair("cosMuMu",cosMuMu_title),make_pair("phi",phi_title))))}; DalitzEff = kTRUE;
  // if you use &mass2Fors you get "ERROR:InputArguments -- RooAbsDataStore::initialize(RelEff_psi2SPi_vs_KPi_B0constr): Data set cannot contain non-fundamental types"
  //effHisto_names[0] = make_pair(make_pair("RelEff_psi2SPi_vs_KPi_B0constr",make_pair(&massVars,make_pair(mKPi_order_relEff,mPsiPi_order_relEff))), make_pair("relEffMasses",make_pair(make_pair("mKPi",massKPi_title),make_pair("mPsiPi",massPsiPi_title)))); DalitzEff = kFALSE;
  effHisto_names[0] = make_pair(make_pair("RelEff_psi2SPi_vs_KPi",make_pair(&massVars,make_pair(mKPi_order_relEff,mPsiPi_order_relEff))), make_pair("relEffMasses",make_pair(make_pair("mKPi",massKPi_title),make_pair("mPsiPi",massPsiPi_title)))); DalitzEff = kFALSE;
  //effHisto_names[0] = make_pair(make_pair("RelEff_cos_Kstar_helicityAngle_vs_KPiSq_varBins",make_pair(&sqDalitz,make_pair(m2KPi_order_relEff,cosKstar_order_relEff))), make_pair("relEffSqDalitz",make_pair(make_pair("m2KPi",mass2KPi_title),make_pair("cosKstar",cosKstar_title)))); DalitzEff = kTRUE;

  pair<RooAbsPdf*, Float_t> effPdf[] = {make_pair(null,0.),make_pair(null,0.)};

  //effFile = 0;
  if (effFile && effFlag)
    for (Int_t iEff=0; iEff <2 ; ++iEff) {
      TString effName = effHisto_names[iEff].second.first;
      TString histName = effHisto_names[iEff].first.first;
      if (tmva) histName.Append("_"+bdtCut_long);
      else histName.Append("_hardCuts_1B0");
      //histName.ReplaceAll("RelEff","RelEffInv"); effName.ReplaceAll("relEff","relInvEff");
      const TH2F* relEffTH2 = (TH2F*)effFile->Get( histName ) ;
      if (!relEffTH2) {
	cout <<"WARNING! No TH2F \"" <<histName <<"\" found in TFile \"" <<effFile->GetName() <<"\".\nSkipping " <<effName <<" correction" <<endl;
	continue; }

      Float_t xMin = 0, xMax = 0; RooBinning* xRooBinning = 0;
      Float_t yMin = 0, yMax = 0; RooBinning* yRooBinning = 0;
      setBinning(relEffTH2,xMin,xMax,yMin,yMax,xRooBinning,yRooBinning);


      RooArgSet* effVars = effHisto_names[iEff].first.second.first ;
      // Set x and y vars
      RooRealVar* x = 0, *y = 0;
      setXY(effVars, x, y);
      x->setRange(xMin,xMax);
      y->setRange(yMin,yMax);


      // with RooHistPDF
      RooDataHist* relEffHist = new RooDataHist(relEffTH2->GetName(), relEffTH2->GetTitle(), *effVars, relEffTH2) ;
      /* // already done in the PROOF macro
	 for (Int_t iBin=0; iBin<relEffHist->numEntries(); ++iBin) { // remove bins whose center is out of Dalitz border
	 *effVars = *(relEffHist->get( iBin ));
	 //cout <<"x = " <<mass2KPi.getVal() <<", y = " <<mass2PsiPi.getVal() <<":\noriginal value = " <<relEffHist->weight(*effVars) <<", corrected value = " <<Dalitz_contour_host(mass2KPi.getVal(), mass2PsiPi.getVal(), kTRUE, psi_nS.Atoi()) * relEffHist->weight(*effVars) <<endl;
	 if (!Dalitz_contour_host(x->getVal(), y->getVal(), DalitzEff, psi_nS.Atoi()))
	 relEffHist->set(*effVars, 0, 0);
	 }
      */
      TString effType = effName; effType.ReplaceAll("relEff","");
      TString first = TString(effType,1); first.ToLower();
      TString efftype = effType; efftype.Remove(0,1); efftype.Prepend(first);
      TString pdfTitle = "#epsilon_{rel}("+efftype+") pdf";
      TString method = "map";
      Int_t xOrder = effHisto_names[iEff].first.second.second.first;
      Int_t yOrder = effHisto_names[iEff].first.second.second.second;
      if (xOrder) method = "interp"+TString::Itoa(xOrder,10); // with RooHistPdf xOrder = yOrder (see previous comment in this commit)
      RooHistPdf* relHistPdf = new RooHistPdf(effName+"PDF_"+method,pdfTitle, *effVars, *relEffHist, xOrder) ; //If last argument is zero, the weight for the bin enclosing the coordinates contained in 'bin' is returned. For higher values, the result is interpolated in the real dimensions of the dataset with an order of interpolation equal to the value provided (more than 10 does not work for Dalitz efficiencies, 9 for masses efficiencies, 10 for angles)
      relHistPdf->setUnitNorm(kTRUE);

      /*
      // with RooNDKeysPDF
      RooRealVar binContent("binContent","binContent",1,0,99);
      effVars->add(binContent);
      RooDataSet* relEffSet = new RooDataSet(relEffTH2->GetName(), relEffTH2->GetTitle(), *effVars, WeightVar("binContent"));
      //RooDataSet* relEffSet = new RooDataSet(relEffTH2->GetName(), relEffTH2->GetTitle(), *effVars);
      for (Int_t i=0; i<=xBins+1; ++i)
      for (Int_t j=0; j<=yBins+1; ++j) {
      x->setVal( relEffTH2->GetXaxis()->GetBinCenter(i) );
      y->setVal( relEffTH2->GetYaxis()->GetBinCenter(j) );
      binContent.setVal( relEffTH2->GetBinContent(i,j) );
      //cout <<"x = " <<x->getVal() <<", y = " <<y->getVal() <<", w = " <<binContent.getVal() <<endl;
      relEffSet->add( *effVars );
      //relEffSet->add( *effVars, binContent.getVal() );
      }
      relEffSet->printArgs(cout); cout <<"; "; relEffSet->printValue(cout); cout <<"; relEffSet->isWeighted() = " <<relEffSet->isWeighted() <<" and relEffSet->sumEntries() = " <<relEffSet->sumEntries() <<endl;
      //TH2F* relEffTH2_fromSet = (TH2F*)relEffSet->createHistogram("relEffTH2_fromSet", *x, Binning(*xRooBinning), YVar(*y, Binning(*yRooBinning)) ) ;
      TH2F* relEffTH2_fromSet = (TH2F*)relEffSet->createHistogram(*x, *y, xBins, yBins, "", "relEffTH2_fromSet") ;
      relEffTH2_fromSet->Draw("LEGO");
      return;
      RooNDKeysPdf* relHistPdf = new RooNDKeysPdf("relHistPdf","relHistPdf", *effVars, *relEffSet, "amdv");
      */

      effPdf[iEff].first = relHistPdf;
      // effPdf[iEff].first = twoDFit(*x, *y, relEffTH2, psi_nS.Atoi(), xOrder, yOrder, effPdf[iEff].second); method = "ROOTfit";
      const Float_t chi2N = effPdf[iEff].second;

      RooAbsPdf* kinematicCheck(0), *effWithKinCheck(0);
      if (effVars->GetName() == massVars_name  ||  effVars->GetName() == mass2Vars_name)
	kinematicCheck = new Dalitz_contour("Dalitz_kinCheck","kinematic check for Dalitz", *x, *y, DalitzEff, psi_nS) ;
      else if (effVars->GetName() == sqDalitz_name)
	kinematicCheck = new sqDalitz_contour("sqDalitz_kinCheck","kinematic check for square Dalitz", *x, *y, DalitzEff, psi_nS.Atoi()) ;
      else if (effVars->GetName() == angleVars_name)
	kinematicCheck = new Angles_contour("angles_kinCheck","kinematic check for angles", *x, *y) ;
      //
      if (kinematicCheck) {
	cout <<"Multiplying " <<relHistPdf->GetTitle() <<" by " <<kinematicCheck->GetTitle() <<endl;
	effWithKinCheck = new RooProdPdf(TString::Format("%s_withKinCheck",relHistPdf->GetName()),TString::Format("%s with kinematic check",relHistPdf->GetTitle()),RooArgSet(*kinematicCheck,*effPdf[iEff].first)) ;
	if (effWithKinCheck)
	  effPdf[iEff].first = effWithKinCheck;
      }

      TString effErrTH2_name = relEffTH2->GetName(); effErrTH2_name.ReplaceAll("Eff","EffErr");
      chi2N_hist(effFile, effErrTH2_name, relEffTH2, effPdf[iEff].first, x, y, method, dir+"/eff", extension);
      //return;
      //relEffTH2->Draw("LEGO");

      plotting(relEffHist, effName, x, y, xRooBinning, yRooBinning, effPdf[iEff].first, pdfTitle, xOrder, effHisto_names[iEff].second.second.first.first, effHisto_names[iEff].second.second.first.second, yOrder, effHisto_names[iEff].second.second.second.first, effHisto_names[iEff].second.second.second.second, chi2N, method, dir+"/eff", extension);
      /*
	x->setVal(4.75); y->setVal(16.);
	cout <<"\neffPdf[iEff].first->getVal(" <<x->getVal() <<", " <<y->getVal() <<") = " <<effPdf[iEff].first->getVal() <<endl;
	cout <<"with relHistPdf->haveUnitNorm() = " <<relHistPdf->haveUnitNorm() <<endl;
	if (relHistPdf->haveUnitNorm()) {
	relHistPdf->setUnitNorm(kFALSE);
	cout <<"\neffPdf[iEff].first->getVal(" <<x->getVal() <<"," <<y->getVal() <<") = " <<effPdf[iEff].first->getVal() <<endl;
	cout <<"with relHistPdf->haveUnitNorm() = " <<relHistPdf->haveUnitNorm() <<endl;
	} else {
	relHistPdf->setUnitNorm(kTRUE);
	cout <<"\neffPdf[iEff].first->getVal(" <<x->.getVal() <<"," <<y->getVal() <<") = " <<effPdf[iEff].first->getVal() <<endl;
	cout <<"with relHistPdf->haveUnitNorm() = " <<relHistPdf->haveUnitNorm() <<endl;
	}
      */

      if (iEff==0) {
 	if (effVars != &massVars) {
          method = "derivedFromFit";
	  /*
	    void deriveMassesPdf(const RooArgSet* massVars, const TString massKPi_name, const TString massPsiPi_name, RooAbsPdf& *pdf, const TString thName, const Int_t xOrder, const TString xName, const TString xTitle, const Int_t yOrder, const TString yName, const TString yTitle, const Float_t chi2N, const TString method, const TString dir, const TString extension) {}

	    TString massesTH_name = "RelEff_psi2SPi_vs_KPi_B0constr";
	    deriveMassesPdf(&massVars, massKPi_name, massPsiPi_name, massesTH_name, xOrder, const TString xName, const TString xTitle, yOrder, const TString yName, const TString yTitle, chi2N, method, dir, extension);
	  */

	  if (kTRUE) { // insert effVars == (m2(K Pi),cos(theta(K*))) check)
	    x = (RooRealVar*)(massVars.find(massKPi_name));
	    y = (RooRealVar*)(massVars.find(massPsiPi_name));

	    effPdf[iEff].first = new sqDalitzToMassesPdf("sqDalitzToMassesPdf","sqDalitzToMassesPdf", *x, *y, effPdf[iEff].first, (RooRealVar*)effVars->find(mass2KPi_name), (RooRealVar*)effVars->find("cosKstar"), massMuMu);

	    TH2F* relEffMassesTH2 = (TH2F*)effFile->Get("RelEff_psi2SPi_vs_KPi_B0constr");
	    if (!relEffMassesTH2) {
	      cout <<"WARNING! No TH2F \"RelEff_psi2SPi_vs_KPi_B0constr\" found in TFile \"" <<effFile->GetName() <<"\".\nSkipping masses efficiency correction" <<endl;
	    } else {
              setBinning(relEffMassesTH2,xMin,xMax,yMin,yMax,xRooBinning,yRooBinning);

	      effName = "relEffMasses";
	      relEffHist = new RooDataHist(relEffMassesTH2->GetName(), relEffMassesTH2->GetTitle(), massVars, relEffMassesTH2) ;
	      plotting(relEffHist, effName, x, y, xRooBinning, yRooBinning, effPdf[iEff].first, pdfTitle, xOrder, effHisto_names[iEff].second.second.first.first, effHisto_names[iEff].second.second.first.second, yOrder, effHisto_names[iEff].second.second.second.first, effHisto_names[iEff].second.second.second.second, 0, method, dir, extension);

	    }
	  } else {

	    // do the same for effVars == (m2(K Pi),m2(J/psi Pi))
	  }
	} // if (effVars != &massVars)

	/*
	// Creating the efficiency as function of masses from the efficiency as function of squared masses, in order to allow the multiplication with sigPDF
	if (DalitzEff) {
	effMasses* massesEffPdf_fromDalitz = new effMasses(TString::Format("%s_fromDalitz",effPdf[iEff].first->GetName()), effPdf[iEff].first->GetTitle(), massKPi, massPsiPi, &mass2KPi, &mass2PsiPi, effPdf[iEff].first);
	//TH1F* massesEffPdf_fromDalitz_mKPTH1 = (TH1F*)massesEffPdf_fromDalitz->createHistogram("massesEffPdf_fromDalitz_mKPTH1", massKPi) ; // not working with message: p.d.f normalization integral is zero or negative
	TH2F* massesEffPdf_fromDalitz_TH2 = (TH2F*)massesEffPdf_fromDalitz->createHistogram("", massKPi, Binning(1*xBins,massKPi_min,massKPi_max), YVar(massPsiPi, Binning(1*yBins,massPsiPi_min,massPsiPi_max)) ) ; massesEffPdf_fromDalitz_TH2->SetName("massesEffPdf_fromDalitz_TH2"); massesEffPdf_fromDalitz_TH2->SetTitleOffset(TH2_offset,"XY");
	massesEffPdf_fromDalitz_TH2->Draw("LEGO");
	gPad->SaveAs(TString::Format("%s/%s%s",dir.Data(),massesEffPdf_fromDalitz_TH2->GetName(),extension.Data()));
	modelWithEff = new RooProdPdf(pdfToCorr_name.Append("__withMassesEff"),TString::Format("(%s)*#epsilon(masses)",pdfToCorrect->GetTitle()),RooArgSet(*pdfToCorrect,*massesEffPdf_fromDalitz)) ;
	}
	*/
	//massKPi.setVal(2.); massPsiPi.setVal(4.);
	//cout <<"massKPi = " <<massKPi.getVal() <<", mass2KPi = " <<mass2KPi.getVal() <<", massPsiPi = " <<massPsiPi.getVal() <<", mass2PsiPi = " <<mass2PsiPi.getVal() <<endl;
	//cout <<"effPdf[iEff].first->getVal() = " <<effPdf[iEff].first->getVal() <<", effPdf[iEff].first->getValV() = " <<effPdf[iEff].first->getValV() <<endl;
	//cout <<"pdfToCorrect->getVal() = " <<pdfToCorrect->getVal() <<", pdfToCorrect->getValV() = " <<pdfToCorrect->getValV() <<endl;
	//cout <<"massKPi = " <<massKPi.getVal() <<", mass2KPi = " <<mass2KPi.getVal() <<", massPsiPi = " <<massPsiPi.getVal() <<", mass2PsiPi = " <<mass2PsiPi.getVal() <<endl;
      } // if (iEff==0)

      //modelWithEff = new RooProdPdf(pdfToCorr_name.Append("__with"+effType+"Eff"),TString::Format("%s * #epsilon("+efftype+")",pdfToCorrect->GetTitle()),RooArgSet(*pdfToCorrect,*effPdf[iEff].first)) ; model = modelWithEff; // replacing sigPdf * eff(M) * eff(A) with sigPdf * [eff(M) * eff(A)] seems to improve the generation time

    } // for (Int_t iEff=0; iEff < 2; ++iEff)
  else
    cout <<"WARNING! TFile \"" <<effFileName <<"\" could not be opened.\nSkipping efficiency correction" <<endl;

  //return;

  RooProdPdf* effModel = 0;
  if (effPdf[0].first && effPdf[1].first)
    effModel = new RooProdPdf(TString::Format("%s_X_%s",effPdf[0].first->GetName(),effPdf[1].first->GetName()),TString::Format("%s * %s",effPdf[0].first->GetTitle(),effPdf[1].first->GetTitle()),RooArgSet(*effPdf[0].first,*effPdf[1].first));

  if ((model != modelWithEff) && effModel) {
    std::cout<<"\nMultiplying " <<effModel->GetTitle() <<" to " <<pdfToCorrect->GetTitle() <<std::endl;
    modelWithEff = new RooProdPdf(TString::Format("%s__with__%s",pdfToCorr_name.Data(),effModel->GetName()),TString::Format("%s * (%s)",pdfToCorrect->GetTitle(),effModel->GetTitle()),*pdfToCorrect,*effModel) ;
    model = modelWithEff;
  }


  // adding background
  RooAddPdf* modelWithBkgHist = 0;
  if (sbsModel) {
    std::cout<<"\nAdding " <<(1 - fixSig.getVal())*100 <<"\% of bkg to pdf" <<std::endl;
    modelWithBkgHist = new RooAddPdf(TString::Format("%s__with__%s",model->GetName(),sbsModel->GetName()),TString::Format("%s + %s",model->GetTitle(),sbsModel->GetTitle()),*model,*sbsModel,fixSig) ;
    model = modelWithBkgHist;
  }
  cout <<"\nmodel->getVal() = " <<model->getVal() <<endl;
  //return;


  Int_t nLegendEntries = 0;

  TString selection = "cutBased";
  if (tmva) selection = bdtCut_long;

  if (generating) {

    if (B0BarFlag)
      B0beauty.setVal(-1.1);
    else
      B0beauty.setVal(+1.1);

    B0beauty.setConstant(kTRUE);

    // Generate toy data from pdf and plot data and p.d.f on frame
    cout <<"\nGenerating " <<nEvents.getVal() <<" events according to " <<model->GetTitle() <<" pdf for " <<model->GetName() <<" with " <<B0beauty.getTitle() <<" = " <<B0beauty.getVal() <<endl;
    timeval genTime;
    gettimeofday(&start, NULL);
    startCPU = times(&startProc);
    //
    TString dataGenName = "Generated_data_from_PDF"; TString dataGen_Name = dataGenName;

    RooDataSet* dataGenPDF = model->generate(kinematicVars, nEvents.getVal(), Verbose(kTRUE), Name(dataGenName)) ; dataGenPDF->SetTitle(dataGenName.ReplaceAll("_"," "));
    //RooDataSet* dataGenPDFB0 = model->generate(kinematicVars_withBeauty, nEvents.getVal(), Verbose(kTRUE), Name(dataGenName)) ; dataGenPDFB0->SetTitle(dataGenName.ReplaceAll("_"," "));

    // RooArgSet genSet = kinematicVars; genSet = kinematicVars_withBeauty;
    // RooDataSet* dataGenPDF = model->generate(genSet, nEvents.getVal(), Verbose(kTRUE), Name(dataGenName)) ; dataGenPDF->SetTitle(dataGenName.ReplaceAll("_"," "));
    //
    stopCPU = times(&stopProc);
    gettimeofday(&stop, NULL);
    timersub(&stop, &start, &genTime);
    Double_t genTimeCPU = (stopCPU - startCPU)*10000;
    Double_t genTimeProc = (stopProc.tms_utime - startProc.tms_utime)*10000 ;

    cout <<"\n" <<nEvents.getVal() <<" events have been genrated in\n" ;
    cout <<"Wallclock time: " << genTime.tv_sec + genTime.tv_usec/1000000.0 << " seconds\n" ;
    cout <<"Total CPU time: " << (genTimeCPU / CLOCKS_PER_SEC) <<" seconds\n" ;
    cout <<"My processes time: " << (genTimeProc / CLOCKS_PER_SEC) << " seconds (differences due to other users' processes on the same CPU)" << endl ;
    //if (nEvents.getVal() > 100000)
    // dataGenPDF->write(TString::Format("%s/%s.txt",datasetsPath.Data(),model->GetName()));

    if (B0BarFlag) {
      dataGenPDF->write(TString::Format("%s/B0bar/%s__%s.txt",datasetsPath.Data(),model->GetName(),selection.Data()));
      //dataGenPDFB0->write(TString::Format("%s/B0Bar/%s__%s__B0Flag.txt",datasetsPath.Data(),model->GetName(),selection.Data()));
    } else {
      dataGenPDF->write(TString::Format("%s/B0/%s__%s.txt",datasetsPath.Data(),model->GetName(),selection.Data()));
      //dataGenPDFB0->write(TString::Format("%s/B0/%s__%s__B0Flag.txt",datasetsPath.Data(),model->GetName(),selection.Data()));
    }

    return;

    dataToFit = dataGenPDF;
  }

  cout <<"\nPlotting data..." <<endl;
  TString plotName = model->GetName();

  Float_t rightMargin = 0.12;
  Float_t cos_limit = 1.02; Float_t cos_margin = 0.02;
  Float_t phi_limit = 3.2; Float_t phi_margin = 0.05;
  if ( totEvents < 50000 ) {
    phi_limit = 3.3; phi_margin = 0.1;
    cos_limit = 1.025; cos_margin = 0.025;
  }
  Int_t cos_bins = 2*cos_limit/cos_margin;
  Int_t phi_bins = 2*phi_limit/phi_margin;

  Float_t KPiMass2_low = 0., KPiMass2_high = 5.; Int_t KPiMass2_bins = 100;
  Float_t MuMuPiMass2_low = 9., MuMuPiMass2_high = 25.; Int_t MuMuPiMass2_bins = 128;

  // Create corresponding m2 dataset + K* helicity angle if absent
  if ( !kinematicVars_m2.contains(cosKstar) )
    kinematicVars_m2.add(cosKstar);

  if (dataToFit) {
    RooDataSet* data_withCos = (RooDataSet*)dataToFit->emptyClone();
    data_withCos->addColumn(cosKstar);

    RooDataSet* data_m2 = new RooDataSet(dataToFit->GetName(), dataToFit->GetTitle(), kinematicVars_m2 ) ;
    for (Int_t iEvent=0; iEvent<dataToFit->numEntries(); ++iEvent) {
      kinematicVars = *(dataToFit->get(iEvent)) ; // this line will propagate the RooRealVars values of the event to all the corresponding RooRealVars in kinematicVars
      //kinematicVars.Print("extras") ; cout <<endl; kinematicVars_m2.Print("extras") ;
      //cout <<"massesEffPdf->printValue(cout) = "; massesEffPdf->printValue(cout); cout <<endl;
      TIterator *massVarsIter = massVars.createIterator() ;
      if (!massVarsIter) {
	cout <<"Cannot create massVars.createIterator()! Please check" <<endl; return;
      }
      RooRealVar* temp = 0;
      while ( (temp = (RooRealVar*)massVarsIter->Next()) ) {
	TString temp_m2Name = temp->GetName(); temp_m2Name.ReplaceAll("mass","mass2"); //cout <<"temp_m2Name = " <<temp_m2Name <<endl; mass2Fors.Print("extras") ;
	RooFormulaVar* temp_m2For = (RooFormulaVar*)mass2Fors.find( temp_m2Name+"For" ) ;
	if (!temp_m2For) {
	  cout <<"temp_m2For = " <<temp_m2For <<"! Please check" <<endl; return;
	} else if ( temp_m2For->dependsOn( *temp ) ) {
	  RooRealVar* temp_m2 = (RooRealVar*)mass2Vars.find(temp_m2Name) ;
	  if (temp_m2)
	    temp_m2->setVal( temp_m2For->getVal() ) ;
	}
      }
      //kinematicVars_m2.Print("extras") ; cout <<endl;
      data_m2->add( kinematicVars_m2 ) ;
      
      //cout <<"\nmassesEffPdf->printValue(cout) = "; massesEffPdf->printValue(cout); cout <<endl;
      
      // set K* helicity angle value
      cosKstar.setVal( cosTheta_FromMasses_host(mass2KPiFor.getVal(), mass2PsiPiFor.getVal(), TMath::Power(massMuMu,2), MBd2, MKaon2, MPion2) );
      data_withCos->add( RooArgSet(kinematicVars,cosKstar) );
    }
    if (dataToFit->numEntries() != data_m2->numEntries()) {
      cout <<"dataToFit->numEntries() (" <<dataToFit->numEntries() <<") != data_m2->numEntries() (" <<data_m2->numEntries() <<")! Please check" <<endl; return;
    }
    if (dataToFit->numEntries() != data_withCos->numEntries()) {
      cout <<"dataToFit->numEntries() (" <<dataToFit->numEntries() <<") != data_withCos->numEntries() (" <<data_withCos->numEntries() <<")! Please check" <<endl; return;
    }
    //return;

    cout <<"\nPlotting angles scatter plot..." <<endl;
    TCanvas* scatter_C = new TCanvas("Angles_scatter_plot","Angles scatter plot",800,600) ; scatter_C->SetRightMargin(rightMargin);
    scatter_C->cd();
    TH2F* scatter = (TH2F*)data_m2->createHistogram("Angles_scatter_plot", cosMuMu, Binning(cos_bins,-cos_limit,cos_limit), YVar(phi, Binning(phi_bins,-phi_limit,phi_limit)) ) ; scatter->SetTitle( TString::Format("Angles scatter plot;%s;%s",cosMuMu.GetTitle(),phi.GetTitle()) ) ;
    gStyle->SetOptStat( 10 ) ;
    scatter->Draw("COLZ");
    scatter_C->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),scatter_C->GetName(),plotName.Data(),extension.Data()));

    cout <<"\nPlotting Dalitz..." <<endl;
    TCanvas* dalitz_C = new TCanvas("Dalitz_C","Dalitz",800,600) ; dalitz_C->SetRightMargin(rightMargin); dalitz_C->cd();
    if ( totEvents < 20000 ) {KPiMass2_bins = 50; MuMuPiMass2_bins = 64;}
    TH2F* dalitz = (TH2F*)data_m2->createHistogram("Dalitz", mass2KPi, Binning(KPiMass2_bins,KPiMass2_low,KPiMass2_high), YVar(mass2PsiPi, Binning(MuMuPiMass2_bins,MuMuPiMass2_low,MuMuPiMass2_high)) ) ; dalitz->SetTitle( TString::Format("Dalitz;%s;%s",mass2KPi.GetTitle(),mass2PsiPi.GetTitle()) ) ;
    gStyle->SetOptStat( 10 ) ;
    dalitz->Draw("COLZ");
    //dalitz->Draw("LEGO");
    dalitz_C->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),dalitz_C->GetTitle(),plotName.Data(),extension.Data()));
    //return;
    cout <<"\nPlotting rectangular Dalitz..." <<endl;
    TCanvas* dalitzRect_C = new TCanvas("DalitzRect_C","Rectangular_Dalitz",800,600) ; dalitzRect_C->SetRightMargin(rightMargin);
    dalitzRect_C->cd();
    TH2F* dalitz_rect = (TH2F*)data_withCos->createHistogram("DalitzRect", massKPi, YVar(cosKstar, Binning(cos_bins,-cos_limit,cos_limit)) ) ; dalitz_rect->SetTitle( TString::Format("Rectangular Dalitz;%s;%s",massKPi.GetTitle(),cosKstar.GetTitle()) ) ;
    gStyle->SetOptStat( 10 ) ;
    dalitz_rect->Draw("COLZ");
    dalitzRect_C->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),dalitzRect_C->GetTitle(),plotName.Data(),extension.Data()));
  }

  vector <RooPlot*> var_frame;
  for (Int_t iVar=0; iVar<kinematicVars.getSize(); ++iVar) {
    RooRealVar var = (RooRealVar&)kinematicVars[sigPDF_varNameTitle[iVar].first];
    var_frame.push_back( var.frame(Title("Projection of "+sigPDF_varNameTitle[iVar].second)) );
    //
    cout <<"\nPlotting " <<var.GetName() <<" ..." <<endl;
    dataToFit->plotOn( var_frame[iVar] ) ;
    TString shortName = var.GetName(); shortName.ReplaceAll("mass","m"); shortName.ReplaceAll("cos","c");
    TCanvas* var_C = new TCanvas( TString::Format("%s_C",var.GetName()), var.GetTitle(), 800,600) ;
    var_C->cd();
    var_frame[iVar]->Draw() ;
    var_C->SaveAs(TString::Format("%s/%s__%s%s",dir.Data(),plotName.Data(),shortName.Data(),extension.Data()));
  }
  //return;
  Int_t fullModelColor = 2; // 2 = kRed
  Int_t bkgColor = fullModelColor;
  TString modelEntry = "Full model";

  cout <<"\nPlotting m(KPi)..." <<endl;
  dataToFit->plotOn(var_frame[0]) ; nLegendEntries++;
  //
  Bool_t fitting = kFALSE; fitting = kTRUE;
  if (!fitting) {
    cout <<"\nPlotting \"" <<model->GetName() <<"\" pdf..." <<endl;
    timeval plotModelTime;
    gettimeofday(&start, NULL);
    startCPU = times(&startProc);
    //
    model->plotOn(var_frame[0],LineColor(fullModelColor),LineStyle(kDashed),Name(modelEntry)) ;
    // 2k events
    // 5h20' with K*(892) + PHSP (1K+1K events)
    // 20' with K*_0(800) + K*_0(1430) + K*_2(1430)
    //
    // 5k events
    // 24' with K*_0(800) + K*_1(892) + K*_1(1410) + K*_0(1430) + K*_2(1430) + K*_3(1780)
    //
    // 20k events and single K* components
    // 220' with K*_0(800) (16'') + K*_1(892) (1035'') + K*_1(1410) (1075'') + K*_0(1430) (15'') + K*_2(1430) (2982'') + K*_3(1780) (6168'')
    //
    if (sigPDF && bkgPDF) {
      model->plotOn(var_frame[0],Components(*bkgPDF),LineColor(bkgColor),LineStyle(kDashed),Name("Bkg")) ;
    }
    //
    stopCPU = times(&stopProc);
    gettimeofday(&stop, NULL);
    timersub(&stop, &start, &plotModelTime);
    Double_t plotModelCPU = (stopCPU - startCPU)*10000;
    Double_t plotModelProc = (stopProc.tms_utime - startProc.tms_utime)*10000 ;
    cout <<"Wallclock time: " << plotModelTime.tv_sec + plotModelTime.tv_usec/1000000.0 << " seconds\n" ;
    cout <<"Total CPU time: " << (plotModelCPU / CLOCKS_PER_SEC) <<" seconds\n" ;
    cout <<"My processes time: " << (plotModelProc / CLOCKS_PER_SEC) << " seconds (differences due to other users' processes on the same CPU)" << endl ;
  }
  nLegendEntries++; // either for generation or fit

  Float_t topRightCorner = 0.9;
  Bool_t plotSingleKstars = kTRUE; plotSingleKstars = kFALSE;
  Float_t yLegLow = topRightCorner -(nLegendEntries+(plotSingleKstars ? nKstars : 1))*0.05 ;
  Float_t xMin = 0.6;
  TLegend* leg = new TLegend(xMin, yLegLow, topRightCorner, topRightCorner); leg->SetFillColor(kWhite);
  leg->AddEntry(dataToFit,"","ep");
  if (!fitting) {
    if (sigPDF && bkgPDF) {
      leg->AddEntry(var_frame[0]->findObject(modelEntry),modelEntry,"l");
      leg->AddEntry(var_frame[0]->findObject("Bkg"),TString::Format("%s (%.1f%%)",bkgPDF->GetTitle(),bkgFrac.getVal()*100),"l");
    } else
      leg->AddEntry(var_frame[0]->findObject(modelEntry),model->GetTitle(),"l");
  }

  // Fitting
  RooFitResult* fitres = 0;
  if (fitting) {
    cout <<"\nFitting ..." <<endl ;
    vector <TString> toSetConst;
    toSetConst.push_back("a892_1_0"); toSetConst.push_back("b892_1_0"); // by convention
    toSetConst.push_back("a892_1_m1"); toSetConst.push_back("b892_1_m1");
    for (Int_t iConst=0; iConst<(Int_t)toSetConst.size(); ++iConst) {
      RooRealVar* toSetConstVar = (RooRealVar*)amplitudeVars.find(toSetConst[iConst]);
      if (toSetConstVar) toSetConstVar->setConstant(kTRUE);
    }

    RooAbsReal* nll = model->createNLL(*dataToFit,Extended(kTRUE),NumCPU(nCPU), Verbose(kTRUE), PrintLevel(3)) ;
    RooArgSet toMinimize(*nll);
    vector <TString> toPenalize;
    //toPenalize.push_back("a892_1_p1"); toPenalize.push_back("b892_1_p1");
    //toPenalize.push_back("a892_1_m1"); toPenalize.push_back("b892_1_m1");
    for (Int_t iPenalty=0; iPenalty<(Int_t)toPenalize.size(); ++iPenalty) {
      RooFormulaVar penalty("penaltyFor_"+toPenalize[iPenalty],"pow(@0 - @1,2)", RooArgSet(*(RooRealVar*)amplitudeVars.find(toPenalize[iPenalty]), RooConstVar(toPenalize[iPenalty]+"_RooConst",toPenalize[iPenalty]+"_RooConst",helJ_map.find(toPenalize[iPenalty])->second.first)));
      toMinimize.add( penalty );
    }
    /*
      RooAddition nllp("nllp","nllp",toMinimize);
      RooMinuit m(nllp); m.setVerbose();
    */
    timeval fitModelTime;
    gettimeofday(&start, NULL);
    startCPU = times(&startProc);
    //
    fitres = model->fitTo(*dataToFit, Hesse(kFALSE), Minos(kFALSE), Save(kTRUE), NumCPU(nCPU), Verbose(kTRUE), PrintLevel(3)) ;
    // 75' with 2k events, 8 Lambda*, 1 NumCPU; 80' with 2k events, 1 K*, 4 NumCPU;
    // with cos(theta_K*) formula in the wrong place:
    // 15h with 2k events, K*(892), 24 CPUs
    //
    //
    //
    // with cos(theta_K*) formula in the right place:
    //
    //
    // with a/b892_1_0 fixed:
    //
    // 1k
    // - 43' for ??  calls (?+?)   with Hesse and 1k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // 2k
    // - 35' for 63  calls (40+12[4]+11[3]) without Hesse and 2k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // - 48' for 94  calls (94+0)  without Hesse and 2k events, 2/4 K*(892) parameters free, 24 CPUs, a and b constrained
    // - 45' for 83  calls (?+?)   with Hesse and 2k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // - 1h  for 118 calls (?+?)   with Hesse and 2k events, 2/4 K*(892) parameters free, 24 CPUs, a and b constrained
    // - 51' for 96  calls (80+16) with Hesse and 2k events, 2/4 K*(892) parameters free, 24 CPUs
    // - 130' for 231 calls (?+?) with Hesse and 2k events, 4/4 K*(892) parameters free, 24 CPUs, a constrained
    // 5k
    // - 40' for 77  calls (63+14) without Hesse and 5k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // 10k
    // - 29' for 51  calls (51+0)  without Hesse and 10k events, 2/4 K*(892) parameters free, 24 CPUs, a and b constrained but b goes at limit (+TMath::Pi)
    // - 27' for 47  calls (47+0)  without Hesse and 10k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // 20k
    // - 24' for 40  calls (30+10) with Hesse and 20k events, 2/4 K*(892) parameters free, 24 CPUs, a and b constrained but b goes at limit (+TMath::Pi)
    // - 34' for 59  calls (49+10) with Hesse and 20k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // - 27' for 49  calls (49+0)  without Hesse and 20k events, 2/4 K*(892) parameters free, 40 CPUs on HPC, a constrained
    // - 28' for 49  calls (49+0)  without Hesse and 20k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // - 19' for 30  calls (30+0)  without Hesse and 20k events, 2/4 K*(892) parameters free, 24 CPUs, a and b constrained but b goes at limit (+TMath::Pi)
    // 60k
    // - 27' for 48  calls (48+0)  without Hesse and 60k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // 100k
    // - 45' for 88  calls and NO  convergence even if plot looks reasonable without Hesse and 100k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // - 56' for 108 calls and NO  convergence even if plot looks reasonable without Hesse and 100k events, 2/4 K*(892) parameters free, 24 CPUs, a and b constrained
    // 200k
    // - 63' for 122 calls and NO  convergence even if plot looks reasonable without Hesse and 200k events, 2/4 K*(892) parameters free, 24 CPUs
    // - 38' for 69  calls (39+7[1]+12[2]+11[3]) without Hesse and 200k events, 2/4 K*(892) parameters free, 24 CPUs, a constrained
    // - 48' for 89  calls and NO  convergence even if plot looks reasonable without Hesse and 200k events, 2/4 K*(892) parameters free, 24 CPUs, a and b constrained
    //
    //
    // with a*[cos(b) + i*sen(b)]:
    // - 64' for 118 calls (?+?) with Hesse and 2k events, 2/4 K*(892) parameters free, 24 CPUs, a and b constrained
    //
    // [1] MIGRAD FAILS TO FIND IMPROVEMENT. MACHINE ACCURACY LIMITS FURTHER IMPROVEMENT. MIGRAD MINIMIZATION HAS CONVERGED.
    // [2] MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.
    // [3] MIGRAD FAILS TO FIND IMPROVEMENT. MIGRAD TERMINATED WITHOUT CONVERGENCE.
    // [4] COVARIANCE MATRIX CALCULATION and ERR MATRIX NOT POS-DEF
    /*
      m.migrad() ;
      m.hesse() ;
      //m.minos() ;
      */
    //
    stopCPU = times(&stopProc);
    gettimeofday(&stop, NULL);
    timersub(&stop, &start, &fitModelTime);
    Double_t fitModelCPU = (stopCPU - startCPU)*10000;
    Double_t fitModelProc = (stopProc.tms_utime - startProc.tms_utime)*10000 ;
    cout <<"\nWallclock time (" <<nCPU <<" CPUs) : " << fitModelTime.tv_sec + fitModelTime.tv_usec/1000000.0 << " seconds\n" ;
    cout <<"Total CPU time (" <<nCPU <<" CPUs) : " << (fitModelCPU / CLOCKS_PER_SEC) <<" seconds\n" ;
    cout <<"My processes time (" <<nCPU <<" CPUs) : " << (fitModelProc / CLOCKS_PER_SEC) << " seconds (differences due to other users' processes on the same CPU)" << endl ;

    fitres->Print("v");
    dataToFit->plotOn(var_frame[0]);
    model->paramOn(var_frame[0], Parameters(amplitudeVars), Layout(xMin,0.99,yLegLow));
    //model->paramOn(var_frame[0], Parameters(RooArgSet(a1600L0S1,b1600L0S1,a1600L1S1,b1600L1S1)), Layout(0.6,0.95,0.9));
    model->plotOn(var_frame[0], LineColor(fullModelColor), Name("4D fit projection")) ;
    leg->AddEntry(var_frame[0]->findObject("4D fit projection"),"4D fit projection","l");
    plotName += "_fit";
    //return;
  }

  //massKPr_frame->SetAxisRange(0,140,"Y") ;

  //Here you can also project other variables.

  if (sigPDF && plotSingleKstars) {
    cout <<"\nIntegrating " <<sigPDF->GetName() <<" over " <<massKPi.GetName() <<"..." <<endl;

    // negligible
    /*
      timeval sigIntegralTime;
      gettimeofday(&start, NULL);
      startCPU = times(&startProc);
      //
      RooAbsReal* sigPDF_integral = sigPDF->createIntegral( kinematicVars ) ;
      Double_t full_signal_integral = sigPDF_integral->getVal(); cout <<"\nFull signal integral = " <<full_signal_integral <<endl;
      //
      stopCPU = times(&stopProc);
      gettimeofday(&stop, NULL);
      timersub(&stop, &start, &sigIntegralTime);
      Double_t sigIntegralCPU = (stopCPU - startCPU)*10000;
      Double_t sigIntegralProc = (stopProc.tms_utime - startProc.tms_utime)*10000 ;
      cout <<"Wallclock time: " << sigIntegralTime.tv_sec + sigIntegralTime.tv_usec/1000000.0 << " seconds\n" ;
      cout <<"Total CPU time: " << (sigIntegralCPU / CLOCKS_PER_SEC) <<" seconds\n" ;
      cout <<"My processes time: " << (sigIntegralProc / CLOCKS_PER_SEC) << " seconds (differences due to other users' processes on the same CPU)" << endl ;
    */
    RooAbsReal* sigPDF_integral = sigPDF->createIntegral( kinematicVars ) ;
    Double_t full_signal_integral = sigPDF_integral->getVal(); cout <<"\nFull signal integral = " <<full_signal_integral <<endl;

    if (nKstars > 1) {
      // 50' with 5 K*
      vector< Double_t > Kstar_integral(nKstars,-1);
      vector<timeval> KstarIntegralTime(nKstars), KstarPlotTime(nKstars);
      vector<Double_t> KstarIntegralCPU(nKstars,0), KstarPlotCPU(nKstars,0);
      vector<Double_t> KstarIntegralProc(nKstars,0), KstarPlotProc(nKstars,0);

      for (Int_t iKstar_S=0; iKstar_S<nKstars; ++iKstar_S) {
	TString Kstar_name = Kstar_spin[iKstar_S].first;

	for (Int_t iVar=0; iVar<(Int_t)varNames.size(); ++iVar) {
	  TString abName = varNames[iVar](1, varNames[iVar].Length());
	  if (helJ_map.find(abName) != helJ_map.end())
	    if ( varNames[iVar].Contains(Kstar_name) ) {
	      pair<Double_t, Double_t> a_b = helJ_map.find(abName)->second ;
	      //cout <<"a_b.first = " <<a_b.first <<", a_b.second = " <<a_b.second <<endl;
	      if ( varNames[iVar].BeginsWith("a") )
		amplitudeVars.setRealValue(varNames[iVar], a_b.first, kTRUE);
	      else if ( varNames[iVar].BeginsWith("b") )
		amplitudeVars.setRealValue(varNames[iVar], a_b.second, kTRUE);
	    }
	    else
	      amplitudeVars.setRealValue(varNames[iVar], 0., kTRUE);
	  else
	    cout <<"Element \"" <<abName <<"\" not found in map helJ_map. Please check map filling." <<endl;
	}

	TString legName = TString::Format("K*_{%s}(%s)", &Kstar_name(Kstar_name.Length() -1), (TString(Kstar_name(0, Kstar_name.Length() -2))).Data());
	cout <<"\n\nIntegrating with " <<legName <<" only..." <<endl;
	// negligible
	/*
	  gettimeofday(&start, NULL);
	  startCPU = times(&startProc);
	  //
	  Kstar_integral[iKstar_S] = (sigPDF->createIntegral(kinematicVars))->getVal() ;
	  //
	  stopCPU = times(&stopProc);
	  gettimeofday(&stop, NULL);
	  timersub(&stop, &start, &KstarIntegralTime[iKstar_S]);
	  KstarIntegralCPU[iKstar_S] = (stopCPU - startCPU)*10000;
	  KstarIntegralProc[iKstar_S] = (stopProc.tms_utime - startProc.tms_utime)*10000 ;
	  cout <<"Wallclock time: " << KstarIntegralTime[iKstar_S].tv_sec + KstarIntegralTime[iKstar_S].tv_usec/1000000.0 << " seconds\n" ;
	  cout <<"Total CPU time: " << (KstarIntegralCPU[iKstar_S] / CLOCKS_PER_SEC) <<" seconds\n" ;
	  cout <<"My processes time: " << (KstarIntegralProc[iKstar_S] / CLOCKS_PER_SEC) << " seconds (differences due to other users' processes on the same CPU)" << endl ;
	  //
	  */
	Kstar_integral[iKstar_S] = (sigPDF->createIntegral(kinematicVars))->getVal() ;

	cout <<"Integral with " <<legName <<" only is: " <<Kstar_integral[iKstar_S] <<endl;
	Double_t fraction = (Kstar_integral[iKstar_S] / full_signal_integral)*(sigFrac.getVal());
	cout <<legName <<" fraction is: " <<fraction*100 <<"%" <<endl;
	cout <<"\nPlotting " <<legName <<" only..." <<endl;

	gettimeofday(&start, NULL);
	startCPU = times(&startProc);
	//
	model->plotOn(var_frame[0],LineColor(iKstar_S + fullModelColor+1), LineStyle(kDashed), Name(Kstar_name), Normalization(fraction,RooAbsReal::Relative)) ;
	//
	stopCPU = times(&stopProc);
	gettimeofday(&stop, NULL);
	timersub(&stop, &start, &KstarPlotTime[iKstar_S]);
	KstarPlotCPU[iKstar_S] = (stopCPU - startCPU)*10000;
	KstarPlotProc[iKstar_S] = (stopProc.tms_utime - startProc.tms_utime)*10000 ;
	cout <<"Wallclock time: " << KstarPlotTime[iKstar_S].tv_sec + KstarPlotTime[iKstar_S].tv_usec/1000000.0 << " seconds\n" ;
	cout <<"Total CPU time: " << (KstarPlotCPU[iKstar_S] / CLOCKS_PER_SEC) <<" seconds\n" ;
	cout <<"My processes time: " << (KstarPlotProc[iKstar_S] / CLOCKS_PER_SEC) << " seconds (differences due to other users' processes on the same CPU)" << endl ;

	leg->AddEntry(var_frame[0]->findObject(Kstar_name),TString::Format("%s (%.1f%%)",legName.Data(),fraction*100),"l");
      }
    } // if (nKstars > 1)
  } // if (sigPDF)

  TCanvas* massKP_C = new TCanvas( TString::Format("%s_C",massKPi.GetName()), massKPi.GetName(), 800,600) ;
  massKP_C->cd();
  var_frame[0]->Draw() ;
  leg->Draw();

  massKP_C->SaveAs(TString::Format("%s/%s%s",dir.Data(),plotName.Data(),extension.Data()));
  gPad->SetLogy();
  massKP_C->SaveAs(TString::Format("%s/%s_logy%s",dir.Data(),plotName.Data(),extension.Data()));

  cout <<"\nEnd of macro!" <<endl;
}
