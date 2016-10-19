//
//  Analysis.C
//
//  A simple AA analysis code.
//  Created by Ivan Heredia de la Cruz on 4/25/16.
//  Developed by Leonardo Cristella
//
// root -l
// .x myPDF.cxx+
// .x Analysis.C
//

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "TMath.h"

#include "RooFitResult.h"
#include <vector>
#include <utility> // std::make_pair
#include "TLegend.h"
#include "RooConstVar.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"

#include <sys/time.h> // for timeval
#include <sys/times.h> // for tms
#include "TSystem.h"

#include "myPDF.h"

using namespace RooFit ;

/*
timeval startFullFitTime, stopFullFullFitTime, totalFullFitTime;
timeval startFullPlotTime, stopFullPlotTime, totalFullPlotTime;
clock_t startFullFitCPU, stopFullFitCPU;
clock_t startFullPlotCPU, stopFullPlotCPU;
tms startFullFitProc, stopFullFitProc;
tms startFullPlotProc, stopFullPlotProc;

vector<timeval> fitTime, plotTime;
vector<clock_t> fitCPU, plotCPU;
vector<tms> fitProc, plotProc;
*/
void Analysis()
{
  SysInfo_t* s = new SysInfo_t();
  gSystem->GetSysInfo(s);
  Int_t nCPU = s->fCpus;

    cout<<TComplex::Exp(TComplex::I()*(TMath::Pi()+0.2))<<endl;
    cout<<TComplex::Exp(TComplex::I()*(-TMath::Pi()+0.2))<<endl;

  /*
  // Lambda*(1600)
  RooRealVar a1600L0S1("a1600L0S1","a1600L0S1",0.64, 0.,1.);
  RooRealVar b1600L0S1("b1600L0S1","b1600L0S1",0, 0.,1.); // see slide 12
  RooRealVar a1600L1S1("a1600L1S1","a1600L1S1",0.64, 0.,1.);
  RooRealVar b1600L1S1("b1600L1S1","b1600L1S1",0, 0.,1.);
  RooRealVar a1600L1S3("a1600L1S3","a1600L1S3",0.64, 0.,1.);
  RooRealVar b1600L1S3("b1600L1S3","b1600L1S3",0, 0.,1.);
  RooRealVar a1600L2S3("a1600L2S3","a1600L2S3",0.64, 0.,1.);
  RooRealVar b1600L2S3("b1600L2S3","b1600L2S3",0, 0.,1.);
  // Lambda*(1670)
  RooRealVar a1670L0S1("a1670L0S1","a1670L0S1",0.0145, 0.,1.);
  RooRealVar b1670L0S1("b1670L0S1","b1670L0S1",0, 0.,1.);
  RooRealVar a1670L1S1("a1670L1S1","a1670L1S1",0.0145, 0.,1.);
  RooRealVar b1670L1S1("b1670L1S1","b1670L1S1",0, 0.,1.);
  RooRealVar a1670L1S3("a1670L1S3","a1670L1S3",0.0145, 0.,1.);
  RooRealVar b1670L1S3("b1670L1S3","b1670L1S3",0, 0.,1.);
  RooRealVar a1670L2S3("a1670L2S3","a1670L2S3",0.0145, 0.,1.);
  RooRealVar b1670L2S3("b1670L2S3","b1670L2S3",0, 0.,1.);
  RooArgSet starResonances(a1600L0S1, b1600L0S1, a1600L1S1, b1600L1S1, a1600L1S3, b1600L1S3, a1600L2S3, b1600L2S3, a1670L0S1, b1670L0S1, a1670L1S1, b1670L1S1, a1670L1S3, b1670L1S3, a1670L2S3, b1670L2S3);
  // Lambda_b -> J/psi Lambda* -> mu+ mu- K- p
  RooRealVar massKPr("massKPr","m(pK^{-}) [GeV]",1.6,1.44,2.52);
  RooRealVar cosLambdaB("cosLambdaB","cosLambdaB",0.5,-1,1); // cosine of the angle between Lambda* and LambdaB in the rest fram of lambdaB
  RooRealVar cosJpsi("cosJpsi","cosJpsi",0.5,-1,1); // cosine of the J/psi helicity angle
  RooRealVar cosLambdaStar("cosLambdaStar","cosLambdaStar",0.5,-1,1); // cosine of the Lambda* helicity angle
  RooRealVar phiMu("phiMu","phi(#mu^{+})",0.25,-TMath::Pi(),TMath::Pi());
  RooRealVar phiK("phiK","phi(K^{+})",0.25,-TMath::Pi(),TMath::Pi());
  RooArgSet kinematcVars(massKPr, cosLambdaB, cosJpsi, cosLambdaStar, phiMu, phiK);

  myPDF sigPDF("signal_pdf","Signal pdf", massKPr, cosLambdaB, cosJpsi, cosLambdaStar, phiMu, phiK,
	       a1600L0S1, b1600L0S1, a1600L1S1, b1600L1S1, a1600L1S3, b1600L1S3, a1600L2S3, b1600L2S3, a1670L0S1, b1670L0S1, a1670L1S1, b1670L1S1, a1670L1S3, b1670L1S3, a1670L2S3, b1670L2S3 ) ;
  */
  /*
  // K*(892)
  RooRealVar a892m1("a892m1","a892m1",0.64, 0.,1.);
  RooRealVar b892m1("b892m1","b892m1",0, -TMath::Pi(), TMath::Pi());
  RooRealVar a892z("a892z","a892z",0.64, 0.,1.);
  RooRealVar b892z("b892z","b892z",0, -TMath::Pi(), TMath::Pi());
  RooRealVar a892p1("a892p1","a892p1",0.64, 0.,1.);
  RooRealVar b892p1("b892p1","b892p1",0, -TMath::Pi(), TMath::Pi());
  */
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
  // Belle B0->J/psi K+ pi- values
  // 5h20' with K*(892) + PHSP (1K+1K events)
  // 20' with K*(800) + K*(1430)0 + K*(1430)2 (2K events)
  TString dataName = "dataset";

  cout <<"Adding K*(892)..." <<endl;
  const Double_t M892 = 0.89581 ; const Double_t G892 = 0.0474; // From PDG neutral only K*(892)
  Kstar_spin.push_back( make_pair("892_1", make_pair(M892,G892) ) ) ;
  helJ_map["892_1_0"] = make_pair(1.,0.); helJ_map["892_1_p1"] = make_pair(0.844,3.14); helJ_map["892_1_m1"] = make_pair(0.196,-1.7);
  dataName += "__k892_1"

  cout <<"Adding K*(800)..." <<endl;
  //const Double_t M800 = 0.682; const Double_t G800 = 0.547; // From PDG
  const Double_t M800 = 0.931; const Double_t G800 = 0.578; // From Belle
  Kstar_spin.push_back( make_pair("800_0", make_pair(M800,G800) ) ) ;
  helJ_map["800_0_0"] = make_pair(1.12,2.3);
  dataName += "__k800_0"

  /*cout <<"Adding K*(1410)..." <<endl;
  const Double_t M1410 = 1.414; const Double_t G1410 = 0.232;
  Kstar_spin.push_back( make_pair("1410_1", make_pair(M1410,G1410) ) ) ;
  helJ_map["1410_1_0"] = make_pair(0.119,0.81); helJ_map["1410_1_p1"] = make_pair(0.123,-1.04); helJ_map["1410_1_m1"] = make_pair(0.036,0.67);
  dataName += "__k1410_1"

  cout <<"Adding K*(1430)_0..." <<endl;
  const Double_t M1430_0 = 1.425; const Double_t G1430_0 = 0.270;
  Kstar_spin.push_back( make_pair("1430_0", make_pair(M1430_0,G1430_0) ) ) ;
  helJ_map["1430_0_0"] = make_pair(0.89,-2.17);
  dataName += "__k1430_0"

  cout <<"Adding K*(1430)_2..." <<endl;
  const Double_t M1430_2 = 1.4324; const Double_t G1430_2 = 0.109; // From PDG neutral only
  Kstar_spin.push_back( make_pair("1430_2", make_pair(M1430_2,G1430_2) ) ) ;
  helJ_map["1430_2_0"] = make_pair(4.66,-0.32); helJ_map["1430_2_p1"] = make_pair(4.65,-3.05); helJ_map["1430_2_m1"] = make_pair(1.26,-1.92);
  dataName += "__k1430_2"*/
  //const Double_t M1780_3 = 1.776; const Double_t G1780_3 = 0.159; // From PDG neutral only
  //Kstar_spin.push_back( make_pair("1780_3", make_pair(M1780_3,G1780_3) ) ) ;
  //helJ_map["1780_3_0"] = make_pair(16.8,-1.43); helJ_map["1780_3_p1"] = make_pair(19.1,2.03); helJ_map["1780_3_m1"] = make_pair(10.2,1.55);
  //dataName += "__k1780_3"

  dataName += ".txt";
  vector< RooRealVar* > amplitudeRooRealVar;
  vector< TString > varNames;
  RooArgSet amplitudeVars("amplitudeVars_set");
  Double_t aMin = -9999.; Double_t aMax = +9999.;
  Double_t bMin = -9999.; Double_t bMax = +9999.;
  TString helJ[] = {"m1","0","p1"} ;

  Int_t nKstar = Kstar_spin.size();
  for (Int_t iKstar_S=0; iKstar_S<nKstar; ++iKstar_S)
    for (Int_t iHelJ=0; iHelJ<3; ++iHelJ) {
      if (Kstar_spin[iKstar_S].first.Contains("_0") && !helJ[iHelJ].EqualTo("0")) continue ;
      TString name = Kstar_spin[iKstar_S].first + "_" + helJ[iHelJ] ;
      if (helJ_map.find(name) != helJ_map.end()) {
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

  for (Int_t iVar=0; iVar<(Int_t)amplitudeRooRealVar.size(); ++iVar) {
    amplitudeVars.add( *amplitudeRooRealVar[iVar] ) ;
  }

  // B^0 -> psi(nS) K* -> mu+ mu- K- pi+
  RooRealVar massKPi("massKPi","m(K^{-}#pi^{+}) [GeV]",1.,0.6,2.2);
  RooRealVar cosMuMu("cosMuMu","cos(J/#psi)",0.,-1,1); // cosine of the psi(nS) helicity angle
  RooRealVar cosKstar("cosKstar","cos(K*)",0.,-1,1); // cosine of the K* helicity angle
  RooRealVar phi("phi","#phi",0.25,-TMath::Pi(),TMath::Pi());
  RooArgSet kinematcVars(massKPi, cosMuMu, cosKstar, phi);
  const Double_t dRadB0 = 5.0; const Double_t dRadKs = 1.5;
  /*
  myPDF sigPDF("signal_pdf","Signal pdf", massKPi, cosMuMu, cosKstar, phi,
	       a892m1, b892m1, a892z, b892z, a892p1, b892p1) ;
  */
  TString psi_nS = "1"; //psi_nS = "2";
  myPDF sigPDF("Kstars_signal","K*s signal", massKPi, cosMuMu, cosKstar, phi,
	       Kstar_spin, varNames, amplitudeVars, psi_nS, dRadB0, dRadKs) ;

  //cout <<"\nsigPDF.getVal() =\n" <<sigPDF.getVal() <<endl; return;

  RooConstVar mBd("mBd", "m(B^{0})", MBd) ;
  RooConstVar mKaon("mKaon", "m(K^{-})", MKaon) ;
  RooConstVar mPion("mPion", "m(#pi^{+})", MPion) ;
  Double_t massMuMu = 0. ;
  if (psi_nS.EqualTo("1")) massMuMu = MJpsi ;
  else if (psi_nS.EqualTo("2")) massMuMu = MPsi2S ;
  else {
    cout <<"psi_nS is neither 1 nor 2, please check it." <<endl;
    return ; }
  RooConstVar mMuMu("mMuMu", "m(#mu^{+}#mu^{-})", massMuMu);
  const Double_t smearing = 0. ;
  RooConstVar smear("smear", "smear", smearing) ;

  // B^{0} -> psi(nS) #pi^{+} K^{-}
  //GooPdf* phaseSpace = new     ThreeBodiesPsiPiK ("BdToMuMuPiK_PHSP", Variable* _x,Variable* _mp,Variable* _m1,Variable* _m2,Variable* _m3);//,Variable* nBk);
  RooAbsPdf* BdToMuMuPiK_PHSP = new RooGenericPdf("BdToMuMuPiK_PHSP","3-body PHSP","sqrt( pow(massKPi,4) + pow(mPion,4) + pow(mKaon,4) - 2*pow(massKPi,2)*pow(mPion,2) - 2*pow(massKPi,2)*pow(mKaon,2) - 2*pow(mPion,2)*pow(mKaon,2) ) * sqrt( pow(mBd,4) + pow(massKPi,4) + pow(mMuMu,4) - 2*pow(mBd,2)*pow(massKPi,2) - 2*pow(mBd,2)*pow(mMuMu,2) - 2*pow(massKPi,2)*pow(mMuMu,2) ) / (massKPi)", RooArgSet(massKPi,mPion,mKaon,mBd,mMuMu)); // variables name used in the formula must be = name of the RooVariables in the list
  //cout <<"\nBdToMuMuPiK_PHSP.getVal() =\n" <<BdToMuMuPiK_PHSP->getVal() <<endl; return;

  RooAbsPdf* bkgPDF = 0;
  //bkgPDF = BdToMuMuPiK_PHSP;
  RooRealVar nBkg("nBkg", "n_{BKG}", 1000, 0, 1E4);
  //RooExtendPdf *extendedBkgPDF = new RooExtendPdf("extendedBkgPDF", "Signal 0 PDF", *bkgPDF, nBkg) ;
  //RooPlot* test_frame = massKPi.frame() ; test_frame->SetTitle( "Projection of "+massKPi.getTitle() ); extendedBkgPDF->plotOn(test_frame) ; test_frame->Draw() ; return;

  RooRealVar nSig("nSig", "n_{SIG}", 2000, 0, 1E4);
  nSig.setVal(10 * nSig.getVal());

  TString modelName = sigPDF.GetName();
  TString modelTitle = TString::Format("%s*(%s)",nSig.GetTitle(),sigPDF.GetTitle());
  RooAbsPdf* model = 0;
  if (bkgPDF) {
    nSig.setVal(nBkg.getVal());
    modelName.Append(TString::Format("_plus_%s",bkgPDF->GetName()));
    modelTitle.Append(TString::Format(" + %s*(%s(",nBkg.GetTitle(),bkgPDF->GetTitle()));
    model = (RooAbsPdf*) new RooAddPdf(modelName,modelTitle,RooArgList(sigPDF,*bkgPDF),RooArgList(nSig,nBkg)) ;
  } else {
    nBkg.setVal( 0 );
    model = (RooAbsPdf*) &sigPDF ;
  }
  RooRealVar nEvents("nEvents","nEvents",nSig.getVal() + nBkg.getVal()) ; nEvents.setConstant(kTRUE);
  RooFormulaVar sigFrac("sigFraction",TString::Format("%s fraction",nSig.GetTitle()),"nSig/nEvents",RooArgSet(nSig,nEvents));
  RooFormulaVar bkgFrac("bkgFraction",TString::Format("%s fraction",nBkg.GetTitle()),"nBkg/nEvents",RooArgSet(nBkg,nEvents));

  RooRealVar* starResonancesMass = 0;
  starResonancesMass = &massKPi;

  RooPlot* massKP_frame = starResonancesMass->frame() ; massKP_frame->SetTitle( "Projection of "+starResonancesMass->getTitle() );
  Int_t nLegendEntries = 0;

  timeval startGenTime, stopGenTime, genTime;
  clock_t startGenCPU, stopGenCPU;
  tms startGenProc, stopGenProc;
  // Generate toy data from pdf and plot data and p.d.f on frame
  cout <<"\nGenerating " <<nEvents.getVal() <<" events according to " <<model->GetName() <<" pdf..." <<endl;

  gettimeofday(&startGenTime, NULL);
  startGenCPU = times(&startGenProc);
  RooDataSet* data = model->generate(kinematcVars, nEvents.getVal()) ;
  stopGenCPU = times(&stopGenProc);
  gettimeofday(&stopGenTime, NULL);
  timersub(&stopGenTime, &startGenTime, &genTime);
  Double_t genTimeCPU = (stopGenCPU - startGenCPU)*10000;
  Double_t procTime = (stopGenProc.tms_utime - startGenProc.tms_utime)*10000 ;
  data->write(dataName.Data());
  cout <<"\n" <<nEvents.getVal() <<" events have been genrated in\n" ;
  cout <<"Wallclock time: " << genTime.tv_sec + genTime.tv_usec/1000000.0 << " seconds\n" ;
  cout <<"Total CPU time: " << (genTimeCPU / CLOCKS_PER_SEC) <<" seconds\n" ;
  cout <<"My processes time: " << (procTime / CLOCKS_PER_SEC) << " seconds (differences due to other processes on the same CPU)" << endl ;

  //return;

  cout <<"\nPlotting data..." <<endl;
  data->plotOn(massKP_frame) ; nLegendEntries++;
  cout <<"\nPlotting " <<model->GetName() <<" pdf..." <<endl;
  Int_t fullModelColor = 2; // 2 = kRed
  Int_t bkgColor = fullModelColor;
  TString modelEntry = "Full model";
  model->plotOn(massKP_frame,LineColor(fullModelColor),Name(modelEntry)) ; nLegendEntries++;
  if (bkgPDF) {
    model->plotOn(massKP_frame,Components(*bkgPDF),LineColor(bkgColor),LineStyle(kDashed),Name("Bkg")) ; nLegendEntries++;
  }
  Float_t topRightCorner = 0.9;
  TLegend* leg = new TLegend(0.6, topRightCorner -(nLegendEntries+nKstar)*0.05, topRightCorner, topRightCorner);
  leg->AddEntry(data,"Data generated","ep");
  if (bkgPDF) {
    leg->AddEntry(massKP_frame->findObject(modelEntry),modelEntry,"l");
    leg->AddEntry(massKP_frame->findObject("Bkg"),TString::Format("%s (%.1f%%)",bkgPDF->GetTitle(),bkgFrac.getVal()*100),"l");
  } else
    leg->AddEntry(massKP_frame->findObject(modelEntry),sigPDF.GetTitle(),"l");


  /*
  RooFitResult* fitres = model->fitTo(*data, Hesse(kTRUE), Minos(kFALSE), Save(kTRUE), NumCPU(nCPU)) ; // 75' with 2000 events, 8 Lambda*, 1 NumCPU; 80' with 2000 events, 1 K*, 4 NumCPU
  fitres->Print("v");
  model->paramOn(massKP_frame);
  //model->paramOn(massKP_frame, Parameters(RooArgSet(a1600L0S1,b1600L0S1,a1600L1S1,b1600L1S1)), Layout(0.6,0.95,0.9));
  model->plotOn(massKP_frame, LineColor(kRed)) ;
  */
  //massKPr_frame->SetAxisRange(0,140,"Y") ;

  //Here you can also project other variables.

  cout <<"\nIntegrating " <<sigPDF.GetName() <<" over " <<starResonancesMass->GetName() <<"..." <<endl;
  RooAbsReal* sigPDF_integral = sigPDF.createIntegral( *starResonancesMass ) ;
  Double_t full_signal_integral = sigPDF_integral->getVal(); cout <<"\nFull signal integral = " <<full_signal_integral <<endl;
  /*
  a1600L1S1.setVal(0.); a1600L1S1.setConstant();
  b1600L1S1.setVal(0.); b1600L1S1.setConstant();
  a1600L1S3.setVal(0.); a1600L1S3.setConstant();
  b1600L1S3.setVal(0.); b1600L1S3.setConstant();
  a1600L2S3.setVal(0.); a1600L2S3.setConstant();
  b1600L2S3.setVal(0.); b1600L2S3.setConstant();
  a1670L0S1.setVal(0.); a1670L0S1.setConstant();
  b1670L0S1.setVal(0.); b1670L0S1.setConstant();
  a1670L1S1.setVal(0.); a1670L1S1.setConstant();
  b1670L1S1.setVal(0.); b1670L1S1.setConstant();
  a1670L1S3.setVal(0.); a1670L1S3.setConstant();
  b1670L1S3.setVal(0.); b1670L1S3.setConstant();
  a1670L2S3.setVal(0.); a1670L2S3.setConstant();
  b1670L2S3.setVal(0.); b1670L2S3.setConstant();
  */
  /*
  a892p1.setVal(0.); a892p1.setConstant();
  b892p1.setVal(0.); b892p1.setConstant();
  a892z.setVal(0.); a892z.setConstant();
  b892z.setVal(0.); b892z.setConstant();

  RooAbsReal* onlyOneResonance = sigPDF.createIntegral(massKPr) ;
  Double_t onlyOneResonance_integral = onlyOneResonance->getVal(); cout <<"\nonlyOneResonance integral = " <<onlyOneResonance_integral <<endl;
  //cout <<"\n1600L0S1 fraction is: " << onlyOneResonance_integral / full_signal_integral * 100 <<"%" <<endl;
  cout <<"\nK*(892)_{-1} fraction is: " << onlyOneResonance_integral / full_signal_integral * 100 <<"%" <<endl;
  */

  if (nKstar > 1) {
    // 50' with 5 K*
    vector< Double_t > Kstar_integral(nKstar,-1);
    for (Int_t iKstar_S=0; iKstar_S<nKstar; ++iKstar_S) {
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
      cosMuMu.setVal(-0.4);
      Kstar_integral[iKstar_S] = (sigPDF.createIntegral(*starResonancesMass))->getVal() ;
      cout <<"CosMuMu : "<<cosMuMu.getVal()<<"Integral with " <<legName <<" only is: " <<Kstar_integral[iKstar_S] <<endl;

      cosMuMu.setVal(0.1);
      Kstar_integral[iKstar_S] = (sigPDF.createIntegral(*starResonancesMass))->getVal() ;
      cout <<"CosMuMu : "<<cosMuMu.getVal()<<"Integral with " <<legName <<" only is: " <<Kstar_integral[iKstar_S] <<endl;

      cosMuMu.setVal(0.3);
      Kstar_integral[iKstar_S] = (sigPDF.createIntegral(*starResonancesMass))->getVal() ;
      cout <<"CosMuMu : "<<cosMuMu.getVal()<<"Integral with " <<legName <<" only is: " <<Kstar_integral[iKstar_S] <<endl;

      RooArgSet observariables(massKPi,cosMuMu,cosKstar,phi);

      cout <<"Total Integral with " <<legName <<" only is: " <<(sigPDF.createIntegral(observariables))->getVal() <<endl;

      Double_t fraction = (Kstar_integral[iKstar_S] / full_signal_integral)*(sigFrac.getVal());
      cout <<legName <<" fraction is: " <<fraction*100 <<"%" <<endl;
      cout <<"\nPlotting " <<legName <<" only..." <<endl;
      model->plotOn(massKP_frame,LineColor(iKstar_S + fullModelColor+1), Name(Kstar_name), Normalization(fraction,RooAbsReal::Relative)) ;
      leg->AddEntry(massKP_frame->findObject(Kstar_name),TString::Format("%s (%.1f%%)",legName.Data(),fraction*100),"l");
    }
  }


  TCanvas* massKP_C = new TCanvas( TString::Format("%s_C",starResonancesMass->GetName()), starResonancesMass->GetName(), 800,600) ;
  massKP_C->cd();
  massKP_frame->Draw() ;
  leg->Draw();

  TString plotName = model->GetName();
  if (psi_nS.EqualTo("1"))
    plotName.Append("_forJpsi");
  else if (psi_nS.EqualTo("2"))
    plotName.Append("_forPsi2S");

  massKP_C->SaveAs(TString::Format("plots/%s_test.png",plotName.Data()));
  gPad->SetLogy();
  massKP_C->SaveAs(TString::Format("plots/%s_logy_test.png",plotName.Data()));
}
