//
//  Analysis.C
//  
//  A simple AA analysis code.
//  Created by Ivan Heredia de la Cruz on 4/25/16.
//  Developed by Leonardo Cristella
//
// root -l
// .x myPDF.cxx+
// .x Analysis.C+
//
// or
// time root -l -b -q run_Analysys.sh > log.txt
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
#include "RooNDKeysPdf.h"

#include <sys/time.h> // for timeval
#include <sys/times.h> // for tms
#include "TSystem.h" // to get number of CPUs
#include "TStyle.h" // to use gStyle
#include <TFile.h>
#include <TNtupleD.h>

#include "myPDF.h"
#include "Dalitz_contour.h"
#include "effMasses.h"
#include "Angles_contour.h"

using namespace RooFit ;

timeval start, stop;
clock_t startCPU, stopCPU;
tms startProc, stopProc;

void Analysis()
{
  SysInfo_t* s = new SysInfo_t();
  gSystem->GetSysInfo(s);
  Int_t nCPU = s->fCpus;
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
  RooArgSet kinematicVars(massKPr, cosLambdaB, cosJpsi, cosLambdaStar, phiMu, phiK);

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
  cout <<"Adding K*(892)..." <<endl;
  //const Double_t M892 = 0.89581 ; const Double_t G892 = 0.0474; // From PDG neutral only K*(892)
  const Double_t M892 = 0.8961 ; const Double_t G892 = 0.0507; // From EvtGen
  Kstar_spin.push_back( make_pair("892_1", make_pair(M892,G892) ) ) ;
  //helJ_map["892_1_0"] = make_pair(1.,0.); helJ_map["892_1_p1"] = make_pair(0.844,3.14); helJ_map["892_1_m1"] = make_pair(0.196,-1.7); // from Belle
  helJ_map["892_1_0"] = make_pair(0.775,0.); helJ_map["892_1_p1"] = make_pair(0.159,1.563); helJ_map["892_1_m1"] = make_pair(0.612,2.712); // from EvtGen
  /*
  cout <<"Adding K*(800)..." <<endl;  
  //const Double_t M800 = 0.682; const Double_t G800 = 0.547; // From PDG
  const Double_t M800 = 0.931; const Double_t G800 = 0.578; // From Belle
  Kstar_spin.push_back( make_pair("800_0", make_pair(M800,G800) ) ) ;
  helJ_map["800_0_0"] = make_pair(1.12,2.3);
  
  cout <<"Adding K*(1410)..." <<endl; 
  const Double_t M1410 = 1.414; const Double_t G1410 = 0.232;
  Kstar_spin.push_back( make_pair("1410_1", make_pair(M1410,G1410) ) ) ;
  helJ_map["1410_1_0"] = make_pair(0.119,0.81); helJ_map["1410_1_p1"] = make_pair(0.123,-1.04); helJ_map["1410_1_m1"] = make_pair(0.036,0.67);
  
  cout <<"Adding K*(1430)_0..." <<endl; 
  const Double_t M1430_0 = 1.425; const Double_t G1430_0 = 0.270;
  Kstar_spin.push_back( make_pair("1430_0", make_pair(M1430_0,G1430_0) ) ) ;
  helJ_map["1430_0_0"] = make_pair(0.89,-2.17);
  //helJ_map["1430_0_0"] = make_pair(1.,0.);
  
  cout <<"Adding K*(1430)_2..." <<endl;
  const Double_t M1430_2 = 1.4324; const Double_t G1430_2 = 0.109; // From PDG neutral only
  Kstar_spin.push_back( make_pair("1430_2", make_pair(M1430_2,G1430_2) ) ) ;
  helJ_map["1430_2_0"] = make_pair(4.66,-0.32); helJ_map["1430_2_p1"] = make_pair(4.65,-3.05); helJ_map["1430_2_m1"] = make_pair(1.26,-1.92);
  
  cout <<"Adding K*(1780)_3..." <<endl;
  const Double_t M1780_3 = 1.776; const Double_t G1780_3 = 0.159; // From PDG neutral only
  Kstar_spin.push_back( make_pair("1780_3", make_pair(M1780_3,G1780_3) ) ) ;
  helJ_map["1780_3_0"] = make_pair(16.8,-1.43); helJ_map["1780_3_p1"] = make_pair(19.1,2.03); helJ_map["1780_3_m1"] = make_pair(10.2,1.55);
  *//*
  cout <<"Adding K*(2380)_5..." <<endl;
  const Double_t M2380_5 = 2.382; const Double_t G2380_5 = 0.178; // From PDG
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
  RooFormulaVar mass2KPiFor(mass2KPi_name+"For","m^{2}(K^{-}#pi^{+}) [GeV^{2}]","pow(massKPi,2)",massKPi);
  RooRealVar mass2KPi(mass2KPi_name,mass2KPiFor.getTitle(),TMath::Power(massKPi.getVal(),2),TMath::Power(massKPi.getMin(),2),TMath::Power(massKPi.getMax(),2));

  TString massPsiPi_name = "massMuMuPi"; TString mass2PsiPi_name = massPsiPi_name; mass2PsiPi_name.ReplaceAll("mass","mass2"); 
  TString massPsiPi_title = "m(#psi#pi^{+})";
  Float_t massPsiPi_min = 3.2, massPsiPi_max = 4.9; 
  RooRealVar massPsiPi(massPsiPi_name,massPsiPi_title+" [GeV]",TMath::Sqrt(23),massPsiPi_min,massPsiPi_max);
  //RooRealVar massPsiPi(massPsiPi_name,massPsiPi_title+" [GeV]",TMath::Sqrt(23),0.,99.9);
  RooFormulaVar mass2PsiPiFor(mass2PsiPi_name+"For","m^{2}(#psi#pi^{+}) [GeV^{2}]","pow(massMuMuPi,2)",massPsiPi);
  RooRealVar mass2PsiPi(mass2PsiPi_name,mass2PsiPiFor.getTitle(),TMath::Power(massPsiPi.getVal(),2),TMath::Power(massPsiPi.getMin(),2),TMath::Power(massPsiPi.getMax(),2));

  RooArgSet massVars(massKPi, massPsiPi);
  RooArgSet mass2Fors(mass2KPiFor, mass2PsiPiFor);
  RooArgSet mass2Vars(mass2KPi, mass2PsiPi);

  TString cosMuMu_title = "cos(#theta_{#psi})";
  RooRealVar cosMuMu("cosMuMu",cosMuMu_title,0.,-1,1); // cosine of the psi(nS) helicity angle
  RooRealVar cosKstar("cosKstar","cos(#theta_{K*})",0.,-1,1); // cosine of the K* helicity angle
  TString phi_title = "#phi";
  RooRealVar phi("phi",phi_title,0.25,-TMath::Pi(),TMath::Pi());
  //RooRealVar phi("phi","#phi",0.25,-2*TMath::Pi(),2*TMath::Pi());
  RooArgSet angleVars(cosMuMu, phi);

  RooArgSet kinematicVars(massVars, angleVars);
  RooArgSet kinematicVars_m2(mass2Vars, angleVars);

  
  // not accessible on cmssusyX
  TString path = "/lustrehome/cristella/work/Z_analysis/exclusive/clean_14ott/original/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/test/sanjay/selector/";
  //TString inputFileName = "MC_K892_JPsi_Bd2MuMuKpi_2p0Sig_4p0to6p0SB.root";
  TString inputFileName = "MC_K892_JPsi_Bd2MuMuKpi_B0massConstraint.root";
  TString fullInputFileName = path+inputFileName ;
  TFile *inputFile = TFile::Open( fullInputFileName ); //inputFile = 0;
  //TFile *inputFile = TFile::Open( inputFileName );
  TString datasetsPath = "datasets";
  RooDataSet *dataToFit = 0;
  inputFile = 0;
  if (!inputFile) {
    cout <<"Warning: unable to open file \"" <<fullInputFileName <<"\"" <<endl;
  } else {
    //TTree *BkgTree = (TTree*)f->Get("BkgTree"); 
    TNtupleD* dataNTuple = (TNtupleD*)inputFile->Get("AAVars");
    //dataNTuple->Print() ;
    RooDataSet *dataGen = new RooDataSet("dataGen","Generated data", dataNTuple, kinematicVars);
    cout <<"\nImported TTree with " <<dataGen->numEntries() <<" entries and the following variables:" <<endl ;
    dataGen->printArgs(cout) ; cout <<"\n" ;
    //cout <<"\n"; kinematicVars.Print("extras") ;
    cout <<"\nImporting dataNTuple..." <<endl;
    RooRealVar B0beauty("B0beauty","B^{0} beauty",0,-2,2);
    RooArgSet kinematicVars_withBeauty(kinematicVars, B0beauty, TString::Format("%s_with%s",kinematicVars.GetName(),B0beauty.GetName())) ;
    RooDataSet *dataGen_B0 = new RooDataSet("dataGen_B0_B0massConstraint","Generated B0 data", dataNTuple, kinematicVars_withBeauty, "B0beauty > 0");
    cout <<"\nImported TTree with " <<dataGen_B0->numEntries() <<" B0" <<endl ;
    dataGen_B0->write(TString::Format("%s/%s.txt",datasetsPath.Data(),dataGen_B0->GetName()));
    RooDataSet *dataGen_B0bar = new RooDataSet("dataGen_B0bar_B0massConstraint","Generated B0bar data", dataNTuple, kinematicVars_withBeauty, "B0beauty < 0");
    cout <<"\nImported TTree with " <<dataGen_B0bar->numEntries() <<" B0bar" <<endl ;
    dataGen_B0bar->write(TString::Format("%s/%s.txt",datasetsPath.Data(),dataGen_B0bar->GetName()));
    //RooPlot* phi_frame = phi.frame(Name(massKPi.getTitle()+"_frame"),Title("Projection of "+massKPi.getTitle())) ; 
    //dataGen_B0->plotOn(phi_frame); dataGen_B0bar->plotOn(phi_frame, LineColor(kRed)); phi_frame->Draw() ; return;

    dataToFit = dataGen_B0;
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
		     Kstar_spin, varNames, amplitudeVars, psi_nS, dRadB0, dRadKs) ;

  if (sigPDF && !nKstars) {
    cout <<"sigPDF set up with no K*! Please check" <<endl;
    return ;
  }
  //cout <<"\nsigPDF->getVal() =\n" <<sigPDF->getVal() <<endl; return;

  
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
  RooConstVar m2Psi("m2Psi", "m^{2}(#mu^{+}#mu^{-})", massMuMu*massMuMu); 
  const Double_t smearing = 0. ; 
  RooConstVar smear("smear", "smear", smearing) ;

  // B^{0} -> psi(nS) #pi^{+} K^{-}
  //RooAbsPdf* BdToPsiPiK_PHSP = new RooGenericPdf("BdToPsiPiK_PHSP","3-body PHSP","sqrt( pow(massKPi,4) + pow(mPion,4) + pow(mKaon,4) - 2*pow(massKPi,2)*pow(mPion,2) - 2*pow(massKPi,2)*pow(mKaon,2) - 2*pow(mPion,2)*pow(mKaon,2) ) * sqrt( pow(mBd,4) + pow(massKPi,4) + pow(mPsi,4) - 2*pow(mBd,2)*pow(massKPi,2) - 2*pow(mBd,2)*pow(mPsi,2) - 2*pow(massKPi,2)*pow(mPsi,2) ) / (massKPi)", RooArgSet(massKPi,mPion,mKaon,mBd,mPsi)); // variables name used in the formula must be = name of the RooVariables in the list
  //cout <<"\nBdToPsiPiK_PHSP.getVal() =\n" <<BdToPsiPiK_PHSP->getVal() <<endl; return;
  RooAbsPdf* BdToPsiPiK_PHSP = new RooGenericPdf("BdToPsiPiK_PHSP","3-body PHSP","sqrt( pow(mass2KPiFor,2) + pow(m2Pion,2) + pow(m2Kaon,2) - 2*mass2KPiFor*m2Pion - 2*mass2KPiFor*m2Kaon - 2*m2Pion*m2Kaon ) * sqrt( pow(m2Bd,2) + pow(mass2KPiFor,2) + pow(m2Psi,2) - 2*m2Bd*mass2KPiFor - 2*m2Bd*m2Psi - 2*mass2KPiFor*m2Psi ) / sqrt(mass2KPiFor)", RooArgSet(mass2KPiFor,m2Pion,m2Kaon,m2Bd,m2Psi)); // variables name used in the formula must be = name of the RooVariables in the list
  //cout <<"\nBdToPsiPiK_PHSP.getVal() =\n" <<BdToPsiPiK_PHSP->getVal() <<endl; return;

  RooAbsPdf* bkgPDF = BdToPsiPiK_PHSP; bkgPDF = 0;

  Double_t totEvents = 2000;
  //totEvents *= 2.5;
  totEvents *= 5;
  //totEvents *= 10;
  //totEvents *= 10; totEvents *= 5;
  //totEvents *= 10; totEvents *= 3;
  //totEvents /= 2;
  RooRealVar nSig("nSig", "n_{SIG}", 0, 0, 1E6);
  //nSig.setVal( 10*nSig.getVal() ); // Does not work on the fly
  RooRealVar nBkg("nBkg", "n_{BKG}", 0, 0, 1E6);
  //RooExtendPdf *extendedBkgPDF = new RooExtendPdf("extendedBkgPDF", "Signal 0 PDF", *bkgPDF, nBkg) ;
  //RooPlot* test_frame = massKPi.frame() ; test_frame->SetTitle( "Projection of "+massKPi.getTitle() ); extendedBkgPDF->plotOn(test_frame) ; test_frame->Draw() ; return;

  TString modelName, modelTitle;
  RooAbsPdf* model = 0;

  if (sigPDF) {
    cout <<"Building " <<sigPDF->GetTitle() <<endl;
    nSig.setVal( totEvents );
    model = (RooAbsPdf*)sigPDF;
    if (bkgPDF) {
      cout <<"Adding " <<bkgPDF->GetTitle() <<endl;
      nSig.setVal( totEvents/2 ); nBkg.setVal( totEvents/2 );
      model = (RooAbsPdf*) new RooAddPdf("","",RooArgList(*sigPDF,*bkgPDF),RooArgList(nSig,nBkg)) ;
      model->SetName( TString::Format("%s__plus__%s",sigPDF->GetName(),bkgPDF->GetName()) );
      model->SetTitle( TString::Format("%s + %s",model->GetTitle(),bkgPDF->GetTitle()) );
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
  Dalitz_contour* kinematicCheck = new Dalitz_contour("kinematicCheck","kinematic check", massKPi, massPsiPi, kFALSE, psi_nS) ;
  //RooProdPdf* modelWithKinCheck = new RooProdPdf(modelName,model->GetTitle(),RooArgSet(*kinematicCheck,*model)) ; model = modelWithKinCheck;


  TString dir = "plots/";
  if (psi_nS.EqualTo("1"))
    dir.Append("JPsi");
  else if (psi_nS.EqualTo("2"))
    dir.Append("psi2S");
  TString extension = ".png"; //extension.Prepend("_test");

  gStyle->SetOptStat( 10 ) ;
  Float_t TH2_offset = 1.6;
 	
  // Masses efficiency
  RooAbsPdf* massesEffPdf = 0;
  TFile *effFile = TFile::Open(path+"officialMC_noPtEtaCuts_JPsi_Bd2MuMuKPi_2p0Sig_4p0to6p0SB.root");
  TString relEffName = "RelEff_psi2SPi_vs_KPi_B0constr_Dalitz"; RooArgSet* massesEffVars = &mass2Vars;
  Bool_t DalitzEff = kTRUE; DalitzEff = kFALSE;
  if (!DalitzEff) {
    relEffName = "RelEff_psi2SPi_vs_KPi_B0constr"; massesEffVars = &massVars;
  }
  //massesEffVars = &mass2Fors; // ERROR:InputArguments -- RooAbsDataStore::initialize(RelEff_psi2SPi_vs_KPi_B0constr): Data set cannot contain non-fundamental types
  TH2F* relEffTH2 = (TH2F*)effFile->Get(relEffName) ;
  //relEffTH2 = 0;
  if (!relEffTH2) {
    cout <<"WARNING! No TH2F \"" <<relEffName <<"\" found in TFile \"" <<effFile->GetName() <<"\".\nExiting" <<endl;
  } else {
    Float_t xMin = relEffTH2->GetXaxis()->GetXmin(), xMax = relEffTH2->GetXaxis()->GetXmax(); Int_t xBins = relEffTH2->GetNbinsX();
    Float_t yMin = relEffTH2->GetYaxis()->GetXmin(), yMax = relEffTH2->GetYaxis()->GetXmax(); Int_t yBins = relEffTH2->GetNbinsY();
    
    // get x and y vars
    RooRealVar* x = (RooRealVar*)(massesEffVars->find(mass2KPi_name));
    RooRealVar* y = (RooRealVar*)(massesEffVars->find(mass2PsiPi_name));
    if (!DalitzEff) {
      x = (RooRealVar*)(massesEffVars->find(massKPi_name));
      y = (RooRealVar*)(massesEffVars->find(massPsiPi_name));
    }

    // with RooHistPDF
    RooDataHist* relEffHist = new RooDataHist(relEffTH2->GetName(), relEffTH2->GetTitle(), *massesEffVars, relEffTH2) ;
    for (Int_t iBin=0; iBin<relEffHist->numEntries(); ++iBin) { // remove bins whose center is out of Dalitz border
      *massesEffVars = *(relEffHist->get( iBin ));
      //cout <<"x = " <<mass2KPi.getVal() <<", y = " <<mass2PsiPi.getVal() <<":\noriginal value = " <<relEffHist->weight(*massesEffVars) <<", corrected value = " <<Dalitz_contour_host(mass2KPi.getVal(), mass2PsiPi.getVal(), kTRUE, psi_nS.Atoi()) * relEffHist->weight(*massesEffVars) <<endl;
      if (!Dalitz_contour_host(x->getVal(), y->getVal(), DalitzEff, psi_nS.Atoi()))
	relEffHist->set(*massesEffVars, 0, 0);
    }

    TString effName = "relEffMasses";
    RooHistPdf* relEffPDF = new RooHistPdf(effName+"PDF","relative #epsilon(masses) pdf", *massesEffVars, *relEffHist, 0) ; //If last argument is zero, the weight for the bin enclosing the coordinates contained in 'bin' is returned. For higher values, the result is interpolated in the real dimensions of the dataset with an order of interpolation equal to the value provided (more than 10 does not work)
    relEffPDF->setUnitNorm(kTRUE);
    Dalitz_contour* kinematicCheck_forEff = new Dalitz_contour("kinematicCheck_forEff","kinematic check", *x, *y, DalitzEff, psi_nS) ;
    RooProdPdf* effWithKinCheck = new RooProdPdf(TString::Format("%s_withKinCheck",relEffPDF->GetName()),TString::Format("%s with kinematic check",relEffPDF->GetTitle()),RooArgSet(*kinematicCheck_forEff,*relEffPDF)) ;
    
    /*
    // with RooNDKeysPDF
    RooRealVar binContent("binContent","binContent",1,0,99);
    massesEffVars->add(binContent);
    RooDataSet* relEffSet = new RooDataSet(relEffTH2->GetName(), relEffTH2->GetTitle(), *massesEffVars, WeightVar("binContent"));
    //RooDataSet* relEffSet = new RooDataSet(relEffTH2->GetName(), relEffTH2->GetTitle(), *massesEffVars);
    for (Int_t i=0; i<=xBins+1; ++i)
      for (Int_t j=0; j<=yBins+1; ++j) {
	x->setVal( relEffTH2->GetXaxis()->GetBinCenter(i) );
	y->setVal( relEffTH2->GetYaxis()->GetBinCenter(j) );
	binContent.setVal( relEffTH2->GetBinContent(i,j) );
	//cout <<"x = " <<x->getVal() <<", y = " <<y->getVal() <<", w = " <<binContent.getVal() <<endl;
	relEffSet->add( *massesEffVars );
	//relEffSet->add( *massesEffVars, binContent.getVal() );
      }
    relEffSet->printArgs(cout); cout <<"; "; relEffSet->printValue(cout); cout <<"; relEffSet->isWeighted() = " <<relEffSet->isWeighted() <<" and relEffSet->sumEntries() = " <<relEffSet->sumEntries() <<endl;
    //TH2F* relEffTH2_fromSet = (TH2F*)relEffSet->createHistogram("relEffTH2_fromSet", *x, Binning(xBins,xMin,xMax), YVar(*y, Binning(yBins,yMin,yMax)) ) ;
    TH2F* relEffTH2_fromSet = (TH2F*)relEffSet->createHistogram(*x, *y, xBins, yBins, "", "relEffTH2_fromSet") ;
    relEffTH2_fromSet->Draw("lego");
    return;
    RooNDKeysPdf* relEffPDF = new RooNDKeysPdf("relEffPDF","relEffPDF", *massesEffVars, *relEffSet, "amdv");
    */

    massesEffPdf = relEffPDF;
    massesEffPdf = effWithKinCheck;

    RooProdPdf* modelWithEff = new RooProdPdf(modelName.Append("__withMassesEff"),TString::Format("(%s)*#epsilon(masses)",model->GetTitle()),RooArgSet(*model,*massesEffPdf)) ;

    //relEffTH2->Draw("lego");
    TH2F* relEffTH2_fromHist = (TH2F*)relEffHist->createHistogram(effName+"TH2_fromHist", *x, Binning(xBins,xMin,xMax), YVar(*y, Binning(yBins,yMin,yMax)) ) ; relEffTH2_fromHist->SetTitle(relEffTH2->GetTitle()); relEffTH2_fromHist->SetTitleOffset(TH2_offset,"XY");
    relEffTH2_fromHist->Draw("lego");
    TH2F* massesEffPdf_TH2 = (TH2F*)massesEffPdf->createHistogram("", *x, Binning(1*xBins,xMin,xMax), YVar(*y, Binning(1*yBins,yMin,yMax)) ) ; massesEffPdf_TH2->SetName(TString::Format("%s_TH2",relEffPDF->GetName())); massesEffPdf_TH2->SetTitle(relEffTH2->GetTitle()); massesEffPdf_TH2->SetTitleOffset(TH2_offset,"XY");
    
    massesEffPdf_TH2->SetLineColor(kRed);
    massesEffPdf_TH2->Draw("surf2");
    //massesEffPdf_TH2->Draw("lego");
    relEffTH2_fromHist->Draw("lego same");
    gPad->SaveAs(TString::Format("%s/%s%s",dir.Data(),massesEffPdf_TH2->GetName(),extension.Data()));

    // Projection X
    //TH1D* relEff_mKP = relEffTH2->ProjectionX("relEff_mKP",1,xBins);
    TH1D* relEff_mKP = relEffTH2_fromHist->ProjectionX("relEff_mKPi"); relEff_mKP->SetTitle(TString::Format("Projection of %s",relEffTH2_fromHist->GetTitle()));
    relEff_mKP->Draw();
    //
    //TH1F* massesEffPdf_mKPHist = (TH1F*)massesEffPdf->createHistogram("massesEffPdf_mKPHist", x, Binning(1*xBins,xMin,xMax) ) ; // not working with message: function value is NAN
    TH1D* massesEffPdf_mKPHist = massesEffPdf_TH2->ProjectionX(); //"massesEffPdf_mKPHist",1,xBins
    massesEffPdf_mKPHist->Scale(relEff_mKP->Integral()); massesEffPdf_mKPHist->SetLineColor(kRed); massesEffPdf_mKPHist->Draw("same");
    gPad->SaveAs(TString::Format("%s/%s%s",dir.Data(),relEff_mKP->GetName(),extension.Data()));

    // Projection Y
    //TH1D* relEff_mPsiP = relEffTH2->ProjectionX("relEff_mPsiP",1,xBins);
    TH1D* relEff_mPsiP = relEffTH2_fromHist->ProjectionY("relEff_mPsiPi"); relEff_mPsiP->SetTitle(TString::Format("Projection of %s",relEffTH2_fromHist->GetTitle()));
    relEff_mPsiP->Draw();
    //
    //TH1F* massesEffPdf_mPsiPiHist = (TH1F*)massesEffPdf->createHistogram("massesEffPdf_mPsiPiHist", y, Binning(1*yBins,yMin,yMax) ) ; massesEffPdf_mPsiPiHist->SetLineColor(kRed); massesEffPdf_mPsiPiHist->Draw();
    TH1D* massesEffPdf_mPsiPHist = massesEffPdf_TH2->ProjectionY(); //"massesEffPdf_mPsiPHist",1,yBins
    massesEffPdf_mPsiPHist->Scale(relEff_mPsiP->Integral()); massesEffPdf_mPsiPHist->SetLineColor(kRed); massesEffPdf_mPsiPHist->Draw("same");
    gPad->SaveAs(TString::Format("%s/%s%s",dir.Data(),relEff_mPsiP->GetName(),extension.Data()));
    //return;
    /*
    x->setVal(4.75); y->setVal(16.);
    cout <<"\nmassesEffPdf->getVal(" <<x->getVal() <<", " <<y->getVal() <<") = " <<massesEffPdf->getVal() <<endl;
    cout <<"with relEffPDF->haveUnitNorm() = " <<relEffPDF->haveUnitNorm() <<endl;
    if (relEffPDF->haveUnitNorm()) {
      relEffPDF->setUnitNorm(kFALSE);
      cout <<"\nmassesEffPdf->getVal(" <<x->getVal() <<"," <<y->getVal() <<") = " <<massesEffPdf->getVal() <<endl;
      cout <<"with relEffPDF->haveUnitNorm() = " <<relEffPDF->haveUnitNorm() <<endl;
    } else {
      relEffPDF->setUnitNorm(kTRUE);
      cout <<"\nmassesEffPdf->getVal(" <<x->.getVal() <<"," <<y->getVal() <<") = " <<massesEffPdf->getVal() <<endl;
      cout <<"with relEffPDF->haveUnitNorm() = " <<relEffPDF->haveUnitNorm() <<endl;
    }
    */

    // Creating the efficiency as function of masses from the efficiency as function of squared masses, in order to allow the multiplication with sigPDF 
    if (DalitzEff) {
      effMasses* massesEffPdf_fromDalitz = new effMasses(TString::Format("%s_fromDalitz",massesEffPdf->GetName()), massesEffPdf->GetTitle(), massKPi, massPsiPi, &mass2KPi, &mass2PsiPi, massesEffPdf);
      //TH1F* massesEffPdf_fromDalitz_mKPTH1 = (TH1F*)massesEffPdf_fromDalitz->createHistogram("massesEffPdf_fromDalitz_mKPTH1", massKPi) ; // not working with message: p.d.f normalization integral is zero or negative
      TH2F* massesEffPdf_fromDalitz_TH2 = (TH2F*)massesEffPdf_fromDalitz->createHistogram("", massKPi, Binning(1*xBins,massKPi_min,massKPi_max), YVar(massPsiPi, Binning(1*yBins,massPsiPi_min,massPsiPi_max)) ) ; massesEffPdf_fromDalitz_TH2->SetName("massesEffPdf_fromDalitz_TH2"); massesEffPdf_fromDalitz_TH2->SetTitleOffset(TH2_offset,"XY");
      massesEffPdf_fromDalitz_TH2->Draw("lego");
      gPad->SaveAs(TString::Format("%s/%s%s",dir.Data(),massesEffPdf_fromDalitz_TH2->GetName(),extension.Data()));
      modelWithEff = new RooProdPdf(modelName.Append("__withMassesEff"),TString::Format("(%s)*#epsilon(masses)",model->GetTitle()),RooArgSet(*model,*massesEffPdf_fromDalitz)) ;
    }

    model = modelWithEff;
  }
  //massKPi.setVal(2.); massPsiPi.setVal(4.);
  //cout <<"massKPi = " <<massKPi.getVal() <<", mass2KPi = " <<mass2KPi.getVal() <<", massPsiPi = " <<massPsiPi.getVal() <<", mass2PsiPi = " <<mass2PsiPi.getVal() <<endl;
  //cout <<"massesEffPdf->getVal() = " <<massesEffPdf->getVal() <<", massesEffPdf->getValV() = " <<massesEffPdf->getValV() <<endl;
  //cout <<"model->getVal() = " <<model->getVal() <<", model->getValV() = " <<model->getValV() <<endl;
  //cout <<"massKPi = " <<massKPi.getVal() <<", mass2KPi = " <<mass2KPi.getVal() <<", massPsiPi = " <<massPsiPi.getVal() <<", mass2PsiPi = " <<mass2PsiPi.getVal() <<endl;
  //return;


  // Angles efficiency
  RooAbsPdf* anglesEffPdf = 0;
  if (!effFile)
    effFile = TFile::Open(path+"officialMC_noPtEtaCuts_JPsi_Bd2MuMuKPi_2p0Sig_4p0to6p0SB.root");
  relEffName = "RelEff_planesAngle_vs_cos_psi2S_helicityAngle";
  relEffTH2 = (TH2F*)effFile->Get(relEffName) ;
  //relEffTH2 = 0;
  if (!relEffTH2) {
    cout <<"WARNING! No TH2F \"" <<relEffName <<"\" found in TFile \"" <<effFile->GetName() <<"\".\nExiting" <<endl;
  } else {
    Float_t xMin = relEffTH2->GetXaxis()->GetXmin(), xMax = relEffTH2->GetXaxis()->GetXmax(); Int_t xBins = relEffTH2->GetNbinsX();
    Float_t yMin = relEffTH2->GetYaxis()->GetXmin(), yMax = relEffTH2->GetYaxis()->GetXmax(); Int_t yBins = relEffTH2->GetNbinsY();
    
    // with RooHistPDF
    cosMuMu.setRange(xMin,xMax); phi.setRange(yMin,yMax);
    RooDataHist* relEffHist = new RooDataHist(relEffTH2->GetName(), relEffTH2->GetTitle(), angleVars, relEffTH2) ;
    /*
    for (Int_t iBin=0; iBin<relEffHist->numEntries(); ++iBin) { // remove bins whose center is out of Dalitz border
      angleVars = *(relEffHist->get( iBin ));
      //cout <<"x = " <<mass2KPi.getVal() <<", y = " <<mass2PsiPi.getVal() <<":\noriginal value = " <<relEffHist->weight(angleVars) <<", corrected value = " <<Dalitz_contour_host(mass2KPi.getVal(), mass2PsiPi.getVal(), kTRUE, psi_nS.Atoi()) * relEffHist->weight(angleVars) <<endl;
      if (!Dalitz_contour_host(x->getVal(), y->getVal(), DalitzEff, psi_nS.Atoi()))
	relEffHist->set(angleVars, 0, 0);
    }
    */
    TString effName = "relEffAngles";
    RooHistPdf* relEffPDF = new RooHistPdf(effName+"PDF","relative #epsilon(angles) pdf", angleVars, *relEffHist, 0) ; //If last argument is zero, the weight for the bin enclosing the coordinates contained in 'bin' is returned. For higher values, the result is interpolated in the real dimensions of the dataset with an order of interpolation equal to the value provided (more than 10 does not work)
    relEffPDF->setUnitNorm(kTRUE);
    Angles_contour* kinematicCheck_forEff = new Angles_contour("kinematicCheck_forEff","kinematic check", cosMuMu, phi) ;
    RooProdPdf* effWithKinCheck = new RooProdPdf(TString::Format("%s_withKinCheck",relEffPDF->GetName()),TString::Format("%s with kinematic check",relEffPDF->GetTitle()),RooArgSet(*kinematicCheck_forEff,*relEffPDF)) ;

    anglesEffPdf = relEffPDF;
    anglesEffPdf = effWithKinCheck;

    RooProdPdf* modelWithEff = new RooProdPdf(modelName.Append("__withAnglesEff"),TString::Format("(%s)*#epsilon(angles)",model->GetTitle()),RooArgSet(*model,*anglesEffPdf)) ;

    //relEffTH2->Draw("lego");
    TH2F* relEffTH2_fromHist = (TH2F*)relEffHist->createHistogram(effName+"TH2_fromHist", cosMuMu, Binning(xBins,xMin,xMax), YVar(phi, Binning(yBins,yMin,yMax)) ) ; relEffTH2_fromHist->SetTitle(relEffTH2->GetTitle()); relEffTH2_fromHist->SetTitleOffset(TH2_offset,"XY");
    TH2F* anglesEffPdf_TH2 = (TH2F*)anglesEffPdf->createHistogram("", cosMuMu, Binning(1*xBins,xMin,xMax), YVar(phi, Binning(1*yBins,yMin,yMax)) ) ;
    anglesEffPdf_TH2->SetName(TString::Format("%s_TH2",relEffPDF->GetName())); anglesEffPdf_TH2->SetTitle(relEffTH2->GetTitle()); anglesEffPdf_TH2->SetTitleOffset(TH2_offset,"XY");
    
    anglesEffPdf_TH2->SetLineColor(kRed);
    anglesEffPdf_TH2->Draw("surf2");
    //anglesEffPdf_TH2->Draw("lego");
    relEffTH2_fromHist->Draw("lego same");
    gPad->SaveAs(TString::Format("%s/%s%s",dir.Data(),anglesEffPdf_TH2->GetName(),extension.Data()));

    // Projection X
    //TH1D* relEff_cosMuMu = relEffTH2->ProjectionX("relEff_cosMuMu",1,xBins);
    TH1D* relEff_cosMuMu = relEffTH2_fromHist->ProjectionX("relEff_cosMuMu"); relEff_cosMuMu->SetTitle(TString::Format("Projection of %s",relEffTH2_fromHist->GetTitle()));
    relEff_cosMuMu->Draw();
    //
    //TH1F* anglesEffPdf_cosMuMuHist = (TH1F*)anglesEffPdf->createHistogram("anglesEffPdf_cosMuMuHist", x, Binning(1*xBins,xMin,xMax) ) ; // not working with message: function value is NAN
    TH1D* anglesEffPdf_cosMuMuHist = anglesEffPdf_TH2->ProjectionX(); //"anglesEffPdf_cosMuMuHist",1,xBins
    anglesEffPdf_cosMuMuHist->SetMinimum(0);
    anglesEffPdf_cosMuMuHist->Scale(relEff_cosMuMu->Integral()); anglesEffPdf_cosMuMuHist->SetLineColor(kRed); anglesEffPdf_cosMuMuHist->Draw("same");
    gPad->SaveAs(TString::Format("%s/%s%s",dir.Data(),relEff_cosMuMu->GetName(),extension.Data()));
    //return;
    // Projection Y
    //TH1D* relEff_phi = relEffTH2->ProjectionX("relEff_phi",1,xBins);
    TH1D* relEff_phi = relEffTH2_fromHist->ProjectionY("relEff_phi"); relEff_phi->SetTitle(TString::Format("Projection of %s",relEffTH2_fromHist->GetTitle()));
    relEff_phi->Draw();
    //
    //TH1F* anglesEffPdf_phiHist = (TH1F*)anglesEffPdf->createHistogram("anglesEffPdf_phiHist", y, Binning(1*yBins,yMin,yMax) ) ; anglesEffPdf_mPsiPiHist->SetLineColor(kRed); anglesEffPdf_mPsiPiHist->Draw();
    TH1D* anglesEffPdf_phiHist = anglesEffPdf_TH2->ProjectionY(); //"anglesEffPdf_phiHist",1,yBins
    anglesEffPdf_phiHist->SetMinimum(0);
    anglesEffPdf_phiHist->Scale(relEff_phi->Integral()); anglesEffPdf_phiHist->SetLineColor(kRed); anglesEffPdf_phiHist->Draw("same");
    gPad->SaveAs(TString::Format("%s/%s%s",dir.Data(),relEff_phi->GetName(),extension.Data()));
    //return;

    model = modelWithEff;
  }
  //return;

  RooPlot* massKP_frame = massKPi.frame() ; massKP_frame->SetTitle("Projection of "+massKPi_title);
  vector <RooPlot*> var_frame;
  for (Int_t iVar=0; iVar<kinematicVars.getSize(); ++iVar) {
    var_frame.push_back( ((RooRealVar&)kinematicVars[sigPDF_varNameTitle[iVar].first]).frame(Title("Projection of "+sigPDF_varNameTitle[iVar].second)) );
  }

  Int_t nLegendEntries = 0;

  // Generate toy data from pdf and plot data and p.d.f on frame
  cout <<"\nGenerating " <<nEvents.getVal() <<" events according to " <<model->GetTitle() <<" pdf for " <<model->GetName() <<endl;
  timeval genTime;
  gettimeofday(&start, NULL);
  startCPU = times(&startProc);
  //
  RooDataSet* dataGenPDF = model->generate(kinematicVars, nEvents.getVal(), Verbose(kTRUE), Name("Generated_data_from_PDF")) ; dataGenPDF->SetTitle("Generated data from PDF");
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
  if (nEvents.getVal() > 100000)
    dataGenPDF->write(TString::Format("%s/%s__EvtGen.txt",datasetsPath.Data(),model->GetName()));
  //return;

  //dataToFit = dataGenPDF;

  // Create corresponding m2 dataset + K* helicity angle if absent
  if ( !kinematicVars_m2.contains(cosKstar) )
    kinematicVars_m2.add(cosKstar);
  RooDataSet* dataGenPDF_withCos = (RooDataSet*)dataGenPDF->emptyClone(); 	
  dataGenPDF_withCos->addColumn(cosKstar);

  RooDataSet* dataGenPDF_m2 = new RooDataSet(dataGenPDF->GetName(), dataGenPDF->GetTitle(), kinematicVars_m2 ) ;
  for (Int_t iEvent=0; iEvent<dataGenPDF->numEntries(); ++iEvent) {
    kinematicVars = *(dataGenPDF->get(iEvent)) ; // this line will propagate the RooRealVars values of the event to all the corresponding RooRealVars in kinematicVars  
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
    //cout <<"\nmassesEffPdf->printValue(cout) = "; massesEffPdf->printValue(cout); cout <<endl;
    
    // set K* helicity angle value
    cosKstar.setVal( cosTheta_FromMasses_host(mass2KPiFor.getVal(), mass2PsiPiFor.getVal(), massMuMu, MBd2, MKaon2, MPion2) );
    
    dataGenPDF_m2->add( kinematicVars_m2 ) ;
    dataGenPDF_withCos->add( RooArgSet(kinematicVars,cosKstar) );     
  }
  if (dataGenPDF->numEntries() != dataGenPDF_m2->numEntries()) {
    cout <<"dataGenPDF->numEntries() (" <<dataGenPDF->numEntries() <<") != dataGenPDF_m2->numEntries() (" <<dataGenPDF_m2->numEntries() <<")! Please check" <<endl; return;
  }
  if (dataGenPDF->numEntries() != dataGenPDF_withCos->numEntries()) {
    cout <<"dataGenPDF->numEntries() (" <<dataGenPDF->numEntries() <<") != dataGenPDF_withCos->numEntries() (" <<dataGenPDF_withCos->numEntries() <<")! Please check" <<endl; return;
  }
  //return;
       
  cout <<"\nPlotting data..." <<endl;
  TString plotName = model->GetName();

  Float_t rightMargin = 0.12;
  cout <<"\nPlotting angles scatter plot..." <<endl;
  TCanvas* scatter_C = new TCanvas("Angles_scatter_plot","Angles scatter plot",800,600) ; scatter_C->SetRightMargin(rightMargin);
  scatter_C->cd();
  Float_t cos_limit = 1.02; Float_t cos_margin = 0.02;
  Float_t phi_limit = 3.2; Float_t phi_margin = 0.05;
  if ( totEvents < 50000 ) {
    phi_limit = 3.3; phi_margin = 0.1;
    cos_limit = 1.025; cos_margin = 0.025;
  }
  Int_t cos_bins = 2*cos_limit/cos_margin;
  Int_t phi_bins = 2*phi_limit/phi_margin;
  TH2F* scatter = (TH2F*)dataGenPDF_m2->createHistogram("Angles_scatter_plot", cosMuMu, Binning(cos_bins,-cos_limit,cos_limit), YVar(phi, Binning(phi_bins,-phi_limit,phi_limit)) ) ; scatter->SetTitle( TString::Format("Angles scatter plot;%s;%s",cosMuMu.GetTitle(),phi.GetTitle()) ) ;
  gStyle->SetOptStat( 10 ) ;
  scatter->Draw("colz"); 	
  scatter_C->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),scatter_C->GetName(),plotName.Data(),extension.Data()));

  cout <<"\nPlotting Dalitz..." <<endl;
  Float_t KPiMass2_low = 0., KPiMass2_high = 5.; Int_t KPiMass2_bins = 100;
  Float_t MuMuPiMass2_low = 9., MuMuPiMass2_high = 25.; Int_t MuMuPiMass2_bins = 128;

  TCanvas* dalitz_C = new TCanvas("Dalitz_C","Dalitz",800,600) ; dalitz_C->SetRightMargin(rightMargin); dalitz_C->cd();
  if ( totEvents < 20000 ) {KPiMass2_bins = 50; MuMuPiMass2_bins = 64;}
  TH2F* dalitz = (TH2F*)dataGenPDF_m2->createHistogram("Dalitz", mass2KPi, Binning(KPiMass2_bins,KPiMass2_low,KPiMass2_high), YVar(mass2PsiPi, Binning(MuMuPiMass2_bins,MuMuPiMass2_low,MuMuPiMass2_high)) ) ; dalitz->SetTitle( TString::Format("Dalitz;%s;%s",mass2KPi.GetTitle(),mass2PsiPi.GetTitle()) ) ;
  gStyle->SetOptStat( 10 ) ;
  dalitz->Draw("colz"); 	
  //dalitz->Draw("lego"); 	
  dalitz_C->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),dalitz_C->GetTitle(),plotName.Data(),extension.Data()));
  //return;
  cout <<"\nPlotting rectangular Dalitz..." <<endl;
  TCanvas* dalitzRect_C = new TCanvas("DalitzRect_C","Rectangular_Dalitz",800,600) ; dalitzRect_C->SetRightMargin(rightMargin);
  dalitzRect_C->cd();
  TH2F* dalitz_rect = (TH2F*)dataGenPDF_withCos->createHistogram("DalitzRect", massKPi, YVar(cosKstar, Binning(cos_bins,-cos_limit,cos_limit)) ) ; dalitz_rect->SetTitle( TString::Format("Rectangular Dalitz;%s;%s",massKPi.GetTitle(),cosKstar.GetTitle()) ) ;
  gStyle->SetOptStat( 10 ) ;
  dalitz_rect->Draw("colz"); 	
  dalitzRect_C->SaveAs(TString::Format("%s/%s_%s%s",dir.Data(),dalitzRect_C->GetTitle(),plotName.Data(),extension.Data()));
 
  cout <<"\nPlotting m(psiPi)..." <<endl;
  dataGenPDF->plotOn( var_frame[2] ) ;
  TCanvas* massPsiP_C = new TCanvas( TString::Format("%s_C",massPsiPi.GetName()), massPsiPi.GetName(), 800,600) ;
  massPsiP_C->cd();
  var_frame[2]->Draw() ;
  massPsiP_C->SaveAs(TString::Format("%s/%s__MuMuPi%s",dir.Data(),plotName.Data(),extension.Data()));
  return;
  Int_t fullModelColor = 2; // 2 = kRed
  Int_t bkgColor = fullModelColor;
  TString modelEntry = "Full model";
  
  cout <<"\nPlotting m(KPi)..." <<endl;
  dataGenPDF->plotOn(var_frame[0]) ; nLegendEntries++;
  //
  Bool_t fitting = kFALSE; //fitting = kTRUE;
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
  leg->AddEntry(dataGenPDF,"","ep");
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
    RooAddition nllp("nllp","nllp",toMinimize);
    RooMinuit m(nllp); m.setVerbose();
    
    timeval fitModelTime;
    gettimeofday(&start, NULL);
    startCPU = times(&startProc);
    //
    fitres = model->fitTo(*dataGenPDF, Hesse(kFALSE), Minos(kFALSE), Save(kTRUE), NumCPU(nCPU), Verbose(kTRUE), PrintLevel(3)) ;
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

