#define psiPrimePiK_MC_cxx
// The class definition in psiPrimePiK_MC.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("psiPrimePiK_MC.C")
// Root > T->Process("psiPrimePiK_MC.C","some options")
// Root > T->Process("psiPrimePiK_MC.C+")
//

#include "psiPrimePiK_MC.h"
#include <TStyle.h>
#include <TPad.h>
#include <TCanvas.h>
// RooFit
//gInterpreter->AddIncludePath("/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms5/include") ;
//gROOT->ProcessLine(".include /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms5/include") ;
// one of the 2 above lines need to be executed interactively before following includes:
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooConstVar.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooAbsReal.h>
#include <RooPlot.h>

#include <TColor.h>
#include <TLine.h>

using namespace RooFit; // needed even without .C+


void psiPrimePiK_MC::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   if (option == "MC" || option == "officialMC") {
     MC = kTRUE ;
     if (option == "officialMC")
       officialMC = kTRUE ;
   } else {
     MC = kFALSE ; 
     officialMC = kFALSE ;  
   }
   // need to reinitialize variables used in Terminate
   muon_mass = 0.1056583715;
   // phi_mass = 1.019455;
   jpsi_mass = 3.0967;
   psi2S_mass = 3.686109 ;
   pionCh_mass = 0.13957018;
   kaonCh_mass = 0.493677 ; 
   Kstar0_mass = 0.8961; // from evt.pdl
   // Bplus_mass = 5.27925;
   B0_mass = 5.27958;
   //
   phi_left = 1.01; phi_right = 1.03;
   Bs_left = 5.3; Bs_right = 5.43;
}

void psiPrimePiK_MC::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   if (option == "MC" || option == "officialMC") {
     MC = kTRUE ;
     if (option == "officialMC")
       officialMC = kTRUE ;
   } else {
     MC = kFALSE ; 
     officialMC = kFALSE ;  
   }
   /////////////////////////////////////////////////////////////////

   newSoftMuonID = kFALSE ; //newSoftMuonID = kTRUE ;
   only1psiPrime = kFALSE ; //  only1psiPrime = kTRUE ;
   Kai_cuts = kFALSE ;      //  Kai_cuts = kTRUE ; // Kai stands for 2011+2012
   HLT_Displaced = kFALSE ; // HLT_Displaced = kTRUE ;
   fixedMassWind = kFALSE ; // fixedMassWind = kTRUE ;
   sideband = kFALSE ;      // sideband = kTRUE ;

   _nt = new TNtupleD("AAVars","AA variables","runNum:evtNum:B0_mass:KPi_mass:Jpsi_mass:JpsiPi_mass:B0_3mom:KPi_3mom:Jpsi_3mom:K_3mom:Pi_3mom:Theta_Kstar:Theta_Jpsi:Phi"); // output ntuple
   if ( MC ) {
     if ( !officialMC )
       OutFile = new TProofOutputFile( "MC_Bd2Psi2SKpi.root" );
     else
       OutFile = new TProofOutputFile( "officialMC_Bd2Psi2SKpi.root" );
   } else {
     if (!HLT_Displaced) {
       if (!Kai_cuts) { // kai stands for 2011+2012
	 if (!newSoftMuonID) {
	   if (!only1psiPrime) {
	     if (!fixedMassWind)
	       OutFile = new TProofOutputFile( "HLT5.root" );
	     else if (fixedMassWind)
	       OutFile = new TProofOutputFile( "HLT5_fixedMassWind.root" );
	   } else if (only1psiPrime) {
	     if (!fixedMassWind)
	       OutFile = new TProofOutputFile( "HLT5_onePsiPrime.root" );
	     else if (fixedMassWind)
	       OutFile = new TProofOutputFile( "HLT5_onePsiPrime_fixedMassWind.root" );
	   }
	 } else if (!only1psiPrime) {
	   //OutFile = new TProofOutputFile( "HLT5_newSoftMuon.root" );
	   OutFile = new TProofOutputFile( "HLT5_newSoftMuon_alsoInPsiPrimeWind.root" );
	 } else {
	   //OutFile = new TProofOutputFile( "HLT5_only1psiPrime.root" );
	   //OutFile = new TProofOutputFile( "HLT5_onePsiPrime.root" );
	   OutFile = new TProofOutputFile( "HLT5_onePsiPrime_check.root" );
	 }
       } else { // include newSoftMuonID  // kai stands for 2011+2012
	 if (!fixedMassWind) {
	   if (!sideband) {
	     //OutFile = new TProofOutputFile( "2011.root" );
	     OutFile = new TProofOutputFile( "kai_onePsiPrime.root" ); // kai stands for 2011+2012
	   } else {
	     OutFile = new TProofOutputFile( "kai_onePsiPrime_sb.root" ); } // kai stands for 2011+2012
	 } else {
	   OutFile = new TProofOutputFile( "fixedMassWind.root" ); }
	 //OutFile = new TProofOutputFile( "fixedMassWind_test.root" ); }
       }
     } else
       if (!fixedMassWind)
	 OutFile = new TProofOutputFile( "hlt3p5_displaced.root" );
       else
	 OutFile = new TProofOutputFile( "hlt3p5_displaced_fixedMassWind.root" );
   }

   /*
  fOut = OutFile->OpenFile("RECREATE");
  SaraOut = SaraFile->OpenFile("RECREATE");
  if ( MC )
    AlexisOut = AlexisFile->OpenFile("RECREATE");
   */

  // Open the files

  if ( !(fOut = OutFile->OpenFile("RECREATE")) )
    Warning("SlaveBegin", "problems opening file: %s/%s", OutFile->GetDir(), OutFile->GetFileName());

  //SaraFile = new TProofOutputFile( "Sara_histos.root" );
  //if ( !(SaraOut = SaraFile->OpenFile("RECREATE")) )
    //Warning("SlaveBegin", "problems opening file: %s/%s", SaraFile->GetDir(), SaraFile->GetFileName());
  SaraOut = 0 ; // disable Sara's histos

  if ( MC ) {
     if ( !officialMC )
       AlexisFile = new TProofOutputFile( "26Jan_histos.root" );
     else
       AlexisFile = new TProofOutputFile( "26Jan_histos.root" );
     //
     if ( !(AlexisOut = AlexisFile->OpenFile("RECREATE")) ) 
       Warning("SlaveBegin", "problems opening file: %s/%s", AlexisFile->GetDir(), AlexisFile->GetFileName());
  } else
    AlexisOut = 0 ;
    
  if ( !officialMC )
    tightCutsFile = new TProofOutputFile( "TightCuts_histos.root" );
  else
    tightCutsFile = new TProofOutputFile( "TightCuts_histos.root" );
  //
  if ( !(tightCutsOut = tightCutsFile->OpenFile("RECREATE")) )
    Warning("SlaveBegin", "problems opening file: %s/%s", tightCutsFile->GetDir(), tightCutsFile->GetFileName());

  //////// HISTOGRAMS
  // put them here:
  //
  Float_t MuMuMass_low = 3.4, MuMuMass_high = 4. ;  Int_t MuMuMass_bins = 120 ;
  Float_t MuMuPiMass_low = 3.8, MuMuPiMass_high = 5.5 ;  Int_t MuMuPiMass_bins = 136 ;
  Float_t MuMuKPiMass_low = 5.15, MuMuKPiMass_high = 5.45 ;  Int_t MuMuKPiMass_bins = 100 ;
  myInclusiveMuMuMass = new TH1F("myInclusiveMuMuMass","myInclusiveMuMuMass;m(#mu^{+}#mu^{-}) [GeV]; Candidates / 5 MeV",MuMuMass_bins,MuMuMass_low,MuMuMass_high); // 200,2.5,4.5
  myInclusiveMuMuMass_noTrigg = new TH1F("myInclusiveMuMuMass_noTrigg","myInclusiveMuMuMass noTrigg",400,2.1,4.1);
  //
  // HLT_0_v*
  myInclusiveMuMuMass_HLT_0_v3 = new TH1F("myInclusiveMuMuMass_HLT_0_v3","myInclusiveMuMuMass_HLT_0_v3",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v3A = new TH1F("myInclusiveMuMuMass_HLT_0_v3A","myInclusiveMuMuMass_HLT_0_v3A",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v3B = new TH1F("myInclusiveMuMuMass_HLT_0_v3B","myInclusiveMuMuMass_HLT_0_v3B",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v4 = new TH1F("myInclusiveMuMuMass_HLT_0_v4","myInclusiveMuMuMass_HLT_0_v4",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v3Bv4 = new TH1F("myInclusiveMuMuMass_HLT_0_v3Bv4","myInclusiveMuMuMass_HLT_0_v3Bv4",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v5 = new TH1F("myInclusiveMuMuMass_HLT_0_v5","myInclusiveMuMuMass_HLT_0_v5",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v6 = new TH1F("myInclusiveMuMuMass_HLT_0_v6","myInclusiveMuMuMass_HLT_0_v6",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v6C = new TH1F("myInclusiveMuMuMass_HLT_0_v6C","myInclusiveMuMuMass_HLT_0_v6C",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v5v6C = new TH1F("myInclusiveMuMuMass_HLT_0_v5v6C","myInclusiveMuMuMass_HLT_0_v5v6C",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_0_v6D = new TH1F("myInclusiveMuMuMass_HLT_0_v6D","myInclusiveMuMuMass_HLT_0_v6D",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  // HLT_5_v*
  myInclusiveMuMuMass_HLT_5_v3 = new TH1F("myInclusiveMuMuMass_HLT_5_v3","myInclusiveMuMuMass_HLT_5_v3",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v3A = new TH1F("myInclusiveMuMuMass_HLT_5_v3A","myInclusiveMuMuMass_HLT_5_v3A",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v3B = new TH1F("myInclusiveMuMuMass_HLT_5_v3B","myInclusiveMuMuMass_HLT_5_v3B",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v4 = new TH1F("myInclusiveMuMuMass_HLT_5_v4","myInclusiveMuMuMass_HLT_5_v4",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v3Bv4 = new TH1F("myInclusiveMuMuMass_HLT_5_v3Bv4","myInclusiveMuMuMass_HLT_5_v3Bv4",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v5 = new TH1F("myInclusiveMuMuMass_HLT_5_v5","myInclusiveMuMuMass_HLT_5_v5",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v6 = new TH1F("myInclusiveMuMuMass_HLT_5_v6","myInclusiveMuMuMass_HLT_5_v6",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v6C = new TH1F("myInclusiveMuMuMass_HLT_5_v6C","myInclusiveMuMuMass_HLT_5_v6C",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v5v6C = new TH1F("myInclusiveMuMuMass_HLT_5_v5v6C","myInclusiveMuMuMass_HLT_5_v5v6C",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_5_v6D = new TH1F("myInclusiveMuMuMass_HLT_5_v6D","myInclusiveMuMuMass_HLT_5_v6D",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  // HLT_7_v*
  myInclusiveMuMuMass_HLT_7_v1 = new TH1F("myInclusiveMuMuMass_HLT_7_v1","myInclusiveMuMuMass_HLT_7_v1",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_7_v2 = new TH1F("myInclusiveMuMuMass_HLT_7_v2","myInclusiveMuMuMass_HLT_7_v2",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_7_v3 = new TH1F("myInclusiveMuMuMass_HLT_7_v3","myInclusiveMuMuMass_HLT_7_v3",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_7_v3C = new TH1F("myInclusiveMuMuMass_HLT_7_v3C","myInclusiveMuMuMass_HLT_7_v3C",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_7_v2v3C = new TH1F("myInclusiveMuMuMass_HLT_7_v2v3C","myInclusiveMuMuMass_HLT_7_v2v3C",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_HLT_7_v3D = new TH1F("myInclusiveMuMuMass_HLT_7_v3D","myInclusiveMuMuMass_HLT_7_v3D",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  // datasets
  myInclusiveMuMuMass_12A = new TH1F("myInclusiveMuMuMass_12A","myInclusiveMuMuMass_12A",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12A_HLT_0_v3 = new TH1F("myInclusiveMuMuMass_12A_HLT_0_v3","myInclusiveMuMuMass_12A_HLT_0_v3",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12A_HLT_5_v3 = new TH1F("myInclusiveMuMuMass_12A_HLT_5_v3","myInclusiveMuMuMass_12A_HLT_5_v3",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12A_HLT_9_v9 = new TH1F("myInclusiveMuMuMass_12A_HLT_9_v9","myInclusiveMuMuMass_12A_HLT_9_v9",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  //
  myInclusiveMuMuMass_12B = new TH1F("myInclusiveMuMuMass_12B","myInclusiveMuMuMass_12B",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12B_HLT_0_v3v4 = new TH1F("myInclusiveMuMuMass_12B_HLT_0_v3v4","myInclusiveMuMuMass_12B_HLT_0_v3v4",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12B_HLT_5_v3v4 = new TH1F("myInclusiveMuMuMass_12B_HLT_5_v3v4","myInclusiveMuMuMass_12B_HLT_5_v3v4",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12B_HLT_7_v1_9_v9 = new TH1F("myInclusiveMuMuMass_12B_HLT_7_v1_9_v9","myInclusiveMuMuMass_12B_HLT_7_v1_9_v9",MuMuMass_bins,MuMuMass_low,MuMuMass_high);

  myInclusiveMuMuMass_12C = new TH1F("myInclusiveMuMuMass_12C","myInclusiveMuMuMass_12C",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12C_HLT_0_v5 = new TH1F("myInclusiveMuMuMass_12C_HLT_0_v5","myInclusiveMuMuMass_12C_HLT_0_v5",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12C_HLT_5_v5 = new TH1F("myInclusiveMuMuMass_12C_HLT_5_v5","myInclusiveMuMuMass_12C_HLT_5_v5",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12C_HLT_7_v2 = new TH1F("myInclusiveMuMuMass_12C_HLT_7_v2","myInclusiveMuMuMass_12C_HLT_7_v2",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12C_HLT_5_v5_or_7_v2 = new TH1F("myInclusiveMuMuMass_12C_HLT_5_v5_or_7_v2","myInclusiveMuMuMass_12C_HLT_5_v5_or_7_v2",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  //
  myInclusiveMuMuMass_12D = new TH1F("myInclusiveMuMuMass_12D","myInclusiveMuMuMass_12D",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12D_HLT_0_v6 = new TH1F("myInclusiveMuMuMass_12D_HLT_0_v6","myInclusiveMuMuMass_12D_HLT_0_v6",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12D_HLT_5_v6 = new TH1F("myInclusiveMuMuMass_12D_HLT_5_v6","myInclusiveMuMuMass_12D_HLT_5_v6",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  myInclusiveMuMuMass_12D_HLT_7_v3 = new TH1F("myInclusiveMuMuMass_12D_HLT_7_v3","myInclusiveMuMuMass_12D_HLT_7_v3",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  //
  myInclusiveMuMuMass_out = new TH1F("myInclusiveMuMuMass_out","myInclusiveMuMuMass_out",MuMuMass_bins,MuMuMass_low,MuMuMass_high);
  //
  myMuMuMass = new TH1F("myMuMuMass", "myMuMuMass;m(#mu^{+}#mu^{-}) [GeV]", 200, 2.5, 4.5);
  //myMuMuMass_test = new TH1F("myMuMuMass_test", "myMuMuMass_test;m(#mu^{+}#mu^{-}) [GeV]", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  myMuMuMass_noMuID = new TH1F("myMuMuMass_noMuID", "myMuMuMass no MuID;m(#mu^{+}#mu^{-}) [GeV]", MuMuMass_bins, MuMuMass_low, MuMuMass_high); 
  myMuMuMass_Mu1ID = new TH1F("myMuMuMass_Mu1ID", "myMuMuMass Mu1ID;m(#mu^{+}#mu^{-}) [GeV]", MuMuMass_bins, MuMuMass_low, MuMuMass_high); 
  myMuMuMass_Mu2ID = new TH1F("myMuMuMass_Mu2ID", "myMuMuMass Mu2ID;m(#mu^{+}#mu^{-}) [GeV]", MuMuMass_bins, MuMuMass_low, MuMuMass_high); 
  myMuMuMass_MuID = new TH1F("myMuMuMass_MuID", "myMuMuMass MuID;m(#mu^{+}#mu^{-}) [GeV]", MuMuMass_bins, MuMuMass_low, MuMuMass_high); 

  myMuMuMass_HLT = new TH1F("myMuMuMass_HLT", "myMuMuMass HLT;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high); // 200, 2.5, 4.5
  myMuMuMass_HLT5_7comp = new TH1F("myMuMuMass_HLT5_7comp", "myMuMuMass HLT5 7comp;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  myMuMuMass_HLT3p5 = new TH1F("myMuMuMass_HLT3p5", "myMuMuMass HLT 3.5 displaced;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  //myMuMuMass_HLT_yMu = new TH1F("myMuMuMass_HLT_yMu", "myMuMuMass HLT for y(#mu^{#pm}) > 2.4;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  //myMuMuMass_HLT_yMu_ypsi = new TH1F("psi2S_HLT_yMu_ypsi", "myMuMuMass HLT for y(#mu^{#pm}, #psi') > 2.4;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  myMuMuMass_HLT_final = new TH1F("myMuMuMass_HLT_final", "myMuMuMass HLT all cuts;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  myMuMuMass_HLT_max2psiPrime = new TH1F("myMuMuMass_HLT_max2psiPrime", "myMuMuMass HLT all cuts (##psi' #leq 2);m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high); // 400, 2.5, 4.5
  myMuMuMass_HLT_2psiPrime = new TH1F("myMuMuMass_HLT_2psiPrime", "myMuMuMass HLT all cuts (##psi' = 2);m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  myMuMuMass_HLT_1psiPrime = new TH1F("myMuMuMass_HLT_1psiPrime", "myMuMuMass HLT all cuts (##psi' = 1);m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  myMuMuMass_HLT_VtxCLDiff_cut = new TH1F("myMuMuMass_HLT_VtxCLDiff_cut", "myMuMuMass HLT all cuts (|1^{st} - 2^{nd}| #psi' vertex CL < 0.3);m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  myMuMuMass_HLT_DR1 = new TH1F("myMuMuMass_HLT_DR1", "myMuMuMass HLT for #DeltaR(#mu^{+}, #mu^{-})<1;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  myMuMuMass_HLT7 = new TH1F("myMuMuMass_HLT7", "myMuMuMass HLT7;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high); // 200, 2.5, 4.5
  //
  myMuMuKMass = new TH1F("myMuMuKMass", "myMuMuKMass;m(#mu^{+}#mu^{-}K^{+}) [GeV]", 200, 0.0, 10.0);
  // Z histos
  myMuMuPiMass = new TH1F("myMuMuPiMass", "myMuMuPiMass;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high); // 128, 3.9, 5.5
  myMuMuPiMass_pr = new TH1F("myMuMuPiMass_pr", "myMuMuPiMass prompt;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_pr_DzCut = new TH1F("myMuMuPiMass_pr_DzCut", "myMuMuPiMass prompt d_{z}(#pi^{+})<0.4;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_pr_Dz01 = new TH1F("myMuMuPiMass_pr_Dz01", "myMuMuPiMass prompt d_{z}(#pi^{+})<1;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_quasipr = new TH1F("myMuMuPiMass_quasipr", "myMuMuPiMass quasi prompt;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_nonpr = new TH1F("myMuMuPiMass_nonpr", "myMuMuPiMass non prompt;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_DeltaR0p1_nonpr5 = new TH1F("myMuMuPiMass_DeltaR0p1_nonpr5", "myMuMuPiMass non prompt 5, #DeltaR(#psi',#pi^{+}) < 0.1;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_DeltaR0p2_nonpr5 = new TH1F("myMuMuPiMass_DeltaR0p2_nonpr5", "myMuMuPiMass non prompt 5, #DeltaR(#psi',#pi^{+}) < 0.2;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_DeltaR0p1_nonpr8 = new TH1F("myMuMuPiMass_DeltaR0p1_nonpr8", "myMuMuPiMass non prompt 8, #DeltaR(#psi',#pi^{+}) < 0.1;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_DeltaR0p2_nonpr8 = new TH1F("myMuMuPiMass_DeltaR0p2_nonpr8", "myMuMuPiMass non prompt 8, #DeltaR(#psi',#pi^{+}) < 0.2;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_nonpr_Dz1inf = new TH1F("myMuMuPiMass_nonpr_Dz1inf", "myMuMuPiMass non prompt d_{z}(#pi^{+})>1;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_barrel = new TH1F("myMuMuPiMass_barrel", "myMuMuPiMass for #psi' in barrel;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_barrel_pr = new TH1F("myMuMuPiMass_barrel_pr", "myMuMuPiMass prompt for #psi' in barrel;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_nonBarrel_pr = new TH1F("myMuMuPiMass_nonBarrel_pr", "myMuMuPiMass prompt for #psi' not in barrel;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_barrel_nonpr = new TH1F("myMuMuPiMass_barrel_nonpr", "myMuMuPiMass non prompt for #psi' in barrel;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass->SetLineWidth(2) ; myMuMuPiMass_pr->SetLineWidth(2) ; myMuMuPiMass_quasipr->SetLineWidth(2) ; myMuMuPiMass_nonpr->SetLineWidth(2) ; myMuMuPiMass_DeltaR0p1_nonpr5->SetLineWidth(2) ; myMuMuPiMass_DeltaR0p2_nonpr5->SetLineWidth(2) ; myMuMuPiMass_DeltaR0p1_nonpr8->SetLineWidth(2) ; myMuMuPiMass_DeltaR0p2_nonpr8->SetLineWidth(2) ;
  myMuMuPiMass_Kai = new TH1F("myMuMuPiMass_Kai", "myMuMuPiMass Kai cuts;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);  myMuMuPiMass_Kai->SetLineWidth(2) ;
  myMuMuPiMass_Kai_HLT7 = new TH1F("myMuMuPiMass_Kai_HLT7", "myMuMuPiMass Kai cuts HLT7;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_HLT5to3p5 = new TH1F("myMuMuPiMass_HLT5to3p5", "myMuMuPiMass HLT 5 to 3.5;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuPiMass_HLT3p5_allCuts = new TH1F("myMuMuPiMass_HLT3p5_allCuts", "myMuMuPiMass HLT 3.5 displaced all cuts reproduced;m(#mu^{+}#mu^{-}#pi^{+}) [GeV];Candidates / 5 MeV", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuKMass_HLT5to3p5 = new TH1F("myMuMuKMass_HLT5to3p5", "myMuMuKMass HLT 5 to 3.5;m(#mu^{+}#mu^{-}K^{+}) [GeV]", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
  myMuMuKMass_HLT3p5_allCuts = new TH1F("myMuMuKMass_HLT3p5_allCuts", "myMuMuKMass HLT 3.5 displaced all cuts reproduced;m(#mu^{+}#mu^{-}K^{+}) [GeV];Candidates / 5 MeV", MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);

  myMuMuKMass_Kai = new TH1F("myMuMuKMass_Kai", "myMuMuKMass Kai cuts;m(#mu^{+}#mu^{-}K^{+}) [GeV]", MuMuPiMass_bins, MuMuPiMass_low, MuMuPiMass_high);
  myMuMuMass_noMuID_KaiBF = new TH1F("myMuMuMass_noMuID_KaiBF", "myMuMuMass no MuID;m(#mu^{+}#mu^{-}) [GeV]", MuMuMass_bins, MuMuMass_low, MuMuMass_high); 
  myMuMuMass_MuID_KaiBF = new TH1F("myMuMuMass_MuID_KaiBF", "myMuMuMass MuID;m(#mu^{+}#mu^{-}) [GeV]", MuMuMass_bins, MuMuMass_low, MuMuMass_high); 
  myMuMuMass_KaiBF = new TH1F("myMuMuMass_KaiBF", "myMuMuMass Kai BF cuts;m(#mu^{+}#mu^{-}) [GeV]", MuMuMass_bins, MuMuMass_low, MuMuMass_high); 
  myMuMuKMass_KaiBF = new TH1F("myMuMuKMass_KaiBF", "myMuMuKMass Kai BF cuts;m(#mu^{+}#mu^{-}K^{+}) [GeV]", 520, 5.15, 5.41); 
  dEdx_vs_p_Bplus = new TH2F("dEdx_vs_p_KFromBplus", "dE/dx vs p for K^{+} from B^{+};p(K^{+}) [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;

  
  // Sara plots
  for (Int_t i = 0; i<=5; ++i)
    Z_CTauBSOverE.push_back( (Float_t)i ) ;
  for (Int_t i = 0; i<=20; ++i)
    mu_pT.push_back( i/2. ) ;
  for (Int_t i = 0; i<=15; ++i)
    trk_pT.push_back( 0.6 +i*0.2 ) ;
  for (Int_t i = 0; i<=10; ++i)
    Z_CL.push_back( i*0.001 ) ;
  for (Int_t i = 1; i<20; ++i)
    Z_CL.push_back( 0.01 +i*0.01 ) ;

  for (Int_t i = 0; i < abs(Z_CTauBSOverE.size()); ++i) {
    TString name = TString::Format("MuMuKMass_ZctauBS_over_E_%.0f", Z_CTauBSOverE[i]) ; name.ReplaceAll(".","p") ;
    TString title = TString::Format("MuMuKMass with c*#tau_{BS}(#mu#mu#pi)/#sigma > %.0f;m(#mu^{+}#mu^{-}K^{+}) [GeV]", Z_CTauBSOverE[i]) ; 
    MuMuKMass_ZCTauBSOverE.push_back( new TH1F(name, title, 130, 5.15, 5.41) ) ; 
  }
  for (Int_t i = 0; i < abs(mu_pT.size()); ++i) {
    TString name1 = TString::Format("MuMuKMass_mu1pT%.1f", mu_pT[i]) ; name1.ReplaceAll(".","p") ;
    TString title1 = TString::Format("MuMuKMass with p_{T}(#mu_{1}) > %.1f;m(#mu^{+}#mu^{-}K^{+}) [GeV]", mu_pT[i]) ; 
    MuMuKMass_mu1Pt.push_back( new TH1F(name1, title1, 130, 5.15, 5.41) ) ; 
    TString name2 = TString::Format("MuMuKMass_mu2pT%.1f", mu_pT[i]) ; name2.ReplaceAll(".","p") ;
    TString title2 = TString::Format("MuMuKMass with p_{T}(#mu_{2}) > %.1f;m(#mu^{+}#mu^{-}K^{+}) [GeV]", mu_pT[i]) ; 
    MuMuKMass_mu2Pt.push_back( new TH1F(name2, title2, 130, 5.15, 5.41) ) ; 
  }
  for (Int_t i = 0; i < abs(trk_pT.size()); ++i) {
    TString name1 = TString::Format("MuMuKMass_trk1pT%.1f", trk_pT[i]) ; name1.ReplaceAll(".","p") ;
    TString title1 = TString::Format("MuMuKMass with p_{T}(trk_{1}) > %.1f;m(#mu^{+}#mu^{-}K^{+}) [GeV]", trk_pT[i]) ; 
    MuMuKMass_trk1Pt.push_back( new TH1F(name1, title1, 130, 5.15, 5.41) ) ; 
    TString name2 = TString::Format("MuMuKMass_trk2pT%.1f", trk_pT[i]) ; name2.ReplaceAll(".","p") ;
    TString title2 = TString::Format("MuMuKMass with p_{T}(trk_{2}) > %.1f;m(#mu^{+}#mu^{-}K^{+}) [GeV]", trk_pT[i]) ; 
    MuMuKMass_trk2Pt.push_back( new TH1F(name2, title2, 130, 5.15, 5.41) ) ; 
  }
  for (Int_t i = 0; i < abs(Z_CL.size()); ++i) {
    TString name = TString::Format("MuMuKMass_ZVtxCL%.3f", Z_CL[i]) ; name.ReplaceAll(".","p") ;
    TString title = TString::Format("MuMuKMass with Z vertex CL > %.3f;m(#mu^{+}#mu^{-}K^{+}) [GeV]", Z_CL[i]) ; 
    MuMuKMass_ZCL.push_back( new TH1F(name, title, 130, 5.15, 5.41) ) ; 
  }
  
  // Nairit plots
  MuMuKMass_Nairit = new TH1F("MuMuKMass_Nairit", "MuMuKMass after Nairit cuts;m(#mu^{+}#mu^{-}K^{+}) [GeV]", 250, 5.15, 5.4) ;
  MuMuKMass_sb = new TH1F("MuMuKMass_sb", "MuMuKMass in sidebands;m(#mu^{+}#mu^{-}K^{+}) [GeV]", 250, 5.15, 5.4) ;
  trk_pT_test_h = new TH1F("K_pT_test", "K^{+} p_{T};p_{T} [GeV]", 150, 0, 15) ;

  // K* plots
  myExclusiveMuMuKPiMass = new TH1F("myExclusiveMuMuKPiMass","MuMuKPiMass with HLT5 and nB0.gt.0 BEFORE selection",100, 4.8, 5.8);
  kaonP_vs_pionP = new TH2F("kaonP_vs_pionP","p(K^{+}) vs p(#pi^{-});p(#pi^{-}) [GeV];p(K^{+}) [GeV]",100,0.,10.,100,0.,10.) ;
  kaonPt_vs_pionPt = new TH2F("kaonPt_vs_pionPt","p_{T}(K^{+}) vs p_{T}(#pi^{-});p_{T}(#pi^{-}) [GeV];p_{T}(K^{+}) [GeV]",100,0.,10.,100,0.,10.) ;
  trk_dxy = new TH1F("trk_Dxy", "trk d_{xy};d_{xy} [cm]", 400, -2, 2) ;
  trk_dz = new TH1F("trk_Dz", "trk d_{z};d_{z} [cm]", 300, -15, 15) ;
  // test
  hmyPsiPKPiMass_B0pTtest = new TH1F("myPsiPKPiMass_B0pTtest","myPsiPKPiMass B0pTtest", 100, 4.8, 5.8);
  hPionFromKstar_p = new TH1F("PionFromKstar_p","#pi^{-} from K*;p [GeV]", 100,0.,10.);
  hKaonFromKstar_p = new TH1F("KaonFromKstar_p","K^{+} from K*;p [GeV]", 100,0.,10.);
  hmyPsiPKPiMass_B0cTautest = new TH1F("myPsiPKPiMass_B0cTautest","myPsiPKPiMass B0cTautest", 100, 4.8, 5.8);
  hmyPsiPKPiMass_deltaRtest = new TH1F("myPsiPKPiMass_deltaRtest","myPsiPKPiMass deltaRtest", 100, 4.8, 5.8);
  //
  hmyPsiPKPiMassBaseSelAlt = new TH1F("myPsiPKPiMassBaseSelAlt","myPsiPKPiMassSel BASELINE ALTERNATIVE", 100, 4.8, 5.8);
  hPsiPrimefromBMassSel = new TH1F("PsiPrimefromBMassSel","PsiPrimefromBMassSel", 100, 3.45, 3.9); 
  hmyPsiPMassSel = new TH1F("myPsiPMassSel","myPsiPMassSel;m(#mu^{+}#mu^{-}) [GeV]", 100, 3.45, 3.9); 
  hmyPsiPKPiMassZoom = new TH1F("myPsiPKPiMassZoom","myPsiPKPiMassSel ZOOM;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high);
  hmyPsiPKPiMassAlt = new TH1F("myPsiPKPiMassAlt","myPsiPKPiMassSel ALTERNATIVE;m(#psi'K^{+}#pi^{-}) [GeV]", 200, 4.7, 6.7); 
  hmyPsiPKPiMassAltZoom = new TH1F("myPsiPKPiMassAltZoom","myPsiPKPiMassSel ALTERNATIVE ZOOM;m(#psi'K^{+}#pi^{-}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high);
  //
  hmyPiPiMass = new TH1F("myPiPiMass","PiPi mass;m(#pi^{+}#pi^{-}) [GeV]", 220, 0., 2.2);
  hmyPiKMass = new TH1F("myPiKMass","PiK mass;m(#pi^{+}K^{-}) [GeV]", 400, 0., 4.);
  hmyKKMass = new TH1F("myKKMass","KK mass;m(K^{+}K^{-}) [GeV]", 400, 0., 4.);
  hmyPsiPPiPiMass = new TH1F("myPsiPPiPiMass","PsiPPiPi mass;m(#psi'#pi^{+}#pi^{-}) [GeV]", 110, 4.7, 5.8);
  hmyPsiPPiKMass = new TH1F("myPsiPPiKMass","PsiPPiK mass;m(#psi'#pi^{+}K^{-}) [GeV]", 110, 4.7, 5.8);
  hmyPsiPPiKMass_zoom = new TH1F("myPsiPPiKMass_zoom","PsiPPiK mass;m(#psi'#pi^{+}K^{-}) [GeV]", 60, 5.15, 5.45);
  hmyPsiPKKMass = new TH1F("myPsiPKKMass","PsiPKK mass;m(#psi'K^{+}K^{-}) [GeV]", 200, 4.7, 6.7);
  hmyPsiPPhiMass = new TH1F("myPsiPPhiMass","PsiPPhi mass;m(#psi'#phi) [GeV]", 200, 4.7, 6.7);
  hmyPsiPPiKMass_fromBs = new TH1F("myPsiPPiKMass_fromBs","PsiPPiK mass from B_{s} region;m(#psi'#pi^{+}K^{-}) [GeV]", 110, 4.7, 5.8);
  hmyPsiPPiKMass_fromBs_zoom = new TH1F("myPsiPPiKMass_fromBs_zoom","PsiPPiK mass from B_{s} region;m(#psi'#pi^{+}K^{-}) [GeV]", 60, 5.15, 5.45);
  hmyPiPiMass_vs_KKMass = new TH2F("myPiPiMass_vs_KKMass","PiPi vs KK mass;m(#pi^{+}#pi^{-});m(K^{+}K^{-}) [GeV]", 220,0.,2., 400,0.,4.);
  hmyPsiPPiPiMass_vs_PiPiMass = new TH2F("myPsiPPiPiMass_vs_PiPiMass","PsiPPiPi vs PiPi mass;m(#psi'#pi^{+}#pi^{-});m(#pi^{+}#pi^{-}) [GeV]", 110,4.7,5.8, 220,0.,2.2);
  hmyPsiPKKMass_vs_KKMass = new TH2F("myPsiPKKMass_vs_KKMass","PsiPKK vs KK mass;m(#psi'K^{+}K^{-});m(K*{+}K^{-}) [GeV]", 200,4.7,6.7, 500,0.,5.);
  //
  psi2SPiSqMass_vs_KPiSqMass_AC = new TH2F("reco_psi2SPi_vs_KPi_dalitz_AC","#psi'#pi^{-} vs K^{+}#pi^{-} reco;m^{2}(K^{+}#pi^{-});m^{2}(#psi'#pi^{-}) [GeV]", 15*2,0.,3., 24*2,13.,25.) ;
  // tight cuts
  // single
  piPt_B0.push_back(0.5) ; piPt_B0.push_back(1.0) ; piPt_B0.push_back(1.5) ; piPt_B0.push_back(2.0) ; piPt_B0.push_back(2.5) ; piPt_B0.push_back(3.0) ; piPt_B0.push_back(3.5) ; piPt_B0.push_back(4.0) ; piPt_B0.push_back(4.5) ; piPt_B0.push_back(5.0) ;
  //KPt_B0.push_back(1.) ; 
  //KPt_B0.push_back( piPt_B0[i] ) ;
  vtxCL_B0.push_back(0.075) ; vtxCL_B0.push_back(0.1) ; vtxCL_B0.push_back(0.125) ; vtxCL_B0.push_back(0.15) ;
  cTau_B0.push_back(4.) ; cTau_B0.push_back(5.) ; cTau_B0.push_back(6.) ; cTau_B0.push_back(7.) ;
  // combined
  trksPt_B0VtxCL_B0CTau.push_back( {1.0, 0.1, 4.} ) ; trksPt_B0VtxCL_B0CTau.push_back( {2.0, 0.15, 5.} ) ; trksPt_B0VtxCL_B0CTau.push_back( {2.0, 0.2, 7.} ) ; trksPt_B0VtxCL_B0CTau.push_back( {3.0, 0.15, 6.} ) ; trksPt_B0VtxCL_B0CTau.push_back( {4.0, 0.2, 6.} ) ; trksPt_B0VtxCL_B0CTau.push_back( {5.0, 0.2, 7.} ) ;
  //
  for (Int_t i = 0; i < abs(piPt_B0.size()); ++i) {
    TString name = TString::Format("myMuMuKPiMass_piPt%.1f", piPt_B0[i] ) ; name.ReplaceAll(".","p") ;
    TString title = TString::Format("myMuMuKPiMass with p_{T}(#pi^{-}) > %.1f GeV;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", piPt_B0[i] ) ; 
    myMuMuKPiMass_piPt.push_back( new TH1F(name, title, MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ) ; 
  }
  //for (Int_t i = 0; i < abs(KPt_B0.size()); ++i) {
  for (Int_t i = 0; i < abs(piPt_B0.size()); ++i) {
    KPt_B0.push_back( piPt_B0[i] ) ;
    TString name = TString::Format("myMuMuKPiMass_kPT_%.1f", KPt_B0[i] ) ; name.ReplaceAll(".","p") ;
    TString title = TString::Format("myMuMuKPiMass with p_{T}(K^{+}) > %.1f GeV;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", KPt_B0[i] ) ; 
    myMuMuKPiMass_KPt.push_back( new TH1F(name, title, MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ) ; 
  }
  //
  for (Int_t i = 0; i < abs(vtxCL_B0.size()); ++i) {
    TString name = TString::Format("myMuMuKPiMass_vtxCL_%.2f", vtxCL_B0[i] ) ; name.ReplaceAll(".","p") ;
    TString title = TString::Format("myMuMuKPiMass with B^{0} vertex CL > %.2f;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", vtxCL_B0[i] ) ; 
    myMuMuKPiMass_vtxCL.push_back( new TH1F(name, title, MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ) ; 
  }
  //
  for (Int_t i = 0; i < abs(cTau_B0.size()); ++i) {
    TString name = TString::Format("myMuMuKPiMass_cTau_%.1f", cTau_B0[i] ) ; name.ReplaceAll(".","p") ;
    TString title = TString::Format("myMuMuKPiMass with c*#tau(B^{0})/#sigma > %.1f;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", cTau_B0[i] ) ; 
    myMuMuKPiMass_cTau.push_back( new TH1F(name, title, MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ) ; 
  }
  //
  for (Int_t i = 0; i < abs(trksPt_B0VtxCL_B0CTau.size()); ++i) {
    TString name = TString::Format("myMuMuKPiMass_trks%.1f_vtxCL%.2f_cTau%.1f", trksPt_B0VtxCL_B0CTau[i].trksPt, trksPt_B0VtxCL_B0CTau[i].B0VtxCL, trksPt_B0VtxCL_B0CTau[i].B0CTau ) ; name.ReplaceAll(".","p") ;
    TString title = TString::Format("myMuMuKPiMass with p_{T}(K^{+},#pi^{-}) > %.1f GeV, B^{0} vertex CL > %.2f, c*#tau(B^{0})/#sigma > %.1f;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", trksPt_B0VtxCL_B0CTau[i].trksPt, trksPt_B0VtxCL_B0CTau[i].B0VtxCL, trksPt_B0VtxCL_B0CTau[i].B0CTau ) ; 
    myMuMuKPiMass_trksAndB0vtxCLAndB0Ctau.push_back( new TH1F(name, title, MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ) ; 
  }
  //
  // PV matching
  priVtx_deltaZ = new TH1F("priVtx_deltaz", "#Deltaz(standard PV, B^{0}_cos(#alpha) PV);#Deltaz", 300, 0, 30) ;
  B0_pointingAngle_PVvsB0LessPV[0] = new TH2F("B0_cosAlpha_PVvsB0less_pT", "B^{0} cos(#alpha) wrt PV vs wrt B0LessPV (best #Sigma(p_{T}^{2}));B^{0} cos(#alpha) wrt PV;B^{0} cos(#alpha) wrt B^{0}LessPV", 1000, 0.995, 1., 1000, 0.995, 1.) ;
  B0_pointingAngle_PVvsB0LessPV[1] = new TH2F("B0_cosAlpha_PVvsB0less_cos", "B^{0} cos(#alpha) wrt PV vs wrt B0LessPV (best cos(#alpha));B^{0} cos(#alpha) wrt PV;B^{0} cos(#alpha) wrt B^{0}LessPV", 1000, 0.995, 1., 1000, 0.995, 1.) ; 
  z_PVvsB0LessPV[0] = new TH2F("z_PVvsB0less_pT", "z(PV) vs wrt z(B0LessPV) (best #Sigma(p_{T}^{2}));z(PV) [cm];z(B0LessPV) [cm]", 100,-25,25, 100,-25,25) ;
  z_PVvsB0LessPV[1] = new TH2F("z_PVvsB0less_cos", "z(PV) vs wrt z(B0LessPV) (best cos(#alpha));z(PV) [cm];z(B0LessPV) [cm]", 100,-25,25, 100,-25,25) ;
  B0_PVCTau_h = new TH1F("B0_PVCTau", "B^{0} flight lenght signif. w.r.t PV;c*#tau_{PV} / #sigma(c*#tau_{PV})", 200,0,40) ;
  B0_CosAlphaCTau_h = new TH1F("B0_CosAlphaCTau", "B^{0} flight lenght signif. w.r.t B^{0}_cos(#alpha) PV;c*#tau_{B^{0}_cos(#alpha)} / #sigma(c*#tau_{B^{0}_cos(#alpha)})", 200,0,40) ;
  B0_PVXCTau_h = new TH1F("B0_PVXCTau", "B^{0} flight lenght signif. w.r.t PVX;c*#tau_{PVX} / #sigma(c*#tau_{PVX})", 200,0,40) ;
  B0_PVX3DCTau_h = new TH1F("B0_PVX3DCTau", "B^{0} flight lenght signif. w.r.t PVX3D;c*#tau_{PVX3D} / #sigma(c*#tau_{PVX3D})", 200,0,40) ;
  priVtx_matched = new TH1I("priVtx_matched", "standard PV matched with B^{0}_cos(#alpha) PV;n(PV matched)", 4, -1.5, 2.5) ;
  //
  hmyPsiPKPiMass = new TH1F("myPsiPKPiMass","RECALCULATED B0Mass with HLT5 and nB0.gt.0 BEFORE selection", 100, 4.8, 5.8); //check
  hmyPsiPKPiMass_KPiExchanged = new TH1F("myPsiPKPiMass_KPiExchanged","RECALCULATED B0Mass (with HLT5 and nB0.gt.0) BEF sel with inverted K-Pi mass hypotheses", 100, 4.8, 5.8); //check
  hmyPsiPKPiMass_bothComb = new TH1F("myPsiPKPiMass_bothComb","myPsiPKPiMass_bothComb", 100, 4.8, 5.8); //check
  //
  hB0VtxCL = new TH1F("B0VtxCL","B0VtxCL", 200, 0., 1.);
  hB0_Pt = new TH1F("B0_Pt","B0_Pt", 100, 0., 30.);
  //
  hB0CTauPVSignif = new TH1F("B0CTauPVSignif","B0CTauPVSignif", 120, 0., 12.);
  hB0CTauPVXSignif = new TH1F("B0CTauPVXSignif","B0CTauPVXSignif", 120, 0., 12.);
  hPsiPrimefromBMass = new TH1F("PsiPrimefromBMass","PsiPrimefromBMass", 100, 3.45, 3.9); 
  hPsiPrime_fromB_Pt = new TH1F("PsiPrime_fromB_Pt","PsiPrime_fromB_Pt", 100, 0., 30.);
  hmyPsiP_Pt = new TH1F("myPsiP_Pt","myPsiP_Pt", 100, 0., 30.);
  hPionfromBpt = new TH1F("PionfromB_pT","#pi^{-} from B^{0};p_{T} [GeV]", 200,0.,10.);
  hKaonfromBpt = new TH1F("KaonfromB_pT","K^{+} from B^{0};p_{T} [GeV]", 200,0.,10.);
  hPionFromKstar_pT = new TH1F("PionFromKstar_pT","#pi^{-} from K*;p_{T} [GeV]", 100, 0., 10.);
  hKaonFromKstar_pT = new TH1F("KaonFromKstar_pT","K^{+} from K*;p_{T} [GeV]", 100, 0., 10.);
  hPionFromKstar_ex_pT = new TH1F("PionFromKstar_ex_pT","#pi^{-} from K* exchanged;p_{T} [GeV]", 100, 0., 10.);
  hKaonFromKstar_ex_pT = new TH1F("KaonFromKstar_ex_pT","K^{+} from K* exchanged;p_{T} [GeV]", 100, 0., 10.);
  hPionFromKstar_sb_pT = new TH1F("PionFromKstar_sb_pT","#pi^{-} from K* sidebands;p_{T} [GeV]", 100, 0., 10.);
  hKaonFromKstar_sb_pT = new TH1F("KaonFromKstar_sb_pT","K^{+} from K* sidebands;p_{T} [GeV]", 100, 0., 10.);
  //
  myExclMuMuKPiMassSel = new TH1F("myExclMuMuKPiMassSel","MuMuKPiMass  AFTER selection",100, 4.8, 5.8);
  hmyPsiPKPiMass_bothCombSel = new TH1F("myPsiPKPiMass_bothCombSel","myPsiPKPiMass_bothCombSel",100, 4.8, 5.8); //check
  psi2SPiSqMass_vs_KPiSqMass_ABaseC = new TH2F("reco_psi2SPi_vs_KPi_dalitz_ABaseC","#psi'#pi^{-} vs K^{+}#pi^{-} reco;m^{2}(K^{+}#pi^{-});m^{2}(#psi'#pi^{-})", 45,0.2,3.8, 56,12.3,36.3) ;
  hmyKPiMass = new TH1F("myKPiMass","K*;m(K^{+}#pi^{-}) [GeV]", 120,0.6,1.8) ;
  hmyKPiMass_ex = new TH1F("myKPiMass_ex","K* exchanged;m(K^{+}#pi^{-}) [GeV]", 120,0.6,1.8) ;
  hmyKPiMass_sb = new TH1F("myKPiMass_sb","myKPiMass sidebands;m(K^{+}#pi^{-}) [GeV]", 200, 0., 2.) ;
  // high purity
  hmyKPiMass_KpTcut = new TH1F("myKPiMass_KpTcut","K* with p_{T}(K^{+}) < 1 GeV;m(K^{+}#pi^{-}) [GeV]", 120,0.6,1.8) ;
  hmyKPiMass_ex_KpTcut = new TH1F("myKPiMass_ex_KpTcut","K* exchanged with p_{T}(K^{+}) < 1 GeV;m(K^{+}#pi^{-}) [GeV]", 120,0.6,1.8) ;
  hmyKPiMass_KpTcut_tightTr = new TH1F("myKPiMass_KpTcut_tightTr","K* with p_{T}(K^{+}) < 1 GeV (tight tracks);m(K^{+}#pi^{-}) [GeV]", 120,0.6,1.8) ;
  hmyKPiMass_ex_KpTcut_tightTr = new TH1F("myKPiMass_ex_KpTcut_tightTr","K* exchanged with p_{T}(K^{+}) < 1 GeV (tight tracks);m(K^{+}#pi^{-}) [GeV]", 120,0.6,1.8) ;
  hmyKPiMass_KpTcut_HP = new TH1F("myKPiMass_KpTcut_HPTr","K* with p_{T}(K^{+}) < 1 GeV (HP tracks);m(K^{+}#pi^{-}) [GeV]", 120,0.6,1.8) ;
  hmyKPiMass_ex_KpTcut_HP = new TH1F("myKPiMass_ex_KpTcut_HPTr","K* exchanged with p_{T}(K^{+}) < 1 GeV (HP tracks);m(K^{+}#pi^{-}) [GeV]", 120,0.6,1.8) ;
  //
  hmyPsiPKPiMass_fromKstar = new TH1F("myPsiPKPiMass_fromKstar","myPsiPKPiMass_fromKstar;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hmyPsiPKPiMass_fromKstar_ex_noRange = new TH1F("myPsiPKPiMass_fromKstar_ex_noRange","myPsiPKPiMass_fromKstar exchanged;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]",160,4.8,5.6);
  hmyPsiPKPiMass_fromKstar_diff = new TH1F("myPsiPKPiMass_fromKstar_ex_diff","myPsiPKPiMass_fromKstar diff;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) - m(#mu^{+}#mu^{-}K^{-}#pi^{+}) [GeV]",80,-0.4,0.4);
  hmyPsiPKPiMass_fromKstar_ex = new TH1F("myPsiPKPiMass_fromKstar_ex","myPsiPKPiMass_fromKstar exchanged;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hmyPsiPKPiMass_fromKstar_sb = new TH1F("myPsiPKPiMass_fromKstar_sb","myPsiPKPiMass_fromKstar sidebands;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hpsiPK_vs_KPi = new TH2F("psiPK_vs_KPi","m(#psi'K^{+}) vs m(K^{+}#pi^{-}) from K*;m(K^{+}#pi^{-}) [GeV];m(#psi'K^{+}) [GeV]", 120,0.6,1.8, 200,3.6,6.) ;
  hpsiPK_vs_KPi_ex = new TH2F("psiPK_vs_KPi_ex","m(#psi'K^{+}) vs m(K^{+}#pi^{-}) from K* exchanged;m(K^{+}#pi^{-}) [GeV];m(#psi'K^{+}) [GeV]", 120,0.6,1.8, 200,3.6,6.);
  hpsiPK_vs_KPi_sb = new TH2F("psiPK_vs_KPi_sb","m(#psi'K^{+}) vs m(K^{+}#pi^{-}) from K* sidebands;m(K^{+}#pi^{-}) [GeV];m(#psi'K^{+}) [GeV]", 130,0.7,2.0, 200,3.6,6.);
  dEdx_vs_pT_Kstar = new TH2F("dEdx_vs_pT_Kstar", "dE/dx vs p_{T} for K* decay;p_{T} [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;
  dEdx_vs_p_Kstar = new TH2F("dEdx_vs_p_Kstar", "dE/dx vs p for K* decay;p [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;
  dEdx_vs_pT_piFromKstar = new TH2F("dEdx_vs_pT_piFromKstar", "dE/dx vs p_{T} for #pi^{-} from K*;p_{T}(#pi^{-}) [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;
  dEdx_vs_pT_KFromKstar = new TH2F("dEdx_vs_pT_KFromKstar", "dE/dx vs p_{T} for K^{+} from K*;p_{T}(K^{+}) [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;
  dEdx_vs_p_piFromKstar = new TH2F("dEdx_vs_p_piFromKstar", "dE/dx vs p for #pi^{-} from K*;p(#pi^{-}) [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;
  dEdx_vs_p_KFromKstar = new TH2F("dEdx_vs_p_KFromKstar", "dE/dx vs p for K^{+} from K*;p(K^{+}) [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;
  dEdx_vs_p_piFromKstar_KpCut = new TH2F("dEdx_vs_p_piFromKstar_KpCut", "dE/dx vs p for #pi^{-} from K* for p(K^{+}) < 1 GeV;p(#pi^{-}) [GeV];dE/dx [MeV/cm]", 120,0,3, 160,0,8) ;
  dEdx_vs_p_KFromKstar_KpCut = new TH2F("dEdx_vs_p_KFromKstar_KpCut", "dE/dx vs p for K^{+} from K* for p(K^{+}) < 1 GeV;p(K^{+}) [GeV];dE/dx [MeV/cm]", 120,0,3, 160,0,8) ;
  dEdxK_vs_dEdxPi = new TH2F("dEdxK_vs_dEdxPi", "dE/dx(K^{+}) vs dE/dx(#pi^{-}) from K*;dE/dx(#pi^{-}) [MeV/cm];dE/dx(K^{+}) [MeV/cm]", 200,0,5, 320, 0, 8) ;
  piFromKstar_stripHits_h = new TH1I("piFromKstar_stripHits", "#pi^{-} strip hits;strip hits", 25, 0, 25) ;
  KFromKstar_stripHits_h = new TH1I("KFromKstar_stripHits", "K^{+} strip hits;strip hits", 25, 0, 25) ;
  piFromKstar_stripHits_vs_y_h = new TH2F("y_vs_piFromKstar_stripHits", "y vs strip hits #pi^{-};rapidity;strip hits", 66,-3.3,3.3, 30,0,30) ;
  KFromKstar_stripHits_vs_y_h = new TH2F("y_vs_KFromKstar_stripHits", "y vs strip hits K^{+};rapidity;strip hits", 66,-3.3,3.3, 30,0,30) ;
  dEdx_vs_pT_KhqFromKstar = new TH2F("dEdx_vs_pT_KhqFromKstar", "dE/dx vs p_{T} for K^{+} from K*;p_{T}(K^{+}) [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;
  dEdx_vs_p_KhqFromKstar = new TH2F("dEdx_vs_p_KhqFromKstar", "dE/dx vs p for K^{+} from K*;p(K^{+}) [GeV];dE/dx [MeV/cm]", 200,0,5, 200, 0, 10) ;
  dEdxPi_pCut = new TH1F("dEdxPi_pCut", "dE/dx(#pi^{-} from K*) for p(K^{+}) < 0.8;dE/dx [MeV/cm]", 80,0,8) ;
  dEdxK_pCut = new TH1F("dEdxK_pCut", "dE/dx(K^{+} from K*) for p(K^{+}) < 0.8;dE/dx [MeV/cm]", 80,0,8) ;
  dEdxPi_pCut_stripCut = new TH1F("dEdxPi_pCut_stripCut", "dE/dx(#pi^{-} from K*) for p(K^{+}) < 0.8 && stripHits >= 17;dE/dx [MeV/cm]", 80,0,8) ;
  dEdxK_pCut_stripCut = new TH1F("dEdxK_pCut_stripCut", "dE/dx(K^{+} from K*) for p(K^{+}) < 0.8 && stripHits >= 17;dE/dx [MeV/cm]", 80,0,8) ;
  Armenteros_Kstar_orig = new TH2F("Armenteros_Kstar_orig", "p_{T} vs #alpha for K*;#alpha = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}));p_{T}(#pi^{-}) = p_{T}(K^{+}) [GeV]", 200,-1,1, 80,0,0.4) ;
  Armenteros_Kstar = new TH2F("Armenteros_Kstar", "p_{T} vs #alpha for K*;#alpha = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}));p_{T}(#pi^{-}) = p_{T}(K^{+}) [GeV]", 200,-1,1, 80,0,0.4) ;
  Armenteros_Kstar_orig_ex = new TH2F("Armenteros_Kstar_orig_ex", "p_{T} vs #alpha for K* exchanged;#alpha = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}));p_{T}(#pi^{-}) = p_{T}(K^{+}) [GeV]", 200,-1,1, 80,0,0.4) ;
  Armenteros_Kstar_orig_sb = new TH2F("Armenteros_Kstar_orig_sb", "p_{T} vs #alpha for K* sidebands;#alpha = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}));p_{T}(#pi^{-}) = p_{T}(K^{+}) [GeV]", 200,-1,1, 160,0,0.8) ;
  alpha_vs_alphaEx = new TH2F("alpha_vs_alphaEx", "#alpha vs #alpha_{ex} for K*;#alpha = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}));#alpha_{ex} = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}))", 200,-1,1, 200,-1,1) ;
  pT_diff = new TH1F("pT_diff","p_{T}(#pi^{-}) - p_{T}(K^{+}) from K*;p_{T}(#pi^{-}) - p_{T}(K^{+}) [GeV]",200,-0.001,0.001) ;
  pT_diff_ex = new TH1F("pT_diff_ex","p_{T}(#pi^{-}) - p_{T}(K^{+}) from K* exchanged;p_{T}(#pi^{-}) - p_{T}(K^{+}) [GeV]",200,-0.001,0.001) ;
  pT_diff_sb = new TH1F("pT_diff_sb","p_{T}(#pi^{-}) - p_{T}(K^{+}) from K* sidebands;p_{T}(#pi^{-}) - p_{T}(K^{+}) [GeV]",200,-0.001,0.001) ;
  //
  
  // my scan
  // pi and/or Z pT cuts
  mode.push_back("all") ; mode.push_back("pr") ; mode.push_back("pr_barrel") ; mode.push_back("quasipr") ; mode.push_back("nonpr") ; mode.push_back("nonpr_barrel") ;
  myMuMuPiMass_piPt.resize(mode.size()) ;  myMuMuPiMass_ZPt.resize(mode.size()) ;  myMuMuPiMass_piAndZPt.resize(mode.size()) ; 

  //pi_pT.push_back(0.6) ; pi_pT.push_back(0.8) ; // almost ineffective
  pi_pT.push_back(1.) ; pi_pT.push_back(2.) ; pi_pT.push_back(3.) ; pi_pT.push_back(4.) ; pi_pT.push_back(5.) ; pi_pT.push_back(7.) ;
  Z_pT.push_back(10.) ; Z_pT.push_back(12.) ; Z_pT.push_back(14.) ; Z_pT.push_back(17.) ; Z_pT.push_back(20.) ; Z_pT.push_back(28.) ; Z_pT.push_back(40.) ;
  piAndZ_pT.push_back( make_pair(1.,15.) ) ; piAndZ_pT.push_back( make_pair(2.,20.) ) ; piAndZ_pT.push_back( make_pair(4.,25.) ) ; piAndZ_pT.push_back( make_pair(7.,30.) ) ; piAndZ_pT.push_back( make_pair(10.,30.) ) ; piAndZ_pT.push_back( make_pair(12.,30.) ) ; piAndZ_pT.push_back( make_pair(15.,30.) ) ; piAndZ_pT.push_back( make_pair(10.,40.) ) ; piAndZ_pT.push_back( make_pair(12.,40.) ) ; piAndZ_pT.push_back( make_pair(12.,50.) ) ;
  
  deltaR_piPsi2S.push_back(1.0) ; deltaR_piPsi2S.push_back(0.8) ; deltaR_piPsi2S.push_back(0.6) ; deltaR_piPsi2S.push_back(0.4) ; deltaR_piPsi2S.push_back(0.2) ; deltaR_piPsi2S.push_back(0.1) ;

  cosAlpha_cut = 0.9 ; cosAlpha_cut2 = 0.99 ;
  scale_factor = 3. ;
  Z_cosPA_min = -1. ; Z_cosPA_max = 1. ;
  //pair<Float_t, Float_t> Z_cosPA_min = make_pair(-1.,0.85) ;
  //pair<Float_t, Float_t> Z_cosPA_max = make_pair(-0.85,1.) ; 
  Int_t Z_cosPA_nbins = 200 ;
  //Int_t Z_cosPA_nbins = 115 ;
  Float_t Z_cosPA_binSize = (Z_cosPA_max-Z_cosPA_min)/Z_cosPA_nbins ;

  for (Int_t i = 0; i < abs(mode.size()); ++i) {
    TString name_Z_cosPA = TString::Format("Z_cosAlpha_%s", mode[i].Data()) ; name_Z_cosPA.ReplaceAll(".","p") ;
    TString title_Z_cosPA = TString::Format("Z cos(#alpha) %s;cos(#alpha)", mode[i].Data()) ; 
    Z_cosPA.push_back( new TH1F(name_Z_cosPA, title_Z_cosPA, Z_cosPA_nbins, Z_cosPA_min, Z_cosPA_max) ) ; // 1000, -1., 1. // 100, 0.999, 1. ;
    //Z_cosPA.push_back( make_pair(new TH1F(name_Z_cosPA+"_low", title_Z_cosPA+" low", Z_cosPA_nbins, Z_cosPA_min.first, Z_cosPA_max.first), 
		       //new TH1F(name_Z_cosPA+"_high", title_Z_cosPA+" high", Z_cosPA_nbins, Z_cosPA_min.second, Z_cosPA_max.second) ) ) ; // 1000, -1., 1. // 100, 0.999, 1. ;
    Z_cosPA[i]->GetXaxis()->SetBinLabel(1, "-1.0") ;
    Z_cosPA[i]->GetXaxis()->SetBinLabel( (-cosAlpha_cut/scale_factor -Z_cosPA_min)/Z_cosPA_binSize, "-0.9") ;
    Z_cosPA[i]->GetXaxis()->SetBinLabel( (cosAlpha_cut/scale_factor -Z_cosPA_min)/Z_cosPA_binSize, "0.9") ;
    Z_cosPA[i]->GetXaxis()->SetBinLabel(Z_cosPA_nbins, "1.0") ;
    Z_cosPA[i]->LabelsOption("h") ;  Z_cosPA[i]->SetLabelSize(0.05,"X") ;
    Z_cosPA[i]->SetLineWidth(2) ;

    TString name_DeltaR = TString::Format("psi2S_DeltaR_pi_%s", mode[i].Data()) ; name_DeltaR.ReplaceAll(".","p") ;
    TString title_DeltaR = TString::Format("#DeltaR(#psi',#pi^{+}) %s;#DeltaR(#psi',#pi^{+})", mode[i].Data()) ; 
    psi2S_DeltaR_pi_h.push_back( new TH1F(name_DeltaR, title_DeltaR, 110, 0, 1.1) ) ; // 160, 0, 1.6
    psi2S_DeltaR_pi_h[i]->SetLineWidth(2) ;

    TString name_DeltaR_vs_cosPA = TString::Format("psi2S_DeltaR_pi_vs_Z_cosAlpha_%s", mode[i].Data()) ; name_DeltaR_vs_cosPA.ReplaceAll(".","p") ;
    TString title_DeltaR_vs_cosPA = TString::Format("#DeltaR(#psi',#pi^{+}) vs Z |cos(#alpha)| %s;Z |cos(#alpha)|;#DeltaR(#psi',#pi^{+})", mode[i].Data()) ; 
    Z_DeltaR_vs_cosPA.push_back( new TH2F(name_DeltaR_vs_cosPA, title_DeltaR_vs_cosPA, 100, cosAlpha_cut2, Z_cosPA_max, 110, 0, 1.1) ) ;
    

    myMuMuPiMass_piPt[i].resize(pi_pT.size()) ;
    myMuMuPiMass_ZPt[i].resize(Z_pT.size()) ;
    myMuMuPiMass_piAndZPt[i].resize(piAndZ_pT.size()) ; 
      
    for (Int_t k = 0; k < abs(deltaR_piPsi2S.size()); ++k) {

      for (Int_t j = 0; j < abs(pi_pT.size()); ++j) {
	//TString nameTitle = TString::Format("myMuMuPiMass_pi%.1f_%s", pi_pT[j], mode[i].Data()) ; 
	TString name = TString::Format("myMuMuPiMass_%s_pi%.1f_DeltaR%.1f", mode[i].Data(), pi_pT[j], deltaR_piPsi2S[k]) ; name.ReplaceAll(".","p") ;
	TString title = TString::Format("myMuMuPiMass %s with p_{T}(#pi^{+}) > %.1f GeV, #DeltaR(#psi',#pi^{+}) < %.1f;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", mode[i].Data(), pi_pT[j], deltaR_piPsi2S[k] ) ; 
	//myMuMuPiMass_piPt[i][j] = new TH1F(nameTitle, nameTitle, MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high);
	myMuMuPiMass_piPt[i][j].push_back( new TH1F(name, title, MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high) ) ; 
	myMuMuPiMass_piPt[i][j][k]->SetLineWidth(2) ;
      }
      
      for (Int_t j = 0; j < abs(Z_pT.size()); ++j) {
	//TString nameTitle = TString::Format("myMuMuPiMass_Z%.1f_%s", Z_pT[j], mode[i].Data()) ; 
	TString name = TString::Format("myMuMuPiMass_%s_Z%.1f_DeltaR%.1f", mode[i].Data(), Z_pT[j], deltaR_piPsi2S[k]) ; name.ReplaceAll(".","p") ;
	TString title = TString::Format("myMuMuPiMass %s with p_{T}(Z^{+}) > %.1f GeV, #DeltaR(#psi',#pi^{+}) < %.1f;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", mode[i].Data(), Z_pT[j], deltaR_piPsi2S[k] ) ; 
	//myMuMuPiMass_ZPt[i][j] = new TH1F(nameTitle, nameTitle, MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high); 
	myMuMuPiMass_ZPt[i][j].push_back( new TH1F(name, title, MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high) ) ; 
	myMuMuPiMass_ZPt[i][j][k]->SetLineWidth(2) ;
      }
      
      for (Int_t j = 0; j < abs(piAndZ_pT.size()); ++j) {
	//TString nameTitle = TString::Format("myMuMuPiMass_Z%.1f_%s", Z_pT[j], mode[i].Data()) ; 
	TString name = TString::Format("myMuMuPiMass_%s_pi%.1f_Z%.1f_DeltaR%.1f", mode[i].Data(), piAndZ_pT[j].first, piAndZ_pT[j].second, deltaR_piPsi2S[k]) ; name.ReplaceAll(".","p") ;
	TString title = TString::Format("myMuMuPiMass %s with p_{T}(#pi^{+},Z^{+}) > (%.1f,%.1f) GeV, #DeltaR(#psi',#pi^{+}) < %.1f;m(#mu^{+}#mu^{-}#pi^{+}) [GeV]", mode[i].Data(), piAndZ_pT[j].first, piAndZ_pT[j].second, deltaR_piPsi2S[k] ) ; 
	//myMuMuPiMass_piAndZPt[i][j] = new TH1F(nameTitle, nameTitle, MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high); 
	myMuMuPiMass_piAndZPt[i][j].push_back( new TH1F(name, title, MuMuPiMass_bins,MuMuPiMass_low,MuMuPiMass_high) ) ; 
	myMuMuPiMass_piAndZPt[i][j][k]->SetLineWidth(2) ;
      }
    } // end DeltaR cuts loop
  } // end mode loop
  psi2S_DeltaR_pi_vs_mumupiMass_nonpr = new TH2F("psi2S_DeltaR_pi_vs_myMuMuPiMass_nonpr", "#DeltaR(#psi',#pi^{+}) vs m(#mu^{+}#mu^{-}#pi^{+});m(#mu^{+}#mu^{-}#pi^{+}) [GeV];#DeltaR(#psi',#pi^{+})", MuMuPiMass_bins, MuMuPiMass_low, MuMuPiMass_high, 110, 0, 1.1) ;

  // psi' histos
  // need to initialize varables from .h here since this method is called before Init method in .h
  rapid_min = -2.4; rapid_max = 2.4; rapid_bins = 24;
  rapid_binSize = (rapid_max - rapid_min)/rapid_bins ;
  psiPrime_rapidBin.resize(rapid_bins) ; psiPrime_2011_rapidBin.resize(rapid_bins) ;
  psiPrimeMass_min = MuMuMass_low ;  psiPrimeMass_max = MuMuMass_high ;  psiPrimeMass_bins = MuMuMass_bins ;
  
  for (Int_t i = 0; i < rapid_bins; i++) {
    Float_t bin_low = rapid_min + (rapid_max - rapid_min)/rapid_bins*i ;
    Float_t bin_high = bin_low + (rapid_max - rapid_min)/rapid_bins ;
    
    TString psiPrimeMass_binName = TString::Format("PsiPrimeMass_bin%d", i+1) ; 
    TString psiPrimeMass_binTitle = TString::Format("#psi' for y in [%.1f, %.1f];m(#mu^{+}#mu^{-}) [GeV]", bin_low, bin_high) ; 
    psiPrime_rapidBin[i] = new TH1F(psiPrimeMass_binName.Data(), psiPrimeMass_binTitle.Data(), psiPrimeMass_bins, psiPrimeMass_min, psiPrimeMass_max);

    TString psiPrimeMass_2011_binName = TString::Format("PsiPrimeMass_2011_bin%d", i+1) ; 
    psiPrime_2011_rapidBin[i] = new TH1F(psiPrimeMass_2011_binName.Data(), psiPrimeMass_binTitle.Data(), psiPrimeMass_bins, psiPrimeMass_min, psiPrimeMass_max);
  }
  // selection variables
  // Z- variables
  nZ_h = new TH1I("Z_multiplicity", "Z multiplicity;Z multiplicity;Events", 150, -0.5, 149.5) ;
  nB0AC_noMassWin_h = new TH1I("B0_multiplicity_AC_anyMass", "B0 multiplicity AC;B^{0} multiplicity;Events", 10, -0.5, 9.5) ;
  nB0AC_h = new TH1I("B0_multiplicity_AC", "B0 multiplicity AC;B^{0} multiplicity;Events", 10, -0.5, 9.5) ;
  nB0AC_signalWin_h = new TH1I("B0_multiplicity_AC_signalWin", "B0 multiplicity AC in signal Win;B^{0} multiplicity;Events", 10, -0.5, 9.5) ;
  nB0AC_purityWin_h = new TH1I("B0_multiplicity_AC_purityWin", "B0 multiplicity AC in purity Win;B^{0} multiplicity;Events", 10, -0.5, 9.5) ;
  hB0Mass_1B0 = new TH1F("B0Mass_1B0", "mass of single B^{0} candidates;m(B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_1B0matched[0][0] = new TH1F("B0Mass_1B0matched", "mass of single matched B^{0} candidates;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_1B0matched[0][1] = new TH1F("B0Mass_1B0matched_dauCh", "mass of single matched (with daughters charge) B^{0} candidates;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_1B0matched[1][0] = new TH1F("B0Mass_1B0signalWin_matched", "mass of single matched B^{0} candidates in signal;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_1B0matched[1][1] = new TH1F("B0Mass_1B0signalWin_matched_dauCh", "mass of single matched (with daughters charge) B^{0} candidates in signal;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0CTau_1B0 = new TH1F("B0CTau_1B0", "c*#tau of single B^{0} candidates;c*#tau(B^{0}_{reco}) [GeV]", 100,0,20) ;
  hB0CTau_1B0matched = new TH1F("B0CTau_1B0matched", "c*#tau of single matched B^{0} candidates;c*#tau(matched B^{0}_{reco}) [GeV]", 100,0,20) ;
  hB0Mass_noTwins_noSignalWinNotTwins = new TH1F("hB0Mass_noTwins_noSignalWinNotTwins", "mass of B^{0} candidates w/o (twins & both B^{0} in signal);m(B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  psi2SPiSqMass_vs_KPiSqMass_1B0[0] = new TH2F("psi2SPi_vs_KPi_dalitz_1B0","#psi'#pi^{-} vs K^{+}#pi^{-} of single B^{0} candidates;m^{2}(K^{+}#pi^{-});m^{2}(#psi'#pi^{-})", 2*15,0.,3., 2*24,13.,25.) ;
  psi2S_helicityAngle = new TH1F("psi2S_helicityAngle","#psi' helicity angle reco;#theta_{#psi'}", 160, -0.03, 3.17) ;
  cos_psi2S_helicityAngle = new TH1F("cos_psi2S_helicityAngle","cosine of #psi' helicity angle reco;cos #theta_{#psi'}", 102, -1.02, 1.02) ;
  cos_psi2S_helicityAngle_rot = new TH1F("cos_psi2S_helicityAngle_rot","cosine of #psi' helicity angle with rotation;cos #theta_{#psi'} rot", 102, -1.02, 1.02) ;
  cos_psi2S_helicityAngle_rotVSnotRot = new TH2F("cos_psi2S_helicityAngle_rotVSnotRot","cosine of #psi' helicity angle rot VS notRot;cos #theta_{#psi'} not rot;cos #theta_{#psi'} rot", 102,-1.02,1.02, 102,-1.02,1.02) ;
  planesAngle = new TH1F("planesAngle","angle between K#pi and #mu#mu planes reco;#phi(#psi',K*)", 160, -3.2, 3.2) ;
  hB0Mass_2B0twin = new TH1F("B0Mass_2B0twin", "mass of 2 B^{0} twin;m(twin B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  // MC analysis
  nMCB0_BT_h = new TH1I("MC_B0_multiplicity_BT", "B^{0} MC multiplicity before trigger;B^{0} multiplicity;Events", 10, -0.5, 9.5) ;
  psi2SPiSqMass_vs_KPiSqMass_BT_gen = new TH2F("gen_psi2SPi_vs_KPi_dalitz_BT","#psi'#pi^{-} vs K^{+}#pi^{-} gen before trigger;m^{2}(K^{+}#pi^{-});m^{2}(#psi'#pi^{-})", 2*15,0.,3., 2*24,13.,25.) ;
  KPiSqMass_vs_psi2SKSqMass_BT_gen = new TH2F("gen_KPi_vs_psi2SK_dalitz_BT","K^{+}#pi^{-} vs #psi'K^{+} gen before trigger;m^{2}(#psi'K^{+});m^{2}(K^{+}#pi^{-})", 2*24,13.,25., 2*15,0.,3.) ;
  nMCB0_h = new TH1I("MC_B0_multiplicity", "B^{0} MC multiplicity;B^{0} multiplicity;Events", 10, -0.5, 9.5) ;
  B0_gen_p_h = new TH1F("B0_gen_p", "gen B^{0} p;p [GeV]", 200, 0, 40) ;
  B0_gen_pT_h = new TH1F("B0_gen_pT", "gen B^{0} p_{T};p_{T} [GeV]", 200, 0, 40) ;
  psi2S_gen_p_h = new TH1F("psi2S_gen_p", "gen #psi' p;p [GeV]", 150, 0, 30) ;
  psi2S_gen_pT_h = new TH1F("psi2S_gen_pT", "gen #psi' p_{T};p_{T} [GeV]", 150, 0, 30) ;
  K_gen_p_h = new TH1F("K_gen_p", "gen K^{-} p;p [GeV]", 150, 0, 15) ;
  K_gen_pT_h = new TH1F("K_gen_pT", "gen K^{-} p_{T};p_{T} [GeV]", 150, 0, 15) ;
  pi_gen_p_h = new TH1F("pi_gen_p", "gen #pi^{+} p;p [GeV]", 150, 0, 15) ;
  pi_gen_pT_h = new TH1F("pi_gen_pT", "gen #pi^{+} p_{T};p_{T} [GeV]", 150, 0, 15) ;
  //
  psi2SPiSqMass_vs_KPiSqMass_gen = new TH2F("gen_psi2SPi_vs_KPi_dalitz","#psi'#pi^{-} vs K^{+}#pi^{-} gen;m^{2}(K^{+}#pi^{-});m^{2}(#psi'#pi^{-})", 2*15,0.,3., 2*24,13.,25.) ;
  psi2S_vs_mu_pT = new TH2F("psi2S_vs_mu_pT","p_{T}(#psi(2S)) vs minimum p_{T}(#mu);p_{T}(#psi(2S));min p_{T}(#mu)",150, 0, 30, 150, 0, 30);
  psi2S_helicityAngle_BT_gen = new TH1F("psi2S_helicityAngle_BT_gen","#psi' helicity angle gen;#theta_{#psi'}", 160, -0.03, 3.17) ;
  cos_psi2S_helicityAngle_BT_gen = new TH1F("cos_psi2S_helicityAngle_BT_gen","cosine of #psi' helicity angle gen;cos #theta_{#psi'}", 102, -1.02, 1.02) ;
  cos_psi2S_helicityAngle_BT_gen_rot = new TH1F("cos_psi2S_helicityAngle_BT_gen_rot","cosine of #psi' helicity angle gen with rotation;cos #theta_{#psi'} rot", 102, -1.02, 1.02) ;
  cos_psi2S_helicityAngle_BT_gen_rotVSnotRot = new TH2F("cos_psi2S_helicityAngle_BT_gen_rotVSnotRot","cosine of #psi' helicity angle gen rot VS notRot;cos #theta_{#psi'} not rot;cos #theta_{#psi'} rot", 102,-1.02,1.02, 102,-1.02,1.02) ;
  cos_psi2S_helicityAngle_vs_psi2SPiMass = new TH2F("cos_psi2S_helicityAngle_vs_psi2SPi_BT_gen","cosine of #psi' helicity angle vs m(#psi'#pi^{-}) gen;m(#psi'#pi^{-});cos #theta_{#psi'}", 2*24,3.8,4.8, 102,-1.02,1.02) ;
  cos_psi2S_helicityAngle_vs_KPiMass = new TH2F("cos_psi2S_helicityAngle_vs_KPi_BT_gen","cosine of #psi' helicity angle vs m(K^{+}#pi^{-}) gen;m(K^{+}#pi^{-});cos #theta_{#psi'}", 2*25,0.58,1.62, 102,-1.02,1.02) ;
  planesAngle_BT_gen = new TH1F("planesAngle_BT_gen","angle between K#pi and #mu#mu planes gen;#phi(#psi',K*)", 160, -3.2, 3.2) ;
  //
  B0CosAlpha_2Dvs3D[0] = new TH2F("PV_2Dvs3D", "B^{0}_cos(#alpha)_3D vs B^{0}_cos(#alpha)_2D wrt standard PV;B^{0} cos(#alpha) 2D;B^{0} cos(#alpha) 3D", 2000, 0.8, 1., 2000, 0.8, 1.) ;
  B0CosAlpha_2Dvs3D[1] = new TH2F("B0CosAlpha2D_vs_B0CosAlpha3D", "best B^{0}_cos(#alpha)_3D vs best B^{0}_cos(#alpha)_2D wrt standard PV;B^{0} cos(#alpha) 2D;B^{0} cos(#alpha) 3D", 2000, 0.8, 1., 2000, 0.8, 1.) ;
  B0CosAlpha_2Dvs3D[2] = new TH2F("B0LessPV_2Dvs3D", "B^{0}_cos(#alpha)_3D vs B^{0}_cos(#alpha)_2D;B^{0} cos(#alpha) 2D wrt B^{0}_lessPV;B^{0} cos(#alpha) 3D", 2000, 0.8, 1., 2000, 0.8, 1.) ;
  B0CosAlpha_2Dvs3D[3] = new TH2F("B0LessB0CosAlpha2D_vs_B0LessB0CosAlpha3D", "best B^{0}_cos(#alpha)_3D vs best B^{0}_cos(#alpha)_2D wrt B^{0}_lessPV;B^{0} cos(#alpha) 2D;B^{0} cos(#alpha) 3D", 2000, 0.8, 1., 2000, 0.8, 1.) ;
  //
  priVtx_z[0] = new TH1F("priVtx_z", "interaction PV z;z [cm]", 300, -30, 30) ;
  priVtx_z[1] = new TH1F("priVtx_B0_z", "standard B^{0} PV z;z [cm]", 300, -30, 30) ;
  priVtx_z[2] = new TH1F("priVtx_B0CosAlpha_z", "best B^{0}_cos(#alpha) PV z;z [cm]", 300, -30, 30) ;
  priVtx_z[3] = new TH1F("priVtx_B0CosAlpha3D_z", "best B^{0}_cos(#alpha)_3D PV z;z [cm]", 300, -30, 30) ;
  priVtx_z[4] = new TH1F("priVtx_B0Less_z", "B^{0}_less PV z;z [cm]", 300, -30, 30) ;
  priVtx_z[5] = new TH1F("priVtx_B0LessB0CosAlpha_z", "best B^{0}_cos(#alpha) B^{0}_lessPV z;z [cm]", 300, -30, 30) ;
  priVtx_z[6] = new TH1F("priVtx_B0LessB0CosAlpha3D_z", "best B^{0}_cos(#alpha)_3D B^{0}_lessPV z;z [cm]", 300, -30, 30) ;
  
  TString xyz[3] = {"x","y","z"}, XYZ[3] = {"X","Y","Z"} ;
  TString zoom[2] = {"","_zoom"} ;
  TString dim[2] = {"2D","3D"} ;
  Float_t binHigh_xyz[2][3] = {{0.1, 0.1, 30.}, {0.02, 0.02, 0.1}} ;
  Int_t nBins[2][3] = {{200, 200, 100}, {100, 100, 200}} ;

  for (Int_t j=0; j<2; j++) {
    for (Int_t i=0; i<3; i++) {
      priVtx_delta[j][i] = new TH1F(TString::Format("priVtx_delta%s%s", XYZ[i].Data(), zoom[j].Data()), TString::Format("#Delta%s(interaction PV, gen B^{0} PV);#Delta%s [cm]", xyz[i].Data(), xyz[i].Data()), nBins[j][i]*2, -binHigh_xyz[j][i]/10, binHigh_xyz[j][i]/10) ;
      priVtx_delta[j][i]->SetLineWidth(2) ; //priVtx_delta[j][i]->SetLineColor(1) ; 
      priVtx_B0_delta[j][i] = new TH1F(TString::Format("priVtx_B0_delta%s%s", XYZ[i].Data(), zoom[j].Data()), TString::Format("#Delta%s(interaction PV, standard PV);#Delta%s [cm]", xyz[i].Data(), xyz[i].Data()), nBins[j][i], -binHigh_xyz[j][i], binHigh_xyz[j][i]) ;
      priVtx_B0_delta[j][i]->SetLineColor(1) ; priVtx_B0_delta[j][i]->SetLineWidth(2) ;
      for (Int_t k=0; k<2; k++) {
	priVtx_B0CosAlpha_delta[k][j][i] = new TH1F(TString::Format("priVtx_B0CosAlpha%s_delta%s%s", dim[k].Data(), XYZ[i].Data(), zoom[j].Data()), TString::Format("#Delta%s(interaction PV, best B^{0}_cos(#alpha_{%s}) PV);#Delta%s [cm]", xyz[i].Data(), dim[k].Data(), xyz[i].Data()), nBins[j][i], -binHigh_xyz[j][i], binHigh_xyz[j][i]) ;
	priVtx_B0CosAlpha_delta[k][j][i]->SetLineColor(2) ; // 2 = red
	priVtx_B0CosAlpha_delta[k][j][i]->SetLineWidth(2) ;
      }
      priVtxB0Less_delta[j][i] = new TH1F(TString::Format("priVtxB0Less_delta%s%s", XYZ[i].Data(), zoom[j].Data()), TString::Format("#Delta%s(interaction PV, B^{0}_lessPV);#Delta%s [cm]", xyz[i].Data(), xyz[i].Data()), nBins[j][i], -binHigh_xyz[j][i], binHigh_xyz[j][i]) ;
      priVtxB0Less_delta[j][i]->SetLineColor(3) ; // 3 = green
      priVtxB0Less_delta[j][i]->SetLineWidth(2) ;
      for (Int_t k=0; k<2; k++) {
	priVtxB0Less_B0CosAlpha_delta[k][j][i] = new TH1F(TString::Format("priVtxB0Less_B0CosAlpha%s_delta%s%s", dim[k].Data(), XYZ[i].Data(), zoom[j].Data()), TString::Format("#Delta%s(interaction PV, best B^{0}_cos(#alpha_{%s}) B^{0}_lessPV);#Delta%s [cm]", xyz[i].Data(), dim[k].Data(), xyz[i].Data()), nBins[j][i], -binHigh_xyz[j][i], binHigh_xyz[j][i]) ;
	priVtxB0Less_B0CosAlpha_delta[k][j][i]->SetLineColor(4) ; // 4 = blue
	priVtxB0Less_B0CosAlpha_delta[k][j][i]->SetLineWidth(2) ;
      }
    }
    dz_pri_priB0Less[j] = new TH1F(TString::Format("priVtx_priVtxB0Less_deltaZ%s", zoom[j].Data()), "#Deltaz(standard PV, best B^{0}_lessPV);#Deltaz [cm]", nBins[j][2], -30, 30) ;
    for (Int_t k=0; k<2; k++) {
      dz_pri_priB0CosAlpha[k][j] = new TH1F(TString::Format("priVtx_priVtxB0CosAlpha%s_deltaZ%s", dim[k].Data(), zoom[j].Data()), TString::Format("#Deltaz(standard PV, best B^{0}_cos(#alpha_{%s}) PV);#Deltaz [cm]", dim[k].Data()), nBins[j][2], -30, 30) ;
      dz_pri_priB0LessB0CosAlpha[k][j] = new TH1F(TString::Format("priVtx_priVtxB0LessB0CosAlpha%s_deltaZ%s", dim[k].Data(), zoom[j].Data()), TString::Format("#Deltaz(standard PV, best B^{0}_cos(#alpha_{%s}) B^{0}_lessPV);#Deltaz [cm]", dim[k].Data()), nBins[j][2], -30, 30) ;
      dz_priB0CosAlpha_priB0LessB0CosAlpha[k][j] = new TH1F(TString::Format("priVtxB0CosAlpha%s_priVtxB0LessB0CosAlpha%s_deltaZ%s", dim[k].Data(), dim[k].Data(), zoom[j].Data()), TString::Format("#Deltaz(best B^{0}_cos(#alpha)_%s PV, best B^{0}_cos(#alpha_{%s}) B^{0}_lessPV);#Deltaz [cm]", dim[k].Data(), dim[k].Data()), nBins[j][2], -30, 30) ;
      dz_priB0Less_priB0LessB0CosAlpha[k][j] = new TH1F(TString::Format("priVtxB0Less_priVtxB0LessB0CosAlpha%s_deltaZ%s", dim[k].Data(), zoom[j].Data()), TString::Format("#Deltaz(best B^{0}_lessPV, best B^{0}_cos(#alpha_{%s}) B^{0}_lessPV);#Deltaz [cm]", dim[k].Data()), nBins[j][2], -30, 30) ;
    }
  }
  priVtx_deltaX_test = new TH1F("priVtx_deltaX_test", "#Deltax(interaction PV, gen B^{0} PV) after #Deltaz matching;#Deltax [cm]", 100, -0.02, 0.02) ;
  priVtx_deltaY_test = new TH1F("priVtx_deltaY_test", "#Deltay(interaction PV, gen B^{0} PV) after #Deltaz matching;#Deltay [cm]", 100, -0.02, 0.02) ;
  //
  B0_gen_mass_h = new TH1F("B0_gen_mass","gen B^{0} mass;m(#mu^{+}#mu^{-}K^{+}#pi^{-}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high);
  piCh_vs_KCh_gen = new TH2F("gen_piCh_vs_KCh","gen K charge vs gen #pi charge;#pi charge;K charge", 5,-2.5,2.5, 5,-2.5,2.5) ;
  hMCDeltaRK_fromTwins = new TH1F("MCDeltaRK_fromTwins", "#DeltaR for kaon from B^{0} twins;#DeltaR(K_{gen},K_{reco})",100, 0, 1);
  hMCPtRecoVsPtTrueK_fromTwins = new TH2F("MCPtRecoVsPtTrueK_fromTwins", "p_{T}(K_{reco}) vs p_{T}(K_{gen}) from B^{0} twins;p_{T}(K_{gen}) [GeV];p_{T}(K_{reco}) [GeV]", 80,0.,40., 80,0.,40.);
  hMCDeltaPtK_fromTwins = new TH1F("MCDeltaPtK_fromTwins", "#Deltap_{T}/p_{T} for K from B^{0} twins;|p_{T}(K^{+}_{gen}) - p_{T}(K^{+}_{reco})| / p_{T}(K^{+}_{gen})", 250,0,0.5) ;
  hMCDeltaPtK_fromTwinsPiMatched = new TH1F("MCDeltaPtK_fromTwinsPiMatched", "#Deltap_{T}/p_{T} for K from B^{0} twins with matched #pi;|p_{T}(K^{+}_{gen}) - p_{T}(K^{+}_{reco})| / p_{T}(K^{+}_{gen})", 250,0,0.5) ;
  hMCpTTrueK_fromTwins = new TH1F("MCpTTrueK_fromTwins", "matched K p_{T} from B^{0} twins;p_{T}(K_{reco}) [GeV]", 150, 0, 15);
  //
  nB0ACInMC_h = new TH1I("B0_multiplicity_ACInMC", "B0 multiplicity AC in MC;B^{0} multiplicity;Events", 10, -0.5, 9.5) ;
  psi2SPiSqMass_vs_KPiSqMass_1B0[1] = new TH2F("psi2SPi_vs_KPi_dalitz_1B0_signalWin","#psi'#pi^{-} vs K^{+}#pi^{-} of single B^{0} candidates in signal;m^{2}(K^{+}#pi^{-});m^{2}(#psi'#pi^{-})", 15,0.,3., 2*24,13.,25.) ;
  hMCDeltaR_1B0 = new TH1F("MCDeltaR_1B0", "#DeltaR for 1 B^{0} candidates;#DeltaR(B^{0}_{gen},B^{0}_{reco})", 200, 0, 0.2) ;
  hMCDeltaPt_1B0 = new TH1F("MCDeltaPt_1B0", "#Deltap_{T}/p_{T} for 1 B^{0} candidates;|p_{T}(B^{0}_{gen}) - p_{T}(B^{0}_{reco})| / p_{T}(B^{0}_{gen})", 100, 0, 1) ;
  hMCTruthMatching_1B0[0][0] = new TH1I("MCTruthMatching_1B0", "MC matching for 1 B^{0} candidate;n(B^{0} matched);Events", 4, -1.5, 2.5) ;
  hMCTruthMatching_1B0[0][1] = new TH1I("MCTruthMatching_1B0_dauCh", "MC matching (with daughters charge) for 1 B^{0} candidate;n(B^{0} matched);Events", 4, -1.5, 2.5) ;
  hMCTruthMatching_1B0[0][2] = new TH1I("MCTruthMatching_1B0fromDau", "MC truth daughters matching for 1 B^{0} candidate;n(B^{0} matched);Events", 4, -1.5, 2.5) ;
  for (Int_t i = 0; i < 2; ++i) {
    TString nameR = TString::Format("MCDeltaR_%dof2B0", i+1) ; 
    TString titleR = TString::Format("#DeltaR for %d of 2 B^{0} candidates;#DeltaR(B^{0}_{gen},B^{0}_{reco})", i+1) ;
    hMCDeltaR_2B0.push_back( new TH1F(nameR, titleR, 200, 0, 0.2) ) ;
    TString nameP = TString::Format("MCDeltaPt_%dof2B0", i+1) ; 
    TString titleP = TString::Format("#Deltap_{T}/p_{T} for %d of 2 B^{0} candidates;|p_{T}(B^{0}_{gen}) - p_{T}(B^{0}_{reco})| / p_{T}(B^{0}_{gen})", i+1) ;
    hMCDeltaPt_2B0.push_back( new TH1F(nameP, titleP, 100, 0, 1) ) ;
  }
  hMCDeltaPt_1vs2_2B0 = new TH2F("MCDeltaPt_1vs2_2B0", "#Deltap_{T}/p_{T} 1^{st} vs 2^{nd} of 2 B^{0} candidates;|p_{T}(B^{0}_{gen}) - p_{T}(1^{st}B^{0}_{reco})| / p_{T}(B^{0}_{gen});|p_{T}(B^{0}_{gen}) - p_{T}(2^{nd}B^{0}_{reco})| / p_{T}(B^{0}_{gen})", 200,0,0.2, 200,0,0.2) ;
  hMCDeltaR_1vs2_2B0 = new TH2F("MCDeltaR_1vs2_2B0", "#DeltaR 1^{st} vs 2^{nd} of 2 B^{0} candidates;#DeltaR(B^{0}_{gen},1^{st}B^{0}_{reco});#DeltaR(B^{0}_{gen},2^{nd}B^{0}_{reco})", 200,0,0.2, 200,0,0.2) ;
  // twins
  hMCDeltaR_1vs2_2B0twin = new TH2F("MCDeltaR_1vs2_2B0twin", "#DeltaR 1^{st} vs 2^{nd} of 2 B^{0} twin candidates;#DeltaR(B^{0}_{gen},1^{st}B^{0}_{reco});#DeltaR(B^{0}_{gen},2^{nd}B^{0}_{reco})", 200,0,0.2, 200,0,0.2) ;
  hMCDeltaPt_1vs2_2B0twin = new TH2F("MCDeltaPt_1vs2_2B0twin", "#Deltap_{T}/p_{T} 1^{st} vs 2^{nd} of 2 B^{0} twin candidates;|p_{T}(B^{0}_{gen}) - p_{T}(1^{st}B^{0}_{reco})| / p_{T}(B^{0}_{gen});|p_{T}(B^{0}_{gen}) - p_{T}(2^{nd}B^{0}_{reco})| / p_{T}(B^{0}_{gen})", 100,0,0.02, 100,0,0.02) ;
  hB0VtxCL_1vs2_2B0twin = new TH2F("B0VtxCL_1vs2_2B0twin", "vertex CL 1^{st} vs 2^{nd} of 2 B^{0} twin candidates;vertex CL(1^{st}B^{0}_{reco});vertex CL(2^{nd}B^{0}_{reco})", 200,0,1, 200,0,1) ;
  hB0Mass_1vs2_2B0twin = new TH2F("B0Mass_1vs2_2B0twin", "mass 1^{st} vs 2^{nd} of 2 B^{0} twin candidates;m(1^{st}B^{0}_{reco}) [GeV];m(2^{nd}B^{0}_{reco}) [GeV]", 60,5.25,5.31, 60,5.25,5.31) ;
  hB0CTau_1vs2_2B0twin = new TH2F("B0CTau_1vs2_2B0twin", "c*#tau 1^{st} vs 2^{nd} of 2 B^{0} twin candidates;c*#tau(1^{st}B^{0}_{reco}) [GeV];c*#tau(2^{nd}B^{0}_{reco}) [GeV]", 100,0,20, 100,0,20) ;
  hAlpha_1vs2_2B0twin = new TH2F("Alpha_1vs2_2B0twin", "#alpha 1^{st} vs 2^{nd} of 2 B^{0} twin candidates;#alpha_{1} = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}));#alpha_{2} = (p_{L}(#pi^{-}) - p_{L}(K^{+}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}))", 200,-1,1, 200,-1,1) ;
  // matched twins
  hMCDeltaR_matchedTwin[1] = new TH1F("MCDeltaR_matchedTwin", "#DeltaR of matched B^{0} between 2 B^{0} twin;#DeltaR(B^{0}_{gen},matched B^{0}_{reco});", 250,0,0.05) ;
  hMCDeltaPt_matchedTwin[1] = new TH1F("MCDeltaPt_matchedTwin", "#Deltap_{T}/p_{T} of matched B^{0} between 2 B^{0} twin;|p_{T}(B^{0}_{gen}) - p_{T}(matched B^{0}_{reco})| / p_{T}(B^{0}_{gen})", 100,0,0.02) ;
  hB0VtxCL_matchedTwin[1] = new TH1F("B0VtxCL_matchedTwin", "vertex CL of matched B^{0} between 2 B^{0} twin;vertex CL(matched B^{0}_{reco})", 200,0,1) ;
  hB0Mass_matchedTwin[1] = new TH1F("B0Mass_matchedTwin", "mass of matched B^{0} between 2 B^{0} twin;m(matched B^{0}_{reco}) [GeV]", 60,5.25,5.31) ;
  hB0CTau_matchedTwin[1] = new TH1F("B0CTau_matchedTwin", "c*#tau of matched B^{0} between 2 B^{0} twin;c*#tau(matched B^{0}_{reco}) [GeV]", 100,0,20) ;
  hAlpha_matchedTwin[1] = new TH1F("Alpha_matchedTwin", "#alpha of matched B^{0} between 2 B^{0} twin;#alpha_{matched} = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}))", 200,-1,1) ;
  // not matched
  hMCDeltaR_matchedTwin[0] = new TH1F("MCDeltaR_notMatchedTwin", "#DeltaR of not matched B^{0} between 2 B^{0} twin;#DeltaR(B^{0}_{gen},not matched B^{0}_{reco});", 250,0,0.05) ;
  hMCDeltaPt_matchedTwin[0] = new TH1F("MCDeltaPt_notMatchedTwin", "#Deltap_{T}/p_{T} of not matched B^{0} between 2 B^{0} twin;|p_{T}(B^{0}_{gen}) - p_{T}(not matched B^{0}_{reco})| / p_{T}(B^{0}_{gen})", 100,0,0.02) ;
  hB0VtxCL_matchedTwin[0] = new TH1F("B0VtxCL_notMatchedTwin", "vertex CL of not matched B^{0} between 2 B^{0} twin;vertex CL(not matched B^{0}_{reco})", 200,0,1) ;
  hB0Mass_matchedTwin[0] = new TH1F("B0Mass_notMatchedTwin", "mass of not matched B^{0} between 2 B^{0} twin;m(not matched B^{0}_{reco}) [GeV]", 60,5.25,5.31) ;
  hB0CTau_matchedTwin[0] = new TH1F("B0CTau_notMatchedTwin", "c*#tau of not matched B^{0} between 2 B^{0} twin;c*#tau(not matched B^{0}_{reco}) [GeV]", 100,0,20) ;
  hAlpha_matchedTwin[0] = new TH1F("Alpha_notMatchedTwin", "#alpha of not matched B^{0} between 2 B^{0} twin;#alpha_{not matched} = (p_{L}(K^{+}) - p_{L}(#pi^{-}))/(p_{L}(K^{+}) + p_{L}(#pi^{-}))", 200,-1,1) ;
  nTwins_h[0] = new TH1I("B0twins_multiplicity_AC", "B0 twins multiplicity AC;B^{0} twins;Events", 4, -1.5, 2.5) ;
  nTwins_h[1] = new TH1I("B0twins_multiplicity_AC_signalWin", "B0 twins multiplicity AC in signal;B^{0} twins;Events", 4, -1.5, 2.5) ;
  hB0Mass_2B0notTwin = new TH1F("B0Mass_2B0notTwin", "mass of 2 B^{0} not twin;m(not twin B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hMCTruthMatching_2B0notTwin[0] = new TH1I("MCTruthMatching_2B0notTwin", "MC matching for 2 not twin B^{0} candidates;n(not twin B^{0} matched);Events", 5, -1.5, 3.5) ;
  hMCTruthMatching_2B0notTwin[1] = new TH1I("MCTruthMatching_2B0notTwin_dauCh", "MC matching (with daughters charge) for 2 not twin B^{0} candidates;n(not twin B^{0} matched);Events", 5, -1.5, 3.5) ;
  hMCTruthMatching_2B0notTwin[2] = new TH1I("MCTruthMatching_2B0notTwin_fromDau", "MC truth daughters matching for 2 not twin B^{0} candidates;n(not twin B^{0} matched);Events", 5, -1.5, 3.5) ;
  hMCDeltaRPi_2B0 = new TH1F("MCDeltaRPi_2B0", "#DeltaR for pions from 2 not twin B^{0} caniddates;#DeltaR(#pi_{gen},#pi_{reco})",200, 0, 2);
  hMCDeltaRK_2B0 = new TH1F("MCDeltaRK_2B0", "#DeltaR for kaons from 2 not twin B^{0} caniddates;#DeltaR(K_{gen},K_{reco})",200, 0, 2);
  //
  hB0Mass_2B0notTwin_matched[0] = new TH1F("B0Mass_2B0notTwin_matched", "mass of matched not twin B^{0};m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_2B0notTwin_matched[1] = new TH1F("B0Mass_2B0notTwin_matched_dauCh", "mass of matched (with daughters charge) not twin B^{0} candidates;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_2B0notTwin_matched[2] = new TH1F("B0Mass_2B0notTwin_matched_fromDau", "mass of truth matched not twin B^{0} candidates;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hDeltaB0Mass_2B0notTwin[0] = new TH2F("DeltaB0Mass_2B0notTwin", "#Deltam(matched) vs #Deltam(not matched) of 2 not twin B^{0} candidates;#Deltam(matched B^{0}_{reco}) [GeV];#Deltam(not matched B^{0}_{reco}) [GeV]", 100,-0.05,0.05, 100,-0.05,0.05) ;
  hMCTruthMatching_2B0twin[0] = new TH1I("MCTruthMatching_2B0twin", "MC matching for 2 twin B^{0} candidates;n(twin B^{0} matched);Events", 5, -1.5, 3.5) ;
  hMCTruthMatching_2B0twin[1] = new TH1I("MCTruthMatching_2B0twin_dauCh", "MC matching (with daughters charge) for 2 twin B^{0} candidates;n(twin B^{0} matched);Events", 5, -1.5, 3.5) ;
  hMCTruthMatching_2B0twin[2] = new TH1I("MCTruthMatching_2B0twin_fromDau", "MC truth daughters matching for 2 twin B^{0} candidates;n(twin B^{0} matched);Events", 5, -1.5, 3.5) ;
  hB0Mass_2B0twin_matched[0] = new TH1F("B0Mass_2B0twin_matched", "mass of matched twin B^{0};m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_2B0twin_matched[1] = new TH1F("B0Mass_2B0twin_matched_dauCh", "mass of matched (with daughters charge) twin B^{0} candidates;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_2B0twin_matched[2] = new TH1F("B0Mass_2B0twin_matched_fromDau", "mass of truth matched twin B^{0} candidates;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hDeltaB0Mass_2B0twin[0] = new TH2F("DeltaB0Mass_2B0twin", "#Deltam(matched) vs #Deltam(not matched) of 2 B^{0} twin candidates;#Deltam(matched B^{0}_{reco}) [GeV];#Deltam(not matched B^{0}_{reco}) [GeV]", 100,-0.2,0.2, 100,-0.2,0.2) ;
  hAlpha_2B0twin[0] = new TH2F("alpha_2B0twin", "#alpha(matched) vs #alpha(not matched) of 2 B^{0} twin candidates;#alpha(matched B^{0}_{reco});#alpha(not matched B^{0}_{reco})", 200,-1,1, 200,-1,1) ;
  hB0Mass_2B0twin_hourglass[0][0] = new TH1F("B0Mass_outside_notMatchedTwin", "mass of not matched B^{0} twin outside;m(not matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_2B0twin_hourglass[0][1] = new TH1F("B0Mass_outside_matchedTwin", "mass of matched B^{0} twin outside;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_2B0twin_hourglass[1][0] = new TH1F("B0Mass_inside_notMatchedTwin", "mass of not matched B^{0} twin inside;m(not matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  hB0Mass_2B0twin_hourglass[1][1] = new TH1F("B0Mass_inside_matchedTwin", "mass of matched B^{0} twin inside;m(matched B^{0}_{reco}) [GeV]", MuMuKPiMass_bins, MuMuKPiMass_low, MuMuKPiMass_high) ;
  //
  hMCTruthMatching_1B0[1][0] = new TH1I("MCTruthMatching_1B0signalWin", "MC matching for 1 signal B^{0} candidate;n(B^{0} matched);Events", 4, -1.5, 2.5) ;
  hMCTruthMatching_1B0[1][1] = new TH1I("MCTruthMatching_1B0signalWin_dauCh", "MC matching (with daughters charge) for 1 signal B^{0} candidate;n(B^{0} matched);Events", 4, -1.5, 2.5) ;
  hMCTruthMatching_1B0[1][2] = new TH1I("MCTruthMatching_1B0signalWin_fromDau", "MC truth daughters matching for 1 signal B^{0} candidate;n(B^{0} matched);Events", 4, -1.5, 2.5) ;
  hMCTruthMatching_2B0signalWin_notTwin = new TH1I("MCTruthMatching_2B0signalWin_notTwin", "MC truth matching for 2 not twin signal B^{0} candidates;n(not twin B^{0} matched);Events", 5, -1.5, 3.5) ;
  hDeltaB0Mass_2B0twin[1] = new TH2F("DeltaB0Mass_2B0twin_signalWin", "#Deltam(matched) vs #Deltam(not matched) of 2 signal B^{0} twin candidates;#Deltam(matched B^{0}_{reco}) [GeV];#Deltam(not matched B^{0}_{reco}) [GeV]", 100,-0.05,0.05, 100,-0.05,0.05) ;
  hAlpha_2B0twin[1] = new TH2F("alpha_2B0twin_signalWin", "#alpha(matched) vs #alpha(not matched) of 2 signal B^{0} twin candidates;#alpha(matched B^{0}_{reco});#alpha(not matched B^{0}_{reco})", 200,-1,1, 200,-1,1) ;
  hDeltaB0Mass_2B0notTwin[1] = new TH2F("DeltaB0Mass_2B0signalWin_notTwin", "#Deltam(matched) vs #Deltam(not matched) of 2 not twin signal B^{0} candidates;#Deltam(matched B^{0}_{reco}) [GeV];#Deltam(not matched B^{0}_{reco}) [GeV]", 100,-0.05,0.05, 100,-0.05,0.05) ;
  // plots with candidate generated variables
  //hMCnGenX = new TH1I ("NGenX", "Nr of generated Cand",5,-0.5,4.5);
  /*
  GenXPlots * GenXPlots_NoSel = new GenXPlots("NoSel");
  AcceptanceMapPlots * AccMap_NoSel = new AcceptanceMapPlots("NoSel");
  AcceptanceMapPlots * AccMap_NoSel_cowboy = new AcceptanceMapPlots("NoSel_cowboy");
  AcceptanceMapPlots * AccMap_NoSel_seagull = new AcceptanceMapPlots("NoSel_seagull");
  AcceptanceMapPlots * AccMap_NoSel_reweight = new AcceptanceMapPlots("NoSel_reweight");
  AcceptanceMapPlots * AccMap_NoSel_reweightplus = new AcceptanceMapPlots("NoSel_reweightplus");
  AcceptanceMapPlots * AccMap_NoSel_reweightminus = new AcceptanceMapPlots("NoSel_reweightminus");
  GenXPlots * GenXPlots_AfterTrig = new GenXPlots("AfterTrig");
  GenXPlots * GenXPlots_WithJPsi = new GenXPlots("WithJPsi");
  AcceptanceMapPlots * AccMap_WithJPsi = new AcceptanceMapPlots("WithJPsi");
  AcceptanceMapPlots * AccMap_WithJPsi_cowboy = new AcceptanceMapPlots("WithJPsi_cowboy");
  AcceptanceMapPlots * AccMap_WithJPsi_seagull = new AcceptanceMapPlots("WithJPsi_seagull");
  AcceptanceMapPlots * AccMap_WithJPsi_reweight = new AcceptanceMapPlots("WithJPsi_reweight");
  AcceptanceMapPlots * AccMap_WithJPsi_reweightplus = new AcceptanceMapPlots("WithJPsi_reweightplus");
  AcceptanceMapPlots * AccMap_WithJPsi_reweightminus = new AcceptanceMapPlots("WithJPsi_reweightminus");
  GenXPlots * GenXPlots_WithJPsiInAcc = new GenXPlots("WithJPsiInAcc");
  AcceptanceMapPlots * AccMap_WithJPsiAndPions = new AcceptanceMapPlots("WithJPsiAndPions");
  GenXPlots * GenXPlots_WithX = new GenXPlots("WithX");
  GenXPlots * GenXPlots_WithXFinal = new GenXPlots("WithXFinal");
  GenXPlots * GenXPlots_NotMatched = new GenXPlots("NotMatched");
  AcceptanceMapPlots * AccMap_WithXFinal = new AcceptanceMapPlots("WithXFinal");
  AcceptanceMapPlots * AccMap_WithXFinalCTau = new AcceptanceMapPlots("WithXFinalCTau");
  AcceptanceMapPlots * AccMap_WithXFinal_cowboy = new AcceptanceMapPlots("WithXFinal_cowboy");
  AcceptanceMapPlots * AccMap_WithXFinal_seagull = new AcceptanceMapPlots("WithXFinal_seagull");
  AcceptanceMapPlots * AccMap_WithXFinal_reweight = new AcceptanceMapPlots("WithXFinal_reweight");
  AcceptanceMapPlots * AccMap_WithXFinal_reweightplus = new AcceptanceMapPlots("WithXFinal_reweightplus");
  AcceptanceMapPlots * AccMap_WithXFinal_reweightminus = new AcceptanceMapPlots("WithXFinal_reweightminus");
  
  // plots with reco variables for True X
  Sidebands * MCTrueX = new Sidebands ("MCTrueX");
  Sidebands * MCTrueXAfterJPsi = new Sidebands ("MCTrueXAfterJPsi");
  Sidebands * MCTrueXAfterQualityTrack = new Sidebands ("MCTrueXAfterQualityTrack");
  Sidebands * MCTrueXAfterMomentumTrack = new Sidebands ("MCTrueXAfterMomentumTrack");
  Sidebands * MCTrueXAfterVertex = new Sidebands ("MCTrueXAfterVertex");
  Sidebands * MCTrueXAfterDeltaR = new Sidebands ("MCTrueXAfterDeltaR");
  Sidebands * MCTrueXFinal = new Sidebands ("MCTrueXFinal");
  */
  SpectrumFromBPlus = new TH1F("SpectrumFromBPlus", "SpectrumFromBPlus",120,3.6,4.0);
  SpectrumFromBPlusCtaued = new TH1F("SpectrumFromBPlusCtaued", "SpectrumFromBPlusCtaued",120,3.6,4.0);
  
  // mumu matching
  hMCDeltaRMuMuMatch = new TH1F("MCDeltaRMuMuMatch", "#DeltaR for #psi' matching;#DeltaR(#psi'_{gen},#psi'_{reco})",100, 0, 0.1);
  hMCPtRecoVsPtTrueMuMuMatch = new TH2F("MCPtRecoVsPtTrueMuMuMatch", "p_{T}(#psi'_{reco}) vs p_{T}(#psi'_{gen});p_{T}(#psi'_{gen}) [GeV];p_{T}(#psi'_{reco}) [GeV]", 80,0.,40., 80,0.,40.);
  hMCpTTrueMuMuMatch = new TH1F("MCpTTrueMuMuMatch", "matched #psi' p_{T};p_{T} [GeV]", 200, 0, 20);
  // pion matching
  hMCDeltaRPiMatch = new TH1F("MCDeltaRPiMatch", "#DeltaR for pion matching;#DeltaR(#pi_{gen},#pi_{reco})",100, 0, 1);
  hMCPtRecoVsPtTruePiMatch = new TH2F("MCPtRecoVsPtTruePiMatch", "p_{T}(#pi_{reco}) vs p_{T}(#pi_{gen});p_{T}(#pi_{gen}) [GeV];p_{T}(#pi_{reco}) [GeV]", 80,0.,40., 80,0.,40.);
  hMCPtTruePiMatch = new TH1F("MCPtTruePiMatch", "matched #pi p_{T};p_{T} [GeV]", 150, 0, 15);
  hMCStripsHitsTruePiMatch = new TH1F("Dan_MCXTrueNrShits", "Nr Strip hits in True Pions",50, 0, 50);
  hMCPixelHitsTruePiMatch = new TH1F("Dan_MCXTrueNrPixel", "Nr Pixel hits in True Pions",10, 0, 10);
  hMCTotalHitsTruePiMatch = new TH1F("Dan_MCXTrueTotalHits", "Total Hits hits in True Pions",60, 0, 60);
  hMCChiTruePiMatch = new TH1F("Dan_MCXTrueChi", "Chi/NDF hits in True Pions",100, 0, 10);
  hMCStripsHitsTrackNotMached = new TH1F("Dan_NoMatched_MCXTrueNrShits", "Nr Strip hits in False Pions",50, 0, 50);
  hMCPixelHitsTrackNotMached = new TH1F("Dan_NoMatched_MCXTrueNrPixel", "Nr Pixel hits in False Pions",10, 0, 10);
  hMCTotalHitsTrackNotMached = new TH1F("Dan_NoMatched_MCXTrueTotalHits", "Total Hits hits in False Pions",60, 0, 60);
  hMCChiTrackNotMached = new TH1F("Dan_NoMatched_MCXTrueChi", "Chi/NDF hits in False Pions",100, 0, 10);
  //
  kaonDr_vs_pionDr_matched = new TH2F("kaonDr_vs_pionDr_matched","#DeltaR(K) vs #DeltaR(#pi);#DeltaR(#pi);#DeltaR(K)", 50,0,0.1, 50,0,0.1) ;
  kaonP_vs_pionP_matched = new TH2F("kaonP_vs_pionP_matched","p(K_{matched}) vs p(#pi_{matched});p(#pi_{matched}) [GeV];p(K_{matched}) [GeV]",200,0,20, 200,0,20) ;
  kaonPt_vs_pionPt_matched = new TH2F("kaonPt_vs_pionPt_matched","p_{T}(K_{matched}) vs p_{T}(#pi_{matched});p_{T}(#pi_{matched}) [GeV];p_{T}(K_{matched}) [GeV]",150,0,15, 150,0,15) ;
  //
  // RECO variables
  nB0_h = new TH1I("B0_multiplicity", "B0 multiplicity;B^{0} multiplicity;Events", 150, -0.5, 149.5) ;
  dEdxMass_h[0] = new TH1F("dEdxMass_h", "dE/dx;dE/dx mass [GeV]", 200, 0, 10) ;
  dEdxMass_h[1] = new TH1F("dEdxMass_sig_h", "dE/dx mass after Nairit cuts;dE/dx mass [GeV]", 200, 0, 10) ;
  dEdxMass_vs_pT = new TH2F("dEdxMass_vs_pT", "dE/dx mass vs p_{T};p_{T} [GeV];dE/dx mass [GeV]", 320, 0, 8, 200, 0, 10) ;
  dEdx_vs_p = new TH2F("dEdx_vs_p", "dE/dx vs p;p [GeV];dE/dx [MeV/cm]", 320, 0, 8, 200, 0, 10) ;
  dEdx_vs_pT[0] = new TH2F("dEdx_vs_pT", "dE/dx vs p_{T};p_{T} [GeV];dE/dx [MeV/cm]", 320, 0, 8, 200, 0, 10) ;
  dEdx_vs_pT[1] = new TH2F("dEdx_vs_pT_sig", "dE/dx vs p_{T} in signal region;p_{T} [GeV];dE/dx [MeV/cm]", 320, 0, 8, 200, 0, 10) ;
  dEdx_vs_pT[2] = new TH2F("dEdx_vs_pT_sb", "dE/dx vs p_{T} in sideband;p_{T} [GeV];dE/dx [MeV/cm]", 320, 0, 8, 200, 0, 10) ;
  //
  B0_Mass_h = new TH1F("B0_mass", "B0 mass", 200, 4, 6) ;
  myB0_Mass_h = new TH1F("myB0_mass", "myB0 mass", 200, 4, 6) ;
  Z_VtxCL_h = new TH1F("Z_vertexCL", "Z vertex CL;vtxCL;Events", 100, 0, 0.1) ;
  Z_CTauBS_h = new TH1F("Z_ctauBS", "Z c*#tau BS;c*#tau_{BS}", 1000, 0.001, 0.1) ;
  Z_CTauBSE_h = new TH1F("Z_ctauBSE", "Z c*#tau BS error;#sigma(c*#tau_{BS})", 1000, 0.001, 0.1) ;
  Z_CTauBS_over_Z_CTauBSE_h = new TH1F("Z_ctauBS_over_E", "Z c*#tau BS over error;c*#tau_{BS} / #sigma(c*#tau_{BS})", 100, 0, 10) ;
  Z_pT_h = new TH1F("Z_pT", "Z p_{T};p_{T} [GeV]", 150, 0, 30) ;
  Z_pointingAngle[0] = new TH1F("Z_cosAlpha_BC", "Z cos(#alpha) before cuts;cos(#alpha)", 100, 0.999, 1.) ; // 840, -1.05, 1.05) ;
  Z_pointingAngle[1] = new TH1F("Z_cosAlpha_sign", "Z cos(#alpha) after cuts;cos(#alpha)", 100, 0.999, 1.) ; // 840, -1.05, 1.05) ;
  Z_pointingAngle[2] = new TH1F("Z_cosAlpha_sb", "Z cos(#alpha);cos(#alpha)", 100, 0.999, 1.) ; // 840, -1.05, 1.05) ;
  Z_vs_pi_pT_h = new TH2F("Z_vs_pi_pT", "Z vs #pi^{+} p_{T};p_{T}(#pi^{+}) [GeV];p_{T}(Z) [GeV]", 200, 0, 40, 400, 0, 80) ;
  //Z_vs_pi_pT_h[0] = new TH2F("Z_vs_pi_pT", "Z vs #pi^{+} p_{T}", 200, 0, 40, 400, 0, 80) ;
  //Z_vs_pi_pT_h[1] = new TH2F("Z_vs_pi_pT_pr", "Z vs #pi^{+} p_{T} prompt", 200, 0, 40, 400, 0, 80) ;
  //Z_vs_pi_pT_h[2] = new TH2F("Z_vs_pi_pT_quasipr", "Z vs #pi^{+} p_{T} quasi prompt", 200, 0, 40, 400, 0, 80) ;
  //Z_vs_pi_pT_h[3] = new TH2F("Z_vs_pi_pT_nonpr", "Z vs #pi^{+} p_{T} non prompt", 200, 0, 40, 400, 0, 80) ;
  // psi(2S) variables
  psi2S_VtxCL_h = new TH1F("psi2S_vertexCL", "#psi' vertex CL;vtxCL(#psi')", 100, 0, 0.1) ;
  psi2S_VtxCL_HLT5_h = new TH1F("psi2S_vertexCL_HLT5", "#psi' vertex CL HLT5;vtxCL(#psi')", 100, 0, 0.1) ;
  psi2S_VtxCL_1_h = new TH1F("psi2S_vertexCL_1", "first #psi' vertex CL;vtxCL(1^{st}#psi')", 250, 0, 1) ;
  psi2S_VtxCL_2_h = new TH1F("psi2S_vertexCL_2", "second #psi' vertex CL;vtxCL(2^{nd}#psi')", 250, 0, 1) ;
  psi2S_VtxCL_diff_h = new TH1F("psi2S_vertexCL_diff", "|1^{st} - 2^{nd}| #psi' vertex CL;vtxCL(1^{st}#psi') - vtxCL(2^{nd}#psi')", 250, 0, 1) ;
  psi2SMuonTrigMatch_h = new TH1I("psi2S_MuonTrigMatch", "#psi' trigger matching", 6, -2.5, 3.5) ;
  psi2S_y_h = new TH1F("psi2S_y", "#psi' y;y", 66, -3.3, 3.3) ;
  psi2S_HLT5_yMu_h = new TH1F("psi2S_HLT5_yMu", "#psi' HLT5 mass for y(#mu^{#pm}) > 2.4;m(#mu^{+}#mu^{-}) [GeV];Candidates / 5 MeV", MuMuMass_bins, MuMuMass_low, MuMuMass_high);
  psi2S_pT_h = new TH1F("psi2S_pT", "#psi' p_{T};p_{T} [GeV]", 200, 0, 20) ;
  psi2S_2011_pT.push_back( 7. ) ; psi2S_2011_pT.push_back( 9. ) ; psi2S_2011_pT.push_back( 11. ) ;
  psi2S_2011v_pT_h.resize(psi2S_2011_pT.size()) ; myInclusiveMuMuMass_2011_HLT.resize(psi2S_2011_pT.size()) ;
  HLT_v.resize(psi2S_2011_pT.size()) ;
  HLT_v[0].push_back( "v1" ) ;  HLT_v[0].push_back( "v2" ) ; HLT_v[0].push_back( "v1v2" ) ;  HLT_v[0].push_back( "v3" ) ;  HLT_v[0].push_back( "v4" ) ;  HLT_v[0].push_back( "v5" ) ;
  HLT_v[1].push_back( "v1" ) ;  HLT_v[1].push_back( "v4" ) ;  HLT_v[1].push_back( "v5" ) ;
  HLT_v[2].push_back( "v1" ) ;  HLT_v[2].push_back( "v4" ) ;  HLT_v[2].push_back( "v5" ) ;

  for (Int_t i = 0; i < abs(psi2S_2011_pT.size()); ++i) {
    TString name = TString::Format("psi2S_2011_HLT_%.0f", psi2S_2011_pT[i]) ; 
    TString title = TString::Format("HLT_Dimuon%.0f_PsiPrime;p_{T} [GeV]", psi2S_2011_pT[i]) ;
    psi2S_2011_pT_h.push_back( new TH1F(name, title, 200, 0, 20) ) ;

    for (Int_t j = 0; j < abs(HLT_v[i].size()); ++j) {
      TString namev = TString::Format("psi2S_2011_HLT_%.0f_%s", psi2S_2011_pT[i], HLT_v[i][j].Data()) ; 
      TString titlev = TString::Format("HLT2011_Dimuon%.0f_PsiPrime_%s;p_{T} [GeV]", psi2S_2011_pT[i], HLT_v[i][j].Data()) ;
      psi2S_2011v_pT_h[i].push_back( new TH1F(namev, titlev, 200, 0, 20) ) ;

      TString nameMass = TString::Format("myInclusiveMuMuMass_2011_HLT_%.0f_%s", psi2S_2011_pT[i], HLT_v[i][j].Data() ) ; 
      TString titleMass = TString::Format("HLT2011_Dimuon%.0f_PsiPrime_%s;m(#mu^{+}#mu^{-}) [GeV]", psi2S_2011_pT[i], HLT_v[i][j].Data() ) ;
      myInclusiveMuMuMass_2011_HLT[i].push_back( new TH1F(nameMass, titleMass, MuMuMass_bins, MuMuMass_low, MuMuMass_high) ) ;
    }
  }

  /*
  psi2S_DeltaR_pi_h = new TH1F("psi2S_DeltaR_pi", "#DeltaR(#psi', #pi^{+})", 160, 0, 1.6) ;
  psi2S_DeltaR_pi_pr_h = new TH1F("psi2S_DeltaR_pi_pr", "#DeltaR(#psi', #pi^{+}) prompt", 160, 0, 1.6) ;
  psi2S_DeltaR_pi_quasipr_h = new TH1F("psi2S_DeltaR_pi_quasipr", "#DeltaR(#psi', #pi^{+}) quasi prompt", 160, 0, 1.6) ;
  psi2S_DeltaR_pi_nonpr_h = new TH1F("psi2S_DeltaR_pi_nonpr", "#DeltaR(#psi', #pi^{+}) non prompt", 160, 0, 1.6) ;
  */
  //psi2S_pointingAngle = new TH1F("psi2S_cosAlpha", "#psi' cos(#alpha)", 210, -1.05, 1.05) ; // empty variable
  //psi2S_pointingAngle_sign = new TH1F("psi2S_cosAlpha_sign", "#psi' cos(#alpha)", 210, -1.05, 1.05) ; // empty variable
  //psi2S_pointingAngle_sb = new TH1F("psi2S_cosAlpha_sb", "#psi' cos(#alpha)", 210, -1.05, 1.05) ; // empty variable
  pT_vs_y_psi2S_h = new TH2F("psi2S_pT_vs_y", "#psi' p_{T} vs y;y;p_{T} [GeV]", 88, -2.2, 2.2, 200, 0, 100) ;
  nPsi2S_h  = new TH1I("psi2S_multiplicity", "#psi' multiplicity;Candidates/Event;Events", 5, 0.5, 5.5) ;
  // muons variables
  mu_px = new TH1F("mu_px", "#mu p_{x};p_{x} [GeV]", 200, -10, 10) ;
  mu_px_psi2S = new TH1F("mu_px_psi2S", "#mu p_{x};p_{x} [GeV]", 200, -10, 10) ;
  fitmu_px = new TH1F("fitmu_px", "fit #mu p_{x};p_{x} [GeV]", 200, -10, 10) ;
  mu1_pixelHits_h = new TH1I("mu1_pixelHits", "#mu_{1} pixel hits;pixel hits", 10, 0, 10) ;
  mu2_pixelHits_h = new TH1I("mu2_pixelHits", "#mu_{2} pixel hits;pixel hits", 10, 0, 10) ;
  mu1_stripHits_h = new TH1I("mu1_stripHits", "#mu_{1} strip hits;strip hits", 25, 0, 25) ;
  mu2_stripHits_h = new TH1I("mu2_stripHits", "#mu_{2} strip hits;strip hits", 25, 0, 25) ;
  mu1_charge_h = new TH1F("mu1_charge", "#mu_{1} charge;charge", 7, -1.5, 2) ;
  mu2_charge_h = new TH1F("mu2_charge", "#mu_{2} charge;charge", 7, -1.5, 2) ;
  mu1_Chi2N_h = new TH1F("mu1_Chi2N", "#mu_{1} #chi^{2}_{N};#chi^{2}_{N}", 40, 0, 20) ;
  mu2_Chi2N_h = new TH1F("mu2_Chi2N", "#mu_{2} #chi^{2}_{N};#chi^{2}_{N}", 40, 0, 20) ;
  mu1Tr_Chi2N_h  = new TH1F("mu1Tr_Chi2N", "#mu_{1} track #chi^{2}_{N};#chi^{2}_{N}", 40, 0, 20) ;
  mu2Tr_Chi2N_h  = new TH1F("mu2Tr_Chi2N", "#mu_{2} track #chi^{2}_{N};#chi^{2}_{N}", 40, 0, 20) ;
  mu1GlTr_Chi2N_h  = new TH1F("mu1GlTr_Chi2N", "#mu_{1} globalTrack #chi^{2}_{N};#chi^{2}_{N}", 40, 0, 20) ;
  mu2GlTr_Chi2N_h  = new TH1F("mu2GlTr_Chi2N", "#mu_{2} globalTrack #chi^{2}_{N};#chi^{2}_{N}", 40, 0, 20) ;
  mu1_Dz_h = new TH1F("mu1_Dz", "#mu_{1} d_{z};d_{z} [cm]", 100, -50, 50) ;
  mu2_Dz_h = new TH1F("mu2_Dz", "#mu_{2} d_{z};d_{z} [cm]", 100, -50, 50) ;
  mu1_Dxy_h = new TH1F("mu1_Dxy", "#mu_{1} d_{xy};d_{xy} [cm]", 200, -10, 10) ;
  mu2_Dxy_h = new TH1F("mu2_Dxy", "#mu_{2} d_{xy};d_{xy} [cm]", 200, -10, 10) ;
  mu1_y_h = new TH1F("mu1_y", "#mu_{1} y;y;events", 66, -3.3, 3.3) ;
  mu2_y_h = new TH1F("mu2_y", "#mu_{2} y;y;events", 66, -3.3, 3.3) ;
  mu1_pT_h = new TH1F("mu1_pT", "#mu_{1} p_{T};p_{T} [GeV];events", 200, 0, 20) ;
  mu2_pT_h = new TH1F("mu2_pT", "#mu_{2} p_{T};p_{T} [GeV];events", 200, 0, 20) ;
  mu1_pT_test_h = new TH1F("mu1_pT_test", "#mu_{1} p_{T};p_{T} [GeV];events", 200, 0, 20) ;
  mu2_pT_test_h = new TH1F("mu2_pT_test", "#mu_{2} p_{T};p_{T} [GeV];events", 200, 0, 20) ;
  mu1_eta_test_h = new TH1F("mu1_eta_test", "#mu_{1} #eta;#eta;events", 72, -2.4, 2.4) ;
  mu2_eta_test_h = new TH1F("mu2_eta_test", "#mu_{2} #eta;#eta;events", 72, -2.4, 2.4) ;
  mu1_DeltaR_mu2_h = new TH1F("mu1_DeltaR_mu2", "#DeltaR(#mu_{1}, #mu_{2})", 160, 0, 1.6) ;
  //mu2_DeltaR_mu1_h = new TH1F("mu2_DeltaR_mu1", "#DeltaR(#mu_{2}, #mu_{1})", 160, 0, 1.6) ;
  //mu1_vs_mu2_pT_h = new TH2F("mu1_vs_mu2_pT", "#mu_{1} vs #mu_{2} p_{T}", 200, 0, 20, 200, 0, 20) ;
  pT_vs_y_mu1_h = new TH2F("mu1_pT_vs_y", "#mu_{1} p_{T} vs y;y;p_{T} [GeV]", 66, -3.3, 3.3, 200, 0, 20) ;
  pT_vs_y_mu2_h = new TH2F("mu2_pT_vs_y", "#mu_{2} p_{T} vs y;y;p_{T} [GeV]", 66, -3.3, 3.3, 200, 0, 20) ;
  pT10_vs_eta_mu_h = new TH2F("mu_pT10_vs_eta", "#mu p_{T} vs #eta;|#eta|;p_{T} [GeV]", 50, 0., 2.5, 100, 0, 10) ;
  pT_vs_eta_mu_h = new TH2F("mu_pT_vs_eta", "#mu p_{T} vs #eta;|#eta|;p_{T} [GeV]", 50, 0., 2.5, 300, 0, 30) ;
  // pi+ variables
  pi_pixelHits_h = new TH1I("pi_pixelHits", "#pi^{+} pixel hits;pixel hits", 10, 0, 10) ;
  pi_stripHits_h = new TH1I("pi_stripHits", "#pi^{+} strip hits;strip hits", 25, 0, 25) ;
  pi_Chi2N_h = new TH1F("pi_Chi2N", "#pi^{+} #chi^{2}_{N};#chi^{2}_{N}", 80, 0, 4) ;
  pi_pT_h = new TH1F("pi_pT", "#pi^{+} p_{T};p_{T} [GeV]", 150, 0, 15) ;
  pi_y_h = new TH1F("pi_y", "#pi^{+} y;y", 66, -3.3, 3.3) ;
  pi_Chi2N_test_h = new TH1F("pi_Chi2N_test", "#pi^{+} #chi^{2}_{N};#chi^{2}_{N}", 80, 0, 4) ;
  pi_pT_test_h = new TH1F("pi_pT_test", "#pi^{+} p_{T};p_{T} [GeV]", 150, 0, 15) ;
  pi_y_test_h = new TH1F("pi_y_test", "#pi^{+} y;y", 66, -3.3, 3.3) ;
  pT_vs_pixelHits_pi_h = new TH2F("pT_vs_pixelHits_pi", "p_{T} vs pixelHits #pi^{+};pixel hits;p_{T} [GeV]", 10, 0, 10, 30, 0, 3) ;
  pT_vs_stripHits_pi_h = new TH2F("pT_vs_stripHits_pi", "p_{T} vs stripHits #pi^{+};strip hits;p_{T} [GeV]", 25, 0, 25, 30, 0, 3) ;
  pT_vs_y_pi_h = new TH2F("pT_vs_y_pi", "p_{T} vs y #pi^{+};y;p_{T} [GeV]", 66, -3.3, 3.3, 150, 0, 15) ;
  pi_Dz_h = new TH1F("pi_Dz", "#pi^{+} d_{z};d_{z} [cm]", 20, -1, 1) ;
  pi_Dxy_h = new TH1F("pi_Dxy", "#pi^{+} d_{xy};d_{xy} [cm]", 200, -1, 1) ;
  pi_Dz_pr_h = new TH1F("pi_Dz_pr", "#pi^{+} d_{z} prompt;d_{z} [cm]", 20, -1, 1) ;
  pi_Dxy_pr_h = new TH1F("pi_Dxy_pr", "#pi^{+} d_{xy} prompt;d_{xy} [cm]", 200, -1, 1) ;
  pi_Dz_nonpr_h = new TH1F("pi_Dz_nonpr", "#pi^{+} d_{z} non prompt;d_{z} [cm]", 20, -1, 1) ; 
  pi_Dxy_nonpr_h = new TH1F("pi_Dxy_nonpr", "#pi^{+} d_{xy} non prompt;d_{xy} [cm]", 200, -1, 1) ;
  
   //////////////////////////////////////////////////////////////////
}

Bool_t psiPrimePiK_MC::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either psiPrimePiK_MC::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  GetEntry(entry);
  //if (entry % 10000 == 0) cout << "I am running " << entry << "th entries out of " << nentries << " total entries" << endl; //not with PROOF
  /*
  cout <<"running Process with MC set to = " <<(MC ? "true" : "false") <<endl ;
  cout <<"running Process with officialMC set to = " <<(officialMC ? "true" : "false") <<endl ;
  cout <<"running Process with oldNtuple set to = " <<(oldNtuple ? "true" : "false") <<endl ;
  cout <<"running Process with priVtxsInfo = " <<(priVtxsInfo ? "true" : "false") <<endl ;
  */
  //////////////////////////////// put the code here

  //////////////////////////////////////// EXCLUSIVE KPi generated
  if ( MC ) {
    nMCB0_BT_h->Fill( nMCB0 ) ;
    for (Int_t myMCB0Idx = 0; myMCB0Idx < abs(nMCB0); myMCB0Idx++) { // only 0 and 1 candidate
      TLorentzVector B0p4_gen, Jp4_gen, muPp4_gen, muMp4_gen, Kp4_gen, Pip4_gen, Kstar_gen ;
      //      
      Float_t mumu_mass = psi2S_mass; 
      // for B+ -> J/psi K+ test or B0 -> J/psi K*
      mumu_mass = jpsi_mass; 
      //
      Float_t K_E_gen = sqrt( pow((*MCkaonPx)[myMCB0Idx], 2) + pow((*MCkaonPy)[myMCB0Idx], 2) + pow((*MCkaonPz)[myMCB0Idx], 2) + pow(kaonCh_mass, 2) ) ;
      Float_t Pi_E_gen = sqrt( pow((*MCpionPx)[myMCB0Idx], 2) + pow((*MCpionPy)[myMCB0Idx], 2) + pow((*MCpionPz)[myMCB0Idx], 2) + pow(pionCh_mass, 2) ) ;
      Float_t muP_E_gen = sqrt( pow((*MCmupPx)[myMCB0Idx], 2) + pow((*MCmupPy)[myMCB0Idx], 2) + pow((*MCmupPz)[myMCB0Idx], 2) + pow(muon_mass, 2) ) ;
      Float_t muM_E_gen = sqrt( pow((*MCmumPx)[myMCB0Idx], 2) + pow((*MCmumPy)[myMCB0Idx], 2) + pow((*MCmumPz)[myMCB0Idx], 2) + pow(muon_mass, 2) ) ;
      //Float_t J_E_gen = sqrt( pow((*MCpsi2SPx)[myMCB0Idx], 2) + pow((*MCpsi2SPy)[myMCB0Idx], 2) + pow((*MCpsi2SPz)[myMCB0Idx], 2) + pow(mumu_mass, 2) ) ;
      //
      muPp4_gen.SetPxPyPzE( (*MCmupPx)[myMCB0Idx], (*MCmupPy)[myMCB0Idx], (*MCmupPz)[myMCB0Idx], muP_E_gen ) ;
      muMp4_gen.SetPxPyPzE( (*MCmumPx)[myMCB0Idx], (*MCmumPy)[myMCB0Idx], (*MCmumPz)[myMCB0Idx], muM_E_gen ) ;
      //Jp4_gen.SetPxPyPzE( (*MCpsi2SPx)[myMCB0Idx], (*MCpsi2SPy)[myMCB0Idx], (*MCpsi2SPz)[myMCB0Idx], J_E_gen ) ; replaced by:
      Jp4_gen = muPp4_gen + muMp4_gen;
      Kp4_gen.SetPxPyPzE( (*MCkaonPx)[myMCB0Idx], (*MCkaonPy)[myMCB0Idx], (*MCkaonPz)[myMCB0Idx], K_E_gen ) ;
      Pip4_gen.SetPxPyPzE( (*MCpionPx)[myMCB0Idx], (*MCpionPy)[myMCB0Idx], (*MCpionPz)[myMCB0Idx], Pi_E_gen ) ;
      B0p4_gen = Jp4_gen + Pip4_gen + Kp4_gen ;
      //B0p4_gen = Jp4_gen + Kp4_gen ; // for B+ -> J/psi K+ test
      Kstar_gen = Pip4_gen + Kp4_gen ;
      /*
      B0_gen_p_h->Fill( B0p4_gen.P() ) ; B0_gen_pT_h->Fill( B0p4_gen.Pt() ) ;
      psi2S_gen_p_h->Fill( Jp4_gen.P() ) ; psi2S_gen_pT_h->Fill( Jp4_gen.Pt() ) ;
      K_gen_p_h->Fill( Kp4_gen.P() ) ; K_gen_pT_h->Fill( Kp4_gen.Pt() ) ;
      pi_gen_p_h->Fill( Pip4_gen.P() ) ; pi_gen_pT_h->Fill( Pip4_gen.Pt() ) ;
      */
      //mumuMass->Fill( (muPp4_gen + muMp4_gen).M() ) ;
      psi2SPiSqMass_vs_KPiSqMass_BT_gen->Fill( (Kp4_gen + Pip4_gen).M2(), (Jp4_gen + Pip4_gen).M2() );
      KPiSqMass_vs_psi2SKSqMass_BT_gen->Fill( (Jp4_gen + Kp4_gen).M2(), (Kp4_gen + Pip4_gen).M2() );

      // gen angle variables
      psi2S_vs_mu_pT->Fill(Jp4_gen.Pt(), TMath::Min(muPp4_gen.Pt(),muMp4_gen.Pt())) ;
      //
      // psi helicity angle
      //Double_t theta_psi_gen = getMumuHelicityAngle(B0p4_gen, Kstar_gen, Jp4_gen, muMp4_gen);
      //Double_t theta_psi_gen = getAleMumuHelicityAngle(B0p4_gen, Kstar_gen, Jp4_gen, muMp4_gen);
      Double_t theta_psi_gen = GetThetaMuMu(B0p4_gen, Jp4_gen, muPp4_gen, beam_energy, mumu_mass, muon_mass) ;
      psi2S_helicityAngle_BT_gen->Fill( theta_psi_gen ) ; 
      cos_psi2S_helicityAngle_BT_gen->Fill( TMath::Cos(theta_psi_gen) ) ;
      Double_t theta_Kstar_gen = GetThetaMuMu(B0p4_gen, Kstar_gen, (*MCpionCh)[myMCB0Idx]>0 ? Pip4_gen : Kp4_gen, beam_energy, Kstar0_mass, (*MCpionCh)[myMCB0Idx]>0 ? pionCh_mass : kaonCh_mass) ;
      Double_t theta_psi_gen_rot = GetThetaMuMu_rot(B0p4_gen, Jp4_gen, muPp4_gen, beam_energy, mumu_mass, muon_mass) ;
      cos_psi2S_helicityAngle_BT_gen_rot->Fill( TMath::Cos(theta_psi_gen_rot) ) ; 
      cos_psi2S_helicityAngle_BT_gen_rotVSnotRot->Fill( TMath::Cos(theta_psi_gen), TMath::Cos(theta_psi_gen_rot) ) ;
      //
      cos_psi2S_helicityAngle_vs_psi2SPiMass->Fill((Jp4_gen + Pip4_gen).M(), TMath::Cos(theta_psi_gen)) ;
      cos_psi2S_helicityAngle_vs_KPiMass->Fill((Kp4_gen + Pip4_gen).M(), TMath::Cos(theta_psi_gen)) ;
      // planes angle
      Double_t phi_gen = getPlanesAngle(B0p4_gen, Kp4_gen, Pip4_gen, muMp4_gen, muPp4_gen) ;
      //Double_t phi_gen = getBellePlanesAngle(B0p4_gen, Kp4_gen, Pip4_gen, muMp4_gen, muPp4_gen) ;
      planesAngle_BT_gen->Fill( phi_gen ) ;

      // respect TNtuple structure: "runNum : evtNum : B0_mass : KPi_mass : Jpsi_mass : JpsiPi_mass : B0_3mom : KPi_3mom : Jpsi_3mom : K_3mom : Pi_3mom : Theta_Kstar : Theta_Jpsi : Phi"
      Double_t NT_var[] = {runNum, evtNum, B0p4_gen.M(), Kstar_gen.M(), Jp4_gen.M(), (Jp4_gen + Pip4_gen).M(), B0p4_gen.P(), Kstar_gen.P(), Jp4_gen.P(), Kp4_gen.P(), Pip4_gen.P(), theta_Kstar_gen, theta_psi_gen, phi_gen};
      _nt->Fill( NT_var );
    }
  }


  Bool_t year2011 = false, year2012 = false;
  /////////////////////// HLT
  // 2011
  Bool_t HLT_2011_7_v1 = false,  HLT_2011_7_v2 = false,  HLT_2011_7_v1v2 = false, HLT_2011_7_v3 = false, HLT_2011_7_v4 = false, HLT_2011_7_v5 = false ;
  Bool_t HLT_2011_9_v1 = false,  HLT_2011_9_v4 = false,  HLT_2011_9_v5 = false ;
  Bool_t HLT_2011_11_v1 = false, HLT_2011_11_v4 = false, HLT_2011_11_v5 = false ;
  //
  Bool_t HLT_2011_Dimuon[] = {false, false, false} ;
  Bool_t HLT_2011 = false ;
  //
  // 2012
  Bool_t HLT_0_v3 = false, HLT_0_v4 = false, HLT_0_v5 = false, HLT_0_v6 = false ;
  Bool_t HLT_5_v3 = false, HLT_5_v4 = false, HLT_5_v5 = false, HLT_5_v6 = false ;
  Bool_t HLT_7_v1 = false, HLT_7_v2 = false, HLT_7_v3 = false ;
  Bool_t HLT_9_v9 = false, HLT_7_v1_9_v9 = false ;
  //
  Bool_t HLT_3p5_v3 = false, HLT_3p5_v4 = false, HLT_3p5_v5 = false, HLT_3p5_v6 = false ;
  //
  Bool_t HLT_5_vAll = false, HLT_7_vAll = false, HLT_3p5_vAll = false ;
  

  for (Int_t i = 0; i != abs(TrigRes->size()); ++i)
    {
      // 2011
      if ( MC || (runNum >= 165088  &&  runNum <= 180252) ) {
	year2011 = true ;
	if (!HLT_Displaced) {
	  if ( TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v1") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_7_v1 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v2") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_7_v2 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v3") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_7_v3 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v4") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_7_v4 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v5") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_7_v5 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon9_PsiPrime_v1") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_9_v1 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon9_PsiPrime_v4") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_9_v4 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon9_PsiPrime_v5") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_9_v5 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon11_PsiPrime_v1") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_11_v1 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon11_PsiPrime_v4") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_11_v4 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon11_PsiPrime_v5") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_2011_11_v5 = true;
	}
      }
      // 2012
      if ( MC || (runNum >= 190645  &&  runNum <= 208686) ) {
	year2012 = true ;
	if (!HLT_Displaced) {
	  //--0
	  if ( TrigNames->at(i).find("HLT_Dimuon0_PsiPrime_v3") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_0_v3 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon0_PsiPrime_v4") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_0_v4 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon0_PsiPrime_v5") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_0_v5 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon0_PsiPrime_v6") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_0_v6 = true;
	  //
	  //-- 5
	  if ( TrigNames->at(i).find("HLT_Dimuon5_PsiPrime_v3") != string::npos  &&  TrigRes->at(i) == 1) HLT_5_v3 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon5_PsiPrime_v4") != string::npos  &&  TrigRes->at(i) == 1) HLT_5_v4 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon5_PsiPrime_v5") != string::npos  &&  TrigRes->at(i) == 1) HLT_5_v5 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon5_PsiPrime_v6") != string::npos  &&  TrigRes->at(i) == 1) HLT_5_v6 = true;
	  //
	  //- 7/9
	  if ( TrigNames->at(i).find("HLT_Dimuon9_PsiPrime_v9") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_9_v9 = true;
	  if ( (TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v1") != string::npos  &&  TrigRes->at(i) == 1 ) ||
	       (TrigNames->at(i).find("HLT_Dimuon9_PsiPrime_v9") != string::npos  &&  TrigRes->at(i) == 1) ) HLT_7_v1_9_v9 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v1") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_7_v1 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v2") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_7_v2 = true;
	  if ( TrigNames->at(i).find("HLT_Dimuon7_PsiPrime_v3") != string::npos  &&  TrigRes->at(i) == 1 ) HLT_7_v3 = true;
	  //
	} else {
	  //- 3.5
	  if ( TrigNames->at(i).find("HLT_DoubleMu3p5_LowMass_Displaced_v3") != string::npos  &&  TrigRes->at(i) == 1) HLT_3p5_v3 = true;
	  if ( TrigNames->at(i).find("HLT_DoubleMu3p5_LowMass_Displaced_v4") != string::npos  &&  TrigRes->at(i) == 1) HLT_3p5_v4 = true;
	  if ( TrigNames->at(i).find("HLT_DoubleMu3p5_LowMass_Displaced_v5") != string::npos  &&  TrigRes->at(i) == 1) HLT_3p5_v5 = true;
	  if ( TrigNames->at(i).find("HLT_DoubleMu3p5_LowMass_Displaced_v6") != string::npos  &&  TrigRes->at(i) == 1) HLT_3p5_v6 = true;
	}
      }
    } // endl loop over TrigRes

  if (HLT_2011_7_v1 || HLT_2011_7_v2) HLT_2011_7_v1v2 = true;
  if (HLT_2011_7_v1 || HLT_2011_7_v2 || HLT_2011_7_v3 || HLT_2011_7_v4 || HLT_2011_7_v5) HLT_2011_Dimuon[0] = true;
  if (HLT_2011_9_v1 || HLT_2011_9_v4 || HLT_2011_9_v5) HLT_2011_Dimuon[1] = true;
  if (HLT_2011_11_v1 || HLT_2011_11_v4 || HLT_2011_11_v5) HLT_2011_Dimuon[2] = true;
  if (Kai_cuts)
    if ( HLT_2011_Dimuon[0] || HLT_2011_Dimuon[1] || HLT_2011_Dimuon[2] ) HLT_2011 = true;
  
  vector< vector<Bool_t> > HLT_2011_Dimuon_v (3) ;
  HLT_2011_Dimuon_v[0].push_back(HLT_2011_7_v1) ; HLT_2011_Dimuon_v[0].push_back(HLT_2011_7_v2) ; HLT_2011_Dimuon_v[0].push_back(HLT_2011_7_v1v2) ; HLT_2011_Dimuon_v[0].push_back(HLT_2011_7_v3) ; HLT_2011_Dimuon_v[0].push_back(HLT_2011_7_v4) ; HLT_2011_Dimuon_v[0].push_back(HLT_2011_7_v5) ; 
  HLT_2011_Dimuon_v[1].push_back(HLT_2011_9_v1) ; HLT_2011_Dimuon_v[1].push_back(HLT_2011_9_v4) ; HLT_2011_Dimuon_v[1].push_back(HLT_2011_9_v5) ; 
  HLT_2011_Dimuon_v[2].push_back(HLT_2011_11_v1) ; HLT_2011_Dimuon_v[2].push_back(HLT_2011_11_v4) ; HLT_2011_Dimuon_v[2].push_back(HLT_2011_11_v5) ; 

  if (HLT_5_v3 || HLT_5_v4 || HLT_5_v5 || HLT_5_v6) HLT_5_vAll = true;
  if (HLT_7_v1 || HLT_7_v2 || HLT_7_v3) HLT_7_vAll = true;
  //
  if (HLT_3p5_v3 || HLT_3p5_v4 || HLT_3p5_v5 || HLT_3p5_v6) HLT_3p5_vAll = true;

  //
  /////////////////////////////////
  for (Int_t myJIdx=0; myJIdx<abs(nMuMu); myJIdx++) 
    myInclusiveMuMuMass_noTrigg->Fill( (*MuMuMass)[myJIdx] ) ;

  //if (!HLT_5_vAll) return kTRUE ;
  
  int muonQual[] = {1,3,4,12};
  Bool_t isNewSoftMuons = kTRUE;

  Float_t fit_mass = 0, fit_sigma = 0 ;
  //Float_t mass_low=0, mass_high=0 ;
  // signal and sidebands definition
  Float_t signal_win = 2. ; // signal_win = 2.5 ;  signal_win = 5. ;  // remember to rename the output file
  Float_t sb_low = 4. ;
  Float_t sb_high = sb_low + 1*signal_win ;
  
  Int_t n_psiPrime = 0 ;
  Float_t first_psi2S_VtxCL = 0., second_psi2S_VtxCL = 0. ;
  Float_t first_psi2S_mass = 0., second_psi2S_mass = 0. ;

  Int_t first_psi2S_index = 0 ;
  Float_t first_psi2S_y = 0 ;

  // with newSoftMuonsID
  /*  // 2011 without pT(psi') > 7 cut
      Float_t mean_2011[] = {3.67833, 3.67738, 3.67772, 3.67811, 3.67949, 3.68135, 3.681468, 3.68224, 3.68249, 3.68266, 3.68298, 3.68272, 3.68209, 3.68246, 3.68168, 3.68114, 3.67979, 3.68023, 3.67879, 3.67795, 3.68050, 3.6757} ; // 22 elements from bin 1 to 23
      Float_t sigma_2011[] = {0.06412, 0.05474, 0.04943, 0.04732, 0.04368, 0.03981, 0.035095, 0.02943, 0.02572, 0.02371, 0.02138, 0.02152, 0.02377, 0.02547, 0.03038, 0.03578, 0.04034, 0.04381, 0.04696, 0.04798, 0.0550, 0.0637} ; // 22 elements from bin 1 to 23
  */
  // 2011 with pT(psi') > 7 cut
  Float_t mean_2011[] = {3.6784, 3.67749, 3.67776, 3.67814, 3.67952, 3.68127, 3.68161, 3.68223, 3.68244, 3.68266, 3.68296, 3.68271, 3.6821, 3.68246, 3.68165, 3.68116, 3.67976, 3.68033, 3.67900, 3.67782, 3.68066, 3.6758} ; // 22 elements from bin 1 to 23
  Float_t sigma_2011[] = {0.0638, 0.05458, 0.04922, 0.04744, 0.04356, 0.03974, 0.03499, 0.02942, 0.02575, 0.02371, 0.02138, 0.02153, 0.0237, 0.02547, 0.03038, 0.03580, 0.04037, 0.04375, 0.04693, 0.04809, 0.05478, 0.0640} ; // 22 elements from bin 1 to 23
  // 2012
  Float_t mean[] = {3.6742, 3.67536, 3.67607, 3.67600, 3.67811, 3.67958, 3.68074, 3.68218, 3.68180, 3.68190, 3.68224, 3.68197, 3.68180, 3.68166, 3.68188, 3.68067, 3.67915, 3.67836, 3.67590, 3.67692, 3.67597, 3.6751} ; // 22 elements from bin 1 to 23
  Float_t sigma[] = {0.0578, 0.05019, 0.04817, 0.04648, 0.04359, 0.03976, 0.03516, 0.02963, 0.02679, 0.02418, 0.02259, 0.02274, 0.02452, 0.02658, 0.02987, 0.03506, 0.03965, 0.04443, 0.04667, 0.04723, 0.05099, 0.0512} ; // 22 elements from bin 1 to 23

  for (Int_t myJIdx = 0; myJIdx < abs(nMuMu); myJIdx++) 
    {
      //int numGoodLooseMuon=(int) (*muIsGoodLooseMuon)[ (*mupIdx)[myJIdx] ] + (int) (*muIsGoodLooseMuon)[ (*mumIdx)[myJIdx] ];
      //int numGoodTightMuon=(int) (*muIsGoodTightMuon)[ (*mupIdx)[myJIdx] ] + (int) (*muIsGoodTightMuon)[ (*mumIdx)[myJIdx] ];   

      Int_t mupidx = (*mu1Idx)[myJIdx] ;
      Int_t mumidx = (*mu2Idx)[myJIdx] ;
    
      TLorentzVector muP, muM ;
      Float_t muP_E = 0., muM_E = 0. ;
      muP_E = sqrt( pow((*muPx)[mupidx], 2) + pow((*muPy)[mupidx], 2) + pow((*muPz)[mupidx], 2) + pow(muon_mass, 2) ) ;    
      muM_E = sqrt( pow((*muPx)[mumidx], 2) + pow((*muPy)[mumidx], 2) + pow((*muPz)[mumidx], 2) + pow(muon_mass, 2) ) ;    
      muP.SetPxPyPzE( (*muPx)[mupidx], (*muPy)[mupidx], (*muPz)[mupidx], muP_E) ;
      muM.SetPxPyPzE( (*muPx)[mumidx], (*muPy)[mumidx], (*muPz)[mumidx], muM_E) ;

      Float_t Jmass = (*MuMuMass)[myJIdx] ;
      Float_t psiPrime_E = sqrt( pow((*MuMuPx)[myJIdx], 2) + pow((*MuMuPy)[myJIdx], 2) + pow((*MuMuPz)[myJIdx], 2) + pow(Jmass, 2)) ;
      TLorentzVector psiPrime4Vec((*MuMuPx)[myJIdx], (*MuMuPy)[myJIdx], (*MuMuPz)[myJIdx], psiPrime_E);
      Float_t psiPrimeRapid = psiPrime4Vec.Rapidity();

      if ( HLT_5_vAll ) {
	myMuMuMass_noMuID->Fill( Jmass );
	if ( 1 // old SoftMuonID
	     &&  ( ((*muQual)[mupidx]) & (1 << muonQual[3]) ) && ( ((*muQual)[mumidx]) & (1 << muonQual[3]) ) // TMOneStationTight
	     &&  ((*muShits)[mupidx] >= 5) && ((*muShits)[mumidx] >= 5) // trackerLayersWithMeasurement
	     &&  ((*muPhits)[mupidx] >= 1) && ((*muPhits)[mumidx] >= 1) // pixelLayersWithMeasurement 
	     //&&  (( (*mu1Chi_MuMu2)[mupidx] / (*mu1NDF_MuMu)[mupidx] ) < 1.8)  &&  (( (*mu2Chi2_MuMu)[mumidx] / (*mu2NDF_MuMu)[mumidx] ) < 1.8) // wrong variables!
	     &&  (fabs((*muDzVtx)[mupidx]) < 30) && (fabs((*muDzVtx)[mumidx]) < 30)
	     &&  (fabs((*muDxyVtx)[mupidx]) < 3) && (fabs((*muDxyVtx)[mumidx]) < 3)
	     ) {
	  if ( ((*muChi2)[mupidx] / (*muNDF)[mupidx] ) < 1.8 )
	    myMuMuMass_Mu1ID->Fill( Jmass );
	  if ( ((*muChi2)[mumidx] / (*muNDF)[mumidx] ) < 1.8 )
	    myMuMuMass_Mu2ID->Fill( Jmass );
	  if ((*muChi2)[mupidx]/(*muNDF)[mupidx] < 3.0  &&  (*muChi2)[mumidx]/(*muNDF)[mumidx] < 3.0)
	    myMuMuMass_MuID->Fill( Jmass );
	}
      }

      if (newSoftMuonID) {
	if ( !( 1
		&&  ( ((*muQual)[mupidx]) & (1 << muonQual[3]) ) && ( ((*muQual)[mumidx]) & (1 << muonQual[3]) ) // TMOneStationTight
		&&  ((*muShits)[mupidx] > 5) && ((*muShits)[mumidx] > 5) // trackerLayersWithMeasurement
		&&  ((*muPhits)[mupidx] > 0) && ((*muPhits)[mumidx] > 0) // pixelLayersWithMeasurement 
		&&  ((*muChi2)[mupidx]/(*muNDF)[mupidx] < 3.0) && ((*muChi2)[mumidx]/(*muNDF)[mumidx] < 3.0)
		&&  (fabs((*muDzVtx)[mupidx]) < 20) && (fabs((*muDzVtx)[mumidx]) < 20)
		&&  (fabs((*muDxyVtx)[mupidx]) < 0.3) && (fabs((*muDxyVtx)[mumidx]) < 0.3)
		) )
	  isNewSoftMuons = kFALSE;
      }
      
      psi2S_VtxCL_h->Fill( (*MuMuVtx_CL)[myJIdx] ) ;
      if (HLT_5_vAll) psi2S_VtxCL_HLT5_h->Fill( (*MuMuVtx_CL)[myJIdx] ) ;

      if ( (*MuMuMuonTrigMatch)[myJIdx] ) {
	psi2SMuonTrigMatch_h->Fill( 1 ) ;
	if (HLT_5_vAll)
	  psi2SMuonTrigMatch_h->Fill( 2 ) ;
      } else {
	psi2SMuonTrigMatch_h->Fill( 0 ) ;
	if (HLT_5_vAll)
	  psi2SMuonTrigMatch_h->Fill( -1 ) ;
      }

      if ( 1
	   && (*MuMuVtx_CL)[myJIdx] > 0.005 // 0.000001 Kai
	   //&& numGoodLooseMuon>=2
	   //&& numGoodTightMuon>=2
	   ) {
	myInclusiveMuMuMass->Fill( Jmass );

	// 2011
	for (Int_t i = 0; i < abs(psi2S_2011_pT.size()); ++i)
	  for (Int_t j = 0; j < abs(HLT_v[i].size()); ++j)
	    if (HLT_2011_Dimuon_v[i][j])
	      myInclusiveMuMuMass_2011_HLT[i][j]->Fill( Jmass );

	// 2012
	// HLT_0_v*
	if (HLT_0_v3) myInclusiveMuMuMass_HLT_0_v3->Fill( Jmass );
	if (HLT_0_v4) myInclusiveMuMuMass_HLT_0_v4->Fill( Jmass );
	if (HLT_0_v5) myInclusiveMuMuMass_HLT_0_v5->Fill( Jmass );
	if (HLT_0_v6) myInclusiveMuMuMass_HLT_0_v6->Fill( Jmass );
	// HLT_5_v*
	if ( MC || (runNum >= 190645 && runNum <= 194712) )
	  if (HLT_5_v3) myInclusiveMuMuMass_HLT_5_v3->Fill( Jmass );
	if ( MC || (runNum >= 194735 && runNum <= 196531) )
	  if (HLT_5_v4) myInclusiveMuMuMass_HLT_5_v4->Fill( Jmass );
	if ( MC || (runNum >= 198049 && runNum <= 199608) )
	  if (HLT_5_v5) myInclusiveMuMuMass_HLT_5_v5->Fill( Jmass );
	if ( MC || (runNum >= 199698 && runNum <= 208686) )
	  if (HLT_5_v6) myInclusiveMuMuMass_HLT_5_v6->Fill( Jmass );
	// HLT_7_v*
	if ( MC || (runNum >= 194735 && runNum <= 196531) )
	  if ( HLT_7_v1 ) myInclusiveMuMuMass_HLT_7_v1->Fill( Jmass );
	if ( MC || (runNum >= 198049 && runNum <= 199608) )
	  if ( HLT_7_v2 ) myInclusiveMuMuMass_HLT_7_v2->Fill( Jmass );
	if ( MC || (runNum >= 199698 && runNum <= 208686) )
	  if ( HLT_7_v3 ) myInclusiveMuMuMass_HLT_7_v3->Fill( Jmass );
	// datasets
	if ( MC || (runNum >= 190456 && runNum <= 193621) ) { // 2012_A
	  myInclusiveMuMuMass_12A->Fill( Jmass );
	  if ( HLT_0_v3 ) { myInclusiveMuMuMass_12A_HLT_0_v3->Fill( Jmass );
	    myInclusiveMuMuMass_HLT_0_v3A->Fill( Jmass ); }
	  if ( HLT_5_v3 ) { myInclusiveMuMuMass_12A_HLT_5_v3->Fill( Jmass );
	    myInclusiveMuMuMass_HLT_5_v3A->Fill( Jmass ); }
	  if ( HLT_9_v9 ) myInclusiveMuMuMass_12A_HLT_9_v9->Fill( Jmass );
	} else if ( MC || (runNum >= 193834 && runNum <= 196531) ) {  // 2012_B
	  myInclusiveMuMuMass_12B->Fill( Jmass );
	  if( HLT_0_v3 ) myInclusiveMuMuMass_HLT_0_v3B->Fill( Jmass ); 
	  if( HLT_5_v3 ) myInclusiveMuMuMass_HLT_5_v3B->Fill( Jmass ); 
	  if( HLT_0_v3 || HLT_0_v4 ) { myInclusiveMuMuMass_12B_HLT_0_v3v4->Fill( Jmass );
	    myInclusiveMuMuMass_HLT_0_v3Bv4->Fill( Jmass ); }
	  if( HLT_5_v3 || HLT_5_v4 ) { myInclusiveMuMuMass_12B_HLT_5_v3v4->Fill( Jmass );
	    myInclusiveMuMuMass_HLT_5_v3Bv4->Fill( Jmass ); }
	  if( HLT_7_v1_9_v9 ) myInclusiveMuMuMass_12B_HLT_7_v1_9_v9->Fill( Jmass );
	} else if ( MC || (runNum >= 198022 && runNum <= 203742) ) {  // 2012_C
	  myInclusiveMuMuMass_12C->Fill( Jmass );
	  if( HLT_0_v5 ) myInclusiveMuMuMass_12C_HLT_0_v5->Fill( Jmass );
	  if( HLT_0_v6 ) myInclusiveMuMuMass_HLT_0_v6C->Fill( Jmass ); 
	  if( HLT_0_v5 || HLT_0_v6 ) myInclusiveMuMuMass_HLT_0_v5v6C->Fill( Jmass );
	  if( HLT_5_v5 ) myInclusiveMuMuMass_12C_HLT_5_v5->Fill( Jmass );
	  if( HLT_5_v5 || HLT_7_v2 ) myInclusiveMuMuMass_12C_HLT_5_v5_or_7_v2->Fill( Jmass );
	  if( HLT_5_v6 ) myInclusiveMuMuMass_HLT_5_v6C->Fill( Jmass ); 
	  if( HLT_5_v5 || HLT_5_v6 ) myInclusiveMuMuMass_HLT_5_v5v6C->Fill( Jmass );
	  if( HLT_7_v2 ) myInclusiveMuMuMass_12C_HLT_7_v2->Fill( Jmass );
	  if( HLT_7_v3 ) myInclusiveMuMuMass_HLT_7_v3C->Fill( Jmass );
	  if( HLT_7_v2 || HLT_7_v3 ) myInclusiveMuMuMass_HLT_7_v2v3C->Fill( Jmass );
	} else if ( MC || (runNum >= 203777 && runNum <= 208686) ) {  // 2012_D
	  myInclusiveMuMuMass_12D->Fill( Jmass );
	  if( HLT_0_v6 ) { myInclusiveMuMuMass_12D_HLT_0_v6->Fill( Jmass );
	    myInclusiveMuMuMass_HLT_0_v6D->Fill( Jmass ); }
	  if( HLT_5_v6 ) { myInclusiveMuMuMass_12D_HLT_5_v6->Fill( Jmass );
	    myInclusiveMuMuMass_HLT_5_v6D->Fill( Jmass ); }
	  if( HLT_7_v3 ) { myInclusiveMuMuMass_12D_HLT_7_v3->Fill( Jmass );
	    myInclusiveMuMuMass_HLT_7_v3D->Fill( Jmass ); }
	} else {
	  myInclusiveMuMuMass_out->Fill( Jmass ); // check ok: it's void
	}
	
	
	for (Int_t i = 0+1; i < rapid_bins-1; i++) // excluding first and last bin
	  if ( psiPrimeRapid > rapid_min + rapid_binSize*i &&
	       psiPrimeRapid < rapid_min + rapid_binSize*(i+1) ) {
	    if (year2011) {
	      //mass_low = mean_2011[i-1] - signal_win*sigma; mass_high = mean_2011[i-1] + signal_win*sigma;
	      fit_mass = mean_2011[i-1] ;  fit_sigma = sigma_2011[i-1] ; }
	    else if (year2012) {
	      //mass_low = mean[i-1] - signal_win*sigma; mass_high = mean[i-1] + signal_win*sigma;
	      fit_mass = mean[i-1] ;  fit_sigma = sigma[i-1] ;
	    }
  	  }
	
	if (1 // muon ID
	    //&& (*goodLooseMuon)[mupidx] && (*goodLooseMuon)[mumidx]
	    //&& ( ((*muQual)[mupidx]) & (1 << muonQual[3]) )  &&  ( ((*muQual)[mumidx]) & (1 << muonQual[3]) )
	    // new soft muons
	    && isNewSoftMuons
	    ) {
	  
	  if (HLT_5_vAll || HLT_2011 || HLT_3p5_vAll) { // trigger required here!
	    myMuMuMass_HLT->Fill( Jmass ) ;
	    if (HLT_3p5_vAll) myMuMuMass_HLT3p5->Fill( Jmass ) ;
	    /*
	      if ( fabs(muP.Rapidity()) < 2.4  &&  fabs(muM.Rapidity()) < 2.4) {
	      myMuMuMass_HLT_yMu_ypsi->Fill( Jmass ) ;
	      if (fabs(psiPrimeRapid) < 2.4) myMuMuMass_HLT_yMu->Fill( Jmass ) ; // no tighter selection
	      }
	    */
	    if ( HLT_5_vAll  /*&&  runNum >= 194735*/ ) myMuMuMass_HLT5_7comp->Fill( Jmass ) ;
	    
	    // psi(2S) variables
	    psi2S_y_h->Fill( psiPrimeRapid ) ;
	    psi2S_pT_h->Fill( psiPrime4Vec.Pt() ) ;
	    for (Int_t i = 0; i < abs(psi2S_2011_pT.size()); ++i)
	      if (HLT_2011_Dimuon[i]) {
		psi2S_2011_pT_h[i]->Fill( psiPrime4Vec.Pt() ) ;
	        for (Int_t j = 0; j < abs(HLT_v[i].size()); ++j)
		  if (HLT_2011_Dimuon_v[i][j])
		    psi2S_2011v_pT_h[i][j]->Fill( psiPrime4Vec.Pt() ) ;
	      }
	    
	    if ( 1 
		 &&  ( ( HLT_5_vAll && psiPrime4Vec.Pt() > 5 )
		       || ( HLT_2011 && psiPrime4Vec.Pt() > 7 ) 
		       || ( HLT_3p5_vAll && psiPrime4Vec.Pt() > 3.5 ) )
		 //&&  fabs(psiPrimeRapid) < 2.4 // useless given the following:
		 &&  fabs(muP.Rapidity()) < 2.4  &&  fabs(muM.Rapidity()) < 2.4
		 ) {
	      // multiplicity stuff
	      if (Jmass > myMuMuMass_HLT_final->GetXaxis()->GetXmin()  &&  Jmass < myMuMuMass_HLT_final->GetXaxis()->GetXmax()) {
		n_psiPrime++ ;
		//myMuMuMass_HLT_final->Fill( Jmass ) ;
	      
		if (n_psiPrime == 1) {
		  first_psi2S_index = myJIdx ;
		  first_psi2S_y = psiPrimeRapid ;
		  first_psi2S_mass = Jmass ;
		  first_psi2S_VtxCL = (*MuMuVtx_CL)[myJIdx] ;
		  psi2S_VtxCL_1_h->Fill( first_psi2S_VtxCL ) ;
		} else if (n_psiPrime == 2) {
		  second_psi2S_mass = Jmass ;
		  second_psi2S_VtxCL = (*MuMuVtx_CL)[myJIdx] ;
		  psi2S_VtxCL_2_h->Fill( second_psi2S_VtxCL ) ;
		}
	      }
	      
	      myMuMuMass_HLT_final->Fill( Jmass ) ;
	      mu1_DeltaR_mu2_h->Fill( muP.DeltaR(muM) ) ;
	      //mu2_DeltaR_mu1_h->Fill( muM.DeltaR(muP) ) ; // same as histogram above
	      if (muP.DeltaR(muM) < 1) myMuMuMass_HLT_DR1->Fill( Jmass ) ;
	      
	      // Hermine plot with fit mass and fit sigma
	      if ( fabs(Jmass - fit_mass) < 3*fit_sigma ) { 
		pT_vs_eta_mu_h->Fill(fabs(muP.PseudoRapidity()), muP.Pt()) ;
		pT_vs_eta_mu_h->Fill(fabs(muM.PseudoRapidity()), muM.Pt()) ;
		pT10_vs_eta_mu_h->Fill(fabs(muP.PseudoRapidity()), muP.Pt()) ;
		pT10_vs_eta_mu_h->Fill(fabs(muM.PseudoRapidity()), muM.Pt()) ;
		pT_vs_y_psi2S_h->Fill(psiPrimeRapid, psiPrime4Vec.Pt()) ;
	      }
	      
	      if (n_psiPrime < 2) { // check
		for (Int_t i = 0; i < rapid_bins; i++) {
		  if ( psiPrimeRapid > rapid_min + rapid_binSize*i &&
		       psiPrimeRapid < rapid_min + rapid_binSize*(i+1) ) { 
		    if (HLT_5_vAll) psiPrime_rapidBin[i]->Fill(Jmass); 
		    if (HLT_2011) psiPrime_2011_rapidBin[i]->Fill(Jmass); 
		  }
		}
	      }
	      
	    } // final cuts
	  } // HLT_*
	  if (HLT_7_vAll) myMuMuMass_HLT7->Fill( Jmass );
	} // muon ID
		
	// mu variables
	mu1_pixelHits_h->Fill( (*muPhits)[mupidx] ) ;  mu2_pixelHits_h->Fill( (*muPhits)[mumidx] ) ; 
	mu1_stripHits_h->Fill( (*muShits)[mupidx] ) ;  mu2_stripHits_h->Fill( (*muShits)[mumidx] ) ;
	mu1_charge_h->Fill( (*muCharge)[mupidx] ) ;    mu2_charge_h->Fill( (*muCharge)[mumidx] ) ;
	mu1_Chi2N_h->Fill( (*mu1Chi_MuMu2)[mupidx] / (*mu1NDF_MuMu)[mupidx] ) ;  mu2_Chi2N_h->Fill( (*mu2Chi2_MuMu)[mumidx] / (*mu2NDF_MuMu)[mumidx] ) ;
	mu1Tr_Chi2N_h->Fill( (*muChi2)[mupidx] / (*muNDF)[mupidx] ) ; mu2Tr_Chi2N_h->Fill( (*muChi2)[mumidx] / (*muNDF)[mumidx] ) ;
	mu1GlTr_Chi2N_h->Fill( (*muGlChi2)[mupidx] / (*muGlNDF)[mupidx] ) ; mu2GlTr_Chi2N_h->Fill( (*muGlChi2)[mumidx] / (*muGlNDF)[mumidx] ) ;
	mu1_Dz_h->Fill( (*muDzVtx)[mupidx] ) ;         mu2_Dz_h->Fill( (*muDzVtx)[mumidx] ) ;
	mu1_Dxy_h->Fill( (*muDxyVtx)[mupidx] ) ;       mu2_Dxy_h->Fill( (*muDxyVtx)[mumidx] ) ;
	
      } // PsiPrimeVtx_CL cut    
    } // end loop on psi'
  
  nPsi2S_h->Fill(n_psiPrime) ;
  
  if (n_psiPrime > 0  &&  n_psiPrime <= 2) {
    myMuMuMass_HLT_max2psiPrime->Fill( first_psi2S_mass ) ;
    if (n_psiPrime == 2) {
      myMuMuMass_HLT_max2psiPrime->Fill( second_psi2S_mass ) ;
      myMuMuMass_HLT_2psiPrime->Fill( first_psi2S_mass ) ;
      myMuMuMass_HLT_2psiPrime->Fill( second_psi2S_mass ) ;
      psi2S_VtxCL_diff_h->Fill( fabs(first_psi2S_VtxCL - second_psi2S_VtxCL) ) ;

      if ( fabs(first_psi2S_VtxCL - second_psi2S_VtxCL) < 0.3) {
	myMuMuMass_HLT_VtxCLDiff_cut->Fill( first_psi2S_mass ) ;
	myMuMuMass_HLT_VtxCLDiff_cut->Fill( second_psi2S_mass ) ;
      } else if ( fabs(first_psi2S_VtxCL - second_psi2S_VtxCL) > 0.3) {
	if (first_psi2S_VtxCL > 0.3)
	  myMuMuMass_HLT_VtxCLDiff_cut->Fill( first_psi2S_mass ) ;
	if (second_psi2S_VtxCL > 0.3)
	  myMuMuMass_HLT_VtxCLDiff_cut->Fill( second_psi2S_mass ) ;
      }
      
    }
  }

  // only 1 psi' since now:
  if (only1psiPrime  &&  n_psiPrime != 1) return kTRUE ;

  if (only1psiPrime) {
    for (Int_t i = 0+1; i < rapid_bins-1; i++) { // excluding first and last bin
      if ( first_psi2S_y > rapid_min + rapid_binSize*i  &&
	   first_psi2S_y < rapid_min + rapid_binSize*(i+1) ) {
	if (!fixedMassWind) 
	  if (!sideband)
	    if ( fabs(first_psi2S_mass - fit_mass) > signal_win*fit_sigma ) return kTRUE ;
	    else myMuMuMass_HLT_1psiPrime->Fill( first_psi2S_mass ) ;
	  else // sideband:
	    if ( fabs(first_psi2S_mass - fit_mass) < sb_low*fit_sigma  ||  fabs(first_psi2S_mass - fit_mass) > sb_high*fit_sigma ) return kTRUE ;
	    else myMuMuMass_HLT_1psiPrime->Fill( first_psi2S_mass ) ;
	else if (fixedMassWind) {
	  //Float_t fixed_sigma = 0.150 ;
	  Float_t fixed_sigma = 0.200 ;
	  if (!sideband)
	    if ( fabs(first_psi2S_mass - psi2S_mass) > fixed_sigma ) return kTRUE ;
	    else myMuMuMass_HLT_1psiPrime->Fill( first_psi2S_mass ) ;
	  else // sideband:
	    if ( fabs(first_psi2S_mass - fit_mass) < sb_low*fixed_sigma  ||  fabs(first_psi2S_mass - fit_mass) > sb_high*fixed_sigma ) return kTRUE ; // check sidebands
	    else myMuMuMass_HLT_1psiPrime->Fill( first_psi2S_mass ) ;
	}
      }
    }
  }

  
  //cout<<"nZ="<<nZ<<", size="<< (*B_cMass).size()<<endl;  
        
  /*
    for (int myZIdx = 0; myZIdx < abs(nZ); myZIdx++) {
    if ( 1
    //&& fabs((*JPsiMass)[myXIdx]-3.0967)<0.1
    && (*MuMuPiVtx_CL)[myZIdx] > 0.1
    && fabs( (*MuMuMass)[ (*JPsiIndex)[myZIdx] ] - psi2S_mass ) < 0.2 
    && sqrt( pow((*B0_fkaonPx)[myXIdx], 2) + pow((*B0_fkaonPy)[myXIdx], 2) ) > 0.9
    && ( (*B0_CTauBS)[ myXIdx] /  (*B0_CTauBSE)[ myXIdx] ) > 3.0
    && sqrt( pow((*MuMuPx)[(*JPsiIndex)[myXIdx]], 2) + pow((*MuMuPy)[(*JPsiIndex)[myXIdx]], 2) ) > 7.0 
	 
	 
    &&  fabs((*ZCosAlphaPV)[myZIdx]) > cosAlpha_cut2 // = 0.99, our is cosAlpha_cut = 0.9, 1k difference
    &&  fabs(Kp4.Eta()) < 2.4
    &&  fabs(mu1.Eta()) < 2.0 && fabs(mu2.Eta()) < 2.0 // 3p5 is 2.1 
    &&  Kp4.Pt() > 1.0
    //
    &&  (*ZCTauBS)[myZIdx] / (*ZCTauBSE)[myZIdx] > 4.0 // our is 3.0: 5k difference
    &&  Jp4.Pt() > 7 //7 our is 5: 1k difference
	 
    )
    {
    myMuMuKMass->Fill( (*B0_Mass)[myXIdx] -  (*MuMuMass)[ (*JPsiIndex)[myXIdx] ] + psi2S_mass );         
    }
    }
  
    //continue;
    */
  
  //( ((*muQual)[ (*MyFourMuonMu1Idx)[myFourMuIdx] ]) & (1 << muonQual[0] ) )

  //nZ_h->Fill(nZ);
  
  nB0_h->Fill( nB0 );
  // if (nMCB0 >= 1) for matching
  for (Int_t myB0Idx = 0; myB0Idx < abs(nB0); myB0Idx++) {
    B0_Mass_h->Fill((*B0Mass)[myB0Idx]) ;

    int myjpsiindex = 0 ;
    if (!only1psiPrime)
      myjpsiindex = (*B0MuMuIdx)[myB0Idx] ;
    else
      myjpsiindex = first_psi2S_index ;
    
    TLorentzVector Jp4, Kp4 ; //Pip4,
    //      
    Float_t J_E = 0. ;
    J_E = sqrt( pow((*MuMuPx)[myjpsiindex], 2) + pow((*MuMuPy)[myjpsiindex], 2) + pow((*MuMuPz)[myjpsiindex], 2) + pow((*MuMuMass)[myjpsiindex], 2) ) ;    
    //
    Jp4.SetPxPyPzE( (*MuMuPx)[myjpsiindex], (*MuMuPy)[myjpsiindex], (*MuMuPz)[myjpsiindex], J_E) ;
    //Pip4.SetPxPyPzE( (*PionPx_MuMuPiK)[myB0Idx], (*PionPy_MuMuPiK)[myB0Idx], (*PionPz_MuMuPiK)[myB0Idx], (*PionE_MuMuPiK)[myB0Idx]) ;
    Kp4.SetPxPyPzE( (*KaonPx_MuMuPiK)[myB0Idx], (*KaonPy_MuMuPiK)[myB0Idx], (*KaonPz_MuMuPiK)[myB0Idx], (*KaonE_MuMuPiK)[myB0Idx]) ;

    for (Int_t i = 0; i < abs(trk_pT.size()); ++i) 
      if ( Kp4.Pt() > trk_pT[i] )
	MuMuKMass_trk2Pt[i]->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;

  }

  // dEdx variables
  const float KData = 2.65945;
  const float CData = 2.79998;
  //
  //const float KMC = 2.34138;
  //const float CMC = 2.62959;
        
  //Mass bands cuts
  const float dEdxCut = 4.1;
  const float ProtonMin = 0.78;
  const float ProtonMax = 1.3;
  const float KaonMin = 0.43;
  const float KaonMax = 0.70;
  const float DeuteronMin = 1.6;
  const float DeuteronMax = 2.15;

  Int_t nB0AC_noMassWin = 0, nB0AC = 0, nB0AC_purityWin = 0, nB0AC_signalWin = 0, nB0ACInMC = 0, nTwins = 0, nTwins_signalWin = 0 ;
  //
  Float_t priVtx_XYZ[5], priVtx_B0CosAlpha_XYZ[2][5], priVtxB0Less_XYZ[5], priVtxB0Less_B0CosAlpha_XYZ[2][5], priVtxGen_XYZ[3], genEvtVtx_XYZ[3] ;
  Float_t pB0_XYZ[3], secVtx_XYZ[3] ;
  for (Int_t i=0; i<5; i++) {
    priVtx_XYZ[i] = 999. ;
    priVtx_B0CosAlpha_XYZ[0][i] = 999. ; priVtx_B0CosAlpha_XYZ[1][i] = 999. ;
    priVtxB0Less_XYZ[i] = 999. ;
    priVtxB0Less_B0CosAlpha_XYZ[0][i] = 999. ; priVtxB0Less_B0CosAlpha_XYZ[1][i] = 999. ;
    if (i<3) {
      priVtxGen_XYZ[i] = 999. ; genEvtVtx_XYZ[i] = 999. ; pB0_XYZ[i] = 999. ; secVtx_XYZ[i] = 999. ;
    }
  }
  //
  Float_t momentumCut = 0.1;
  Float_t DeltaR_1B0[2] = {99,99}, DeltaPt_1B0[2] = {99,99} ;
  Int_t chargePi_1B0[2] = {0,0}, chargeK_1B0[2] = {0,0} ; 
  Bool_t MuMuMatched_1B0[2] = {kFALSE,kFALSE}, pionMatched_1B0[2] = {kFALSE,kFALSE}, kaonMatched_1B0[2] = {kFALSE,kFALSE} ;

  // twin variables
  Float_t DeltaR_2B0[2][2] = {{99,99},{99,99}}, DeltaPt_2B0[2][2] = {{99,99},{99,99}}, B0VtxCL_2B0[2] = {-1,-1}, B0Mass_2B0[2][2] = { {-1,-1},{-1,-1} }, B0CTau_2B0[2][2] = { {-1,-1},{-1,-1} }, alpha_2B0[2] = {-2,-2} ;
  Double_t DrPi_2B0[2] = {99,99}, DrK_2B0[2] = {99,99} ;
  Int_t chargePi_2B0[2][2] = { {0,0},{0,0} }, chargeK_2B0[2][2] = { {0,0},{0,0} } ;
  Float_t B0Mass_gen = 0 ;
  TLorentzVector PsiPp4_2B0[2][2], Pip4_2B0[2][2], Kp4_2B0[2][2] ;
  Double_t theta_psi = 10., cosTheta_psi = 10., theta_psi_rot = 10., cosTheta_psi_rot = 10., phi = 10. ; 
  Int_t PiCharge_2B0[2][2] = { {0,0},{0,0} }, KCharge_2B0[2][2] = { {0,0},{0,0} } ; 
  //Bool_t MuMuMatched[2][2] = {{kFALSE,kFALSE},{kFALSE,kFALSE}}, pionMatched[2][2] = {{kFALSE,kFALSE},{kFALSE,kFALSE}}, kaonMatched[2][2] = {{kFALSE,kFALSE},{kFALSE,kFALSE}} ;
  Bool_t MuMuMatched[2][2] = { {kFALSE,kFALSE},{kFALSE,kFALSE} }, pionMatched[2][2] = { {kFALSE,kFALSE},{kFALSE,kFALSE} }, kaonMatched[2][2] = { {kFALSE,kFALSE},{kFALSE,kFALSE} } ;
  Float_t twins_var[2][4] = { {0,0,0,0},{0,0,0,0} } ;


  //for (Int_t myZIdx = 0; myZIdx < abs(nZ); myZIdx++) {
  for (Int_t myB0Idx = 0; myB0Idx < abs(nB0); myB0Idx++) {  
    int myjpsiindex = 0 ;
    if (!only1psiPrime)
      myjpsiindex = (*B0MuMuIdx)[myB0Idx] ;
    else
      myjpsiindex = first_psi2S_index ;
    //if (myB0Idx == 0) myMuMuMass_test->Fill( (*MuMuMass)[myjpsiindex] );
    int mymupidx = (*mu1Idx)[myjpsiindex] ;
    int mymumidx = (*mu2Idx)[myjpsiindex] ;

    mu_px->Fill( (*muPx)[mymupidx] ) ;  //mu_px->Fill( (*muPx)[mymumidx] ) ;
    mu_px_psi2S->Fill( (*mu1Px_MuMu)[mymupidx] ) ; //mu_px_psi2S->Fill( (*mumPx)[mymumidx] ) ;
    fitmu_px->Fill( (*Muon1Px_MuMuPiK)[myjpsiindex] ) ; //fitmu_px->Fill( (*fitmu2Px)[myjpsiindex] ) ;
    TLorentzVector mu1, mu2, Jp4, Bcp4, Zp4, Pip4, Kp4 ;
    //      
    Float_t J_E = 0., mu1_E = 0., mu2_E = 0., K_E = 0. ; // BcE = 0.,
    // 
    mu1_E = sqrt( pow((*muPx)[mymupidx], 2) + pow((*muPy)[mymupidx], 2) + pow((*muPz)[mymupidx], 2) + pow(muon_mass, 2) ) ;    
    mu2_E = sqrt( pow((*muPx)[mymumidx], 2) + pow((*muPy)[mymumidx], 2) + pow((*muPz)[mymumidx], 2) + pow(muon_mass, 2) ) ;    
    J_E = sqrt( pow((*MuMuPx)[myjpsiindex], 2) + pow((*MuMuPy)[myjpsiindex], 2) + pow((*MuMuPz)[myjpsiindex], 2) + pow((*MuMuMass)[myjpsiindex], 2) ) ;
    //
    mu1.SetPxPyPzE( (*muPx)[mymupidx], (*muPy)[mymupidx], (*muPz)[mymupidx], mu1_E) ;
    mu2.SetPxPyPzE( (*muPx)[mymumidx], (*muPy)[mymumidx], (*muPz)[mymumidx], mu2_E) ;
    Jp4.SetPxPyPzE( (*MuMuPx)[myjpsiindex], (*MuMuPy)[myjpsiindex], (*MuMuPz)[myjpsiindex], J_E) ;

    Int_t mypiIdx = (*B0PionIdx)[myB0Idx] ;
    //cout<<"mypiIdx="<<(*muPhits)[mymumidx]<<endl;

    //Pip4.SetPxPyPzE( (*PionPx_MuMuPi)[myB0Idx], (*PionPy_MuMuPi)[myB0Idx], (*PionPz_MuMuPi)[myB0Idx], (*PionE_MuMuPi)[myB0Idx]) ; // refitted
    Pip4.SetPxPyPzE( (*trackPx)[mypiIdx], (*trackPy)[mypiIdx], (*trackPz)[mypiIdx], (*trackEnergy)[mypiIdx]) ; // original

    //K_E = sqrt( pow((*PionPx_MuMuPi)[myB0Idx], 2) + pow((*PionPy_MuMuPi)[myB0Idx], 2) + pow((*PionPz_MuMuPi)[myB0Idx], 2) + pow(kaonCh_mass, 2) ) ; // refitted
    //Kp4.SetPxPyPzE( (*PionPx_MuMuPi)[myB0Idx], (*PionPy_MuMuPi)[myB0Idx], (*PionPz_MuMuPi)[myB0Idx], K_E) ; // refitted
    K_E = sqrt( pow((*trackPx)[mypiIdx], 2) + pow((*trackPy)[mypiIdx], 2) + pow((*trackPz)[mypiIdx], 2) + pow(kaonCh_mass, 2) ) ; // original   
    Kp4.SetPxPyPzE( (*trackPx)[mypiIdx], (*trackPy)[mypiIdx], (*trackPz)[mypiIdx], K_E) ; // original

    //          BcE = sqrt( pow((*xPx)[myB0Idx],2) + pow((*xPy)[myB0Idx],2) + pow((*xPz)[myB0Idx],2) + pow((*xMass)[myB0Idx],2) ) ;      
    //Bcp4.SetPxPyPzE((*xPx)[myB0Idx],(*xPy)[myB0Idx],(*xPz)[myB0Idx],BcE);
    Bcp4 = Jp4 + Kp4 ;
    Zp4 = Jp4 + Pip4 ;

    //Float_t Z_Pt = sqrt( pow((*MuMuPiPx)[myB0Idx],2) + pow((*MuMuPiPy)[myB0Idx],2) ) ;
    Float_t Z_Pt = Zp4.Pt() ;

    // dEdx computation
    Int_t kaon_pID = -1;
    if ( (*trackShits)[mypiIdx] < 4  ||
	 (*trackChi2)[mypiIdx] / (*trackNDF)[mypiIdx] > (1.5*((*trackShits)[mypiIdx]) + 4*sqrt(1.5*((*trackShits)[mypiIdx])))  ||
	 (*tr_dedx)[mypiIdx] < 1e-6 ) {
      kaon_pID = -1;
    }
    //cout << "no found"<<endl;}
    else {
      float dEdxMass = ComputedEdxMass(mypiIdx, Pip4.Pt(), KData, CData);
      //cout << dEdxMass << " Mass " << (*tr_dedxMass)[mypiIdx]; << endl;

      dEdxMass_h[0]->Fill(dEdxMass) ;
      dEdxMass_vs_pT->Fill(Pip4.Pt(), dEdxMass) ;
      dEdx_vs_pT[0]->Fill(Pip4.Pt(), ((*tr_dedx)[mypiIdx])) ;
      dEdx_vs_p->Fill(Pip4.P(), ((*tr_dedx)[mypiIdx])) ;

      float dEdxMin = dEdxCut;        
      //float theo_dEdx = ComputeTheodEdx(mypiIdx, Pip4.Pt(), KData, CData);
      //cout << theo_dEdx << " Mass " << (*tr_theo)[mypiIdx];;
      // if( (*tr_dedx)[mypiIdx] > dEdxMin)
      // {
      kaon_pID = 0;
      if (dEdxMass > ProtonMin  &&  dEdxMass < ProtonMax)
	{
	  kaon_pID = 2212;
	}
      else if(dEdxMass > KaonMin && dEdxMass < KaonMax)
	{
	  kaon_pID = 310;
	}
      else if(dEdxMass > DeuteronMin && dEdxMass < DeuteronMax)
	{
	  kaon_pID = 1000;
	}
      //}
      else if( (*tr_dedx)[mypiIdx] < (dEdxMin - 1.)  &&  (*tr_dedx)[mypiIdx] > (dEdxMin - 2.)  &&  Pip4.Pt() < 1.1  &&  dEdxMass < 0.2  &&  dEdxMass > 0.05)
	{
	  kaon_pID = 135;
	} 
    }

    // to use with a 2nd track
    /*
    int kaon_mID = -1;
    if ( pimShits < 4 || pimchi2NDF > (1.5*pimShits + 4*sqrt(1.5*pimShits)) || pim_dedx < 1e-6) {
      kaon_mID = -1;
      //cout << "no found"<<endl;}
    }
    else {
      float dEdxMass = ComputedEdxMass(mykIdx, Kp4.Pt(), KData, CData);
      //cout << dEdxMass << " Mass " << pip_dedxMass << endl;
      float dEdxMin = dEdxCut;
      
      //float theo_dEdx = ComputeTheodEdx(pipNr, pipPt, KData, CData);
      //cout << theo_dEdx << " Mass " << pip_theo;
      //  if(pim_dedx > dEdxMin)
      // {
      kaon_mID = 0;
      if(dEdxMass > ProtonMin && dEdxMass < ProtonMax)
      {
	kaon_mID = 2212;
      }
      else if(dEdxMass > KaonMin && dEdxMass < KaonMax)
	{
	  kaon_mID = 310;
	}
      else if(dEdxMass > DeuteronMin && dEdxMass < DeuteronMax)
	{
	  kaon_mID = 1000;
	}
      // }
      else if(pim_dedx < (dEdxMin - 1.) && pim_dedx > (dEdxMin - 2.) && pimPt < 1.1 && dEdxMass < 0.2 && dEdxMass > 0.05)
	{
	  kaon_mID = 135;
	}
    }
	
    int KaonPairID;
    if (kaon_pID == 310  &&  kaon_mID == 310) {
      KaonPairID = 1 ;
      //Phi_identified->Fill(KKMass2);
    } else if (kaon_pID == 310  ||  kaon_mID == 310) {
      //Phi_Oneidentified->Fill(KKMass2);
      KaonPairID =2;
    } else {
      KaonPairID =0;
      //Phi_Notidentified->Fill(KKMass2);
    }
    */
    
    //myB0_Mass_h->Fill( (Jp4 + Kp4 + Pip4).M() ) ; // meaningless because Kp4 = Pip4

    Z_pT_h->Fill( Z_Pt ) ;
    //Z_VtxCL_h->Fill( (*MuMuPiVtx_CL)[myB0Idx] ) ;
    //Z_CTauBS_h->Fill( (*ZCTauBS)[myB0Idx] ) ;
    //Z_CTauBSE_h->Fill( (*ZCTauBSE)[myB0Idx] ) ;
    //Z_CTauBS_over_Z_CTauBSE_h->Fill( (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] ) ; moved after cuts
    //Z_pointingAngle[0]->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;
    //Float_t coeff = (*ZCosAlphaPV)[myB0Idx] < 0 ? -cosAlpha_cut/scale_factor + ((*ZCosAlphaPV)[myB0Idx] +cosAlpha_cut)*( (Z_cosPA_max-cosAlpha_cut/scale_factor)/(Z_cosPA_max-cosAlpha_cut) ) : +cosAlpha_cut/scale_factor + ((*ZCosAlphaPV)[myB0Idx] -cosAlpha_cut)*( (Z_cosPA_max-cosAlpha_cut/scale_factor)/(Z_cosPA_max-cosAlpha_cut) ) ;
    // mu variables
    mu1_y_h->Fill( mu1.Rapidity() ) ; mu2_y_h->Fill( mu2.Rapidity() ) ;
    mu1_pT_h->Fill( mu1.Pt() ) ; mu2_pT_h->Fill( mu2.Pt() ) ;
    //mu1_vs_mu2_pT_h->Fill( mu1.Pt(), mu2.Pt() ) ;
    pT_vs_y_mu1_h->Fill( mu1.Rapidity(), mu1.Pt()) ; 
    pT_vs_y_mu2_h->Fill( mu2.Rapidity(), mu2.Pt()) ; 
    //psi2S_pointingAngle->Fill( (*JCosAlpha)[myjpsiindex] ) ; 
    // pi+ variables
    pi_pixelHits_h->Fill( (*trackPhits)[mypiIdx] ) ;
    pi_stripHits_h->Fill( (*trackShits)[mypiIdx] ) ;
    pi_Chi2N_h->Fill( (*trackChi2)[mypiIdx] / (*trackNDF)[mypiIdx] ) ;
    pi_pT_h->Fill( Pip4.Pt() ) ;
    pi_y_h->Fill( Pip4.Rapidity() ) ;
      
    // Kai cuts
    if ( 1
	 //&&  (*MuMuPiVtx_CL)[myB0Idx] > 0.1 // 0.1 Kai
	 &&  fabs(Pip4.Eta()) < 1.5
	 &&  fabs(Jp4.Eta()) < 1.5
	 &&  fabs( Jp4.M() - psi2S_mass ) < 0.06  // we applied a moving window on the onePsiPrime selection
	 //&&  (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] < 3.0
	 //&&  (*ZCTauBSE)[myB0Idx] < 0.01 // only 150 events cut with the other cuts already applied
	 &&  Pip4.Pt() > 7.0
	 //&&  muon ID
	 &&  mu1.Pt() > 7.0  &&  mu2.Pt() > 7.0
	 ) {
      myMuMuPiMass_Kai->Fill( Zp4.M() - Jp4.M() + psi2S_mass ) ;
      if (HLT_7_vAll  &&  Jp4.Pt() > 7.) myMuMuPiMass_Kai_HLT7->Fill( Zp4.M() - Jp4.M() + psi2S_mass ) ;
    }


    if ( 1  // reproducing Kai cuts for B+ -> psi(2S) phi K+
	 &&  HLT_5_vAll
	 ) {
      myMuMuMass_noMuID_KaiBF->Fill( Jp4.M() ) ;
      if ( 1
	   //&&  (*goodLooseMuon)[mymupidx] && (*goodLooseMuon)[mymumidx]
	   &&  ( ((*muQual)[mymupidx]) & (1 << muonQual[3]) )  &&  ( ((*muQual)[mymumidx]) & (1 << muonQual[3]) )
	   &&  ( ( (*muChi2)[mymupidx] / (*muNDF)[mymupidx] ) <  1.8 )  &&  ( ( ((*muChi2)[mymumidx]) / (*muNDF)[mymumidx] ) < 1.8 )
	   &&  fabs((*muDzVtx)[mymupidx]) < 30.0  &&  fabs((*muDzVtx)[mymumidx]) < 30.0 
	   &&  fabs((*muDxyVtx)[mymupidx]) < 3.0  &&  fabs((*muDxyVtx)[mymumidx]) < 3.0 
	   &&  (*muPhits)[mymupidx] >= 1 && (*muPhits)[mymumidx] >= 1 
	   &&  (*muShits)[mymupidx] >= 5 && (*muShits)[mymumidx] >= 5 
	   ) {
	myMuMuMass_MuID_KaiBF->Fill( Jp4.M() ) ;
	if ( 1
	     //&&  fabs((*ZCosAlphaPV)[myB0Idx]) > cosAlpha_cut2 // = 0.99, our is cosAlpha_cut = 0.9, 1k difference
	     &&  fabs(Kp4.Eta()) < 2.4
	     &&  fabs(mu1.Eta()) < 2.0 && fabs(mu2.Eta()) < 2.0 // 3p5 is 2.1 
	     &&  Kp4.Pt() > 1.0
	     //
	     //&&  (*MuMuPiVtx_CL)[myB0Idx] > 0.1 // our is 0.01: 10k difference
	     //&&  (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] > 4.0 // our is 3.0: 5k difference
	     &&  Jp4.Pt() > 7 //7 our is 5: 1k difference
	     ) { 
	  if (!only1psiPrime) {
	    if (fixedMassWind) {
	      //Float_t fixed_sigma = 0.150 ;
	      Float_t fixed_sigma = 0.200 ;
	      if ( fabs(Jp4.M() - psi2S_mass) < fixed_sigma ) {
		//myMuMuMass_KaiBF->Fill( Jp4.M() ) ;
		myMuMuKMass_KaiBF->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;

		if ( ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass) > 5.268  &&  ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass) < 5.288 )
		  myMuMuMass_KaiBF->Fill( Jp4.M() ) ;
	      }
	    } else if ( fabs( Jp4.M() - psi2S_mass ) < 0.12 ) {
	      myMuMuKMass_KaiBF->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	    }
	  } else if (only1psiPrime) {
	    myMuMuKMass_KaiBF->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	      
	    if ( ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass) > 5.268  &&  ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass) < 5.288 ) 
	      myMuMuMass_KaiBF->Fill( Jp4.M() ) ;
	  }
	} // remaining cuts
      } // muID
    } // HLT_5
      
      //////////////////////////////////////// Sara plots
    if ( HLT_5_vAll ) { // trigger
      if ( 1 // muonID
	   &&  ( ((*muQual)[mymupidx]) & (1 << muonQual[3]) )  &&  ( ((*muQual)[mymumidx]) & (1 << muonQual[3]) )
	   &&  ( ( (*muChi2)[mymupidx] / (*muNDF)[mymupidx] ) <  1.8 )  &&  ( ( ((*muChi2)[mymumidx]) / (*muNDF)[mymumidx] ) < 1.8 )
	   &&  fabs((*muDzVtx)[mymupidx]) < 30.0  &&  fabs((*muDzVtx)[mymumidx]) < 30.0 
	   &&  fabs((*muDxyVtx)[mymupidx]) < 3.0  &&  fabs((*muDxyVtx)[mymumidx]) < 3.0 
	   &&  (*muPhits)[mymupidx] >= 1 && (*muPhits)[mymumidx] >= 1 
	   &&  (*muShits)[mymupidx] >= 5 && (*muShits)[mymumidx] >= 5 
	   ) {
	/*
	  for (Int_t i = 0; i < abs(Z_CTauBSOverE.size()); ++i) {
	  if ( (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] > Z_CTauBSOverE[i])
	  MuMuKMass_ZCTauBSOverE[i]->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	  }*/
	for (Int_t i = 0; i < abs(mu_pT.size()); ++i) {
	  if ( mu1.Pt() > mu_pT[i] )
	    MuMuKMass_mu1Pt[i]->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	  if ( mu2.Pt() > mu_pT[i] )
	    MuMuKMass_mu2Pt[i]->Fill(( Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	}
	for (Int_t i = 0; i < abs(trk_pT.size()); ++i) {
	  if ( Pip4.Pt() > trk_pT[i] )
	    MuMuKMass_trk1Pt[i]->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	  //if ( > trk_pT[i] ) // moved in the B0 loop above
	  //MuMuKMass_trk2Pt[i]->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	}
	/*
	  for (Int_t i = 0; i < abs(Z_CL.size()); ++i)
	  if ( (*MuMuPiVtx_CL)[myB0Idx] > Z_CL[i] )
	  MuMuKMass_ZCL[i]->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	*/
      } // muonID
    } // trigger
      
      //////////////////////////////////////// Nairit cuts
    if ( HLT_5_vAll ) { // trigger

      if ( 1 // muonID
	   &&  ( ((*muQual)[mymupidx]) & (1 << muonQual[3]) )  &&  ( ((*muQual)[mymumidx]) & (1 << muonQual[3]) )
	   &&  ( ( (*muChi2)[mymupidx] / (*muNDF)[mymupidx] ) < 3 )  &&  ( ( ((*muChi2)[mymumidx]) / (*muNDF)[mymumidx] ) < 3 )
	   &&  fabs((*muDzVtx)[mymupidx]) < 20.0  &&  fabs((*muDzVtx)[mymumidx]) < 20.0 
	   &&  fabs((*muDxyVtx)[mymupidx]) < 0.3  &&  fabs((*muDxyVtx)[mymumidx]) < 0.3 
	   &&  (*muPhits)[mymupidx] > 0 && (*muPhits)[mymumidx] > 0 
	   &&  (*muShits)[mymupidx] > 5 && (*muShits)[mymumidx] > 5 
	   ) {

	Bool_t mu1_PtSel = false;
	Bool_t mu2_PtSel = false;
	if ((((fabs( mu1.Eta() ) < 1.2) && (mu1.Pt()>4.))) || ((mu1.Eta()>=1.2 || mu1.Eta() <=-1.2) && (mu1.Pt()>3.3))) mu1_PtSel = true;
	if ((((fabs( mu2.Eta() ) < 1.2) && (mu2.Pt()>4.))) || ((mu2.Eta()>=1.2 || mu2.Eta() <=-1.2) && (mu2.Pt()>3.3))) mu2_PtSel = true;

	if ( 1
	     && Jp4.Pt() > 5.05 // unchanged      
	     && fabs(mu1.Eta()) < 2.0 && fabs(mu2.Eta()) < 2.0 // unchanged     
	     && fabs( Jp4.M() - psi2S_mass ) < 0.12  // same as && fabs((*MuMuMass)[myPsiP_fromZ_Index]-psi2SPDGmass) < 0.12  unchanged
	     && ( (*MuMuVtx_CL)[myjpsiindex] > 0.02 ) // optimised value
	     && mu1.Pt() > 3.3  &&  mu2.Pt() > 3.3 // this cut is ok with muons at any eta; later we ask >4.0 for central muons
	     //
	     && mu1_PtSel && mu2_PtSel
	     ) {
	  if ( (((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass > 5.15)  &&  ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass < 5.2))  || 
	       (((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass > 5.35)  &&  ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass < 5.4) )) {
	    MuMuKMass_sb->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	    dEdx_vs_pT[2]->Fill(Kp4.Pt(), ((*tr_dedx)[mypiIdx])) ;	 
	  }

	  if ( 1 // B+ cuts
		 //&& fabs((*ZCosAlphaPV)[myB0Idx]) > 0.998 // optimised value
		 //&& (*MuMuPiVtx_CL)[myB0Idx] > 0.425 // optimised value
		 //&& (*MuMuPiVtx_CL)[myB0Idx] > 0.1
		 //&& ((*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx]) > 3.5 //8.5 // optimised value (just over 90 % purity)
	       && Kp4.Pt() > 0.4 //0.7
	       && fabs(Kp4.Eta()) < 2.4 // unchanged
	       //                 
	       && ((*trackChi2)[mypiIdx]/(*trackNDF)[mypiIdx]) < 5.0
	       && (*trackShits)[mypiIdx] >= 7
	       && (*trackPhits)[mypiIdx] > 0
	       && fabs(Kp4.Eta()) < 2.4
	       && Jp4.DeltaR(Kp4) < 0.4
	       && Bcp4.Pt() > 5 //7
	       ) {
	    MuMuKMass_Nairit->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
	    dEdx_vs_p_Bplus->Fill(Kp4.P(), ((*tr_dedx)[mypiIdx])) ;

	    if ( ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass > 5.25)  &&  ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass < 5.31) ) {
	      trk_pT_test_h->Fill( Kp4.Pt() ) ;
	      dEdx_vs_pT[1]->Fill(Kp4.Pt(), ((*tr_dedx)[mypiIdx])) ;	 
	      dEdxMass_h[1]->Fill( ComputedEdxMass(mypiIdx, Kp4.Pt(), KData, CData) ) ;
	    }
	      
	  }     
	}
	  
      } // muonID
    } // trigger


    //////////////////////////////////////// EXCLUSIVE KPi    
    if ( HLT_5_vAll ) 
      {
	//for (Int_t myB0Idx=0; myB0Idx<abs(nB0); myB0Idx++) {  // already opened above    
	// set indices
	Int_t psip_index = (*B0MuMuIdx)[myB0Idx];
	Int_t mu1_fromB0_index = (*mu1Idx)[psip_index];
	Int_t mu2_fromB0_index = (*mu2Idx)[psip_index];
	Int_t pi_orig_Index = (*B0PionIdx)[myB0Idx];
	Int_t ka_orig_Index = (*B0KaonIdx)[myB0Idx];
	//
	// Muon 1 & 2
	// this is the 4-momentum of the original muons
	TLorentzVector mu1_p4, mu2_p4;
	Float_t mu1E = 0., mu2E = 0.;
	//
	mu1E = sqrt( pow((*muPx)[mu1_fromB0_index],2) + pow((*muPy)[mu1_fromB0_index],2) + pow((*muPz)[mu1_fromB0_index],2) + pow(muon_mass,2) );
	mu2E = sqrt( pow((*muPx)[mu2_fromB0_index],2) + pow((*muPy)[mu2_fromB0_index],2) + pow((*muPz)[mu2_fromB0_index],2) + pow(muon_mass,2) );
	mu1_p4.SetPxPyPzE( (*muPx)[mu1_fromB0_index], (*muPy)[mu1_fromB0_index], (*muPz)[mu1_fromB0_index], mu1E);
	mu2_p4.SetPxPyPzE( (*muPx)[mu2_fromB0_index], (*muPy)[mu2_fromB0_index], (*muPz)[mu2_fromB0_index], mu2E);
	//
	TLorentzVector PsiPp4;
	Float_t PsiP_E = 0.;
	PsiP_E = sqrt( pow((*MuMuPx)[psip_index],2) + pow((*MuMuPy)[psip_index], 2) + pow((*MuMuPz)[psip_index], 2) + pow((*MuMuMass)[psip_index], 2) ) ;
	PsiPp4.SetPxPyPzE( (*MuMuPx)[psip_index], (*MuMuPy)[psip_index], (*MuMuPz)[psip_index], PsiP_E) ;
	hmyPsiP_Pt->Fill(PsiPp4.Pt());
	//
	TLorentzVector pip4, pip4_exchanged;
	pip4.SetPxPyPzE( (*PionPx_MuMuPiK)[myB0Idx], (*PionPy_MuMuPiK)[myB0Idx], (*PionPz_MuMuPiK)[myB0Idx], (*PionE_MuMuPiK)[myB0Idx]) ;
	pip4_exchanged.SetVectM(pip4.Vect(), kaonCh_mass) ;
	      
	TLorentzVector kp4, kp4_exchanged;
	kp4.SetPxPyPzE( (*KaonPx_MuMuPiK)[myB0Idx], (*KaonPy_MuMuPiK)[myB0Idx], (*KaonPz_MuMuPiK)[myB0Idx], (*KaonE_MuMuPiK)[myB0Idx]) ;
	kp4_exchanged.SetVectM(kp4.Vect(), pionCh_mass) ;
	      
	TLorentzVector Pip4_orig, Pip4_orig_ex, Kp4_orig, Kp4_orig_ex;
	Pip4_orig.SetPxPyPzE( (*trackPx)[pi_orig_Index], (*trackPy)[pi_orig_Index], (*trackPz)[pi_orig_Index], (*trackEnergy)[pi_orig_Index]) ; // original
	Pip4_orig_ex.SetVectM(Pip4_orig.Vect(), kaonCh_mass) ;

	K_E = sqrt( pow((*trackPx)[ka_orig_Index],2) + pow((*trackPy)[ka_orig_Index],2) + pow((*trackPz)[ka_orig_Index],2) + pow(kaonCh_mass,2) ) ; 
	Kp4_orig.SetPxPyPzE( (*trackPx)[ka_orig_Index], (*trackPy)[ka_orig_Index], (*trackPz)[ka_orig_Index], K_E) ; 
	//Kp4_orig.SetPxPyPzE( (*trackPx)[ka_orig_Index], (*trackPy)[ka_orig_Index], (*trackPz)[ka_orig_Index], (*trackEnergy)[ka_orig_Index]) ; // original but trackEnergy is filled with the mass assumpion of ?
	Kp4_orig_ex.SetVectM(Kp4_orig.Vect(), pionCh_mass) ;
	      
	TLorentzVector Kstarp4, Kstarp4_orig, Kstarp4_orig_ex, myB0p4, myB0p4_orig, myB0p4_KPiExchanged;
	Kstarp4 = pip4 + kp4 ;
	Kstarp4_orig = Pip4_orig + Kp4_orig ;
	Kstarp4_orig_ex = Pip4_orig_ex + Kp4_orig_ex ;
	myB0p4 = PsiPp4 + pip4 + kp4;
	myB0p4_orig = PsiPp4 + Pip4_orig + Kp4_orig;
	// alternatively:
	Float_t myB0MassAlt = (PsiPp4 + Pip4_orig + Kp4_orig).Mag() - PsiPp4.M() + psi2S_mass;
	Float_t myB0MassAlt_ex = (PsiPp4 + Pip4_orig_ex + Kp4_orig_ex).Mag() - PsiPp4.M() + psi2S_mass;
	//
	myB0p4_KPiExchanged = PsiPp4 + pip4_exchanged + kp4_exchanged;
	//
	// Armenteros plot
	// original tracks
	Float_t pi_KstarPl_orig = ( (Pip4_orig.Vect()).Dot(Kstarp4_orig.Vect()) ) / Kstarp4_orig.Rho()  ; 
	Float_t K_KstarPl_orig = ( (Kp4_orig.Vect()).Dot(Kstarp4_orig.Vect()) ) / Kstarp4_orig.Rho() ;  
	Float_t alpha_orig = (pi_KstarPl_orig - K_KstarPl_orig)/(pi_KstarPl_orig + K_KstarPl_orig) ; 
	Float_t pi_KstarPt_orig = sqrt( pow(Pip4_orig.Rho(),2) - pow(pi_KstarPl_orig,2) ) ;
	Float_t K_KstarPt_orig = sqrt( pow(Kp4_orig.Rho(),2) - pow(K_KstarPl_orig,2) ) ;
	// refitted tracks
	Float_t pi_KstarPl = ( (pip4.Vect()).Dot(Kstarp4.Vect()) ) / Kstarp4.Rho()  ; 
	Float_t K_KstarPl = ( (kp4.Vect()).Dot(Kstarp4.Vect()) ) / Kstarp4.Rho() ;  
	Float_t alpha = (pi_KstarPl - K_KstarPl)/(pi_KstarPl + K_KstarPl) ; 
	Float_t pi_KstarPt = sqrt( pow(pip4.Rho(),2) - pow(pi_KstarPl,2) ) ;
	//Float_t K_KstarPt = sqrt( pow(kp4.Rho(),2) - pow(K_KstarPl,2) ) ;
	// for exchanged particles
	Float_t pi_ex_KstarPl = ( (Pip4_orig_ex.Vect()).Dot(Kstarp4_orig_ex.Vect()) ) / Kstarp4_orig_ex.Rho()  ; 
	Float_t K_ex_KstarPl = ( (Kp4_orig_ex.Vect()).Dot(Kstarp4_orig_ex.Vect()) ) / Kstarp4_orig_ex.Rho() ;  
	//Float_t alpha_orig_ex = (pi_ex_KstarPl - K_ex_KstarPl)/(pi_ex_KstarPl + K_ex_KstarPl) ; 
	Float_t alpha_orig_ex = (K_ex_KstarPl - pi_ex_KstarPl)/(pi_ex_KstarPl + K_ex_KstarPl) ; 
	Float_t pi_ex_KstarPt = sqrt( pow(Pip4_orig_ex.Rho(),2) - pow(pi_ex_KstarPl,2) ) ;
	Float_t K_ex_KstarPt = sqrt( pow(Kp4_orig_ex.Rho(),2) - pow(K_ex_KstarPl,2) ) ;
	//
	// reco angle variables
	//
	// psi helicity angle
	//theta_psi = getMumuHelicityAngle(myB0p4_orig, Kstarp4_orig, PsiPp4, (*muCharge)[mu1_fromB0_index]<0 ? mu1_p4 : mu2_p4) ;
	//theta_psi = getAleMumuHelicityAngle(myB0p4_orig, Kstarp4_orig, PsiPp4, (*muCharge)[mu1_fromB0_index]<0 ? mu1_p4 : mu2_p4) ;
	theta_psi = GetThetaMuMu(myB0p4_orig, PsiPp4, (*muCharge)[mu1_fromB0_index]>0 ? mu1_p4 : mu2_p4, beam_energy, psi2S_mass, muon_mass) ;
	cosTheta_psi = TMath::Cos(theta_psi) ;
	theta_psi_rot = GetThetaMuMu_rot(myB0p4_orig, PsiPp4, (*muCharge)[mu1_fromB0_index]>0 ? mu1_p4 : mu2_p4, beam_energy, psi2S_mass, muon_mass) ;
	cosTheta_psi_rot = TMath::Cos(theta_psi_rot) ;
	// planes angle
	phi = getPlanesAngle(myB0p4_orig, Kp4_orig, Pip4_orig, (*muCharge)[mu1_fromB0_index]<0 ? mu1_p4 : mu2_p4, (*muCharge)[mu1_fromB0_index]<0 ? mu2_p4 : mu1_p4) ;
	//phi = getBellePlanesAngle(myB0p4_orig, Kp4_orig, Pip4_orig, (*muCharge)[mu1_fromB0_index]<0 ? mu1_p4 : mu2_p4, (*muCharge)[mu1_fromB0_index]<0 ? mu2_p4 : mu1_p4) ;

	myExclusiveMuMuKPiMass->Fill( (*B0Mass)[myB0Idx] ) ; // B0mass before any selection
	hmyPsiPKPiMass->Fill(myB0p4.M());                    // recalculated with 4-vectors
	hmyPsiPKPiMass_KPiExchanged->Fill(myB0p4_KPiExchanged.M()); // recalculated exchanging mass hypotheses
	// putting together for 2 tracks both the recalculated combinations K-Pi & Pi-K
	hmyPsiPKPiMass_bothComb->Fill(myB0p4.M());
	hmyPsiPKPiMass_bothComb->Fill(myB0p4_KPiExchanged.M());
	//
	hPsiPrimefromBMass->Fill( (*MuMuMass)[(*B0MuMuIdx)[myB0Idx]] );
	Float_t PsiPrime_fromB0_Pt = sqrt( pow((*MuMuPx)[(*B0MuMuIdx)[myB0Idx]],2) + pow((*MuMuPy)[(*B0MuMuIdx)[myB0Idx]],2) ) ;
	//
	hB0VtxCL->Fill( (*B0Vtx_CL)[myB0Idx] );
	//
	Float_t B0_Pt = sqrt( pow((*B0Px)[myB0Idx],2) + pow((*B0Py)[myB0Idx],2) ) ;
	hB0_Pt->Fill( B0_Pt );
	//
	Float_t pion_fromB0_Pt = sqrt( pow((*PionPx_MuMuPiK)[myB0Idx],2) + pow((*PionPy_MuMuPiK)[myB0Idx],2) ) ;
	Float_t kaon_fromB0_Pt = sqrt( pow((*KaonPx_MuMuPiK)[myB0Idx],2) + pow((*KaonPy_MuMuPiK)[myB0Idx],2) ) ;
	hPionfromBpt->Fill( pion_fromB0_Pt );
	hKaonfromBpt->Fill( kaon_fromB0_Pt );
	//
	////////////////////////////////
	Bool_t newsoftID = kFALSE ;
	Bool_t PV_matching_cuts = kFALSE, BO_selectionCuts = kFALSE ;
	Bool_t baselineB0Sel = kFALSE, baselinePVMatching = kFALSE ;
	//
	// NEW SOFT MUON 2013
	if ( 1 && 
	     ((*muQual)[mu1_fromB0_index] & (1 << muonQual[3])) && ((*muQual)[mu2_fromB0_index] & (1 << muonQual[3])) && //TMOneSTationTight
	     ((*muChi2)[mu1_fromB0_index]/(*muNDF)[mu1_fromB0_index])<3. && ((*muChi2)[mu2_fromB0_index]/(*muNDF)[mu2_fromB0_index])<3. && 
	     (*muShits)[mu1_fromB0_index]>5 && (*muShits)[mu2_fromB0_index]>5 &&
	     (*muPhits)[mu1_fromB0_index]>0 && (*muPhits)[mu2_fromB0_index]>0 &&                      // new soft muons 2013 (was 1 in old soft muons)
	     fabs((*muDzVtx)[mu1_fromB0_index])<20. && fabs((*muDzVtx)[mu2_fromB0_index])<20. &&        // new soft muons 2013 (was 30 in old soft muons)
	     fabs((*muDxyVtx)[mu1_fromB0_index])<0.3 && fabs((*muDxyVtx)[mu2_fromB0_index])<0.3         // new soft muons 2013 (was 3 in old soft muons)
	     )
	  {
	    newsoftID = true;
	  }
	////////////////////////////////
	Bool_t Mu1_PtSel = false;
	Bool_t Mu2_PtSel = false;
	if ( ( ( (fabs(mu1_p4.Eta()) < 1.2 ) && (mu1_p4.Pt()>4.)) ) || ( ( mu1_p4.Eta() >= 1.2 || mu1_p4.Eta() <=- 1.2) && (mu1_p4.Pt()>3.3) ) ) Mu1_PtSel = true;
	if ( ( ( (fabs(mu2_p4.Eta()) < 1.2 ) && (mu2_p4.Pt()>4.)) ) || ( ( mu2_p4.Eta() >= 1.2 || mu2_p4.Eta() <=- 1.2) && (mu2_p4.Pt()>3.3) ) ) Mu2_PtSel = true;
	////////////////////////////////
	//
	if ( 1
	     && ( (*B0Vtx_CL)[myB0Idx] > 0.2 )
	     && ( (*B0CosAlphaPV)[myB0Idx] > 0.95 )
	     //&& ( (*B0CTauBS)[myB0Idx]/(*B0CTauBSE)[myB0Idx] > 7 )
	     )
	  PV_matching_cuts = kTRUE;
	//
	if ( 1 
	     && ( (*B0Vtx_CL)[myB0Idx] > 0.05 )
	     && ( (*B0CosAlphaPV)[myB0Idx] > 0.99 ) // 0.9
	     //&& ( (*B0CosAlphaPV)[myB0Idx] > 0.99 ) && ( (*B0CosAlpha3DPV)[myB0Idx] > 0.99 )
	     && ( (*B0CTauPV)[myB0Idx]/(*B0CTauPVE)[myB0Idx] > 3.5 ) // 7
	     )
	  BO_selectionCuts = kTRUE ;
	//
	if ( 1 && newsoftID 
	     && ( fabs( (*MuMuMass)[psip_index] - psi2S_mass ) < 0.12 )
	     && ( (*MuMuVtx_CL)[psip_index] > 0.05 )
	     //&& ( pion_fromB0_Pt > 0.5 ) // 0.45
	     //&& ( kaon_fromB0_Pt > 0.5 ) // 0.45
	     && ( B0_Pt > 7.)
	     && ( PsiPrime_fromB0_Pt > 5.)
	     //&& ( mu1_p4.Pt() > 3.3 && mu2_p4.Pt() > 3.3)
	     && Mu1_PtSel && Mu2_PtSel
	     )
	  {
	    if ( BO_selectionCuts ) {
	      hmyPsiPKPiMassBaseSelAlt->Fill( myB0MassAlt );
	      myExclMuMuKPiMassSel->Fill( (*B0Mass)[myB0Idx] ) ;
	      //
	      hmyPsiPKPiMass_bothCombSel->Fill(myB0p4.M());
	      hmyPsiPKPiMass_bothCombSel->Fill(myB0p4_KPiExchanged.M());
	      
	      psi2SPiSqMass_vs_KPiSqMass_ABaseC->Fill( (kp4 + pip4).M2(), (PsiPp4 + pip4).M2() );
	      
	      baselineB0Sel = kTRUE;
	    }
	    if ( PV_matching_cuts )
	      baselinePVMatching = kTRUE;
	  }
	//
	if ( 1 
	     && B0_Pt > 8.
	     //&& pion_fromB0_Pt > 0.7
	     && PsiPp4.DeltaR(Pip4_orig) < 0.5 //1.0
	     && PsiPp4.DeltaR(Kp4_orig) < 0.5 //1.0
	     //
	     && ((*trackChi2)[pi_orig_Index]/(*trackNDF)[pi_orig_Index]) < 5.0
	     && ((*trackChi2)[ka_orig_Index]/(*trackNDF)[ka_orig_Index]) < 5.0
	     && (*trackShits)[pi_orig_Index] >= 10
	     && (*trackShits)[ka_orig_Index] >= 10
	     && (*trackPhits)[pi_orig_Index] > 0
	     && (*trackPhits)[ka_orig_Index] > 0
	     && Pip4_orig.Eta() < 2.4
	     && Kp4_orig.Eta() < 2.4
	     )
	  {
	    //cout <<"before PV matching" <<endl ;
	    // PV matching
	    if ( baselinePVMatching
		 && Pip4_orig.Pt() > 2. && Kp4_orig.Pt() > 2.
		 && fabs(myB0MassAlt - B0_massFit) < 3*B0_sigmaFit
		 ) { 
	      if ( PriVtx_B0CosAlpha_Z )
		priVtx_deltaZ->Fill( fabs((*PriVtx_B0CosAlpha_Z)[myB0Idx] - priVtx_Z) );
	      if ( B0CosAlphaPV && B0CosAlphaB0LessPV )
		B0_pointingAngle_PVvsB0LessPV[0]->Fill( (*B0CosAlphaPV)[myB0Idx], (*B0CosAlphaB0LessPV)[myB0Idx] ) ;
	      if ( B0CosAlphaPVCosAlpha && B0CosAlphaB0LessPVCosAlpha )
		B0_pointingAngle_PVvsB0LessPV[1]->Fill( (*B0CosAlphaPVCosAlpha)[myB0Idx], (*B0CosAlphaB0LessPVCosAlpha)[myB0Idx] ) ;
	      if ( PriVtxB0Less_Z )
		z_PVvsB0LessPV[0]->Fill( priVtx_Z, (*PriVtxB0Less_Z)[myB0Idx] ) ;
	      if ( PriVtx_B0CosAlpha_Z && PriVtxB0Less_B0CosAlpha_Z )
		z_PVvsB0LessPV[1]->Fill( (*PriVtx_B0CosAlpha_Z)[myB0Idx], (*PriVtxB0Less_B0CosAlpha_Z)[myB0Idx] ) ;
	      //
	      B0_PVCTau_h->Fill( (*B0CTauPV)[myB0Idx] / (*B0CTauPVE)[myB0Idx] ) ;
	      B0_CosAlphaCTau_h->Fill( (*B0CTauPVCosAlpha)[myB0Idx] / (*B0CTauPVCosAlphaE)[myB0Idx] ) ;
	      B0_PVXCTau_h->Fill( (*B0CTauPVX)[myB0Idx] / (*B0CTauPVXE)[myB0Idx] ) ;
	      //B0_PVX3DCTau_h->Fill( (*B0CTauPVX_3D)[myB0Idx] / (*B0CTauPVX_3D_err)[myB0Idx] ) ;

	      if ( (priVtx_X == (*PriVtx_B0CosAlpha_X)[myB0Idx]) &&
		   (priVtx_Y == (*PriVtx_B0CosAlpha_Y)[myB0Idx]) &&
		   (priVtx_Z == (*PriVtx_B0CosAlpha_Z)[myB0Idx]) )
		priVtx_matched->Fill( 1 ) ;
	      else 
		priVtx_matched->Fill( 0 ) ;
	      /*
		if ( ((*PriVtxB0Corr_X)[myB0Idx] == (*PriVtx_X)) &&
		     ((*PriVtxB0Corr_Y)[myB0Idx] == (*PriVtx_Y)) &&
		     ((*PriVtxB0Corr_Z)[myB0Idx] == (*PriVtx_Z)) )
		  priVtx_matched->Fill( 1 ) ;
		else 
		  priVtx_matched->Fill( 0 ) ;
	      */
	    }
	    //
	    if ( baselineB0Sel ) {
	      if (PsiPp4.DeltaR(pip4) < 0.4  &&  PsiPp4.DeltaR(kp4) < 0.4)
		hmyPsiPKPiMass_deltaRtest->Fill( myB0MassAlt );
	      if ((*B0CTauPV)[myB0Idx]/(*B0CTauPVE)[myB0Idx] > 11)
		hmyPsiPKPiMass_B0cTautest->Fill( myB0MassAlt );
	      if (B0_Pt > 15.) {
		hmyPsiPKPiMass_B0pTtest->Fill( myB0MassAlt );
		if ( (Pip4_orig + Kp4_orig).M() > 0.85  &&  (Pip4_orig + Kp4_orig).M() < 0.95 )
		  if ( fabs(myB0MassAlt - 5.280808) < 2*0.00816 ) {
		    hPionFromKstar_p->Fill( Pip4_orig.P() );
		    hKaonFromKstar_p->Fill( Kp4_orig.P() );
		  }
	      }
	      nB0AC_noMassWin++ ;
	      
	      kaonP_vs_pionP->Fill( Pip4_orig.P(),Kp4_orig.P() ) ;
	      kaonPt_vs_pionPt->Fill( Pip4_orig.Pt(),Kp4_orig.Pt() ) ;
	      trk_dz->Fill((*trackDzVtx)[pi_orig_Index]) ; trk_dz->Fill((*trackDzVtx)[ka_orig_Index]) ;
	      trk_dxy->Fill((*trackDxyVtx)[pi_orig_Index]) ; trk_dxy->Fill((*trackDxyVtx)[ka_orig_Index]) ;
	      hmyPsiPMassSel->Fill( PsiPp4.M() ); hPsiPrimefromBMassSel->Fill( (*MuMuMass)[(*B0MuMuIdx)[myB0Idx]] ) ;
	      hmyPsiPKPiMassZoom->Fill( myB0p4_orig.M() );
	      hmyPsiPKPiMassAlt->Fill( myB0MassAlt );
	      hmyPsiPKPiMassAltZoom->Fill( myB0MassAlt );
	      // look for reflections M(KK), M(KPi), M(PiPi), M(PsiKK), M(PsiPiPi)
	      hmyPiPiMass->Fill( (Pip4_orig + Kp4_orig_ex).M() ) ;
	      hmyPiKMass->Fill( (Pip4_orig + Kp4_orig).M() ) ;
	      hmyKKMass->Fill( (Pip4_orig_ex + Kp4_orig).M() ) ;
	      hmyPsiPPiPiMass->Fill( (PsiPp4 + Pip4_orig + Kp4_orig_ex).M() ) ;
	      hmyPsiPPiKMass->Fill( (PsiPp4 + Pip4_orig + Kp4_orig).M() ) ; hmyPsiPPiKMass_zoom->Fill( (PsiPp4 + Pip4_orig + Kp4_orig).M() ) ;
	      hmyPsiPKKMass->Fill( (PsiPp4 + Pip4_orig_ex + Kp4_orig).M() ) ;
	      if ( (Pip4_orig_ex + Kp4_orig).M() > phi_left  &&  (Pip4_orig_ex + Kp4_orig).M() < phi_right) {
		hmyPsiPPhiMass->Fill( (PsiPp4 + Pip4_orig_ex + Kp4_orig).M() ) ;
		if ((PsiPp4 + Pip4_orig_ex + Kp4_orig).M() > Bs_left  &&  (PsiPp4 + Pip4_orig_ex + Kp4_orig).M() < Bs_right) { // to tune
		  hmyPsiPPiKMass_fromBs->Fill( (PsiPp4 + Pip4_orig + Kp4_orig).M() ); hmyPsiPPiKMass_fromBs_zoom->Fill( (PsiPp4 + Pip4_orig + Kp4_orig).M() ); }
	      }
	      hmyPiPiMass_vs_KKMass->Fill( (Pip4_orig + Kp4_orig_ex).M(), (Pip4_orig_ex + Kp4_orig).M() ) ;
	      hmyPsiPPiPiMass_vs_PiPiMass->Fill( (PsiPp4 + Pip4_orig + Kp4_orig_ex).M(), (Pip4_orig + Kp4_orig_ex).M() ) ;
	      hmyPsiPKKMass_vs_KKMass->Fill( (PsiPp4 + Pip4_orig_ex + Kp4_orig).M(), (Pip4_orig_ex + Kp4_orig).M() ) ;
	      //
	    
	      // tight cuts
	      for (Int_t i = 0; i < abs(piPt_B0.size()); ++i)
		if (Pip4_orig.Pt() > piPt_B0[i]) 
		  myMuMuKPiMass_piPt[i]->Fill( myB0MassAlt ) ;
	      for (Int_t i = 0; i < abs(KPt_B0.size()); ++i)
		if (Kp4_orig.Pt() > KPt_B0[i]) 
		  myMuMuKPiMass_KPt[i]->Fill( myB0MassAlt ) ;
	      for (Int_t i = 0; i < abs(vtxCL_B0.size()); ++i)
		if ((*B0Vtx_CL)[myB0Idx] > vtxCL_B0[i]) 
		  myMuMuKPiMass_vtxCL[i]->Fill( myB0MassAlt ) ;
	      for (Int_t i = 0; i < abs(cTau_B0.size()); ++i)
		if ((*B0CTauPV)[myB0Idx]/(*B0CTauPVE)[myB0Idx] > cTau_B0[i]) 
		  myMuMuKPiMass_cTau[i]->Fill( myB0MassAlt ) ;
	      // combined
	      for (Int_t i = 0; i < abs(trksPt_B0VtxCL_B0CTau.size()); ++i)
		if ( (Pip4_orig.Pt() > trksPt_B0VtxCL_B0CTau[i].trksPt) && (Kp4_orig.Pt() > trksPt_B0VtxCL_B0CTau[i].trksPt) &&
		     ((*B0Vtx_CL)[myB0Idx] > trksPt_B0VtxCL_B0CTau[i].B0VtxCL) && ((*B0CTauPV)[myB0Idx]/(*B0CTauPVE)[myB0Idx] > trksPt_B0VtxCL_B0CTau[i].B0CTau) ) 
		  myMuMuKPiMass_trksAndB0vtxCLAndB0Ctau[i]->Fill( myB0MassAlt ) ;
	      
	      ///////////////////////////////////////////////////////////////////
	      //cout <<"before nB0AC++" <<endl ;
	      if ( myB0MassAlt > 5.15  &&  myB0MassAlt < 5.45 ) {
		nB0AC++ ;

		psi2SPiSqMass_vs_KPiSqMass_AC->Fill( (kp4 + pip4).M2(), (PsiPp4 + pip4).M2() );

		if ( fabs(myB0MassAlt - B0_massFit) < 2*B0_sigmaFit ) 
		  nB0AC_purityWin++ ;
		
		// MC matching
		//if ( MC ) { // moved below
		Float_t minDeltaRMuMu = 0.05, minDeltaRPi = 0.5, minDeltaRK = 0.5;
		//Float_t minDeltaRMuMu = 0.5, minDeltaRPi = 0.5, minDeltaRK = 0.5; // same result
		//Float_t minDeltaRMuMu = 0.05, minDeltaRPi = 0.9, minDeltaRK = 0.9;  
		//int myMuMuMatched[] = {-1,-1}, myPiMatched[] = {-1,-1}, myKMatched[] = {-1,-1};
		
		vector<int> truePiCharge, trueKCharge ;
		vector<TLorentzVector> truePi4Vec, trueK4Vec, trueMuMu4Vec, trueB04Vec ;   

		TLorentzVector B0p4_gen, Jp4_gen, muPp4_gen, muMp4_gen, Kp4_gen, Pip4_gen, Kstar_gen ;
		
		if ( MC ) {
		  nMCB0_h->Fill( nMCB0 ) ;
		  for (Int_t myMCB0Idx = 0; myMCB0Idx < abs(nMCB0); myMCB0Idx++) { // only 0 and 1 candidate
		    priVtxGen_XYZ[0] = (*PriVtxGen_X)[myMCB0Idx] ; priVtxGen_XYZ[1] = (*PriVtxGen_Y)[myMCB0Idx] ; priVtxGen_XYZ[2] = (*PriVtxGen_Z)[myMCB0Idx] ;
		    genEvtVtx_XYZ[0] = (*genEvtVtx_X)[myMCB0Idx] ; genEvtVtx_XYZ[1] = (*genEvtVtx_Y)[myMCB0Idx] ; genEvtVtx_XYZ[2] = (*genEvtVtx_Z)[myMCB0Idx] ;
		    //      
		    Float_t J_E_gen = 0, K_E_gen = 0, Pi_E_gen = 0, muP_E_gen = 0, muM_E_gen = 0;
		    J_E_gen = sqrt( pow((*MCpsi2SPx)[myMCB0Idx], 2) + pow((*MCpsi2SPy)[myMCB0Idx], 2) + pow((*MCpsi2SPz)[myMCB0Idx], 2) + pow(psi2S_mass, 2) ) ;
		    K_E_gen = sqrt( pow((*MCkaonPx)[myMCB0Idx], 2) + pow((*MCkaonPy)[myMCB0Idx], 2) + pow((*MCkaonPz)[myMCB0Idx], 2) + pow(kaonCh_mass, 2) ) ;
		    Pi_E_gen = sqrt( pow((*MCpionPx)[myMCB0Idx], 2) + pow((*MCpionPy)[myMCB0Idx], 2) + pow((*MCpionPz)[myMCB0Idx], 2) + pow(pionCh_mass, 2) ) ;
		    muP_E_gen = sqrt( pow((*MCmupPx)[myMCB0Idx], 2) + pow((*MCmupPy)[myMCB0Idx], 2) + pow((*MCmupPz)[myMCB0Idx], 2) + pow(muon_mass, 2) ) ;
		    muM_E_gen = sqrt( pow((*MCmumPx)[myMCB0Idx], 2) + pow((*MCmumPy)[myMCB0Idx], 2) + pow((*MCmumPz)[myMCB0Idx], 2) + pow(muon_mass, 2) ) ;
		    //
		    Jp4_gen.SetPxPyPzE( (*MCpsi2SPx)[myMCB0Idx], (*MCpsi2SPy)[myMCB0Idx], (*MCpsi2SPz)[myMCB0Idx], J_E_gen ) ;
		    muPp4_gen.SetPxPyPzE( (*MCmupPx)[myMCB0Idx], (*MCmupPy)[myMCB0Idx], (*MCmupPz)[myMCB0Idx], muP_E_gen ) ;
		    muMp4_gen.SetPxPyPzE( (*MCmumPx)[myMCB0Idx], (*MCmumPy)[myMCB0Idx], (*MCmumPz)[myMCB0Idx], muM_E_gen ) ;
		    Kp4_gen.SetPxPyPzE( (*MCkaonPx)[myMCB0Idx], (*MCkaonPy)[myMCB0Idx], (*MCkaonPz)[myMCB0Idx], K_E_gen ) ;
		    Pip4_gen.SetPxPyPzE( (*MCpionPx)[myMCB0Idx], (*MCpionPy)[myMCB0Idx], (*MCpionPz)[myMCB0Idx], Pi_E_gen ) ;
		    B0p4_gen = Jp4_gen + Pip4_gen + Kp4_gen ;
		    Kstar_gen = Pip4_gen + Kp4_gen ;

		    // NO GEN FILL INSIDE RECO LOOP!
		    /*
		    B0_gen_p_h->Fill( B0p4_gen.P() ) ; B0_gen_pT_h->Fill( B0p4_gen.Pt() ) ;
		    psi2S_gen_p_h->Fill( Jp4_gen.P() ) ; psi2S_gen_pT_h->Fill( Jp4_gen.Pt() ) ;
		    K_gen_p_h->Fill( Kp4_gen.P() ) ; K_gen_pT_h->Fill( Kp4_gen.Pt() ) ;
		    pi_gen_p_h->Fill( Pip4_gen.P() ) ; pi_gen_pT_h->Fill( Pip4_gen.Pt() ) ;

		    psi2SPiSqMass_vs_KPiSqMass_gen->Fill( (Kp4_gen + Pip4_gen).M2(), (Jp4_gen + Pip4_gen).M2() );
		    */
		    truePi4Vec.push_back( Pip4_gen ); trueK4Vec.push_back( Kp4_gen ); trueMuMu4Vec.push_back( Jp4_gen ); trueB04Vec.push_back( B0p4_gen );
		    truePiCharge.push_back( (*MCpionCh)[myMCB0Idx] ); trueKCharge.push_back( (*MCkaonCh)[myMCB0Idx] ); 
		  }
		  if (nMCB0 > 0)
		    B0Mass_gen = B0p4_gen.M() ;
		}
		
		if ( !MC || (MC && nMCB0>0) ) {

		  if ( (nB0AC == 1 ) ) {
		    
		    pB0_XYZ[0] = (*B0Px)[myB0Idx] ; pB0_XYZ[1] = (*B0Py)[myB0Idx] ; pB0_XYZ[2] = (*B0Pz)[myB0Idx] ;
		    B0CTau_2B0[0][0] = (*B0CTauPV)[myB0Idx]/(*B0CTauPVE)[myB0Idx] ;
		    secVtx_XYZ[0] = (*B0DecayVtx_X)[myB0Idx] ; secVtx_XYZ[1] = (*B0DecayVtx_Y)[myB0Idx] ; secVtx_XYZ[2] = (*B0DecayVtx_Z)[myB0Idx] ;
		    //if ( priVtx_Z )
		    priVtx_XYZ[0] = priVtx_X; priVtx_XYZ[1] = priVtx_Y; priVtx_XYZ[2] = priVtx_Z;
		    priVtx_XYZ[3] = (*B0CosAlphaPV)[myB0Idx] ; //priVtx_XYZ[4] = (*B0CosAlpha3DPV)[myB0Idx] ;

		    cout <<"before PriVtx stuff" <<endl ;
		    if ( PriVtx_B0CosAlpha_Z ) {
		      priVtx_B0CosAlpha_XYZ[0][0] = (*PriVtx_B0CosAlpha_X)[myB0Idx]; priVtx_B0CosAlpha_XYZ[0][1] = (*PriVtx_B0CosAlpha_Y)[myB0Idx]; priVtx_B0CosAlpha_XYZ[0][2] = (*PriVtx_B0CosAlpha_Z)[myB0Idx] ; 
		      priVtx_B0CosAlpha_XYZ[0][3] = (*B0CosAlphaPVCosAlpha)[myB0Idx]; priVtx_B0CosAlpha_XYZ[0][4] = 0; // to add (*B0CosAlphaPVCosAlpha3D)[myB0Idx];
		    }
		    if ( PriVtx_B0CosAlpha3D_Z ) {
		      priVtx_B0CosAlpha_XYZ[1][0] = (*PriVtx_B0CosAlpha3D_X)[myB0Idx]; priVtx_B0CosAlpha_XYZ[1][1] = (*PriVtx_B0CosAlpha3D_Y)[myB0Idx]; priVtx_B0CosAlpha_XYZ[1][2] = (*PriVtx_B0CosAlpha3D_Z)[myB0Idx] ; 
		      priVtx_B0CosAlpha_XYZ[1][3] = 0; priVtx_B0CosAlpha_XYZ[1][4] = (*B0CosAlphaPVCosAlpha3D)[myB0Idx]; // to add (*B0CosAlphaPVCosAlpha2D)[myB0Idx];
		    }
		    if ( PriVtxB0Less_Z ) {
		      priVtxB0Less_XYZ[0] = (*PriVtxB0Less_X)[myB0Idx]; priVtxB0Less_XYZ[1] = (*PriVtxB0Less_Y)[myB0Idx]; priVtxB0Less_XYZ[2] = (*PriVtxB0Less_Z)[myB0Idx]; 
		      priVtxB0Less_XYZ[3] = (*B0CosAlphaB0LessPV)[myB0Idx] ; priVtxB0Less_XYZ[4] = 0; // to add (*B0CosAlpha3DB0LessPV)[myB0Idx] ;
		    }
		    if ( PriVtxB0Less_B0CosAlpha_Z ) {
		      priVtxB0Less_B0CosAlpha_XYZ[0][0] = (*PriVtxB0Less_B0CosAlpha_X)[myB0Idx]; priVtxB0Less_B0CosAlpha_XYZ[0][1] = (*PriVtxB0Less_B0CosAlpha_Y)[myB0Idx]; priVtxB0Less_B0CosAlpha_XYZ[0][2] = (*PriVtxB0Less_B0CosAlpha_Z)[myB0Idx]; 
		      priVtxB0Less_B0CosAlpha_XYZ[0][3] = (*B0CosAlphaB0LessPVCosAlpha)[myB0Idx]; priVtxB0Less_B0CosAlpha_XYZ[0][4] = 0; // to add (*B0CosAlphaB0LessPVCosAlpha3D)[myB0Idx];
		    }
		    if ( PriVtxB0Less_B0CosAlpha3D_Z ) {
		      priVtxB0Less_B0CosAlpha_XYZ[1][0] = (*PriVtxB0Less_B0CosAlpha3D_X)[myB0Idx]; priVtxB0Less_B0CosAlpha_XYZ[1][1] = (*PriVtxB0Less_B0CosAlpha3D_Y)[myB0Idx]; priVtxB0Less_B0CosAlpha_XYZ[1][2] = (*PriVtxB0Less_B0CosAlpha3D_Z)[myB0Idx]; 
		      priVtxB0Less_B0CosAlpha_XYZ[1][3] = 0; priVtxB0Less_B0CosAlpha_XYZ[1][4] = (*B0CosAlphaB0LessPVCosAlpha3D)[myB0Idx]; // to add (*B0CosAlphaB0LessPVCosAlpha2D)[myB0Idx];
		    }
		       
 
		    B0Mass_2B0[0][nB0AC-1] = myB0MassAlt ;
		    PsiPp4_2B0[0][nB0AC-1] = PsiPp4; Pip4_2B0[0][nB0AC-1] = Pip4_orig; Kp4_2B0[0][nB0AC-1] = Kp4_orig ;
		    PiCharge_2B0[0][nB0AC-1] = (*trackCharge)[pi_orig_Index] ; KCharge_2B0[0][nB0AC-1] = (*trackCharge)[ka_orig_Index] ;
		    //
		    twins_var[0][0] = PsiPp4.M(); twins_var[0][1] = PsiPp4.Pt(); twins_var[0][2] = Kp4_orig.Pt(); twins_var[0][3] = Pip4_orig.Pt();

		    if ( MC ) {
		      DeltaR_1B0[0] = myB0p4_orig.DeltaR(trueB04Vec[0]) ;
		      DeltaPt_1B0[0] = fabs(myB0p4_orig.Pt() - trueB04Vec[0].Pt()) / trueB04Vec[0].Pt() ;
		      
		      Double_t DrMuMu = trueMuMu4Vec[0].DeltaR(PsiPp4) ; // DeltaR between true psi' and reco psi'
		      Double_t MuMuPtDiff = fabs(trueMuMu4Vec[0].Pt() - PsiPp4.Pt())/trueMuMu4Vec[0].Pt() ;
		      Double_t DrPi = truePi4Vec[0].DeltaR(Pip4_orig) ; // DeltaR between true pi and reco pi
		      Double_t piPtDiff = fabs(truePi4Vec[0].Pt() - Pip4_orig.Pt())/truePi4Vec[0].Pt() ;
		      Double_t DrK = trueK4Vec[0].DeltaR(Kp4_orig) ; // DeltaR between true K and reco K
		      Double_t kPtDiff = fabs(trueK4Vec[0].Pt() - Kp4_orig.Pt())/trueK4Vec[0].Pt() ;
		      
		      chargePi_1B0[0] = truePiCharge[0] * PiCharge_2B0[0][nB0AC-1] ;
		      chargeK_1B0[0] = trueKCharge[0] * KCharge_2B0[0][nB0AC-1] ;
		      
		      // check truth matching
		      if ( (DrMuMu < minDeltaRMuMu) && (MuMuPtDiff < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			MuMuMatched_1B0[0] = kTRUE ;
		      
		      if ( chargePi_1B0[0] > 0 ) 
			if ( (DrPi < minDeltaRPi) && (piPtDiff < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			  pionMatched_1B0[0] = kTRUE ;
		      
		      if ( chargeK_1B0[0] > 0 ) 
			if ( (DrK < minDeltaRK) && (kPtDiff < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			  kaonMatched_1B0[0] = kTRUE ;
		      
		    }
		  }
		  else if ( nB0AC == 2 ) {
		    B0Mass_2B0[0][nB0AC-1] = myB0MassAlt ;
		    PsiPp4_2B0[0][nB0AC-1] = PsiPp4; Pip4_2B0[0][nB0AC-1] = Pip4_orig; Kp4_2B0[0][nB0AC-1] = Kp4_orig ;
		    PiCharge_2B0[0][nB0AC-1] = (*trackCharge)[pi_orig_Index] ; KCharge_2B0[0][nB0AC-1] = (*trackCharge)[ka_orig_Index] ;
		    //
		    twins_var[1][0] = PsiPp4.M(); twins_var[1][1] = PsiPp4.Pt(); twins_var[1][2] = Kp4_orig.Pt(); twins_var[1][3] = Pip4_orig.Pt();
		    if (twins_var[0][0] == twins_var[1][0]) // psi(2S) mass
		      if (twins_var[0][1] == twins_var[1][1]) // psi(2S) pT
			if (twins_var[0][2] == twins_var[1][3]) // track pT
			  if (twins_var[0][3] == twins_var[1][2]) // track pT
			    nTwins++ ;
		    //
		    if ( MC ) {
		      DeltaR_2B0[0][nB0AC-2] = DeltaR_1B0[0] ;
		      DeltaPt_2B0[0][nB0AC-2] = DeltaPt_1B0[0] ;
		      DeltaR_2B0[0][nB0AC-1] = myB0p4_orig.DeltaR(trueB04Vec[0]) ; 
		      DeltaPt_2B0[0][nB0AC-1] = fabs(myB0p4_orig.Pt() - trueB04Vec[0].Pt()) / trueB04Vec[0].Pt() ; 
		      
		      Double_t DrMuMu[2] = { trueMuMu4Vec[0].DeltaR(PsiPp4_2B0[0][0]), trueMuMu4Vec[0].DeltaR(PsiPp4_2B0[0][1]) } ; // DeltaR between true psi' and reco psi'
		      Double_t MuMuPtDiff[2] = { fabs(trueMuMu4Vec[0].Pt() - PsiPp4_2B0[0][0].Pt())/trueMuMu4Vec[0].Pt(), fabs(trueMuMu4Vec[0].Pt() - PsiPp4_2B0[0][1].Pt())/trueMuMu4Vec[0].Pt() } ;
		      Double_t piPtDiff[2] = { fabs(truePi4Vec[0].Pt() - Pip4_2B0[0][0].Pt())/truePi4Vec[0].Pt(), fabs(truePi4Vec[0].Pt() - Pip4_2B0[0][1].Pt())/truePi4Vec[0].Pt() } ;
		      Double_t kPtDiff[2] = { fabs(trueK4Vec[0].Pt() - Kp4_2B0[0][0].Pt())/trueK4Vec[0].Pt(), fabs(trueK4Vec[0].Pt() - Kp4_2B0[0][1].Pt())/trueK4Vec[0].Pt() } ;
		    
		      for (Int_t i = 0; i < nB0AC; i++) { 
			DrPi_2B0[i] = truePi4Vec[0].DeltaR(Pip4_2B0[0][i]) ; // DeltaR between true psi' and reco psi'
			DrK_2B0[i] = trueK4Vec[0].DeltaR(Kp4_2B0[0][i]) ; // DeltaR between true psi' and reco psi'
			chargePi_2B0[0][i] = truePiCharge[0] * PiCharge_2B0[0][i] ;
			chargeK_2B0[0][i] = trueKCharge[0] * KCharge_2B0[0][i] ;

			if ( (DrMuMu[i] < minDeltaRMuMu) && (MuMuPtDiff[i] < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			  MuMuMatched[0][i] = kTRUE ;
		      
			if ( chargePi_2B0[0][i] > 0 ) 
			  if ( (DrPi_2B0[i] < minDeltaRPi) && (piPtDiff[i] < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			    pionMatched[0][i] = kTRUE ;
		      
			if ( chargeK_2B0[0][i] > 0 ) 
			  if ( (DrK_2B0[i] < minDeltaRK) && (kPtDiff[i] < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			    kaonMatched[0][i] = kTRUE ;
		      }
		    }
		  } // if ( nB0AC == 2 ) {
		} // if ( !MC || (MC && nMCB0>0) ) {
	      
		cout <<"before nB0AC_signalWin++" <<endl ;
		if ( fabs(myB0MassAlt - B0_massFit) < 3*B0_sigmaFit ) {
		  nB0AC_signalWin++ ;
		
		  if ( MC ) {
		    if ( nMCB0 > 0 ) {
		      nB0ACInMC++ ;
		      B0_gen_mass_h->Fill( trueB04Vec[0].M() ) ;		  
		      piCh_vs_KCh_gen->Fill(truePiCharge[0], trueKCharge[0]) ;
		    }
		  } 

		  if ( (!MC && (nB0AC_signalWin == 1)) || nB0ACInMC == 1 ) {

		    B0Mass_2B0[1][nB0ACInMC-1] = myB0MassAlt ;
		    B0VtxCL_2B0[0] = (*B0Vtx_CL)[myB0Idx] ;
		    B0CTau_2B0[1][0] = (*B0CTauPV)[myB0Idx]/(*B0CTauPVE)[myB0Idx] ;
		    alpha_2B0[0] = alpha_orig;
		    PsiPp4_2B0[1][0] = PsiPp4; Pip4_2B0[1][0] = Pip4_orig; Kp4_2B0[1][0] = Kp4_orig ;
		    PiCharge_2B0[1][0] = (*trackCharge)[pi_orig_Index] ; KCharge_2B0[1][0] = (*trackCharge)[ka_orig_Index] ;
		    //
		    twins_var[0][0] = PsiPp4.M(); twins_var[0][1] = PsiPp4.Pt(); twins_var[0][2] = Kp4_orig.Pt(); twins_var[0][3] = Pip4_orig.Pt();

		    if ( MC ) { // truth matching
		      if ( nMCB0 < 2 ) {   // skip events with more than 1 generated B0		  

			DeltaR_1B0[1] = myB0p4_orig.DeltaR(trueB04Vec[0]) ;
			DeltaPt_1B0[1] = fabs(myB0p4_orig.Pt() - trueB04Vec[0].Pt()) / trueB04Vec[0].Pt() ;
			chargePi_1B0[1] = truePiCharge[0] * PiCharge_2B0[1][0] ;
			chargeK_1B0[1] = trueKCharge[0] * KCharge_2B0[1][0] ;
			chargePi_2B0[1][0] = chargePi_1B0[1] ;
			chargeK_2B0[1][0] = chargeK_1B0[1] ;
			  
			Double_t DrMuMu = trueMuMu4Vec[0].DeltaR(PsiPp4) ; // DeltaR between true psi' and reco psi'
			Double_t MuMuPtDiff = fabs(trueMuMu4Vec[0].Pt() - PsiPp4.Pt())/trueMuMu4Vec[0].Pt() ;
			Double_t DrPi = truePi4Vec[0].DeltaR(Pip4_orig) ; // DeltaR between true psi' and reco psi'
			Double_t piPtDiff = fabs(truePi4Vec[0].Pt() - Pip4_orig.Pt())/truePi4Vec[0].Pt() ;
			Double_t DrK = trueK4Vec[0].DeltaR(Kp4_orig) ; // DeltaR between true psi' and reco psi'
			Double_t kPtDiff = fabs(trueK4Vec[0].Pt() - Kp4_orig.Pt())/trueK4Vec[0].Pt() ;
		      
			// check truth matching
			if ( (DrMuMu < minDeltaRMuMu) && (MuMuPtDiff < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			  MuMuMatched_1B0[1] = kTRUE ;
			  
			if ( chargePi_1B0[1] > 0 ) 
			  if ( (DrPi < minDeltaRPi) && (piPtDiff < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			    pionMatched_1B0[1] = kTRUE ;
			  
			if ( chargeK_1B0[1] > 0 ) 
			  if ( (DrK < minDeltaRK) && (kPtDiff < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			    kaonMatched_1B0[1] = kTRUE ;

		      }
		    }
		  
		  } // if ( (!MC && (nB0AC_signalWin == 1)) || nB0ACInMC == 1 )
		  else if ( (!MC && (nB0AC_signalWin == 2)) || nB0ACInMC == 2 ) {
		    B0Mass_2B0[1][nB0ACInMC-1] = myB0MassAlt ;
		    B0VtxCL_2B0[1] = (*B0Vtx_CL)[myB0Idx] ;
		    B0CTau_2B0[1][1] = (*B0CTauPV)[myB0Idx]/(*B0CTauPVE)[myB0Idx] ;
		    alpha_2B0[1] = alpha_orig;
		    PsiPp4_2B0[1][1] = PsiPp4; Pip4_2B0[1][1] = Pip4_orig; Kp4_2B0[1][1] = Kp4_orig ;
		    PiCharge_2B0[1][1] = (*trackCharge)[pi_orig_Index] ; KCharge_2B0[1][1] = (*trackCharge)[ka_orig_Index] ;
		    //
		    twins_var[1][0] = PsiPp4.M(); twins_var[1][1] = PsiPp4.Pt(); twins_var[1][2] = Kp4_orig.Pt(); twins_var[1][3] = Pip4_orig.Pt();
		    if (twins_var[0][0] == twins_var[1][0]) // psi(2S) mass
		      if (twins_var[0][1] == twins_var[1][1]) // psi(2S) pT
			if (twins_var[0][2] == twins_var[1][3]) // track pT
			  if (twins_var[0][3] == twins_var[1][2]) // track pT
			    nTwins_signalWin++ ;
		
		    if ( MC ) { // truth matching
		      if ( nMCB0 < 2 ) {   // skip events with more than 1 generated B0		  
			DeltaR_2B0[1][nB0ACInMC-2] = DeltaR_1B0[1] ;
			DeltaPt_2B0[1][nB0ACInMC-2] = DeltaPt_1B0[1] ;
			DeltaR_2B0[1][nB0ACInMC-1] = myB0p4_orig.DeltaR(trueB04Vec[0]) ; 
			DeltaPt_2B0[1][nB0ACInMC-1] = fabs(myB0p4_orig.Pt() - trueB04Vec[0].Pt()) / trueB04Vec[0].Pt() ; 
			chargePi_2B0[1][1] = truePiCharge[0] * PiCharge_2B0[1][1] ;
			chargeK_2B0[1][1] = trueKCharge[0] * KCharge_2B0[1][1] ;
		      
			Double_t DrMuMu[2] = { trueMuMu4Vec[0].DeltaR(PsiPp4_2B0[1][0]), trueMuMu4Vec[0].DeltaR(PsiPp4_2B0[1][1]) } ; // DeltaR between true psi' and reco psi'
			Double_t MuMuPtDiff[2] = { fabs(trueMuMu4Vec[0].Pt() - PsiPp4_2B0[1][0].Pt())/trueMuMu4Vec[0].Pt(), fabs(trueMuMu4Vec[0].Pt() - PsiPp4_2B0[1][1].Pt())/trueMuMu4Vec[0].Pt() } ;
			Double_t DrPi[2] = { truePi4Vec[0].DeltaR(Pip4_2B0[1][0]), truePi4Vec[0].DeltaR(Pip4_2B0[1][1]) } ; // DeltaR between true psi' and reco psi'
			Double_t piPtDiff[2] = { fabs(truePi4Vec[0].Pt() - Pip4_2B0[1][0].Pt())/truePi4Vec[0].Pt(), fabs(truePi4Vec[0].Pt() - Pip4_2B0[1][1].Pt())/truePi4Vec[0].Pt() } ;
			Double_t DrK[2] = { trueK4Vec[0].DeltaR(Kp4_2B0[1][0]), trueK4Vec[0].DeltaR(Kp4_2B0[1][1]) } ; // DeltaR between true psi' and reco psi'
			Double_t kPtDiff[2] = { fabs(trueK4Vec[0].Pt() - Kp4_2B0[1][0].Pt())/trueK4Vec[0].Pt(), fabs(trueK4Vec[0].Pt() - Kp4_2B0[1][1].Pt())/trueK4Vec[0].Pt() } ;
		      
			//if ( nTwins_signalWin == 1 ) { // let's match any 2 B0 candidates

			// check truth matching
			for (Int_t i = 0; i < nB0AC_signalWin+1; i++) { // if ( (!MC && (nB0AC_signalWin == 2)) || nB0ACInMC == 2 )

			  if ( (DrMuMu[i] < minDeltaRMuMu) && (MuMuPtDiff[i] < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			    MuMuMatched[1][i] = kTRUE ;
			  
			  if ( chargePi_2B0[1][i] > 0 ) 
			    if ( (DrPi[i] < minDeltaRPi) && (piPtDiff[i] < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			      pionMatched[1][i] = kTRUE ;
			  
			  if ( chargeK_2B0[1][i] > 0 ) 
			    if ( (DrK[i] < minDeltaRK) && (kPtDiff[i] < momentumCut) ) // matching based on DeltaR and below momentumCut relative difference in pT
			      kaonMatched[1][i] = kTRUE ;
			  
			  // fill histograms
			  if ( MuMuMatched[1][i] ) { // always true
			    /*
			      if ( pionMatched[1][i] ) { // matched pions
			      hMCDeltaRPiMatch->Fill( minDeltaRPi[pit] );
			      hMCPtRecoVsPtTruePiMatch->Fill(truePiPt, matchedPiPt);
			      hMCPtTruePiMatch->Fill( matchedPiPt );
			      if ( matchedStrip > 0 ) {
			      hMCStripsHitsTruePiMatch->Fill( matchedStrip );  
			      hMCPixelHitsTruePiMatch->Fill( matchedPixel ); 
			      hMCTotalHitsTruePiMatch->Fill( matchedStrip+matchedPixel );  
			      hMCChiTruePiMatch->Fill( matchedTrackChiNDF ); 
			      }
			      } else { // not matched pions
			      //for (UInt_t it = 0; it < trackPx->size() ; it ++) {
			      //if ((abs(it) != myPiMatched[0])) {
			      //if ((abs(it) != myPiMatched[0]) && (it != myPiMatched[1])) {
			      hMCStripsHitsTrackNotMached->Fill( (*trackShits)[it] );  
			      hMCPixelHitsTrackNotMached->Fill( (*trackPhits)[it] );  
			      hMCTotalHitsTrackNotMached->Fill( (*trackPhits)[it] + (*trackShits)[it] ); 
			      hMCChiTrackNotMached->Fill(1.0*(*trackChi2)[it] / 1.0*(*trackNDF)[it]);
			      }
			    */
			    // matched kaons
			    hMCDeltaRK_fromTwins->Fill( DrK[i] );
			    hMCPtRecoVsPtTrueK_fromTwins->Fill( trueK4Vec[0].Pt(), Kp4_2B0[1][i].Pt());
			    hMCDeltaPtK_fromTwins->Fill( kPtDiff[i] ) ;
			    if ( pionMatched[1][i] )
			      hMCDeltaPtK_fromTwinsPiMatched->Fill( kPtDiff[i] ) ;
			    hMCpTTrueK_fromTwins->Fill( Kp4_2B0[1][i].Pt() );
			    
			    /*
			      weight = GetWeightPtShape(xGenPt,histo_weight);
			      weight_plus = GetWeightPtShape(xGenPt,histo_weightplus);
			      weight_minus = GetWeightPtShape(xGenPt,histo_weightminus);
			    */
			    
			    if ( MuMuMatched[1][i] && pionMatched[1][i] && kaonMatched[1][i] ) {
			      kaonDr_vs_pionDr_matched->Fill( minDeltaRPi, minDeltaRK ) ;
			      kaonP_vs_pionP_matched->Fill( Pip4_2B0[1][i].P(), Kp4_2B0[1][i].P() ) ;
			      kaonPt_vs_pionPt_matched->Fill( Pip4_2B0[1][i].Pt(), Kp4_2B0[1][i].Pt() ) ;
			      //hMCDeltaRB0Match->Fill( Kp4_orig.DeltaR(trueK4Vec[kit]) ) ;
			    }
			    
			  } // if ( MuMuMatched[1][i] )
			} // for (Int_t i=0; i<2; i++)
			//} // if (nTwins_signalWin == 1)
		      } // if (nMCB0 < 2)
		    } // if MC

		  } // if (nB0ACInMC == 2)
		  //} // if (nMCB0 > 0) removed to use data too
		

		  ///////////////////////////////////////////////////////////////////
		  /*
		  // Daniele code //
		  //float weight = 1., weight_plus = 1., weight_minus =1.;
		  //float XTruePdl = -10.;
		  //int TrueMother = 0;
	    
		  int nGenX = 0;
		  float xGenEta=999, xGenPhi=999, xGenPt=999, xGenP=999, xGenY=999;
		  float psi2SGenY=999, psi2SGenPt=999;
		  float kTruePt = 999, kTrueRap =999;
		  double GenDeltaPhiMu = 999.;
		  vector<TLorentzVector> trueMu4Vec, truePi4Vec;
		  TLorentzVector xGen4Mom, psi2SGen4Mom, kTrueVec;
		  vector<int> truePiCharge;
		  vector<bool> trueMuInAcc;
		  int myPiMatched[] = {-1,-1};
		  int myKMatched = -1;
		  float minDeltaR[] = {0.5, 0.5};
		  float minDrK = 0.5;
    
		  //float partMass = (partNum == psi2S_PDGId) ? psi2S_mass : X3872PDG;
		  //float minMassRange = (partNum == psi2S_PDGId) ? 3.656 : 3.827;                // 3.855 : 3.665; Sideband
		  //float maxMassRange = (partNum == psi2S_PDGId) ? 3.716 : 3.917;               // 3.885  : 3.707; Sidebands

		  int knumber = 0;
    
		  /////////////// Get generated particles MC Info ////////////////////////////////
		  for (UInt_t it = 0; it < nMCAll ; it++) { // loop on gen particles
		  int pId = (*MCPdgIdAll)[it]; 
		  // cout << " it "<< it << " pdi " << pId << endl;
		  if (pId == B0_Id) {
		  nGenX++;
		  //TLorentzVector xGen4Mom;
		  xGen4Mom.SetXYZM((*MCPx)[it],(*MCPy)[it],(*MCPz)[it], B0_mass);
		  //xGenEta = xGen4Mom.Eta(); xGenPhi = xGen4Mom.Phi(); xGenPt = xGen4Mom.Pt(); xGenP = xGen4Mom.P(); xGenY =  xGen4Mom.Rapidity();
		  //XTruePdl = 10*(*ppdlTrue)[it]; //?
		  //TrueMother = (*MC_Real_Mom)[it]; //?
		  } 
		  else if(abs(pId) == muon_Id) { // muon
		  TLorentzVector tempMu4Vec;
		  tempMu4Vec.SetXYZM((*MCPx)[it],(*MCPy)[it],(*MCPz)[it], muon_mass);
		  trueMu4Vec.push_back(tempMu4Vec);
		  bool inMuAcc = isMuonInAccept(fabs(tempMu4Vec.Eta()), tempMu4Vec.Pt(), tempMu4Vec.P());
		  trueMuInAcc.push_back(inMuAcc);
		  }
		  else if(abs(pId) == pion_Id) {       // pion
		  TLorentzVector tempPi4Vec;
		  tempPi4Vec.SetXYZM((*MCPx)[it],(*MCPy)[it],(*MCPz)[it], pionCh_mass);
		  truePi4Vec.push_back( tempPi4Vec );
		  int charge = (pId == pion_Id)? 1 : -1;
		  truePiCharge.push_back(charge);
		  }
		  else if(abs(pId) == psi2S_Id) {   // psi(2S)
		  psi2SGen4Mom.SetXYZM((*MCPx)[it],(*MCPy)[it],(*MCPz)[it], psi2S_mass);
		  //psi2SGenPt = psi2SGen4Mom.Pt(); psi2SGenY = psi2SGen4Mom.Rapidity();
		  }
		  else if(abs(pId) == Kaon_Id) {   // K
		  knumber ++;
		  kTrueVec.SetXYZM((*MCPx)[it],(*MCPy)[it],(*MCPz)[it], kaonCh_mass);
		  //kTruePt = kTrueVec.Pt(); kTrueRap = kTrueVec.Rapidity();
		  }
	
		  } // end loop on gen particles
      
		  //GenDeltaPhiMu = trueMu4Vec[0].DeltaPhi(trueMu4Vec[1]);
      
		  if (nGenX < 2) {   // skip events with more than one generated X (X = B0 in our case)

		  UInt_t pionNumber_min = 1 ;
		  if (truePiCharge.size() >= pionNumber_min) { // skip events without pionNumber_min generated pions

		  hMCnGenX->Fill( nGenX );
	
		  // Match Pions Gen-Reco done with original pions, before fit
		  float momCut = 0.1; // no more than 10% difference in momentum
		  for (UInt_t pit = 0;  pit < pionNumber_min; pit ++) {
		  int matchedPixel = -1;
		  int matchedStrip = -1;
		  int matchedTrackChiNDF = 0;
		  float truePiP = (truePi4Vec[pit]).P();
		  float matchedPiP = 999.;
		  for (UInt_t it = 0; it < trackPx->size() ; it ++) { // loop on reco track
	    
		  if( (truePiCharge[pit] * (*trackCharge)[it]) < 0) continue;  // check same charge
		  TLorentzVector recPi4Mom;
		  recPi4Mom.SetPxPyPzE((*trackPx)[it],(*trackPy)[it],(*trackPz)[it],(*trackEnergy)[it]);
		  float dr = recPi4Mom.DeltaR(truePi4Vec[pit]); // DeltaR between true pion and reco pion
		  float recoPiP = recPi4Mom.P();
		  float piPdiff = fabs(truePiP - recoPiP) / truePiP;
		  if( (dr < minDeltaR[pit]) && (piPdiff < momCut)) { // matching based on DeltaR and below 10% difference in P
		  minDeltaR[pit] = dr;
		  myPiMatched[pit] = it;
		  matchedPiP = recoPiP;
		  matchedPixel = (*trackPhits)[it];
		  matchedStrip = (*trackShits)[it];
		  matchedTrackChiNDF = ((*trackChi2)[it] / (*trackNDF)[it]);
		  }
	    
		  }
		  hMCDeltaRPiMatch->Fill( minDeltaR[pit] );
		  hMCPRecoVsPTruePiMatch->Fill(truePiP, matchedPiP);
		  if ( matchedStrip > 0 ) {
		  hMCStripsHitsTruePiMatch->Fill( matchedStrip );  
		  hMCPixelHitsTruePiMatch->Fill( matchedPixel ); 
		  hMCTotalHitsTruePiMatch->Fill( matchedStrip+matchedPixel );  
		  hMCChiTruePiMatch->Fill( matchedTrackChiNDF ); 
		  }
		  }
	
		  // not matched pions
		  for (UInt_t it = 0; it < trackPx->size() ; it ++) {
		  //if ((abs(it) != myPiMatched[0]) && (it != myPiMatched[1])) {
		  if ((abs(it) != myPiMatched[0])) {
		  hMCStripsHitsTrackNotMached->Fill( (*trackShits)[it] );  
		  hMCPixelHitsTrackNotMached->Fill( (*trackPhits)[it] );  
		  hMCTotalHitsTrackNotMached->Fill( (*trackPhits)[it] + (*trackShits)[it] ); 
		  hMCChiTrackNotMached->Fill(1.0*(*trackChi2)[it] / 1.0*(*trackNDF)[it]);
		  }
		  }
	
		  if (knumber > 0) {
		  int matchedPixelK = -1, matchedStripK = -1;
		  float matchedKP;
		  //minDrK = 0.5;
		  for (UInt_t it = 0; it < trackPx->size() ; it ++) { // loop on reco track
	    
		  TLorentzVector recoKMom;
		  recoKMom.SetPxPyPzE((*trackPx)[it],(*trackPy)[it],(*trackPz)[it],(*trackEnergy)[it]);
		  float drk = recoKMom.DeltaR(kTrueVec); // DeltaR between true kaon and reco kaon
		  float recokP = recoKMom.P();
		  float kPdiff = fabs(kTrueVec.P() - recokP) / kTrueVec.P();
		  if( (drk < minDrK) &&  (kPdiff < momCut)) { // matching based on DeltaR and below 10% difference in P
		  minDrK = drk;
		  myKMatched = it;
		  matchedKP = recokP;
		  matchedPixelK = (*trackPhits)[it];
		  matchedStripK = (*trackShits)[it];
		  //matchedTrackChiNDF = ((*TrackChi2)[it] / (*TrackNDF)[it]);
		  }
	    
		  }
	  
		  }    
	
		  //weight = GetWeightPtShape(xGenPt,histo_weight);
		  //weight_plus = GetWeightPtShape(xGenPt,histo_weightplus);
		  //weight_minus = GetWeightPtShape(xGenPt,histo_weightminus);
	
		  } // if (truePiCharge.size() >= pionNumber_min)
		  } // if (nGenX < 2)
		  */ // Daniele code //

		  //} // if (MC) // removed to use data too
		} // if ( fabs(myB0MassAlt - B0mass) < 3*sigma ) {
	      } // if ( myB0MassAlt > 5.15  &&  myB0MassAlt < 5.45 ) {	    
	      ////////////////////////////////////////////////////////////////////////////////////
	    
	      cout <<"before K*" <<endl ;
	      // K* 
	      if (myB0MassAlt > 5.25  &&  myB0MassAlt < 5.3) {
		hmyKPiMass->Fill( (Pip4_orig + Kp4_orig).M() ) ;
		hmyKPiMass_ex->Fill( (Pip4_orig_ex + Kp4_orig_ex).M() ) ;
		if (Kp4_orig.Pt() < 1.) {
		  hmyKPiMass_KpTcut->Fill( (Pip4_orig + Kp4_orig).M() ) ;
		  hmyKPiMass_ex_KpTcut->Fill( (Pip4_orig_ex + Kp4_orig_ex).M() ) ;
		  // high purity
		  if ((*TrackTight)[pi_orig_Index] && (*TrackTight)[pi_orig_Index]) {

		    hmyKPiMass_KpTcut_tightTr->Fill( (Pip4_orig + Kp4_orig).M() ) ;
		    hmyKPiMass_ex_KpTcut_tightTr->Fill( (Pip4_orig_ex + Kp4_orig_ex).M() ) ;
		    if ((*TrackHighPurity)[pi_orig_Index] && (*TrackHighPurity)[ka_orig_Index]) {

		      hmyKPiMass_KpTcut_HP->Fill( (Pip4_orig + Kp4_orig).M() ) ;
		      hmyKPiMass_ex_KpTcut_HP->Fill( (Pip4_orig_ex + Kp4_orig_ex).M() ) ;
		    }
		  }
		}
	      } // B0 mass window (5.25,5.3)
	      if ( (Pip4_orig + Kp4_orig).M() > 0.85  &&  (Pip4_orig + Kp4_orig).M() < 0.95 ) { //  for right pi-K masses assignment

		hmyPsiPKPiMass_fromKstar->Fill( myB0MassAlt ) ;
		if (myB0MassAlt > hmyPsiPKPiMass_fromKstar->GetXaxis()->GetXmin()  &&  myB0MassAlt < hmyPsiPKPiMass_fromKstar->GetXaxis()->GetXmax()) {
		  hmyPsiPKPiMass_fromKstar_ex_noRange->Fill( myB0MassAlt_ex ) ;
		  hmyPsiPKPiMass_fromKstar_diff->Fill(myB0MassAlt - myB0MassAlt_ex) ;
		}
		hpsiPK_vs_KPi->Fill( (Pip4_orig + Kp4_orig).M(), (PsiPp4 + Kp4_orig).M() ) ;
		if ( fabs(myB0MassAlt - 5.280808) < 2*0.00816 ) {
			
		  //hPionFromKstar_pT->Fill( pion_fromB0_Pt );
		  //hKaonFromKstar_pT->Fill( kaon_fromB0_Pt );
		  hPionFromKstar_pT->Fill( Pip4_orig.Pt() );
		  hKaonFromKstar_pT->Fill( Kp4_orig.Pt() );
		  //
		  dEdx_vs_pT_Kstar->Fill( Pip4_orig.Pt(), ((*tr_dedx)[pi_orig_Index]) ) ; 
		  dEdx_vs_pT_Kstar->Fill( Kp4_orig.Pt(), ((*tr_dedx)[ka_orig_Index]) ) ; 
		  dEdx_vs_p_Kstar->Fill( Pip4_orig.P(), ((*tr_dedx)[pi_orig_Index]) ) ; 
		  dEdx_vs_p_Kstar->Fill( Kp4_orig.P(), ((*tr_dedx)[ka_orig_Index]) ) ; 
		  dEdx_vs_pT_piFromKstar->Fill( Pip4_orig.Pt(), ((*tr_dedx)[pi_orig_Index]) ) ; 
		  dEdx_vs_pT_KFromKstar->Fill( Kp4_orig.Pt(), ((*tr_dedx)[ka_orig_Index]) ) ; 
		  dEdx_vs_p_piFromKstar->Fill( Pip4_orig.P(), ((*tr_dedx)[pi_orig_Index]) ) ; 
		  dEdx_vs_p_KFromKstar->Fill( Kp4_orig.P(), ((*tr_dedx)[ka_orig_Index]) ) ; 
		  if (Kp4_orig.P() < 1.) {

		    dEdx_vs_p_piFromKstar_KpCut->Fill( Pip4_orig.P(), ((*tr_dedx)[pi_orig_Index]) ) ; 
		    dEdx_vs_p_KFromKstar_KpCut->Fill( Kp4_orig.P(), ((*tr_dedx)[ka_orig_Index]) ) ; 
		    if ( Kp4_orig.P() < 0.8 ) {
		      dEdxK_pCut->Fill( (*tr_dedx)[ka_orig_Index] ) ;			  
		      dEdxPi_pCut->Fill( (*tr_dedx)[pi_orig_Index] ) ;
		      if ( (*trackShits)[ka_orig_Index] >= 17) {
			dEdxK_pCut_stripCut->Fill( (*tr_dedx)[ka_orig_Index] ) ;			  
			dEdxPi_pCut_stripCut->Fill( (*tr_dedx)[pi_orig_Index] ) ;
		      }
		      if ( Kp4_orig.P() < 0.6) 
			dEdxK_vs_dEdxPi->Fill((*tr_dedx)[pi_orig_Index], (*tr_dedx)[ka_orig_Index]) ;
			    
		    }
		  }
			
		  piFromKstar_stripHits_h->Fill( (*trackShits)[pi_orig_Index] ) ;
		  piFromKstar_stripHits_vs_y_h->Fill( Pip4_orig.Rapidity(), (Float_t)((*trackShits)[pi_orig_Index]) ) ;
		  KFromKstar_stripHits_h->Fill( (*trackShits)[ka_orig_Index] ) ;
		  KFromKstar_stripHits_vs_y_h->Fill( Kp4_orig.Rapidity(), (Float_t)((*trackShits)[ka_orig_Index]) ) ;
		  if ( (*trackShits)[ka_orig_Index] >= 15 ) {
			  
		    dEdx_vs_pT_KhqFromKstar->Fill( Kp4_orig.Pt(), ((*tr_dedx)[ka_orig_Index]) ) ; 
		    dEdx_vs_p_KhqFromKstar->Fill( Kp4_orig.P(), ((*tr_dedx)[ka_orig_Index]) ) ; 
		  }
			
		  Armenteros_Kstar_orig->Fill(alpha_orig, pi_KstarPt_orig) ;
		  Armenteros_Kstar->Fill(alpha, pi_KstarPt) ;
		  pT_diff->Fill(pi_KstarPt_orig - K_KstarPt_orig) ;
		  alpha_vs_alphaEx->Fill(alpha_orig, alpha_orig_ex) ;
		} // B0 2*sigma mass window
	      }// right pi-K masses assignment
	      else if ( (Pip4_orig + Kp4_orig).M() < 0.8  ||  (Pip4_orig + Kp4_orig).M() > 1. ) { // sideband for right pi-K masses assignment
		
		hmyPsiPKPiMass_fromKstar_sb->Fill( myB0MassAlt ) ;
		hpsiPK_vs_KPi_sb->Fill( (Pip4_orig + Kp4_orig).M(), (PsiPp4 + Kp4_orig).M()) ;
		if ( fabs(myB0MassAlt - 5.28057) < 2*0.01227 ) {
			
		  //hPionFromKstar_sb_pT->Fill( pion_fromB0_Pt );
		  //hKaonFromKstar_sb_pT->Fill( kaon_fromB0_Pt );
		  hPionFromKstar_sb_pT->Fill( Pip4_orig.Pt() );
		  hKaonFromKstar_sb_pT->Fill( Kp4_orig.Pt() );
		  //
		  Armenteros_Kstar_orig_sb->Fill(alpha_orig,pi_KstarPt_orig) ;
		  pT_diff_sb->Fill(pi_KstarPt_orig - K_KstarPt_orig) ;
		}		    
	      }
	      //
	      // for exchanged pi-K masses
	      if ( (Pip4_orig_ex + Kp4_orig_ex).M() > 0.85  &&  (Pip4_orig_ex + Kp4_orig_ex).M() < 0.95 ) {
		      
		hmyPsiPKPiMass_fromKstar_ex->Fill( myB0MassAlt_ex ) ;
		hpsiPK_vs_KPi_ex->Fill( (Pip4_orig_ex + Kp4_orig_ex).M(), (PsiPp4 + Kp4_orig_ex).M()) ;
		if ( fabs(myB0MassAlt_ex - 5.277776) < 2*0.01161 ) {
		    
		  //hPionFromKstar_ex_pT->Fill( pion_fromB0_Pt );
		  //hKaonFromKstar_ex_pT->Fill( kaon_fromB0_Pt );
		  hPionFromKstar_ex_pT->Fill( Pip4_orig_ex.Pt() );
		  hKaonFromKstar_ex_pT->Fill( Kp4_orig_ex.Pt() );
		  //
		  Armenteros_Kstar_orig_ex->Fill(alpha_orig_ex,pi_ex_KstarPt) ;
		  pT_diff_ex->Fill(pi_ex_KstarPt - K_ex_KstarPt) ;
		}
	      }// for exchanged pi-K masses
	    //} // B0 mass window (5.25,5.3)
	    } // if ( baselineB0Sel ) {
	  } // cuts
	//} // nB0 loop // closed below
      } // if ( HLT_5_vAll && nB0>0 )
      
    /*
    //////////////////////////////////////// my scan
    if ( 
      1
      //         && (*trackPhits)[mypiIdx]>=2 && (*trackShits)[mypiIdx]>=7
      //         && (*muPhits)[mymupidx]>=2 && (*muShits)[mymupidx]>=7
      //         && (*muPhits)[mymumidx]>=2 && (*muShits)[mymumidx]>=7
      //         &&  (((*muQual)[ mymupidx ]) & (1 << muonQual[0] )) &&  (((*muQual)[ mymumidx ]) & (1 << muonQual[0] ))
      //
      // look in the Analyzer for def of goodTightMuon since some cuts are commented out
      //&& (*goodTightMuon)[mymupidx] && (*goodTightMuon)[mymumidx]
      //
      
      //// my soft muons
      //&& ( ((*muQual)[ mymupidx ]) & (1 << muonQual[3]) )  &&  ( ((*muQual)[ mymumidx ]) & (1 << muonQual[3]) )
      ////&& trackerLayersWithMeasurement() > 5 // absent variable
      //&& (*muPhits)[mymupidx] > 1  &&  (*muPhits)[mymumidx] > 1
      //&& (*muChi2)[mymupidx]/(*muNDF)[mymupidx] < 1.8  &&  (*muChi2)[mymumidx]/(*muNDF)[mymumidx] < 1.8
      //&& fabs((*muDxyVtx)[mymupidx]) < 3  &&  fabs((*muDxyVtx)[mymumidx]) < 3
      //&& fabs((*muDzVtx)[mymupidx]) < 30  &&  fabs((*muDzVtx)[mymumidx]) < 30
      
      //
      // 2012 soft muons
      //&& (*goodLooseMuon)[mymupidx] && (*goodLooseMuon)[mymumidx] // not implemented in our rootuples
      && ( ((*muQual)[ mymupidx ]) & (1 << muonQual[3]) )  &&  ( ((*muQual)[ mymumidx ]) & (1 << muonQual[3]) )
      && (*muChi2)[mymupidx]/(*muNDF)[mymupidx] < 1.8  &&  (*muChi2)[mymumidx]/(*muNDF)[mymumidx] < 1.8
      // new soft muons
      && isNewSoftMuons
      //
      //&& sqrt ( pow((*muPx)[mymupidx],2) + pow((*muPy)[mymupidx],2) ) > 7.0 // mu+
      //&& sqrt ( pow((*muPx)[mymumidx],2) + pow((*muPy)[mymumidx],2) ) > 7.0 // mu-
      //&& (*MuMuPiVtx_CL)[myB0Idx] > 0.1
      //&& sqrt( pow((*PionPx_MuMuPi)[myB0Idx],2) + pow((*PionPy_MuMuPi)[myB0Idx],2) ) > 7.0
      //&& ((*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] ) < 3
      //
      //         && ((*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] ) >3.0
      //
      //         && (*ZCTauBS)[myB0Idx] <0.01
      //         && (*ZCTauBS)[myB0Idx] >=0.006
      //
      // && (*ZCTauBSE)[myB0Idx] < 0.01
      //
      //         && sqrt( pow((*MuMuPx)[(*JPsiIndex)[myB0Idx]], 2) + pow((*MuMuPy)[(*JPsiIndex)[myB0Idx]], 2) ) > 7.5 
      //
      // && fabs(Jp4.Eta())<1.5 && fabs(Pip4.Eta())<1.5
      //         && Jp4.DeltaR(Pip4)<0.7
      //         && Z_Pt > 10.0
      
      // my cuts
      && (*MuMuVtx_CL)[myjpsiindex] > 0.005 // Lucia: 0.01  // not needed with one_psi2S already required
      //&& (*muPhits)[mymupidx] >= 1 && (*muPhits)[mymumidx] >= 1  // moved above
      //&& (*muShits)[mymupidx] >= 8 && (*muShits)[mymumidx] >= 8  // moved above
      //&& (*muChi2)[mymupidx] < 7 && (*muChi2)[mymumidx] < 7 // to check  // moved above
      && fabs(Jp4.Rapidity()) < 2.4 // not needed with one_psi2S already required
      )
      {  
      // sigma fit parameters
      
      //Float_t mean[] = {3.6757, 3.67518, 3.67583, 3.67592, 3.67805, 3.67982, 3.68091, 3.68223, 3.68185, 3.68189, 3.68222, 3.68191, 3.68170, 3.68166, 3.68195, 3.68057, 3.67929, 3.67864, 3.67586, 3.67667, 3.67571, 3.6757} ; // 22 elements
      ////Float_t p[] = {0.02281, 0.05011, 1.637, -0.005282} ;
      ////Float_t sigma = p[0] + (p[1]-p[0]) * (1 - TMath::Cos(p[2]*Jp4.Rapidity() - p[3]))/2 ; // mean fit is not good 
      //Float_t sigma[] = {0.0584, 0.05042, 0.04848, 0.04695, 0.04376, 0.04010, 0.03562, 0.03003, 0.02712, 0.02454, 0.02288, 0.02304, 0.02486, 0.02698, 0.03027, 0.03554, 0.04006, 0.04485, 0.04672, 0.04754, 0.05142, 0.0524} ; // 22 elements
      
      
      //// with newSoftMuonsID
      //Float_t mean[] = {3.6742, 3.67536, 3.67607, 3.67600, 3.67811, 3.67958, 3.68074, 3.68218, 3.68180, 3.68190, 3.68224, 3.68197, 3.68180, 3.68166, 3.68188, 3.68067, 3.67915, 3.67836, 3.67590, 3.67692, 3.67597, 3.6751} ; // 22 elements
      //Float_t sigma[] = {0.0578, 0.05019, 0.04817, 0.04648, 0.04359, 0.03976, 0.03516, 0.02963, 0.02679, 0.02418, 0.02259, 0.02274, 0.02452, 0.02658, 0.02987, 0.03506, 0.03965, 0.04443, 0.04667, 0.04723, 0.05099, 0.0512} ; // 22 elements
      
      //Float_t fit_mass=0, fit_sigma=0 ;
      ////Float_t mass_low=0, mass_high=0 ;
      //// signal and sidebands definition
      //Float_t signal_win=2, sb_low=4.5, sb_high=6.5 ;
      
      //for (Int_t i = 0+1; i < rapid_bins-1; i++) { // excluding first and last bin
      //  if ( Jp4.Rapidity() > rapid_min + rapid_binSize*i &&
      //	 Jp4.Rapidity() < rapid_min + rapid_binSize*(i+1) ) { 
      //    //mass_low = mean[i-1] - signal_win*sigma; mass_high = mean[i-1] + signal_win*sigma;
      //    fit_mass = mean[i-1] ; fit_sigma = sigma[i-1] ;
      //  }
      //}
      
      //if (Jp4.Pt() > 5) // not needed with one_psi2S already required
      if ( ( HLT_5_vAll  &&  Jp4.Pt() > 5 )
      || ( HLT_2011  &&  Jp4.Pt() > 7 ) 
      || ( HLT_3p5_vAll  &&  Jp4.Pt() > 3.5 ) ) // not needed with one_psi2S already required
	    
      { 
      if ( HLT_5_vAll ) { // not needed with onePsiPrime
      if ( fabs(mu1.Rapidity()) < 2.4  &&  fabs(mu2.Rapidity()) < 2.4 )
      psi2S_HLT5_yMu_h->Fill( Jp4.M() ) ;
      }
      
      if ( 1  // baseline cuts
      //&& fabs( Jp4.M() - psi2S_mass ) < 0.06 
      //&& (HLT_5_vAll) // not needed with onePsiPrime
      //      &&(*ZCosAlphaBS)[myB0Idx]>0.95
      // my cuts
      && (*MuMuPiVtx_CL)[myB0Idx] > 0.01 // Lucia: 0.01
      && (*ZCosAlphaPV)[myB0Idx] > 0.9925
      && fabs((*ZCosAlphaPV)[myB0Idx]) > cosAlpha_cut // = 26° = |Lucia| = 0.9
      && (*ZCosAlphaPV)[myB0Idx] > 0.99985 // w/o muID: 0.9993 optimized
      && fabs(Pip4.Rapidity()) < 2.4
      && Pip4.Pt() > 0.7 //1
      && Jp4.DeltaR(Pip4) < 1.0
      && (*trackPhits)[mypiIdx] >= 2 && (*trackShits)[mypiIdx] >= 8 
      ) 
      {
      myMuMuMass->Fill( Jp4.M() );
      // pi+ variables
      pi_Chi2N_test_h->Fill( (*trackChi2)[mypiIdx] / (*trackNDF)[mypiIdx] ) ;
      pi_pT_test_h->Fill( Pip4.Pt() ) ;
      pi_y_test_h->Fill( Pip4.Rapidity() ) ;
      pT_vs_pixelHits_pi_h->Fill( (Float_t)((*trackPhits)[mypiIdx]), Pip4.Pt() ) ;
      pT_vs_stripHits_pi_h->Fill( (Float_t)((*trackShits)[mypiIdx]), Pip4.Pt() ) ;
      pT_vs_y_pi_h->Fill( Pip4.Rapidity(), Pip4.Pt() ) ;
      pi_Dz_h->Fill((*trackDzVtx)[mypiIdx]) ; pi_Dxy_h->Fill((*trackDxyVtx)[mypiIdx]) ;
      
      Float_t Z_pT_baseline = 0.; // check
      
      //if (!fixedMassWind) Z_pT_baseline = 8. ; // = Lucia 
      //else Z_pT_baseline = 0. ;
      
      if ( Z_Pt > Z_pT_baseline ) { 
      // psi' mass window
      //if ( Jp4.M() > mass_low && Jp4.M() < mass_high ) { // fit mass
      //if ( fabs(Jp4.M() - fit_mass) < signal_win*sigma ) { // fit mass and interpolated sigma
      // moved above
      // if ( fabs(Jp4.M() - fit_mass) < signal_win*fit_sigma ) { // fit mass and fit sigma
      
      psi2S_DeltaR_pi_h[0]->Fill( Jp4.DeltaR(Pip4) ) ; 
      if ((*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] > 3.0)
      psi2S_DeltaR_pi_vs_mumupiMass_nonpr->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass, Jp4.DeltaR(Pip4) ) ;
      //psi2S_pointingAngle_sign->Fill( (*JCosAlpha)[myjpsiindex] ) ;  
      Z_cosPA[0]->Fill( coeff ) ;
      Z_cosPA[0].first->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;  Z_cosPA[0].second->Fill( (*ZCosAlphaPV)[myB0Idx] ) ; 
      Z_pointingAngle[1]->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;
      Z_DeltaR_vs_cosPA[0]->Fill( fabs((*ZCosAlphaPV)[myB0Idx]), Jp4.DeltaR(Pip4) ) ;
      if ( fabs((*ZCosAlphaPV)[myB0Idx]) > cosAlpha_cut2 ) {
      Z_vs_pi_pT_h->Fill( Pip4.Pt(), Z_Pt) ;
      Z_CTauBS_over_Z_CTauBSE_h->Fill( (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] ) ;
      }
      myMuMuPiMass->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ;
      
      
      ///////////////////////// START SCANNING /////////////////////////
      
      for (Int_t k = 0; k < abs(deltaR_piPsi2S.size()); ++k)
        if (Jp4.DeltaR(Pip4) < deltaR_piPsi2S[k]) {
	  for (Int_t j = 0; j < abs(pi_pT.size()); ++j)
	    if (Pip4.Pt() > pi_pT[j]) 
	      myMuMuPiMass_piPt[0][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=0 -> all
	    for (Int_t j = 0; j < abs(Z_pT.size()); ++j) 
	      if (Z_Pt > Z_pT[j]) 
	        myMuMuPiMass_ZPt[0][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=0 -> all
		  for (Int_t j = 0; j < abs(piAndZ_pT.size()); ++j)
		    if (Pip4.Pt() > piAndZ_pT[j].first  &&  Z_Pt > piAndZ_pT[j].second) 
		      myMuMuPiMass_piAndZPt[0][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=0 -> all ;
	}
		    
      if ( (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] < 1.8 ) { // prompt i=1
        pi_Dz_pr_h->Fill((*trackDzVtx)[mypiIdx]) ; pi_Dxy_pr_h->Fill((*trackDxyVtx)[mypiIdx]) ;
    psi2S_DeltaR_pi_h[1]->Fill( Jp4.DeltaR(Pip4) ) ; psi2S_DeltaR_pi_h[1]->SetLineColor( kBlack ) ;
    Z_cosPA[1]->Fill( coeff ) ; Z_cosPA[1]->SetLineColor( kBlack ) ;
    Z_cosPA[1].first->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;  Z_cosPA[1].second->Fill( (*ZCosAlphaPV)[myB0Idx] ) ; 
    Z_DeltaR_vs_cosPA[1]->Fill( fabs((*ZCosAlphaPV)[myB0Idx]), Jp4.DeltaR(Pip4) ) ;
    if ( fabs((*ZCosAlphaPV)[myB0Idx]) > cosAlpha_cut2 ) { 
    myMuMuPiMass_pr->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; myMuMuPiMass_pr->SetLineColor( kBlack ) ;
    for (Int_t k = 0; k < abs(deltaR_piPsi2S.size()); ++k)
    if (Jp4.DeltaR(Pip4) < deltaR_piPsi2S[k]) {
    for (Int_t j = 0; j < abs(pi_pT.size()); ++j) 
    if (Pip4.Pt() > pi_pT[j]) {
    myMuMuPiMass_piPt[1][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=1 -> prompt
    myMuMuPiMass_piPt[1][j][k]->SetLineColor( kBlack ) ; } // i=1 -> prompt
    for (Int_t j = 0; j < abs(Z_pT.size()); ++j) 
    if (Z_Pt > Z_pT[j]) {
    myMuMuPiMass_ZPt[1][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=1 -> prompt
    myMuMuPiMass_ZPt[1][j][k]->SetLineColor( kBlack ) ; } // i=1 -> prompt
    for (Int_t j = 0; j < abs(piAndZ_pT.size()); ++j)
    if (Pip4.Pt() > piAndZ_pT[j].first  &&  Z_Pt > piAndZ_pT[j].second) {
    myMuMuPiMass_piAndZPt[1][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=1 -> prompt ;
    myMuMuPiMass_piAndZPt[1][j][k]->SetLineColor( kBlack ) ; } // i=1 -> prompt
    }
		      
    if ( fabs((*trackDzVtx)[mypiIdx]) < 1) { // pi Dz < 1
    myMuMuPiMass_pr_Dz01->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; 
    if ( fabs((*trackDzVtx)[mypiIdx]) < 0.4) // pi Dz < 0.4
    myMuMuPiMass_pr_DzCut->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; }
    }

    if ( fabs(Jp4.Rapidity()) < 1 ) { // prompt barrel i=2
    psi2S_DeltaR_pi_h[2]->Fill( Jp4.DeltaR(Pip4) ) ; psi2S_DeltaR_pi_h[2]->SetLineColor( kBlack ) ;
    Z_cosPA[2]->Fill( coeff ) ; Z_cosPA[2]->SetLineColor( kBlack ) ;
    Z_cosPA[2].first->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;  Z_cosPA[2].second->Fill( (*ZCosAlphaPV)[myB0Idx] ) ; 
    Z_DeltaR_vs_cosPA[2]->Fill( fabs((*ZCosAlphaPV)[myB0Idx]), Jp4.DeltaR(Pip4) ) ;
    if ( fabs((*ZCosAlphaPV)[myB0Idx]) > cosAlpha_cut2 ) {
    myMuMuPiMass_barrel_pr->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; myMuMuPiMass_barrel_pr->SetLineColor( kBlack ) ;
    for (Int_t k = 0; k < abs(deltaR_piPsi2S.size()); ++k) 
    if (Jp4.DeltaR(Pip4) < deltaR_piPsi2S[k]) {
    for (Int_t j = 0; j < abs(pi_pT.size()); ++j) 
    if (Pip4.Pt() > pi_pT[j]) {
    myMuMuPiMass_piPt[2][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=2 -> prompt barrel
    myMuMuPiMass_piPt[2][j][k]->SetLineColor( kBlack ) ; } // i=1 -> prompt barrel
    for (Int_t j = 0; j < abs(Z_pT.size()); ++j) 
    if (Z_Pt > Z_pT[j]) {
    myMuMuPiMass_ZPt[2][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=2 -> prompt barrel
    myMuMuPiMass_ZPt[2][j][k]->SetLineColor( kBlack ) ; } // i=2 -> prompt barrel
    for (Int_t j = 0; j < abs(piAndZ_pT.size()); ++j)
    if (Pip4.Pt() > piAndZ_pT[j].first  &&  Z_Pt > piAndZ_pT[j].second) {
    myMuMuPiMass_piAndZPt[2][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=2 -> prompt barrel;
    myMuMuPiMass_piAndZPt[2][j][k]->SetLineColor( kBlack ) ; } // i=2 -> prompt barrel
    }
    } 
    } else // non barrel
    if ( fabs((*ZCosAlphaPV)[myB0Idx]) > cosAlpha_cut2 ) {
    myMuMuPiMass_nonBarrel_pr->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ;
    myMuMuPiMass_nonBarrel_pr->SetLineColor( kBlack ) ; }
		      
    } else if ( (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] < 3.0) { // quasi prompt i=3
    psi2S_DeltaR_pi_h[3]->Fill( Jp4.DeltaR(Pip4) ) ;
    Z_cosPA[3]->Fill( coeff ) ; 
    Z_cosPA[3].first->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;  Z_cosPA[3].second->Fill( (*ZCosAlphaPV)[myB0Idx] ) ; 
    Z_DeltaR_vs_cosPA[3]->Fill( fabs((*ZCosAlphaPV)[myB0Idx]), Jp4.DeltaR(Pip4) ) ;
    if ( (*ZCosAlphaPV)[myB0Idx] > cosAlpha_cut2 ) {
    myMuMuPiMass_quasipr->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ;
    for (Int_t k = 0; k < abs(deltaR_piPsi2S.size()); ++k) 
    if (Jp4.DeltaR(Pip4) < deltaR_piPsi2S[k]) {
    for (Int_t j = 0; j < abs(pi_pT.size()); ++j) 
    if (Pip4.Pt() > pi_pT[j]) 
    myMuMuPiMass_piPt[3][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=3 -> quasi prompt
    for (Int_t j = 0; j<abs(Z_pT.size()); ++j) 
    if (Z_Pt > Z_pT[j]) 
    myMuMuPiMass_ZPt[3][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=3 -> quasi prompt
    for (Int_t j = 0; j < abs(piAndZ_pT.size()); ++j)
    if (Pip4.Pt() > piAndZ_pT[j].first  &&  Z_Pt > piAndZ_pT[j].second) 
    myMuMuPiMass_piAndZPt[3][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=3 -> quasi prompt ;
    }
    }

    } else { // non prompt i=4
    pi_Dz_nonpr_h->Fill((*trackDzVtx)[mypiIdx]) ; pi_Dxy_nonpr_h->Fill((*trackDxyVtx)[mypiIdx]) ;
    psi2S_DeltaR_pi_h[4]->Fill( Jp4.DeltaR(Pip4) ) ; psi2S_DeltaR_pi_h[4]->SetLineColor( kRed ) ;
    Z_cosPA[4]->Fill( coeff ) ; Z_cosPA[4]->SetLineColor( kRed ) ;
    Z_cosPA[4].first->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;  Z_cosPA[4].second->Fill( (*ZCosAlphaPV)[myB0Idx] ) ; 
    Z_DeltaR_vs_cosPA[4]->Fill( fabs((*ZCosAlphaPV)[myB0Idx]), Jp4.DeltaR(Pip4) ) ;


    if ( 1  // reproducing HLT_DoubleMu3p5_LowMass_Displaced cuts
    &&  fabs(mu1.Eta()) < 2.1 && fabs(mu2.Eta()) < 2.1
    &&  fabs(mu1.Eta() - mu2.Eta()) < 2.2
    &&  (*MuMuVtx_CL)[myjpsiindex] > 0.10
    &&  Jp4.Pt() > 7 
    &&  fabs((*ZCosAlphaPV)[myB0Idx]) > cosAlpha_cut // already asked
    &&  (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] > 3.0 // already in non prompt selection
    ) { 
    Float_t mu_pT_cut = 3.5 ; //mu_pT_cut = 4.0 ; mu_pT_cut = 4.5 ;
    if (HLT_5_vAll) {
    myMuMuKMass_Kai->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass) ;
    if (mu1.Pt() > mu_pT_cut  &&  mu2.Pt() > mu_pT_cut) {
    myMuMuPiMass_HLT5to3p5->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ;
    myMuMuKMass_HLT5to3p5->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ; 
    if ( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass  >  5.24  &&
    (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass  <  5.32 ) {
    mu1_pT_test_h->Fill( mu1.Pt() ) ;  mu2_pT_test_h->Fill( mu2.Pt() ) ;
    mu1_eta_test_h->Fill( mu1.Eta() ) ;  mu2_eta_test_h->Fill( mu2.Eta() ) ; }
    }
    }
    else if (HLT_3p5_vAll) {
    if (mu1.Pt() > mu_pT_cut  &&  mu2.Pt() > mu_pT_cut) {
    myMuMuPiMass_HLT3p5_allCuts->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ;
    myMuMuKMass_HLT3p5_allCuts->Fill( (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass ) ;
    if ((Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass  >  5.24  &&
    (Jp4 + Kp4).Mag() - Jp4.M() + psi2S_mass  <  5.32 ) {
    mu1_pT_test_h->Fill( mu1.Pt() ) ;  mu2_pT_test_h->Fill( mu2.Pt() ) ;
    mu1_eta_test_h->Fill( mu1.Eta() ) ;  mu2_eta_test_h->Fill( mu2.Eta() ) ; }
    }
    }
    }

    if ( (*ZCosAlphaPV)[myB0Idx] > cosAlpha_cut2 ) {
    myMuMuPiMass_nonpr->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; myMuMuPiMass_nonpr->SetLineColor( kRed ) ;
    if ( (*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] > 5.0  &&  Jp4.DeltaR(Pip4) < 0.2
    &&  Pip4.Pt() > 12  &&  Z_Pt > 40 ) {
    myMuMuPiMass_DeltaR0p2_nonpr5->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; myMuMuPiMass_DeltaR0p2_nonpr5->SetLineColor( kRed ) ; 
    if (Jp4.DeltaR(Pip4) < 0.1) {
    myMuMuPiMass_DeltaR0p1_nonpr5->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; myMuMuPiMass_DeltaR0p1_nonpr5->SetLineColor( kRed ) ; }
    if ((*ZCTauBS)[myB0Idx] / (*ZCTauBSE)[myB0Idx] > 8.0) {
    myMuMuPiMass_DeltaR0p2_nonpr8->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; myMuMuPiMass_DeltaR0p2_nonpr8->SetLineColor( kRed ) ; 
    if (Jp4.DeltaR(Pip4) < 0.1) {
    myMuMuPiMass_DeltaR0p1_nonpr8->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; myMuMuPiMass_DeltaR0p1_nonpr8->SetLineColor( kRed ) ; }
    }
    }
    for (Int_t k = 0; k < abs(deltaR_piPsi2S.size()); ++k)
    if (Jp4.DeltaR(Pip4) < deltaR_piPsi2S[k]) {
    for (Int_t j = 0; j < abs(pi_pT.size()); ++j) 
    if (Pip4.Pt() > pi_pT[j]) {
    myMuMuPiMass_piPt[4][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=4 -> non prompt
    myMuMuPiMass_piPt[4][j][k]->SetLineColor( kRed ) ; } // i=4 -> non prompt
    for (Int_t j = 0; j<abs(Z_pT.size()); ++j) 
    if (Z_Pt > Z_pT[j]) {
    myMuMuPiMass_ZPt[4][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=4 -> non prompt
    myMuMuPiMass_ZPt[4][j][k]->SetLineColor( kRed ) ; } // i=4 -> non prompt
    for (Int_t j = 0; j < abs(piAndZ_pT.size()); ++j)
    if (Pip4.Pt() > piAndZ_pT[j].first  &&  Z_Pt > piAndZ_pT[j].second) {
    myMuMuPiMass_piAndZPt[4][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // 4=1 -> non prompt ;
    myMuMuPiMass_piAndZPt[4][j][k]->SetLineColor( kRed ) ; } // i=4 -> non prompt
    }
    }
		      
    if ( fabs(Jp4.Rapidity()) < 1 ) { // non prompt barrel i=5
    psi2S_DeltaR_pi_h[5]->Fill( Jp4.DeltaR(Pip4) ) ; psi2S_DeltaR_pi_h[5]->SetLineColor( kRed ) ;
    Z_cosPA[5]->Fill( coeff ) ; Z_cosPA[5]->SetLineColor( kRed ) ;
    Z_cosPA[5].first->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;  Z_cosPA[5].second->Fill( (*ZCosAlphaPV)[myB0Idx] ) ; 
    Z_DeltaR_vs_cosPA[5]->Fill( fabs((*ZCosAlphaPV)[myB0Idx]), Jp4.DeltaR(Pip4) ) ;
    if ( (*ZCosAlphaPV)[myB0Idx] > cosAlpha_cut2 ) {
    myMuMuPiMass_barrel_nonpr->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ;
    for (Int_t k = 0; k < abs(deltaR_piPsi2S.size()); ++k) 
    if (Jp4.DeltaR(Pip4) < deltaR_piPsi2S[k]) {
    for (Int_t j = 0; j < abs(pi_pT.size()); ++j) 
    if (Pip4.Pt() > pi_pT[j]) {
    myMuMuPiMass_piPt[5][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=5 -> non prompt barrel
    myMuMuPiMass_piPt[5][j][k]->SetLineColor( kRed ) ; } // i=5 -> non prompt barrel
    for (Int_t j = 0; j < abs(Z_pT.size()); ++j) 
    if (Z_Pt > Z_pT[j]) {
    myMuMuPiMass_ZPt[5][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=5 -> non prompt barrel
    myMuMuPiMass_ZPt[5][j][k]->SetLineColor( kRed ) ; } // i=5 -> non prompt barrel
    for (Int_t j = 0; j < abs(piAndZ_pT.size()); ++j)
    if (Pip4.Pt() > piAndZ_pT[j].first  &&  Z_Pt > piAndZ_pT[j].second) { 
    myMuMuPiMass_piAndZPt[5][j][k]->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; // i=5 -> non prompt barrel
    myMuMuPiMass_piAndZPt[5][j][k]->SetLineColor( kRed ) ; } // i=5 -> non prompt barrel
    }
						
    if ( fabs((*trackDzVtx)[mypiIdx]) > 1) // pi Dz > 1
    myMuMuPiMass_nonpr_Dz1inf->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ; 
    }
    }
    }
    if ( fabs(Jp4.Rapidity()) < 1 ) // barrel
    myMuMuPiMass_barrel->Fill( (*MuMuPiMass)[myB0Idx] - Jp4.M() + psi2S_mass ) ;
    //} else if ( fabs(Jp4.M() - fit_mass) > sb_low*sigma && fabs(Jp4.M() - fit_mass) < sb_high*sigma ) { // fit mass and interpolated sigma
    // } else if ( fabs(Jp4.M() - fit_mass) > sb_low*fit_sigma && fabs(Jp4.M() - fit_mass) < sb_high*fit_sigma ) { // fit mass and fit sigma
    ////psi2S_pointingAngle_sb->Fill( (*JCosAlpha)[myjpsiindex] ) ; 
    //Z_pointingAngle[2]->Fill( (*ZCosAlphaPV)[myB0Idx] ) ; // moved below
    //}  // psi' signal_window
    } 
    } else Z_pointingAngle[2]->Fill( (*ZCosAlphaPV)[myB0Idx] ) ;
    } 
    //
    //cout<<(*B_cMass)[myB0Idx]   <<endl;
    //
    }*/
  } // end loop over Z/B0

  nB0AC_noMassWin_h->Fill( nB0AC_noMassWin ) ;
  nB0AC_h->Fill( nB0AC ) ;  
  nB0AC_signalWin_h->Fill( nB0AC_signalWin ) ;
  nB0AC_purityWin_h->Fill( nB0AC_purityWin ) ;

  Float_t minDeltaRB0 = 0.1, minDeltaPtB0 = 0.02 ; // (0.1,0.02)
  if (nB0AC == 1) {
    hB0Mass_1B0->Fill( B0Mass_2B0[0][0] ) ;
    hB0CTau_1B0->Fill( B0CTau_2B0[0][0] ) ;
    hB0Mass_noTwins_noSignalWinNotTwins->Fill( B0Mass_2B0[0][0] ) ;
    if ( fabs(B0Mass_2B0[0][0] - B0_massFit) < 3*B0_sigmaFit ) { 
      psi2SPiSqMass_vs_KPiSqMass_1B0[0]->Fill( (Kp4_2B0[0][0] + Pip4_2B0[0][0]).M2(), (PsiPp4_2B0[0][0] + Pip4_2B0[0][0]).M2() ) ;
      psi2S_helicityAngle->Fill( theta_psi ) ;
      //cos_psi2S_helicityAngle->Fill( TMath::Cos(theta_psi) ) ;
      cos_psi2S_helicityAngle->Fill( cosTheta_psi ) ;
      cos_psi2S_helicityAngle_rot->Fill( cosTheta_psi_rot ) ;
      cos_psi2S_helicityAngle_rotVSnotRot->Fill( cosTheta_psi, cosTheta_psi_rot ) ;
      //
      planesAngle->Fill( phi ) ;
    }

    if ( MC ) {
      if ( (DeltaR_1B0[0] < minDeltaRB0) && (DeltaPt_1B0[0] < minDeltaPtB0) ) {
	hB0Mass_1B0matched[0][0]->Fill( B0Mass_2B0[0][0] ) ;
	hB0CTau_1B0matched->Fill( B0CTau_2B0[0][0] ) ; 

	hMCTruthMatching_1B0[0][0]->Fill(1) ;
	//if ( (chargePi_1B0[0] > 0) && (chargeK_1B0[0] > 0) ) {
	if ( MuMuMatched_1B0[0] && pionMatched_1B0[0] && kaonMatched_1B0[0] ) {
	  hMCTruthMatching_1B0[0][1]->Fill(1) ;
	  hB0Mass_1B0matched[0][1]->Fill( B0Mass_2B0[0][0] ) ;
	  
	  if ( fabs(B0Mass_2B0[0][0] - B0_massFit) < 3*B0_sigmaFit ) {
	    // primary vertex comparison
	    priVtx_z[0]->Fill( genEvtVtx_XYZ[2] ) ;
	    for (Int_t j=0; j<2; j++) 
	      for (Int_t i=0; i<3; i++) 
		priVtx_delta[j][i]->Fill( genEvtVtx_XYZ[i] - priVtxGen_XYZ[i] ) ;

	    if ( fabs(genEvtVtx_XYZ[2] - priVtxGen_XYZ[2]) < 0.0001) {
	      if ( fabs(priVtx_XYZ[2] - genEvtVtx_XYZ[2]) < 0.1 ) 
		if ( fabs(priVtx_XYZ[2] - priVtxB0Less_XYZ[2]) > 0.1 ) {
		  cout <<"\nB0CTau_2B0[0][0] = " <<B0CTau_2B0[0][0] <<" ; " ;
		  for (Int_t i=0; i<3; i++) {
		    cout <<"\ngenEvtVtx_XYZ[" <<i <<"] = " <<genEvtVtx_XYZ[i] <<" ; " ;
		    cout <<"pB0_XYZ[" <<i <<"] = " <<pB0_XYZ[i] <<" ; " ;
		    cout <<"secVtx_XYZ[" <<i <<"] = " <<secVtx_XYZ[i] <<" ; " ;
		    cout <<"priVtx_XYZ[" <<i <<"] = " <<priVtx_XYZ[i] <<" ; " ;
		    cout <<"priVtx_B0CosAlpha2D_XYZ[" <<i <<"] = " <<priVtx_B0CosAlpha_XYZ[0][i] <<" ; " ;
		    cout <<"priVtx_B0CosAlpha3D_XYZ[" <<i <<"] = " <<priVtx_B0CosAlpha_XYZ[1][i] <<" ; " ;
		    cout <<"priVtxB0Less_XYZ[" <<i <<"] = " <<priVtxB0Less_XYZ[i] <<" ; " ;
		    cout <<"priVtxB0Less_B0CosAlpha2D_XYZ[" <<i <<"] = " <<priVtxB0Less_B0CosAlpha_XYZ[0][i] <<" ; " ;
		    cout <<"priVtxB0Less_B0CosAlpha3D_XYZ[" <<i <<"] = " <<priVtxB0Less_B0CosAlpha_XYZ[1][i] <<" ;\n" ;
		  }
		  cout <<endl ;
		}
	      priVtx_deltaX_test->Fill( genEvtVtx_XYZ[0] - priVtxGen_XYZ[0] ) ;
	      priVtx_deltaY_test->Fill( genEvtVtx_XYZ[1] - priVtxGen_XYZ[1] ) ;
	      //
	      // B0 cos(alpha)
	      B0CosAlpha_2Dvs3D[0]->Fill(priVtx_XYZ[3], priVtx_XYZ[4]) ;
	      B0CosAlpha_2Dvs3D[1]->Fill(priVtx_B0CosAlpha_XYZ[0][3], priVtx_B0CosAlpha_XYZ[1][4]) ;
	      //B0CosAlpha_2Dvs3D[2]->Fill(priVtxB0Less_XYZ[3], priVtxB0Less_XYZ[4]) ;
	      B0CosAlpha_2Dvs3D[3]->Fill(priVtxB0Less_B0CosAlpha_XYZ[0][3], priVtxB0Less_B0CosAlpha_XYZ[1][4]) ;
	      //
	      priVtx_z[1]->Fill( priVtx_XYZ[2] ) ;
	      priVtx_z[2]->Fill( priVtx_B0CosAlpha_XYZ[0][2] ) ;
	      priVtx_z[3]->Fill( priVtx_B0CosAlpha_XYZ[1][2] ) ;
	      priVtx_z[4]->Fill( priVtxB0Less_XYZ[2] ) ;
	      priVtx_z[5]->Fill( priVtxB0Less_B0CosAlpha_XYZ[0][2] ) ;
	      priVtx_z[6]->Fill( priVtxB0Less_B0CosAlpha_XYZ[1][2] ) ;
	      //
	      for (Int_t j=0; j<2; j++) {
		for (Int_t i=0; i<3; i++) {
		  priVtx_B0_delta[j][i]->Fill( genEvtVtx_XYZ[i] - priVtx_XYZ[i] ) ;
		  for (Int_t k=0; k<2; k++)
		    priVtx_B0CosAlpha_delta[k][j][i]->Fill( genEvtVtx_XYZ[i] - priVtx_B0CosAlpha_XYZ[k][i] );
		  priVtxB0Less_delta[j][i]->Fill( genEvtVtx_XYZ[i] - priVtxB0Less_XYZ[i] ) ;
		  for (Int_t k=0; k<2; k++)
		    priVtxB0Less_B0CosAlpha_delta[k][j][i]->Fill( genEvtVtx_XYZ[i] - priVtxB0Less_B0CosAlpha_XYZ[k][i] ) ;
		}   
		dz_pri_priB0Less[j]->Fill( priVtx_XYZ[2] - priVtxB0Less_XYZ[2] ) ;
		for (Int_t k=0; k<2; k++) {
		  dz_pri_priB0CosAlpha[k][j]->Fill( priVtx_XYZ[2] - priVtx_B0CosAlpha_XYZ[k][2] ) ;
		  dz_pri_priB0LessB0CosAlpha[k][j]->Fill( priVtx_XYZ[2] - priVtxB0Less_B0CosAlpha_XYZ[k][2] ) ;
		  dz_priB0CosAlpha_priB0LessB0CosAlpha[k][j]->Fill( priVtx_B0CosAlpha_XYZ[k][2] - priVtxB0Less_B0CosAlpha_XYZ[k][2] ) ;
		  dz_priB0Less_priB0LessB0CosAlpha[k][j]->Fill( priVtxB0Less_XYZ[2] - priVtxB0Less_B0CosAlpha_XYZ[k][2] ) ;
		}
	      }
	    }
	  }

	} else
	  hMCTruthMatching_1B0[0][1]->Fill(0) ;
      } else {
	hMCTruthMatching_1B0[0][0]->Fill(0) ;
	hMCTruthMatching_1B0[0][1]->Fill(0) ;
      }
      //
      if ( MuMuMatched_1B0[0] && pionMatched_1B0[0] && kaonMatched_1B0[0] )
	hMCTruthMatching_1B0[0][2]->Fill(1) ;
      else
	hMCTruthMatching_1B0[0][2]->Fill(0) ;
    }
  } else if (nB0AC == 2) {
    // twins histograms
    nTwins_h[0]->Fill( nTwins ) ;
    // no twin
    if (nTwins == 0) { 
      hB0Mass_2B0notTwin->Fill( B0Mass_2B0[0][0] ) ; hB0Mass_2B0notTwin->Fill( B0Mass_2B0[0][1] ) ; 
      if ( MC ) {
	if ( (DeltaR_2B0[0][0] < minDeltaRB0 && DeltaPt_2B0[0][0] < minDeltaPtB0) && (DeltaR_2B0[0][1] < minDeltaRB0 && DeltaPt_2B0[0][1] < minDeltaPtB0) ) {
	  hMCTruthMatching_2B0notTwin[0]->Fill( 2 ) ;
	  for (Int_t i=0; i<nTwins+2; i++) {
	    hB0Mass_2B0notTwin_matched[0]->Fill( B0Mass_2B0[0][i] ) ; 
	    hMCDeltaRPi_2B0->Fill( DrPi_2B0[i] );
	    hMCDeltaRK_2B0->Fill( DrK_2B0[i] );
	    //if ( (chargePi_2B0[0][i] > 0) && (chargeK_2B0[0][i] > 0) ) {
	    if ( MuMuMatched[0][i] && pionMatched[0][i] && kaonMatched[0][i] ) {
	      hMCTruthMatching_2B0notTwin[1]->Fill( 2 ) ;
	      hB0Mass_2B0notTwin_matched[1]->Fill( B0Mass_2B0[0][i] ) ; 
	    } else
	      hMCTruthMatching_2B0notTwin[1]->Fill( 0 ) ;
	  }	 
	} else if ( (DeltaR_2B0[0][0] < minDeltaRB0 && DeltaPt_2B0[0][0] < minDeltaPtB0) || (DeltaR_2B0[0][1] < minDeltaRB0 && DeltaPt_2B0[0][1] < minDeltaPtB0) ) {
	  hMCTruthMatching_2B0notTwin[0]->Fill( 1 ) ;
	  for (Int_t i=0; i<nTwins+2; i++)
	    if ( (DeltaR_2B0[0][i] < minDeltaRB0) && (DeltaPt_2B0[0][i] < minDeltaPtB0) ) {
	      hB0Mass_2B0notTwin_matched[0]->Fill( B0Mass_2B0[0][i] ) ;
	      hMCDeltaRPi_2B0->Fill( DrPi_2B0[i] );
	      hMCDeltaRK_2B0->Fill( DrK_2B0[i] );
	      //if ( (chargePi_2B0[0][i] > 0) && (chargeK_2B0[0][i] > 0) ) {
	      if ( MuMuMatched[0][i] && pionMatched[0][i] && kaonMatched[0][i] ) {
		hMCTruthMatching_2B0notTwin[1]->Fill( 1 ) ;
		hB0Mass_2B0notTwin_matched[1]->Fill( B0Mass_2B0[0][i] ) ; 
	      } else
		hMCTruthMatching_2B0notTwin[1]->Fill( 0 ) ;
	    }  
	} else {
	  hMCTruthMatching_2B0notTwin[0]->Fill( 0 ) ;
	  hMCTruthMatching_2B0notTwin[1]->Fill( 0 ) ;
	}
	//
	if ( (MuMuMatched[0][0] && pionMatched[0][0] && kaonMatched[0][0]) && (MuMuMatched[0][1] && pionMatched[0][1] && kaonMatched[0][1]) )
	  hMCTruthMatching_2B0notTwin[2]->Fill( 2 ) ;
	else if ( (MuMuMatched[0][0] && pionMatched[0][0] && kaonMatched[0][0]) || (MuMuMatched[0][1] && pionMatched[0][1] && kaonMatched[0][1]) )
	  hMCTruthMatching_2B0notTwin[2]->Fill( 1 ) ;
	else
	  hMCTruthMatching_2B0notTwin[2]->Fill( 0 ) ;

	for (Int_t i=0; i<nTwins+2; i++) {
	  if ( MuMuMatched[0][i] && pionMatched[0][i] && kaonMatched[0][i] ) {
	    hB0Mass_2B0notTwin_matched[2]->Fill( B0Mass_2B0[0][i] ) ;
	    if (i==0) 
	      hDeltaB0Mass_2B0notTwin[0]->Fill( B0Mass_2B0[0][i] - B0Mass_gen, B0Mass_2B0[0][1] - B0Mass_gen) ;
	    else
	      hDeltaB0Mass_2B0notTwin[0]->Fill( B0Mass_2B0[0][i] - B0Mass_gen, B0Mass_2B0[0][0] - B0Mass_gen) ;
	  }
	}
      }
      //
      if ( (fabs(B0Mass_2B0[0][0] - B0_massFit) > 3*B0_sigmaFit)  ||  (fabs(B0Mass_2B0[0][1] - B0_massFit) > 3*B0_sigmaFit) ) {
	hB0Mass_noTwins_noSignalWinNotTwins->Fill( B0Mass_2B0[0][0] ) ;
	hB0Mass_noTwins_noSignalWinNotTwins->Fill( B0Mass_2B0[0][1] ) ;
      }
    } else { ////////// twin //////////
      hB0Mass_2B0twin->Fill( B0Mass_2B0[0][0] ) ; hB0Mass_2B0twin->Fill( B0Mass_2B0[0][1] ) ;
      // MC matching
      if ( MC ) {
	if ( (DeltaR_2B0[0][0] < minDeltaRB0 && DeltaPt_2B0[0][0] < minDeltaPtB0) && (DeltaR_2B0[0][1] < minDeltaRB0 && DeltaPt_2B0[0][1] < minDeltaPtB0) ) {
	  hMCTruthMatching_2B0twin[0]->Fill( 2 ) ;
	  for (Int_t i=0; i<nTwins+1; i++) {
	    hB0Mass_2B0twin_matched[0]->Fill( B0Mass_2B0[0][i] ) ;
	    //if ( (chargePi_2B0[0][i] > 0) && (chargeK_2B0[0][i] > 0) ) {
	    if ( MuMuMatched[0][i] && pionMatched[0][i] && kaonMatched[0][i] ) {
	      hMCTruthMatching_2B0twin[1]->Fill( 2 ) ;
	      hB0Mass_2B0twin_matched[1]->Fill( B0Mass_2B0[0][i] ) ;
	      //
	      if (i==0) {
		hDeltaB0Mass_2B0twin[0]->Fill( B0Mass_2B0[0][i] - B0Mass_gen, B0Mass_2B0[0][1] - B0Mass_gen ) ;
		hAlpha_2B0twin[0]->Fill( alpha_2B0[i], alpha_2B0[1] ) ;
		if ( fabs(B0Mass_2B0[0][i] - B0Mass_gen) < fabs(B0Mass_2B0[0][1] - B0Mass_gen) ) { // inside
		  hB0Mass_2B0twin_hourglass[1][1]->Fill( B0Mass_2B0[0][i] ) ; // matched
		  hB0Mass_2B0twin_hourglass[1][0]->Fill( B0Mass_2B0[0][1] ) ; // unmatched
		} else { // outside
		  hB0Mass_2B0twin_hourglass[0][1]->Fill( B0Mass_2B0[0][i] ) ; // matched
		  hB0Mass_2B0twin_hourglass[0][0]->Fill( B0Mass_2B0[0][1] ) ; // unmatched
		}
	      } else {
		hDeltaB0Mass_2B0twin[0]->Fill( B0Mass_2B0[0][i] - B0Mass_gen, B0Mass_2B0[0][0] - B0Mass_gen ) ;
		hAlpha_2B0twin[0]->Fill( alpha_2B0[i], alpha_2B0[0] ) ;	  
		if ( fabs(B0Mass_2B0[0][i] - B0Mass_gen) < fabs(B0Mass_2B0[0][0] - B0Mass_gen) ) { // inside
		  hB0Mass_2B0twin_hourglass[1][1]->Fill( B0Mass_2B0[0][i] ) ; // matched
		  hB0Mass_2B0twin_hourglass[1][0]->Fill( B0Mass_2B0[0][0] ) ; // unmatched
		} else { // outside
		  hB0Mass_2B0twin_hourglass[0][1]->Fill( B0Mass_2B0[0][i] ) ; // matched
		  hB0Mass_2B0twin_hourglass[0][0]->Fill( B0Mass_2B0[0][0] ) ; // unmatched
		}
	      }
	    } else
	      hMCTruthMatching_2B0twin[1]->Fill( 0 ) ;
	  }
	} else if ( (DeltaR_2B0[0][0] < minDeltaRB0 && DeltaPt_2B0[0][0] < minDeltaPtB0) || (DeltaR_2B0[0][1] < minDeltaRB0 && DeltaPt_2B0[0][1] < minDeltaPtB0) ) {
	  hMCTruthMatching_2B0twin[0]->Fill( 1 ) ;
	  for (Int_t i=0; i<nTwins+1; i++) 
	    if ( (DeltaR_2B0[0][i] < minDeltaRB0) && (DeltaPt_2B0[0][i] < minDeltaPtB0) ) {
	      hB0Mass_2B0twin_matched[0]->Fill( B0Mass_2B0[0][i] ) ; 
	      //if ( (chargePi_2B0[0][i] > 0) && (chargeK_2B0[0][i] > 0) ) {
	      if ( MuMuMatched[0][i] && pionMatched[0][i] && kaonMatched[0][i] ) {
		hMCTruthMatching_2B0twin[1]->Fill( 1 ) ;
		hB0Mass_2B0twin_matched[1]->Fill( B0Mass_2B0[0][i] ) ;
		//
		if (i==0) {
		  hDeltaB0Mass_2B0twin[0]->Fill( B0Mass_2B0[0][i] - B0Mass_gen, B0Mass_2B0[0][1] - B0Mass_gen ) ;
		  hAlpha_2B0twin[0]->Fill( alpha_2B0[i], alpha_2B0[1] ) ;
		  if ( fabs(B0Mass_2B0[0][i] - B0Mass_gen) < fabs(B0Mass_2B0[0][1] - B0Mass_gen) ) { // inside
		    hB0Mass_2B0twin_hourglass[1][1]->Fill( B0Mass_2B0[0][i] ) ; // matched
		    hB0Mass_2B0twin_hourglass[1][0]->Fill( B0Mass_2B0[0][1] ) ; // unmatched
		  } else { // outside
		    hB0Mass_2B0twin_hourglass[0][1]->Fill( B0Mass_2B0[0][i] ) ; // matched
		    hB0Mass_2B0twin_hourglass[0][0]->Fill( B0Mass_2B0[0][1] ) ; // unmatched
		  }
		} else {
		  hDeltaB0Mass_2B0twin[0]->Fill( B0Mass_2B0[0][i] - B0Mass_gen, B0Mass_2B0[0][0] - B0Mass_gen ) ;
		  hAlpha_2B0twin[0]->Fill( alpha_2B0[i], alpha_2B0[0] ) ;	  
		  if ( fabs(B0Mass_2B0[0][i] - B0Mass_gen) < fabs(B0Mass_2B0[0][0] - B0Mass_gen) ) { // inside
		    hB0Mass_2B0twin_hourglass[1][1]->Fill( B0Mass_2B0[0][i] ) ; // matched
		    hB0Mass_2B0twin_hourglass[1][0]->Fill( B0Mass_2B0[0][0] ) ; // unmatched
		  } else { // outside
		    hB0Mass_2B0twin_hourglass[0][1]->Fill( B0Mass_2B0[0][i] ) ; // matched
		    hB0Mass_2B0twin_hourglass[0][0]->Fill( B0Mass_2B0[0][0] ) ; // unmatched
		  }
		}
	      } else
		hMCTruthMatching_2B0twin[1]->Fill( 0 ) ;
	    }
	} else {
	  hMCTruthMatching_2B0twin[0]->Fill( 0 ) ;
	  hMCTruthMatching_2B0twin[1]->Fill( 0 ) ;
	}
	//
	if ( (MuMuMatched[0][0] && pionMatched[0][0] && kaonMatched[0][0]) && (MuMuMatched[0][1] && pionMatched[0][1] && kaonMatched[0][1]) )
	  hMCTruthMatching_2B0twin[2]->Fill( 2 ) ;
	else if ( (MuMuMatched[0][0] && pionMatched[0][0] && kaonMatched[0][0]) || (MuMuMatched[0][1] && pionMatched[0][1] && kaonMatched[0][1]) )
	  hMCTruthMatching_2B0twin[2]->Fill( 1 ) ;
	else
	  hMCTruthMatching_2B0twin[2]->Fill( 0 ) ;

      } // MC

    } // twin
  } // if (nB0AC == 2)

  // in signal win
  nB0ACInMC_h->Fill( nB0ACInMC ) ;
  //if ( (!MC && (nB0AC_signalWin == nB0AC)) || (nB0ACInMC == nB0AC) ) { // ?
  if ( (!MC && (nB0AC_signalWin > 0)) || (nB0ACInMC > 0) ) { // nB0ACInMC > 0 only if (MC && nB0AC_signalWin>0) 
    // only 1 B0
    if ( (!MC && (nB0AC_signalWin == 1)) || (nB0ACInMC == 1) ) 
      psi2SPiSqMass_vs_KPiSqMass_1B0[1]->Fill( (Kp4_2B0[1][0] + Pip4_2B0[1][0]).M2(), (PsiPp4_2B0[1][0] + Pip4_2B0[1][0]).M2() ) ;

    if ( nB0ACInMC == 1 ) { // nB0ACInMC > 0 only if (MC && nB0AC_signalWin>0)
      hMCDeltaR_1B0->Fill( DeltaR_1B0[1] ) ;
      hMCDeltaPt_1B0->Fill( DeltaPt_1B0[1] ) ;

      if ( (DeltaR_1B0[1] < minDeltaRB0) && (DeltaPt_1B0[1] < minDeltaPtB0) ) {
	hMCTruthMatching_1B0[1][0]->Fill(1) ;
	hB0Mass_1B0matched[1][0]->Fill( B0Mass_2B0[1][0] ) ;

	//if ( (chargePi_1B0[1] > 0) && (chargeK_1B0[1] > 0) ) { 
	if ( MuMuMatched_1B0[1] && pionMatched_1B0[1] && kaonMatched_1B0[1] ) {
	  hMCTruthMatching_1B0[1][1]->Fill(1) ;
	  hB0Mass_1B0matched[1][1]->Fill( B0Mass_2B0[1][0] ) ;	    
	}
	else
	  hMCTruthMatching_1B0[1][1]->Fill(0) ;
      } else {
	hMCTruthMatching_1B0[1][0]->Fill(0) ;
	hMCTruthMatching_1B0[1][1]->Fill(0) ;
      }
      //
      if ( MuMuMatched_1B0[1] && pionMatched_1B0[1] && kaonMatched_1B0[1] )
	hMCTruthMatching_1B0[1][2]->Fill(1) ;
      else
	hMCTruthMatching_1B0[1][2]->Fill(0) ;
    }
    else if ( (!MC && (nB0AC_signalWin == 2)) || (nB0ACInMC == 2) ) { // nB0ACInMC > 0 only if (MC && nB0AC_signalWin>0)
      if ( MC ) {
	hMCDeltaR_2B0[0]->Fill( DeltaR_2B0[1][0] ) ;
	hMCDeltaPt_2B0[0]->Fill( DeltaPt_2B0[1][0] ) ; 
	hMCDeltaR_2B0[1]->Fill( DeltaR_2B0[1][1] ) ; 
	hMCDeltaPt_2B0[1]->Fill( DeltaPt_2B0[1][1] ) ; 
	hMCDeltaR_1vs2_2B0->Fill( DeltaR_2B0[1][0], DeltaR_2B0[1][1] ) ;
	hMCDeltaPt_1vs2_2B0->Fill( DeltaPt_2B0[1][0], DeltaPt_2B0[1][1] ) ;
      }
      nTwins_h[1]->Fill( nTwins_signalWin ) ;
      if (nTwins_signalWin == 0) { // MC matching
	if ( MC ) { 
	  if ( (DeltaR_2B0[1][0] < minDeltaRB0 && DeltaPt_2B0[1][0] < minDeltaPtB0) && (DeltaR_2B0[1][1] < minDeltaRB0 && DeltaPt_2B0[1][1] < minDeltaPtB0) ) 
	    hMCTruthMatching_2B0signalWin_notTwin->Fill( 2 ) ;
	  else if ( (DeltaR_2B0[1][0] < minDeltaRB0 && DeltaPt_2B0[1][0] < minDeltaPtB0) || (DeltaR_2B0[1][1] < minDeltaRB0 && DeltaPt_2B0[1][1] < minDeltaPtB0) ) 
	    hMCTruthMatching_2B0signalWin_notTwin->Fill( 1 ) ;
	  else 
	    hMCTruthMatching_2B0signalWin_notTwin->Fill( 0 ) ;
	  //
	  /*
	    if ( (DeltaR_2B0[1][0] < minDeltaRB0) && (DeltaPt_2B0[1][0] < minDeltaPtB0) )
	    hDeltaB0Mass_2B0notTwin[1]->Fill( B0Mass_2B0[1][0] - B0Mass_gen, B0Mass_2B0[1][1] - B0Mass_gen) ;
	    if ( (DeltaR_2B0[1][1] < minDeltaRB0) && (DeltaPt_2B0[1][1] < minDeltaPtB0) ) 
	    hDeltaB0Mass_2B0notTwin[1]->Fill( B0Mass_2B0[1][1] - B0Mass_gen, B0Mass_2B0[1][0] - B0Mass_gen) ;
	  */
	  if ( MuMuMatched[1][0] && pionMatched[1][0] && kaonMatched[1][0] )
	    hDeltaB0Mass_2B0notTwin[1]->Fill( B0Mass_2B0[1][0] - B0Mass_gen, B0Mass_2B0[1][1] - B0Mass_gen) ;
	  if ( MuMuMatched[1][1] && pionMatched[1][1] && kaonMatched[1][1] )
	    hDeltaB0Mass_2B0notTwin[1]->Fill( B0Mass_2B0[1][1] - B0Mass_gen, B0Mass_2B0[1][0] - B0Mass_gen) ;
	}
      } else {
	if ( MC ) { 
	  hMCDeltaR_1vs2_2B0twin->Fill( DeltaR_2B0[1][0], DeltaR_2B0[1][1] ) ;
	  hMCDeltaPt_1vs2_2B0twin->Fill( DeltaPt_2B0[1][0], DeltaPt_2B0[1][1] ) ;
	}
	hB0VtxCL_1vs2_2B0twin->Fill( B0VtxCL_2B0[0], B0VtxCL_2B0[1] ) ;
	hB0Mass_1vs2_2B0twin->Fill( B0Mass_2B0[1][0], B0Mass_2B0[1][1] ) ;
	hB0CTau_1vs2_2B0twin->Fill( B0CTau_2B0[1][0], B0CTau_2B0[1][1] ) ;
	hAlpha_1vs2_2B0twin->Fill( alpha_2B0[0], alpha_2B0[1] ) ;
	
	if ( MC )
	  for (Int_t i=0; i<=nTwins_signalWin; i++)
	    if ( MuMuMatched[1][i] && pionMatched[1][i] && kaonMatched[1][i] ) {
	      hMCDeltaR_matchedTwin[1]->Fill( DeltaR_2B0[1][i] ) ;
	      hMCDeltaPt_matchedTwin[1]->Fill( DeltaPt_2B0[1][i] ) ;
	      hB0VtxCL_matchedTwin[1]->Fill( B0VtxCL_2B0[i] ) ;
	      hB0Mass_matchedTwin[1]->Fill( B0Mass_2B0[1][i] ) ;
	      hB0CTau_matchedTwin[1]->Fill( B0CTau_2B0[1][i] ) ;
	      hAlpha_matchedTwin[1]->Fill( alpha_2B0[i] ) ;
	      if (i==0) {
		hDeltaB0Mass_2B0twin[1]->Fill( B0Mass_2B0[1][i] - B0Mass_gen, B0Mass_2B0[1][1] - B0Mass_gen) ;
		hAlpha_2B0twin[1]->Fill( alpha_2B0[i], alpha_2B0[1] ) ;
	      } else {
		hDeltaB0Mass_2B0twin[1]->Fill( B0Mass_2B0[1][i] - B0Mass_gen, B0Mass_2B0[1][0] - B0Mass_gen) ;
		hAlpha_2B0twin[1]->Fill( alpha_2B0[i], alpha_2B0[0] ) ;	  
	      }
	    } else {
	      hMCDeltaR_matchedTwin[0]->Fill( DeltaR_2B0[1][i] ) ;
	      hMCDeltaPt_matchedTwin[0]->Fill( DeltaPt_2B0[1][i] ) ;
	      hB0VtxCL_matchedTwin[0]->Fill( B0VtxCL_2B0[i] ) ;
	      hB0Mass_matchedTwin[0]->Fill( B0Mass_2B0[1][i] ) ;
	      hB0CTau_matchedTwin[0]->Fill( B0CTau_2B0[1][i] ) ;
	      hAlpha_matchedTwin[0]->Fill( alpha_2B0[i] ) ;
	    }
      } // if (nTwins_signalWin > 0)
    } // if (nB0ACInMC == 2)
  }
  ///////////////////////////////////////////////////

  return kTRUE;
}

float psiPrimePiK_MC::ComputedEdxMass(int itrack, float TrackPt, float K, float C) 
{
  if(C == 0.) C = K;
  //cout << TrackPt << "  "<< (*trackPz)[itrack]<<endl;
  float P = sqrt(TrackPt*TrackPt + (*trackPz)[itrack]*(*trackPz)[itrack]);
  //cout << P << endl;
  float I = (*tr_dedx)[itrack];
  //cout << I << endl;
  
  float result;
  if ( (I-C)>0 )
    result = sqrt((I - C) / K) * P;
  else
    result = 0 ;
  //cout << " Result "<< result  <<endl;
  return result;
}


void psiPrimePiK_MC::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  cout <<"\nSlaveTerminate:" <<endl ;

  TDirectory *savedir = gDirectory;
  if (fOut) 
    {      
      fOut->cd();
      //
      _nt->Write();
      //// Write histograms in output file
      gStyle->SetOptStat(111111) ;
      //
      // hname->Write();
      //
      myMuMuMass->Write() ;  myMuMuMass_noMuID->Write() ; myMuMuMass_MuID->Write() ; myMuMuMass_Mu1ID->Write() ; myMuMuMass_Mu2ID->Write() ; 
      myMuMuMass_HLT->Write() ; myMuMuMass_HLT7->Write() ; myMuMuMass_HLT5_7comp->Write() ; myMuMuMass_HLT3p5->Write() ; //myMuMuMass_test->Write() ;  //myMuMuMass_HLT_yMu->Write() ;  myMuMuMass_HLT_yMu_ypsi->Write() ;
      myMuMuMass_HLT_final->Write() ;  myMuMuMass_HLT_DR1->Write() ;  myMuMuMass_HLT_VtxCLDiff_cut->Write() ;
      myMuMuMass_HLT_max2psiPrime->Write() ;  myMuMuMass_HLT_2psiPrime->Write() ;  myMuMuMass_HLT_1psiPrime->Write() ;
      myMuMuKMass->Write() ;
      myMuMuPiMass_Kai->Write() ; myMuMuPiMass_Kai_HLT7->Write() ; myMuMuKMass_Kai->Write() ;
      myMuMuMass_KaiBF->Write() ; myMuMuKMass_KaiBF->Write() ; dEdx_vs_p_Bplus->Write() ;
      myMuMuMass_noMuID_KaiBF->Write() ; myMuMuMass_MuID_KaiBF->Write() ;
      myMuMuPiMass_HLT5to3p5->Write() ; myMuMuPiMass_HLT3p5_allCuts->Write() ; myMuMuKMass_HLT5to3p5->Write() ; myMuMuKMass_HLT3p5_allCuts->Write() ;
      //
      //
      myInclusiveMuMuMass->Write() ; myInclusiveMuMuMass_noTrigg->Write() ;
      // HLT_0_v*
      myInclusiveMuMuMass_HLT_0_v3->Write() ; myInclusiveMuMuMass_HLT_0_v4->Write() ; myInclusiveMuMuMass_HLT_0_v5->Write() ; myInclusiveMuMuMass_HLT_0_v6->Write() ;
      myInclusiveMuMuMass_HLT_0_v3A->Write() ; myInclusiveMuMuMass_HLT_0_v3B->Write() ;
      myInclusiveMuMuMass_HLT_0_v6C->Write() ; myInclusiveMuMuMass_HLT_0_v6D->Write() ;
      myInclusiveMuMuMass_HLT_0_v3Bv4->Write() ; myInclusiveMuMuMass_HLT_0_v5v6C->Write() ;
      // HLT_5_v*
      myInclusiveMuMuMass_HLT_5_v3->Write() ; myInclusiveMuMuMass_HLT_5_v4->Write() ; myInclusiveMuMuMass_HLT_5_v5->Write() ; myInclusiveMuMuMass_HLT_5_v6->Write() ;
      myInclusiveMuMuMass_HLT_5_v3A->Write() ; myInclusiveMuMuMass_HLT_5_v3B->Write() ;
      myInclusiveMuMuMass_HLT_5_v6C->Write() ; myInclusiveMuMuMass_HLT_5_v6D->Write() ;
      myInclusiveMuMuMass_HLT_5_v3Bv4->Write() ; myInclusiveMuMuMass_HLT_5_v5v6C->Write() ;
      // HLT_7_v*
      myInclusiveMuMuMass_HLT_7_v1->Write() ; myInclusiveMuMuMass_HLT_7_v2->Write() ; myInclusiveMuMuMass_HLT_7_v3->Write() ;
      myInclusiveMuMuMass_HLT_7_v3C->Write() ; myInclusiveMuMuMass_HLT_7_v2v3C->Write() ; myInclusiveMuMuMass_HLT_7_v3D->Write() ;
      // datasets 
      myInclusiveMuMuMass_12A->Write() ;
      myInclusiveMuMuMass_12A_HLT_0_v3->Write() ; myInclusiveMuMuMass_12A_HLT_5_v3->Write() ; myInclusiveMuMuMass_12A_HLT_9_v9->Write() ;
      //
      myInclusiveMuMuMass_12B->Write() ;
      myInclusiveMuMuMass_12B_HLT_0_v3v4->Write() ; myInclusiveMuMuMass_12B_HLT_5_v3v4->Write() ; myInclusiveMuMuMass_12B_HLT_7_v1_9_v9->Write() ;
      //
      myInclusiveMuMuMass_12C->Write() ;
      myInclusiveMuMuMass_12C_HLT_0_v5->Write() ; myInclusiveMuMuMass_12C_HLT_5_v5->Write() ; myInclusiveMuMuMass_12C_HLT_7_v2->Write() ;
      myInclusiveMuMuMass_12C_HLT_5_v5_or_7_v2->Write() ;
      //
      myInclusiveMuMuMass_12D->Write() ;
      myInclusiveMuMuMass_12D_HLT_0_v6->Write() ; myInclusiveMuMuMass_12D_HLT_5_v6->Write() ; myInclusiveMuMuMass_12D_HLT_7_v3->Write() ;
      myInclusiveMuMuMass_out->Write() ;
      
      // Nairit plots
      MuMuKMass_Nairit->Write() ; MuMuKMass_sb->Write() ;
      trk_pT_test_h->Write() ;
      for (Int_t i=1; i<3; i++)
	dEdx_vs_pT[i]->Write() ;
      dEdxMass_h[1]->Write() ; 
      //
      // EXCLUSIVE KPi
      myExclusiveMuMuKPiMass->Write() ;
      hmyPsiPKPiMassBaseSelAlt->Write() ;  
      kaonP_vs_pionP->Write() ; kaonPt_vs_pionPt->Write() ; trk_dxy->Write() ; trk_dz->Write() ;
      hmyPsiPKPiMass_B0pTtest->Write() ; hPionFromKstar_p->Write() ; hKaonFromKstar_p->Write() ;
      hmyPsiPKPiMass_deltaRtest->Write() ; hmyPsiPKPiMass_B0cTautest->Write() ;
      hmyPsiPKPiMass->Write() ; hmyPsiPKPiMass_KPiExchanged->Write() ; hmyPsiPKPiMass_bothComb->Write() ;
      hB0VtxCL->Write() ; hB0_Pt->Write() ;
      hB0CTauPVSignif->Write() ; hB0CTauPVXSignif->Write() ;
      hPsiPrimefromBMass->Write() ; hPsiPrime_fromB_Pt->Write() ; hmyPsiP_Pt->Write() ;
      hPionfromBpt->Write() ; hKaonfromBpt->Write() ;
      myExclMuMuKPiMassSel->Write() ; hmyPsiPKPiMass_bothCombSel->Write() ;
      psi2SPiSqMass_vs_KPiSqMass_ABaseC->Write() ;
      hmyPsiPMassSel->Write() ; hPsiPrimefromBMassSel->Write() ;
      hmyPsiPKPiMassZoom->Write() ; hmyPsiPKPiMassAlt->Write() ; hmyPsiPKPiMassAltZoom->Write() ; 
      hmyPiPiMass->Write() ; hmyPiKMass->Write() ; hmyKKMass->Write() ; hmyPsiPPiPiMass->Write() ; hmyPsiPPiKMass->Write() ; hmyPsiPPiKMass_zoom->Write() ; hmyPsiPKKMass->Write() ; hmyPsiPPhiMass->Write() ; hmyPsiPPiKMass_fromBs->Write() ; hmyPsiPPiKMass_fromBs_zoom->Write() ;
      hmyPiPiMass_vs_KKMass->Write() ; hmyPsiPPiPiMass_vs_PiPiMass->Write() ; hmyPsiPKKMass_vs_KKMass->Write() ;
      psi2SPiSqMass_vs_KPiSqMass_AC->Write() ;
      //
      // MC analysis
      nMCB0_BT_h->Write() ;
      psi2SPiSqMass_vs_KPiSqMass_BT_gen->Write() ; KPiSqMass_vs_psi2SKSqMass_BT_gen->Write() ;
      nB0AC_noMassWin_h->Write() ; nB0AC_h->Write() ; nB0AC_signalWin_h->Write() ; nB0AC_purityWin_h->Write() ; 
      hB0Mass_1B0->Write() ; hB0Mass_noTwins_noSignalWinNotTwins->Write() ; 
      for (Int_t i=0; i<2; i++)
	psi2SPiSqMass_vs_KPiSqMass_1B0[i]->Write() ;
      psi2S_helicityAngle->Write() ; cos_psi2S_helicityAngle->Write() ; cos_psi2S_helicityAngle_rot->Write() ; cos_psi2S_helicityAngle_rotVSnotRot->Write() ; planesAngle->Write() ;
      hB0CTau_1B0->Write() ; hB0CTau_1B0matched->Write() ;
      nMCB0_h->Write() ;
      B0_gen_p_h->Write() ; B0_gen_pT_h->Write() ; psi2S_gen_p_h->Write() ; psi2S_gen_pT_h->Write() ; K_gen_p_h->Write() ; K_gen_pT_h->Write() ; pi_gen_p_h->Write() ; pi_gen_pT_h->Write() ;
      psi2SPiSqMass_vs_KPiSqMass_gen->Write() ;
      psi2S_vs_mu_pT->Write() ;
      psi2S_helicityAngle_BT_gen->Write() ; cos_psi2S_helicityAngle_BT_gen->Write() ; cos_psi2S_helicityAngle_BT_gen_rot->Write() ; cos_psi2S_helicityAngle_BT_gen_rotVSnotRot->Write() ;
      cos_psi2S_helicityAngle_vs_psi2SPiMass->Write() ; cos_psi2S_helicityAngle_vs_KPiMass->Write() ; planesAngle_BT_gen->Write() ;
      B0_gen_mass_h->Write() ;
      piCh_vs_KCh_gen->Write() ;
      hMCDeltaRK_fromTwins->Write() ; hMCPtRecoVsPtTrueK_fromTwins->Write() ; hMCDeltaPtK_fromTwins->Write() ; hMCDeltaPtK_fromTwinsPiMatched->Write() ; hMCpTTrueK_fromTwins->Write() ;
      kaonDr_vs_pionDr_matched->Write() ; kaonP_vs_pionP_matched->Write() ; kaonPt_vs_pionPt_matched->Write() ;
      nB0ACInMC_h->Write() ;
      hMCDeltaR_1B0->Write() ; hMCDeltaPt_1B0->Write() ;
      for (Int_t i=0; i<3; i++) {
	hMCTruthMatching_1B0[0][i]->Write() ;
	if ( i < 2 )
	  hB0Mass_1B0matched[0][i]->Write() ;
      }      
      //
      for (Int_t i=0; i<4; i++) 
	B0CosAlpha_2Dvs3D[i]->Write() ;
      for (Int_t i=0; i<7; i++) 
	priVtx_z[i]->Write() ;
      for (Int_t j=0; j<2; j++) {
	for (Int_t i=0; i<3; i++) {
	  priVtx_delta[j][i]->Write() ;
	  priVtx_B0_delta[j][i]->Write() ;
	  for (Int_t k=0; k<2; k++)
	    priVtx_B0CosAlpha_delta[k][j][i]->Write() ;
	  priVtxB0Less_delta[j][i]->Write() ;
	  for (Int_t k=0; k<2; k++)
	    priVtxB0Less_B0CosAlpha_delta[k][j][i]->Write() ;
	}
	dz_pri_priB0Less[j]->Write() ;
	for (Int_t k=0; k<2; k++) {
	  dz_pri_priB0CosAlpha[k][j]->Write() ;  dz_pri_priB0LessB0CosAlpha[k][j]->Write() ; dz_priB0CosAlpha_priB0LessB0CosAlpha[k][j]->Write() ; dz_priB0Less_priB0LessB0CosAlpha[k][j]->Write() ;
	}      
      }
      priVtx_deltaX_test->Write() ; priVtx_deltaY_test->Write() ;
      //
      hDeltaB0Mass_2B0twin[0]->Write() ; hAlpha_2B0twin[0]->Write() ;
      for (Int_t i=0; i<2; i++) 
	for (Int_t j=0; j<2; j++) 
	  hB0Mass_2B0twin_hourglass[i][j]->Write() ;

      for (Int_t i = 0; i < abs(hMCDeltaR_2B0.size()); i++) {
	hMCDeltaR_2B0[i]->Write() ;
	hMCDeltaPt_2B0[i]->Write() ;
      }
      hMCDeltaPt_1vs2_2B0->Write() ; hMCDeltaR_1vs2_2B0->Write() ;
      hMCDeltaR_1vs2_2B0twin->Write() ; hMCDeltaPt_1vs2_2B0twin->Write() ; hB0VtxCL_1vs2_2B0twin->Write() ; hB0Mass_1vs2_2B0twin->Write() ; hB0CTau_1vs2_2B0twin->Write() ; hAlpha_1vs2_2B0twin->Write() ;
      for (Int_t i=0; i<2; i++) {
	hMCDeltaR_matchedTwin[i]->Write() ; hMCDeltaPt_matchedTwin[i]->Write() ; hB0VtxCL_matchedTwin[i]->Write() ; hB0Mass_matchedTwin[i]->Write() ; hB0CTau_matchedTwin[i]->Write();  hAlpha_matchedTwin[i]->Write() ;
      }
      nTwins_h[0]->Write() ;
      hB0Mass_2B0notTwin->Write() ;
      for (Int_t i=0; i<3; i++) {
	hMCTruthMatching_2B0notTwin[i]->Write() ;
	hB0Mass_2B0notTwin_matched[i]->Write() ;
	hMCTruthMatching_2B0twin[i]->Write() ;
	hB0Mass_2B0twin_matched[i]->Write() ;
      }
      hDeltaB0Mass_2B0twin[1]->Write() ; hAlpha_2B0twin[1]->Write() ;
      hMCDeltaRPi_2B0->Write() ; hMCDeltaRK_2B0->Write() ;
      nTwins_h[1]->Write() ;
      for (Int_t i=0; i<3; i++) {
	hMCTruthMatching_1B0[1][i]->Write() ;
	if ( i < 2 )
	  hB0Mass_1B0matched[1][i]->Write() ;
      }
      hDeltaB0Mass_2B0notTwin[1]->Write() ; hMCTruthMatching_2B0signalWin_notTwin->Write() ; hDeltaB0Mass_2B0notTwin[0]->Write() ;
      hB0Mass_2B0twin->Write() ; 
      // plots with candidate generated variables
      //hMCnGenX->Write() ; SpectrumFromBPlus->Write() ;  SpectrumFromBPlusCtaued->Write() ;
      // psi(2S) matching
      hMCDeltaRMuMuMatch->Write() ; hMCPtRecoVsPtTrueMuMuMatch->Write() ; hMCpTTrueMuMuMatch->Write() ;
      // pion matching
      hMCDeltaRPiMatch->Write() ; hMCPtRecoVsPtTruePiMatch->Write() ; hMCPtTruePiMatch->Write() ;
      hMCStripsHitsTruePiMatch->Write() ; hMCPixelHitsTruePiMatch->Write() ; hMCTotalHitsTruePiMatch->Write() ; hMCChiTruePiMatch->Write() ;
      //
      kaonDr_vs_pionDr_matched->Write() ; kaonP_vs_pionP_matched->Write() ; kaonPt_vs_pionPt_matched->Write() ;
      //
      // K*
      hmyPsiPKPiMass_fromKstar->Write() ; hmyPsiPKPiMass_fromKstar_ex_noRange->Write() ; hmyPsiPKPiMass_fromKstar_diff->Write() ; hmyPsiPKPiMass_fromKstar_ex->Write() ; hmyPsiPKPiMass_fromKstar_sb->Write() ; 
      hmyKPiMass->Write() ; hmyKPiMass_ex->Write() ; hmyKPiMass_sb->Write() ;
      hmyKPiMass_KpTcut->Write() ; hmyKPiMass_ex_KpTcut->Write() ; 
      hmyKPiMass_KpTcut_tightTr->Write() ; hmyKPiMass_ex_KpTcut_tightTr->Write() ; hmyKPiMass_KpTcut_HP->Write() ; hmyKPiMass_ex_KpTcut_HP->Write() ;

      hpsiPK_vs_KPi->Write() ; hpsiPK_vs_KPi_ex->Write() ; hpsiPK_vs_KPi_sb->Write() ;
      hPionFromKstar_pT->Write() ; hKaonFromKstar_pT->Write() ; hPionFromKstar_ex_pT->Write() ; hKaonFromKstar_ex_pT->Write() ; hPionFromKstar_sb_pT->Write() ; hKaonFromKstar_sb_pT->Write() ;
      dEdx_vs_pT_Kstar->Write() ; dEdx_vs_p_Kstar->Write() ; 
      dEdx_vs_pT_piFromKstar->Write() ; dEdx_vs_pT_KFromKstar->Write() ;
      dEdx_vs_p_piFromKstar->Write() ; dEdx_vs_p_KFromKstar->Write() ; dEdx_vs_p_piFromKstar_KpCut->Write() ; dEdx_vs_p_KFromKstar_KpCut->Write() ;
      dEdxK_vs_dEdxPi->Write() ; 
      piFromKstar_stripHits_h->Write() ; KFromKstar_stripHits_h->Write() ; piFromKstar_stripHits_vs_y_h->Write() ; KFromKstar_stripHits_vs_y_h->Write() ;
      dEdx_vs_pT_KhqFromKstar->Write() ;
      dEdx_vs_p_KhqFromKstar->Write() ; dEdxPi_pCut->Write() ; dEdxK_pCut->Write() ; dEdxPi_pCut_stripCut->Write() ; dEdxK_pCut_stripCut->Write() ;
      Armenteros_Kstar_orig->Write() ; Armenteros_Kstar->Write() ; Armenteros_Kstar_orig_ex->Write() ; Armenteros_Kstar_orig_sb->Write() ; alpha_vs_alphaEx->Write() ; 
      pT_diff->Write() ; pT_diff_ex->Write() ; pT_diff_sb->Write() ;
      //
      /*
      // my scan
      myMuMuPiMass->Write() ; myMuMuPiMass_pr->Write() ; myMuMuPiMass_quasipr->Write() ; myMuMuPiMass_nonpr->Write() ;
      myMuMuPiMass_barrel->Write() ;  myMuMuPiMass_barrel_pr->Write() ;  myMuMuPiMass_nonBarrel_pr->Write() ;  myMuMuPiMass_barrel_nonpr->Write() ;
      myMuMuPiMass_pr_DzCut->Write() ; myMuMuPiMass_pr_Dz01->Write() ; myMuMuPiMass_nonpr_Dz1inf->Write() ;
      myMuMuPiMass_DeltaR0p1_nonpr5->Write() ; myMuMuPiMass_DeltaR0p2_nonpr5->Write() ; myMuMuPiMass_DeltaR0p1_nonpr8->Write() ; myMuMuPiMass_DeltaR0p2_nonpr8->Write() ; 
      for (Int_t i = 0; i < abs(mode.size()); ++i) {
	psi2S_DeltaR_pi_h[i]->Write() ;
	//Z_cosPA[i].first->Write() ;  Z_cosPA[i].second->Write() ;
	Z_cosPA[i]->Write() ;
	Z_DeltaR_vs_cosPA[i]->Write() ;
	for (Int_t k = 0; k < abs(deltaR_piPsi2S.size()); ++k) {
	  for (Int_t j = 0; j < abs(pi_pT.size()); ++j) 
	    myMuMuPiMass_piPt[i][j][k]->Write() ;
	  for (Int_t j = 0; j < abs(Z_pT.size()); ++j) 
	    myMuMuPiMass_ZPt[i][j][k]->Write() ;
	  for (Int_t j = 0; j < abs(piAndZ_pT.size()); ++j) 
	    myMuMuPiMass_piAndZPt[i][j][k]->Write() ;
	}
      }
      */
      psi2S_DeltaR_pi_vs_mumupiMass_nonpr->Write() ;
      
      for (Int_t i = 0; i < rapid_bins; i++) {
	psiPrime_rapidBin[i]->Write() ;
	//psiPrime_2011_rapidBin[i]->Write() ; 
      }
      //
      // reco variables
      nZ_h->Write() ; nB0_h->Write() ;
      dEdxMass_h[0]->Write() ; dEdxMass_vs_pT->Write() ; dEdx_vs_pT[0]->Write() ; dEdx_vs_p->Write() ;
      B0_Mass_h ->Write() ; myB0_Mass_h ->Write() ;
      // selection variables
      // Z+ variables
      Z_VtxCL_h->Write() ;
      Z_CTauBS_h->Write() ;
      Z_CTauBSE_h->Write() ;
      Z_CTauBS_over_Z_CTauBSE_h->Write() ;
      Z_pT_h->Write() ;
      for (Int_t i=0; i<3; ++i)
	Z_pointingAngle[i]->Write() ;
      //for (Int_t i = 0; i < 4; ++i)
      Z_vs_pi_pT_h->Write() ;
      // psi(2S) variables
      psi2S_VtxCL_h->Write() ;  psi2S_VtxCL_HLT5_h->Write() ;
      psi2S_VtxCL_1_h->Write() ;  psi2S_VtxCL_2_h->Write() ;  psi2S_VtxCL_diff_h->Write() ;
      psi2SMuonTrigMatch_h->Write() ;
      psi2S_y_h->Write() ;
      psi2S_HLT5_yMu_h->Write() ; 
      psi2S_pT_h->Write() ;
      /*
      for (Int_t i = 0; i < abs(psi2S_2011_pT.size()); ++i) {
	psi2S_2011_pT_h[i]->Write() ;
	for (Int_t j = 0; j < abs(HLT_v[i].size()); ++j) {
	  psi2S_2011v_pT_h[i][j]->Write() ;
	  myInclusiveMuMuMass_2011_HLT[i][j]->Write() ;
	}
	}*/

      //psi2S_DeltaR_pi_h->Write() ; psi2S_DeltaR_pi_pr_h->Write() ; psi2S_DeltaR_pi_quasipr_h->Write() ; psi2S_DeltaR_pi_nonpr_h->Write() ; 
      //psi2S_pointingAngle->Write() ; psi2S_pointingAngle_sign->Write() ; psi2S_pointingAngle_sb->Write() ;
      pT_vs_y_psi2S_h->Write() ;
      nPsi2S_h->Write() ;
      // mu variables
      mu_px_psi2S->Write() ; mu_px->Write() ; fitmu_px->Write() ;
      mu1_pixelHits_h->Write() ; mu2_pixelHits_h->Write() ;
      mu1_stripHits_h->Write() ; mu2_stripHits_h->Write() ;
      mu1_charge_h->Write() ; mu2_charge_h->Write() ;
      mu1_Chi2N_h->Write() ; mu2_Chi2N_h->Write() ;
      mu1Tr_Chi2N_h->Write() ; mu2Tr_Chi2N_h->Write() ;
      mu1GlTr_Chi2N_h->Write() ; mu2GlTr_Chi2N_h->Write() ;
      mu1_Dxy_h->Write() ; mu2_Dxy_h->Write() ; 
      mu1_Dz_h->Write() ; mu2_Dz_h->Write() ; 
      mu1_y_h->Write() ; mu2_y_h->Write() ;
      mu1_pT_h->Write() ; mu2_pT_h->Write() ;
      mu1_pT_test_h->Write() ; mu2_pT_test_h->Write() ;
      mu1_eta_test_h->Write() ; mu2_eta_test_h->Write() ;
      mu1_DeltaR_mu2_h->Write() ; //mu2_DeltaR_mu1_h->Write() ; 
      //mu1_vs_mu2_pT_h->Write() ;
      pT_vs_y_mu1_h->Write() ; pT_vs_y_mu2_h->Write() ;
      pT_vs_eta_mu_h->Write() ; pT10_vs_eta_mu_h->Write() ;
      // pi+ variables
      pi_pixelHits_h->Write() ;
      pi_stripHits_h->Write() ;
      pi_Chi2N_h->Write() ; pi_Chi2N_test_h->Write() ;
      pi_y_h->Write() ; pi_y_test_h->Write() ;
      pi_pT_h->Write() ; pi_pT_test_h->Write() ;
      pT_vs_pixelHits_pi_h->Write() ; pT_vs_stripHits_pi_h->Write() ;
      pT_vs_y_pi_h->Write() ;
      pi_Dz_h->Write() ; pi_Dxy_h->Write() ; pi_Dz_pr_h->Write() ; pi_Dxy_pr_h->Write() ; pi_Dz_nonpr_h->Write() ; pi_Dxy_nonpr_h->Write() ;
      //.............
      //
      OutFile->Print();
      fOutput->Add(OutFile); 
      myInclusiveMuMuMass->SetDirectory(0); // any hname
      gDirectory = savedir;
      fOut->Close();
      //
    }
  
  if (SaraOut) 
    {      
      SaraOut->cd();
      //// Write histograms in output file
      gStyle->SetOptStat(111111) ;
      //
      // hname->Write();
      //
      // Sara plots
      for (Int_t i = 0; i < abs(Z_CTauBSOverE.size()); ++i) 
	MuMuKMass_ZCTauBSOverE[i]->Write() ;
      for (Int_t i = 0; i < abs(mu_pT.size()); ++i) {
        MuMuKMass_mu1Pt[i]->Write() ;
	MuMuKMass_mu2Pt[i]->Write() ;
      }
      for (Int_t i = 0; i < abs(trk_pT.size()); ++i) {
        MuMuKMass_trk1Pt[i]->Write() ;
	MuMuKMass_trk2Pt[i]->Write() ;
      }
      for (Int_t i = 0; i < abs(Z_CL.size()); ++i) 
        MuMuKMass_ZCL[i]->Write() ;
      //.............
      //
      SaraFile->Print();
      fOutput->Add(SaraFile); 
      MuMuKMass_ZCL[0]->SetDirectory(0); // any hname
      gDirectory = savedir;
      SaraOut->Close();
      //
    }

  if ( AlexisOut ) 
    {      
      AlexisOut->cd();
      //// Write histograms in output file
      gStyle->SetOptStat(111111) ;
      //
      // hname->Write();
      //
      // Alexis plots
      nMCB0_BT_h->Write() ;
      psi2SPiSqMass_vs_KPiSqMass_BT_gen->Write() ; KPiSqMass_vs_psi2SKSqMass_BT_gen->Write() ;
      nB0_h->Write() ; nB0AC_noMassWin_h->Write() ; nB0AC_h->Write() ; nB0AC_signalWin_h->Write() ; nB0AC_purityWin_h->Write() ; nB0ACInMC_h->Write() ;
      hmyPsiPMassSel->Write() ; hmyPsiPKPiMassZoom->Write() ;  hmyPsiPKPiMassAlt->Write() ; hmyPsiPKPiMassAltZoom->Write() ;
      hmyPiPiMass->Write() ; hmyPiKMass->Write() ; hmyKKMass->Write() ; hmyPsiPPiPiMass->Write() ; hmyPsiPPiKMass->Write() ; hmyPsiPPiKMass_zoom->Write() ; hmyPsiPKKMass->Write() ; hmyPsiPPhiMass->Write() ; hmyPsiPPiKMass_fromBs->Write() ; hmyPsiPPiKMass_fromBs_zoom->Write() ;
      hmyPiPiMass_vs_KKMass->Write() ; hmyPsiPPiPiMass_vs_PiPiMass->Write() ; hmyPsiPKKMass_vs_KKMass->Write() ;
      psi2SPiSqMass_vs_KPiSqMass_AC->Write() ;
      // 
      hB0Mass_1B0->Write() ; hB0Mass_noTwins_noSignalWinNotTwins->Write() ;
      for (Int_t i=0; i<2; i++)
	psi2SPiSqMass_vs_KPiSqMass_1B0[i]->Write() ; 
      psi2S_helicityAngle->Write() ; cos_psi2S_helicityAngle->Write() ; cos_psi2S_helicityAngle_rot->Write() ; cos_psi2S_helicityAngle_rotVSnotRot->Write() ; planesAngle->Write() ;
      hB0CTau_1B0->Write() ; hB0CTau_1B0matched->Write() ;
      hMCDeltaR_1B0->Write() ; hMCDeltaPt_1B0->Write() ; 
      for (Int_t i=0; i<3; i++) {
	hMCTruthMatching_1B0[0][i]->Write() ;
	if ( i < 2 )
	  hB0Mass_1B0matched[0][i]->Write() ;
      }      
      //
      for (Int_t i=0; i<4; i++) 
	B0CosAlpha_2Dvs3D[i]->Write() ;
      for (Int_t i=0; i<7; i++) 
	priVtx_z[i]->Write() ;
      for (Int_t j=0; j<2; j++) { 
	for (Int_t i=0; i<3; i++) {
	  priVtx_delta[j][i]->Write() ;
	  priVtx_B0_delta[j][i]->Write() ;
	  for (Int_t k=0; k<2; k++)  
	    priVtx_B0CosAlpha_delta[k][j][i]->Write() ;
	  priVtxB0Less_delta[j][i]->Write() ;
	  for (Int_t k=0; k<2; k++)  
	    priVtxB0Less_B0CosAlpha_delta[k][j][i]->Write() ;
	}
	dz_pri_priB0Less[j]->Write() ;
	for (Int_t k=0; k<2; k++) {
	  dz_pri_priB0CosAlpha[k][j]->Write() ;  dz_pri_priB0LessB0CosAlpha[k][j]->Write() ; dz_priB0CosAlpha_priB0LessB0CosAlpha[k][j]->Write() ; dz_priB0Less_priB0LessB0CosAlpha[k][j]->Write() ;
	}      
      }
      priVtx_deltaX_test->Write() ; priVtx_deltaY_test->Write() ;
      //
      hDeltaB0Mass_2B0twin[0]->Write() ; hAlpha_2B0twin[0]->Write() ;
      for (Int_t i=0; i<2; i++) 
	for (Int_t j=0; j<2; j++) 
	  hB0Mass_2B0twin_hourglass[i][j]->Write() ;
	    
      for (Int_t i = 0; i < abs(hMCDeltaR_2B0.size()); i++) {
	hMCDeltaR_2B0[i]->Write() ;
	hMCDeltaPt_2B0[i]->Write() ;
      }
      hMCDeltaPt_1vs2_2B0->Write() ; hMCDeltaR_1vs2_2B0->Write() ;
      hMCDeltaR_1vs2_2B0twin->Write() ; hMCDeltaPt_1vs2_2B0twin->Write() ; hB0VtxCL_1vs2_2B0twin->Write() ; hB0Mass_1vs2_2B0twin->Write() ; hB0CTau_1vs2_2B0twin->Write() ; hAlpha_1vs2_2B0twin->Write() ;
      psi2SPiSqMass_vs_KPiSqMass_gen->Write() ;
      psi2S_vs_mu_pT->Write() ;
      psi2S_helicityAngle_BT_gen->Write() ; cos_psi2S_helicityAngle_BT_gen->Write() ; cos_psi2S_helicityAngle_BT_gen_rot->Write() ; cos_psi2S_helicityAngle_BT_gen_rotVSnotRot->Write() ;
      cos_psi2S_helicityAngle_vs_psi2SPiMass->Write() ; cos_psi2S_helicityAngle_vs_KPiMass->Write() ; planesAngle_BT_gen->Write() ;
      B0_gen_mass_h->Write() ;
      piCh_vs_KCh_gen->Write() ; 
      hMCDeltaRK_fromTwins->Write() ; hMCPtRecoVsPtTrueK_fromTwins->Write() ; hMCDeltaPtK_fromTwins->Write() ; hMCDeltaPtK_fromTwinsPiMatched->Write() ; hMCpTTrueK_fromTwins->Write() ;
      kaonDr_vs_pionDr_matched->Write() ; kaonP_vs_pionP_matched->Write() ; kaonPt_vs_pionPt_matched->Write() ;
      for (Int_t i=0; i<2; i++) {
	hMCDeltaR_matchedTwin[i]->Write() ; hMCDeltaPt_matchedTwin[i]->Write() ; hB0VtxCL_matchedTwin[i]->Write() ; hB0Mass_matchedTwin[i]->Write() ; hB0CTau_matchedTwin[i]->Write(); hAlpha_matchedTwin[i]->Write() ;
      }
      nTwins_h[0]->Write() ; 
      hDeltaB0Mass_2B0notTwin[1]->Write() ; hB0Mass_2B0notTwin->Write() ;
      hB0Mass_2B0twin->Write() ;
      for (Int_t i=0; i<3; i++) {
	hMCTruthMatching_2B0notTwin[i]->Write() ;
	hB0Mass_2B0notTwin_matched[i]->Write() ;
	hMCTruthMatching_2B0twin[i]->Write() ;
	hB0Mass_2B0twin_matched[i]->Write() ;
      }
      nTwins_h[1]->Write() ;
      for (Int_t i=0; i<3; i++) {
	hMCTruthMatching_1B0[1][i]->Write() ;
	if ( i < 2 )
	  hB0Mass_1B0matched[1][i]->Write() ;
      }
      hMCTruthMatching_2B0signalWin_notTwin->Write() ; hDeltaB0Mass_2B0notTwin[0]->Write() ; hMCDeltaRPi_2B0->Write() ; hMCDeltaRK_2B0->Write() ;
      hDeltaB0Mass_2B0twin[1]->Write() ; hAlpha_2B0twin[1]->Write() ;
      //.............
      //
      AlexisFile->Print();
      fOutput->Add(AlexisFile); 
      hmyPsiPMassSel->SetDirectory(0); // any hname
      gDirectory = savedir;
      AlexisOut->Close();
      //
    }

  if ( tightCutsOut ) 
    {      
      tightCutsOut->cd();
      //
      hmyPsiPMassSel->Write() ; hmyPsiPKPiMassZoom->Write() ; hmyPsiPKPiMassAlt->Write() ; hmyPsiPKPiMassAltZoom->Write() ; 
      psi2SPiSqMass_vs_KPiSqMass_AC->Write() ;
      nB0AC_h->Write() ; nB0AC_signalWin_h->Write() ; nB0AC_purityWin_h->Write() ;
      for (Int_t i = 0; i < abs(piPt_B0.size()); ++i)
	myMuMuKPiMass_piPt[i]->Write() ;
      for (Int_t i = 0; i < abs(KPt_B0.size()); ++i)
	myMuMuKPiMass_KPt[i]->Write() ;
      for (Int_t i = 0; i < abs(vtxCL_B0.size()); ++i)
	myMuMuKPiMass_vtxCL[i]->Write() ;
      for (Int_t i = 0; i < abs(cTau_B0.size()); ++i)
	myMuMuKPiMass_cTau[i]->Write() ;
      for (Int_t i = 0; i < abs(trksPt_B0VtxCL_B0CTau.size()); ++i)
	myMuMuKPiMass_trksAndB0vtxCLAndB0Ctau[i]->Write() ;

      priVtx_deltaZ->Write() ;
      B0_pointingAngle_PVvsB0LessPV[0]->Write() ; B0_pointingAngle_PVvsB0LessPV[1]->Write() ;
      z_PVvsB0LessPV[0]->Write() ; z_PVvsB0LessPV[1]->Write() ;
      B0_PVCTau_h->Write() ; B0_CosAlphaCTau_h->Write() ; B0_PVXCTau_h->Write() ; B0_PVX3DCTau_h->Write() ;
      priVtx_matched->Write() ;
      //.............
      //
      tightCutsFile->Print();
      fOutput->Add(tightCutsFile); 
      hmyPsiPMassSel->SetDirectory(0); // any hname
      gDirectory = savedir;
      tightCutsOut->Close();
      //
    }
  //
}

void psiPrimePiK_MC::Terminate()
{
  // The Terminate() function is the last function to be called during a query. It always runs on the client, it can be used to present the results graphically or save the results to file.
  cout <<"\nTerminate:" <<endl ;

  TString dir = "";

  TString option = GetOption();
  if (option == "officialMC") {
    //dir = "officialMC" ;
    //dir = "officialMC/0p9" ; 
    dir = "officialMC_extended" ;
  }
  if ( option.Contains("MC") )
    dir.Prepend("mc/");
  else
   dir.Prepend("data/");

  gSystem->mkdir(TString::Format("plots/%s", dir.Data()), true);
  
  TString fileName = "26Jan_histos.root"; fileName = "HLT5.root";
  if ( (OutFile = dynamic_cast<TProofOutputFile*>(fOutput->FindObject( fileName.Data() ))) ) {

    TString outputFile(OutFile->GetOutputFileName());
    TString outputName(OutFile->GetName()); outputName += ".root";
    Printf("\nOutput file found: %s", outputFile.Data());
    
    // Open the file
    AlexisOut = TFile::Open(outputFile);
    if (AlexisOut) {
      Printf("Managed to open file: %s\n", outputFile.Data());
      ////////// Read histograms //////////
      // 1 B0
      if ( (hMCTruthMatching_1B0[0][0] = (TH1I*) AlexisOut->Get("MCTruthMatching_1B0")) ) {
	cout <<"MC matching efficiency for only 1 B0 = " <<( 1. - hMCTruthMatching_1B0[0][0]->GetBinContent(hMCTruthMatching_1B0[0][0]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_1B0[0][0]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_1B0\" not found");
      if ( (hMCTruthMatching_1B0[0][1] = (TH1I*) AlexisOut->Get("MCTruthMatching_1B0_dauCh")) ) {
	cout <<"MC matching (with daughters charge) efficiency for only 1 B0 = " <<( 1. - hMCTruthMatching_1B0[0][1]->GetBinContent(hMCTruthMatching_1B0[0][1]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_1B0[0][1]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_1B0_dauCh\" not found");
      //
      if ( (hMCTruthMatching_1B0[0][2] = (TH1I*) AlexisOut->Get("MCTruthMatching_1B0fromDau")) ) {
	cout <<"MC truth daughters matching efficiency for only 1 B0 = " <<( 1. - hMCTruthMatching_1B0[0][2]->GetBinContent(hMCTruthMatching_1B0[0][2]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_1B0[0][2]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_1B0fromDau\" not found");
      //

      // 2 B0 not twin
      if ( (hMCTruthMatching_2B0notTwin[0] = (TH1I*) AlexisOut->Get("MCTruthMatching_2B0notTwin")) ) {
	cout <<"MC matching efficiency for 2 B0 not twin = " <<( 1. - hMCTruthMatching_2B0notTwin[0]->GetBinContent(hMCTruthMatching_2B0notTwin[0]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_2B0notTwin[0]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_2B0notTwin\" not found");
      if ( (hMCTruthMatching_2B0notTwin[1] = (TH1I*) AlexisOut->Get("MCTruthMatching_2B0notTwin_dauCh")) ) {
	cout <<"MC matching (with daughters charge) efficiency for 2 B0 not twin = " <<( 1. - hMCTruthMatching_2B0notTwin[1]->GetBinContent(hMCTruthMatching_2B0notTwin[1]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_2B0notTwin[1]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_2B0notTwin_dauCh\" not found");
      //
      if ( (hMCTruthMatching_2B0notTwin[2] = (TH1I*) AlexisOut->Get("MCTruthMatching_2B0notTwin_fromDau")) ) {
	cout <<"MC truth daughters matching efficiency for 2 B0 not twin = " <<( 1. - hMCTruthMatching_2B0notTwin[2]->GetBinContent(hMCTruthMatching_2B0notTwin[2]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_2B0notTwin[2]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_2B0notTwin_fromDau\" not found");
      //
      
      // 2 B0 twin 
      if ( (hMCTruthMatching_2B0twin[0] = (TH1I*) AlexisOut->Get("MCTruthMatching_2B0twin")) ) {
	cout <<"MC matching efficiency for 2 B0 twin = " <<( 1. - hMCTruthMatching_2B0twin[0]->GetBinContent(hMCTruthMatching_2B0twin[0]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_2B0twin[0]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_2B0twin\" not found");
      if ( (hMCTruthMatching_2B0twin[1] = (TH1I*) AlexisOut->Get("MCTruthMatching_2B0twin_dauCh")) ) {
	cout <<"MC matching (with daughters charge) efficiency for 2 B0 twin = " <<( 1. - hMCTruthMatching_2B0twin[1]->GetBinContent(hMCTruthMatching_2B0twin[1]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_2B0twin[1]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_2B0twin_dauCh\" not found");
      //
      if ( (hMCTruthMatching_2B0twin[2] = (TH1I*) AlexisOut->Get("MCTruthMatching_2B0twin_fromDau")) ) {
	cout <<"MC truth daughters matching efficiency for 2 B0 twin = " <<( 1. - hMCTruthMatching_2B0twin[2]->GetBinContent(hMCTruthMatching_2B0twin[2]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_2B0twin[2]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_2B0twin_fromDau\" not found");
      //
      
      // signal win
      cout <<"\n\nMC matching in signal mass window:\n" <<endl ;
      if ( (hMCTruthMatching_1B0[1][0] = (TH1I*) AlexisOut->Get("MCTruthMatching_1B0signalWin")) ) {
	cout <<"MC matching efficiency for only 1 B0 in signal = " <<( 1. - hMCTruthMatching_1B0[1][0]->GetBinContent(hMCTruthMatching_1B0[1][0]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_1B0[1][0]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_1B0signalWin\" not found");
      if ( (hMCTruthMatching_1B0[1][1] = (TH1I*) AlexisOut->Get("MCTruthMatching_1B0signalWin_dauCh")) ) {
	cout <<"MC matching (with daughters charge) efficiency for only 1 B0 in signal = " <<( 1. - hMCTruthMatching_1B0[1][1]->GetBinContent(hMCTruthMatching_1B0[1][1]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_1B0[1][1]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_1B0signalWin_dauCh\" not found");
      if ( (hMCTruthMatching_1B0[1][2] = (TH1I*) AlexisOut->Get("MCTruthMatching_1B0signalWin_fromDau")) ) {
	cout <<"MC truth daughters matching efficiency for only 1 B0 in signal = " <<( 1. - hMCTruthMatching_1B0[1][2]->GetBinContent(hMCTruthMatching_1B0[1][2]->GetXaxis()->FindBin(0.)) / hMCTruthMatching_1B0[1][2]->GetEntries() )*100 <<"%" <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_1B0signalWin_fromDau\" not found");
      //
      if ( (hMCTruthMatching_2B0signalWin_notTwin = (TH1I*) AlexisOut->Get("MCTruthMatching_2B0signalWin_notTwin")) ) {
	cout <<"\nMC truth matching efficiency for 2 B0 not twin in signal = " <<( 1 - hMCTruthMatching_2B0signalWin_notTwin->GetBinContent(hMCTruthMatching_2B0signalWin_notTwin->GetXaxis()->FindBin(0.)) / hMCTruthMatching_2B0signalWin_notTwin->GetEntries() )*100 <<"%"  <<endl ;
      } else Warning("Terminate", "histogram \"MCTruthMatching_2B0notTwin\" not found");
      //
      ////////// Plot histograms //////////
      //
      vector< pair<TString, TString> > name_TH ;
      gStyle->SetOptStat( 10 ) ;

      name_TH.push_back( make_pair("B0_multiplicity_AC","TH1I") ); name_TH.push_back( make_pair("B0twins_multiplicity_AC","TH1I") );
      name_TH.push_back( make_pair("myPsiPKPiMassAlt","TH1F") );
      name_TH.push_back( make_pair("myPiPiMass","TH1F") ); name_TH.push_back( make_pair("myPiKMass","TH1F") ); name_TH.push_back( make_pair("myKKMass","TH1F") );
      name_TH.push_back( make_pair("myPsiPPiPiMass","TH1F") ); name_TH.push_back( make_pair("myPsiPPiKMass","TH1F") ); name_TH.push_back( make_pair("myPsiPPiKMass_zoom","TH1F") ); name_TH.push_back( make_pair("myPsiPKKMass","TH1F") ); name_TH.push_back( make_pair("myPsiPPhiMass","TH1F") );
      name_TH.push_back( make_pair("myPsiPPiKMass_fromBs","TH1F") ); name_TH.push_back( make_pair("myPsiPPiKMass_fromBs_zoom","TH1F") ); name_TH.push_back( make_pair("myPiPiMass_vs_KKMass","TH2F") ); name_TH.push_back( make_pair("myPsiPPiPiMass_vs_PiPiMass","TH2F") ); name_TH.push_back( make_pair("myPsiPKKMass_vs_KKMass","TH2F") );

      for (Int_t i = 0; i < abs(name_TH.size()); ++i) {
	TH1I* histo1I = 0; TH1F* histo1F = 0; TH2F* histo2F = 0;
	if ( name_TH[i].second.EqualTo("TH1I") )
	  if ( (histo1I = (TH1I*) AlexisOut->Get( name_TH[i].first.Data() )) ) {
	    histo1I->Draw() ;
	    histo1I->Draw("text30same") ; gPad->SetLogy() ;
	    gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), histo1I->GetName()) );
	  } else Warning("Terminate", TString::Format("histogram \"%s\" not found", name_TH[i].first.Data() )) ;
	else if ( name_TH[i].second.EqualTo("TH1F") )
	  if ( (histo1F = (TH1F*) AlexisOut->Get( name_TH[i].first.Data() )) ) {
	    histo1F->Draw() ; gPad->SetLogy(0) ;
	    gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), histo1F->GetName()) );
	  } else Warning("Terminate", TString::Format("histogram \"%s\" not found", name_TH[i].first.Data() )) ;
	else if ( name_TH[i].second.EqualTo("TH2F") )
	  if ( (histo2F = (TH2F*) AlexisOut->Get( name_TH[i].first.Data() )) ) {
	    histo2F->Draw() ;
	    histo2F->Draw("colz") ; gPad->SetLogy(0) ;
	    gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), histo2F->GetName()) );
	  } else Warning("Terminate", TString::Format("histogram \"%s\" not found", name_TH[i].first.Data() )) ;
	else Warning("Terminate", "histogram type is neither \"TH1I\" nor \"TH1F\" nor \"TH2F\"!") ;
      }

      TCanvas *Bs_reflection = new TCanvas("Bs_reflection","B_{s} reflection", 2*600, 2*400); Bs_reflection->Divide(2,2);
      TString BsReflection_hist[] = {"myKKMass","myPsiPPhiMass","myPsiPPiKMass","myPsiPPiKMass_zoom"};
      Float_t lineX[2][2] = {{phi_left,phi_right}, {Bs_left,Bs_right}};
      TH1F* histo[4] = {0,0,0,0};
      for (Int_t i=0; i<4; ++i)
	if ( (histo[i] = (TH1F*) AlexisOut->Get(BsReflection_hist[i])) ) {	 
	  Bs_reflection->cd(i+1);
	  histo[i]->Draw(); Bs_reflection->Update();
	  if (i<2) {
	    TLine* left = new TLine(lineX[i][0],gPad->GetUymin(),lineX[i][0],gPad->GetUymax()); left->SetLineColor(kRed); left->SetLineWidth(2);  left->Draw();
	    TLine* right = new TLine(lineX[i][1],gPad->GetUymin(),lineX[i][1],gPad->GetUymax()); right->SetLineColor(kRed); right->SetLineWidth(2); right->Draw();
	    if (i==0) { gPad->SetLogx(); histo[i]->GetXaxis()->SetRangeUser(0.9,histo[i]->GetXaxis()->GetXmax()) ; }
	  }
	  else if (i==2) { 
	      TH1F* histoSame = (TH1F*) AlexisOut->Get("myPsiPPiKMass_fromBs") ;
	      if (histoSame) {
		histoSame->SetFillColor(kRed); histoSame->Draw("same");
		//gPad->Update(); gPad->BuildLegend();
	      } else Warning("Terminate", "histogram \"myPsiPPiKMass_fromBs\" not found") ;

	  } else if (i==3) { 
	    histo[i]->SetMinimum(0);
	    TH1F* histoSame = (TH1F*) AlexisOut->Get("myPsiPPiKMass_fromBs_zoom") ;
	    if (histoSame) {
	      histoSame->SetFillColor(kRed); histoSame->Draw("same"); //gPad->BuildLegend();
	    } else Warning("Terminate", "histogram \"myPsiPPiKMass_fromBs_zoom\" not found") ;
	  }
	} else Warning("Terminate", TString::Format("histogram \"%s\" not found", BsReflection_hist[i].Data() )) ;

      Bs_reflection->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), Bs_reflection->GetName()) );
       
      /*
      if ( (hmyPsiPKPiMassSelAltZoom = (TH1F*) AlexisOut->Get("myPsiPKPiMassSelAltZoom")) ) {

	Double_t xMin = hmyPsiPKPiMassSelAltZoom->GetXaxis()->GetXmin();
	Double_t xMax = hmyPsiPKPiMassSelAltZoom->GetXaxis()->GetXmax();
	Int_t nBins = hmyPsiPKPiMassSelAltZoom->GetNbinsX();
	Double_t bin_width = (xMax-xMin)/nBins ;
  
	RooRealVar xVar("xVar", hmyPsiPKPiMassSelAltZoom->GetXaxis()->GetTitle(), xMin, xMax) ;  xVar.setBins(nBins) ;  
	RooDataHist *RooPsiPKPiMassSelAltZoom = new RooDataHist(hmyPsiPKPiMassSelAltZoom->GetName(), hmyPsiPKPiMassSelAltZoom->GetTitle(), RooArgSet(xVar), Import(*hmyPsiPKPiMassSelAltZoom, kFALSE)) ;
	//RooConstVar nEntries("Entries", "Total number of entries", hmyPsiPKPiMassSelAltZoom->Integral()) ;
  
	Double_t fitMax = xMax, fitMin = xMin ;

	RooRealVar mean("m", "mean", 5.28, 5.25, 5.31) ;
	RooRealVar sigma1("#sigma_{1}", "sigma1", 0.01, 0.001, 0.02);
	//RooRealVar sigma2("#sigma_{2}", "sigma2", 0.02, 0.01, 0.2);
	RooGaussian sigPDF("sigPDF", "Signal component", xVar, mean, sigma1); sigma1.SetNameTitle("#sigma","sigma");
	//RooGaussian sig1PDF("sig1PDF", "Signal component 1", xVar, mean, sigma1);
	//RooGaussian sig2PDF("sig2PDF", "Signal component 2", xVar, mean, sigma2);
	//RooRealVar ratio("sig1frac", "fraction of component 1 in signal", 0.1, 0, 0.2);
	//RooAddPdf sigPDF("sigPDF", "Signal component", RooArgList(sig1PDF, sig2PDF), ratio);

	RooRealVar nSig("N_{Sig}", "Number of signal candidates", 2.5e+3, 1e+2, 1e+4);

	RooRealVar c1("c_{1}", "c1", +0.02, -10, 10);
	RooRealVar c2("c_{2}", "c2", -0.001, -10, 10);
	RooChebychev bkgPDF("bkgPDF", "bkgPDF", xVar, RooArgSet(c1,c2));
	RooRealVar nBkg("nBkg", "Bkg component", 1e+3, 1e+2, 1e+4);     
  
	RooAddPdf *totalPDF = new RooAddPdf("totalPDF", "totalPDF", RooArgList(sigPDF, bkgPDF), RooArgList(nSig, nBkg));  // with 5 sigma
	cout <<"\n\n---------------  Fitting histo \"" <<RooPsiPKPiMassSelAltZoom->GetTitle() <<"\"  ---------------\n" <<endl;
  
	totalPDF->fitTo(*RooPsiPKPiMassSelAltZoom, Range(fitMin, fitMax), Extended(kTRUE));

	RooPlot* xframe = xVar.frame();
	//xframe->SetAxisRange(fitMin, xMax) ;
	xframe->SetTitle( RooPsiPKPiMassSelAltZoom->GetTitle() ) ;
	xframe->SetYTitle( TString::Format("Candidates / %.1f MeV/c^{2}", bin_width*1000) ) ;
  
	RooPsiPKPiMassSelAltZoom->plotOn(xframe);
	totalPDF->plotOn(xframe);
	totalPDF->plotOn(xframe, Components(RooArgSet(sigPDF)), LineColor(kRed));
	//totalPDF->plotOn(xframe, Components(RooArgSet(sig1PDF)), LineColor(kOrange));
	//totalPDF->plotOn(xframe, Components(RooArgSet(sig2PDF)), LineColor(kMagenta));

	totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig,mean,sigma1)), ShowConstants(kTRUE), Layout(0.55,0.9,0.9)) ;
	//totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig,mean,sigma1,sigma2)), Layout(0.55,0.9,0.9)) ;

	xframe->Draw();
	//gPad->SaveAs( TString::Format("plots/%s/B0_mass_1B0_fit.png", dir.Data()) );

      } else Warning("Terminate", "histogram \"myPsiPKPiMassSelAltZoom\" not found");
      //
      
      // 1 B0
      if ( (hB0Mass_1B0 = (TH1F*) AlexisOut->Get("B0Mass_1B0")) ) {
	if ( (hB0Mass_1B0matched[0][0] = (TH1F*) AlexisOut->Get("B0Mass_1B0matched")) ) {
	  TH1F* diff_hist = (TH1F*) hB0Mass_1B0->Clone("diff_hist"); diff_hist->Sumw2();
	  diff_hist->Add(hB0Mass_1B0matched[0][0], -1);
	  hB0Mass_1B0->Draw() ;

	  if ( (hB0Mass_1B0matched[0][1] = (TH1F*) AlexisOut->Get("B0Mass_1B0matched_dauCh")) ) {
	    TH1F* diff_hist2 = (TH1F*) hB0Mass_1B0->Clone("diff_hist2"); diff_hist2->Sumw2();
	    diff_hist2->Add(hB0Mass_1B0matched[0][1], -1);
	    diff_hist2->SetFillColor( kYellow ) ; diff_hist2->Draw("hist same") ;
	    //diff_hist->SetFillColor( kRed ) ; diff_hist->Draw("hist same") ;
	    hB0Mass_1B0->SetLineColor( kBlack ) ; hB0Mass_1B0->SetLineWidth( 2 ) ; hB0Mass_1B0->Draw("same") ;
	    gPad->SaveAs( TString::Format("plots/%s/B0mass_1B0.png", dir.Data()) );

	  } else Warning("Terminate", "histogram \"B0Mass_1B0matched_dauCh\" not found");
	} else Warning("Terminate", "histogram \"B0Mass_1B0matched\" not found");
	
	TCanvas *Bpeak_C = new TCanvas("Bpeak_C","Bpeak_C", 800, 600) ;  

	Double_t xMin = hB0Mass_1B0->GetXaxis()->GetXmin();
	Double_t xMax = hB0Mass_1B0->GetXaxis()->GetXmax();
	Int_t nBins = hB0Mass_1B0->GetNbinsX();
	Double_t bin_width = (xMax-xMin)/nBins ;
  
	RooRealVar xVar("xVar", hB0Mass_1B0->GetXaxis()->GetTitle(), xMin, xMax) ;  xVar.setBins(nBins) ;  
	RooDataHist *RooB0Mass_1B0 = new RooDataHist(hB0Mass_1B0->GetName(), hB0Mass_1B0->GetTitle(), RooArgSet(xVar), Import(*hB0Mass_1B0, kFALSE)) ;
	RooConstVar nEntries("Entries", "Total number of entries", hB0Mass_1B0->Integral()) ;
	Double_t fitMax = xMax, fitMin = xMin ;

	RooRealVar mean("m", "mean", 5.28, 5.25, 5.31) ;
	RooRealVar sigma1("#sigma_{1}", "sigma1", 0.01, 0.001, 0.02);
	//RooRealVar sigma2("#sigma_{2}", "sigma2", 0.02, 0.01, 0.2);
	
	RooGaussian sigPDF("sigPDF", "Signal component", xVar, mean, sigma1); sigma1.SetNameTitle("#sigma","sigma");
	//RooGaussian sig1PDF("sig1PDF", "Signal component 1", xVar, mean, sigma1);
	//RooGaussian sig2PDF("sig2PDF", "Signal component 2", xVar, mean, sigma2);
	//RooRealVar ratio("sig1frac", "fraction of component 1 in signal", 0.1, 0, 0.2);
	//RooAddPdf sigPDF("sigPDF", "Signal component", RooArgList(sig1PDF, sig2PDF), ratio);
	
	RooRealVar nSig("N_{Sig}", "Number of signal candidates", 2.5e+3, 1e+2, 1e+5);

	RooRealVar c1("c_{1}", "c1", +0.02, -10, 10);
	RooRealVar c2("c_{2}", "c2", -0.001, -10, 10);
	RooChebychev bkgPDF("bkgPDF", "bkgPDF", xVar, RooArgSet(c1,c2));
	RooRealVar nBkg("nBkg", "Bkg component", 1e+3, 1e+2, 1e+4);     
  
	RooAddPdf *totalPDF = new RooAddPdf("totalPDF", "totalPDF", RooArgList(sigPDF, bkgPDF), RooArgList(nSig, nBkg));  // with 5 sigma
	cout <<"\n\n---------------  Fitting histo \"" <<RooB0Mass_1B0->GetTitle() <<"\"  ---------------\n" <<endl;
  
	totalPDF->fitTo(*RooB0Mass_1B0, Range(fitMin, fitMax), Extended(kTRUE));
  
	RooPlot* xframe = xVar.frame();
	//xframe->SetAxisRange(fitMin, xMax) ;
	xframe->SetTitle( RooB0Mass_1B0->GetTitle() ) ;
	xframe->SetYTitle( TString::Format("Candidates / %.1f MeV/c^{2}", bin_width*1000) ) ;
  
	RooB0Mass_1B0->plotOn(xframe);
	totalPDF->plotOn(xframe);
	totalPDF->plotOn(xframe, Components(RooArgSet(sigPDF)), LineColor(kRed));
	//totalPDF->plotOn(xframe, Components(RooArgSet(sig1PDF)), LineColor(kOrange));
	//totalPDF->plotOn(xframe, Components(RooArgSet(sig2PDF)), LineColor(kMagenta));

	totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig,mean,sigma1)), ShowConstants(kTRUE), Layout(0.55,0.9,0.9)) ;
	//totalPDF->paramOn(xframe, Parameters(RooArgSet(nSig,mean,sigma1,sigma2)), Layout(0.55,0.9,0.9)) ;

	xframe->SetTitleOffset(1.3, "XY");
	xframe->Draw();
	Bpeak_C->SaveAs( TString::Format("plots/%s/B0_mass_1B0_fit.png", dir.Data()) );

	// calculating purity
	xVar.setRange("signal", mean.getVal() -2*(sigma1.getVal()), mean.getVal() +2*(sigma1.getVal())) ;
	RooAbsReal* tot_sigRange = totalPDF->createIntegral(xVar, NormSet(xVar), Range("signal")) ;
	RooAbsReal* sig_sigRange = sigPDF.createIntegral(xVar, NormSet(xVar), Range("signal")) ;

	Double_t tot_sigRange_val = tot_sigRange->getVal() ; // compute integral
	cout <<"\nIntegral of total PDF in signal region = " << nEntries.getVal()*tot_sigRange_val << endl;
	cout <<"\nIntegral of signal PDF in signal region = " << nSig.getVal()*(sig_sigRange->getVal()) << endl;
	cout <<"\nBackground events in signal region = " << nEntries.getVal()*tot_sigRange_val << " - " << nSig.getVal()*(sig_sigRange->getVal()) << " = " << nEntries.getVal()*tot_sigRange_val - nSig.getVal()*(sig_sigRange->getVal()) << endl;
	
      } else Warning("Terminate", "histogram \"B0Mass_1B0\" not found");
      //
      
      // 2 B0 not twin
      if ( (hB0Mass_2B0notTwin = (TH1F*) AlexisOut->Get("B0Mass_2B0notTwin")) ) {
	if ( (hB0Mass_2B0notTwin_matched[0] = (TH1F*) AlexisOut->Get("B0Mass_2B0notTwin_matched")) ) {
	  TH1F* diff_hist = (TH1F*) hB0Mass_2B0notTwin->Clone("diff_hist"); diff_hist->Sumw2();
	  diff_hist->Add(hB0Mass_2B0notTwin_matched[0], -1);
	  hB0Mass_2B0notTwin->Draw() ;

	  if ( (hB0Mass_2B0notTwin_matched[1] = (TH1F*) AlexisOut->Get("B0Mass_2B0notTwin_matched_dauCh")) ) {
	    TH1F* diff_hist2 = (TH1F*) hB0Mass_2B0notTwin->Clone("diff_hist2"); diff_hist2->Sumw2();
	    diff_hist2->Add(hB0Mass_2B0notTwin_matched[1], -1);
	    diff_hist2->SetFillColor( kYellow ) ; diff_hist2->Draw("hist same") ;
	    //diff_hist->SetFillColor( kRed ) ; diff_hist->Draw("hist same") ;
	    hB0Mass_2B0notTwin->SetLineColor( kBlack ) ; hB0Mass_2B0notTwin->SetLineWidth( 2 ) ; hB0Mass_2B0notTwin->Draw("same") ;
	    gPad->SaveAs( TString::Format("plots/%s/B0mass_2B0notTwin.png", dir.Data()) );

	  } else Warning("Terminate", "histogram \"B0Mass_2B0notTwin_matched_dauCh\" not found");
	} else Warning("Terminate", "histogram \"B0Mass_2B0notTwin_matched\" not found");
      } else Warning("Terminate", "histogram \"B0Mass_2B0notTwin\" not found");
      //
      
      // 2 B0 twin
      if ( (hB0Mass_2B0twin = (TH1F*) AlexisOut->Get("B0Mass_2B0twin")) ) {
	if ( (hB0Mass_2B0twin_matched[0] = (TH1F*) AlexisOut->Get("B0Mass_2B0twin_matched")) ) {
	  TH1F* diff_hist = (TH1F*) hB0Mass_2B0twin->Clone("diff_hist"); diff_hist->Sumw2();
	  diff_hist->Add(hB0Mass_2B0twin_matched[0], -1);
	  hB0Mass_2B0twin->Draw() ;

	  if ( (hB0Mass_2B0twin_matched[1] = (TH1F*) AlexisOut->Get("B0Mass_2B0twin_matched_dauCh")) ) {
	    TH1F* diff_hist2 = (TH1F*) hB0Mass_2B0twin->Clone("diff_hist2"); diff_hist2->Sumw2();
	    diff_hist2->Add(hB0Mass_2B0twin_matched[1], -1);
	    diff_hist2->SetFillColor( kYellow ) ; diff_hist2->Draw("hist same") ;
	    //diff_hist->SetFillColor( kRed ) ; diff_hist->Draw("hist same") ;
	    hB0Mass_2B0twin->SetLineColor( kBlack ) ; hB0Mass_2B0twin->SetLineWidth( 2 ) ; hB0Mass_2B0twin->Draw("same") ;
	    gPad->SaveAs( TString::Format("plots/%s/B0mass_2B0twin.png", dir.Data()) );

	  } else Warning("Terminate", "histogram \"B0Mass_2B0twin_matched_dauCh\" not found");
	} else Warning("Terminate", "histogram \"B0Mass_2B0twin_matched\" not found");
      } else Warning("Terminate", "histogram \"B0Mass_2B0twin\" not found");
      //
      
      // hourglass
      cout <<"\nHourglass plot:" <<endl ;

      if ( ( hDeltaB0Mass_2B0twin[0] = (TH2F*) AlexisOut->Get("DeltaB0Mass_2B0twin")) ) {
	hDeltaB0Mass_2B0twin[0]->Draw("colz") ;
	gPad->SaveAs( TString::Format("plots/%s/Hourglass.png", dir.Data()) ) ;
      } else Warning("Terminate", "histogram \"DeltaB0Mass_2B0twin\" not found");

      if ( ( hDeltaB0Mass_2B0twin[1] = (TH2F*) AlexisOut->Get("DeltaB0Mass_2B0twin_signalWin")) ) {
	hDeltaB0Mass_2B0twin[1]->Draw("colz") ;
	gPad->SaveAs( TString::Format("plots/%s/Hourglass_signalWin.png", dir.Data()) ) ;
      } else Warning("Terminate", "histogram \"DeltaB0Mass_2B0twin_signalWin\" not found");

      if ( ( hB0Mass_2B0twin_hourglass[1][1] = (TH1F*) AlexisOut->Get("B0Mass_inside_matchedTwin")) ) {
	hB0Mass_2B0twin_hourglass[1][1]->SetLineWidth(2) ;
	hB0Mass_2B0twin_hourglass[1][1]->SetLineColor( kBlack ) ;
	hB0Mass_2B0twin_hourglass[1][1]->SetStats(0) ;
	hB0Mass_2B0twin_hourglass[1][1]->Draw() ;

	if ( ( hB0Mass_2B0twin_hourglass[1][0] = (TH1F*) AlexisOut->Get("B0Mass_inside_notMatchedTwin")) ) {
	  //hB0Mass_2B0twin_hourglass[1][0]->SetLineColor(kYellow) ; //hB0Mass_2B0twin_hourglass[1][0]->SetLineWidth(2) ;
	  hB0Mass_2B0twin_hourglass[1][0]->SetFillColor(kYellow) ; 
	  hB0Mass_2B0twin_hourglass[1][0]->Draw("same") ; 

	  if ( ( hB0Mass_2B0twin_hourglass[0][1] = (TH1F*) AlexisOut->Get("B0Mass_outside_matchedTwin")) ) {
	    hB0Mass_2B0twin_hourglass[0][1]->SetLineWidth(2) ;
	    hB0Mass_2B0twin_hourglass[0][1]->SetLineColor(kBlue) ;
	    hB0Mass_2B0twin_hourglass[0][1]->Draw("same") ;

	    if ( ( hB0Mass_2B0twin_hourglass[0][0] = (TH1F*) AlexisOut->Get("B0Mass_outside_notMatchedTwin")) ) {
	      hB0Mass_2B0twin_hourglass[0][0]->SetLineWidth(2) ;
	      hB0Mass_2B0twin_hourglass[0][0]->SetLineColor(kRed) ;
	      hB0Mass_2B0twin_hourglass[0][0]->Draw("same") ;

	    } else Warning("Terminate", "histogram \"B0Mass_outside_notMatchedTwin\" not found");
	  } else Warning("Terminate", "histogram \"B0Mass_outside_matchedTwin\" not found");
	} else Warning("Terminate", "histogram \"B0Mass_inside_notMatchedTwin\" not found");

	gPad->BuildLegend() ; hB0Mass_2B0twin_hourglass[1][1]->Draw("same") ;
	gPad->RedrawAxis() ;
	gPad->SaveAs( TString::Format("plots/%s/B0_mass_2B0twin_hourglass.png", dir.Data()) ) ;
	
	gPad->SetLogy() ;
	gPad->SaveAs( TString::Format("plots/%s/B0_mass_2B0twin_hourglass_logy.png", dir.Data()) ) ;

      } else Warning("Terminate", "histogram \"B0Mass_inside_matchedTwin\" not found");
      */
      
      cout <<"\nDalitz plot:" <<endl ;
      //
      TString nameGen_squared = "cos_psi2S_helicityAngle_vs_psi2SPi_BT_gen";
      TH2F* dalitzGen_squared = 0 ;
      if ( (dalitzGen_squared = (TH2F*) AlexisOut->Get( nameGen_squared.Data() )) ) {
	gStyle->SetOptStat( 10 ) ;
	
	dalitzGen_squared->Draw("colz") ;
	gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), nameGen_squared.Data()) );
      } else Warning("Terminate", TString::Format("histogram %s not found", nameGen_squared.Data()));
      //
      nameGen_squared = "cos_psi2S_helicityAngle_vs_KPi_BT_gen";
      if ( (dalitzGen_squared = (TH2F*) AlexisOut->Get( nameGen_squared.Data() )) ) {
	gStyle->SetOptStat( 10 ) ;

	dalitzGen_squared->Draw("colz") ;
	gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), nameGen_squared.Data()) );
      } else Warning("Terminate", TString::Format("histogram %s not found", nameGen_squared.Data()));
      //
      //TString nameGen = "gen_psi2SPi_vs_KPi_dalitz";
      TString nameGen = "gen_psi2SPi_vs_KPi_dalitz_BT";
      //nameGen = "gen_KPi_vs_psi2SK_dalitz_BT";
      TH2F* dalitzGen = 0 ;
      if ( (dalitzGen = (TH2F*) AlexisOut->Get( nameGen.Data() )) ) {

	gStyle->SetOptStat( 10 ) ;
	
	//dalitzGen->Draw("lego") ;
	dalitzGen->Draw("colz") ;

	// Dalitz contour
	Int_t m = 19800;
	Double_t x[19800], m23_max[19800], m23_min[19800];
	Double_t E2[19800], E3[19800];
	
	Double_t m_mother = B0_mass; 
	//Double_t m_dau1 = kaonCh_mass, m_dau2 = pionCh_mass, m_dau3 = psi2S_mass ; 
	Double_t m_dau1 = kaonCh_mass, m_dau2 = pionCh_mass, m_dau3 = psi2S_mass ; 

	Double_t m12_min = (m_dau1+m_dau2)*(m_dau1+m_dau2);
	Double_t m12_max = (m_mother-m_dau3)*(m_mother-m_dau3);
	Double_t step = (m12_max - m12_min)/(m-1);
	
	x[0] = m12_min + 0.0001;
	
	for (Int_t k=1; k<m; k++ ) 
	  x[k] = x[k-1] + step;
	
	Int_t n = 19799;
	for (Int_t i=0; i<n; i++) {
	  E2[i] = (x[i] - m_dau1*m_dau1 + m_dau2*m_dau2)/(2*sqrt(x[i]));
	  E3[i] = (m_mother*m_mother - x[i] - m_dau3*m_dau3)/(2*sqrt(x[i]));
	  m23_min[i] = (E2[i]+E3[i])*(E2[i]+E3[i]) - TMath::Power((sqrt(E2[i]*E2[i] - m_dau2*m_dau2) + sqrt(E3[i]*E3[i] - m_dau3*m_dau3)),2);
	  m23_max[i] = (E2[i]+E3[i])*(E2[i]+E3[i]) - TMath::Power((sqrt(E2[i]*E2[i] - m_dau2*m_dau2) - sqrt(E3[i]*E3[i] - m_dau3*m_dau3)),2);
	}
	
	TGraph *cont_up = new TGraph(n,x,m23_min); cont_up->SetLineWidth(3);
	TGraph *cont_down = new TGraph(n,x,m23_max); cont_down->SetLineWidth(3);
	cont_up->Draw("lsame"); 
	cont_down->Draw("lsame"); 

	gPad->SaveAs( TString::Format("plots/%s/dalitz_gen.png", dir.Data()) );
	//dalitzGen->Rebin2D();
	//gPad->SaveAs( TString::Format("plots/%s/dalitz_gen_rebinned.png", dir.Data()) );

	//TString nameReco = "reco_psi2SPi_vs_KPi_dalitz" ;
	TString nameReco = "psi2SPi_vs_KPi_dalitz_1B0" ;
	TH2F* dalitzReco = 0 ;
	//if ( (psi2SPiSqMass_vs_KPiSqMass_AC = (TH2F*) AlexisOut->Get( nameReco.Data() )) ) {
	if ( (dalitzReco = (TH2F*) AlexisOut->Get( nameReco.Data() )) ) {
	  dalitzReco->Draw("colz") ;
	  cont_up->Draw("lsame"); 
	  cont_down->Draw("lsame"); 

	  gPad->SaveAs( TString::Format("plots/%s/dalitz_reco.png", dir.Data()) );
	  //dalitzReco->Rebin2D();
	  //gPad->SaveAs( TString::Format("plots/%s/dalitz_reco_rebinned.png", dir.Data()) );

	  cout <<"\nCreating efficiency Dalitz" <<endl ;
	  TH2D *hAbsEff, *hAbsEffErr, *hRelEff;
	  TH1D *hAbsEffDistr, *hGenDistr, *hRecoDistr;
	  
	  Float_t xLow = dalitzGen->GetXaxis()->GetXmin() ;
	  Float_t xHigh = dalitzGen->GetXaxis()->GetXmax() ;
	  Int_t nX = dalitzGen->GetNbinsX() ;
	  Float_t yLow = dalitzGen->GetYaxis()->GetXmin() ;
	  Float_t yHigh = dalitzGen->GetYaxis()->GetXmax() ;
	  Int_t nY = dalitzGen->GetNbinsY() ;
	  hAbsEff = new TH2D ("AbsEff", "AbsEff", nX, xLow, xHigh, nY, yLow, yHigh);
	  hAbsEff->SetXTitle( dalitzGen->GetXaxis()->GetTitle() );
	  hAbsEff->SetYTitle( dalitzGen->GetYaxis()->GetTitle() );
	  hRelEff = new TH2D ("RelEff", "RelEff", nX, xLow, xHigh, nY, yLow, yHigh);
	  hRelEff->SetXTitle( dalitzGen->GetXaxis()->GetTitle() );
	  hRelEff->SetYTitle( dalitzGen->GetYaxis()->GetTitle() );
	  hAbsEffErr = new TH2D ("AbsEffErr", "AbsEffErr", nX, xLow, xHigh, nY, yLow, yHigh);
	  hAbsEffErr->SetXTitle( dalitzGen->GetXaxis()->GetTitle() );
	  hAbsEffErr->SetYTitle( dalitzGen->GetYaxis()->GetTitle() );
	  hAbsEffDistr = new TH1D ("hAbsEffDistr","Efficiency distribution", 100, -0.001, 0.006);
	  hGenDistr = new TH1D ("hGenDistr","Bins' Gen entries distribution", 100, 0., 25000);
	  hRecoDistr = new TH1D ("hRecoDistr","Bins' Reco entries distribution", 110, 0., 55.);
	 	  
	  Int_t MyPalette[100];
	  Double_t red[]   = { 1.0, 0.0, 1.0, 1.0 };
	  Double_t green[] = { 1.0, 0.0, 0.0, 1.0 };
	  Double_t blue[]  = { 1.0, 1.0, 0.0, 0.0 };
	  Double_t stop[] = { 0.0, 0.5, 0.9, 1.0 };
	  Int_t FI = TColor::CreateGradientColorTable(4, stop, red, green, blue, 100);
	  for (int i=0; i<100; i++)
	    MyPalette[i] = FI+i;
	  
	  gStyle->SetPalette(100, MyPalette);
	  
	  // 32 righe & 40 colonne
	  Double_t genBinCont[nX][nY], recoBinCont[nX][nY];
	  Double_t absEff[nX][nY], absEffErr[nX][nY];
	  Double_t absEffSingleValue, absEffErrSingleValue;
	  Double_t genSingleValue, recoSingleValue;
	  //
	  Int_t count=0;
	  Double_t sumAbsEff = 0, sumAbsEffErr = 0.;
	  //
	  for (Int_t i=0; i<nX; i++)
	    {
	      for (Int_t j=0; j<nY; j++)
		{
		  genBinCont[i][j] = dalitzGen->GetBinContent(i,j);
		  recoBinCont[i][j] = dalitzReco->GetBinContent(i,j);
		  //cout << "riga j=" << j << " e colonna i=" << i << " -> no. candidati generati = " << genBinCont[j][i] << " & no. ricostruiti = " << recoBinCont[j][i] << endl;
		  //
		  genSingleValue = genBinCont[i][j];
		  recoSingleValue = recoBinCont[i][j];
		  hGenDistr->Fill(genSingleValue);
		  hRecoDistr->Fill(recoSingleValue);
		  //
		  //if ( genBinCont[i][j] > 500. && recoBinCont[i][j] <= genBinCont[i][j] ) // fix 500 with sth not hardcoded
		  if ( recoBinCont[i][j] <= genBinCont[i][j] ) 
		    {
		      absEff[i][j] = recoBinCont[i][j] / genBinCont[i][j];
		      absEffSingleValue = absEff[i][j];
		      hAbsEffDistr->Fill(absEffSingleValue);
		      absEffErr[i][j] = (1 / genBinCont[i][j] ) * (sqrt( genBinCont[i][j] * absEff[i][j] * (1 - absEff[i][j]) ) );
		      absEffErrSingleValue = absEffErr[i][j];
		      //
		      // calculate weighted mean: two possible intervals: ]0,0.00225[ or ]0,0.005[
		      //if (absEffSingleValue > 0. && absEffSingleValue < 0.01) //0.005
		      if (absEffSingleValue > 0.) //0.005
			{
			  sumAbsEff = sumAbsEff + (absEffSingleValue/(absEffErrSingleValue * absEffErrSingleValue)) ;
			  sumAbsEffErr = sumAbsEffErr + (1./(absEffErrSingleValue * absEffErrSingleValue));
			  count++ ;
			}
		      //
		    } else {
		    absEff[i][j] = 0.;
		    absEffErr[i][j] = 0.;
		    absEffSingleValue = -0.0005;
		    hAbsEffDistr->Fill(absEffSingleValue);
		  }
		  //cout << "raw-j=" << j << " e col-i=" << i << " -> no. gen-cand = " << genBinCont[j][i] << " & no. reco-cand = " << recoBinCont[j][i] << " & eff= " << absEff[j][i] << endl;
		  //
		  hAbsEff->SetBinContent(i,j,absEff[i][j]);
		  hAbsEffErr->SetBinContent(i,j,absEffErr[i][j]);
		  //
		} // for (Int_t j=0; j<nY; j++)
	    } // for (Int_t i=0; i<nX; i++)
	  //
	  cout <<"count = " <<count << " with weighted mean = " << sumAbsEff/sumAbsEffErr << " and weighted mean uncertainty = " << 1./sqrt(sumAbsEffErr) << endl; // about 10-3 with error 10-5
	  //
	  hAbsEff->Draw("colz"); hAbsEff->SetStats(kFALSE);
	  cont_up->Draw("lsame"); 
	  cont_down->Draw("lsame"); 
	  gPad->SaveAs( TString::Format("plots/%s/dalitz_absEff.png", dir.Data()) );

	  hAbsEffErr->Draw("colz"); hAbsEffErr->SetStats(kFALSE);
	  cont_up->Draw("lsame"); 
	  cont_down->Draw("lsame"); 
	  gPad->SaveAs( TString::Format("plots/%s/dalitz_absEffErr.png", dir.Data()) );

	  //////////////////////////////////////// RELATIVE EFFICIENCY CALCULATION
	  //   
	  Double_t relEff[nX][nY];
	  //Double_t relEffErr[nX][nY];
	  //
	  Double_t weightedMeanEff = sumAbsEff/sumAbsEffErr;
	  // Double_t weightedMeanEffErr = 1./sqrt(sumAbsEffErr);
	  //
	  for (Int_t i=0; i<nX; i++)
	    {
	      for (Int_t j=0; j<nY; j++)
		{
		  //genBinCont[i][j] = dalitzGen->GetBinContent(i,j);
		  //recoBinCont[i][j] = dalitzReco->GetBinContent(i,j);	   
		  //if ( genBinCont[i][j] > 500. && recoBinCont[i][j] <= genBinCont[i][j] ) 
		  if ( recoBinCont[i][j] <= genBinCont[i][j] )
		    {
		      //absEff[i][j] = recoBinCont[i][j] / genBinCont[i][j];
		      //if (absEff[i][j] > 0. && absEff[i][j] < 0.05) // 0. && 0.005
		      if (absEff[i][j] > 0.)
			{
			  relEff[i][j] = absEff[i][j]/weightedMeanEff;
			}
		      else 
			relEff[i][j] = 0. ;
		    }
		  else
		    relEff[i][j] = 0. ;
		  hRelEff->SetBinContent(i,j,relEff[i][j]);
		}
	    }
	  //
	  hRelEff->Draw("colz"); hRelEff->SetStats(kFALSE);
	  hRelEff->SetMinimum(0.); hRelEff->SetMaximum(2.2);
	  cont_up->Draw("lsame"); cont_down->Draw("lsame"); 
	  gPad->SaveAs( TString::Format("plots/%s/dalitz_relEff.png", dir.Data()) );

	  hAbsEffDistr->Draw() ;
	} else Warning("Terminate", TString::Format("histogram \"%s\" not found", nameReco.Data()));

      } else Warning("Terminate", TString::Format("histogram \"%s\" not found", nameGen.Data()));
     
      cout <<"\nAngles plot:" <<endl ;
      // psi helicity angle
      TH1F* helicityAngleGen = 0, *cosHelicityAngleGen = 0 ;
      TString name = "psi2S_helicityAngle_BT_gen";
      TString nameCos = "cos_psi2S_helicityAngle_BT_gen";
      if ( (helicityAngleGen = (TH1F*) AlexisOut->Get( name.Data() )) ) {
	gStyle->SetOptStat( 10 ) ;
	helicityAngleGen->Draw() ;
	gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), helicityAngleGen->GetName()) );

	if ( (cosHelicityAngleGen = (TH1F*) AlexisOut->Get( nameCos.Data() )) ) {
	  cosHelicityAngleGen->Draw() ;
	  gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), cosHelicityAngleGen->GetName()) );
	  
	  TH1F* helicityAngleReco = 0, *cosHelicityAngleReco = 0 ;
	  TString nameReco = "psi2S_helicityAngle" ;
	  TString nameRecoCos = "cos_psi2S_helicityAngle" ;
	  if ( (helicityAngleReco = (TH1F*) AlexisOut->Get( nameReco.Data() )) ) {
	    helicityAngleReco->Draw() ;
	    gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), helicityAngleReco->GetName()) );

	    if ( (cosHelicityAngleReco = (TH1F*) AlexisOut->Get( nameRecoCos.Data() )) ) {
	      cosHelicityAngleReco->Draw() ;
	      gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), cosHelicityAngleReco->GetName()) );

	      TH1F* cosHelicityAngleEff = (TH1F*) cosHelicityAngleReco->Clone("cos_psi2S_helicityAngle_Eff"); cosHelicityAngleEff->SetTitle("efficiency of #psi' helicity angle cosine") ;
	      cosHelicityAngleEff->Divide( cosHelicityAngleReco, cosHelicityAngleGen, 1,1,"B") ;
	      cosHelicityAngleEff->Draw() ;
	      gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), cosHelicityAngleEff->GetName()) );

	    } else Warning("Terminate", TString::Format("histogram \"%s\" not found", nameRecoCos.Data()));
	  } else Warning("Terminate", TString::Format("histogram \"%s\" not found", nameReco.Data()));
	
	} else Warning("Terminate", TString::Format("histogram \"%s\" not found", nameCos.Data()));
      } else Warning("Terminate", TString::Format("histogram \"%s\" not found", name.Data()));
      //
      // planes angle
      TH1F* planesAngleGen = 0 ;
      name = "planesAngle_BT_gen";
      if ( (planesAngleGen = (TH1F*) AlexisOut->Get( name.Data() )) ) {

	gStyle->SetOptStat( 10 ) ;

	planesAngleGen->Draw() ;
	gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), planesAngleGen->GetName()) );

	TH1F* planesAngleReco = 0 ;
	TString nameReco = "planesAngle" ;
	if ( (planesAngleReco = (TH1F*) AlexisOut->Get( nameReco.Data() )) ) {

	  planesAngleReco->Draw() ;
	  gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), planesAngleReco->GetName()) );

	  TH1F* planesAngleEff = (TH1F*) planesAngleReco->Clone("planesAngle_Eff"); planesAngleEff->SetTitle("efficiency of #phi angle") ;
	  planesAngleEff->Divide( planesAngleReco, planesAngleGen, 1,1,"B") ;
	  planesAngleEff->Draw() ;
	  gPad->SaveAs( TString::Format("plots/%s/%s.png", dir.Data(), planesAngleEff->GetName()) );

	} else Warning("Terminate", TString::Format("histogram \"%s\" not found", nameReco.Data()));
	
      } else Warning("Terminate", TString::Format("histogram \"%s\" not found", name.Data()));
      
    } else Error("Terminate", "could not open file: %s", outputFile.Data()) ;
    
  }
  else Warning("Terminate", "file \"26Jan_histos.root\" not found");

  cout <<endl ;
  
  // Primary vertex plots
  /*
  cout <<"\nPrimary vertex plots:" <<endl ;
  gSystem->mkdir(TString::Format("plots/%s/priVtx", dir.Data()), true);

  if ( (OutFile = dynamic_cast<TProofOutputFile*>(fOutput->FindObject("officialMC_Bd2Psi2SKpi.root"))) ) {
    TString outputFile(OutFile->GetOutputFileName());
    TString outputName(OutFile->GetName()); outputName += ".root";
    Printf("\nOutput file found: %s", outputFile.Data());   
    // Open the file
    fOut = TFile::Open(outputFile);

    gStyle->SetOptStat( 10 ) ;
    // primary vertex
    //TCanvas *priVtx_C = new TCanvas("priVtx_C","priVtx_C", 800, 600) ; // useless
    TString type[7] = {"","_B0","_B0CosAlpha","_B0CosAlpha3D","_B0Less","_B0LessB0CosAlpha","_B0LessB0CosAlpha3D"} ;
    gPad->SetLogy(0) ;
    for (Int_t i=0; i<7; i++) {
      TString name = TString::Format("priVtx%s_z", type[i].Data()) ;
      if ( (priVtx_z[i] = (TH1F*) fOut->Get( name.Data() )) ) {
	priVtx_z[i]->Draw() ;
	gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtx_z[i]->GetName()));
      } else Warning("Terminate", TString::Format("histogram \"%s\" not found", name.Data() ));
    }

    TString xyz[3] = {"x","y","z"}, XYZ[3] = {"X","Y","Z"} ;
    TString zoom[2] = {"","_zoom"} ;
    TString dim[2] = {"2D","3D"} ;
    for (Int_t j=0; j<2; j++) { 
      for (Int_t i=0; i<3; i++) {
	// priVtx_delta
	TString name = TString::Format("priVtx_delta%s%s", XYZ[i].Data(), zoom[j].Data()) ;
	if ( (priVtx_delta[j][i] = (TH1F*) fOut->Get( name.Data() )) ) {
	  priVtx_delta[j][i]->Draw() ; gPad->SetLogy() ;
	  if (i==2 && j==1)
	    priVtx_delta[j][i]->GetXaxis()->SetRangeUser(-0.0002, 0.0002) ;
	  gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtx_delta[j][i]->GetName()));
	} else Warning("Terminate", TString::Format("histogram \"%s\" not found", name.Data()));

	// priVtx_B0_delta
	name = TString::Format("priVtx_B0_delta%s%s", XYZ[i].Data(), zoom[j].Data()) ;
	if ( (priVtx_B0_delta[j][i] = (TH1F*) fOut->Get( name.Data() )) ) {
	  priVtx_B0_delta[j][i]->Draw() ; gPad->SetLogy() ;
	  gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtx_B0_delta[j][i]->GetName()));
	} else Warning("Terminate", TString::Format("histogram \"%s\" not found", name.Data()));

	// priVtx_B0Less_delta
	name = TString::Format("priVtxB0Less_delta%s%s", XYZ[i].Data(), zoom[j].Data()) ;
	if ( (priVtxB0Less_delta[j][i] = (TH1F*) fOut->Get( name.Data() )) ) {
	  priVtxB0Less_delta[j][i]->Draw() ; gPad->SetLogy() ;
	  gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtxB0Less_delta[j][i]->GetName()));
	} else Warning("Terminate", TString::Format("histogram \"%s\" not found", name.Data()));

	// dim = 2D, 3D
	for (Int_t k=0; k<2; k++) { 
	  // priVtx_B0CosAlpha_delta
	  name = TString::Format("priVtx_B0CosAlpha%s_delta%s%s", dim[k].Data(), XYZ[i].Data(), zoom[j].Data()) ;
	  if ( (priVtx_B0CosAlpha_delta[k][j][i] = (TH1F*) fOut->Get( name.Data() )) ) {
	    priVtx_B0CosAlpha_delta[k][j][i]->Draw() ; gPad->SetLogy() ;
	    gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtx_B0CosAlpha_delta[k][j][i]->GetName()));
	  } else Warning("Terminate", TString::Format("histogram \"%s\" not found", name.Data()));
	  
	  // priVtx_B0Less_B0CosAlpha_delta
	  name = TString::Format("priVtxB0Less_B0CosAlpha%s_delta%s%s", dim[k].Data(), XYZ[i].Data(), zoom[j].Data()) ;
	  if ( (priVtxB0Less_B0CosAlpha_delta[k][j][i] = (TH1F*) fOut->Get( name.Data() )) ) {
	    priVtxB0Less_B0CosAlpha_delta[k][j][i]->Draw() ; gPad->SetLogy() ;
	    gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtxB0Less_B0CosAlpha_delta[k][j][i]->GetName()));
	  } else Warning("Terminate", TString::Format("histogram \"%s\" not found", name.Data()));
	
	  // same
	  priVtx_B0_delta[j][i]->Draw() ;
	  priVtx_B0CosAlpha_delta[k][j][i]->Draw("same") ;
	  priVtxB0Less_delta[j][i]->Draw("same") ;
	  priVtxB0Less_B0CosAlpha_delta[k][j][i]->Draw("same") ;
	  gPad->SaveAs( TString::Format("plots/%s/priVtx/priVtx_all%s_delta%s%s.png", dir.Data(), dim[k].Data(), XYZ[i].Data(), zoom[j].Data()));
	}

      }
    }

    TH1F* priVtx_diff[1][3] = {{0,0,0}}, *priVtx_diffZ = 0, *priVtx3D_diff[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}} ;
    Int_t isZoom = 1 ;
    Int_t is3D = 1 ;
    TString diff[3] = {"B0","B0Less","B0Less_B0CosAlpha"} ;

    for (Int_t i=0; i<3; i++) {
      priVtx_diff[0][i] = (TH1F*) priVtxB0Less_delta[isZoom][i]->Clone( TString::Format("priVtx_delta%s_diff", XYZ[i].Data()) ) ; priVtx_diff[0][i]->Sumw2();
      priVtx_diff[0][i]->Add(priVtx_B0_delta[isZoom][i], -1) ;
      priVtx_diff[0][i]->SetTitle( TString::Format("#Delta%s(B^{0}_less PV) - #Delta%s(standard B^{0} PV);#Delta%s [cm]", xyz[i].Data(), xyz[i].Data(), xyz[i].Data()) ) ;
      priVtx_diff[0][i]->SetFillColor( priVtxB0Less_delta[isZoom][i]->GetLineColor() ) ;
      priVtx_diff[0][i]->Draw("hist") ;
      gPad->SetLogy(0) ;
      gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtx_diff[0][i]->GetName()));

      // 3D
      for (Int_t j=0; j<3; j++) {
	priVtx3D_diff[j][i] = (TH1F*) priVtx_B0CosAlpha_delta[is3D][isZoom][i]->Clone( TString::Format("priVtx_B0CosAlpha%s_delta%s_diff_%s", dim[is3D].Data(), XYZ[i].Data(), diff[j].Data()) ) ; priVtx3D_diff[j][i]->Sumw2();
	if (j == 0) {
	  priVtx3D_diff[j][i]->Add(priVtx_B0_delta[isZoom][i], -1) ;
	  priVtx3D_diff[j][i]->SetTitle( TString::Format("#Delta%s(B^{0}_cos(#alpha_{%s}) PV) - #Delta%s(standard B^{0} PV);#Delta%s [cm]", xyz[i].Data(), dim[is3D].Data(), xyz[i].Data(), xyz[i].Data()) ) ;
	  priVtx3D_diff[j][i]->SetFillColor( priVtx_B0_delta[isZoom][i]->GetLineColor() ) ;
	} else if (j == 1) {
	  priVtx3D_diff[j][i]->Add(priVtxB0Less_delta[isZoom][i], -1) ;
	  priVtx3D_diff[j][i]->SetTitle( TString::Format("#Delta%s(B^{0}_cos(#alpha_{%s}) PV) - #Delta%s(B^{0}_less PV);#Delta%s [cm]", xyz[i].Data(), dim[is3D].Data(), xyz[i].Data(), xyz[i].Data()) ) ;
	  priVtx3D_diff[j][i]->SetFillColor( priVtxB0Less_delta[isZoom][i]->GetLineColor() ) ;
	} else if (j == 2) {
	  priVtx3D_diff[j][i]->Add(priVtxB0Less_B0CosAlpha_delta[is3D-1][isZoom][i], -1) ;
	  priVtx3D_diff[j][i]->SetTitle( TString::Format("#Delta%s(B^{0}_cos(#alpha_{%s}) PV) - #Delta%s(B^{0}_cos(#alpha_{%s}) B^{0}_less PV);#Delta%s [cm]", xyz[i].Data(), dim[is3D].Data(), xyz[i].Data(), dim[is3D-1].Data(), xyz[i].Data()) ) ;
	  priVtx3D_diff[j][i]->SetFillColor( priVtxB0Less_B0CosAlpha_delta[is3D-1][isZoom][i]->GetLineColor() ) ;
	}	
	priVtx3D_diff[j][i]->Draw("hist") ;
	gPad->SetLogy(0) ;
	gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtx3D_diff[j][i]->GetName()));
      }    
    }
    isZoom = 0 ;
    priVtx_diffZ = (TH1F*) priVtxB0Less_delta[isZoom][2]->Clone( TString::Format("priVtx_delta%s_diff_noZoom", XYZ[2].Data()) ) ; priVtx_diffZ->Sumw2() ;
    priVtx_diffZ->Add(priVtx_B0_delta[isZoom][2], -1) ;
    priVtx_diffZ->SetTitle( TString::Format("#Delta%s(B^{0}_less PV) - #Delta%s(standard B^{0} PV);#Delta%s [cm]", xyz[2].Data(), xyz[2].Data(), xyz[2].Data()) ) ;
    priVtx_diffZ->SetFillColor( priVtxB0Less_delta[isZoom][2]->GetLineColor() ) ;
    priVtx_diffZ->Draw("hist") ;
    gPad->SaveAs( TString::Format("plots/%s/priVtx/%s.png", dir.Data(), priVtx_diffZ->GetName()) );

  } else Warning("Terminate", "file \"officialMC_Bd2Psi2SKpi.root\" not found");
  */
  //gROOT->ProcessLine(".x b_peak/b_peakFit.C++");
}


int psiPrimePiK_MC::getXVarType(const double jpsivar, vector<double> vectbin) {
  for(unsigned int i=0;i<vectbin.size()-1;i++) {
    if(jpsivar > vectbin[i] && jpsivar < vectbin[i+1]) return i+1;
  }
  
  return -999;
}

//bool psiPrimePiK_MC::isMuonInAccept(float eta, float pt, float p) {
bool psiPrimePiK_MC::isMuonInAccept(float eta, float pt) {
  // *USE* muon kinematical cuts (eta dependent momentum / pT cuts )
  return (
	  ( (eta < 1.2) && (pt > 4.0) ) ||
	  ( (eta > 1.2) && (eta < 2.4) && (pt > 3.3) )
	  );
  
}

float psiPrimePiK_MC::GetWeightPtShape( float pt, TH1F* histo ) {
  
  float weight = 1.;
  for (int i=1; i<=histo->GetNbinsX() ;i++){
    float lowbin = histo->GetBinLowEdge(i);
    float upbin = lowbin + histo->GetBinWidth(i);
    //cout<<" pt="<<pt<<" lowbin="<<lowbin<<" upbin="<<upbin<<endl;
    if (  pt >= lowbin && pt<upbin ){
      weight = histo->GetBinContent(i);
      break;
    } 
  }
  return weight;
}

double psiPrimePiK_MC::getMumuHelicityAngle(TLorentzVector B0, TLorentzVector Kstar, TLorentzVector mumu, TLorentzVector muM) {

  TLorentzVector mumuBoosted = mumu, muMBoosted = muM, KstarBoosted = Kstar ;
  mumuBoosted.Boost( -B0.BoostVector() ) ;
  muMBoosted.Boost( -B0.BoostVector() ) ; muMBoosted.Boost( -mumuBoosted.BoostVector() ) ;
  KstarBoosted.Boost( -B0.BoostVector() ) ; KstarBoosted.Boost( -mumuBoosted.BoostVector() ) ;

  //return muMBoosted.Angle( KstarBoosted.Vect() );    
  return muMBoosted.Angle( mumuBoosted.Vect() );    
}

double psiPrimePiK_MC::getBelleMumuHelicityAngle(TLorentzVector B0, TLorentzVector Kstar, TLorentzVector mumu, TLorentzVector muM) {
		    
  TLorentzVector mumuBoosted = mumu, muMBoosted = muM, KstarBoosted = Kstar ;
  mumuBoosted.Boost( -B0.BoostVector() ) ;
  muMBoosted.Boost( -B0.BoostVector() ) ; muMBoosted.Boost( -mumuBoosted.BoostVector() ) ;
  TVector3 muM3 = muMBoosted.Vect() ;
  KstarBoosted.Boost( -B0.BoostVector() ) ; KstarBoosted.Boost( -mumuBoosted.BoostVector() ) ;
  TVector3 Kstar3 = KstarBoosted.Vect() ; 

  return TMath::ACos( -(muM3.Dot(Kstar3)) / (muM3.Mag()*Kstar3.Mag()) ) ;    
}

double psiPrimePiK_MC::getAleMumuHelicityAngle(TLorentzVector mumu, TLorentzVector muM) {
  
  //INPUT:
  // JPsi 4-momentum and muon 4-momentum
  //OUTPUT:
  // muon cos(theta) in helicity frame
  // muon phi in helicity frame
  //
  
  TLorentzVector mumuBoosted = mumu, muMBoosted = muM ;
  //JPsi momentum in lab frame
  TVector3 bMumu = mumuBoosted.BoostVector() ;
  TLorentzVector beam1(0., 0.,  4000., 4000.) ; // beam momentum in lab frame
  TLorentzVector beam2(0., 0., -4000., 4000.) ; // beam momentum in lab frame
  beam1.Boost( -bMumu ) ; // beam momentum in JPsi rest frame
  beam2.Boost( -bMumu ) ; // beam momentum in JPsi rest frame
  TVector3 beam1_dir = beam1.Vect().Unit();  // beam direction in JPsi rest frame
  TVector3 beam2_dir = beam2.Vect().Unit();  // beam direction in JPsi rest frame

  TVector3 Y = beam1_dir.Cross( beam2_dir ).Unit(); // the production plane normal
  TVector3 Z = mumuBoosted.Vect().Unit();         // JPsi direction in lab frame
  TVector3 X = Y.Cross(Z).Unit();         // completes the right-handed coordinate

  muMBoosted.Boost( -bMumu );     // muon- momentum in JPsi rest frame

  TRotation rotation;
  rotation.SetToIdentity();
  rotation.RotateAxes(X,Y,Z);
  rotation.Invert(); // transforms coordinates from the "xyz" frame to the new frame
  TVector3 mu_rotated = muMBoosted.Vect();
  mu_rotated.Transform(rotation);
  
  return mu_rotated.CosTheta() ;
  //return TMath::ACos( mu_rotated.CosTheta() );
  //muphi_hx = mu_rotated.Phi()>0. ? mu_rotated.Phi() : mu_rotated.Phi()+2.*TMath::Pi(); // radians
}

double psiPrimePiK_MC::getPlanesAngle(TLorentzVector B0, TLorentzVector K, TLorentzVector Pi, TLorentzVector muM, TLorentzVector muP) {

  TLorentzVector KBoosted = K, PiBoosted = Pi, muMBoosted = muM, muPBoosted = muP ;
  muPBoosted.Boost( -B0.BoostVector() ) ;
  muMBoosted.Boost( -B0.BoostVector() ) ;
  KBoosted.Boost( -B0.BoostVector() ) ;
  PiBoosted.Boost( -B0.BoostVector() ) ;
  TVector3 MuMuPlane_normal = muPBoosted.Vect().Cross( muMBoosted.Vect() ) ;
  TVector3 KPiPlane_normal = KBoosted.Vect().Cross( PiBoosted.Vect() ) ;
  Double_t angle = 10. ;
  /*
          TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
  */
  if (MuMuPlane_normal.Cross(KPiPlane_normal).Dot( -B0.Vect() ) > 0.0)
    angle = MuMuPlane_normal.Angle(KPiPlane_normal);
  else
    angle = -MuMuPlane_normal.Angle(KPiPlane_normal);

  return angle;
}

double psiPrimePiK_MC::getBellePlanesAngle(TLorentzVector B0, TLorentzVector K, TLorentzVector Pi, TLorentzVector muM, TLorentzVector muP) {

  TLorentzVector KBoosted = K, PiBoosted = Pi, muMBoosted = muM, muPBoosted = muP ;
  K.Boost( -B0.BoostVector() ) ; TVector3 K3 = K.Vect() ;
  PiBoosted.Boost( -B0.BoostVector() ) ; TVector3 Pi3 = PiBoosted.Vect() ;
  TVector3 Kstar3 = K3 + Pi3 ;
  muPBoosted.Boost( -B0.BoostVector() ) ; TVector3 muP3 = muPBoosted.Vect() ;
  muMBoosted.Boost( -B0.BoostVector() ) ; TVector3 muM3 = muMBoosted.Vect() ;
  TVector3 MuMu3 = muP3 + muM3 ;

  TVector3 aKplus = K3 - (Kstar3.Dot(K3)/Kstar3.Mag2())*Kstar3 ;
  TVector3 aMuplus = muP3 - (muP3.Dot(MuMu3)/MuMu3.Mag2())*MuMu3 ;

  return aKplus.Angle( aMuplus );
}

float psiPrimePiK_MC::GetThetaMuMu(TLorentzVector BVec, TLorentzVector JPsiVec, TLorentzVector MuPlusVec, float BeamEnergy, float JPsiPDG , float muonPDG) {

  TVector3 JPsiInBFrame, MuInBFrame, MuInJPsiFrame, MuInJPsiFromBFrame;
  TLorentzVector JPsiInBFrameTLVec, MuInBFrameTLVec;
  /*
  // Alessandra
  // get the momentum of the X in the in the B rest-frame : XInBFrame
  GetMomentumInMotherFrame( BVec , XVec , BeamEnergy, XInBFrame);
  XInBFrameTLVec.SetPtEtaPhiM(XInBFrame.Perp() , XInBFrame.Eta(),  XInBFrame.Phi() ,XCandPDG);
  // get the momentum of the J/psi int he X rest-frame with two steps in cascade:
  // 1) step 1 : apply B boost in lab
  GetMomentumInMotherFrame( BVec , JPsiVec, BeamEnergy, JPsiInBFrame);
  JPsiInBFrameTLVec.SetPtEtaPhiM( JPsiInBFrame.Perp() , JPsiInBFrame.Eta(),JPsiInBFrame.Phi() ,JPsiPDG );
  // 2) step 2: apply X boost in B rest-frame
  GetMomentumInMotherFrame( XInBFrameTLVec, JPsiInBFrameTLVec, BeamEnergy, JPsiInXFromBFrame);
  JPsiInXFromBFrameTLVec.SetPtEtaPhiM(JPsiInXFromBFrame.Perp() ,JPsiInXFromBFrame.Eta(),JPsiInXFromBFrame.Phi() ,JPsiPDG );
  */

  // B0 - >J/psi K pi
  
  // get the momentum of the J/psi in the in the B rest-frame : JPsiInBFrame
  GetMomentumInMotherFrame( BVec , JPsiVec , BeamEnergy, JPsiInBFrame);
  JPsiInBFrameTLVec.SetPtEtaPhiM(JPsiInBFrame.Perp() , JPsiInBFrame.Eta(),  JPsiInBFrame.Phi() , JPsiPDG);

  // get momentum of the mu+ in JPsi rest-frame
  GetMomentumInMotherFrame( BVec , MuPlusVec, BeamEnergy, MuInBFrame); // B boost
  MuInBFrameTLVec.SetPtEtaPhiM(MuInBFrame.Perp() , MuInBFrame.Eta(),  MuInBFrame.Phi() ,muonPDG);
  //GetMomentumInMotherFrame(XInBFrameTLVec, MuInBFrameTLVec, BeamEnergy, MuInXFromBFrame); // X in B r.f. boost
  //MuInXFromBFrameTLVec.SetPtEtaPhiM(MuInXFromBFrame.Perp() ,MuInXFromBFrame.Eta(), MuInXFromBFrame.Phi() ,muonPDG);
  //GetMomentumInMotherFrame( JPsiInXFromBFrameTLVec, MuInXFromBFrameTLVec, BeamEnergy, MuInJPsiFromXFromBFrame);
  GetMomentumInMotherFrame( JPsiInBFrameTLVec, MuInBFrameTLVec, BeamEnergy, MuInJPsiFromBFrame);

  float thetaJPsi = MuInJPsiFromBFrame.Angle(JPsiInBFrame);
  
  /*
  // lab frame
  GetMomentumInMotherFrame( JPsiVec , MuPlusVec, BeamEnergy, MuInJPsiFrame); // B boost
  float thetaJPsi = MuInJPsiFrame.Angle(JPsiVec.Vect());
  */
  return thetaJPsi;

}

void psiPrimePiK_MC::GetMomentumInMotherFrame( TLorentzVector Mother,TLorentzVector Particle, double BeamEnergy , TVector3 &Particle_rotated){
  //Mother momentum in lab frame
  TVector3 bMother = Mother.BoostVector();
  TLorentzVector beam1(0., 0.,  BeamEnergy, BeamEnergy); // beam momentum in lab frame
  TLorentzVector beam2(0., 0., -BeamEnergy, BeamEnergy); // beam momentum in lab frame
  beam1.Boost(-bMother);      // beam momentum in JPsi rest frame
  beam2.Boost(-bMother);     // beam momentum in JPsi rest frame
  TVector3 beam1_dir = beam1.Vect().Unit();  // beam direction in Mother rest frame
  TVector3 beam2_dir = beam2.Vect().Unit();  // beam direction in Mother rest frame

  TVector3 Y = beam1_dir.Cross( beam2_dir ).Unit(); // the production plane normal
  TVector3 Z = Mother.Vect().Unit();         // Mother direction in lab frame
  TVector3 X = Y.Cross(Z).Unit();         // completes the right-handed coordinate

  Particle.Boost(-bMother);     // Particle momentum in Mother rest frame


  TRotation rotation;
  rotation.SetToIdentity();
  rotation.RotateAxes(X,Y,Z);
  rotation.Invert(); // transforms coordinates from the "xyz" frame to the new frame
  //TVector3 Particle_rotated = Particle.Vect(); // io: particle coordinates in the rotated frame
  Particle_rotated = Particle.Vect();
  //Particle_rotated.Transform(rotation);

}

float psiPrimePiK_MC::GetThetaMuMu_rot(TLorentzVector BVec, TLorentzVector JPsiVec, TLorentzVector MuPlusVec, float BeamEnergy, float JPsiPDG , float muonPDG) {

  TVector3 JPsiInBFrame, MuInBFrame, MuInJPsiFrame, MuInJPsiFromBFrame;
  TLorentzVector JPsiInBFrameTLVec, MuInBFrameTLVec;
  
  // B0 - >J/psi K pi
  // get the momentum of the J/psi in the in the B rest-frame : JPsiInBFrame
  GetMomentumInMotherFrame_rot( BVec , JPsiVec , BeamEnergy, JPsiInBFrame);
  JPsiInBFrameTLVec.SetPtEtaPhiM(JPsiInBFrame.Perp() , JPsiInBFrame.Eta(),  JPsiInBFrame.Phi() , JPsiPDG);

  // get momentum of the mu+ in JPsi rest-frame
  GetMomentumInMotherFrame_rot( BVec , MuPlusVec, BeamEnergy, MuInBFrame); // B boost
  MuInBFrameTLVec.SetPtEtaPhiM(MuInBFrame.Perp() , MuInBFrame.Eta(),  MuInBFrame.Phi() ,muonPDG);
    
  GetMomentumInMotherFrame_rot( JPsiInBFrameTLVec, MuInBFrameTLVec, BeamEnergy, MuInJPsiFromBFrame);
  float thetaJPsi = MuInJPsiFromBFrame.Angle(JPsiInBFrame);
  
  return thetaJPsi;

}

void psiPrimePiK_MC::GetMomentumInMotherFrame_rot( TLorentzVector Mother,TLorentzVector Particle, double BeamEnergy , TVector3 &Particle_rotated){
  //Mother momentum in lab frame
  TVector3 bMother = Mother.BoostVector();
  TLorentzVector beam1(0., 0.,  BeamEnergy, BeamEnergy); // beam momentum in lab frame
  TLorentzVector beam2(0., 0., -BeamEnergy, BeamEnergy); // beam momentum in lab frame
  beam1.Boost(-bMother);      // beam momentum in JPsi rest frame
  beam2.Boost(-bMother);     // beam momentum in JPsi rest frame
  TVector3 beam1_dir = beam1.Vect().Unit();  // beam direction in Mother rest frame
  TVector3 beam2_dir = beam2.Vect().Unit();  // beam direction in Mother rest frame

  TVector3 Y = beam1_dir.Cross( beam2_dir ).Unit(); // the production plane normal
  TVector3 Z = Mother.Vect().Unit();         // Mother direction in lab frame
  TVector3 X = Y.Cross(Z).Unit();         // completes the right-handed coordinate

  Particle.Boost(-bMother);     // Particle momentum in Mother rest frame


  TRotation rotation;
  rotation.SetToIdentity();
  rotation.RotateAxes(X,Y,Z);
  rotation.Invert(); // transforms coordinates from the "xyz" frame to the new frame
  //TVector3 Particle_rotated = Particle.Vect(); // io: particle coordinates in the rotated frame
  Particle_rotated = Particle.Vect();
  Particle_rotated.Transform(rotation);

}

//rsync -vut --existing psiPrimePiK_MC.C cristella@cmssusy.ba.infn.it:/cmshome/cristella/work/Z_analysis/exclusive/clean_14ott/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/test/sanjay/selector/psiPrimePiK_MC.C
