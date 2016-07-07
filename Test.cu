//
//  Analysis.cu
//
//  AA analysis code in GooFit
//  Created by Ivan Heredia de la Cruz on 4/25/16.
//  Developed by
//
// root -l
// .x myPDF.cxx+
// .x Analysis.C
//

//RooFit
//#include "RooGlobalFunc.h"
//#include "RooRealVar.h"
//#include "RooDataSet.h"
//GooFit
#include "Variable.hh"
#include "ThreeBodiesPsiPiK.hh"

#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "TMath.h"

#include <vector>
#include <utility> // std::make_pair
#include "TLegend.h"


#include <sys/time.h> // for timeval
#include <sys/times.h> // for tms
#include "TSystem.h"

#include "myPDF.h"

//ROOFIT
//using namespace RooFit ;

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

fptype phaseSpaceFunction(fptype x,fptype mP,fptype m12,fptype m1,fptype m2,fptype m3)
{

  fptype function = sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x);

  return function;

}


void printinstruction(){

  std::cerr << "======= Instructions \n"
  			<< "\t-h,--help \t\t\t\ Show this help message\n"
  			<< "\t-n <iter> \t\t\t\t Specify the number of points to use for the kdtree\n"
  			<< "\t-r <path> <stack> <index> \t\t\t\t Read MC Toys in <stack> with <index> previously generated in <path>\n"
  			<< "\t-x <x> \t\t\t\t Select sigma single peak cuts\n"
        << "\t-y <y> \t\t\t\t Select sigma sided peack cuts\n"
        << "\t-z <z> \t\t\t\t Select sigma side cuts\n"
        //<< "\t-e <number of events> \t\t\t Select sigma cuts\n"
        <<std::endl;

}


int main(int argc, char** argv)
{
  SysInfo_t* s = new SysInfo_t();
  gSystem->GetSysInfo(s);
  Int_t nCPU = s->fCpus;

  unsigned int iter = 1000;
  unsigned int events = 3617;

  if(argc<=1)
  {
    printinstruction();
    return 0;
  }

  for (int i = 1; i < argc; ++i)
	{
		std::string arg = argv[i];
    if ((arg == "-h") || (arg == "--help"))
		{
			printinstruction();
			return 0;
		}
    else if (arg == "-n")
		{
			if (i + 1 < argc) // Make sure we aren't at the end of argv!
			{
				i++;
				std::istringstream ss(argv[i]);
				if (!(ss >> iter))
				{
					std::cerr << "Invalid number " << argv[i] << '\n';
					exit(1);

				}
			}
		}
    else if (arg == "-e")
		{
			if (i + 1 < argc) // Make sure we aren't at the end of argv!
			{
				i++;
				std::istringstream ss(argv[i]);
				if (!(ss >> events))
				{
					std::cerr << "Invalid number " << argv[i] << '\n';
					exit(1);
				}
			}
		}
  }


  //vector< pair<TString, pair< pair<fptype, fptype>, pair<fptype, fptype> > > > Kstar_spin;
  vector< pair<TString, pair<const fptype, const fptype> > > Kstar_spin;
  map< TString, pair<fptype, fptype> > helJ_map;
  // Belle B0->J/psi K+ pi- values
  // 5h20' with K*(892) + PHSP (1K+1K events)
  // 20' with K*(800) + K*(1430)0 + K*(1430)2 (2K events)

  cout <<"Adding K*(892)..." <<endl;
  const fptype M892 = 0.89581 ; const fptype G892 = 0.0474; // From PDG neutral only K*(892)
  Kstar_spin.push_back( make_pair("892_1", make_pair(M892,G892) ) ) ;
  helJ_map["892_1_0"] = make_pair(1.,0.); helJ_map["892_1_p1"] = make_pair(0.844,3.14); helJ_map["892_1_m1"] = make_pair(0.196,-1.7);

  cout <<"Adding K*(800)..." <<endl;
  //const fptype M800 = 0.682; const fptype G800 = 0.547; // From PDG
  const fptype M800 = 0.931; const fptype G800 = 0.578; // From Belle
  Kstar_spin.push_back( make_pair("800_0", make_pair(M800,G800) ) ) ;
  helJ_map["800_0_0"] = make_pair(1.12,2.3);

  cout <<"Adding K*(1410)..." <<endl;
  const fptype M1410 = 1.414; const fptype G1410 = 0.232;
  Kstar_spin.push_back( make_pair("1410_1", make_pair(M1410,G1410) ) ) ;
  helJ_map["1410_1_0"] = make_pair(0.119,0.81); helJ_map["1410_1_p1"] = make_pair(0.123,-1.04); helJ_map["1410_1_m1"] = make_pair(0.036,0.67);

  cout <<"Adding K*(1430)_0..." <<endl;
  const fptype M1430_0 = 1.425; const fptype G1430_0 = 0.270;
  Kstar_spin.push_back( make_pair("1430_0", make_pair(M1430_0,G1430_0) ) ) ;
  helJ_map["1430_0_0"] = make_pair(0.89,-2.17);

  cout <<"Adding K*(1430)_2..." <<endl;
  const fptype M1430_2 = 1.4324; const fptype G1430_2 = 0.109; // From PDG neutral only
  Kstar_spin.push_back( make_pair("1430_2", make_pair(M1430_2,G1430_2) ) ) ;
  helJ_map["1430_2_0"] = make_pair(4.66,-0.32); helJ_map["1430_2_p1"] = make_pair(4.65,-3.05); helJ_map["1430_2_m1"] = make_pair(1.26,-1.92);

  //const fptype M1780_3 = 1.776; const fptype G1780_3 = 0.159; // From PDG neutral only
  //Kstar_spin.push_back( make_pair("1780_3", make_pair(M1780_3,G1780_3) ) ) ;
  //helJ_map["1780_3_0"] = make_pair(16.8,-1.43); helJ_map["1780_3_p1"] = make_pair(19.1,2.03); helJ_map["1780_3_m1"] = make_pair(10.2,1.55);

  //ROOFIT
  //vector< RooRealVar* > amplitudeGooVar;
  //GooFit
  vector <Variable*> amplitudeGooVar
  vector< TString > varNames;
  RooArgSet amplitudeVars("amplitudeVars_set");
  fptype aMin = -9999.; fptype aMax = +9999.;
  fptype bMin = -9999.; fptype bMax = +9999.;
  TString helJ[] = {"m1","0","p1"} ;

  Int_t nKstar = Kstar_spin.size();
  for (Int_t iKstar_S=0; iKstar_S<nKstar; ++iKstar_S)
    for (Int_t iHelJ=0; iHelJ<3; ++iHelJ) {
      if (Kstar_spin[iKstar_S].first.Contains("_0") && !helJ[iHelJ].EqualTo("0")) continue ;
      TString name = Kstar_spin[iKstar_S].first + "_" + helJ[iHelJ] ;
      if (helJ_map.find(name) != helJ_map.end()) {
	pair<fptype, fptype> a_b = helJ_map.find(name)->second ;
	TString aName = "a"+name; TString bName = "b"+name;
	//a_b.first = aIvan; a_b.second = bIvan;
  new Variable("Mass",18.0,0.001,14.0,50.0);

	amplitudeGooVar.push_back( new Variable(aName,aName, a_b.first, aMin, aMax) ) ; varNames.push_back( aName ) ;
	amplitudeGooVar.push_back( new Variable(bName,bName, a_b.second, bMin, bMax) ); varNames.push_back( bName ) ;
      }
      else {
	cout <<"Element \"" <<name <<"\" not found in map helJ_map, please check map filling." <<endl;
	return ;
      }
    }

  for (Int_t iVar=0; iVar<(Int_t)amplitudeGooVar.size(); ++iVar) {
    amplitudeVars.add( *amplitudeGooVar[iVar] ) ;
  }

  // B^0 -> psi(nS) K* -> mu+ mu- K- pi+

  //RooFit
  //RooRealVar massKPi("massKPi","m(K^{-}#pi^{+}) [GeV]",1.,0.6,2.2);
  //RooRealVar cosMuMu("cosMuMu","cos(J/#psi)",0.,-1,1); // cosine of the psi(nS) helicity angle
  //RooRealVar cosKstar("cosKstar","cos(K*)",0.,-1,1); // cosine of the K* helicity angle
  //RooRealVar phi("phi","#phi",0.25,-TMath::Pi(),TMath::Pi());
  //RooArgSet kinematcVars(massKPi, cosMuMu, cosKstar, phi);

  //GooFit
  Variable massKPi("massKPi","m(K^{-}#pi^{+}) [GeV]",1.,0.6,2.2);
  Variable cosMuMu("cosMuMu","cos(J/#psi)",0.,-1,1); // cosine of the psi(nS) helicity angle
  Variable cosKstar("cosKstar","cos(K*)",0.,-1,1); // cosine of the K* helicity angle
  Variable phi("phi","#phi",0.25,-TMath::Pi(),TMath::Pi());

  const fptype dRadB0 = 5.0; const fptype dRadKs = 1.5;
  /*
  myPDF sigPDF("signal_pdf","Signal pdf", massKPi, cosMuMu, cosKstar, phi,
	       a892m1, b892m1, a892z, b892z, a892p1, b892p1) ;
  */
  TString psi_nS = "1"; //psi_nS = "2";


  //GooFit
  Variable mBd("mBd", "m(B^{0})", MBd) ;
  Variable mKaon("mKaon", "m(K^{-})", MKaon) ;
  Variable mPion("mPion", "m(#pi^{+})", MPion) ;
  fptype massMuMu = 0. ;
  if (psi_nS.EqualTo("1")) massMuMu = MJpsi ;
  else if (psi_nS.EqualTo("2")) massMuMu = MPsi2S ;
  else {
    cout <<"psi_nS is neither 1 nor 2, please check it." <<endl;
    return ; }
  //RooConstVar mMuMu("mMuMu", "m(#mu^{+}#mu^{-})", massMuMu);
  Variable mMuMu("mMuMu", "m(#mu^{+}#mu^{-})", massMuMu);
  const fptype smearing = 0. ;
  //RooConstVar smear("smear", "smear", smearing) ;
  Variable smear("smear", "smear", smearing) ;

  TH1F* dataHisto = new TH1F("data","data",50,massKPi.lowerlimit,massKPi.upperlimit);
  // B^{0} -> psi(nS) #pi^{+} K^{-}
  GooPdf* phaseSpace = new ThreeBodiesPsiPiK (std::n,&massKPi,&mBd,,mPion,mKaon);
  //cout <<"\nBdToMuMuPiK_PHSP.getVal() =\n" <<BdToMuMuPiK_PHSP->getVal() <<endl; return;

  fptype roll=0.0;
	fptype func=0.0;

  for (int j = 0; j < events; ++j) {
  			massKPi.value = donram.Uniform(massKPi.upperlimit-massKPi.lowerlimit)+massKPi.lowerlimit;
  			func = background(massKPi.value);
  			roll = donram.Uniform(100);
  			if (roll > func) {
  				--j;
  				continue; }

  			if ((massKPi.value < massKPi.lowerlimit) || (massKPi.value > massKPi.upperlimit)) {
  				--j;
  				continue;}

  			dataHisto->Fill(massKPi.value);

  	}
 TCanvas canvas("canvas","canvas",1000,1000);
 dataHisto->Draw();

}
