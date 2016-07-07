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
#include "ThreeBodiesPsiPiKPdf.hh"
#include "FitManager.hh"
#include "UnbinnedDataSet.hh"
#include "BinnedDataSet.hh"

#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"

#include <vector>
#include <sstream>
#include <utility> // std::make_pair
#include "TLegend.h"


#include <sys/time.h> // for timeval
#include <sys/times.h> // for tms
#include "TSystem.h"

#define BINS 50


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

const fptype MLb = 5.61951;
const fptype MBd = 5.27961;
const fptype MPsi2S = 3.686109;
const fptype MJpsi = 3.096916;
const fptype MProton = 0.938272046;
const fptype MKaon = 0.493677;
const fptype MPion = 0.13957018;


fptype phaseSpaceFunction(fptype x,fptype mP,fptype m1,fptype m2,fptype m3)
{

  fptype function = sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x);

  return function;

}


void printinstruction(){

  std::cerr << "======= Instructions \n"
  			<< "\t-h,--help \t\t\t\t Show this help message\n"
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

  // B^0 -> psi(nS) K* -> mu+ mu- K- pi+

  //RooFit
  //RooRealVar massKPi("massKPi","m(K^{-}#pi^{+}) [GeV]",1.,0.6,2.2);
  //RooRealVar cosMuMu("cosMuMu","cos(J/#psi)",0.,-1,1); // cosine of the psi(nS) helicity angle
  //RooRealVar cosKstar("cosKstar","cos(K*)",0.,-1,1); // cosine of the K* helicity angle
  //RooRealVar phi("phi","#phi",0.25,-TMath::Pi(),TMath::Pi());
  //RooArgSet kinematcVars(massKPi, cosMuMu, cosKstar, phi);

  //GooFit
  Variable massKPi("massKPi",1.,0.7,2.1); massKPi->numbins=BINS;
  BinnedDataSet dataSet(&massKPi);
  Variable cosMuMu("cosMuMu",0.,-1,1); // cosine of the psi(nS) helicity angle
  Variable cosKstar("cosKstar",0.,-1,1); // cosine of the K* helicity angle
  Variable phi("phi",0.25,-TMath::Pi(),TMath::Pi());

  const fptype dRadB0 = 5.0; const fptype dRadKs = 1.5;
  /*
  myPDF sigPDF("signal_pdf","Signal pdf", massKPi, cosMuMu, cosKstar, phi,
	       a892m1, b892m1, a892z, b892z, a892p1, b892p1) ;
  */
  TString psi_nS = "1"; //psi_nS = "2";


  //GooFit
  Variable mBd("mBd", MBd) ;
  Variable mKaon("mKaon", MKaon) ;
  Variable mPion("mPion", MPion) ;
  fptype massMuMu = 0. ;
  if (psi_nS.EqualTo("1")) massMuMu = MJpsi ;
  else if (psi_nS.EqualTo("2")) massMuMu = MPsi2S ;
  else {
    cout <<"psi_nS is neither 1 nor 2, please check it." <<endl;
    return 1; }
  //RooConstVar mMuMu("mMuMu", "m(#mu^{+}#mu^{-})", massMuMu);
  Variable mMuMu("mMuMu", massMuMu);
  const fptype smearing = 0. ;
  //RooConstVar smear("smear", "smear", smearing) ;
  Variable smear("smear",smearing) ;

  TH1F* dataHisto = new TH1F("data","data",BINS,massKPi.lowerlimit,massKPi.upperlimit);
  TH1F pdfBkgHist ("bkg","bkg",BINS,massKPi.lowerlimit,massKPi.upperlimit);
  // B^{0} -> psi(nS) #pi^{+} K^{-}
  GooPdf* phaseSpace = new ThreeBodiesPsiPiK ("phasespace",&massKPi,&mBd,&mPion,&mKaon,&mMuMu);
  //cout <<"\nBdToMuMuPiK_PHSP.getVal() =\n" <<BdToMuMuPiK_PHSP->getVal() <<endl; return;

  fptype roll=0.0;
	fptype func=0.0;

  long int ms; struct timeval tp;

  gettimeofday(&tp,NULL);
  ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
  TRandom ranGen(ms);

  for (int j = 0; j < events; ++j) {



  			massKPi.value = ranGen.Uniform(massKPi.upperlimit-massKPi.lowerlimit)+massKPi.lowerlimit;
  			func = phaseSpaceFunction(massKPi.value,MBd,MPion,MKaon,massMuMu);
  			roll = ranGen.Uniform(100);
  			if (roll > func) {
  				--j;
  				continue; }

  			if ((massKPi.value < massKPi.lowerlimit) || (massKPi.value > massKPi.upperlimit)) {
  				--j;
  				continue;}

  			dataHisto->Fill(massKPi.value);

  	}


  for (int i = 0; i < BINS; i++) {
    dataSet.setBinContent(i-1,dataHisto->GetBinContent(i));
  }

  phaseSpace->setData(&dataSet);
	phaseSpace->setFitControl(new BinnedNllFit());

	FitManager fitterNull(phaseSpace);
	fitterNull.fit();
	fitterNull.getMinuitValues();

	vector<double> ValsFondo;

	phaseSpace->evaluateAtPoints(&massKPi,ValsFondo);
	fptype totalFondo=0.0;
	for(int k=0;k<BINS;k++){

    pdfBkgHist.SetBinContent(k+1,ValsFondo[k]);
		totalFondo += ValsFondo[k];

  }

	for(int k=0;k<BINS;k++){
		fptype valTot = pdfBkgHist.GetBinContent(k+1);
		valTot /= totalFondo;
		valTot *= events;
		pdfBkgHist.SetBinContent(k+1, valTot);
		//cout<<" "<<pdfBkgHist.GetBinContent(k+1)<<endl;
		}

		pdfBkgHist.SetFillStyle(3002);
		pdfBkgHist.SetFillColor(kGreen);

		double likeliHoodNull = 0.0;
		double likeliHoodSignal = 0.0;


 TCanvas canvas("canvas","canvas",1000,1000);
 dataHisto->Draw();
 pdfBkgHist.Draw("same");
 canvas.SaveAs("./plots/test.png");
 return 0;

}
