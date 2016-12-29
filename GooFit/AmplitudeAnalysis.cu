#include "Variable.hh"
#include "ThreeBodiesPsiPiKPdf.hh"
#include "FitManager.hh"
#include "UnbinnedDataSet.hh"
#include "BinnedDataSet.hh"
#include "InterHistPdf.hh"
#include "FlatHistoPdf.hh"
#include "AddPdf.hh"
#include "ProdPdf.hh"
#include "MatrixPdf.hh"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TH1.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TAttLine.h"

#include "../utilities.h"
// #include "../Angles_contour.h"
// #include "../Dalitz_contour.h"
// #include "../effMasses.h"

#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <sstream>
#include <utility> // std::make_pair
#include <fstream>
#include "TLegend.h"

#include <sys/time.h> // for timeval
#include <sys/times.h> // for tms
#include <iostream>
#include "TSystem.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

//#define CUDADEBUGGING 50


timeval startTime, stopTime, totalTime;
clock_t startC, stopC;
tms startProc, stopProc;

const fptype M892 = 0.89581 ; const fptype G892 = 0.0474; // From PDG charged only K*(892)
const fptype M892e = 0.8961 ; const fptype G892e = 0.0507; // From EvtGen
const fptype M1410 = 1.414; const fptype G1410 = 0.232; // K*1410
const fptype M800 = 0.931; const fptype G800 = 0.578; // K*800
const fptype M1430_0 = 1.425; const fptype G1430_0 = 0.270; // K*1430_0
const fptype M1430_2 = 1.4324; const fptype G1430_2 = 0.109; // K*1430_2
const fptype M1780_3 = 1.776; const fptype G1780_3 = 0.159; // K*1780_3

const fptype TMATH_PI = TMath::Pi();

std::string migrad("MIGRAD"); std::string m("M");
std::string hesse("HESSE");   std::string h("H");
std::string minos("MINOS");   std::string n("N");

fptype phaseSpaceFunction(fptype x,fptype mP,fptype m1,fptype m2,fptype m3)
{
  fptype function = isnan(sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x)) ? 0 : (sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x));

  return function;
}

void debug(int line) {
  #ifdef CUDADEBUGGING
  std::cout <<"Debugging on line " <<line <<std::endl;
  #endif
}

int parDotSpin (fptype dotSpin) {
  int result = static_cast<int>(floor((dotSpin - floor(dotSpin))*10.+.1));
  return result;
}

std::string doubleToStr (fptype dbl) {
  std::ostringstream strs;
  strs << dbl;
  return strs.str();
}

Int_t compBins = 50, pdfBins = 50, dataBins = 100;
//Int_t compBins = 100, pdfBins = 100, dataBins = 100;

void addHelAmplStat(TPaveText *fitStat, TString hel, Variable* a, Variable* b) {
  TString a_value = TString::Format("a_{%s} = %.2f",hel.Data(),a->value) ;
  TString a_error = TString::Format("#pm %.2f",a->error) ;
  if (a->fixed) a_error = "fixed";

  TString b_value = TString::Format("b_{%s} = %.2f",hel.Data(),b->value) ;
  TString b_error = TString::Format("#pm %.2f",b->error) ;
  if (b->fixed) b_error = "fixed";

  fitStat->AddText(TString::Format("%s %s, %s %s",a_value.Data(),a_error.Data(),b_value.Data(),b_error.Data()));
}

void printinstruction() {

  std::cerr << "======= Instructions \n"
  	        << "\t-h,--help \t\t Show this help message\n"
            << "\t-evtGen \t\t\t Select EvtGen dataset and p.d.f. parameters\n"
            << "\t-eff \t\t\t Add efficiency calculations\n"
	          << "\t-n <events> \t\t Specify the number of events to use\n"
  	        << "\t-r <path> \t\t Read Generated Events from txt in <path>\n"
            << "\t-algos <algo1algo2...>\t\t\t Select the mimimisation algos in the order they should \n be performed (MIGRAD at least once) ["<<m<<" for MIGRAD, "<<h<<" for HESSE, "<<n<<" for MINOS]\n (e.g -algo "<<h<<m<<h<<n<<" for HESSE MIGRAD HESSE MINOS - default: MIGRAD only)"
  	        << "\t-b1 <b1> \t\t Select binning for MassKPi (for normalisation & integration, default: 40)\n"
            << "\t-b2 <b2> \t\t Select binning for CosMuMu (for normalisation & integration, default: 40)\n"
            << "\t-b3 <b3> \t\t Select binning for massPsiPi (for normalisation & integration, default: 40)\n"
            << "\t-b4 <b4> \t\t Select binning for Phi (for normalisation & integration, default: 40)\n"
            << "\t-p1 <p> \t\t Select p.d.f. plotting binning finenness (default: " <<pdfBins <<") for MassKPi \n"
            << "\t-p2 <p> \t\t Select p.d.f. plotting binning finenness (default: " <<pdfBins <<") for CosMuMu \n"
            << "\t-p3 <p> \t\t Select p.d.f. plotting binning finenness (default: " <<pdfBins <<") for MassPsiPi \n"
            << "\t-p4 <p> \t\t Select p.d.f. plotting binning finenness (default: " <<pdfBins <<") for Phi \n"
            << "\t-d1 <p> \t\t Select dataset binning (default: " <<dataBins <<") for MassKPi \n"
            << "\t-d2 <p> \t\t Select dataset binning (default: " <<dataBins <<") for CosMuMu \n"
            << "\t-d3 <p> \t\t Select dataset binning (default: " <<dataBins <<") for MassPsiPi \n"
            << "\t-d4 <p> \t\t Select dataset binning (default: " <<dataBins <<") for Phi \n"
            << "\t-Bb <Bb> \t\t Select bound limits for b parameter (default: 9999)\n"
            << "\t-k800 \t\t\t Add K*_0(800) to p.d.f.\n"
            << "\t-k892 \t\t\t Add K*_1(892) to p.d.f.\n"
            << "\t-k1410 \t\t\t Add K*_1(1410) to p.d.f.\n"
            << "\t-k1430_0 \t\t Add K*_0(1430) to p.d.f.\n"
            << "\t-k1430_2 \t\t Add K*_2(1430) to p.d.f.\n"
            << "\t-k1780 \t\t\t Add K*_3(1780) to p.d.f.\n"
            << "\t-Bkg \t\t\t Add Phase Space Background to p.d.f.\n"
            << std::endl;
}


int main(int argc, char** argv) {
  debug(__LINE__);

  char bufferstring[1024];

  unsigned int events = 10000;
  unsigned int nKstars = 0;

  unsigned int bin1 = compBins, bin2 = compBins, bin3 = compBins, bin4 = compBins;
  unsigned int datapoints1 = dataBins, datapoints2 = dataBins, datapoints3 = dataBins, datapoints4 = dataBins;
  unsigned int plottingfine1 = pdfBins, plottingfine2 = pdfBins, plottingfine3 = pdfBins, plottingfine4 = pdfBins;

  fptype aMax = +9999.;
  fptype bMax = +9999.;

  bool k892Star = false;
  bool k800Star = false;
  bool k1410Star = false;
  bool k1430Star0 = false;
  bool k1430Star2 = false;
  bool k1780Star = false;

  bool bkgPhaseSpace = false;
  bool effPdfProd = false;

  bool evtGen = false;

  //bool hesse = false;

  std::vector<std::string> algos;
  algos.push_back(migrad);

  TString datasetName = "Kstars";
  std::string underscores = "__";
  TString plotsDir = "./plots";
  std::vector< std::string> kStarNames;

  TH2F* relEffTH2 = 0;

    if (argc<=1)
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
  	      if (!(ss >> events))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-b1")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> bin1))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-b2")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> bin2))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-b3")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> bin3))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-b4")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> bin4))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-p1")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> plottingfine1))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-p2")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> plottingfine2))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-p3")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> plottingfine3))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-p4")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> plottingfine4))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-Bb")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> bMax))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	  else
  	    bMax = TMATH_PI;
  	}
        else if (arg == "-k892")
  	{
  	  k892Star = true;
  	  ++nKstars;
  	}
        else if (arg == "-k800")
  	{
  	  k800Star = true;
  	  ++nKstars;
  	}
        else if (arg == "-k1410")
  	{
  	  k1410Star = true;
  	  ++nKstars;
  	}
        else if (arg == "-k1430_0")
  	{
  	  k1430Star0 = true;
  	  ++nKstars;
  	}
        else if (arg == "-k1430_2")
  	{
  	  k1430Star2 = true;
  	  ++nKstars;
  	}
        else if (arg == "-k1780")
  	{
  	  k1780Star = true;
  	  ++nKstars;
  	}
    else if (arg == "-Bkg")
    {
      bkgPhaseSpace = true;
    }
    else if (arg == "-eff")
    {
      effPdfProd = true;
    }
        else if (arg == "-d1")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> datapoints1))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-d2")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> datapoints2))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-d3")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> datapoints3))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
        else if (arg == "-d4")
  	{
  	  if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;
  	      std::istringstream ss(argv[i]);
  	      if (!(ss >> datapoints4))
  		{
  		  std::cerr << "Invalid number " << argv[i] << '\n';
  		  exit(1);
  		}
  	    }
  	}
    //     else if (arg == "-H")
  	// {
  	//   hesse = true;
  	// }
        else if (arg == "-algos")
    {
      algos.clear();
      if (i + 1 < argc) // Make sure we aren't at the end of argv!
  	    {
  	      i++;

  	      std::string algosInput= argv[i];
          std::size_t found = algosInput.find(m);

          if (found == std::string::npos)
          {
            std::cout << "Minimisation algorithms invalid input : MIGRAD to be called at least once \n";
            exit(1);
          }
          std::cout << "- Minimisation algorithms sequence : "<<std::endl;

          for (std::string::size_type l = 0; l < algosInput.length(); ++l)
          {
            std::string::value_type algo = algosInput[l];

            if(algo==m)
            {
              algos.push_back(migrad);
              std::cout<<"  - MIGRAD "<<std::endl;
            }
            else if(algo==h)
            {
              algos.push_back(hesse);
              std::cout<<"  - HESSE "<<std::endl;
            }
            else if(algo==n)
            {
              algos.push_back(minos);
              std::cout<<"  - MINOS "<<std::endl;
            }
            else std:: cout<<"  - \""<<algo<<"\" invalid input, ignored "<<std::endl;
          }

  	    }
      }


    }

      if (bin1 > plottingfine1)
        cout <<"WARNING! Bins for normalisation & integration (" <<bin1 <<") are more than bins for p.d.f. plotting (" <<plottingfine1 <<")\n" <<endl;
      if (bin2 > plottingfine2)
        cout <<"WARNING! Bins for normalisation & integration (" <<bin2 <<") are more than bins for p.d.f. plotting (" <<plottingfine2 <<")\n" <<endl;
      if (bin3 > plottingfine3)
        cout <<"WARNING! Bins for normalisation & integration (" <<bin3 <<") are more than bins for p.d.f. plotting (" <<plottingfine3 <<")\n" <<endl;
      if (bin4 > plottingfine4)
        cout <<"WARNING! Bins for normalisation & integration (" <<bin4 <<") are more than bins for p.d.f. plotting (" <<plottingfine4 <<")\n" <<endl;

  TString plotsName = "";
  TString extension = "eps"; extension = "png";

  if (!nKstars) {
    cout <<"No K* selected (K892,K800,K1410,K1430) please see instructions below" <<endl;
    printinstruction();
    return 1;
  } else {
    cout <<"- Performing Amplitude Analysis fit with\n  " <<nKstars <<" K*(s) on\n  " <<events <<" events, using\n  " <<bin1 <<" bins for normalisation & integration and\n  " <<plottingfine1 <<" bins for p.d.f. plotting" <<endl;
    if (nKstars < 2) {
      datasetName = "Kstar";
      underscores = "_"; }

    if (k892Star) {
      cout <<"  - K*(892)" <<endl;
      datasetName.Append(underscores+"892_1"); plotsName.Append("__892_1");
      kStarNames.push_back("K*_{1}(892)");
    }
    if (k800Star) {
      cout <<"  - K*(800)" <<endl;
      datasetName.Append(underscores+"800_0"); plotsName.Append("__800_0");
      kStarNames.push_back("K*_{0}(800)");
    }
    if (k1410Star) {
      cout <<"  - K*(1410)" <<endl;
      datasetName.Append(underscores+"1410_1"); plotsName.Append("__1410_1");
      kStarNames.push_back("K*_{1}(1410)");
    }
    if (k1430Star0) {
      cout <<"  - K*(1430_0)" <<endl;
      datasetName.Append(underscores+"1430_0"); plotsName.Append("__1430_0");
      kStarNames.push_back("K*_{0}(1430)");
    }
    if (k1430Star2) {
      cout <<"  - K*(1430_2)" <<endl;
      datasetName.Append(underscores+"1430_2"); plotsName.Append("__1430_2");
      kStarNames.push_back("K*_{2}(1430)");}
    if (k1780Star) {
      cout <<"  - K*(1780_3)" <<endl;
      datasetName.Append(underscores+"1780_3"); plotsName.Append("__1780_3");
      kStarNames.push_back("K*_{3}(1780)");}
    if (bkgPhaseSpace) {
      cout <<"  - Three Bodies Phase-space background" <<endl;
      datasetName.Append("__plus__BdToPsiPiK_PHSP"); plotsName.Append("_PHSP");
    }
  }

  fptype aMin = -aMax;
  fptype bMin = -bMax;

  debug(__LINE__);

  //CANVAS
  TCanvas* canvas = new TCanvas("","",2000,1200);

  Variable* dRadB0  = new Variable("dRadB0",5.0);
  Variable* dRadKs  = new Variable("dRadKs",1.5);
  Variable* psi_nS  = new Variable("psi_nS",1.0);

  //std::vector<Variable* > amplitudeGooVars;
  //std::vector<Variable*> KParams;

  //GooFit
  Variable mBd("mBd", 5.27961) ;
  Variable mKaon("mKaon", 0.493677) ;
  Variable mPion("mPion", 0.13957018) ;

  fptype massMuMu = 0. ;
  if (psi_nS->value==1.0) massMuMu = 3.096916 ;
  else if (psi_nS->value==2.0) massMuMu = 3.686109 ;
  else {
    cout <<"psi_nS is neither 1 nor 2, please check it." <<endl;
    return 1; }
  Variable mMuMu("mMuMu", massMuMu);
  const fptype smearing = 0. ;
  Variable smear("smear",smearing) ;
  debug(__LINE__);

  //TH1F* dataHisto = new TH1F("data","data",BINS,massKPi.lowerlimit,massKPi.upperlimit);
  //TH1F pdfBkgHist ("bkg","bkg",BINS,massKPi.lowerlimit,massKPi.upperlimit);
  // B^{0} -> psi(nS) #pi^{+} K^{-}

  //cout <<"\nBdToMuMuPiK_PHSP.getVal() =\n" <<BdToMuMuPiK_PHSP->getVal() <<endl; return;

  //fptype roll=0.0;
  //fptype func=0.0;

  long int ms; struct timeval tp;

  gettimeofday(&tp,NULL);
  ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
  TRandom ranGen(ms);
  debug(__LINE__);
  /*
    for (int j = 0; j < events; ++j) {



    massKPi.value = ranGen.Uniform(massKPi.upperlimit - massKPi.lowerlimit) + massKPi.lowerlimit;
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


    for (int i = 1; i < BINS; i++) {
    dataSet.setBinContent(i-1,dataHisto->GetBinContent(i));
    }

    debug(__LINE__);
    phaseSpace->setData(&dataSet);
    phaseSpace->setFitControl(new BinnedNllFit());

    FitManager fitterNull(phaseSpace);
    fitterNull.fit();
    fitterNull.getMinuitValues();

    vector<double> ValsFondo;
    debug(__LINE__);
    phaseSpace->evaluateAtPoints(&massKPi,ValsFondo);
    fptype totalFondo=0.0;
    for (int k=0;k<BINS;k++) {

    pdfBkgHist.SetBinContent(k+1,ValsFondo[k]);
    totalFondo += ValsFondo[k];
    std::cout << ValsFondo[k]<<std::endl;

    }
    debug(__LINE__);
    for (int k=0;k<BINS;k++) {
    fptype valTot = pdfBkgHist.GetBinContent(k+1);
    valTot /= totalFondo;
    valTot *= events;
    pdfBkgHist.SetBinContent(k+1, valTot);
    cout <<" " <<pdfBkgHist.GetBinContent(k+1)<<endl;
    }
    debug(__LINE__);
    pdfBkgHist.SetFillStyle(3002);
    pdfBkgHist.SetFillColor(kGreen);

    double likeliHoodNull = 0.0;
    double likeliHoodSignal = 0.0;

    debug(__LINE__);
    TCanvas canvas("canvas","canvas",1000,1000);

    dataHisto->Draw();
    pdfBkgHist.Draw("same");

    canvas.SaveAs("./plots/test.png");
  */

  //Defining minimums and maximums
  fptype massKPi_min = 0.6, massKPi_max = 2.2;
  fptype massPsiPi_min = 3.2, massPsiPi_max = 4.9;


  TString massKPi_name = "massKPi", cosMuMu_name = "cosMuMu", massPsiPi_name = "massPsiPi", phi_name = "phi";
  TString massKPi_eff_name = "massKPiEff", massPsiPi_eff_name = "massPsiPiEff";
  Variable* massKPi = new Variable(massKPi_name.Data(),1.,massKPi_min,massKPi_max); massKPi->numbins = bin1;
  Variable* massKPiEff = new Variable(massKPi_eff_name.Data(),1.,0.6,2.2); massKPiEff->numbins = bin1;
  //Variable* massKPi = new Variable(massKPi_name.Data(),1.,0.6,1.67); massKPi->numbins = bin1;
  Variable* massPsiPi = new Variable(massPsiPi_name.Data(),TMath::Sqrt(23),massPsiPi_min,massPsiPi_max); massPsiPi->numbins = bin3;
  Variable* massPsiPiEff = new Variable(massPsiPi_eff_name.Data(),TMath::Sqrt(23),3.2,4.9); massPsiPiEff->numbins = bin3;
  // cosine of the psi(nS) helicity angle
  Variable* cosMuMu = new Variable(cosMuMu_name.Data(),0.,-1,1); cosMuMu->numbins = bin2;
  // cosine of the K* helicity angle
  //Variable* massPsiPi = new Variable(massPsiPi_name.Data(),0.,-1,1); massPsiPi->numbins = bin3;
  // angle between decay planes
  Variable* phi = new Variable(phi_name.Data(),0.25,-TMATH_PI,TMATH_PI); phi->numbins = bin4;

  //fptype ratio = ((fptype)(plottingfine))/((fptype)massKPi->numbins);
  fptype ratioMKPi = ((fptype)(plottingfine1))/((fptype)datapoints1);
  fptype ratioCosMuMu = ((fptype)(plottingfine2))/((fptype)datapoints2);
  fptype ratioMassPsiPi = ((fptype)(plottingfine3))/((fptype)datapoints3);
  fptype ratioPhi = ((fptype)(plottingfine4))/((fptype)datapoints4);

  std::vector<Variable*> obserVariables;
  obserVariables.push_back(massKPi);
  obserVariables.push_back(massPsiPi);
  obserVariables.push_back(cosMuMu);
  obserVariables.push_back(phi);

  // std::vector<Variable*> obserMasses;
  // obserMasses.push_back(massKPi);
  // obserMasses.push_back(massPsiPi);


  std::vector<Variable*> Masses, Gammas, Spins, as, bs;

  if (k892Star) {
    cout <<"\nAdding K*(892) ..." <<endl;

    if (!evtGen) {
      Masses.push_back(new Variable("K_892_Mass_0",M892));
      Gammas.push_back(new Variable("K_892_Gamma_0",G892));
      Spins.push_back(new Variable("K_892_Spin_0",1.0));
      as.push_back(new Variable("a_K_892_0",1.0));//,aMin,aMax) );
      bs.push_back(new Variable("b_K_892_0",0.0));//,bMin,bMax) );
      as.push_back(new Variable("a_K_892_p1",0.844,aMin,aMax) );
      bs.push_back(new Variable("b_K_892_p1",3.14,bMin,bMax) );
      as.push_back(new Variable("a_K_892_m1",0.196,aMin,aMax));
      bs.push_back(new Variable("b_K_892_m1",-1.7,bMin,bMax));
    } else {
      Masses.push_back(new Variable("K_892_Mass_0",M892e));
      Gammas.push_back(new Variable("K_892_Gamma_0",G892e));
      Spins.push_back(new Variable("K_892_Spin_0",1.0));
      // EvtGen
      as.push_back(new Variable("a_K_892_0",0.775));
      //bs.push_back(new Variable("b_K_892_0",0.0));
      //as.push_back(new Variable("a_K_892_0",0.775,0.50,0.8));
      bs.push_back(new Variable("b_K_892_0",0.0));
      as.push_back(new Variable("a_K_892_p1",0.159,0.14,0.17) );
      bs.push_back(new Variable("b_K_892_p1",1.563,1.4,1.57) );
      as.push_back(new Variable("a_K_892_m1",0.612,0.50,0.63));
      bs.push_back(new Variable("b_K_892_m1",2.712,1.0,2.73));
    }

  }

  if (k800Star) {
    cout <<"Adding K*(800) ..." <<endl;

    Masses.push_back(new Variable("K_800_Mass_0",M800));
    Gammas.push_back(new Variable("K_800_Gamma_0",G800));
    Spins.push_back(new Variable("K_800_Spin_0",0.0));
    as.push_back(new Variable("a_K_800_0",1.12,aMin,aMax) );
    bs.push_back(new Variable("b_K_800_0",2.3,bMin,bMax) );
  }

  if (k1410Star) {
    cout <<"Adding K*(1410) ..." <<endl;

    Masses.push_back(new Variable("K_1410_Mass_0",M1410));
    Gammas.push_back(new Variable("K_1410_Gamma_0",G1410));
    Spins.push_back(new Variable("K_1410_Spin_0",1.0));
    as.push_back(new Variable("a_K_1410_0",0.119,aMin,aMax) );
    bs.push_back(new Variable("b_K_1410_0",0.81,bMin,bMax) );

    //as.push_back(new Variable("a_K_1410_0",0.844));
    //bs.push_back(new Variable("b_K_1410_0",3.14,bMin,bMax));

    as.push_back(new Variable("a_K_1410_p1",0.123,aMin,aMax) );
    bs.push_back(new Variable("b_K_1410_p1",-1.04,bMin,bMax) );

    as.push_back(new Variable("a_K_1410_m1",0.036,aMin,aMax));
    bs.push_back(new Variable("b_K_1410_m1",0.67,bMin,bMax));
  }

  if (k1430Star0) {
    cout <<"Adding K*(1430_0) ..." <<endl;

    Masses.push_back(new Variable("K_1430_0_Mass_0",M1430_0));
    Gammas.push_back(new Variable("K_1430_0_Gamma_0",G1430_0));
    Spins.push_back(new Variable("K_1430_0_Spin_0",0.0));
    as.push_back(new Variable("a_K_1430_0_0",0.89,aMin,aMax) );
    bs.push_back(new Variable("b_K_1430_0_0",-2.17,bMin,bMax) );
  }

  if (k1430Star2) {
    cout <<"Adding K*(1430_2) ..." <<endl;

    Masses.push_back(new Variable("K_1430_2_Mass_0",M1430_2));
    Gammas.push_back(new Variable("K_1430_2_Gamma_0",G1430_2));
    Spins.push_back(new Variable("K_1430_2_Spin_0",2.0));
    as.push_back(new Variable("a_K_1430_2_0",4.66,aMin,aMax) );
    bs.push_back(new Variable("b_K_1430_2_0",-0.32,bMin,bMax) );

    //as.push_back(new Variable("a_K_1430_2_0",0.844));
    //bs.push_back(new Variable("b_K_1430_2_0",3.14,bMin,bMax));

    as.push_back(new Variable("a_K_1430_2_p1",4.65,aMin,aMax) );
    bs.push_back(new Variable("b_K_1430_2_p1",-3.05,bMin,bMax) );

    as.push_back(new Variable("a_K_1430_2_m1",1.26,aMin,aMax));
    bs.push_back(new Variable("b_K_1430_2_m1",-1.92,bMin,bMax));
  }

  if (k1780Star) {
    cout <<"Adding K*(1780)_3 ..." <<endl;

    Masses.push_back(new Variable("K_1780_3_Mass_0",M1780_3));
    Gammas.push_back(new Variable("K_1780_3_Gamma_0",G1780_3));
    Spins.push_back(new Variable("K_1780_3_Spin_0",3.0));
    as.push_back(new Variable("a_K_1780_3_0",16.8,aMin,aMax) );
    bs.push_back(new Variable("b_K_1780_3_0",-1.43,bMin,bMax) );

    //as.push_back(new Variable("a_K_1780_3_0",0.844));
    //bs.push_back(new Variable("b_K_1780_3_0",3.14,bMin,bMax));

    as.push_back(new Variable("a_K_1780_3_p1",19.1,aMin,aMax) );
    bs.push_back(new Variable("b_K_1780_3_p1",2.03,bMin,bMax) );

    as.push_back(new Variable("a_K_1780_3_m1",10.2,aMin,aMax));
    bs.push_back(new Variable("b_K_1780_3_m1",1.55,bMin,bMax));
  }

  Int_t nHelAmps = as.size();

  //DATASET

  UnbinnedDataSet dataset(obserVariables);
  //UnbinnedDataSet massesDataset(obserMasses);

  std::cout<<"Dataset : "<<std::endl;
  // std::cout<<" - dataset with "<<dataset.getNumBins()<<" bins "<<std::endl;
  // std::cout<<" - massesDataset with "<<massesDataset.getNumBins()<<" bins "<<std::endl;
  //std::cout<<" - efficiencyDatasetMasses with "<<efficiencyDatasetMasses->getNumBins()<<" bins "<<std::endl;

  //TString massKPi_title = "m(K^{-}#pi^{+})",  cosMuMu_title = "cos(#theta_{J/#psi})",  massPsiPi_title = "cos(#theta_{K*})",  phi_title = "#phi";
  TString massKPi_title = "m(K^{-}#pi^{+})",  cosMuMu_title = "cos(#theta_{J/#psi})",  massPsiPi_title = "m(J/#psi#pi^{+})",  phi_title = "#phi";
  TH1F massKPiHisto(massKPi_name+"_Histo", TString::Format("%s;%s [GeV]",massKPi_name.Data(),massKPi_title.Data()), datapoints1, massKPi->lowerlimit, massKPi->upperlimit); massKPiHisto.SetLineColor(kBlack); massKPiHisto.SetMarkerColor(kBlack);
  TH1F cosMuMuHisto(cosMuMu_name+"_Histo", TString::Format("%s;%s",cosMuMu_name.Data(),cosMuMu_title.Data()), datapoints2, cosMuMu->lowerlimit, cosMuMu->upperlimit); cosMuMuHisto.SetLineColor(kBlack); cosMuMuHisto.SetMarkerColor(kBlack);
  TH1F massPsiPiHisto(massPsiPi_name+"_Histo", massPsiPi_name+";"+massPsiPi_title+" [GeV]", datapoints3, massPsiPi->lowerlimit, massPsiPi->upperlimit); massPsiPiHisto.SetLineColor(kBlack); massPsiPiHisto.SetMarkerColor(kBlack);
  TH1F phiHisto(phi_name+"_Histo", phi_name+";"+phi_title, datapoints4, phi->lowerlimit, phi->upperlimit); phiHisto.SetLineColor(kBlack); phiHisto.SetMarkerColor(kBlack);

  //datasetName = "dataGen_B0"; //datasetName = "dataGen_B0bar";
  //datasetName.Append("_B0massConstraint");
  if (datasetName.Contains("dataGen_B0")) plotsDir.Append("/B0");
  else if (datasetName.Contains("dataGen_B0bar")) plotsDir.Append("/B0bar");

  if(evtGen) datasetName.Append("__EvtGen");
  //datasetName.Append("_mPhi");
  datasetName.Append(".txt");
  TString fullDatasetName = "./datasets/"+datasetName;
  fullDatasetName = "../datasets/"+datasetName;
  ifstream dataTxt(fullDatasetName.Data());
  Int_t totEvents = 0;
  if ( !(dataTxt.good()) ) {
    std::cout <<"No valid input at : " <<fullDatasetName <<" provided.\nReturning." <<std::endl;
    return 1;
  } else {
    totEvents = std::count(std::istreambuf_iterator<char>(dataTxt), std::istreambuf_iterator<char>(), '\n');
    if (events > totEvents) {
      cout <<"\nWARNING! The number of events requested is " <<events <<" but " <<fullDatasetName <<" contains only " <<totEvents <<" events." <<endl;
      events = totEvents;
    }
  }

  fptype var1, var2, var3, var4;

  if (dataTxt.good()) {
    Int_t evt=0;
    cout <<"\n- Reading " <<events <<" out of " <<totEvents <<" events from " <<datasetName <<" and filling variables histograms" <<endl;
    dataTxt.clear(); dataTxt.seekg (0, ios::beg);
    while( (evt < events)  &&  (dataTxt >> var1 >> var2 >> var3 >> var4) ) {
      evt++;
      massKPi->value = var1;
      massPsiPi->value = var2;
      cosMuMu->value = var3;
      phi->value = var4;

      //std::cout << massKPi->value << " - " <<cosMuMu->value << " - " << massPsiPi->value << " - " << phi->value << " - " << std::endl;
      if (Dalitz_contour_host(massKPi->value, massPsiPi->value, kFALSE, (Int_t)psi_nS->value) ) {
          	dataset.addEvent();
            //if(massesDataset.getNumEvents()==0) massesDataset.addEvent();
          	massKPiHisto.Fill(massKPi->value);
          	cosMuMuHisto.Fill(cosMuMu->value);
          	massPsiPiHisto.Fill(massPsiPi->value);
          	phiHisto.Fill(phi->value);
      }

      dataTxt.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }
  dataTxt.close();
  std::cout <<"Added " <<dataset.getNumEvents() <<" events within Dalitz border to GooFit dataset" <<std::endl;
  if (dataset.getNumEvents() < 1) {
    cout <<"No events added from "  <<fullDatasetName <<"\nReturning." <<endl;
    return 0;
  }

  //Efficiencies

  GooPdf* efficiencyHist;
  BinnedDataSet* efficiencyDatasetMasses;
  //BinnedDataSet efficiencyDatasetAngles(obserVariables,"efficiency Dataset Angles");

  if(effPdfProd)
  {
    int iVar1 = 0;
    int iVar2 = 0;
    int outcounter = 0;

    int holdBinVar1, holdBinVar2;

    TFile *effFile = TFile::Open("./effFiles/officialMC_noPtEtaCuts_JPsi_Bd2MuMuKPi_2p0Sig_4p0to6p0SB.root");
    TString relEffName = "RelEff_psi2SPi_vs_KPi_B0constr";
    //TString relEffName = "RelEff_planesAngle_vs_cos_psi2S_helicityAngle";
    relEffTH2 = (TH2F*)effFile->Get(relEffName) ;

    if(!(relEffTH2)){
      std::cout<<"Efficiency TH2 named \'"<<relEffName<<"\' NOT FOUND in found in TFile \'" <<effFile->GetName() <<"\'.\nReturning."<<std::endl;
      return -1;
    }

    std::cout<<"Efficiency TH2 read with bin x = "<<relEffTH2->GetNbinsX()<<" and bin y = "<<relEffTH2->GetNbinsY()<<std::endl;

    holdBinVar1 = obserVariables[iVar1]->numbins;
    holdBinVar2 = obserVariables[iVar2]->numbins;

    obserVariables[iVar1]->numbins = relEffTH2->GetNbinsX();
    obserVariables[iVar2]->numbins = relEffTH2->GetNbinsY();

    efficiencyDatasetMasses = new BinnedDataSet(obserVariables,"efficiency Dataset Masses");

    // for (int j = 0; j < relEffTH2->GetNbinsX(); ++j)
    //   for (int i = 0; i < relEffTH2->GetNbinsY(); ++i)
    //     if(relEffTH2->GetBinContent(i,j)!=0.0)
    //       std::cout<<"Bin center x :"<<relEffTH2->GetXaxis()->GetBinCenter(i)<<" y :"<<relEffTH2->GetYaxis()->GetBinCenter(j)<<" - content : "<<relEffTH2->GetBinContent(i,j)<<std::endl;


    for (int j = 0; j < efficiencyDatasetMasses->getNumBins(); ++j) {
      efficiencyDatasetMasses->setBinContent(j,0.0);
    }

    for (int j = 0; j < efficiencyDatasetMasses->getNumBins(); ++j) {

      efficiencyDatasetMasses->setBinContent(j,relEffTH2->GetBinContent(relEffTH2->FindBin(efficiencyDatasetMasses->getBinCenter(massKPi,j),efficiencyDatasetMasses->getBinCenter(massPsiPi,j))));
      // std::cout<<"Histo content at massKpi : "<<efficiencyDatasetMasses->getBinCenter(massKPi,j)<<" and massPsiPi : " <<efficiencyDatasetMasses->getBinCenter(massPsiPi,j)<<" is = "<<relEffTH2->GetBinContent(relEffTH2->FindBin(efficiencyDatasetMasses->getBinCenter(massKPi,j),efficiencyDatasetMasses->getBinCenter(massPsiPi,j)))<<std::endl;
      // std::cout<<"Binned dataset content : "<<efficiencyDatasetMasses->getBinContent(j)<<"at massKpi : "<<efficiencyDatasetMasses->getBinCenter(massKPi,j)<<" and massPsiPi : " <<efficiencyDatasetMasses->getBinCenter(massPsiPi,j)<<std::endl;
    }

    // for (int j = 0; j < efficiencyDatasetMasses->getNumBins(); ++j) {
    //   // std::cout<<"Histo content at massKpi : "<<efficiencyDatasetMasses->getBinCenter(massKPi,j)<<" and massPsiPi : " <<efficiencyDatasetMasses->getBinCenter(massPsiPi,j)<<" is = "<<relEffTH2->GetBinContent(relEffTH2->FindBin(efficiencyDatasetMasses->getBinCenter(massKPi,j),efficiencyDatasetMasses->getBinCenter(massPsiPi,j)))<<std::endl;
    //   // std::cout<<"Binned dataset content : "<<efficiencyDatasetMasses->getBinContent(j)<<std::endl;
    //   if(efficiencyDatasetMasses->getBinContent(j)!=0.0)
    //     std::cout<<"Dataset bin center x:"<<efficiencyDatasetMasses->getBinCenter(massKPi,j)<<" y : "<<efficiencyDatasetMasses->getBinCenter(massPsiPi,j)<<" - content : "<<efficiencyDatasetMasses->getBinContent(j)<<std::endl;
    //   //if(efficiencyDatasetMasses->getBinContent(j)!=0.0 && j>10000) j=efficiencyDatasetMasses->getNumBins()+1;
    // }

    // for (int j = 0; j < efficiencyDatasetMasses->getNumBins(); ++j) { // remove bins whose center is out of Dalitz border
    //
    //       if (!Dalitz_contour_host(efficiencyDatasetMasses->getBinCenter(obserVariables[iVar1],j), efficiencyDatasetMasses->getBinCenter(obserVariables[iVar2],j), kFALSE, (Int_t)psi_nS->value))
    //       {
    //          efficiencyDatasetMasses->setBinContent(j,0.0);
    //          ++outcounter;
    //        }
    //         //  cout <<"MassKPi = " <<efficiencyDatasetMasses->getBinCenter(massKPi,j) <<", massPsiPi = " <<efficiencyDatasetMasses->getBinCenter(massPsiPi,j) <<":\n original value = " <<std::endl; //relEffHist->weight(*massesEffVars) <<", corrected value = " <<Dalitz_contour_host(mass2KPi.getVal(), mass2PsiPi.getVal(), kTRUE, psi_nS.Atoi()) * relEffHist->weight(*massesEffVars) <<endl;
    //         //  cout<<j<<" OUT!"<<std::endl;
    //     }
        std::cout<<"Out : "<<outcounter<<" on "<<efficiencyDatasetMasses->getNumBins()<<std::endl;

        efficiencyHist = new FlatHistoPdf ("EfficiencyPdf",efficiencyDatasetMasses,obserVariables);

        // efficiencyHist = new InterHistPdf("EfficiencyPdf",efficiencyDatasetMasses,Masses,obserVariables);

        efficiencyHist->setData(&dataset);



        // for (int j = 0; j < dataset.getNumEvents(); ++j) {
        //
        //   massKPi->value = dataset.getValue(massKPi,j);
        //   massPsiPi->value = dataset.getValue(massPsiPi,j);
        //   cosMuMu->value = dataset.getValue(cosMuMu,j);
        //   phi->value = dataset.getValue(phi,j);
        //
        //    std::cout<<"Histo content at massKpi : "<<massKPi->value<<" and massPsiPi : " <<massPsiPi->value<<" is = "<<relEffTH2->GetBinContent(relEffTH2->FindBin(massKPi->value,massPsiPi->value))<<std::endl;
        //    std::cout<<"Dataset content :"<<efficiencyDatasetMasses->getBinContent(efficiencyDatasetMasses->getBinNumber())<<std::endl;
        //    std::cout<<"Pdf value :"<<efficiencyHist->getValue()<<std::endl;
        // }

        obserVariables[iVar1]->numbins = holdBinVar1;
        obserVariables[iVar2]->numbins = holdBinVar2;
        //return 0;
  }





  //PDFs
  GooPdf* totalPdf;

  vector<PdfBase*> pdfComponents;
  vector<Variable*> pdfYield;

  std::string p = "phasespace";

  Variable* sFrac = new Variable("sFrac",0.5,0.,1.0);

  GooPdf* matrix = new MatrixPdf("Kstars_signal", massKPi, cosMuMu, massPsiPi, phi,Masses,Gammas,Spins,as,bs,psi_nS,dRadB0,dRadKs);
  GooPdf* phaseSpace = new ThreeBodiesPsiPiK("phasespace",massKPi,cosMuMu,massPsiPi,phi,&mBd,&mPion,&mKaon,&mMuMu);
  GooPdf* sumPdf = new AddPdf("Kstars_signal + PhaseSpace", sFrac, matrix,phaseSpace);

  if (bkgPhaseSpace)
  {

    if(effPdfProd)
    {
      pdfComponents.push_back(sumPdf);
      pdfComponents.push_back(efficiencyHist);
      totalPdf = new ProdPdf("Kstars_signal + PhaseSpace + Efficiency",pdfComponents);

    }
    else
      totalPdf = sumPdf;
  } else
  {
    if(effPdfProd)
    {
      pdfComponents.push_back(matrix);
      pdfComponents.push_back(efficiencyHist);

      totalPdf = new ProdPdf("Kstars_signal + Efficiency",pdfComponents);
    }
    else
      totalPdf = matrix;
  }

  pdfComponents.clear();
  pdfYield.clear();

  obserVariables[0]->numbins = bin1;
  obserVariables[1]->numbins = bin2;

  totalPdf->setData(&dataset);
  //total->setData(&dataset);

  cout <<"\n- Fitting ..." <<endl;
  //FitManager fitter(total);
  // FitManager* fitter = 0;
  // if (!hesse)
  //   fitter = new FitManager(totalPdf);
  // else
  //   fitter = new FitManager(totalPdf,hesse);

  FitManager* fitter = new FitManager(totalPdf);

  totalPdf->setFitControl(new UnbinnedNllFit());


  gettimeofday(&startTime, NULL);
  startC = times(&startProc);
  //
  fitter->fitOrdered(algos);
  fitter->getMinuitValues();
  //
  stopC = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  fptype fitClocks = (stopC - startC)*10000.;

  // Bring phases within [-TMath::Pi,+TMath::Pi]
  fptype period = 2*TMATH_PI;
  for (int i = 0; i < nHelAmps; i++) {
    while (fabs(bs[i]->value) > TMATH_PI)
      bs[i]->value += bs[i]->value > 0 ? -period : period ;
  }
  gettimeofday(&startTime, NULL);
  startC = times(&startProc);
  //
  UnbinnedDataSet plottingGridData(obserVariables);
  std::vector<UnbinnedDataSet> compData;
  std::vector<std::vector<fptype> > pdfTotalValues;
  std::vector<std::vector<fptype> > pdfTotalSigValues;
  std::vector<std::vector<fptype> > pdfTotalBkgValues;
  /*std::vector<std::vector<std::vector<fptype> >  pdfCompValues;*/
  std::vector<std::vector<fptype> > pdfCompValues;

  for (int id = 0; id < nKstars; ++id) {
    compData.push_back( UnbinnedDataSet(obserVariables) );
  }

  std::vector<fptype> fractions;

  std::vector<fptype> compEvents;

  std::vector<fptype> mkpTotalProjection, mkpTotalSigProjection, mkpTotalBkgProjection;
  std::vector<fptype> cosMuMuTotalProjection, cosMuMuTotalSigProjection, cosMuMuTotalBkgProjection;
  std::vector<fptype> massPsiPiTotalProjection, massPsiPiTotalSigProjection, massPsiPiTotalBkgProjection;
  std::vector<fptype> phiTotalProjection, phiTotalSigProjection, phiTotalBkgProjection;

  massKPi->numbins = plottingfine1;
  cosMuMu->numbins = plottingfine2;
  massPsiPi->numbins = plottingfine3;
  phi->numbins = plottingfine4;

  fptype pointsMKPiXTot[massKPi->numbins];
  fptype pointsMKPiYTot[massKPi->numbins],pointsMKPiYTotSig[massKPi->numbins],pointsMKPiYTotBkg[massKPi->numbins];

  fptype pointsCosMuMuXTot[cosMuMu->numbins];
  fptype pointsCosMuMuYTot[cosMuMu->numbins];

  fptype pointsmassPsiPiXTot[massPsiPi->numbins];
  fptype pointsmassPsiPiYTot[massPsiPi->numbins];

  fptype pointsPhiXTot[phi->numbins];
  fptype pointsPhiYTot[phi->numbins];

  //Total pdf projection histos
  TH1F projMKPiHisto("projMKPiHisto", "projMKPiHisto",massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit);
  TH1F projCosMuMuHisto("projCosMuMuHisto", "projCosMuMuHisto",cosMuMu->numbins, cosMuMu->lowerlimit, cosMuMu->upperlimit);
  TH1F projmassPsiPiHisto("projmassPsiPiHisto", "projmassPsiPiHisto",massPsiPi->numbins, massPsiPi->lowerlimit, massPsiPi->upperlimit);
  TH1F projPhiHisto("projPhiHisto", "projPhiHisto",phi->numbins, phi->lowerlimit, phi->upperlimit);
  //Signal pdf projection histos
  TH1F projMKPiHistoSig("projMKPiHistoSignal", "projMKPiHistoSignal",massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit);
  TH1F projCosMuMuHistoSig("projCosMuMuHistoSignal", "projCosMuMuHistoSignal",cosMuMu->numbins, cosMuMu->lowerlimit, cosMuMu->upperlimit);
  TH1F projmassPsiPiHistoSig("projmassPsiPiHistoSignal", "projmassPsiPiHistoSignal",massPsiPi->numbins, massPsiPi->lowerlimit, massPsiPi->upperlimit);
  TH1F projPhiHistoSig("projPhiHistoSignal", "projPhiHistoSignal",phi->numbins, phi->lowerlimit, phi->upperlimit);
  //Background pdf projection histos
  TH1F projMKPiHistoBkg("projMKPiHistoBkg", "projMKPiHistoBkg",massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit);
  TH1F projCosMuMuHistoBkg("projCosMuMuHistoBkg", "projCosMuMuHistoBkg",cosMuMu->numbins, cosMuMu->lowerlimit, cosMuMu->upperlimit);
  TH1F projmassPsiPiHistoBkg("projmassPsiPiHistoBkg", "projmassPsiPiHistoBkg",massPsiPi->numbins, massPsiPi->lowerlimit, massPsiPi->upperlimit);
  TH1F projPhiHistoBkg("projPhiHistoBkg", "projPhiHistoBkg",phi->numbins, phi->lowerlimit, phi->upperlimit);


  for (int i = 0; i < massKPi->numbins; ++i)
    {
      mkpTotalProjection.push_back(0.0);
      mkpTotalSigProjection.push_back(0.0);
      mkpTotalBkgProjection.push_back(0.0);
    }

  for (int i = 0; i < cosMuMu->numbins; ++i)
    {
      cosMuMuTotalProjection.push_back(0.0);
      cosMuMuTotalSigProjection.push_back(0.0);
      cosMuMuTotalBkgProjection.push_back(0.0);
    }

  for (int i = 0; i < massPsiPi->numbins; ++i)
    {
      massPsiPiTotalProjection.push_back(0.0);
      massPsiPiTotalSigProjection.push_back(0.0);
      massPsiPiTotalBkgProjection.push_back(0.0);
    }

  for (int i = 0; i < phi->numbins; ++i)
    {
      phiTotalProjection.push_back(0.0);
      phiTotalSigProjection.push_back(0.0);
      phiTotalBkgProjection.push_back(0.0);
    }

  fptype sum = 0.0;
  fptype sumSig = 0.0;
  fptype sumBkg = 0.0;

  std::cout <<"\n- Starting plotting cycle" ;
  std::cout <<"\n- Plotting generated dataset" <<std::endl;
  for (int k = 0; k < phi->numbins; ++k) {
    phi->value = phi->lowerlimit + (phi->upperlimit - phi->lowerlimit)*(k + 0.5) / phi->numbins;
    //std::cout <<"Phi : " << k <<std::endl;
    for (int j = 0; j < cosMuMu->numbins; ++j) {
      cosMuMu->value = cosMuMu->lowerlimit + (cosMuMu->upperlimit - cosMuMu->lowerlimit)*(j + 0.5) / cosMuMu->numbins;
      //std::cout <<"CosMu : " << j <<std::endl;
      for (int a = 0; a < massPsiPi->numbins; ++a) {
	massPsiPi->value = massPsiPi->lowerlimit + (massPsiPi->upperlimit - massPsiPi->lowerlimit)*(a + 0.5) / massPsiPi->numbins;
	//std::cout <<"CosK : " << a <<std::endl;
	for (int i = 0; i < massKPi->numbins; ++i) {
	  //std::vector<std::vector<fptype> > tempValues;
	  //UnbinnedDataSet tempData(obserVariables);
	  massKPi->value = massKPi->lowerlimit + (massKPi->upperlimit - massKPi->lowerlimit)*(i + 0.5) / massKPi->numbins;
	  //std::cout <<"MKP : " << i <<std::endl;

	  /*tempData.addEvent();
	    matrix->setData(&tempData);
	    matrix->getCompProbsAtDataPoints(tempValues);

	    std::cout <<massKPi->value<<" ";
	    std::cout <<cosMuMu->value<<" ";
	    std::cout <<massPsiPi->value<<" ";
	    std::cout <<phi->value<<" ";
	    std::cout <<tempValues[0][0]<<std::endl;*/

	  //mkpTotalProjection[i]+=tempValues[0][0];
	  //sum +=tempValues[0][0];

	  plottingGridData.addEvent();
	  /*for (size_t ii = 0; ii < compData.size(); ++ii) {
	    compData[ii].addEvent();
	    }*/
	}
      }
    }
  }
  //
  stopC = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  fptype dataSetClocks = (stopC - startC)*10000.;

  ////////////////////////////////////////////////////////////////////////////////
  ///// TOTAL PDF PLOT
  ////////////////////////////////////////////////////////////////////////////////

  Float_t xMax = 0.95, yMax = 0.9;
  Float_t legLeft = 0.6, legWidth = 0.15;
  TLegend *legPlot = new TLegend(legLeft, 0.6, legLeft+legWidth, yMax); // 0.6 will be replaced later
  TPaveText *fitStat = new TPaveText(legPlot->GetX2(), 0.4, xMax, yMax, "NDC");

  std::cout <<"\n- Evaluating the total p.d.f." <<std::endl;
  totalPdf->setData(&plottingGridData);

  gettimeofday(&startTime, NULL);
  startC = times(&startProc);
  //
  totalPdf->getCompProbsAtDataPoints(pdfTotalValues);
  //std::cout <<" Vector size : " <<pdfTotalValues[0].size()<<std::endl;
  //std::cout <<" Vector proj : " <<pdfTotalValues[0].size()/massKPi->numbins<<std::endl;
  for (int k = 0; k < pdfTotalValues[0].size(); k++) {
    //std::cout <<mkpTotalProjection[k]*events/sum<<std::endl;
    sum += pdfTotalValues[0][k];
      if(bkgPhaseSpace)
      {
        sumSig += pdfTotalValues[1][k];
        sumBkg += pdfTotalValues[2][k];
      }
  }
  //
  stopC = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  fptype sumClocks = (stopC - startC)*10000.;

  gettimeofday(&startTime, NULL);
  startC = times(&startProc);
  //
  std::cout <<"\n[ Total Pdf sum : " <<sum<<" ] " <<std::endl;
  for (int k = 0; k<pdfTotalValues[0].size(); ++k) {
    //std::cout <<mkpTotalProjection[k]*events/sum<<std::endl;
    pdfTotalValues[0][k] /= sum;
    pdfTotalValues[0][k] *= events;
    if(bkgPhaseSpace)
    {
      pdfTotalValues[1][k] /= sumSig;
      pdfTotalValues[1][k] *= (events*sFrac->value);

      pdfTotalValues[2][k] /= sumBkg;
      pdfTotalValues[2][k] *= (events*(1.0-sFrac->value));
    }
  }
  //
  stopC = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  startC = times(&startProc);
  fptype normClocks = (stopC - startC)*10000.;

  //////////////////////////////////////////////////////////////////////
  //PROJECTING PDF ON THE FOUR VARIABLES (mkpi,phi,cosMuMu,cosK)
  //////////////////////////////////////////////////////////////////////
  //   Pdf evaluation vector (pdfTotalValues) structure :
  //
  //   es.
  //   b1 = 4  mKPi
  //   b2 = 5  massPsiPi
  //   b3 = 3  cosMuMu
  //   b4 = 3  phi
  //   ==================================================
  //   array ind   coordinates where the pdf is evaluated
  //   ==================================================
  //   0                phi0 cosMu0 cosk0 mkp0
  //   1                phi0 cosMu0 cosk0 mkp1
  //   2                phi0 cosMu0 cosk0 mkp2
  //   3                phi0 cosMu0 cosk0 mkp3
  //   4(b1)            phi0 cosMu0 cosk1 mkp0
  //   5                phi0 cosMu0 cosk1 mkp1
  //   . . .
  //   19               phi0 cosMu0 cosk4 mkp3
  //   20(b1*b2)        phi0 cosMu1 cosk0 mkp0
  //   . . .
  //   59               phi0 cosMu2 cosk4 mkp3
  //   60(b1*b2*b3)     phi1 cosMu0 cosk0 mkp0
  //   . . .            phi0 cosMu0 cosk0 mkp3
  //   179              phi2 cosMu2 cosk4 mkp3

  int notMPKBins = pdfTotalValues[0].size()/massKPi->numbins;
  int notCosMuMuBins = pdfTotalValues[0].size()/cosMuMu->numbins;
  int notmassPsiPiBins = pdfTotalValues[0].size()/massPsiPi->numbins;
  int notPhiBins = pdfTotalValues[0].size()/phi->numbins;

  //Mass K Pi
  for (int j = 0; j < massKPi->numbins; ++j) {
    for (int i = 0; i < notMPKBins; ++i) {
      mkpTotalProjection[j] += pdfTotalValues[0][j  +  i * massKPi->numbins];
      if(bkgPhaseSpace)
      {
        mkpTotalSigProjection[j] += pdfTotalValues[1][j  +  i * massKPi->numbins];
        mkpTotalBkgProjection[j] += pdfTotalValues[2][j  +  i * massKPi->numbins];
      }
    }
  }

  //Cos Mu Mu
  for (int j = 0; j < massPsiPi->numbins; ++j) {
    for (int k = 0; k < phi->numbins * cosMuMu->numbins; ++k) {
      for (int i = 0; i < massKPi->numbins; ++i) {
        massPsiPiTotalProjection[j] += pdfTotalValues[0][i  +  k * massKPi->numbins * massPsiPi->numbins  +  j * massKPi->numbins];
        // if(bkgPhaseSpace)
        // {
        //   massPsiPiTotalSigProjection[j] += pdfTotalValues[1][i  +  k * massKPi->numbins * massPsiPi->numbins  +  j * massKPi->numbins];
        //   massPsiPiTotalBkgProjection[j] += pdfTotalValues[2][i  +  k * massKPi->numbins * massPsiPi->numbins  +  j * massKPi->numbins];
        // }
      }
    }
  }

  // //Cos Mu Mu
  // for (int j = 0; j < cosMuMu->numbins; ++j) {
  //  for (int k = 0; k < phi->numbins*cosMuMu->numbins; ++k) {
  //    for (int i = 0; i < massKPi->numbins; ++i) {
  //      cosMuMuTotalProjection[j] += pdfTotalValues[0][i+k*massKPi->numbins*massPsiPi->numbins+j*massKPi->numbins];
  //    }
  //  }
  // }

  // //Cos K Star
  // for (int j = 0; j < massPsiPi->numbins; ++j) {
  //  for (int k = 0; k < phi->numbins; ++k) {
  //    for (int i = 0; i < massKPi->numbins*cosMuMu->numbins; ++i) {
  //      massPsiPiTotalProjection[j] += pdfTotalValues[0][i+j*massKPi->numbins*massPsiPi->numbins+k*massKPi->numbins*cosMuMu->numbins*massPsiPi->numbins];
  //    }
  //  }
  // }

  //Cos K Star
  for (int j = 0; j < cosMuMu->numbins; ++j) {
    for (int k = 0; k < phi->numbins; ++k) {
      for (int i = 0; i < massKPi->numbins*cosMuMu->numbins; ++i) {
        cosMuMuTotalProjection[j] += pdfTotalValues[0][i  +  j * massKPi->numbins * massPsiPi->numbins  +  k * massKPi->numbins * cosMuMu->numbins * massPsiPi->numbins];
        // if(bkgPhaseSpace)
        // {
        //   cosMuMuTotalSigProjection[j] += pdfTotalValues[1][i  +  j * massKPi->numbins * massPsiPi->numbins  +  k * massKPi->numbins * cosMuMu->numbins * massPsiPi->numbins];
        //   cosMuMuTotalBkgProjection[j] += pdfTotalValues[2][i  +  j * massKPi->numbins * massPsiPi->numbins  +  k * massKPi->numbins * cosMuMu->numbins * massPsiPi->numbins];
        // }
      }
    }
  }

  //Phi
  for (int j = 0; j < phi->numbins; ++j) {
    for (int k = 0; k < massPsiPi->numbins * massKPi->numbins * cosMuMu->numbins; ++k) {
        phiTotalProjection[j] += pdfTotalValues[0][k  +  j * massKPi->numbins * cosMuMu->numbins * massPsiPi->numbins];
        // if(bkgPhaseSpace)
        // {
        //   phiTotalSigProjection[j] += pdfTotalValues[1][k  +  j * massKPi->numbins * cosMuMu->numbins * massPsiPi->numbins];
        //   phiTotalBkgProjection[j] += pdfTotalValues[2][k  +  j * massKPi->numbins * cosMuMu->numbins * massPsiPi->numbins];
        // }
    }
  }

  //////////////////////////////////////////////////////////////////////
  //Timing
  stopC = times(&stopProc);
  gettimeofday(&stopTime, NULL);

  fptype projClocks = (stopC - startC)*10000.;

  //////////////////////////////////////////////////////////////////////
  //Fillling projection histograms

  for (int j = 0; j < massKPi->numbins; ++j) {
    projMKPiHisto.SetBinContent(j+1,mkpTotalProjection[j]);
    projMKPiHistoSig.SetBinContent(j+1,mkpTotalSigProjection[j]);
    projMKPiHistoBkg.SetBinContent(j+1,mkpTotalBkgProjection[j]);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  projMKPiHisto.Scale(ratioMKPi);
  projMKPiHistoSig.Scale(ratioMKPi);
  projMKPiHistoBkg.Scale(ratioMKPi);

  for (int j = 0; j < cosMuMu->numbins; ++j) {
    projCosMuMuHisto.SetBinContent(j+1,cosMuMuTotalProjection[j]);
    // projCosMuMuHistoSig.SetBinContent(j+1,cosMuMuTotalSigProjection[j]);
    // projCosMuMuHistoBkg.SetBinContent(j+1,cosMuMuTotalBkgProjection[j]);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  projCosMuMuHisto.Scale(ratioCosMuMu);
  // projCosMuMuHistoSig.Scale(ratioCosMuMu);
  // projCosMuMuHistoBkg.Scale(ratioCosMuMu);

  for (int j = 0; j < massPsiPi->numbins; ++j) {
    projmassPsiPiHisto.SetBinContent(j+1,massPsiPiTotalProjection[j]);
    // projmassPsiPiHistoSig.SetBinContent(j+1,massPsiPiTotalSigProjection[j]);
    // projmassPsiPiHistoBkg.SetBinContent(j+1,massPsiPiTotalBkgProjection[j]);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  projmassPsiPiHisto.Scale(ratioMassPsiPi);
  // projmassPsiPiHistoSig.Scale(ratioMassPsiPi);
  // projmassPsiPiHistoBkg.Scale(ratioMassPsiPi);

  for (int j = 0; j < phi->numbins; ++j) {
    projPhiHisto.SetBinContent(j+1,phiTotalProjection[j]);
    // projPhiHistoSig.SetBinContent(j+1,phiTotalSigProjection[j]);
    // projPhiHistoBkg.SetBinContent(j+1,phiTotalBkgProjection[j]);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  projPhiHisto.Scale(ratioPhi);
  // projPhiHistoSig.Scale(ratioPhi);
  // projPhiHistoBkg.Scale(ratioPhi);

  //////////////////////////////////////////////////////////////////////
  //Fillling projection histograms and TGraphs

  for (int j = 0; j < massKPi->numbins; ++j) {
    pointsMKPiXTot[j] = projMKPiHisto.GetBinCenter(j+1);
    pointsMKPiYTot[j] = projMKPiHisto.GetBinContent(j+1);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  for (int j = 0; j < cosMuMu->numbins; ++j) {
    pointsCosMuMuXTot[j] = projCosMuMuHisto.GetBinCenter(j+1);
    pointsCosMuMuYTot[j] = projCosMuMuHisto.GetBinContent(j+1);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  for (int j = 0; j < massPsiPi->numbins; ++j) {
    pointsmassPsiPiXTot[j] = projmassPsiPiHisto.GetBinCenter(j+1);
    pointsmassPsiPiYTot[j] = projmassPsiPiHisto.GetBinContent(j+1);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  for (int j = 0; j < phi->numbins; ++j) {
    pointsPhiXTot[j] = projPhiHisto.GetBinCenter(j+1);
    pointsPhiYTot[j] = projPhiHisto.GetBinContent(j+1);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  projPhiHisto.Scale(ratioPhi);

  //////////////////////////////////////////////////////////////////////
  //Setting Graphs & MultiGraphs

  TMultiGraph* multiGraphMKPi = new TMultiGraph(massKPi_name+"_MultiGraph", TString::Format("%s;%s",massKPiHisto.GetTitle(),(massKPiHisto.GetXaxis())->GetTitle()));
  TMultiGraph* multiGraphCosMuMu = new TMultiGraph(cosMuMu_name+"_MultiGraph", TString::Format("%s;%s",cosMuMuHisto.GetTitle(),cosMuMuHisto.GetXaxis()->GetTitle()));
  TMultiGraph* multiGraphmassPsiPi = new TMultiGraph(massPsiPi_name+"_MultiGraph", TString::Format("%s;%s",massPsiPiHisto.GetTitle(),massPsiPiHisto.GetXaxis()->GetTitle()));
  TMultiGraph* multiGraphPhi = new TMultiGraph(phi_name+"_MultiGraph", TString::Format("%s;%s",phiHisto.GetTitle(),phiHisto.GetXaxis()->GetTitle()));

  TGraph signalTotalPlotMKPi(massKPi->numbins,pointsMKPiXTot,pointsMKPiYTot);
  TGraph signalTotalSigPlotMKPi(massKPi->numbins,pointsMKPiXTot,pointsMKPiYTotSig);
  TGraph signalTotalBkgPlotMKPi(massKPi->numbins,pointsMKPiXTot,pointsMKPiYTotBkg);

  TGraph signalTotalPlotCosMuMu(cosMuMu->numbins,pointsCosMuMuXTot,pointsCosMuMuYTot);
  // TGraph signalTotalSigPlotCosMuMu(cosMuMu->numbins,pointsCosMuMuXTot,pointsCosMuMuYTotSig);
  // TGraph signalTotalBkgPlotCosMuMu(cosMuMu->numbins,pointsCosMuMuXTot,pointsCosMuMuYTotBkg);

  TGraph signalTotalPlotmassPsiPi(massPsiPi->numbins,pointsmassPsiPiXTot,pointsmassPsiPiYTot);
  // TGraph signalTotalSigPlotmassPsiPi(massPsiPi->numbins,pointsmassPsiPiXTot,pointsmassPsiPiYTotSig);
  // TGraph signalTotalBkgPlotmassPsiPi(massPsiPi->numbins,pointsmassPsiPiXTot,pointsmassPsiPiYTotBkg);

  TGraph signalTotalPlotPhi(phi->numbins,pointsPhiXTot,pointsPhiYTot);
  // TGraph signalTotalSigPlotPhi(phi->numbins,pointsPhiXTot,pointsPhiYTotSig);
  // TGraph signalTotalBkgPlotPhi(phi->numbins,pointsPhiXTot,pointsPhiYTotBkg);

  signalTotalPlotMKPi.SetLineColor(kRed); signalTotalPlotMKPi.SetLineWidth(2);
  signalTotalBkgPlotMKPi.SetLineColor(kRed); signalTotalBkgPlotMKPi.SetLineWidth(2);signalTotalBkgPlotMKPi.SetLineStyle(kDashDotted);
  signalTotalBkgPlotMKPi.SetLineColor(kRed); signalTotalBkgPlotMKPi.SetLineWidth(2);signalTotalBkgPlotMKPi.SetLineStyle(kDashed);

  signalTotalPlotCosMuMu.SetLineColor(kRed); signalTotalPlotCosMuMu.SetLineWidth(2);
  // signalTotalSigPlotCosMuMu.SetLineColor(kRed); signalTotalSigPlotCosMuMu.SetLineWidth(2); signalTotalSigPlotCosMuMu.SetLineStyle(kDashDotted);
  // signalTotalBkgPlotCosMuMu.SetLineColor(kRed); signalTotalBkgPlotCosMuMu.SetLineWidth(2); signalTotalSigPlotCosMuMu.SetLineStyle(kDashed);

  signalTotalPlotmassPsiPi.SetLineColor(kRed); signalTotalPlotmassPsiPi.SetLineWidth(2);
  // signalTotalSigPlotmassPsiPi.SetLineColor(kRed); signalTotalSigPlotmassPsiPi.SetLineWidth(2); signalTotalSigPlotmassPsiPi.SetLineStyle(kDashDotted);
  // signalTotalBkgPlotmassPsiPi.SetLineColor(kRed); signalTotalBkgPlotmassPsiPi.SetLineWidth(2); signalTotalBkgPlotmassPsiPi.SetLineStyle(kDashed);

  signalTotalPlotPhi.SetLineColor(kRed); signalTotalPlotPhi.SetLineWidth(2);
  // signalTotalSigPlotPhi.SetLineColor(kRed); signalTotalSigPlotPhi.SetLineWidth(2); signalTotalSigPlotPhi.SetLineStyle(kDashDotted);
  // signalTotalBkgPlotPhi.SetLineColor(kRed); signalTotalBkgPlotPhi.SetLineWidth(2); signalTotalBkgPlotPhi.SetLineStyle(kDashed);

  fptype totalIntegral = totalPdf->normalise();
  fptype compsIntegral = 0.0;
  std::cout <<"\nTotal Normalisation Factor = " <<totalIntegral <<std::endl;

  int kCounter = 0;
  Int_t nStatEntries = 0;
  Int_t amplitudeCounter = 0;
  for (size_t u=0; u<nKstars; ++u) {
    fitStat->AddText(TString::Format("\n------------------  %s  ------------------", kStarNames[kCounter].c_str()));
    ((TText*)fitStat->GetListOfLines()->Last())->SetTextColor(kCounter+3);
    addHelAmplStat(fitStat, "0", as[amplitudeCounter], bs[amplitudeCounter]); ++amplitudeCounter;
    nStatEntries +=2 ;

    if (Spins[u]->value > 0.) {
      addHelAmplStat(fitStat, "+1", as[amplitudeCounter], bs[amplitudeCounter]); ++amplitudeCounter;
      addHelAmplStat(fitStat, "-1", as[amplitudeCounter], bs[amplitudeCounter]); ++amplitudeCounter;
      nStatEntries +=2 ;
    }

    ++kCounter;
  }

  fitStat->SetTextAlign(12);
  fitStat->SetShadowColor(0);
  fitStat->SetFillColor(0);

  totalPdf->clearCurrentFit();

  legPlot->AddEntry(&massKPiHisto, "Generated data", "lpe");
  legPlot->AddEntry(&signalTotalPlotMKPi, "Total fit", "l");
  if(bkgPhaseSpace)
  {
    legPlot->AddEntry(&signalTotalBkgPlotMKPi, "Phase space only", "l");
    legPlot->AddEntry(&signalTotalSigPlotMKPi, "K* signal only", "l");
  }
  //multiGraphMKPi->Add(&signalTotalPlot,"L");
  ////////////////////////////////////////////////////////////////////////////////
  ///// COMPONENTS PDF PLOT
  ////////////////////////////////////////////////////////////////////////////////

  std::vector<fptype> originalAs;
  std::vector<fptype> originalBs;

  for (int i = 0; i < nHelAmps; i++) {
    originalAs.push_back(as[i]->value);
    originalBs.push_back(bs[i]->value);
  }

  kCounter = 0;

  std::vector<Variable*> MassesPlot;
  std::vector<Variable*> GammasPlot;
  std::vector<Variable*> SpinsPlot;
  std::vector<Variable*> asPlot;
  std::vector<Variable*> bsPlot;

  std::vector<TH1F*> compHistosMKPi;
  std::vector<TH1F*> compHistosCosMuMu;
  std::vector<TH1F*> compHistosmassPsiPi;
  std::vector<TH1F*> compHistosPhi;

  Bool_t plotSingleKstars = kTRUE; //plotSingleKstars = kFALSE;

  int lastAmplitude = 0;

  for (int i = 0; i < nKstars; ++i) {
    MassesPlot.push_back(Masses[i]);
    GammasPlot.push_back(Gammas[i]);
    SpinsPlot.push_back(Spins[i]);
  }

  for (int k = 0; k < nKstars; ++k) {

    ////////////////////////////////////////////////////////////////////////////////
    //Initialising projection vectors

    std::vector<fptype> mkpCompProjection;
    std::vector<fptype> cosMuMuCompProjection;
    std::vector<fptype> massPsiPiCompProjection;
    std::vector<fptype> phiCompProjection;

    for (int i = 0; i < massKPi->numbins; ++i) {
      mkpCompProjection.push_back(0.0);
    }

    for (int i = 0; i < cosMuMu->numbins; ++i) {
      cosMuMuCompProjection.push_back(0.0);
    }

    for (int i = 0; i < massPsiPi->numbins; ++i) {
      massPsiPiCompProjection.push_back(0.0);
    }

    for (int i = 0; i < phi->numbins; ++i) {
      phiCompProjection.push_back(0.0);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Pushing histogram for each projection
    sprintf(bufferstring,"comp_%d_plotHisto_MKPi",kCounter);
    compHistosMKPi.push_back(new TH1F(bufferstring,bufferstring,massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit));
    sprintf(bufferstring,"comp_%d_plotHisto_CMM",kCounter);
    compHistosCosMuMu.push_back(new TH1F(bufferstring,bufferstring,cosMuMu->numbins, cosMuMu->lowerlimit, cosMuMu->upperlimit));
    sprintf(bufferstring,"comp_%d_plotHisto_CKS",kCounter);
    compHistosmassPsiPi.push_back(new TH1F(bufferstring,bufferstring,massPsiPi->numbins, massPsiPi->lowerlimit, massPsiPi->upperlimit));
    sprintf(bufferstring,"comp_%d_plotHisto_PHI",kCounter);
    compHistosPhi.push_back(new TH1F(bufferstring,bufferstring,phi->numbins, phi->lowerlimit, phi->upperlimit));


    cout <<"\n- Plotting " <<kStarNames[kCounter] <<" component by setting all other components to zero" <<endl;
    sum = 0.0;

    ////////////////////////////////////////////////////////////////////////////////
    // Setting other components to zero and fixing all useful parameters
    for (int l = 0; l < nHelAmps; ++l) {
      as[l]->fixed = true;
      bs[l]->fixed = true;
    }

    // std::cout<<" Plotting KStars - Mass : "<<Masses[k]->value<<" Spin : "<<Spins[k]->value<<"Last Amplitude : "<<lastAmplitude<<std::endl;
    //For Spin = 0.0 only one component
    if (Spins[k]->value==0.0) {

      as[lastAmplitude]->fixed;
      bs[lastAmplitude]->fixed;
      // std::cout<<" - Amplitude vector pushing: "<<lastAmplitude<<" index ";
      asPlot.push_back(as[lastAmplitude]);
      bsPlot.push_back(bs[lastAmplitude]);

      for (int j = 0; j < nHelAmps; ++j) {
        if (j!=lastAmplitude) {
          // std::cout<<" putting zero: "<<j<<" index "<<std::endl;
      	  asPlot.push_back(new Variable("zero_a",0.0));
      	  bsPlot.push_back(new Variable("zero_b",0.0));
        }
      }
      ++lastAmplitude;
    } else {
      // For Spin != 0 three components
      for (int i = lastAmplitude; i <= lastAmplitude+2; ++i) {
        // std::cout<<" - Amplitude vector pushing: "<<i<<" index ";
        as[i]->fixed;
        bs[i]->fixed;
        asPlot.push_back(as[i]);
        bsPlot.push_back(bs[i]);
      }
      for (int d = 0; d < nHelAmps; ++d) {
	if (d!=lastAmplitude && d!=lastAmplitude+1 && d!=lastAmplitude+2) {
	  // std::cout<<" putting zero: "<<d<<" index "<<std::endl;
	  asPlot.push_back(new Variable("zero_a",0.0));
	  bsPlot.push_back(new Variable("zero_b",0.0));
	}}
      lastAmplitude+=2;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Normalising, integrating and evaluating the single component pdf

    sprintf(bufferstring,"Kstars_signal_plot_%d",kCounter);
    GooPdf* matrixPlot = new MatrixPdf(bufferstring, massKPi, cosMuMu, massPsiPi, phi,MassesPlot,GammasPlot,SpinsPlot,asPlot,bsPlot,psi_nS,dRadB0,dRadKs);
    matrixPlot->setData(&plottingGridData);

    matrixPlot->copyParams();
    compsIntegral = matrixPlot->normalise();

    fractions.push_back(compsIntegral/totalIntegral);

    cout <<"  Component " <<kStarNames[kCounter]<<" normalisation factor : " <<matrixPlot->normalise()<<" (fraction: " <<compsIntegral/totalIntegral*100.0<<"%)" <<endl;

    matrixPlot->getCompProbsAtDataPoints(pdfCompValues);

    matrixPlot->clearCurrentFit();

    for (int k = 0; k<pdfCompValues[0].size();k++) {
      //std::cout <<" Bin : " << k << " pdf : " << pdfCompValues[0][k] <<std::endl;
      sum += pdfCompValues[0][k];
    }

    for (int k = 0; k<pdfCompValues[0].size();k++) {
      pdfCompValues[0][k] /=sum;
      pdfCompValues[0][k] *= events;
      pdfCompValues[0][k] *= (compsIntegral/totalIntegral);
      //compHistosMKPi[kCounter]->SetBinContent(k,pdfCompValues[0][k]);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //Filling single components projections histos

    //MassKPi
    for (int j = 0; j < massKPi->numbins; ++j) {
      for (int i = 0; i < notMPKBins; ++i) {
	mkpCompProjection[j] += pdfCompValues[0][j+i*massKPi->numbins];
      }
      compHistosMKPi[kCounter]->SetBinContent(j+1,mkpCompProjection[j]);
    }

    compHistosMKPi[kCounter]->Scale(ratioMKPi);

    //Cos K Star
    for (int j = 0; j < cosMuMu->numbins; ++j) {
      for (int k = 0; k < phi->numbins; ++k) {
	for (int i = 0; i < massKPi->numbins*cosMuMu->numbins; ++i) {
	  cosMuMuCompProjection[j] += pdfCompValues[0][i+j*massKPi->numbins*massPsiPi->numbins+k*massKPi->numbins*cosMuMu->numbins*massPsiPi->numbins];
	}
      }
      compHistosCosMuMu[kCounter]->SetBinContent(j+1,cosMuMuCompProjection[j]);
    }

    compHistosCosMuMu[kCounter]->Scale(ratioCosMuMu);

    //Cos Mu Mu
    for (int j = 0; j < massPsiPi->numbins; ++j) {
      for (int k = 0; k < phi->numbins*cosMuMu->numbins; ++k) {
	for (int i = 0; i < massKPi->numbins; ++i) {
	  massPsiPiCompProjection[j] += pdfCompValues[0][i+k*massKPi->numbins*massPsiPi->numbins+j*massKPi->numbins];
	}
      }
      compHistosmassPsiPi[kCounter]->SetBinContent(j+1,massPsiPiCompProjection[j]);
    }

    compHistosmassPsiPi[kCounter]->Scale(ratioMassPsiPi);

    //Phi
    for (int j = 0; j < phi->numbins; ++j) {
      for (int k = 0; k < massPsiPi->numbins*massKPi->numbins*cosMuMu->numbins; ++k) {
	phiCompProjection[j] += pdfCompValues[0][k+j*massKPi->numbins*cosMuMu->numbins*massPsiPi->numbins];
      }
      compHistosPhi[kCounter]->SetBinContent(j+1,phiCompProjection[j]);
    }

    compHistosPhi[kCounter]->Scale(ratioPhi);

    if (nKstars > 1  &&  plotSingleKstars) {

      // Kstas components points array for each projection
      fptype pointsXCompMKPi[massKPi->numbins];
      fptype pointsYCompMKPi[massKPi->numbins];

      fptype pointsXCompCosMuMu[cosMuMu->numbins];
      fptype pointsYCompCosMuMu[cosMuMu->numbins];

      fptype pointsXCompmassPsiPi[massPsiPi->numbins];
      fptype pointsYCompmassPsiPi[massPsiPi->numbins];

      fptype pointsXCompPhi[phi->numbins];
      fptype pointsYCompPhi[phi->numbins];

      // Filling vectors for components projections graphs
      for (int k=0; k < massKPi->numbins; k++) {
        pointsXCompMKPi[k] = compHistosMKPi[kCounter]->GetBinCenter(k+1);
	      pointsYCompMKPi[k] = compHistosMKPi[kCounter]->GetBinContent(k+1);
      }
      for (int k=0;k<cosMuMu->numbins;k++) {
      	pointsXCompCosMuMu[k] = compHistosCosMuMu[kCounter]->GetBinCenter(k+1);
      	pointsYCompCosMuMu[k] = compHistosCosMuMu[kCounter]->GetBinContent(k+1);
      }
      for (int k=0;k<massPsiPi->numbins;k++) {
      	pointsXCompmassPsiPi[k] = compHistosmassPsiPi[kCounter]->GetBinCenter(k+1);
      	pointsYCompmassPsiPi[k] = compHistosmassPsiPi[kCounter]->GetBinContent(k+1);
      }
      for (int k=0;k<phi->numbins;k++) {
      	pointsXCompPhi[k] = compHistosPhi[kCounter]->GetBinCenter(k+1);
      	pointsYCompPhi[k] = compHistosPhi[kCounter]->GetBinContent(k+1);
      }

      // Filling components projections graphs
      TGraph* signalCompPlotMKPi = new TGraph(massKPi->numbins, pointsXCompMKPi, pointsYCompMKPi);
      signalCompPlotMKPi->SetLineColor(kCounter+3);
      signalCompPlotMKPi->SetLineWidth(2);
      signalCompPlotMKPi->SetLineStyle(kDashed);
      signalCompPlotMKPi->GetXaxis()->SetTitle("m(K#Pi)");
      sprintf(bufferstring,"Events / (%.3f)",(massKPi->upperlimit - massKPi->lowerlimit)/massKPi->numbins);
      signalCompPlotMKPi->GetYaxis()->SetTitle(bufferstring);
      //signalCompPlotMKPi->Draw("");
      multiGraphMKPi->Add(signalCompPlotMKPi,"L");

      //
      TGraph* signalCompPlotCosMuMu = new TGraph(cosMuMu->numbins,pointsXCompCosMuMu,pointsYCompCosMuMu);
      signalCompPlotCosMuMu->SetLineColor(kCounter+3);
      signalCompPlotCosMuMu->SetLineWidth(2);
      signalCompPlotCosMuMu->SetLineStyle(kDashed);
      signalCompPlotCosMuMu->GetXaxis()->SetTitle("m(K#Pi)");
      sprintf(bufferstring,"Events / (%.3f)",(cosMuMu->upperlimit- cosMuMu->lowerlimit)/cosMuMu->numbins);
      signalCompPlotCosMuMu->GetYaxis()->SetTitle(bufferstring);
      //signalCompPlot->Draw("");
      multiGraphCosMuMu->Add(signalCompPlotCosMuMu,"L");

      TGraph* signalCompPlotmassPsiPi = new TGraph(massPsiPi->numbins,pointsXCompmassPsiPi,pointsYCompmassPsiPi);
      signalCompPlotmassPsiPi->SetLineColor(kCounter+3);
      signalCompPlotmassPsiPi->SetLineWidth(2);
      signalCompPlotmassPsiPi->SetLineStyle(kDashed);
      signalCompPlotmassPsiPi->GetXaxis()->SetTitle("m(K#Pi)");
      sprintf(bufferstring,"Events / (%.3f)",(massPsiPi->upperlimit- massPsiPi->lowerlimit)/massPsiPi->numbins);
      signalCompPlotmassPsiPi->GetYaxis()->SetTitle(bufferstring);
      //signalCompPlot->Draw("");
      multiGraphmassPsiPi->Add(signalCompPlotmassPsiPi,"L");

      TGraph* signalCompPlotPhi = new TGraph(phi->numbins,pointsXCompPhi,pointsYCompPhi);
      signalCompPlotPhi->SetLineColor(kCounter+3);
      signalCompPlotPhi->SetLineWidth(2);
      signalCompPlotPhi->SetLineStyle(kDashed);
      signalCompPlotPhi->GetXaxis()->SetTitle("m(K#Pi)");
      sprintf(bufferstring,"Events / (%.3f)",(phi->upperlimit- phi->lowerlimit)/phi->numbins);
      signalCompPlotPhi->GetYaxis()->SetTitle(bufferstring);
      //signalCompPlot->Draw("");
      multiGraphPhi->Add(signalCompPlotPhi,"L");

      sprintf(bufferstring,"%s (%.2f %)",kStarNames[kCounter].c_str(),compsIntegral/totalIntegral*100.);
      legPlot->AddEntry(signalCompPlotMKPi,bufferstring,"l");

    } // if (plotSingleKstars)

    /*massKPiHisto.Draw("");
      signalCompPlot->Draw("sameL");

      sprintf(bufferstring,"plots/plot%d.eps",kCounter);
      canvas->SetLogy(1);
      canvas->SaveAs(bufferstring);
      canvas->Clear();*/
    ++kCounter;

    asPlot.clear();
    bsPlot.clear();
    pdfCompValues.clear();

  } // for (int k = 0; k < nHelAmps; ++k) {

  //Adding single points to plot better and total plots

  fptype pointsX[2], pointsY[2];
  pointsX[0] = massKPi->lowerlimit; pointsX[1] = massKPi->upperlimit;
  pointsY[0] = 0.01; pointsY[1] = massKPiHisto.GetMaximum();
  TGraph* pointsMKP = new TGraph(2,pointsX,pointsY);
  if(bkgPhaseSpace)
  {
    multiGraphMKPi->Add(&signalTotalBkgPlotMKPi,"L");
    //multiGraphMKPi->Add(&signalTotalSigPlotMKPi,"L");
  }
  multiGraphMKPi->Add(&signalTotalPlotMKPi,"L");
  multiGraphMKPi->Add(pointsMKP,"P");

  pointsX[0] = massPsiPi->lowerlimit; pointsX[1] = massPsiPi->upperlimit;
  pointsY[0] = massPsiPiHisto.GetMinimum();pointsY[1] = massPsiPiHisto.GetMaximum();
  TGraph* pointsCKS = new TGraph(2,pointsX,pointsY);
  // if(bkgPhaseSpace)
  // {
  //   multiGraphmassPsiPi->Add(&signalTotalBkgPlotmassPsiPi,"L");
  //   multiGraphmassPsiPi->Add(&signalTotalSigPlotmassPsiPi,"L");
  // }
  multiGraphmassPsiPi->Add(&signalTotalPlotmassPsiPi,"L");
  multiGraphmassPsiPi->Add(pointsCKS,"P");

  pointsX[0] = cosMuMu->lowerlimit; pointsX[1] = cosMuMu->upperlimit;
  pointsY[0] = cosMuMuHisto.GetMinimum(); pointsY[1] = cosMuMuHisto.GetMaximum();
  TGraph* pointsCMM = new TGraph(2,pointsX,pointsY);
  // if(bkgPhaseSpace)
  // {
  //   multiGraphCosMuMu->Add(&signalTotalBkgPlotCosMuMu,"L");
  //   multiGraphCosMuMu->Add(&signalTotalSigPlotCosMuMu,"L");
  // }
  multiGraphCosMuMu->Add(&signalTotalPlotCosMuMu,"L");
  multiGraphCosMuMu->Add(pointsCMM,"P");

  pointsX[0] = phi->lowerlimit; pointsX[1] = phi->upperlimit;
  pointsY[0] = phiHisto.GetMinimum(); pointsY[1] = phiHisto.GetMaximum();
  TGraph* pointsPHI = new TGraph(2,pointsX,pointsY);
  // multiGraphPhi->Add(&signalTotalBkgPlotPhi,"L");
  // multiGraphPhi->Add(&signalTotalSigPlotPhi,"L");
  multiGraphPhi->Add(&signalTotalPlotPhi,"L");
  multiGraphPhi->Add(pointsPHI,"P");

  ////////////////////////////////////////////////////////////////////////////////
  // PLOTTING

  legPlot->SetY1( yMax - 0.05*(legPlot->GetNRows()) ) ;
  fitStat->SetY1( yMax - 0.03*nStatEntries ) ;
  //Mass K Pi
  multiGraphMKPi->Draw("AL");
  massKPiHisto.Draw("Esame");
  legPlot->Draw(); fitStat->Draw();

  // first Logy(1) and after Logy(0), viceversa does not work
  canvas->SetLogy(1);
  canvas->SaveAs(TString::Format("%s/%s%s__logy.%s",plotsDir.Data(),massKPi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->SetLogy(0);
  canvas->SaveAs(TString::Format("%s/%s%s.%s",plotsDir.Data(),massKPi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->Clear();
  ////////////////////////////////////////////////////////////////////////////////

  //CosMuMu
  multiGraphCosMuMu->Draw("AL");
  cosMuMuHisto.Draw("Esame");
  // it's enough on the m(KPi) plot
  //legPlot->Draw(); //fitStat->Draw();

  canvas->SetLogy(1);
  canvas->SaveAs(TString::Format("%s/%s%s__logy.%s",plotsDir.Data(),cosMuMu_name.Data(),plotsName.Data(),extension.Data()));
  canvas->SetLogy(0);
  canvas->SaveAs(TString::Format("%s/%s%s.%s",plotsDir.Data(),cosMuMu_name.Data(),plotsName.Data(),extension.Data()));
  canvas->Clear();
  ////////////////////////////////////////////////////////////////////////////////

  //massPsiPi
  multiGraphmassPsiPi->Draw("AL");
  massPsiPiHisto.Draw("Esame");
  // it's enough on the m(KPi) plot
  //legPlot->Draw(); //fitStat->Draw();

  canvas->SetLogy(1);
  canvas->SaveAs(TString::Format("%s/%s%s__logy.%s",plotsDir.Data(),massPsiPi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->SetLogy(0);
  canvas->SaveAs(TString::Format("%s/%s%s.%s",plotsDir.Data(),massPsiPi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->Clear();
  ////////////////////////////////////////////////////////////////////////////////

  //Phi
  multiGraphPhi->Draw("AL");
  phiHisto.Draw("Esame");
  // it's enough on the m(KPi) plot
  //legPlot->Draw(); //fitStat->Draw();

  canvas->SetLogy(1);
  canvas->SaveAs(TString::Format("%s/%s%s__logy.%s",plotsDir.Data(),phi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->SetLogy(0);
  canvas->SaveAs(TString::Format("%s/%s%s.%s",plotsDir.Data(),phi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->Clear();
  ////////////////////////////////////////////////////////////////////////////////

  cout <<endl;
  cout <<"~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~" <<endl;
  cout << "PDF fitting time:       " << (fitClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout << "Data plotting time:     " << (dataSetClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout << "PDF sum time:           " << (sumClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout << "PDF normalization time: " << (normClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout << "PDF projection time:    " << (projClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout <<"~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~" <<endl;


  /*
    matrix->getCompProbsAtDataPoints(pdfTotalValuesNorm,events);
    for (int j = 0; j < massKPi->numbins; ++j) {
    projMKPiHisto.SetBinContent(j+1,mkpTotalProjection[j]);
    std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;

    }*/



  /*
    UnbinnedDataSet tempData(obserVariables);

    std::vector<std::vector<fptype> > tempValues;

    massKPi->value = 1.0;
    massPsiPi->value = 0.5;
    cosMuMu->value = 0.5;
    phi->value = 0.25;

    tempData.addEvent();
    matrix->getCompProbsAtDataPoints(tempValues);

    std::cout << "Pdf value : " <<tempValues[0][0]<<std::endl;
  */

  return 0;

}
