#include "Variable.hh"
#include "ThreeBodiesPsiPiKPdf.hh"
#include "FitManager.hh"
#include "UnbinnedDataSet.hh"
#include "BinnedDataSet.hh"
#include "AddPdf.hh"
#include "MatrixPdf.hh"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TMultiGraph.h"
#include "TPaveText.h"

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
const fptype M1410 = 1.414; const fptype G1410 = 0.232; // K*1410
const fptype M800 = 0.931; const fptype G800 = 0.578; // K*800
const fptype M1430_0 = 1.425; const fptype G1430_0 = 0.270; // K*1430_0
const fptype M1430_2 = 1.4324; const fptype G1430_2 = 0.109; // K*1430_2
const fptype M1780_3 = 1.776; const fptype G1780_3 = 0.159; // K*1780_3

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

void printinstruction() {
  std::cerr << "======= Instructions \n"
  	    << "\t-h,--help \t\t Show this help message\n"
	    << "\t-n <events> \t\t Specify the number of events to use\n"
  	    << "\t-r <path> \t\t Read Generated Events from txt in <path>\n"
            << "\t-H \t\t\t Run HESSE (after MIGRAD, defautl: MIGRAD only)\n"
  	    << "\t-b1 <b1> \t\t Select binning for MassKPi (for normalisation & integration, default: 40)\n"
            << "\t-b2 <b2> \t\t Select binning for CosMuMu (for normalisation & integration, default: 40)\n"
            << "\t-b3 <b3> \t\t Select binning for massPsiPi (for normalisation & integration, default: 40)\n"
            << "\t-b4 <b4> \t\t Select binning for Phi (for normalisation & integration, default: 40)\n"
            << "\t-p1 <p> \t\t Select plotting binning finenness (default: 50) for MassKPi \n"
            << "\t-p2 <p> \t\t Select plotting binning finenness (default: 50) for CosMuMu \n"
            //<< "\t-p3 <p> \t\t Select plotting binning finenness (default: 50) for massPsiPi \n"
            << "\t-p3 <p> \t\t Select plotting binning finenness (default: 50) for MassPsiPi \n"
            << "\t-p4 <p> \t\t Select plotting binning finenness (default: 50) for Phi \n"
            << "\t-d1 <p> \t\t Select dataset binning (default: 100) for MassKPi \n"
            << "\t-d2 <p> \t\t Select dataset binning (default: 100) for CosMuMu \n"
            << "\t-d3 <p> \t\t Select dataset binning (default: 100) for MassPsiPi \n"
            << "\t-d4 <p> \t\t Select dataset binning (default: 100) for Phi \n"
            << "\t-Bb <Bb> \t\t Select bound limits for b parameter (default: 9999)\n"
            << "\t-k800 \t\t\t Add K*_0(800) to p.d.f.\n"
            << "\t-k892 \t\t\t Add K*_1(892) to p.d.f.\n"
            << "\t-k1410 \t\t\t Add K*_1(1410) to p.d.f.\n"
            << "\t-k1430_0 \t\t Add K*_0(1430) to p.d.f.\n"
            << "\t-k1430_2 \t\t Add K*_2(1430) to p.d.f.\n"
            << "\t-k1780 \t\t\t Add K*_3(1780) to p.d.f.\n"
            << std::endl;
}


int main(int argc, char** argv) {
  debug(__LINE__);

  char bufferstring[1024];

  unsigned int events = 10000;
  unsigned int nOfKstar = 0;

  unsigned int bin1 = 40;
  unsigned int bin2 = 40;
  unsigned int bin3 = 40;
  unsigned int bin4 = 40;

  unsigned int datapoints1 = 100;
  unsigned int datapoints2 = 100;
  unsigned int datapoints3 = 100;
  unsigned int datapoints4 = 100;

  unsigned int plottingfine1 = 50;
  unsigned int plottingfine2 = 50;
  unsigned int plottingfine3 = 50;
  unsigned int plottingfine4 = 50;

  fptype aMax = +9999.;
  fptype bMax = +9999.;

  bool k892Star = false;
  bool k800Star = false;
  bool k1410Star = false;
  bool k1430Star0 = false;
  bool k1430Star2 = false;
  bool k1780Star = false;

  bool hesse = false;

  std::string datasetName = "Kstars";
  TString plotsDir = "./plots/";
  std::vector< std::string> kStarNames;

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
	    bMax = TMath::Pi();
	}
      else if (arg == "-k892")
	{
	  k892Star = true;
	  ++nOfKstar;
	}
      else if (arg == "-k800")
	{
	  k800Star = true;
	  ++nOfKstar;
	}
      else if (arg == "-k1410")
	{
	  k1410Star = true;
	  ++nOfKstar;
	}
      else if (arg == "-k1430_0")
	{
	  k1430Star0 = true;
	  ++nOfKstar;
	}
      else if (arg == "-k1430_2")
	{
	  k1430Star2 = true;
	  ++nOfKstar;
	}
      else if (arg == "-k1780")
	{
	  k1780Star = true;
	  ++nOfKstar;
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
      else if (arg == "-H")
	{
	  hesse = true;
	}

    }

  TString plotsName = "";
  TString extension = "eps"; extension = "png";

  if (!nOfKstar) {
    cout <<"No K* selected (K892,K800,K1410,K1430) please see instructions below" <<endl;
    printinstruction();
    return 1;
  } else {
    cout <<"Starting Amplitude Analysis fit with " <<nOfKstar <<" K*s:" <<endl;
    if (nOfKstar < 2)
      datasetName = "Kstar";

    if (k892Star) {
      cout <<"- K*(892)"<<endl;
      datasetName += "__892_1"; plotsName.Append("__892_1");
      kStarNames.push_back("K*_{1}(892)");
    }
    if (k800Star) {
      cout<<"- K*(800)"<<endl;
      datasetName += "__800_0"; plotsName.Append("__800_0");
      kStarNames.push_back("K*_{0}(800)");
    }
    if (k1410Star) {
      cout<<"- K*(1410)"<<endl;
      datasetName += "__1410_1"; plotsName.Append("__1410_1");
      kStarNames.push_back("K*_{1}(1410)");
    }
    if (k1430Star0) {
      cout<<"- K*(1430_0)"<<endl;
      datasetName += "__1430_0"; plotsName.Append("__1430_0");
      kStarNames.push_back("K*_{0}(1430)");
    }
    if (k1430Star2) {
      cout<<"- K*(1430_2)"<<endl;
      datasetName += "__1430_2"; plotsName.Append("__1430_2");
      kStarNames.push_back("K*_{2}(1430)");}
    if (k1780Star) {
      cout<<"- K*(1780_3)"<<endl;
      datasetName += "__1780_3"; plotsName.Append("__1780_3");
      kStarNames.push_back("K*_{3}(1780)");}
  }

  datasetName += ".txt";

  fptype aMin = -aMax;
  fptype bMin = -bMax;

  debug(__LINE__);

  //CANVAS
  TCanvas* canvas = new TCanvas("","",2000,1200);

  //START TIME//
  gettimeofday(&startTime, NULL);
  startC = times(&startProc);

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
    std::cout<< ValsFondo[k]<<std::endl;

    }
    debug(__LINE__);
    for (int k=0;k<BINS;k++) {
    fptype valTot = pdfBkgHist.GetBinContent(k+1);
    valTot /= totalFondo;
    valTot *= events;
    pdfBkgHist.SetBinContent(k+1, valTot);
    cout<<" "<<pdfBkgHist.GetBinContent(k+1)<<endl;
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

  TString massKPi_name = "massKPi", cosMuMu_name = "cosMuMu", massPsiPi_name = "massPsiPi", phi_name = "phi";
  Variable* massKPi = new Variable(massKPi_name.Data(),1.,0.6,2.2); massKPi->numbins = bin1;
  Variable* massPsiPi = new Variable(massPsiPi_name.Data(),TMath::Sqrt(23),3.2,4.9); massPsiPi->numbins = bin3;
  // cosine of the psi(nS) helicity angle
  Variable* cosMuMu = new Variable(cosMuMu_name.Data(),0.,-1,1); cosMuMu->numbins = bin2;
  // cosine of the K* helicity angle
  //Variable* massPsiPi = new Variable(massPsiPi_name.Data(),0.,-1,1); massPsiPi->numbins = bin3;
  // angle between decay planes
  Variable* phi = new Variable(phi_name.Data(),0.25,-TMath::Pi(),TMath::Pi()); phi->numbins = bin4;

  GooPdf* phaseSpace = new ThreeBodiesPsiPiK("phasespace",massKPi,&mBd,&mPion,&mKaon,&mMuMu);

  //fptype ratio = ((fptype)(plottingfine))/((fptype)massKPi->numbins);
  fptype ratioMKPi = ((fptype)(plottingfine1))/((fptype)datapoints1);
  fptype ratioCosMuMu = ((fptype)(plottingfine2))/((fptype)datapoints2);
  fptype ratioMassPsiPi = ((fptype)(plottingfine3))/((fptype)datapoints3);
  fptype ratioPhi = ((fptype)(plottingfine4))/((fptype)datapoints4);

  std::vector<Variable*> obserVariables;

  obserVariables.push_back(massKPi);
  obserVariables.push_back(cosMuMu);
  obserVariables.push_back(massPsiPi);
  obserVariables.push_back(phi);

  std::vector<Variable*> Masses, Gammas, Spins, as, bs;

  if (k892Star) {
    cout <<"\nAdding K*(892) ..." <<endl;

    Masses.push_back(new Variable("K_892_Mass_0",M892));
    Gammas.push_back(new Variable("K_892_Gamma_0",G892));
    Spins.push_back(new Variable("K_892_Spin_0",1.0));
    as.push_back(new Variable("a_K_892_0",1.0));//,aMin,aMax) );
    bs.push_back(new Variable("b_K_892_0",0.0));//,bMin,bMax) );

    //as.push_back(new Variable("a_K_892_0",0.844));
    //bs.push_back(new Variable("b_K_892_0",3.14,bMin,bMax));

    Masses.push_back(new Variable("K_892_Mass_p1",M892) );
    Gammas.push_back(new Variable("K_892_Gamma_p1",G892) );
    Spins.push_back(new Variable("K_892_Spin_p1",1.0) );
    as.push_back(new Variable("a_K_892_p1",0.844,aMin,aMax) );
    bs.push_back(new Variable("b_K_892_p1",3.14,bMin,bMax) );

    Masses.push_back(new Variable("K_892_Mass_m1",M892) );
    Gammas.push_back(new Variable("K_892_Gamma_m1",G892) );
    Spins.push_back(new Variable("K_892_Spin_m1",1.0));
    as.push_back(new Variable("a_K_892_m1",0.196,aMin,aMax));
    bs.push_back(new Variable("b_K_892_m1",-1.7,bMin,bMax));
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

    Masses.push_back(new Variable("K_1410_Mass_p1",M1410) );
    Gammas.push_back(new Variable("K_1410_Gamma_p1",G1410) );
    Spins.push_back(new Variable("K_1410_Spin_p1",1.0) );
    as.push_back(new Variable("a_K_1410_p1",0.123,aMin,aMax) );
    bs.push_back(new Variable("b_K_1410_p1",-1.04,bMin,bMax) );

    Masses.push_back(new Variable("K_1410_Mass_m1",M1410) );
    Gammas.push_back(new Variable("K_1410_Gamma_m1",G1410) );
    Spins.push_back(new Variable("K_1410_Spin_m1",1.0));
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

    Masses.push_back(new Variable("K_1430_2_Mass_p1",M1430_2) );
    Gammas.push_back(new Variable("K_1430_2_Gamma_p1",G1430_2) );
    Spins.push_back(new Variable("K_1430_2_Spin_p1",2.0) );
    as.push_back(new Variable("a_K_1430_2_p1",4.65,aMin,aMax) );
    bs.push_back(new Variable("b_K_1430_2_p1",-3.05,bMin,bMax) );

    Masses.push_back(new Variable("K_1430_2_Mass_m1",M1430_2) );
    Gammas.push_back(new Variable("K_1430_2_Gamma_m1",G1430_2) );
    Spins.push_back(new Variable("K_1430_2_Spin_m1",2.0));
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

    Masses.push_back(new Variable("K_1780_3_Mass_p1",M1780_3) );
    Gammas.push_back(new Variable("K_1780_3_Gamma_p1",G1780_3) );
    Spins.push_back(new Variable("K_1780_3_Spin_p1",3.0) );
    as.push_back(new Variable("a_K_1780_3_p1",19.1,aMin,aMax) );
    bs.push_back(new Variable("b_K_1780_3_p1",2.03,bMin,bMax) );

    Masses.push_back(new Variable("K_1780_3_Mass_m1",M1780_3) );
    Gammas.push_back(new Variable("K_1780_3_Gamma_m1",G1780_3) );
    Spins.push_back(new Variable("K_1780_3_Spin_m1",3.0));
    as.push_back(new Variable("a_K_1780_3_m1",10.2,aMin,aMax));
    bs.push_back(new Variable("b_K_1780_3_m1",1.55,bMin,bMax));
  }


  GooPdf* matrix = new MatrixPdf("Kstars_signal", massKPi, cosMuMu, massPsiPi, phi,Masses,Gammas,Spins,as,bs,psi_nS,dRadB0,dRadKs);

  UnbinnedDataSet dataset(obserVariables);

  TString massKPi_title = "m(K^{-}#pi^{+})",  cosMuMu_title = "cos(#theta_{J/#psi})",  massPsiPi_title = "cos(#theta_{K*})",  phi_title = "#phi";
  TH1F massKPiHisto(massKPi_name+"_Histo", TString::Format("%s;%s [GeV]",massKPi_name.Data(),massKPi_title.Data()), datapoints1, massKPi->lowerlimit, massKPi->upperlimit); massKPiHisto.SetLineColor(kBlack); massKPiHisto.SetMarkerColor(kBlack);
  TH1F cosMuMuHisto(cosMuMu_name+"_Histo", TString::Format("%s;%s",cosMuMu_name.Data(),cosMuMu_title.Data()), datapoints2, cosMuMu->lowerlimit, cosMuMu->upperlimit); cosMuMuHisto.SetLineColor(kBlack); cosMuMuHisto.SetMarkerColor(kBlack);
  TH1F massPsiPiHisto(massPsiPi_name+"_Histo", massPsiPi_name+";"+massPsiPi_title, datapoints3, massPsiPi->lowerlimit, massPsiPi->upperlimit); massPsiPiHisto.SetLineColor(kBlack); massPsiPiHisto.SetMarkerColor(kBlack);
  TH1F phiHisto(phi_name+"_Histo", phi_name+";"+phi_title, datapoints4, phi->lowerlimit, phi->upperlimit); phiHisto.SetLineColor(kBlack); phiHisto.SetMarkerColor(kBlack);

  //ifstream dataTxt("dataset.txt");
  ifstream dataTxt(("datasets/"+datasetName).c_str());
  Int_t totEvents = 0;
  if ( !(dataTxt.good()) ) {
    std::cout <<"No valid input at : "<<datasetName <<" provided. Returning." <<std::endl;
    return 1;
  } else
    totEvents = std::count(std::istreambuf_iterator<char>(dataTxt), std::istreambuf_iterator<char>(), '\n');

  fptype var1, var2, var3, var4;

  if (dataTxt.good()) {
    Int_t evt=0;
    cout <<"\nReading " <<events <<" out of " <<totEvents <<" events from " <<datasetName <<" and filling variables histograms" <<endl;
    while( (dataTxt >> var1 >> var2 >> var3 >> var4)  &&  (evt < events)) {
      evt++;
      massKPi->value = var1;
      massPsiPi->value = var2;
      cosMuMu->value = var3;
      phi->value = var4;

      //std::cout<< massKPi->value << " - " <<cosMuMu->value << " - " << massPsiPi->value << " - " << phi->value << " - " << std::endl;
      //std::cout<< massKPi->value << " " <<cosMuMu->value << " " << massPsiPi->value << " " << phi->value << " " << std::endl;

      dataset.addEvent();
      massKPiHisto.Fill(massKPi->value);
      cosMuMuHisto.Fill(cosMuMu->value);
      massPsiPiHisto.Fill(massPsiPi->value);
      phiHisto.Fill(phi->value);
    }
  }
  dataTxt.close();
  return 0;
  matrix->setData(&dataset);
  //total->setData(&dataset);

  FitManager fitter(matrix, hesse);
  //FitManager fitter(total);
  cout <<"\nFitting ..." <<endl;
  fitter.fit();
  fitter.getMinuitValues();

  stopC = times(&stopProc);
  gettimeofday(&stopTime, NULL);

  fptype fitClocks = (stopC - startC)*10000.;

  gettimeofday(&startTime, NULL);
  startC = times(&startProc);
  //
  UnbinnedDataSet currData(obserVariables);
  std::vector<UnbinnedDataSet> compData;
  std::vector<std::vector<fptype> > pdfTotalValues;
  /*std::vector<std::vector<std::vector<fptype> >  pdfCompValues;*/
  std::vector<std::vector<fptype> > pdfCompValues;

  for (int id = 0; id < nOfKstar; ++id) {
    compData.push_back( UnbinnedDataSet(obserVariables) );
  }

  std::vector<fptype> fractions;

  std::vector<fptype> compEvents;


  std::vector<fptype> mkpTotalProjection;
  std::vector<fptype> cosMuMuTotalProjection;
  std::vector<fptype> massPsiPiTotalProjection;
  std::vector<fptype> phiTotalProjection;

  massKPi->numbins = plottingfine1;
  cosMuMu->numbins = plottingfine2;
  massPsiPi->numbins = plottingfine3;
  phi->numbins = plottingfine4;

  fptype pointsMKPiXTot[massKPi->numbins];
  fptype pointsMKPiYTot[massKPi->numbins];

  fptype pointsCosMuMuXTot[cosMuMu->numbins];
  fptype pointsCosMuMuYTot[cosMuMu->numbins];

  fptype pointsmassPsiPiXTot[massPsiPi->numbins];
  fptype pointsmassPsiPiYTot[massPsiPi->numbins];

  fptype pointsPhiXTot[phi->numbins];
  fptype pointsPhiYTot[phi->numbins];

  TH1F projMKPiHisto("projMKPiHisto", "projMKPiHisto",massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit);
  TH1F projCosMuMuHisto ("projCosMuMuHisto", "projCosMuMuHisto",cosMuMu->numbins, cosMuMu->lowerlimit, cosMuMu->upperlimit);
  TH1F projmassPsiPiHisto("projmassPsiPiHisto", "projmassPsiPiHisto",massPsiPi->numbins, massPsiPi->lowerlimit, massPsiPi->upperlimit);
  TH1F projPhiHisto("projPhiHisto", "projPhiHisto",phi->numbins, phi->lowerlimit, phi->upperlimit);


  for (int i = 0; i < massKPi->numbins; ++i) {
    mkpTotalProjection.push_back(0.0);
  }

  for (int i = 0; i < cosMuMu->numbins; ++i) {
    cosMuMuTotalProjection.push_back(0.0);
  }

  for (int i = 0; i < massPsiPi->numbins; ++i) {
    massPsiPiTotalProjection.push_back(0.0);
  }

  for (int i = 0; i < phi->numbins; ++i) {
    phiTotalProjection.push_back(0.0);
  }

  fptype sum = 0.0;

  std::cout<<"\n- Starting plotting cycle "<<std::endl;
  std::cout<<"- Plotting dataset generation "<<std::endl;
  for (int k = 0; k < phi->numbins; ++k) {
    phi->value = phi->lowerlimit + (phi->upperlimit - phi->lowerlimit)*(k + 0.5) / phi->numbins;
    //std::cout<<"Phi : "<< k <<std::endl;
    for (int j = 0; j < cosMuMu->numbins; ++j) {
      cosMuMu->value = cosMuMu->lowerlimit + (cosMuMu->upperlimit - cosMuMu->lowerlimit)*(j + 0.5) / cosMuMu->numbins;
      //std::cout<<"CosMu : "<< j <<std::endl;
      for (int a = 0; a < massPsiPi->numbins; ++a) {
	massPsiPi->value = massPsiPi->lowerlimit + (massPsiPi->upperlimit - massPsiPi->lowerlimit)*(a + 0.5) / massPsiPi->numbins;
	//std::cout<<"CosK : "<< a <<std::endl;
	for (int i = 0; i < massKPi->numbins; ++i) {
	  //std::vector<std::vector<fptype> > tempValues;
	  //UnbinnedDataSet tempData(obserVariables);
	  massKPi->value = massKPi->lowerlimit + (massKPi->upperlimit - massKPi->lowerlimit)*(i + 0.5) / massKPi->numbins;
	  //std::cout<<"MKP : "<< i <<std::endl;

	  /*tempData.addEvent();
	    matrix->setData(&tempData);
	    matrix->getCompProbsAtDataPoints(tempValues);

	    std::cout<<massKPi->value<<" ";
	    std::cout<<cosMuMu->value<<" ";
	    std::cout<<massPsiPi->value<<" ";
	    std::cout<<phi->value<<" ";
	    std::cout<<tempValues[0][0]<<std::endl;*/

	  //mkpTotalProjection[i]+=tempValues[0][0];
	  //sum +=tempValues[0][0];

	  currData.addEvent();
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
  TLegend *legPlot = new TLegend(legLeft, 0.6, legLeft+legWidth, yMax);
  TPaveText *fitStat = new TPaveText(legPlot->GetX2(), 0.4, xMax, yMax, "NDC");

  std::cout<<"\n- Evaluating the total p.d.f."<<std::endl;
  matrix->setData(&currData);

  gettimeofday(&startTime, NULL);
  startC = times(&startProc);
  //
  matrix->getCompProbsAtDataPoints(pdfTotalValues);
  //std::cout<<" Vector size : "<<pdfTotalValues[0].size()<<std::endl;
  //std::cout<<" Vector proj : "<<pdfTotalValues[0].size()/massKPi->numbins<<std::endl;
  for (int k = 0; k < pdfTotalValues[0].size(); k++) {
    //std::cout<<mkpTotalProjection[k]*events/sum<<std::endl;
    sum += pdfTotalValues[0][k];
  }
  //
  stopC = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  fptype sumClocks = (stopC - startC)*10000.;

  gettimeofday(&startTime, NULL);
  startC = times(&startProc);
  //
  std::cout<<"\n[ Total Pdf sum : "<<sum<<" ] "<<std::endl;
  for (int k = 0; k<pdfTotalValues[0].size(); ++k) {
    //std::cout<<mkpTotalProjection[k]*events/sum<<std::endl;
    pdfTotalValues[0][k] /= sum;
    pdfTotalValues[0][k] *= events;
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
    }
  }

  //Cos Mu Mu
  for (int j = 0; j < massPsiPi->numbins; ++j) {
    for (int k = 0; k < phi->numbins * cosMuMu->numbins; ++k) {
      for (int i = 0; i < massKPi->numbins; ++i) {
	massPsiPiTotalProjection[j] += pdfTotalValues[0][i  +  k * massKPi->numbins * massPsiPi->numbins  +  j * massKPi->numbins];
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
      }
    }
  }

  //Phi
  for (int j = 0; j < phi->numbins; ++j) {
    for (int k = 0; k < massPsiPi->numbins * massKPi->numbins * cosMuMu->numbins; ++k) {
      phiTotalProjection[j] += pdfTotalValues[0][k  +  j * massKPi->numbins * cosMuMu->numbins * massPsiPi->numbins];
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
    pointsMKPiXTot[j] = projMKPiHisto.GetBinCenter(j+1);
    pointsMKPiYTot[j] = mkpTotalProjection[j];
    //std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;
  }

  projMKPiHisto.Scale(ratioMKPi);

  for (int j = 0; j < cosMuMu->numbins; ++j) {
    projCosMuMuHisto.SetBinContent(j+1,cosMuMuTotalProjection[j]);
    //std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;
  }

  projCosMuMuHisto.Scale(ratioCosMuMu);

  for (int j = 0; j < massPsiPi->numbins; ++j) {
    projmassPsiPiHisto.SetBinContent(j+1,massPsiPiTotalProjection[j]);
    //std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;
  }

  projmassPsiPiHisto.Scale(ratioMassPsiPi);

  for (int j = 0; j < phi->numbins; ++j) {
    projPhiHisto.SetBinContent(j+1,phiTotalProjection[j]);
    //std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;
  }

  projPhiHisto.Scale(ratioPhi);

  //////////////////////////////////////////////////////////////////////
  //Fillling projection histograms and TGraphs

  for (int j = 0; j < massKPi->numbins; ++j) {
    pointsMKPiXTot[j] = projMKPiHisto.GetBinCenter(j+1);
    pointsMKPiYTot[j] = projMKPiHisto.GetBinContent(j+1);
    //std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;
  }

  for (int j = 0; j < cosMuMu->numbins; ++j) {
    pointsCosMuMuXTot[j] = projCosMuMuHisto.GetBinCenter(j+1);
    pointsCosMuMuYTot[j] = projCosMuMuHisto.GetBinContent(j+1);
    //std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;
  }

  for (int j = 0; j < massPsiPi->numbins; ++j) {
    pointsmassPsiPiXTot[j] = projmassPsiPiHisto.GetBinCenter(j+1);
    pointsmassPsiPiYTot[j] = projmassPsiPiHisto.GetBinContent(j+1);
    //std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;
  }

  for (int j = 0; j < phi->numbins; ++j) {
    pointsPhiXTot[j] = projPhiHisto.GetBinCenter(j+1);
    pointsPhiYTot[j] = projPhiHisto.GetBinContent(j+1);
    //std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;
  }

  projPhiHisto.Scale(ratioPhi);

  //////////////////////////////////////////////////////////////////////
  //Setting Graphs & MultiGraphs

  TMultiGraph* multiGraphMKPi = new TMultiGraph(massKPi_name+"_MultiGraph", TString::Format("%s;%s",massKPiHisto.GetTitle(),(massKPiHisto.GetXaxis())->GetTitle()));
  TMultiGraph* multiGraphCosMuMu = new TMultiGraph(cosMuMu_name+"_MultiGraph", TString::Format("%s;%s",cosMuMuHisto.GetTitle(),cosMuMuHisto.GetXaxis()->GetTitle()));
  TMultiGraph* multiGraphmassPsiPi = new TMultiGraph(massPsiPi_name+"_MultiGraph", TString::Format("%s;%s",massPsiPiHisto.GetTitle(),massPsiPiHisto.GetXaxis()->GetTitle()));
  TMultiGraph* multiGraphPhi = new TMultiGraph(phi_name+"_MultiGraph", TString::Format("%s;%s",phiHisto.GetTitle(),phiHisto.GetXaxis()->GetTitle()));

  TGraph signalTotalPlotMKPi(massKPi->numbins,pointsMKPiXTot,pointsMKPiYTot);
  TGraph signalTotalPlotCosMuMu(cosMuMu->numbins,pointsCosMuMuXTot,pointsCosMuMuYTot);
  TGraph signalTotalPlotmassPsiPi(massPsiPi->numbins,pointsmassPsiPiXTot,pointsmassPsiPiYTot);
  TGraph signalTotalPlotPhi(phi->numbins,pointsPhiXTot,pointsPhiYTot);

  signalTotalPlotMKPi.SetLineColor(kRed); signalTotalPlotMKPi.SetLineWidth(2);
  signalTotalPlotCosMuMu.SetLineColor(kRed); signalTotalPlotCosMuMu.SetLineWidth(2);
  signalTotalPlotmassPsiPi.SetLineColor(kRed); signalTotalPlotmassPsiPi.SetLineWidth(2);
  signalTotalPlotPhi.SetLineColor(kRed); signalTotalPlotPhi.SetLineWidth(2);

  fptype totalIntegral = matrix->normalise();
  fptype compsIntegral = 0.0;
  std::cout <<"\nTotal Normalisation Factor = " <<totalIntegral <<std::endl;

  int kCounter = 0;

  for (size_t u=0; u<as.size(); ++u) {

    if (Spins[u]->value==0.0) {
      fitStat->AddText(TString::Format("\n--------------- %s ---------------", kStarNames[kCounter].c_str()));
      fitStat->AddText(TString::Format("a_{0} = %.2f #pm %.2f, b_{0} = %.2f #pm %.2f",as[u]->value,as[u]->error,bs[u]->value,bs[u]->error));
    } else {
      fitStat->AddText(TString::Format("\n--------------- %s ---------------", kStarNames[kCounter].c_str()));
      fitStat->AddText(TString::Format("a_{0} = %.2f #pm %.2f, b_{0} = %.2f #pm %.2f",as[u]->value,as[u]->error,bs[u]->value,bs[u]->error));
      fitStat->AddText(TString::Format("a_{p1} = %.2f #pm %.2f, b_{p1} = %.2f #pm %.2f",as[u+1]->value,as[u+1]->error,bs[u+1]->value,bs[u+1]->error));
      fitStat->AddText(TString::Format("a_{m1} = %.2f #pm %.2f, b_{m1} = %.2f #pm %.2f",as[u+2]->value,as[u+2]->error,bs[u+2]->value,bs[u+2]->error));
      u+=2;
    }

    ++kCounter;
  }

  fitStat->SetTextAlign(12);
  fitStat->SetShadowColor(0);
  fitStat->SetFillColor(0);

  matrix->clearCurrentFit();

  legPlot->AddEntry(&massKPiHisto, "Generated data", "lpe");
  legPlot->AddEntry(&signalTotalPlotMKPi, "Total fit", "l");

  //multiGraphMKPi->Add(&signalTotalPlot,"L");
  ////////////////////////////////////////////////////////////////////////////////
  ///// COMPONENTS PDF PLOT
  ////////////////////////////////////////////////////////////////////////////////

  std::vector<fptype> originalAs;
  std::vector<fptype> originalBs;

  for (int i = 0; i < (int)as.size(); i++) {
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

  for (int k = 0; k < (int)as.size(); ++k) {

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
    //Components points array for each projection

    fptype pointsXCompMKPi[massKPi->numbins];
    fptype pointsYCompMKPi[massKPi->numbins];

    fptype pointsXCompCosMuMu[cosMuMu->numbins];
    fptype pointsYCompCosMuMu[cosMuMu->numbins];

    fptype pointsXCompmassPsiPi[massPsiPi->numbins];
    fptype pointsYCompmassPsiPi[massPsiPi->numbins];

    fptype pointsXCompPhi[phi->numbins];
    fptype pointsYCompPhi[phi->numbins];

    ////////////////////////////////////////////////////////////////////////////////
    //Pushing histograms for each component for each component

    sprintf(bufferstring,"comp_%d_plotHisto_MKPi",kCounter);
    compHistosMKPi.push_back(new TH1F(bufferstring,bufferstring,massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit));
    sprintf(bufferstring,"comp_%d_plotHisto_CMM",kCounter);
    compHistosCosMuMu.push_back(new TH1F(bufferstring,bufferstring,cosMuMu->numbins, cosMuMu->lowerlimit, cosMuMu->upperlimit));
    sprintf(bufferstring,"comp_%d_plotHisto_CKS",kCounter);
    compHistosmassPsiPi.push_back(new TH1F(bufferstring,bufferstring,massPsiPi->numbins, massPsiPi->lowerlimit, massPsiPi->upperlimit));
    sprintf(bufferstring,"comp_%d_plotHisto_PHI",kCounter);
    compHistosPhi.push_back(new TH1F(bufferstring,bufferstring,phi->numbins, phi->lowerlimit, phi->upperlimit));

    cout <<"\n- Plotting component " <<kStarNames[kCounter] <<" all other components to zero" <<endl;
    sum = 0.0;

    ////////////////////////////////////////////////////////////////////////////////
    //Setting other components to zero and fixing al usefull parameters

    for (int l = 0; l < as.size(); ++l) {
      as[l]->fixed = true;
      bs[l]->fixed = true;
    }

    //For Spin = 0.0 only one component
    if (Spins[k]->value==0.0) {

      as[k]->fixed;
      bs[k]->fixed;

      MassesPlot.push_back(Masses[k]);
      GammasPlot.push_back(Gammas[k]);
      SpinsPlot.push_back(Spins[k]);
      asPlot.push_back(as[k]);
      bsPlot.push_back(bs[k]);


      for (int j = 0; j < (int)as.size(); ++j) {
        if (j!=k) {

          MassesPlot.push_back(Masses[j]);
          GammasPlot.push_back(Gammas[j]);
          SpinsPlot.push_back(Spins[j]);
          asPlot.push_back(new Variable("zero_a",0.0));
          bsPlot.push_back(new Variable("zero_b",0.0));


        }}

    }else{
      //For Spin != 0.0 three components
      for (int i = k; i <= k+2; ++i) {


        as[i]->fixed;
        bs[i]->fixed;

        MassesPlot.push_back(Masses[i]);
        GammasPlot.push_back(Gammas[i]);
        SpinsPlot.push_back(Spins[i]);
        asPlot.push_back(as[i]);
        bsPlot.push_back(bs[i]);

      }
      for (int d = 0; d < as.size(); ++d) {
        if (d!=k && d!=k+1 && d!=k+2) {

          MassesPlot.push_back(Masses[d]);
          GammasPlot.push_back(Gammas[d]);
          SpinsPlot.push_back(Spins[d]);
          asPlot.push_back(new Variable("zero_a",0.0));
          bsPlot.push_back(new Variable("zero_b",0.0));


        }}
      k+=2;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //Normalising, integrating and evaluating the single component pdf

    sprintf(bufferstring,"Kstars_signal_plot_%d",kCounter);
    GooPdf* matrixPlot = new MatrixPdf(bufferstring, massKPi, cosMuMu, massPsiPi, phi,MassesPlot,GammasPlot,SpinsPlot,asPlot,bsPlot,psi_nS,dRadB0,dRadKs);
    matrixPlot->setData(&currData);

    matrixPlot->copyParams();
    compsIntegral = matrixPlot->normalise();

    fractions.push_back(compsIntegral/totalIntegral);

    cout<<"  Component "<<kStarNames[kCounter]<<" normalisation factor : "<<matrixPlot->normalise()<<" (fraction: "<<compsIntegral/totalIntegral*100.0<<"%)"<<endl;

    matrixPlot->getCompProbsAtDataPoints(pdfCompValues);

    matrixPlot->clearCurrentFit();

    for (int k = 0; k<pdfCompValues[0].size();k++) {
      //std::cout<<" Bin : "<< k << " pdf : " << pdfCompValues[0][k] <<std::endl;
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

    ////////////////////////////////////////////////////////////////////////////////
    //Filling vectors

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

    ////////////////////////////////////////////////////////////////////////////////
    //Filling components projections graphs

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


    /*massKPiHisto.Draw("");
      signalCompPlot->Draw("sameL");

      sprintf(bufferstring,"plots/plot%d.eps",kCounter);
      canvas->SetLogy(1);
      canvas->SaveAs(bufferstring);
      canvas->Clear();*/
    ++kCounter;

    MassesPlot.clear();
    GammasPlot.clear();
    SpinsPlot.clear();
    asPlot.clear();
    bsPlot.clear();
    pdfCompValues.clear();

  }

  //Adding single points to plot better and total plots

  fptype pointsX[2], pointsY[2];
  pointsX[0] = massKPi->lowerlimit; pointsX[1] = massKPi->upperlimit;
  pointsY[0] = 0.01; pointsY[1] = massKPiHisto.GetMaximum();
  TGraph* pointsMKP = new TGraph(2,pointsX,pointsY);
  multiGraphMKPi->Add(&signalTotalPlotMKPi,"L");
  multiGraphMKPi->Add(pointsMKP,"P");

  pointsX[0] = massPsiPi->lowerlimit; pointsX[1] = massPsiPi->upperlimit;
  pointsY[0] = massPsiPiHisto.GetMinimum();pointsY[1] = massPsiPiHisto.GetMaximum();
  TGraph* pointsCKS = new TGraph(2,pointsX,pointsY);
  multiGraphmassPsiPi->Add(&signalTotalPlotmassPsiPi,"L");
  multiGraphmassPsiPi->Add(pointsCKS,"P");

  pointsX[0] = cosMuMu->lowerlimit; pointsX[1] = cosMuMu->upperlimit;
  pointsY[0] = cosMuMuHisto.GetMinimum(); pointsY[1] = cosMuMuHisto.GetMaximum();
  TGraph* pointsCMM = new TGraph(2,pointsX,pointsY);
  multiGraphCosMuMu->Add(&signalTotalPlotCosMuMu,"L");
  multiGraphCosMuMu->Add(pointsCMM,"P");

  pointsX[0] = phi->lowerlimit; pointsX[1] = phi->upperlimit;
  pointsY[0] = phiHisto.GetMinimum(); pointsY[1] = phiHisto.GetMaximum();
  TGraph* pointsPHI = new TGraph(2,pointsX,pointsY);
  multiGraphPhi->Add(&signalTotalPlotPhi,"L");
  multiGraphPhi->Add(pointsPHI,"P");

  ////////////////////////////////////////////////////////////////////////////////
  //PLOTTING


  //Mass K Pi
  multiGraphMKPi->Draw("AL");
  massKPiHisto.Draw("same");
  legPlot->Draw(); fitStat->Draw();
  canvas->SetLogy(1);

  canvas->SaveAs(TString::Format("%s/%s%s.%s",plotsDir.Data(),massKPi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->Clear();
  ////////////////////////////////////////////////////////////////////////////////

  //CosMuMu
  multiGraphCosMuMu->Draw("AL");
  cosMuMuHisto.Draw("same");
  // it's enough on the m(KPi) plot
  //legPlot->Draw(); //fitStat->Draw();
  canvas->SetLogy(1);

  canvas->SaveAs(TString::Format("%s/%s%s.%s",plotsDir.Data(),cosMuMu_name.Data(),plotsName.Data(),extension.Data()));
  canvas->Clear();
  ////////////////////////////////////////////////////////////////////////////////

  //massPsiPi
  multiGraphmassPsiPi->Draw("AL");
  massPsiPiHisto.Draw("same");
  // it's enough on the m(KPi) plot
  //legPlot->Draw(); //fitStat->Draw();
  canvas->SetLogy(1);

  canvas->SaveAs(TString::Format("%s/%s%s.%s",plotsDir.Data(),massPsiPi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->Clear();
  ////////////////////////////////////////////////////////////////////////////////

  //Phi
  multiGraphPhi->Draw("AL");
  phiHisto.Draw("same");
  // it's enough on the m(KPi) plot
  //legPlot->Draw(); //fitStat->Draw();
  canvas->SetLogy(1);

  canvas->SaveAs(TString::Format("%s/%s%s.%s",plotsDir.Data(),phi_name.Data(),plotsName.Data(),extension.Data()));
  canvas->Clear();
  ////////////////////////////////////////////////////////////////////////////////

  cout<<endl;
  cout<<"~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~"<<endl;
  cout << "PDF fitting time:       " << (fitClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout << "Data plotting time:     " << (dataSetClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout << "PDF sum time:           " << (sumClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout << "PDF normalization time: " << (normClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout << "PDF projection time:    " << (projClocks / CLOCKS_PER_SEC) << " s" << endl ;
  cout<<"~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~"<<endl;


  /*
    matrix->getCompProbsAtDataPoints(pdfTotalValuesNorm,events);
    for (int j = 0; j < massKPi->numbins; ++j) {
    projMKPiHisto.SetBinContent(j+1,mkpTotalProjection[j]);
    std::cout<<" Bin "<<j<<" center = "<<projMKPiHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;

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

    std::cout<< "Pdf value : "<<tempValues[0][0]<<std::endl;
  */

  return 0;

}
