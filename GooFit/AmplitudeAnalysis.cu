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
#include <sstream>
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

void debug(int line){

  #ifdef CUDADEBUGGING
  std::cout<<"Debugging on line "<<line<<std::endl;
  #endif

}

int parDotSpin (fptype dotSpin){

  int result = static_cast<int>(floor((dotSpin - floor(dotSpin))*10.+.1));
  return result;
}

std::string doubleToStr (fptype dbl){

  std::ostringstream strs;
  strs << dbl;
  return strs.str();
}

void printinstruction(){

  std::cerr << "======= Instructions \n"
  			<< "\t-h,--help \t\t\t\t Show this help message\n"
  			<< "\t-n  <events> \t\t\t\t Specify the number of events\n"
  			<< "\t-r  <path> \t\t\t\t Read Generated Events from txt in <path>\n"
  			<< "\t-b1 <b1> \t\t\t\t Select binning for MassKPi (for normalisation, default : 40)\n"
        << "\t-b2 <b2> \t\t\t\t Select binning for CosMuMu (for normalisation, default : 40)\n"
        << "\t-b3 <b3> \t\t\t\t Select binning for CosKStar (for normalisation, default : 40)\n"
        << "\t-b4 <b4> \t\t\t\t Select binning for Phi (for normalisation, default : 40)\n"
        << "\t-p1  <p> \t\t\t\t Select plotting binning finenness (default : 30) for MassKPi \n"
        << "\t-p2  <p> \t\t\t\t Select plotting binning finenness (default : 30) for CosMuMu \n"
        << "\t-p3  <p> \t\t\t\t Select plotting binning finenness (default : 30) for CosKStar \n"
        << "\t-p4  <p> \t\t\t\t Select plotting binning finenness (default : 30) for Phi \n"
        << "\t-Bb <Bb> \t\t\t\t Select bound limits for b parameter (default : 9999)\n"
        << "\t-k800 \t\t\t\t Add K*(800) to p.d.f.\n"
        << "\t-k892 \t\t\t\t Add K*(892) to p.d.f.\n"
        << "\t-k1410 \t\t\t\t Add K*(1410) to p.d.f.\n"
        << "\t-k1430_0 \t\t\t\t Add K*(1430_0) to p.d.f.\n"
        << "\t-k1430_2 \t\t\t\t Add K*(1430_2) to p.d.f.\n"
        << "\t-k1780 \t\t\t\t Add K*(1430_3) to p.d.f.\n"
        //<< "\t-e <number of events> \t\t\t Select sigma cuts\n"
        <<std::endl;

}


int main(int argc, char** argv)
{
  debug(__LINE__);

  char bufferstring[1024];

  unsigned int events = 10000;
  unsigned int nOfKstar = 0;

  unsigned int bin1 = 40;
  unsigned int bin2 = 40;
  unsigned int bin3 = 40;
  unsigned int bin4 = 40;

  unsigned int datapoints = 100;

  unsigned int plottingfine1 = 30;
  unsigned int plottingfine2 = 30;
  unsigned int plottingfine3 = 30;
  unsigned int plottingfine4 = 30;

  fptype aMax = +9999.;
  fptype bMax = +9999.;

  bool k892Star = false;
  bool k800Star = false;
  bool k1410Star = false;
  bool k1430Star0 = false;
  bool k1430Star2 = false;
  bool k1780Star = false;

  std::string datasetname = "dataset";
  std::string plotsname = "./plots/plot";
  std::vector< std::string> kStarNames;

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
      else if (arg == "-dP")
  		{
        if (i + 1 < argc) // Make sure we aren't at the end of argv!
  			{
  				i++;
  				std::istringstream ss(argv[i]);
  				if (!(ss >> datapoints))
  				{
  					std::cerr << "Invalid number " << argv[i] << '\n';
  					exit(1);
  				}
  			}
  		}


  }

  if(!(k1430Star0 || k1430Star2 || k1410Star || k800Star ||k892Star || k1780Star)){

    cout<<"No K* selected (K892,K800,K1410,K1430) please see instructions below"<<endl;
    printinstruction();
    return 1;

  }else{

    cout<<" Starting Amplitude Analysis Fit With : "<<nOfKstar<<" K* "<<endl;
    if(k892Star) {
      cout<<" - K*(892)"<<endl;
      datasetname +="k892";
      plotsname +="k892";
      kStarNames.push_back("K*(892)");
    }
    if(k800Star) {
      cout<<" - K*(800)"<<endl;
      datasetname +="k800";
      plotsname +="k800";
      kStarNames.push_back("K*(800)");
    }
    if(k1410Star) {
      cout<<" - K*(1410)"<<endl;
      datasetname +="k1410";
      plotsname +="k1410";
      kStarNames.push_back("K*(k1410)");
    }
    if(k1430Star0) {
      cout<<" - K*(1430_0)"<<endl;
      datasetname +="k14300";
      plotsname +="k14300";
      kStarNames.push_back("K*(1430_0)");
    }
    if(k1430Star2) {
      cout<<" - K*(1430_2)"<<endl;
      datasetname +="k14302";
      plotsname +="k14302";
      kStarNames.push_back("K*(1430_2)");}
    if(k1780Star) {
      cout<<" - K*(1780_3)"<<endl;
      datasetname +="k17803";
      plotsname +="k17803";
      kStarNames.push_back("K*(k1780_3)");}
}

  datasetname += ".txt";

  fptype aMin = -aMax;
  fptype bMin = -bMax;

  debug(__LINE__);

  //CANVAS
  TCanvas* canvas = new TCanvas("","",1600,1200);

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
	for(int k=0;k<BINS;k++){

    pdfBkgHist.SetBinContent(k+1,ValsFondo[k]);
		totalFondo += ValsFondo[k];
    std::cout<< ValsFondo[k]<<std::endl;

  }
  debug(__LINE__);
	for(int k=0;k<BINS;k++){
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

 Variable* massKPi = new Variable("massKPi",1.,0.6,2.2);
 Variable* cosMuMu = new Variable("cosMuMu",0.,-1,1); // cosine of the psi(nS) helicity angle
 Variable* cosKstar = new Variable("cosKstar",0.,-1,1); // cosine of the K* helicity angle
 Variable* phi = new Variable("phi",0.25,-TMath::Pi(),TMath::Pi());
 massKPi->numbins = bin1;
 cosMuMu->numbins = bin2;
 cosKstar->numbins = bin3;
 phi->numbins = bin4;

 GooPdf* phaseSpace = new ThreeBodiesPsiPiK("phasespace",massKPi,&mBd,&mPion,&mKaon,&mMuMu);

 //fptype ratio = ((fptype)(plottingfine))/((fptype)massKPi->numbins);
 fptype ratio = ((fptype)(plottingfine1))/((fptype)datapoints);
 std::vector<Variable* > obserVariables;

 obserVariables.push_back(massKPi);
 obserVariables.push_back(cosMuMu);
 obserVariables.push_back(cosKstar);
 obserVariables.push_back(phi);

 std::vector<Variable*> Masses;
 std::vector<Variable*> Gammas;
 std::vector<Variable*> Spins;
 std::vector<Variable*> as;
 std::vector<Variable*> bs;

 if(k892Star){

   cout <<"Adding K*(892)..." <<endl;

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

 if(k800Star){

   cout <<"Adding K*(800)..." <<endl;

   Masses.push_back(new Variable("K_800_Mass_0",M800));
   Gammas.push_back(new Variable("K_800_Gamma_0",G800));
   Spins.push_back(new Variable("K_800_Spin_0",0.0));
   as.push_back(new Variable("a_K_800_0",1.12,aMin,aMax) );
   bs.push_back(new Variable("b_K_800_0",2.3,bMin,bMax) );

 }

 if(k1410Star){

   cout <<"Adding K*(1410)..." <<endl;

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

 if(k1430Star0){

   cout <<"Adding K*(1430_0)..." <<endl;

   Masses.push_back(new Variable("K_1430_0_Mass_0",M1430_0));
   Gammas.push_back(new Variable("K_1430_0_Gamma_0",G1430_0));
   Spins.push_back(new Variable("K_1430_0_Spin_0",0.0));
   as.push_back(new Variable("a_K_1430_0_0",0.89,aMin,aMax) );
   bs.push_back(new Variable("b_K_1430_0_0",-2.17,bMin,bMax) );

 }

 if(k1430Star2){

   cout <<"Adding K*(1430_2)..." <<endl;

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
 if(k1780Star){

   cout <<"Adding K*(1780)_3..." <<endl;

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


 GooPdf* matrix = new MatrixPdf("Kstars_signal", massKPi, cosMuMu, cosKstar, phi,Masses,Gammas,Spins,as,bs,psi_nS,dRadB0,dRadKs);

 UnbinnedDataSet dataset(obserVariables);

 TH1F mKPHisto("mKPHisto", "mKPHisto",datapoints, massKPi->lowerlimit, massKPi->upperlimit);

 //ifstream dataTxt("dataset.txt");
 ifstream dataTxt(datasetname.c_str());
 if(!(dataTxt.good())){
   std::cout<<" No valid input at : "<<datasetname<<" provided. Returning."<<std::endl;
   return 1;
 }

 fptype var1, var2, var3, var4;

 if(dataTxt.good()){
    while(dataTxt>>var1>>var2>>var3>>var4){

      massKPi->value = var1;
      cosMuMu->value = var2;
      cosKstar->value = var3;
      phi->value = var4;

      //std::cout<< massKPi->value << " - " <<cosMuMu->value << " - " << cosKstar->value << " - " << phi->value << " - " << std::endl;
      //std::cout<< massKPi->value << " " <<cosMuMu->value << " " << cosKstar->value << " " << phi->value << " " << std::endl;

      dataset.addEvent();
      mKPHisto.Fill(massKPi->value);

    }

  }
  dataTxt.close();

  matrix->setData(&dataset);
  //total->setData(&dataset);

	FitManager fitter(matrix);
  //FitManager fitter(total);
  fitter.fit();
	fitter.getMinuitValues();

  stopC = times(&stopProc);
  gettimeofday(&stopTime, NULL);

  fptype fitClocks = (stopC - startC)*10000.;

  gettimeofday(&startTime, NULL);
  startC = times(&startProc);

  UnbinnedDataSet currData(obserVariables);
  std::vector<UnbinnedDataSet> compData;
  std::vector<std::vector<fptype> > pdfTotalValues;
  /*std::vector<std::vector<std::vector<fptype> >  pdfCompValues;*/
  std::vector<std::vector<fptype> > pdfCompValues;

  for (int id = 0; id < nOfKstar; ++id) {
    compData.push_back( UnbinnedDataSet(obserVariables));
  }

  std::vector<fptype> fractions;

  std::vector<fptype> compEvents;


  std::vector<fptype> mkpTotalProjection;

  massKPi->numbins = plottingfine1;
  cosMuMu->numbins = plottingfine2;
  cosKstar->numbins = plottingfine3;
  phi->numbins = plottingfine4;

  fptype pointsXTot[massKPi->numbins];
	fptype pointsYTot[massKPi->numbins];

  TH1F plotHisto("plotHisto", "plotHisto",massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit);


  for (size_t i = 0; i < massKPi->numbins; i++) {
    mkpTotalProjection.push_back(0.0);
  }

  fptype sum = 0.0;

  std::cout<<"- Starting plotting cycle "<<std::endl;
  std::cout<<"- Plotting dataset generation "<<std::endl;
  for (int k = 0; k < phi->numbins; ++k) {
      phi->value = phi->lowerlimit + (phi->upperlimit - phi->lowerlimit)*(k + 0.5) / phi->numbins;
      //std::cout<<"Phi : "<< k <<std::endl;
      for (int j = 0; j < cosMuMu->numbins; ++j) {
          cosMuMu->value = cosMuMu->lowerlimit + (cosMuMu->upperlimit - cosMuMu->lowerlimit)*(j + 0.5) / cosMuMu->numbins;
          //std::cout<<"CosMu : "<< j <<std::endl;
          for (int a = 0; a < cosKstar->numbins; ++a) {
              cosKstar->value = cosKstar->lowerlimit + (cosKstar->upperlimit - cosKstar->lowerlimit)*(a + 0.5) / cosKstar->numbins;
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
                  std::cout<<cosKstar->value<<" ";
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

 stopC = times(&stopProc);
 gettimeofday(&stopTime, NULL);

 fptype dataSetClocks = (stopC - startC)*10000.;

 gettimeofday(&startTime, NULL);
 startC = times(&startProc);
 ////////////////////////////////////////////////////////////////////////////////
 ///// TOTAL PDF PLOT
 ////////////////////////////////////////////////////////////////////////////////

 TLegend *legPlot = new TLegend(0.65,0.80,0.95,0.99,"Legend");
 TPaveText *pt = new TPaveText(0.65, .45, .95, .80, "NDC");

 std::cout<<"- Evaluating the total P.d.f."<<std::endl;
 matrix->setData(&currData);
 matrix->getCompProbsAtDataPoints(pdfTotalValues);

 //std::cout<<" Vector size : "<<pdfTotalValues[0].size()<<std::endl;
 //std::cout<<" Vector proj : "<<pdfTotalValues[0].size()/massKPi->numbins<<std::endl;

 int notMPKBins = pdfTotalValues[0].size()/massKPi->numbins;

 for (int k = 0; k<pdfTotalValues[0].size();k++){
   //std::cout<<mkpTotalProjection[k]*events/sum<<std::endl;
   sum += pdfTotalValues[0][k];
 }

 stopC = times(&stopProc);
 gettimeofday(&stopTime, NULL);

 fptype sumClocks = (stopC - startC)*10000.;

 gettimeofday(&startTime, NULL);
 startC = times(&startProc);

 std::cout<<"[ Total Pdf sum : "<<sum<<" ] "<<std::endl;

 for (int k = 0; k<pdfTotalValues[0].size();k++){
   //std::cout<<mkpTotalProjection[k]*events/sum<<std::endl;
   pdfTotalValues[0][k] /=sum;
   pdfTotalValues[0][k] *= 10000;
 }

 stopC = times(&stopProc);
 gettimeofday(&stopTime, NULL);
 startC = times(&startProc);

 fptype normClocks = (stopC - startC)*10000.;


 for (int j = 0; j < massKPi->numbins; ++j) {
   for (int i = 0; i < notMPKBins; ++i) {
     mkpTotalProjection[j] += pdfTotalValues[0][j+i*massKPi->numbins];
   }
 }

 stopC = times(&stopProc);
 gettimeofday(&stopTime, NULL);

 fptype projClocks = (stopC - startC)*10000.;

  for (int j = 0; j < massKPi->numbins; ++j) {
    plotHisto.SetBinContent(j+1,mkpTotalProjection[j]);
    //std::cout<<" Bin "<<j<<" center = "<<plotHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;

  }

  plotHisto.Scale(ratio);

  for(int k=0;k<massKPi->numbins;k++){
    pointsXTot[k] = plotHisto.GetBinCenter(k+1);
    pointsYTot[k] = plotHisto.GetBinContent(k+1);
  }

  TMultiGraph* multiGraph = new TMultiGraph();
  multiGraph->SetTitle("K* plot");

  TGraph signalTotalPlot(massKPi->numbins,pointsXTot,pointsYTot);

  signalTotalPlot.SetLineColor(kRed);
  signalTotalPlot.SetLineWidth(2);

  fptype totalIntegral = matrix->normalise();
  fptype compsIntegral = 0.0;
  std::cout<<"Total Normalisation Factor = "<<totalIntegral<<std::endl;

  int kCounter = 0;

  for (size_t u = 0; u < as.size(); ++u) {

      if(Spins[u]->value==0.0){

        pt->AddText(TString::Format(" ---- %s", kStarNames[kCounter].c_str()));
        pt->AddText(TString::Format(" a_{0} = %.3f #pm %.3f b_{0} = %.3f #pm %.3f",as[u]->value,as[u]->error,bs[u]->value,bs[u]->error));

      }else{

        pt->AddText(TString::Format(" ---- %s", kStarNames[kCounter].c_str()));
        pt->AddText(TString::Format(" a_{0} = %.3f #pm %.3f b_{0} = %.3f #pm %.3f",as[u]->value,as[u]->error,bs[u]->value,bs[u]->error));
        pt->AddText(TString::Format(" a_{p1} = %.3f #pm %.3f b_{p1} = %.3f #pm %.3f",as[u+1]->value,as[u+1]->error,bs[u+1]->value,bs[u+1]->error));
        pt->AddText(TString::Format(" a_{m1} = %.3f #pm %.3f b_{m1} = %.3f #pm %.3f",as[u+2]->value,as[u+2]->error,bs[u+2]->value,bs[u+2]->error));

        u+=2;
      }

      ++kCounter;
  }

  pt->SetTextAlign(12);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);

  matrix->clearCurrentFit();

  legPlot->AddEntry(&mKPHisto,"Generated Data","f");
  legPlot->AddEntry(&signalTotalPlot,"Total Fit","l");

  //multiGraph->Add(&signalTotalPlot,"L");
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

  std::vector<TH1F*> compHistos;

  for (int k = 0; k < (int)as.size(); ++k) {

    std::vector<fptype> mkpCompProjection;

    for (size_t i = 0; i < massKPi->numbins; i++) {
      mkpCompProjection.push_back(0.0);
    }

    fptype pointsXComp[massKPi->numbins];
  	fptype pointsYComp[massKPi->numbins];

    //TH1F tempPlotHisto("tempPlotHisto", "tempPlotHisto",massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit);
    sprintf(bufferstring,"comp_%d_plotHisto",kCounter);
    compHistos.push_back(new TH1F(bufferstring,bufferstring,massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit));


    cout<<"- Plotting component "<<kStarNames[kCounter]<<" all other components to zero"<<endl;
    sum = 0.0;

/*
    for (size_t i = 0; i < massKPi->numbins; i++) {
      mkpCompProjection[k].push_back(0.0);
    }*/

    ///// OTHER COMPONENTS TO ZERO

    for (int l = 0; l < as.size(); ++l) {
      as[l]->fixed = true;
      bs[l]->fixed = true;
    }

    if(Spins[k]->value==0.0){

      as[k]->fixed;
      bs[k]->fixed;

      MassesPlot.push_back(Masses[k]);
      GammasPlot.push_back(Gammas[k]);
      SpinsPlot.push_back(Spins[k]);
      asPlot.push_back(as[k]);
      bsPlot.push_back(bs[k]);


      for (int j = 0; j < (int)as.size(); ++j) {
        if(j!=k){

          MassesPlot.push_back(Masses[j]);
          GammasPlot.push_back(Gammas[j]);
          SpinsPlot.push_back(Spins[j]);
          asPlot.push_back(new Variable("zero_a",0.0));
          bsPlot.push_back(new Variable("zero_b",0.0));


        }}

    }else{
      for (int i = k; i <= k+2; ++i) {


        as[i]->fixed;
        bs[i]->fixed;

        MassesPlot.push_back(Masses[i]);
        GammasPlot.push_back(Gammas[i]);
        SpinsPlot.push_back(Spins[i]);
        asPlot.push_back(as[i]);
        bsPlot.push_back(bs[i]);

      }
      for(int d = 0; d < as.size(); ++d) {
        if(d!=k && d!=k+1 && d!=k+2) {

          MassesPlot.push_back(Masses[d]);
          GammasPlot.push_back(Gammas[d]);
          SpinsPlot.push_back(Spins[d]);
          asPlot.push_back(new Variable("zero_a",0.0));
          bsPlot.push_back(new Variable("zero_b",0.0));


        }}
        k+=2;
      }

      sprintf(bufferstring,"Kstars_signal_plot_%d",kCounter);
      GooPdf* matrixPlot = new MatrixPdf(bufferstring, massKPi, cosMuMu, cosKstar, phi,MassesPlot,GammasPlot,SpinsPlot,asPlot,bsPlot,psi_nS,dRadB0,dRadKs);
      matrixPlot->setData(&currData);


      //FitManager fitterPlot(matrixPlot);
      //fitterPlot.fit();
      //fitterPlot.getMinuitValues();
      matrixPlot->copyParams();
      compsIntegral = matrixPlot->normalise();

      fractions.push_back(compsIntegral/totalIntegral);

      cout<<"  Component "<<kStarNames[kCounter]<<" normalisation factor : "<<matrixPlot->normalise()<<" (fraction: "<<compsIntegral/totalIntegral*100.0<<"%)"<<endl;

      matrixPlot->getCompProbsAtDataPoints(pdfCompValues);

      matrixPlot->clearCurrentFit();

      for (int k = 0; k<pdfCompValues[0].size();k++){
        //std::cout<<" Bin : "<< k << " pdf : " << pdfCompValues[0][k] <<std::endl;
        sum += pdfCompValues[0][k];
      }

      for (int k = 0; k<pdfCompValues[0].size();k++){
        pdfCompValues[0][k] /=sum;
        pdfCompValues[0][k] *= events;
        pdfCompValues[0][k] *= (compsIntegral/totalIntegral);
        //compHistos[kCounter]->SetBinContent(k,pdfCompValues[0][k]);
      }

      for (int j = 0; j < massKPi->numbins; ++j) {
        for (int i = 0; i < notMPKBins; ++i) {
          mkpCompProjection[j] += pdfCompValues[0][j+i*massKPi->numbins];
        }
        compHistos[kCounter]->SetBinContent(j+1,mkpCompProjection[j]);
      }

      compHistos[kCounter]->Scale(ratio);
      pdfCompValues.clear();

      for(int k=0;k<massKPi->numbins;k++){

        pointsXComp[k] = compHistos[kCounter]->GetBinCenter(k+1);
        pointsYComp[k] = compHistos[kCounter]->GetBinContent(k+1);

        //std::cout<<" X : "<< pointsXComp[k] << " Y : " << pointsYComp[k] <<std::endl;
      }

      TGraph* signalCompPlot = new TGraph(massKPi->numbins,pointsXComp,pointsYComp);

      signalCompPlot->SetLineColor(kCounter+3);
      signalCompPlot->SetLineWidth(2);
      signalCompPlot->SetLineStyle(kDashed);
      signalCompPlot->GetXaxis()->SetTitle("m(K#Pi)");
      sprintf(bufferstring,"Events / (%.3f)",(massKPi->upperlimit- massKPi->lowerlimit)/massKPi->numbins);
      signalCompPlot->GetYaxis()->SetTitle(bufferstring);
      //signalCompPlot->Draw("");
      multiGraph->Add(signalCompPlot,"L");

      sprintf(bufferstring,"%s (%.2f %)",kStarNames[kCounter].c_str(),compsIntegral/totalIntegral*100.0);
      legPlot->AddEntry(signalCompPlot,bufferstring,"l");


      /*mKPHisto.Draw("");
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

  multiGraph->Add(&signalTotalPlot,"L");

  multiGraph->Draw("AL");
  mKPHisto.Draw("same");
  //signalTotalPlot.Draw("sameL");
  legPlot->Draw();
  pt->Draw();

  canvas->SetLogy(1);

  plotsname += ".eps";
  canvas->SaveAs(plotsname.c_str());
  //canvas->SaveAs("plots/plot.png");

  cout<<endl;
  cout<<"~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~"<<endl;
  cout << "Fit time:       " << (fitClocks / CLOCKS_PER_SEC) << endl ;
  cout << "Data time:      " << (dataSetClocks / CLOCKS_PER_SEC) << endl ;
  cout << "Sum time:       " << (sumClocks / CLOCKS_PER_SEC) << endl ;
  cout << "Norm time:      " << (normClocks / CLOCKS_PER_SEC) << endl ;
  cout << "Proj time:      " << (projClocks / CLOCKS_PER_SEC) << endl ;
  cout<<"~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~"<<endl;


/*
  matrix->getCompProbsAtDataPoints(pdfTotalValuesNorm,events);
  for (int j = 0; j < massKPi->numbins; ++j) {
    plotHisto.SetBinContent(j+1,mkpTotalProjection[j]);
    std::cout<<" Bin "<<j<<" center = "<<plotHisto.GetBinCenter(j+1)<<" : "<<mkpTotalProjection[j]<<std::endl;

  }*/



/*
 UnbinnedDataSet tempData(obserVariables);

 std::vector<std::vector<fptype> > tempValues;

 massKPi->value = 1.0;
 cosKstar->value = 0.5;
 cosMuMu->value = 0.5;
 phi->value = 0.25;

 tempData.addEvent();
 matrix->getCompProbsAtDataPoints(tempValues);

 std::cout<< "Pdf value : "<<tempValues[0][0]<<std::endl;*/

 return 0;

}
