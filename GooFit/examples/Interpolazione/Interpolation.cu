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

#include "TNtupleD.h"
#include "TTree.h"
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
#include "TGraph2D.h"

#include "BiDimHistoPdf.hh"

#include "utilities.h"
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


int main(int argc, char** argv)

{

  TCanvas* canvas = new TCanvas("","",2000,1200);

  TString path = "/lustre/home/adrianodif/RootFiles/Z4430/";
  TString fileName = "TMVApp_data_withBDTCutAt0p00_JPsi_2p0Sig_6p0to9p0SB.root";

  TFile *file = TFile::Open(path+fileName);

  if (!(file)) {
    std::cout<<"File named \'"<<file->GetName()<<"\' NOT FOUND.\nReturning."<<std::endl;
    return -1;
  }

  TString thName = "psi2SPi_vs_KPi_masses_sbs_BDT";

  TH2F* th = (TH2F*)file->Get(thName);
  //TH2F* thOut = (TH2F*)file->Get(thName);

  if (!(th)) {
    std::cout<<"TH2 named \'"<<thName<<"\' NOT FOUND in found in TFile \'" <<file->GetName() <<"\'.\nReturning."<<std::endl;
    return -1;
  }

  th->Draw();

  fptype massKPi_min = 0.6, massKPi_max = 2.2;
  fptype massPsiPi_min = 3.2, massPsiPi_max = 4.9;

  TString massKPi_name = "massKPi",massPsiPi_name = "massPsiPi";

  Variable* massKPi = new Variable(massKPi_name.Data(),1.,massKPi_min,massKPi_max); massKPi->numbins = th->GetNbinsX();
  Variable* massPsiPi = new Variable(massPsiPi_name.Data(),TMath::Sqrt(23),massPsiPi_min,massPsiPi_max); massPsiPi->numbins = th->GetNbinsY();

  std::vector<Variable*> obserVariables;
  obserVariables.push_back(massKPi);
  obserVariables.push_back(massPsiPi);

  std::vector<Variable*> useless;
  useless.push_back(new Variable("zeroa",0.0));
  useless.push_back(new Variable("zerob",0.0));

  BinnedDataSet* biDataset = new BinnedDataSet(obserVariables,"Bkg Dataset Masses");

  for (int j = 0; j < biDataset->getNumBins(); ++j)
    biDataset->setBinContent(j,0.0);

  // FILLING DATASET WITH HISTOGRAM
  for (int j = 0; j < biDataset->getNumBins(); ++j) {
    biDataset->setBinContent(j, th->GetBinContent(th->FindBin(biDataset->getBinCenter(massKPi,j),biDataset->getBinCenter(massPsiPi,j))));
  }

  for (int j = 0; j < biDataset->getNumBins(); ++j) {

    if (biDataset->getBinContent(j) != th->GetBinContent(th->FindBin(biDataset->getBinCenter(massKPi,j),biDataset->getBinCenter(massPsiPi,j))))
      std::cout<<"Not the same "<<std::endl;

  }

  std::cout<<"Here"<<std::endl;
  //canvas->SaveAs("th2Test.png");

  UnbinnedDataSet plottingGridData(obserVariables);

  int events = biDataset->getNumEvents();

  GooPdf* biDimPdf = new BiDimHistoPdf ("EfficiencyPdf",biDataset,obserVariables);
  std::cout<<"Here"<<std::endl;
  massKPi->value = 1.4;
  massPsiPi->value = 4.0;

  int j = biDataset->getBinNumber();

  fptype massKPiCenter = biDataset->getBinCenter(massKPi,j);
  fptype massPsiPiCenter = biDataset->getBinCenter(massPsiPi,j);

  massKPi->value = massKPiCenter;
  massPsiPi->value =massPsiPiCenter;

  // biDimPdf->setData(&plottingGridData);
  biDimPdf->setData(biDataset);

  std::cout<<"InterPolation : "<<biDimPdf->getValue()<<"Histo : "<<th->GetBinContent(th->FindBin(biDataset->getBinCenter(massKPi,j),biDataset->getBinCenter(massPsiPi,j)))<<std::endl;

  std::vector<fptype> pdfTotalValues;
  std::vector<std::vector<fptype> > pdfTotalValuesPre;

  // biDimPdf->getCompProbsAtDataPoints(pdfTotalValuesPre);
  //
  // std::cout<<pdfTotalValuesPre[0].size()<<std::endl;
  // std::cout<<pdfTotalValuesPre.size()<<std::endl;
  // for (int k = 0; k<pdfTotalValuesPre[0].size();k++)
  //   std::cout <<" Bin : " << k << " pdf : " << pdfTotalValuesPre[0][k] <<std::endl;


  int plottingfine1 = massKPi->numbins;
  int plottingfine2 = massPsiPi->numbins;

  for (int a = 0; a < massPsiPi->numbins; ++a) {
         massPsiPi->value = massPsiPi->lowerlimit + (massPsiPi->upperlimit - massPsiPi->lowerlimit)*(a + 0.5) / massPsiPi->numbins;
         for (int i = 0; i < massKPi->numbins; ++i) {
         massKPi->value = massKPi->lowerlimit + (massKPi->upperlimit - massKPi->lowerlimit)*(i + 0.5) / massKPi->numbins;
         plottingGridData.addEvent();
         pdfTotalValues.push_back(biDimPdf->getValue());

      }
  }



  fptype ratioMKPi = ((fptype)(plottingfine2))/((fptype)massKPi->numbins);
  fptype ratioMassPsiPi = ((fptype)(plottingfine1))/((fptype)massPsiPi->numbins);
  std::cout<<"Ratio "<<ratioMKPi<<std::endl;
  massKPi->numbins = plottingfine1;
  massPsiPi->numbins = plottingfine2;

  TH1F projMKPiHisto("projMKPiHisto", "projMKPiHisto",massKPi->numbins, massKPi->lowerlimit, massKPi->upperlimit);
  TH1F projmassPsiPiHisto("projmassPsiPiHisto", "projmassPsiPiHisto",massPsiPi->numbins, massPsiPi->lowerlimit, massPsiPi->upperlimit);

  fptype pointsmassPsiPiXTot[massPsiPi->numbins];
  fptype pointsmassPsiPiYTot[massPsiPi->numbins];

  fptype pointsMKPiXTot[massKPi->numbins];
  fptype pointsMKPiYTot[massKPi->numbins];

  //biDimPdf->getCompProbsAtDataPoints(pdfTotalValues);

  // std::cout<<pdfTotalValues[0].size()<<std::endl;
  // std::cout<<pdfTotalValues.size()<<std::endl;

  fptype sum = 0.0;

  for (int k = 0; k<pdfTotalValues.size();k++) {
    //std::cout <<" Bin : " << k << " pdf : " << pdfCompValues[0][k] <<std::endl;
    sum += pdfTotalValues[k];
  }

  for (int k = 0; k<pdfTotalValues.size();k++) {
    pdfTotalValues[k] /=sum;
    pdfTotalValues[k] *= (fptype) events;
    //compHistosMKPi[kCounter]->SetBinContent(k,pdfCompValues[0][k]);
  }

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

  // int notMPKBins = pdfTotalValues.size()/massKPi->numbins;


  std::vector<fptype> mkpTotalProjection;
  std::vector<fptype> massPsiPiTotalProjection;

  for (int j = 0; j < biDataset->getNumBins(); ++j) {
    massKPi->value = biDataset->getBinCenter(massKPi,j);
    massPsiPi->value = biDataset->getBinCenter(massPsiPi,j);
    // std::cout<<"Bin : "<<j<<" MassKPi : "<<biDataset->getBinCenter(massKPi,j)<<" MassPsiPi : "<<biDataset->getBinCenter(massPsiPi,j)<<" - con : "<<th->GetBinContent(th->FindBin(biDataset->getBinCenter(massKPi,j),biDataset->getBinCenter(massPsiPi,j)))<<" - pdf : "<<biDimPdf->getValue()<<std::endl;

  }

  //
  //
  // for (int i = 0; i < massKPi->numbins; ++i)
  //     mkpTotalProjection.push_back(0.0);
  //
  // for (int i = 0; i < massPsiPi->numbins; ++i)
  //     massPsiPiTotalProjection.push_back(0.0);
  //
  // for (int j = 0; j < massKPi->numbins; ++j)
  //   for (int i = 0; i < massPsiPi->numbins; ++i)
  //     mkpTotalProjection[j] += pdfTotalValues[j  +  i * massKPi->numbins];
  //
  //
  // for (int j = 0; j < massPsiPi->numbins; ++j)
  //     for (int i = 0; i < massKPi->numbins; ++i)
  //       massPsiPiTotalProjection[j] += pdfTotalValues[i  +  j * massKPi->numbins];
  //
  // for (int j = 0; j < massKPi->numbins; ++j)
  //         projMKPiHisto.SetBinContent(j+1,mkpTotalProjection[j]);
  //
  // for (int j = 0; j < massPsiPi->numbins; ++j)
  //           projmassPsiPiHisto.SetBinContent(j+1,massPsiPiTotalProjection[j]);
  //
  // projMKPiHisto.Scale(ratioMKPi);
  // projmassPsiPiHisto.Scale(ratioMassPsiPi);

  for (int j = 0; j < massKPi->numbins; ++j) {
    pointsMKPiXTot[j] = projMKPiHisto.GetBinCenter(j+1);
    fptype yValue = 0.0;
    massKPi->value = projMKPiHisto.GetBinCenter(j+1);
    for (int i = 0; i < massPsiPi->numbins; ++i) {
      massPsiPi->value = projmassPsiPiHisto.GetBinCenter(i+1);
      yValue += biDimPdf->getValue();
    }
    projMKPiHisto.SetBinContent(j+1,yValue);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  // projMKPiHisto.Scale(1.0/projMKPiHisto.Integral());
  // projMKPiHisto.Scale(th->ProjectionX()->GetEntries());
  //projMKPiHisto.Scale(ratioMKPi);

  for (int j = 0; j < massKPi->numbins; ++j) {
    //pointsmassPsiPiXTot[j] = projmassPsiPiHisto.GetBinCenter(j+1);
    pointsMKPiYTot[j] = projMKPiHisto.GetBinContent(j+1);
    //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  }

  // for (int j = 0; j < massPsiPi->numbins; ++j) {
  //   pointsmassPsiPiXTot[j] = projmassPsiPiHisto.GetBinCenter(j+1);
  //   pointsmassPsiPiYTot[j] = projmassPsiPiHisto.GetBinContent(j+1);
  //   //std::cout <<" Bin " <<j<<" center = " <<projMKPiHisto.GetBinCenter(j+1)<<" : " <<mkpTotalProjection[j]<<std::endl;
  // }

  TGraph interPlotMKPi(massKPi->numbins,pointsMKPiXTot,pointsMKPiYTot);
  // TGraph interPlotmassPsiPi(massPsiPi->numbins,pointsmassPsiPiXTot,pointsmassPsiPiYTot);

  interPlotMKPi.SetLineColor(kRed); interPlotMKPi.SetLineWidth(2);
  // interPlotmassPsiPi.SetLineColor(kRed); interPlotmassPsiPi.SetLineWidth(2);


  th->ProjectionX()->Draw();
  interPlotMKPi.Draw("sameL");
  canvas->SaveAs("x.png");

  canvas->Clear();
  interPlotMKPi.Draw("AL");
  canvas->SaveAs("x1.png");

  canvas->Clear();
  projMKPiHisto.Draw();
  canvas->SaveAs("x2.png");

  // th->ProjectionY()->Draw();
  // interPlotmassPsiPi.Draw("sameL");
  // canvas->SaveAs("y.png");

  return 0;

}
