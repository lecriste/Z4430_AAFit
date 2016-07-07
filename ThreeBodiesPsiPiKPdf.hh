#ifndef THREE_PDF_PSIPIK_HH
#define THREE_PDF_PSIPIK_HH

#include "GooPdf.hh"
#include "TH1F.h"
#include "DataSet.hh"
#include "BinnedDataSet.hh"
#include "UnbinnedDataSet.hh"

class ThreeBodiesPsiPiK : public GooPdf {
public:

    ThreeBodiesPsiPiK (std::string n, Variable* _x,Variable* _mp,Variable* _m1,Variable* _m2,Variable* _m3);//,Variable* nBk);
  __host__ fptype integrate (fptype lo, fptype hi) const;
  __host__ virtual bool hasAnalyticIntegral () const {return false;}
  __global__ void fillHisto (TH1F* histo,unsigned int nevents);
  __global__ void fillDataset (BinnedDataSet dataset,unsigned int nevents);
  __global__ void fillDataset (UnbinnedDataSet dataset,unsigned int nevents);


private:

};

#endif
