#include "ThreeBodiesPsiPiKPdf.hh"
#include "TRandom.h"

EXEC_TARGET fptype device_ThreeBodiesPsiPiK (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[indices[2 + indices[0]]];

  fptype mP = p[indices[1]];
  fptype m1 = p[indices[2]];
  fptype m2 = p[indices[3]];
  fptype m3 = p[indices[4]];

  fptype ret = isnan(sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x)) ? 0 : (sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x));


return ret;
}

EXEC_TARGET fptype device_ThreeBodiesPsiPiK_Point (fptype* point, fptype* p, unsigned int* indices) {
  fptype x = point[0];

  fptype mP = p[indices[1]];
  fptype m1 = p[indices[2]];
  fptype m2 = p[indices[3]];
  fptype m3 = p[indices[4]];

  fptype ret = isnan(sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x)) ? 0 : (sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x));

  return ret;
}

/*
EXEC_TARGET fptype device_Three_Bin (fptype* evt, fptype* p, unsigned int* indices) {
  fptype xmid = evt[indices[2 + indices[0]]];
  fptype xnex = evt[indices[2 + indices[0]]+3];

  fptype bin = xnex-xmid;
  fptype step = bin; //Integration steps
  int n = 4;
  step /= n;

  //STEPS FOR INTEGRATION
  fptype Xs[5];

  for (int k=0;k<=n;k++){

    Xs[k]=(-2+k)*step+xmid;

  }

  //NEWTON-COTES COEFFICIENTS
  fptype cnw = 2.0;
  cnw /=45.0;
  fptype b[5] = {7.0,32.0,12.0,32.0,7.0};

  fptype a[5];

  for (int k=0;k<=n;k++){
        a[k] = b[k]*cnw;
  }

  //FUNCTION IN Xi
  fptype fX[5];
  for (int k=0;k<=n;k++){
    fptype x;
    x=Xs[k];

    fX[k] = (x<1.015)||(x>1.7)? 0 : isnan((sqrt( pow(x+3.0967,4) + pow(3.0967,4) + pow(1.01946,4) - 2*pow(x+3.0967,2)*pow(3.0967,2) - 2*pow(3.0967,2)*pow(1.01946,2) - 2*pow(x+3.0967,2)*pow(1.01946,2) ) * sqrt( pow(5.279,4) + pow(x+3.0967,4) + pow(0.493677,4) - 2*pow(5.279,2)*pow(x+3.0967,2) - 2*pow(5.279,2)*pow(0.493677,2) - 2*pow(x+3.0967,2)*pow(0.493677,2) ) / (x+3.0967)))? 0 : (sqrt( pow(x+3.0967,4) + pow(3.0967,4) + pow(1.01946,4) - 2*pow(x+3.0967,2)*pow(3.0967,2) - 2*pow(3.0967,2)*pow(1.01946,2) - 2*pow(x+3.0967,2)*pow(1.01946,2) ) * sqrt( pow(5.279,4) + pow(x+3.0967,4) + pow(0.493677,4) - 2*pow(5.279,2)*pow(x+3.0967,2) - 2*pow(5.279,2)*pow(0.493677,2) - 2*pow(x+3.0967,2)*pow(0.493677,2) ) / (x+3.0967));

  }

  //INTEGRAL
  fptype integralF = 0;

   for (int k=0;k<=n;k++){
        integralF += a[k]*fX[k];
   }
    integralF *= step;

return integralF;

}
*/

MEM_DEVICE device_function_ptr ptr_to_ThreeBodiesPsiPiK = device_Three;
//MEM_DEVICE device_function_ptr ptr_to_ThreeBodiesPsiPiK_Bin = device_Three_Bin;
MEM_DEVICE device_function_ptr ptr_to_ThreeBodiesPsiPiK_Point = device_ThreeBodiesPsiPiK_Point;


__host__ ThreeBodiesPsiPiK::ThreeBodiesPsiPiK (std::string n, Variable* _x,Variable* _mp,Variable* _m1,Variable* _m2,Variable* _m3)
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;

  GET_FUNCTION_ADDR(ptr_to_ThreeBodiesPsiPiK);
  //GET_INTEGRAL_ADDR(ptr_to_Three_Bin);
  GET_ATPOINTS_ADDR(ptr_to_ThreeBodiesPsiPiK_Point);
  initialise(pindices);
}

__host__ fptype ThreeBodiesPsiPiK::integrate (fptype lo, fptype hi) const {

  int loind=0, hiind=0;
  fptype max =2.0;
  fptype min =1.0;

  if(lo<=1.0){
  if(hi>=2) return 5.06186;
  else {
  while(hi<max){
  max=PUNTI[255-hiind];
  hiind++;
  }
  return AREA[255-hiind];
  }}

  if(hi>=2.0){
  if(lo<=1) return 5.06186;
  else {
  while(lo>min){
  min=PUNTI[loind];
  loind++;
  }
  return AREA[loind];
  }}

  while(lo>min){
  min=PUNTI[loind];
  loind++;
  }

  while(hi<max){
  max=PUNTI[255-hiind];
  hiind++;
  }

  return AREA[255-hiind]-AREA[loind];
}
/*
__global__ void fillHisto (TH1F* histo,unsigned int nevents,fptype fmax){


    fptype roll,func,x;
    fptype xmin = histo->GetBinLowEdge(1);
    fptype xmax = histo->GetBinLowEdge(histo->GetNbinsX()+1);
    long int ms; struct timeval tp;

    for (int j = 0; j < nevents; ++j) {

		    gettimeofday(&tp,NULL);
        ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
        TRandom3 ranGen(ms);

  			x = donram.Uniform(xmax-xmins)+xmin;
  			func = background(xvar->value);
  			roll = donram.Uniform(fmax);

  			if (roll > func) {
  				--j;
  				continue; }

  			if ((x < xmin) || (x > xmax)) {
  				--j;
  				continue;}

  			histo->Fill(x);

  		}


}
*/
