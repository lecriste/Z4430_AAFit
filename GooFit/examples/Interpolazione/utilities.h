Double_t denom2_for_cosTheta_FromMasses(const Double_t sameSideM2, const Double_t psi_nSM2, const Double_t motherM2, const Double_t refM2, const Double_t otherM2) {
  return ((1./4.)*pow(motherM2 - psi_nSM2 + sameSideM2,2) - sameSideM2*motherM2) * ((1./4.)*pow(sameSideM2 - otherM2 + refM2,2) - sameSideM2*refM2) ;
}

// The same function is defined in /lustrehome/cristella/work/MyGooFit/PDFs/MatrixPdf.hh, since "calling a __host__ function from a __device__ function is not allowed"; remember to propagate any change:
Double_t cosTheta_FromMasses_host(const Double_t sameSideM2, const Double_t oppositeSideM2, const Double_t psi_nSM2, const Double_t motherM2, const Double_t refM2, const Double_t otherM2) {

  Double_t num = (sameSideM2/2)*(motherM2 + refM2 - oppositeSideM2) - (1./4.)*(motherM2 - psi_nSM2 + sameSideM2)*(sameSideM2 - otherM2 + refM2) ;
  Double_t denom2 = denom2_for_cosTheta_FromMasses(sameSideM2, psi_nSM2, motherM2, refM2, otherM2) ;

  return num / TMath::Sqrt(denom2) ;
}
/*
TString cosTheta_FromMasses_TFormula_host(const TString sameSideM2, const TString oppositeSideM2, const Double_t psi_nSM2, const Double_t motherM2, const Double_t refM2, const Double_t otherM2) {

  RooConstVar m2Psi_nS("m2Psi_nS","m2Psi_nS",psi_nSM2); const char* m2Psi_nSName = m2Psi_nS.GetName();
  RooConstVar m2Mother("m2Mother","m2Mother",motherM2); const char* m2MotherName = m2Mother.GetName();
  RooConstVar m2Ref("m2Ref","m2Ref",refM2);             const char* m2RefName = m2Ref.GetName();
  RooConstVar m2Other("m2Other","m2Other",otherM2);     const char* m2OtherName = m2Other.GetName();

  const char* m2SameSide = sameSideM2.Data();
  const char* m2OppositeSide = oppositeSideM2.Data();

  TString num = TString::Format("(%s/2)*(%s + %s - %s) - (1./4.)*(%s - %s + %s)*(%s - %s + %s)", m2SameSide,m2MotherName,m2RefName,m2OppositeSide, m2MotherName,m2Psi_nSName,m2SameSide,m2SameSide,m2OtherName,m2RefName);
  TString denom = TString::Format();
  Double_t denom2 = ((1./4.)*pow(motherM2 - psi_nSM2 + sameSideM2,2) - sameSideM2*motherM2) * ((1./4.)*pow(sameSideM2 - otherM2 + refM2,2) - sameSideM2*refM2) ;

  return num / TMath::Sqrt(denom2) ;
}
*/

Bool_t Dalitz_contour_host(const Double_t mKP, const Double_t mPsiP, const Bool_t massSquared, const Int_t psi_nS) {

  Double_t MPsi_nS = 0; 
  if (psi_nS == 1) 
    MPsi_nS = MJpsi;
  else if (psi_nS == 2)
    MPsi_nS = MPsi2S;
  else {
    cout <<"psi_nS = " <<psi_nS <<" not allowed at the moment." <<endl;
    return kFALSE;
  }

  Double_t mKP_1 = mKP;
  Double_t mPsiP_1 = mPsiP;
  if (massSquared) {
    mKP_1 = TMath::Sqrt( mKP );
    mPsiP_1 = TMath::Sqrt( mPsiP );
  }

  if ((mKP_1 < MKaon + MPion) || (mKP_1 > MBd - MPsi_nS) || (mPsiP_1 < MPsi_nS + MPion) || (mPsiP_1 > MBd - MKaon))
    return kFALSE;
  else { // Dalitz border from PDG KINEMATICS 43.4.3.1. 
    Float_t E_P = (mPsiP_1*mPsiP_1 - MJpsi2 + MPion2)/(2*mPsiP_1) ;
    Float_t E_K = (MBd2 - mPsiP_1*mPsiP_1 - MKaon2)/(2*mPsiP_1) ;
    Float_t E_PpE_K_2 = TMath::Power((E_P + E_K),2);
    Float_t sqrt_E_P2mMP2 = TMath::Sqrt(E_P*E_P - MPion2);
    Float_t sqrt_E_K2mMK2 = TMath::Sqrt(E_K*E_K - MKaon2);
    Float_t mKP2_min = E_PpE_K_2 - TMath::Power(sqrt_E_P2mMP2 + sqrt_E_K2mMK2,2);
    Float_t mKP2_max = E_PpE_K_2 - TMath::Power(sqrt_E_P2mMP2 - sqrt_E_K2mMK2,2);
    if ((mKP_1*mKP_1 < mKP2_min) || (mKP_1*mKP_1 > mKP2_max))
      return kFALSE;
  }

  return kTRUE ;

}
