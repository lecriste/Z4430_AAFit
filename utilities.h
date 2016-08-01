Double_t cosTheta_FromMasses(const Double_t sameSideM2, const Double_t oppositeSideM2, const Double_t psi_nSM2, const Double_t motherM2, const Double_t refM2, const Double_t otherM2) {

  Double_t num = (sameSideM2/2)*(motherM2 + refM2 - oppositeSideM2) - (1./4.)*(motherM2 - psi_nSM2 + sameSideM2)*(sameSideM2 - otherM2 + refM2) ;
  Double_t denom2 = ((1./4.)*pow(motherM2 - psi_nSM2 + sameSideM2,2) - sameSideM2*motherM2) * ((1./4.)*pow(sameSideM2 - otherM2 + refM2,2) - sameSideM2*refM2) ;

  return num / TMath::Sqrt(denom2) ;
}
