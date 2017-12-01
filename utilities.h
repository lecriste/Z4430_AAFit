#include "TMath.h"
#include <TLorentzVector.h>
//#include "constants.h"
#include "stdio.h"
#include <iostream>
using namespace std;

Double_t denom2_for_cosTheta_FromMasses(const Double_t sameSideM2, const Double_t psi_nSM2, const Double_t motherM2, const Double_t refM2, const Double_t otherM2) {
  return ((1./4.)*pow(motherM2 - psi_nSM2 + sameSideM2,2) - sameSideM2*motherM2) * ((1./4.)*pow(sameSideM2 - otherM2 + refM2,2) - sameSideM2*refM2) ;
}

// The same function is defined in /lustrehome/cristella/work/MyGooFit/PDFs/MatrixPdf.hh, since "calling a __host__ function from a __device__ function is not allowed"; remember to propagate any change:
Double_t cosTheta_FromMasses_host(const Double_t sameSideM2, const Double_t oppositeSideM2, const Double_t psi_nSM2, const Double_t motherM2, const Double_t refM2, const Double_t otherM2) {

  Double_t num = (sameSideM2/2)*(motherM2 + refM2 - oppositeSideM2) - (1./4.)*(motherM2 - psi_nSM2 + sameSideM2)*(sameSideM2 - otherM2 + refM2) ;
  Double_t denom2 = denom2_for_cosTheta_FromMasses(sameSideM2, psi_nSM2, motherM2, refM2, otherM2) ;

  return num / TMath::Sqrt(denom2) ;
}

//================ costheta_helicity ===========================
// This is a more general form of the above formula valid for any form of decay as follows :
// Mom decays into Dau1 and Dau2, Dau1 decays into GDau1 and GDau2
// The angle returned is the angle between GDau1 in Dau1 rest frame and initial direction of Dau1

Double_t costhetaHel_host(const Double_t m2Mom, const Double_t m2Dau1, const Double_t m2GDau1, const Double_t m2GDau2, const Double_t m2Dau2, const Double_t m2Dau2GDau2) {

    Double_t num      = 0.5*m2Dau1*(m2Mom-m2Dau2GDau2+m2GDau1)-0.25*(m2Mom-m2Dau2+m2Dau1)*(m2GDau1-m2GDau2+m2Dau1);
    Double_t denom2   = (0.25*(pow((m2Mom+m2Dau1-m2Dau2),2))-(m2Mom*m2Dau1))*(0.25*(pow((m2GDau1+m2Dau1-m2GDau2),2))-(m2GDau1*m2Dau1));
    Double_t denom    = TMath::Sqrt(denom2);
    
    return num/denom;
}
//================ costheta_helicity ===========================



//================ Decay Momentum ====================
// Momentum in 2-particle decay : m0->m1+m2

Double_t sq_calc_host(const Double_t x, const Double_t y, const Double_t z)
{
    return pow(x,2)+pow(y,2)+pow(z,2)-2.0*x*y-2.0*x*z-2.0*y*z; // corr sign
}
Double_t dec2mm_host (const Double_t m0, const Double_t m1, const Double_t m2)
{
    Double_t rootterm = sq_calc_host(m0*m0,m1*m1,m2*m2);
    if (rootterm >= 0) {
        return TMath::Sqrt(rootterm)/(2.0*m0);
    }
    else {
        cout <<"WARNING! In \"dec2mm\" function: rootterm (" <<rootterm <<") < 0 for m0 = " << m0 <<" going to m1 = " << m1 << " + m2 = "<< m2 <<" -> returning 0" <<endl;
        return 0.0;
    }
}
//================ Decay Momentum ====================




//================ Alpha =============================
Double_t alpha(const Double_t theta, const Double_t phi, const Double_t m2kpi, const Double_t m2jpsipi, const Double_t mPsi_nS )
{
    Double_t m2Psi = mPsi_nS*mPsi_nS;
    
    Double_t kmom = dec2mm_host(sqrt(m2kpi),MKaon,MPion);
    Double_t costh_k = costhetaHel_host(MBd2,m2kpi,MKaon2,MPion2,m2Psi,m2jpsipi);
    
    TLorentzVector K;
    Double_t pkx = kmom*sin(acos(costh_k));
    Double_t pky = 0.0;
    Double_t pkz = kmom*costh_k;
    Double_t Ek = sqrt(MKaon2+kmom*kmom);
    K.SetPxPyPzE(pkx,pky,pkz,Ek);
    TLorentzVector Pi;
    Double_t ppix = -kmom*sin(acos(costh_k));
    Double_t ppiy = 0.0;
    Double_t ppiz = -kmom*costh_k;
    Double_t Epi = sqrt(MPion2+kmom*kmom);
    Pi.SetPxPyPzE(ppix,ppiy,ppiz,Epi);
    
    // Jpsi mom = K* mom in B0 rest frame
    Double_t jpsimom = dec2mm_host(MBd,mPsi_nS,sqrt(m2kpi));
    TLorentzVector J_b0;
    J_b0.SetPxPyPzE(0.0,0.0,-jpsimom,sqrt(m2Psi+jpsimom*jpsimom));
    TLorentzVector Kstar_b0;
    Kstar_b0.SetPxPyPzE(0.0,0.0,jpsimom,sqrt(m2kpi+jpsimom*jpsimom));
    
    // boosting K* to Jpsi rest frame
    TLorentzVector Kstar_jpsi = Kstar_b0;
    Kstar_jpsi.Boost( -J_b0.BoostVector() );
    
    // boosting Jpsi to K* rest frame
    TLorentzVector J_Kstar = J_b0;
    J_Kstar.Boost( -Kstar_b0.BoostVector() );
    
    // boosting Pi in K* rest frame to Jpsi rest frame
    TLorentzVector Pi_jpsi = Pi;
    Pi_jpsi.Boost( -J_Kstar.BoostVector() );
    
    // Muon 4 momenta in Jpsi rest frame
    Double_t mumom = dec2mm_host(mPsi_nS,MMuon,MMuon);
    TLorentzVector muP;
    Double_t pmuPx = mumom*sin(theta)*cos(phi);
    Double_t pmuPy = -mumom*sin(theta)*sin(phi);
    Double_t pmuPz = -mumom*cos(theta);
    Double_t EmuP = sqrt(MMuon*MMuon+mumom*mumom);
    muP.SetPxPyPzE(pmuPx,pmuPy,pmuPz,EmuP);
    
    TVector3	MuPPiPlane	=	muP.Vect().Cross(Pi_jpsi.Vect()); //muP.Vect().Cross(Pi_jpsi.Vect());
    TVector3	MuPKstPlane		=	muP.Vect().Cross(Kstar_jpsi.Vect()); //muP.Vect().Cross(Kstar_jpsi.Vect());
    Double_t alph;
    if	(	MuPPiPlane.Cross(MuPKstPlane).Dot(-J_b0.Vect())	>	0.0	)
        alph	=	MuPPiPlane.Angle(MuPKstPlane);
    else
        alph	=	-MuPPiPlane.Angle(MuPKstPlane);
    
    return alph;
    
    
}
//================ Alpha =============================

//================ Theta tilde =======================
Double_t costhetatilde(const Double_t theta, const Double_t phi, const Double_t m2kpi, const Double_t m2jpsipi, const Double_t mPsi_nS)
{
    Double_t m2Psi = mPsi_nS*mPsi_nS;
    
    // K momentum in B0 frame
    Double_t pk_B0 = dec2mm_host(MBd,sqrt(m2jpsipi),MKaon);
    TLorentzVector K_B0;
    K_B0.SetPxPyPzE(0.0,0.0,pk_B0,sqrt(MKaon2+pk_B0*pk_B0));
    
    TLorentzVector Zc_B0;
    Zc_B0.SetPxPyPzE(0.0,0.0,-pk_B0,sqrt(m2jpsipi+pk_B0*pk_B0));
    
    // Boost K to Zc rest frame
    TLorentzVector K_Zc_old = K_B0;
    K_Zc_old.Boost( -Zc_B0.BoostVector() );
    
    // 4 momenta in Zc rest frame (with a different coordinate system)
    Double_t thetaz = acos(  costhetaHel_host(MBd2,m2jpsipi,m2Psi,MPion2,MKaon2,m2kpi)  );
    Double_t pk = K_Zc_old.Pz();
    Double_t Ek = sqrt(MKaon2+pk*pk);
    TLorentzVector K_Zc;
    K_Zc.SetPxPyPzE(pk*sin(thetaz),0.0,pk*cos(thetaz),Ek);
    
    Double_t ppi = dec2mm_host(sqrt(m2jpsipi),mPsi_nS,MPion);
    
    Double_t Epi = sqrt(MPion2+ppi*ppi);
    TLorentzVector Pi_Zc;
    Pi_Zc.SetPxPyPzE(0.0,0.0,ppi,Epi);
    
    Double_t Epsi = sqrt(m2Psi+ppi*ppi);
    TLorentzVector Jpsi_Zc;
    Jpsi_Zc.SetPxPyPzE(0.0,0.0,-ppi,Epsi);
    
    // Boost everything to Jpsi frame
    TLorentzVector K_jpsi = K_Zc;
    K_jpsi.Boost( -Jpsi_Zc.BoostVector() );
    
    TLorentzVector Pi_jpsi = Pi_Zc;
    Pi_jpsi.Boost( -Jpsi_Zc.BoostVector() );
    
    // Muon momenta in Jpsi rest frame
    Double_t pmu = dec2mm_host(mPsi_nS,MMuon,MMuon);
    Double_t Emu = sqrt(MMuon*MMuon + pmu*pmu);
    
    Double_t denom = sqrt( (0.25*pow((MBd2-m2kpi+m2Psi),2)-MBd2*m2Psi)*(0.25*m2Psi*m2Psi-MMuon*MMuon*m2Psi) );
    Double_t m2kpimu = 0.5*( MBd2+m2kpi+2.0*MMuon*MMuon-m2Psi-4.0*cos(theta)*denom/m2Psi );
    
    Double_t Ek_jpsi = K_jpsi.E();
    Double_t pkx = K_jpsi.Px();
    Double_t pkz = K_jpsi.Pz();
    
    Double_t Epi_jpsi = Pi_jpsi.E();
    Double_t ppiz = Pi_jpsi.Pz();
    
    Double_t a = pow((Ek_jpsi+Emu+Epi_jpsi),2) - m2kpimu - (pkx*pkx + pkz*pkz +2.0*pkz*ppiz + pmu*pmu + ppiz*ppiz );
    Double_t b = 2.0*pmu*(pkz+ppiz);
    Double_t c = 2.0*pkx*pmu*cos(phi);
    Double_t discr = c*c*(b*b+c*c-a*a);
    
    
    TLorentzVector PKPi = K_jpsi + Pi_jpsi;
    Double_t kpi_angle = PKPi.Vect().Angle(Pi_jpsi.Vect());
    
    
    if ( discr>=0.0 && fabs( TMath::Cos(theta) )<0.98 ) {
        
        Double_t sinth1 = -(a*c*c-b*sqrt(discr))/(b*b*c+c*c*c);
        Double_t costh1 = ( a*b + sqrt(discr) )/(b*b+c*c);
        Double_t th1 = TMath::ACos(costh1);
        
        Double_t sinth2 = -(a*c*c+b*sqrt(discr))/(b*b*c+c*c*c);
        Double_t costh2 = ( a*b - sqrt(discr) )/(b*b+c*c);
        Double_t th2 = TMath::ACos(costh2);
        
        Double_t costh_big, costh_small;
        
        if (th1>th2) {
            costh_big = costh1;
            costh_small = costh2;
        }
        else {
            costh_big = costh2;
            costh_small = costh1;
        }
        
        /*
        if ( fabs(phi)>=1.570796 ) {
            return costh_small;
        }
        else {
            return costh_big;
        }
        */
        
        if (sinth1 <=0.0) {return costh1; }
        else {return costh2;}
        
    }
    else {
        if (TMath::Cos(theta) >= 0.0 ){
            return TMath::Cos(kpi_angle);
        }
        else {
            return -TMath::Cos(kpi_angle);
        }
    }
    
}
//================ Theta tilde =======================


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

Bool_t angles_contour_host(const Double_t cos, const Double_t phi) {

  if ( (fabs(cos) > 1) || (fabs(phi) > TMath::Pi()) )
    return kFALSE;
  else  
    return kTRUE ;

}

