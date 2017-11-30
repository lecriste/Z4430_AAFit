#ifndef Z4430_CONSTANTS_H
#define Z4430_CONSTANTS_H

//const Double_t TMATH_PI = TMath::Pi(); // error: identifier "TMATH_PI" is undefined in device code (because GPU does not know about TMath)
const Double_t TMATH_PI = 3.14159265358979323846;

const Double_t MLb = 5.61951;
const Double_t MBd = 5.27961;
const Double_t MPsi2S = 3.686109;
const Double_t MJpsi = 3.096916;
const Double_t MProton = 0.938272046;
const Double_t MKaon = 0.493677;
const Double_t MPion = 0.13957018;

const Double_t M892 = 0.89581 ; const Double_t G892 = 0.0474; // From PDG neutral only K*(892)
//const Double_t M892 = 0.8961 ; const Double_t G892 = 0.0507; // From EvtGen
const Double_t M800 = 0.682; const Double_t G800 = 0.547; // From PDG
// K*800 Belle values: M = 0.946, G = 0.736 ?
//const Double_t M800 = 0.931; const Double_t G800 = 0.578; // From Belle
const Double_t M1410 = 1.414; const Double_t G1410 = 0.232;
const Double_t M1430_0 = 1.425; const Double_t G1430_0 = 0.270;
const Double_t M1430_2 = 1.4324; const Double_t G1430_2 = 0.109; // From PDG neutral only
const Double_t M1780_3 = 1.776; const Double_t G1780_3 = 0.159; // From PDG neutral only
const Double_t M2380_5 = 2.382; const Double_t G2380_5 = 0.178; // From PDG

const Double_t MLb2 = MLb*MLb;
const Double_t MLb4 = MLb2*MLb2;
const Double_t MBd2 = MBd*MBd;
const Double_t MBd4 = MBd2*MBd2;
const Double_t MPsi2S2 = MPsi2S*MPsi2S;
const Double_t MPsi2S4 = MPsi2S2*MPsi2S2;
const Double_t MJpsi2 = MJpsi*MJpsi;
const Double_t MJpsi4 = MJpsi2*MJpsi2;
const Double_t MProton2 = MProton*MProton;
const Double_t MProton4 = MProton2*MProton2;
const Double_t MKaon2 = MKaon*MKaon;
const Double_t MKaon4 = MKaon2*MKaon2;
const Double_t MPion2 = MPion*MPion;
const Double_t MPion4 = MPion2*MPion2;

#endif
