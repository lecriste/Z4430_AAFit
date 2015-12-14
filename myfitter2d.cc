#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
//
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooDecay.h>
#include <RooGaussModel.h>
#include <RooAddModel.h>
#include <RooPlot.h>
//
#include "myloop.h"
//#include "plotDressing2D.h"

using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1
#define DO_MINOS            kTRUE
// 0 - w/o DISPLAY
// 1 - w/  DISPLAY
#define DISPLAY             1

#define mumuPi_MIN          3.5
#define mumuPi_MAX          6.0
#define piK_MIN             0.5
#define piK_MAX             6.0
#define MASS_PEAK           BP_MASS
#define muon_massPDG        0.1056583715
#define pionCh_massPDG      0.13957018
#define kaonCh_massPDG      0.493677
#define beam_energy         4000
#define SOURCE              "BdToJpsiPPbar_18Mar_MuMuPPbarPAT_ntpl.root"

// angles variables
double getPlanesAngle(TLorentzVector B0, TLorentzVector K, TLorentzVector Pi, TLorentzVector muM, TLorentzVector muP);
float GetThetaMuMu(TLorentzVector BVec, TLorentzVector JPsiVec, TLorentzVector MuPlusVec, float BeamEnergy, float JPsiPDG , float muonPDG);
void GetMomentumInMotherFrame(TLorentzVector Mother, TLorentzVector Particle, double BeamEnergy , TVector3 &Particle_rotated);


void myfitter2d()
{
    // define variables needed in the fit
    RooRealVar mumuPi_mass("mumuPi_mass","mumuPi mass",mumuPi_MIN,mumuPi_MAX); 
    RooRealVar piK_mass("piK_mass","piK mass",piK_MIN,piK_MAX); 
    //RooRealVar cterr("cterr","cterr",0.0001,0.008);
    
    //outout
    TFile *fout = new TFile("myfitter2d.root","recreate");     // output file    
    //TNtupleD *_nt = new TNtupleD("_nt","_nt","mumuPi_mass:piK_mass:cterr"); // output ntuple
    TNtupleD *_nt = new TNtupleD("_nt","_nt","mumuPi_mass:piK_mass"); // output ntuple
    
    // input 
    TFile *fin = new TFile(SOURCE); 
    TTree *tin = (TTree*)fin->Get("mkcands/Z_data");
    if (!tin) cout <<"\nNo TTree \"mkcands/Z_data\" found in file \"" <<fin->GetName() <<"\"" <<endl;
    else cout <<"\nTTree \"mkcands/Z_data\" found in file \"" <<fin->GetName() <<"\"" <<endl;
   
    // setting up rootuple for reading
    MCmupPx = 0; MCmupPy = 0; MCmupPz = 0;
    MCmumPx = 0; MCmumPy = 0; MCmumPz = 0;
    MCpionPx = 0; MCpionPy = 0; MCpionPz = 0;
    MCkaonPx = 0; MCkaonPy = 0; MCkaonPz = 0;
    
    ReducedBranches br;
    cout <<"\nbefore setbranchadd" <<endl;
    br.setbranchadd(tin);
    
    // reading rootuple
    //for (int evt=0; evt < tin->GetEntries(); evt++) 
    for (int evt=0; evt < 10000; evt++) 
       {
	 if (evt%1000 == 0)
	   cout <<"\nGetting entry number " <<evt <<endl ;
	 tin->GetEntry(evt);
	 // 
	 // cuts to select events/cands
	 /*
	   if (br.hltbook[HLT_Dimuon16_Jpsi_v1]!=1) continue;
	   if (br.vtxprob<=0.15) continue;
	   if (br.tk1pt<=2.0) continue;
	 */
	 if (br.nMCB0 <= 0) continue;
	 //
	 //cout <<"\nbefore accessing branches" <<endl;
	 TLorentzVector muP; muP.SetXYZM((*br.MCmupPx)[0], (*br.MCmupPy)[0], (*br.MCmupPz)[0], muon_massPDG);
	 TLorentzVector muM; muM.SetXYZM((*br.MCmumPx)[0], (*br.MCmumPy)[0], (*br.MCmumPz)[0], muon_massPDG);
	 TLorentzVector pion; pion.SetXYZM((*br.MCpionPx)[0], (*br.MCpionPy)[0], (*br.MCpionPz)[0], pionCh_massPDG);
	 TLorentzVector kaon; kaon.SetXYZM((*br.MCkaonPx)[0], (*br.MCkaonPy)[0], (*br.MCkaonPz)[0], kaonCh_massPDG);
	 //cout <<"\nafter accessing branches" <<endl;
	 //
	 TLorentzVector mumu = muP + muM;
	 TLorentzVector mumuPi = mumu + pion;
	 TLorentzVector piK = pion + kaon;
	 TLorentzVector mumuPiK = mumu + piK;
	 //cout <<"\nbefore angles" <<endl;
	 Double_t theta_mumu = GetThetaMuMu(mumuPiK, mumu, muP, beam_energy, mumu.M(), muP.M()) ;
	 Double_t phi_gen = getPlanesAngle(mumuPiK, kaon, pion, muM, muP) ;
	 //cout <<"\nafter angles" <<endl;
	 //
	 // filling the variables in the output ntuple
	 double var[4];
	 var[0] = mumuPi.M();
	 var[1] = piK.M();
	 var[2] = theta_mumu;
	 var[3] = phi_gen;
	
	 //cout <<"\nBefore filling output ntuple" <<endl; 
	 _nt->Fill(var);
	 //cout <<"\nAfter filling output ntuple" <<endl; 
       }
    //
    fin->Close();
    //
    // the dataset contains only the 3 variables of interest
    //RooDataSet *data = new RooDataSet("data","data",_nt,RooArgSet(mumuPi_mass,piK_mass,cterr));
    RooDataSet *data = new RooDataSet("data","data",_nt,RooArgSet(mumuPi_mass,piK_mass));
    //
    /////////////////////////////////////////////////////
    //
    // initialization
    //
    /*
    double n_signal_initial = data->sumEntries(TString::Format("abs(mass-%g)<0.015",MASS_PEAK))
    - data->sumEntries(TString::Format("abs(mass-%g)<0.030&&abs(mass-%g)>0.015",MASS_PEAK,MASS_PEAK));
    //
    double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
    //
    //-----------------------------------------------------------------
    //
    // signal PDF 
    //===========
    //
    // double gaussian for the signal in mass
    //
    RooRealVar m_mean("m_mean","m_mean",MASS_PEAK,mumuPi_MIN,mumuPi_MAX);
    RooRealVar m_sigma1("m_sigma1","m_sigma1",0.016,0.001,0.045);
    RooRealVar m_sigma2("m_sigma2","m_sigma2",0.035,0.001,0.090);
    RooRealVar m_fraction("m_fraction","m_fraction",0.5);
    //
    RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",mass,m_mean,m_sigma1);
    RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",mass,m_mean,m_sigma2);
    //
    RooAddPdf pdf_m_signal("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
    //
    // exponential convoluted with gaussian resolution for the signal in ct
    //
    RooRealVar res_sig_mean("res_sig_mean","res_sig_mean",0.0,-1.,1.);
    RooRealVar res_sig_sigma("res_sig_sigma","res_sig_sigma",1.0,0.3,2.0);
    //
    RooGaussModel res_signal("res_signal","res_signal",ct,res_sig_mean,res_sig_sigma,cterr);
    //
    RooRealVar ctau("ctau","ctau",0.04911,0.010,0.090);
    RooDecay pdf_t_signal("pdf_t_signal","pdf_t_signal",ct,ctau,res_signal,RooDecay::SingleSided);
    //
    // bidimensional signal pdf
    //
    RooProdPdf pdf_signal("pdf_signal","pdf_signal", RooArgSet(pdf_m_signal, pdf_t_signal));
    //
    //
    // combinatorial background PDF (prompt or non-prompt J/psi + random track)
    //=========================================================================
    //
    // exponential for the combinatorial background in mass
    //
    RooRealVar m_par1("m_par1","m_par1",-0.3,-2.,+2.);
    RooExponential pdf_m_combinatorial("pdf_m_combinatorial","pdf_m_combinatorial",mass,m_par1);
    //
    // exponential convoluted with gaussian resolution for the non-prompt background in ct
    //
    RooRealVar ctau_nonprompt("ctau_nonprompt","ctau_nonprompt",0.0500, 0.0010, 0.1000);
    RooDecay pdf_t_nonprompt("pdf_t_nonprompt","pdf_t_nonprompt",ct,ctau_nonprompt,res_signal,RooDecay::SingleSided);
    //
    // Sum of gaussian resolution function (res_signal) for prompt background in ct and the previous exponential for NP-bkg
    //
    RooRealVar prompt_fraction("prompt_fraction","prompt_fraction",0.5,0.0,1.0); 
    //
    RooAddPdf pdf_t_combinatorial("pdf_t_combinatorial","pdf_t_combinatorial",RooArgList(res_signal,pdf_t_nonprompt),RooArgList(prompt_fraction));
    //
    // bidimensional combinatorial-bkg pdf
    //
    RooProdPdf pdf_combinatorial("pdf_combinatorial","pdf_combinatorial",RooArgSet(pdf_m_combinatorial,pdf_t_combinatorial));
    //
    //-----------------------------------------------------------------
    //
    // B->J/psi+track+X background PDF
    //================================
    // 
    // single gaussian for the physical background in mass
    //
    RooRealVar m_jpsix_mean("m_jpsix_mean","m_jpsix_mean",5.1,5.0,5.3);
    RooRealVar m_jpsix_sigma("m_jpsix_sigma","m_jpsix_sigma",0.05,0.01,0.10);
    RooGaussian pdf_m_jpsix("pdf_m_jpsix","pdf_m_jpsix",mass,m_jpsix_mean,m_jpsix_sigma);
    //
    // exponential convoluted with gaussian resolution for the physical background in ct
    //
    RooRealVar ctau_jpsix("ctau_jpsix","ctau_jpsix",0.0500, 0.0010, 0.1000);
    RooDecay pdf_t_jpsix("pdf_t_jpsix","pdf_t_jpsix",ct,ctau_jpsix,res_signal,RooDecay::SingleSided); 
    //
    // bidimensional physical-bkg pdf
    //
    RooProdPdf pdf_jpsix("pdf_jpsix","pdf_jpsix",RooArgSet(pdf_m_jpsix, pdf_t_jpsix));
    //
    //-----------------------------------------------------------------
    // 
    // FULL MODEL (SIGNAL + 2 BKGS)
    //
    // define coefficients for addition of the 3 bidimentional pdfs
    //
    RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
    RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
    RooRealVar n_jpsix("n_jpsix","n_jpsix",200.,0.,data->sumEntries());
    //
    RooAddPdf model("model","model",
                    RooArgList(pdf_signal, pdf_combinatorial, pdf_jpsix),
                    RooArgList(n_signal, n_combinatorial, n_jpsix));
    //
    /////////////////////////////////////////////////////////////////////
    //
    // finally go for fitting !
    //
    model.fitTo(*data,Minos(DO_MINOS),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));
    //
    // go to display plots with fits superimposed on data distributions
#if DISPLAY
    //
    // Display mass plots
    //-------------------
    //
    TCanvas *c1 = canvasDressing("c1");
    //
    RooPlot* frame_m = mass.frame();
    //
    TH1D* histo_data_m = (TH1D*)data->createHistogram("histo_data_m", mass, Binning(50,mass.getMin(),mass.getMax()));
    //
    histo_data_m->Sumw2(false);
    histo_data_m->SetBinErrorOption(TH1::kPoisson);
    histo_data_m->SetMarkerStyle(20);
    histo_data_m->SetMarkerSize(0.8);
    histo_data_m->SetLineColor(kBlack);
    for (int i=1; i<=50; i++)
        if (histo_data_m->GetBinContent(i)==0) histo_data_m->SetBinError(i,0.);
    //
    data->plotOn(frame_m,Binning(50),Invisible());
    model.plotOn(frame_m,Precision(5E-4));
    model.plotOn(frame_m,Precision(5E-4),Components(pdf_signal),LineColor(kRed),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(2), VLines(), DrawOption("F"));
    model.plotOn(frame_m,Precision(5E-4),Components(pdf_combinatorial),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_m,Precision(5E-4),Components(pdf_jpsix),LineColor(kBlack),LineWidth(2),LineStyle(7));
    //
    frame_m->SetTitle("");
    frame_m->GetXaxis()->SetTitle("M_{J/#psi K^{#pm}} [GeV]");
    frame_m->GetXaxis()->SetLabelFont(42);
    frame_m->GetXaxis()->SetLabelOffset(0.01);
    frame_m->GetXaxis()->SetTitleSize(0.06);
    frame_m->GetXaxis()->SetTitleOffset(1.09);
    frame_m->GetXaxis()->SetLabelFont(42);
    frame_m->GetXaxis()->SetLabelSize(0.055);
    frame_m->GetXaxis()->SetTitleFont(42);
    frame_m->GetYaxis()->SetTitle("Events / 20 MeV");
    frame_m->GetYaxis()->SetLabelFont(42);
    frame_m->GetYaxis()->SetLabelOffset(0.01);
    frame_m->GetYaxis()->SetTitleOffset(1.14);
    frame_m->GetYaxis()->SetTitleSize(0.06);
    frame_m->GetYaxis()->SetTitleFont(42);
    frame_m->GetYaxis()->SetLabelFont(42);
    frame_m->GetYaxis()->SetLabelSize(0.055);
    //
    frame_m->Draw();
    histo_data_m->Draw("Esame");
    LegendMassProj();
    //
    c1->SaveAs("massPlot_13TeV.png");
    //
    //c1->Delete();
    //
    /////////////////////////////////////////////////////
    //
    // Display c*proper-time plots 
    //----------------------------
    //
    TCanvas *c2 = canvasDressing("c2");
    //
    RooPlot* frame_t = ct.frame();
    //
    TH1D* histo_data_t = (TH1D*)data->createHistogram("histo_data_t", ct, Binning(120,ct.getMin(),ct.getMax()));
    //
    histo_data_t->Sumw2(false);
    histo_data_t->SetBinErrorOption(TH1::kPoisson);
    histo_data_t->SetMarkerStyle(20);
    histo_data_t->SetMarkerSize(0.8);
    histo_data_t->SetLineColor(kBlack);
    for (int i=1; i<=120; i++)
        if (histo_data_t->GetBinContent(i)==0) histo_data_t->SetBinError(i,0.);
    //
    data->plotOn(frame_t,Binning(120),Invisible());
    model.plotOn(frame_t,Precision(5E-4),ProjWData(cterr,*data));
    model.plotOn(frame_t,Precision(5E-4),ProjWData(cterr,*data),Components(pdf_signal),LineColor(kRed),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(2), VLines(), DrawOption("F"));
    model.plotOn(frame_t,Precision(5E-4),ProjWData(cterr,*data),Components(pdf_combinatorial),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_t,Precision(5E-4),ProjWData(cterr,*data),Components(pdf_jpsix),LineColor(kBlack),LineWidth(2),LineStyle(7));
    //
    frame_t->SetTitle("");
    frame_t->GetXaxis()->SetTitle("ct [cm]");
    frame_t->GetXaxis()->SetLabelFont(42);
    frame_t->GetXaxis()->SetLabelOffset(0.01);
    frame_t->GetXaxis()->SetTitleSize(0.06);
    frame_t->GetXaxis()->SetTitleOffset(1.09);
    frame_t->GetXaxis()->SetLabelFont(42);
    frame_t->GetXaxis()->SetLabelSize(0.055);
    frame_t->GetXaxis()->SetTitleFont(42);
    frame_t->GetYaxis()->SetTitle("Events / 25 #mum");
    frame_t->GetYaxis()->SetLabelFont(42);
    frame_t->GetYaxis()->SetLabelOffset(0.01);
    frame_t->GetYaxis()->SetTitleOffset(1.14);
    frame_t->GetYaxis()->SetTitleSize(0.06);
    frame_t->GetYaxis()->SetTitleFont(42);
    frame_t->GetYaxis()->SetLabelFont(42);
    frame_t->GetYaxis()->SetLabelSize(0.055);
    //
    frame_t->GetYaxis()->SetRangeUser(0.5,histo_data_t->GetMaximum()*2.);
    //
    frame_t->Draw();
    histo_data_t->Draw("Esame");
    LegendLifetimeProj();
    //
    c2->SetLogy();
    //
    c2->SaveAs("ctPlot_13TeV.png");
    //
    c1->Delete();
    c2->Delete();
    //
#endif
    //
    */
    fout->Write();
    fout->Close();
}

// angles variables
double getPlanesAngle(TLorentzVector B0, TLorentzVector K, TLorentzVector Pi, TLorentzVector muM, TLorentzVector muP) {

  TLorentzVector KBoosted = K, PiBoosted = Pi, muMBoosted = muM, muPBoosted = muP ;
  muPBoosted.Boost( -B0.BoostVector() ) ;
  muMBoosted.Boost( -B0.BoostVector() ) ;
  KBoosted.Boost( -B0.BoostVector() ) ;
  PiBoosted.Boost( -B0.BoostVector() ) ;
  TVector3 MuMuPlane_normal = muPBoosted.Vect().Cross( muMBoosted.Vect() ) ;
  TVector3 KPiPlane_normal = KBoosted.Vect().Cross( PiBoosted.Vect() ) ;
  Double_t angle = 10. ;
  /*
          TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
  */
  if (MuMuPlane_normal.Cross(KPiPlane_normal).Dot( -B0.Vect() ) > 0.0)
    angle = MuMuPlane_normal.Angle(KPiPlane_normal);
  else
    angle = -MuMuPlane_normal.Angle(KPiPlane_normal);

  return angle;
}

float GetThetaMuMu(TLorentzVector BVec, TLorentzVector JPsiVec, TLorentzVector MuPlusVec, float BeamEnergy, float JPsiPDG , float muonPDG) {

  TVector3 JPsiInBFrame, MuInBFrame, MuInJPsiFrame, MuInJPsiFromBFrame;
  TLorentzVector JPsiInBFrameTLVec, MuInBFrameTLVec;
  /*
  // Alessandra
  // get the momentum of the X in the in the B rest-frame : XInBFrame
  GetMomentumInMotherFrame( BVec , XVec , BeamEnergy, XInBFrame);
  XInBFrameTLVec.SetPtEtaPhiM(XInBFrame.Perp() , XInBFrame.Eta(),  XInBFrame.Phi() ,XCandPDG);
  // get the momentum of the J/psi int he X rest-frame with two steps in cascade:
  // 1) step 1 : apply B boost in lab
  GetMomentumInMotherFrame( BVec , JPsiVec, BeamEnergy, JPsiInBFrame);
  JPsiInBFrameTLVec.SetPtEtaPhiM( JPsiInBFrame.Perp() , JPsiInBFrame.Eta(),JPsiInBFrame.Phi() ,JPsiPDG );
  // 2) step 2: apply X boost in B rest-frame
  GetMomentumInMotherFrame( XInBFrameTLVec, JPsiInBFrameTLVec, BeamEnergy, JPsiInXFromBFrame);
  JPsiInXFromBFrameTLVec.SetPtEtaPhiM(JPsiInXFromBFrame.Perp() ,JPsiInXFromBFrame.Eta(),JPsiInXFromBFrame.Phi() ,JPsiPDG );
  */

  // B0 - >J/psi K pi
  
  // get the momentum of the J/psi in the in the B rest-frame : JPsiInBFrame
  GetMomentumInMotherFrame( BVec , JPsiVec , BeamEnergy, JPsiInBFrame);
  JPsiInBFrameTLVec.SetPtEtaPhiM(JPsiInBFrame.Perp() , JPsiInBFrame.Eta(),  JPsiInBFrame.Phi() , JPsiPDG);

  // get momentum of the mu+ in JPsi rest-frame
  GetMomentumInMotherFrame( BVec , MuPlusVec, BeamEnergy, MuInBFrame); // B boost
  MuInBFrameTLVec.SetPtEtaPhiM(MuInBFrame.Perp() , MuInBFrame.Eta(),  MuInBFrame.Phi() ,muonPDG);
  //GetMomentumInMotherFrame(XInBFrameTLVec, MuInBFrameTLVec, BeamEnergy, MuInXFromBFrame); // X in B r.f. boost
  //MuInXFromBFrameTLVec.SetPtEtaPhiM(MuInXFromBFrame.Perp() ,MuInXFromBFrame.Eta(), MuInXFromBFrame.Phi() ,muonPDG);
  //GetMomentumInMotherFrame( JPsiInXFromBFrameTLVec, MuInXFromBFrameTLVec, BeamEnergy, MuInJPsiFromXFromBFrame);
  GetMomentumInMotherFrame( JPsiInBFrameTLVec, MuInBFrameTLVec, BeamEnergy, MuInJPsiFromBFrame);

  float thetaJPsi = MuInJPsiFromBFrame.Angle(JPsiInBFrame);
  
  /*
  // lab frame
  GetMomentumInMotherFrame( JPsiVec , MuPlusVec, BeamEnergy, MuInJPsiFrame); // B boost
  float thetaJPsi = MuInJPsiFrame.Angle(JPsiVec.Vect());
  */
  return thetaJPsi;

}

void GetMomentumInMotherFrame( TLorentzVector Mother,TLorentzVector Particle, double BeamEnergy , TVector3 &Particle_rotated){
  //Mother momentum in lab frame
  TVector3 bMother = Mother.BoostVector();
  TLorentzVector beam1(0., 0.,  BeamEnergy, BeamEnergy); // beam momentum in lab frame
  TLorentzVector beam2(0., 0., -BeamEnergy, BeamEnergy); // beam momentum in lab frame
  beam1.Boost(-bMother);      // beam momentum in JPsi rest frame
  beam2.Boost(-bMother);     // beam momentum in JPsi rest frame
  TVector3 beam1_dir = beam1.Vect().Unit();  // beam direction in Mother rest frame
  TVector3 beam2_dir = beam2.Vect().Unit();  // beam direction in Mother rest frame

  TVector3 Y = beam1_dir.Cross( beam2_dir ).Unit(); // the production plane normal
  TVector3 Z = Mother.Vect().Unit();         // Mother direction in lab frame
  TVector3 X = Y.Cross(Z).Unit();         // completes the right-handed coordinate

  Particle.Boost(-bMother);     // Particle momentum in Mother rest frame


  TRotation rotation;
  rotation.SetToIdentity();
  rotation.RotateAxes(X,Y,Z);
  rotation.Invert(); // transforms coordinates from the "xyz" frame to the new frame
  //TVector3 Particle_rotated = Particle.Vect(); // io: particle coordinates in the rotated frame
  Particle_rotated = Particle.Vect();
  //Particle_rotated.Transform(rotation);

}
