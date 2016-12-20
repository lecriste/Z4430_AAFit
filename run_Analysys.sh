{
gSystem->Exec("rm Analysis_C*");

//TString RooFit_includePath = TROOT::GetMacroPath();
//gInterpreter->AddIncludePath("/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms5/include");
gROOT->ProcessLine(".x myPDF.cxx+");
gROOT->ProcessLine(".x Dalitz_contour.cxx+");
gROOT->ProcessLine(".x effMasses.cxx+");
gROOT->ProcessLine(".x Angles_contour.cxx+");
gROOT->ProcessLine(".x sqDalitzToMassesPdf.cxx+");
gROOT->ProcessLine(".x sqDalitz_contour.cxx+");
gROOT->ProcessLine(".x Analysis.C+");
}
