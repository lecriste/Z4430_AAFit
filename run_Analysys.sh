{
gInterpreter->AddIncludePath("/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms5/include");
gROOT->ProcessLine(".x myPDF.cxx+");
gROOT->ProcessLine(".x Dalitz_contour.cxx+");
gROOT->ProcessLine(".x Analysis.C+");
}
