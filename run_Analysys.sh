{
gROOT->ProcessLine(".x myPDF.cxx+");
gInterpreter->AddIncludePath("/cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/roofit/5.34.09-cms5/include");
gROOT->ProcessLine(".x /cmshome/cristella/work/Z_analysis/AA_fit/Analysis.C+");
}
