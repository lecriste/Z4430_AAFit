#Amplitude Analysis GooFit

REQUIREMENTS
============

- GooFit with the matrix p.d.f. - I have it working in this git repo : https://github.com/AdrianoDee/MyGooFit
- Root version 5.34/28 - surely works with other version, I tested with this
- ... to be completed

INSTRUCTIONS
=============

=To Build
- Setup GooFit
- Edit the Makefile so that the $(GOOFITDIRECTORY) point to the directory where GooFit is installed (I'd add it to gitignore to avoid conflicts)
- ... to be completed

=To Run
- need a .txt dataset written from RooDataSet in the form of e.g. datasetK800k1410.txt for k800 and k1410 k* taken into account
- ./AmplitudeAnalysis with no argument prints help instructions
- ... to be completed

NOTES FOR FURTHER UPDATES
=========================
- dataset generation within GooFit 
- plotting stuff ported to GPU
- adding Hesse DONE
- cleaning mass,gamma,spin inputs in MatrixPdf (constant variables repetition)
- parameters plot DONE
- projection plots on all the four variables DONE
