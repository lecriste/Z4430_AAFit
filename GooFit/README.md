#Amplitude Analysis GooFit

REQUIREMENTS
============

- GooFit with the matrix p.d.f. - I have it working in this git repo : https://github.com/AdrianoDee/MyGooFit
- Root version 5.34/28 - surely works with other version, I tested with this
- ... to be completed

INSTRUCTIONS
=============

= To Build

- Setup GooFit

- Edit the Makefile so that the $(GOOFITDIRECTORY) point to the directory where GooFit is installed (I'd add it to gitignore to avoid conflicts)

- ... to be completed

= To Run

- if "-r" option :need a .txt dataset written from RooDataSet in the executable path in the form of e.g. datasetK800k1410.txt for k800 and k1410 k* taken into account; unless path provided. The txt dataset must be in the form:

      0.709429 0.561751 -0.984957 -1.25312
      0.767005 0.756151 -0.493894 -2.07288
      0.883276 -0.22407 0.0196103 -0.176145
      1.1992 0.144528 0.881793 0.781775
      0.93699 0.904563 -0.0396737 0.297077
      0.846126 0.685711 0.880613 -0.62682

in order : MassKPi, CosJPsi, CosKStar,Phi

- ./AmplitudeAnalysis with no argument prints help instructions

- ... to be completed

NOTES FOR FURTHER UPDATES
=========================
- cleaning mass,gamma,spin inputs in MatrixPdf (constant variables repetition)
- plotting stuff ported to GPU
- perform studies to chose the right integration step
- adding Hesse DONE
- parameters plot DONE
- projection plots on all the four variables DONE
- dataset generation within GooFit DONE (but still CPU-side!)
