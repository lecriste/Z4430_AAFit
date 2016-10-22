CUDA_PATH       ?= /usr/local/cuda-7.5
CXX= nvcc -ccbin g++
LD=g++
OutPutOpt = -o

CXXFLAGS     = -O3 -arch=sm_35 -g
CUDALIBDIR=lib64
CXXFLAGS    += -Xcompiler -fPIC -Xlinker -fPIE
#CXXFLAGS    += -Xcompiler -Xlinker
#CXXFLAGS    += -XLinker -fPIE

UNAME=$(shell uname)
ifeq ($(UNAME), Darwin)
CUDALIBDIR=lib64
CXXFLAGS+=-m64
else
# Adds library symbols to debug info 
# - not strictly required, nice to have.
# Turn it off if it gives you trouble. 
#LDFLAGS += -rdynamic 
endif

##CXXFLAGS    += -Xcompiler -fopenmp -lgomp
LDFLAGS     += -fPIC -fPIE #-fopenmp #-lgomp
#LD     += -XLinker -fPIE
ifneq ($(TARGET_OMP),)
CXXFLAGS    += -Xcompiler -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP #-lgomp
LDFLAGS     += -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP #-lgomp
endif

ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda-7.5
#CUDALOCATION = /cmshome/adrianodif/CUDAs/cuda-6.0
#CUDALOCATION = /cmshome/adrianodif/CUDAs/cuda-6.5
#CUDALOCATION = /usr/local/cuda-6.5

endif
CUDAHEADERS = $(CUDALOCATION)/include/
GOOFITDIR = $(GOOFITDIRECTORY)

INCLUDES += -I$(CUDALOCATION)/samples/common/inc/ -I$(CUDAHEADERS) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff/ -I$(GOOFITDIR)/rootstuff/Math -I$(GOOFITDIR)/PDFs/
LIBS += -L$(CUDALOCATION)/$(CUDALIBDIR) -lcuda -lcudart -L$(GOOFITDIR)/rootstuff/ -lRootUtils

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I$(ROOTSYS)/include/ -I$(ROOTSYS)/include/Math
ROOT_LIBS     = -L$(ROOTSYS)/lib/ -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit -lFoam # -lRooFit -lRooFitCore  -lRooStats -lRootAuth  
WRKDIR = $(GOOFITDIR)/wrkdir/

THRUSTO         = $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o
ROOTRIPDIR      = $(GOOFITDIR)/rootstuff/
ROOTUTILLIB     = $(ROOTRIPDIR)/libRootUtils.so
PROGRAMS        = AmplitudeAnalysis 

.SUFFIXES: 

all:    $(PROGRAMS)

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

AmplitudeAnalysis:	AmplitudeAnalysis.o
			$(LD) $(LDFLAGS) $^ $(THRUSTO) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
			@echo "$@ done"


clean:
			@rm -f *.o core $(PROGRAMS) 
