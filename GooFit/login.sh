
#########
#ROOT

#source /opt/exp_soft/root/root_v5.34.07.Linux-slc5_amd64-gcc4.3/bin/thisroot.sh
#source $HOME/Tools/root/bin/thisroot.sh
source /lustre/home/adrianodif/Tools/root/bin/thisroot.sh
export PATH=$PATH:$HOME/Tools/root/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib:/lustre/home/adrianodif/Git/MyGooFit/rootstuff/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export GOOFITDIRECTORY=/lustre/home/adrianodif/Git/MyGooFit

##########

##################################################################
##CUDA 7.5
##################################################################
export PATH=$PATH:/usr/local/cuda-7.5/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-7.5/lib64
