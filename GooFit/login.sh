#########
#ROOT

#source /opt/exp_soft/root/root_v5.34.07.Linux-slc5_amd64-gcc4.3/bin/thisroot.sh
#ROOT_path=$HOME
# if you want to use Adriano's local ROOT
#ROOT_path=/lustre/home/adrianodif/
#source $ROOT_path/Tools/root/bin/thisroot.sh
# not needed apparently
#export PATH=$PATH:$HOME/Tools/root/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
#export GOOFITDIRECTORY=$ROOT_path/Git/MyGooFit
export GOOFITDIRECTORY=/lustrehome/cristella/work/MyGooFit
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib:/$ROOT_path/Git/MyGooFit/rootstuff/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib:/$GOOFITDIRECTORY/rootstuff/

##########

##################################################################
##CUDA 7.5
##################################################################
export PATH=$PATH:/usr/local/cuda-7.5/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-7.5/lib64
