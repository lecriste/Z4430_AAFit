########
#GIT

function gitStart(){

eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_git

}

alias gpull='git pull && gmake'

#########
#ROOT

#source /opt/exp_soft/root/root_v5.34.07.Linux-slc5_amd64-gcc4.3/bin/thisroot.sh 
source $HOME/Tools/root/bin/thisroot.sh
export PATH=$PATH:$HOME/Tools/root/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib:/lustre/home/adrianodif/Git/MyGooFit/rootstuff/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export GOOFITDIRECTORY=/lustre/home/adrianodif/Git/MyGooFit

##########

#module load cuda

##################################################################
##CUDA 7.5
##################################################################
export PATH=$PATH:/usr/local/cuda-7.5/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-7.5/lib64

##################################################################
#MPS
##################################################################

NGPUS=1

function ClearAllMPS (){
        echo "=============================================================="
        echo "Clearing MPS from all GPUs"
        echo "=============================================================="
        for ((i=0; i<$NGPUS; i++))
        do
        echo "--GPU " $i "Cleared"
        #export CUDA_MPS_PIPE_DIRECTORY=/tmp/mps_$i
        echo "quit" | nvidia-cuda-mps-control
        rm -r /tmp/mps_$i
        rm -r /tmp/mps_log_$i
        done


}

function SetUpMPS(){
        echo "=============================================================="
        echo "Setting Up  MPS from all GPUs"
        echo "=============================================================="

        for ((i=0; i<$NGPUS; i++))
        do
        echo "--GPU " $i "MPS Set Up"
         mkdir /tmp/mps_$i
         mkdir /tmp/mps_log_$i
        export CUDA_VISIBLE_DEVICES=$i
        export CUDA_MPS_PIPE_DIRECTORY=/tmp/mps_$i
        export CUDA_MPS_LOG_DIRECTORY=/tmp/mps_log_$i
        nvidia-cuda-mps-control -d
        done
}

#alias FourMuonsToys='/lustre/home/adrianodif/Git/FourMuonsToys/FourMuonsToys/FourMuonToys'

function BoundFourSeq(){

        echo "=============================================================="
        echo "Running FourMuons Process for $1 times - One CPU bound"
        echo "=============================================================="
        for ((i=0;i<$1;i++))
        do
        (taskset -c $i ./FourMuonsToys -n $2 -x -y -z&);
        sleep 0.5;
        done
         echo "=============================================================="
        echo "FourMuonsToys $2 Process for $1 times DONE"
        echo "=============================================================="

}

alias killemall='kill `ps -o pid= -N T`'

