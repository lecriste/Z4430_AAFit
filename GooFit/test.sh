#PBS -N testAA
#PBS -M adriano.diflorior@ba.infn.it
#PBS -S /bin/bash
#PBS -d /lustre/home/adrianodif/Git/FourMuonsToys/FourMuonsToys/
#PBS -e /lustre/home/adrianodif/testAA.err
#PBS -o /lustre/home/adrianodif/testAA.out
#PBS -q bigmpi2@sauron.recas.ba.infn.it
#PBS -l gpus=1 
#PBS -l ncpus=1

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
let HOSTNAME = "$(hostname)"
echo HOST: "$(hostname)"
echo -----------------------------------------------------

module load cuda/cuda-75 

#CUDA + GOOFIT
source /lustre/home/adrianodif/login.sh

date

cd /lustre/home/adrianodif/Git/Z4430_AAFitDif

gmake

time ./Test -n 10000

date


wait



wait
nvidia-smi
wait
