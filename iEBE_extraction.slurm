#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=0-02:00:00
#SBATCH --array=1-100
#SBATCH --output=iEBE_extraction.o%A.out

FOLDER=centrality_50_60/vis_20/
JOB=$SLURM_ARRAY_TASK_ID

cd /store/user/palmerjh/Results/${FOLDER}
unzip job-$JOB.zip 2> /dev/null

cd job-$JOB

/scratch/palmerjh/iEBE/EBE-Node/EbeCollector/unpackDatabase.py particles.db

cd /home/palmerjh/CMSSW_7_5_0_pre5/src/

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472
eval `scramv1 runtime -sh`

cp /home/palmerjh/makeTree.C ./iEBE_trees/${FOLDER}
cd ./iEBE_trees/${FOLDER}

root -q -b "makeTree.C($JOB)" &> iEBE_extraction_${JOB}.log

exit 0
