#!/bin/sh
#
# This is a simple example of a SGE PE mpich2 batch script
#
# request Bourne shell as shell for job
#$ -S /bin/zsh
#
# Your job name 
#$ -N tt_seed
#
# Use current working directory
#$ -cwd -o /nfs/dust/cms/user/volinar/tools/Test_AREA/MCFM-8.0/Bin -e /nfs/dust/cms/user/volinar/tools/Test_AREA/MCFM-8.0/Bin
#
# Join stdout and stderr
#$ -j y
#
# pe request for distributed configuration. Set your number of processors here. 
# -pe mpich2 4
#
# Currently only sld5 support
#$ -l os=sld6
#
# Disk space
#$ -l h_vmem=4G
#
# Run time 
#$ -l h_rt=72:00:00 
#
# Exclude broken nodes (STATUS of 9.11.14):
#$ -l hostname="!(bird212.desy.de|bird170.desy.de)"
#
# Check job parameters, but do not submit job
# -w v
#
# The following is for reporting only. It is not really needed
# to run the job. It will show up in your output file.

#export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export LHAPDFSYS=/afs/desy.de/user/v/volinar/Top_PHYS/Main
export PATH=${PATH}:${LHAPDFSYS}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LHAPDFSYS}/lib
export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
echo " +++++++++++++++++++++++++++++++ + +  + + + + + +  + + + + "
localSetupROOT --rootVersion=5.34.25-x86_64-slc6-gcc48-opt

#CDPATH=.:$HOME:
PRINTER=cmscp2
LPDEST=$PRINTER
export LPDEST PRINTER
export LANGUAGE="C	"
export LC_ALL="C"

#THIS IS IMPORTNAT FOR YOU
module use -a /afs/desy.de/group/cms/modulefiles/
module avail
alias la='ls -a'

############################### ROOT AT CVMFS ###############################
#   !!!!   ROOT version should be the same as gcc   !!!!
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh
gcc --version
#available ROOT versions at /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/ then /x86_64-slc6-gcc'**'-opt
cd /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.19/x86_64-slc6-gcc48-opt/root/
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.19/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh
#cd /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/
#source  /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
alias root='root -l'
cd -

###################################  Compiler ###############################
export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/qcdnum/bin:$PATH
#############################################################################


export RIVET_REF_PATH=/afs/desy.de/user/v/volinar/Top_PHYS/Rivet/ATLAS_ttlep


export PATH=/nfs/dust/cms/user/volinar/tools/Test_AREA/Main/bin/:$PATH
export LD_LIBRARY_PATH=/nfs/dust/cms/user/volinar/tools/Test_AREA/Main/lib/:$LD_LIBRARY_PATH
export LDFLAGS=-L/nfs/dust/cms/user/volinar/tools/Test_AREA/Main/lib/:$LDFLAGS
###export LHAPATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/lhapdf/share/LHAPDF/PDFsets/:$LHAPATH
export MCFMBRIDGE_LDFLAGS="`mcfmbridge-config --ldflags`"
export LHAPATH=/nfs/dust/cms/user/volinar/tools/Test_AREA/Main/share/lhapdf/PDFsets/:$LHAPATH
#
# Wminus_run.sh is operating file for W^- production with process 6 at NAF
#

date

SAVEDIR=/nfs/dust/cms/user/volinar/tools/Test_AREA/MCFM-8.0/Bin
mkdir -p ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed

cp -r /nfs/dust/cms/user/volinar/tools/Test_AREA/MCFM-8.0/ $TMPDIR
cd $TMPDIR/MCFM-8.0/Bin/

echo "******************************"
date
echo $PWD
ls
echo "******************************"

./mcfm input.DAT >log_1run
cp grid-*   ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
cp *.top ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
cp log_*    ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
echo "******************************"
date
echo $PWD
ls
echo "******************************"


echo "TMPDIR = $TMPDIR"
echo "JOB_SCRIPT = $JOB_SCRIPT"
echo "SGE_JOB_SPOOL_DIR = $SGE_JOB_SPOOL_DIR"

echo "i'm still in $PWD"


echo "#############"
ls
echo "#############"

cd ${SAVEDIR}
mv input.DAT ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
mv Proc_ProcNum_Order_Scales_seed.o* ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed	
echo "+++++++++++++"
ls
echo "-------------"

exit 0
