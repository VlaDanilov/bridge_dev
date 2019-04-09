#!/bin/zsh
#
# This is a simple example of a SGE PE mpich2 batch script
#
# request Bourne shell as shell for job
#$ -S /bin/zsh
#
# Your job name 
# -N tt_seed
#
# Use current working directory
# -cwd -o /nfs/dust/cms/user/volinar/tools/Test_AREA/MCFM-8.0/Bin -e /nfs/dust/cms/user/volinar/tools/Test_AREA/MCFM-8.0/Bin
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
#$ -l h_vmem=8G
#
# Run time 
#$ -l h_rt=23:00:00 
#
# Exclude broken nodes (STATUS of 9.11.14):
#$ -l hostname="!(bird212.desy.de|bird170.desy.de)"
#
# Check job parameters, but do not submit job
# -w v
#
# The following is for reporting only. It is not really needed
# to run the job. It will show up in your output file.

#export LHAPDFSYS=/afs/desy.de/user/v/volinar/Top_PHYS/Main
#export PATH=${PATH}:${LHAPDFSYS}/bin
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LHAPDFSYS}/lib
#export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK

#module use -a /afs/desy.de/group/cms/modulefiles/
#module avail
#alias la='ls -a'

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
#export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/qcdnum/bin:$PATH
#############################################################################


#export RIVET_REF_PATH=/afs/desy.de/user/v/volinar/Top_PHYS/Rivet/ATLAS_ttlep


export PATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/bin/:$PATH
export LD_LIBRARY_PATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/lib/:$LD_LIBRARY_PATH
export LDFLAGS=-L/afs/desy.de/user/v/volinar/APPLgrid_development/Main/lib/:$LDFLAGS
###export LHAPATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/lhapdf/share/LHAPDF/PDFsets/:$LHAPATH
export MCFMBRIDGE_LDFLAGS="`mcfmbridge-config --ldflags`"
#export LHAPATH=/nfs/dust/cms/user/volinar/tools/Test_AREA/Main/share/lhapdf/PDFsets/:$LHAPATH
export LHAPATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/share/LHAPDF/$LHAPATH
#
# Wminus_run.sh is operating file for W^- production with process 6 at NAF
#

SAVEDIR=/afs/desy.de/user/v/volinar/APPLgrid_development/BIRD_results
MCFMDIR=/afs/desy.de/user/v/volinar/APPLgrid_development/MCFM-6.8
#mkdir -p ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed

cp -r $MCFMDIR ${TMP}
echo "TMPDIR is $TMP"
echo "Switching to $TMP/MCFM-8.0/Bin"
cd ${TMP}/MCFM-8.0/Bin/

echo "******************************"
echo "Start at $date"
echo "In $PWD"
echo ""
echo "Current setup:"
echo "PATH=$PATH"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "LDFLAGS=$LDFLAGS"
echo "MCFMBRIDGE_LDFLAGS=`mcfmbridge-config --ldflags`"
echo "LHAPATH=$LHAPATH"
echo "******************************"
        rm grid*
        rm *.top
        rm *.grid
        rm *.C
        rm grid-*
echo "******************************"
echo "*****STARTING FIRST RUN*******"
echo "******************************"
./mcfm input.DAT 
echo "******************************"
echo "*****STARTING SECOND RUN******"
echo "******************************"
./mcfm input.DAT 
cp grid-*   ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
cp *.top ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
cp NONEXISTEN_FILE ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
echo "******************************"
echo "Work is done at $date"
echo "In $PWD"
ls
echo "******************************"
echo ""

echo "TMPDIR = $TMP"
echo "JOB_SCRIPT = $JOB_SCRIPT"
echo "SGE_JOB_SPOOL_DIR = $SGE_JOB_SPOOL_DIR"

echo "i'm still in $PWD"


echo "#############"
ls
echo "#############"

cd ${MCFMDIR}/Bin
mv input.DAT ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
echo "+++++++++++++"
ls
echo "-------------"

exit 0
