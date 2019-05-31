#!/bin/zsh
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

export PATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/bin/:$PATH
export LD_LIBRARY_PATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/lib/:$LD_LIBRARY_PATH
export LDFLAGS=-L/afs/desy.de/user/v/volinar/APPLgrid_development/Main/lib/:$LDFLAGS
export MCFMBRIDGE_LDFLAGS="`mcfmbridge-config --ldflags`"
export LHAPATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/share/lhapdf/PDFsets/:$LHAPATH
#
# Wminus_run.sh is operating file for W^- production with process 6 at NAF
#
date
source /cvmfs/cms.cern.ch/cmsset_default.sh

MCFM=MCFM-8.0

SAVEDIR=/afs/desy.de/user/v/volinar/APPLgrid_development/BIRD_results/${MCFM}
MCFMDIR=/afs/desy.de/user/v/volinar/APPLgrid_development/${MCFM}/
#mkdir -p ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed

cp -r $MCFMDIR ${TMP}
echo "TMPDIR is $TMP"
echo "I'm in ${PWD}"
cd ${TMP}/${MCFM}/Bin/     #   <-----  MCFM VERSION MUST AGREE !!!!

echo "******************************"
date
echo $PWD
ls

rm grid*
rm *.top
rm *.grid
rm *.C

./mcfm input.DAT 
./mcfm input.DAT 

cp grid-*          ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
cp *.top           ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed
cp NONEXISTEN_FILE ${SAVEDIR}/Process_ProcNum/DynamicScale_Scales/Order/seed_numberoftheseed

echo "******************************"
date

echo $PWD
ls
echo "******************************"


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
