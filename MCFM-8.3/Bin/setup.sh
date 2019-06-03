	
#export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
#export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK 
#source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

#export AtlasSetup=/afs/cern.ch/atlas/software/dist/AtlasSetup



#alias asetup='source $AtlasSetup/scripts/asetup.sh'
#lsetup asetup
#echo $External

#asetup 16.6.0.1,32,slc5,opt,gcc43,AtlasProduction  --testarea $PWD
#asetup 19.2.5.1,here, cvmfsonly
#asetup '19.2.4.10.2','here, MCProd'
#asetup 17.0.3,here


#localSetup
#localSetupROOT --rootVersion=5.34.25-x86_64-slc6-gcc48-opt

##############################  CERN  ################################
#kinit vldanilo@CERN.CH
######################################################################

#CDPATH=.:$HOME:
PRINTER=cmscp2
LPDEST=$PRINTER
export LPDEST PRINTER
export LANGUAGE="C	"
export LC_ALL="C"

#export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH:/products/dcache/lib	
# cernlib



#THIS IS IMPORTNAT FOR YOU
module use -a /afs/desy.de/group/cms/modulefiles/
#module avail
#module load root5
#AND ONE OF THE LINES BELOW - PLEASE CHOOSE
#module load cmssw/slc6_amd64_gcc493
#module load cmssw/slc6_amd64_gcc
#module load root/5.34
#-------------------------------------------------------------------------------------------------------------------------------------------------------
#alias ls='ls -ltrF --time-style=+"%d.%m.%Y %H:%M" --color=auto -F'
alias la='ls -a'
#alias crabini='source /cvmfs/cms.cern.ch/crab3/crab.sh; crab --version'
#alias proxyini='voms-proxy-init -verify -voms cms:/cms/dcms -valid 96:00'
#alias analysini='cd ~/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src;cmsenv;cd Analysis/MssmHbb/bin'

############################### ROOT AT CVMFS ###############################
#   !!!!   ROOT version should be the same as gcc   !!!!
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/5.1.0/x86_64-slc6/setup.sh
gcc --version
#available ROOT versions at /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/ then /x86_64-slc6-gcc'**'-opt
cd /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.19/x86_64-slc6-gcc48-opt/root
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.19/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh
cd -
#cd /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/
#source  /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
alias root='root -l'
alias ll='ls -altr'


# NOT SURE ABOUT NEXT TWO LINES...
#source /cvmfs/cms.cern.ch/cmsset_default.sh                                                                         
#export SCRAM_ARCH=slc6_amd64_gcc530
#

echo ">>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<"
echo ">>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<"
echo ">>>>>>>>>>>>>>>APPLGRID DEVELOPMENT<<<<<<<<<<<<<<"
echo ">>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<"
echo ">>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<"

############################# xFitter part  ################################# 
#export CURRENTDIR=/nfs/dust/cms/user/volinar/tools/xFitter
#export version=`cat version`
#export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/xfitter-1.2.2/bin:$PATH
#export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/hoppet/bin:$PATH
#export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/applgrid/bin:$PATH
#export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/lhapdf/bin:$PATH
#export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/apfel/bin:$PATH
#export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/mela/bin:$PATH
#export LD_LIBRARY_PATH=$CURRENTDIR/deps/hoppet/lib/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$CURRENTDIR/deps/lhapdf/lib/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$CURRENTDIR/deps/applgrid/lib/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$CURRENTDIR/deps/apfel/lib/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$CURRENTDIR/deps/mela/lib/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$CURRENTDIR/deps/qcdnum/lib/:$LD_LIBRARY_PATH
#############################################################################

###################################  Compiler ###############################
#. /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh
#cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc46-opt/root
#. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh
#cd -
#export PATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/qcdnum/bin:$PATH
#############################################################################

################### MCFM-6.8 + mcfm-bridge ################################## 
#export PATH=/afs/desy.de/user/v/volinar/Top_PHYS/Main/bin/:$PATH
#export LD_LIBRARY_PATH=/afs/desy.de/user/v/volinar/Top_PHYS/Main/lib/:$LD_LIBRARY_PATH
#export LDFLAGS=-L$HOME/Top_PHYS/Main/lib/:$LDFLAGS
#export MCFMBRIDGE_LDFLAGS="`mcfmbridge-config --ldflags`"
#export LHAPATH=/afs/desy.de/user/v/volinar/Top_PHYS/Main/share/lhapdf/PDFsets/:$LHAPATH
#########################  RIVET  ##########################################

#export RIVET_REF_PATH=/afs/desy.de/user/v/volinar/Top_PHYS/Rivet/ATLAS_ttlep

#################### TESTING DIRECTORY ######################################


#export PATH=/afs/desy.de/user/v/volinar/xxl/TESTING_AREA/testMain/bin/:$PATH
#export LD_LIBRARY_PATH=/afs/desy.de/user/v/volinar/xxl/TESTING_AREA/testMain/lib/:$LD_LIBRARY_PATH
#export LDFLAGS=-L/afs/desy.de/user/v/volinar/xxl/TESTING_AREA/testMain/lib/:$LDFLAGS
#export MCFMBRIDGE_LDFLAGS="`mcfmbridge-config --ldflags`"


export PATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/bin/:$PATH
export LD_LIBRARY_PATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/lib/:$LD_LIBRARY_PATH
export LDFLAGS=-L/afs/desy.de/user/v/volinar/APPLgrid_development/Main/lib/:$LDFLAGS
###export LHAPATH=/nfs/dust/cms/user/volinar/tools/xFitter/deps/lhapdf/share/LHAPDF/PDFsets/:$LHAPATH
export MCFMBRIDGE_LDFLAGS="`mcfmbridge-config --ldflags`"

export LHAPATH=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/share/LHAPDF/:$LHAPATH
#export LHAPATH=/nfs/dust/cms/user/volinar/tools/Test_AREA/Main/share/lhapdf/PDFsets/:$LHAPATH

#export PATH=/nfs/dust/cms/user/volinar/tools/DELETE_ME/Main/bin/:$PATH
#export LD_LIBRARY_PATH=/nfs/dust/cms/user/volinar/tools/DELETE_ME/Main/lib/:$LD_LIBRARY_PATH
#export LDFLAGS=-L/nfs/dust/cms/user/volinar/tools/DELETE_ME/Main/lib/:$LDFLAGS
#export MCFMBRIDGE_LDFLAGS="`mcfmbridge-config --ldflags`"

alias Go2nfs='/nfs/dust/cms/user/volinar/tools/Test_AREA/'
alias Bird2nfs='/nfs/dust/cms/user/volinar/BIRD/MCFM-8.0/'
alias condor_q='condor_q | grep volinar'
