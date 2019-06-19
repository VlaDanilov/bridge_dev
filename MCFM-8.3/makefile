# Flag for compiling with MPI (YES) or not (anything else)
USEMPI = NO

# Set this to NATIVE/PDFLIB/LHAPDF
#   NATIVE -- internal routines
#   LHAPDF -- Les Houches library
PDFROUTINES = LHAPDF

# Replace this with the location of Cernlib on your system (if desired)
CERNLIB     =
# Replace this with the location of LHAPDF on your system (if desired)
LHAPDFLIB=/afs/desy.de/user/v/volinar/APPLgrid_development/Main/lib 

LDFLAGS += -u book_gridptr $(shell mcfmbridge-config --ldflags)

FASTJETINC  = $(shell fastjet-config --cxxflags)
FASTJETLIB  = $(shell fastjet-config --libs)


# Set this to NO/YES/FROOT
#   NO  -- no n-tuple output or unweighting is possible
#   YES -- n-tuples, unweighting available - requires CERNLIB libraries
#   FROOT -- n-tuples, unweighting available -- requires ROOT libraries;
#            output uses the FROOT package of P. Nadolsky.
NTUPLES = NO

OSXGFORTRAN =
OSXGPP =

ifeq ($(origin FC),environment)
    $(info Inheriting FC, CXX from environment)
    $(info FC = $(FC))
    $(info CXX = $(CXX))
else
	ifeq ($(USEMPI),YES)
	  FC = mpifort
	  CXX = mpic++
	  # for OpenMPI
	  export OMPI_CXX=g++
	  export OMPI_FC=gfortran
	else
	  FC = gfortran
	  CXX = g++
	endif

	ifeq ($(OSXGFORTRAN),)
	else
		FC=$(OSXGFORTRAN)
	endif

	ifeq ($(OSXGPP),)
	else
		CXX=$(OSXGPP)
	endif
endif



ifeq ($(USEMPI),YES)
  MAIN = mcfm_mpi.o
  MPIDUMMY=.
else
  MAIN = mcfm_omp.o
  MPIDUMMY=$(PWD)/src/mpidummy
endif

MCFMHOME        = $(PWD)
SOURCEDIR       = $(MCFMHOME)/src
OBJNAME = obj
OBJDIR		= $(MCFMHOME)/$(OBJNAME)
OUTPUT_OPTION = -o $(OBJDIR)/$@
BIN		= $(MCFMHOME)/Bin
INCPATH  	= $(SOURCEDIR)/Inc
TENSORREDDIR	= $(MCFMHOME)/TensorReduction
PVDIR		= $(TENSORREDDIR)/pv
PVEXTDIR	= $(TENSORREDDIR)/pvext
RECURDIR	= $(TENSORREDDIR)/recur
OVDIR		= $(TENSORREDDIR)/ov
MOD_QCDLOOP=$(MCFMHOME)/qcdloop-2.0.2
VPATH		= $(DIRS):$(MOD_QCDLOOP):$(INCPATH):$(TENSORREDDIR)/Include:$(SOURCEDIR)/mpidummy

# Flags for compilation
DEVELOPMENT = -Wall -Wextra -Wno-unused-variable -Wno-compare-reals -Wno-unused-parameter -Wno-unused-dummy-argument -Wno-unused-label -Wuninitialized -Wno-error=conversion -Wno-tabs -O0
# do not use -Ofast with gcc, this enables unsafe math operations (-ffast-math) and will break things!
PRODUCTION = -finit-local-zero -O2 -mtune=native
FFLAGS 	= $(PRODUCTION) -fno-f2c -ffixed-line-length-none -fopenmp \
 -I$(INCPATH) -I$(MPIDUMMY) -I$(TENSORREDDIR)/Include -L$(PVDIR) -L$(PVEXTDIR) -J$(OBJNAME) -I$(MOD_QCDLOOP)


# If using FROOT package for ROOT ntuples, first specify C++ compiler:
CXXFLAGS=$(CXXFLAGS0) -fopenmp -std=c++11 -Wall $(DROOT)
# ROOTLIBS and ROOTINCLUDE are locations of ROOT libraries and header files.
# Find the ROOT libraries and include directory automatically by running root-config 
# or specify them manually by providing values for ROOTLIBS and ROOTINCLUDE below
DMYROOT= -DMYROOT
ROOTLIBS     := $(shell root-config --libs 2> /dev/null)
ROOTINCLUDE  := $(shell root-config --cflags 2> /dev/null)


DIRS	=	$(MCFMHOME):\
		$(OBJDIR):\
		$(SOURCEDIR)/User:$(SOURCEDIR)/Procdep:$(SOURCEDIR)/Vol:\
		$(SOURCEDIR)/Need:$(SOURCEDIR)/Lib:$(SOURCEDIR)/Phase:\
		$(SOURCEDIR)/Parton:$(SOURCEDIR)/Integrate:\
		$(SOURCEDIR)/Spinor:$(SOURCEDIR)/Spinor_Weyl:\
		$(SOURCEDIR)/Wbb:$(SOURCEDIR)/Zbb:\
		$(SOURCEDIR)/WHbbar:$(SOURCEDIR)/ZHbbar:\
		$(SOURCEDIR)/WW:$(SOURCEDIR)/WZ:$(SOURCEDIR)/ZZ:\
		$(SOURCEDIR)/Diboson: \
		$(SOURCEDIR)/Wgajet:\
		$(SOURCEDIR)/Top:$(SOURCEDIR)/Topdk:$(SOURCEDIR)/Singletop:\
		$(SOURCEDIR)/Singletop2:\
		$(SOURCEDIR)/helas:\
		$(SOURCEDIR)/TopH:$(SOURCEDIR)/TopZ:$(SOURCEDIR)/TopW:\
		$(SOURCEDIR)/HWW:$(SOURCEDIR)/HZZ:$(SOURCEDIR)/Tau:\
		$(SOURCEDIR)/Httbar:$(SOURCEDIR)/Hjetmass\
		$(SOURCEDIR)/Hjetmass/sushi\
		$(SOURCEDIR)/Hjetmass/split\
		$(SOURCEDIR)/W:$(SOURCEDIR)/Z:\
		$(SOURCEDIR)/qpW:\
		$(SOURCEDIR)/W1jet:$(SOURCEDIR)/Z1jet:\
		$(SOURCEDIR)/Wcjet:$(SOURCEDIR)/Wbfrmc:$(SOURCEDIR)/Wbjet:\
		$(SOURCEDIR)/W2jet:$(SOURCEDIR)/W2jetvirt:\
		$(SOURCEDIR)/Wbbm:$(SOURCEDIR)/Zbbm:\
		$(SOURCEDIR)/Wgam:$(SOURCEDIR)/Wgamgam:$(SOURCEDIR)/Zgam:\
                $(SOURCEDIR)/Z2jet:$(SOURCEDIR)/Zb:\
		$(SOURCEDIR)/bbHiggs:$(SOURCEDIR)/Wt:\
                $(SOURCEDIR)/qqH:$(SOURCEDIR)/qqHWW:$(SOURCEDIR)/qqHZZ:\
		$(SOURCEDIR)/ZQjet:\
                $(SOURCEDIR)/ggH:$(SOURCEDIR)/ggHg:\
                $(SOURCEDIR)/ggHggreal:$(SOURCEDIR)/ggHggvirt:\
                $(SOURCEDIR)/H4pCode:\
                $(SOURCEDIR)/WHWW:$(SOURCEDIR)/ZHWW:\
		$(SOURCEDIR)/WHZZ:$(SOURCEDIR)/ZHZZ:\
		$(SOURCEDIR)/WHgaga:$(SOURCEDIR)/ZHgaga:$(SOURCEDIR)/qqHgaga:\
		$(SOURCEDIR)/Stopb:$(SOURCEDIR)/Stopjet:\
		$(SOURCEDIR)/epem3j:\
		$(SOURCEDIR)/WpWp2j:$(SOURCEDIR)/F90\
		$(SOURCEDIR)/qqZtt:$(SOURCEDIR)/ggHgaga:\
		$(SOURCEDIR)/ggHZga:$(SOURCEDIR)/HH:\
		$(SOURCEDIR)/Gamgam:$(SOURCEDIR)/Dirgam:\
		$(SOURCEDIR)/TopdkBSY:$(SOURCEDIR)/Topdecay:\
		$(SOURCEDIR)/Frag:$(SOURCEDIR)/Zgamgam:$(SOURCEDIR)/Zgamjet:\
		$(SOURCEDIR)/Zaj:\
		$(SOURCEDIR)/SingletopH:$(SOURCEDIR)/SingletopZ:\
		$(SOURCEDIR)/WpmZbj:$(SOURCEDIR)/DM:$(SOURCEDIR)/Monojet:\
                $(SOURCEDIR)/Vec_DM:$(SOURCEDIR)/GG_DM:\
                $(SOURCEDIR)/Ax_DM:$(SOURCEDIR)/Scal_DM:\
                $(SOURCEDIR)/PS_DM:$(SOURCEDIR)/Trigam:\
                $(SOURCEDIR)/VV:$(SOURCEDIR)/Fourgam:\
                $(SOURCEDIR)/SCET:$(SOURCEDIR)/TDHPL:\
                $(SOURCEDIR)/SCET0j:\
                $(SOURCEDIR)/Trigam/Mad:\
                $(SOURCEDIR)/pwgplots:$(SOURCEDIR)/Multichan:$(SOURCEDIR)/TopH/Store/Mad\
                $(SOURCEDIR)/UTools:$(SOURCEDIR)/WBFZZ\
                $(SOURCEDIR)/WBFWW:$(SOURCEDIR)/WBFWpWp:$(SOURCEDIR)/WBFWZ\
                $(SOURCEDIR)/WH1jet:$(SOURCEDIR)/ZH1jet:$(SOURCEDIR)/QT:$(SOURCEDIR)/Mad\
                $(SOURCEDIR)/ZEW:$(SOURCEDIR)/TopEW:$(SOURCEDIR)/JetsEW\
				$(SOURCEDIR)/Mods


# -----------------------------------------------------------------------------
# Specify the object files. 

EWFILES = \
mcfmccdb0.o \
qqb_z_ew_sudakov.o \
z_ew_box.o \
sudakov_DY.o \
scalar_int.o \
charge_op.o \
qqb_z_ew_exact.o \
self_VV_new.o \
ggdilep.o

EWFILES += \
ggQQb_ew_oneloop.o \
qqbQQb_ew_oneloop.o \
qqb_QQb_ew_sudakov.o \
qqb_QQb_mix.o \
qqb_QQb_mix_g.o \
qqb_QQb_mix_gs.o \
qqb_QQb_mix_v.o \
qqb_QQb_mix_z.o \

EWFILES += \
single.o \
dijet_gg_ew.o \
dijet_qqb_ew_bx.o \
dijet_qqb_ew_tree.o \
dijet_qqb_ew_vert.o \
qqb_twojet_ew.o \
qqb_twojet_ew_sudakov.o \
qqb_twojet_ew_exact.o \
twojet_ew_tree.o \
qqb_twojet_mix.o \
qqb_twojet_mix_g.o \
qqb_twojet_mix_gs.o \
qqb_twojet_mix_gqsub.o \
qqb_twojet_mix_gqsub2.o \
qqb_twojet_mix_z.o

WH1JETFILES = \
WHqqbgg.o \
qqb_WH1jet.o \
qqb_WH1jet_g.o \
qqb_WH1jet_gs.o \
qqb_WH1jet_gvec.o \
qqb_WH1jet_v.o \
qqb_WH1jet_z.o \
qqb_ZH1jet.o \
qqb_ZH1jet_g.o \
qqb_ZH1jet_gs.o \
qqb_ZH_gvec.o \
qqb_ZH1jet_v.o \
qqb_ZH1jet_z.o \
A51_VH.o\
A52_VH.o\
virt5_VH.o\
virt5_VHtop.o\
A5NLO_VH.o\
A5NLO_VHtop.o\
qqbgWH_topamp.o\
virt5_ZHtop.o\
A5NLO_ZHtop.o\
qqbgZH_topamp.o\
genclust_kt_ab.o\
combine_ab.o\

SCETFILES = \
ap.o \
auxfunctions.o \
gen3taucut.o \
gen4taucut.o \
genparton.o \
genpartonVj.o \
genVgataucut.o \
genVgajtaucut.o \
genVH.o \
genVHjjtaucut.o \
genVHjtaucut.o \
genVjtaucut.o \
hardgg.o \
hardqq.o \
hplog.o \
I1gg.o \
I1qq.o \
I2gg.o \
I2qq.o \
plus.o \
qqcoeff.o \
setupscet.o \
xbeam1bis.o \
xbeam2bis.o \
xcdil.o \
xspenz.o

SCET0jFILES = \
assemble.o \
gamgamampsq.o \
gamgamampsq_new.o \
lumxmsq_gaga.o \
lumxmsq_h.o \
lumxmsq_h_Zga.o \
lumxmsq_w.o \
lumxmsq_wh.o \
lumxmsq_z.o \
lumxmsq_zh.o \
qqb_gamgam_vbis.o \
softggbis.o \
softqqbis.o \
dxpdf_dfridr.o \
midpoint_devpdf.o \
tau0_powcorr_gg.o \
tau0_powcorr_qa.o

SCET1jFILES = \
ampgggpmm.o \
ampgggppp.o \
assemblejet.o \
averageoverW.o \
averageoverZ.o \
computeIijm.o \
gg_hgbis.o \
gg_hg_vbis.o \
Hampgggsq.o \
hel_H.o \
hel_V.o \
I0JSTW.o \
I1JSTW.o \
jet.o \
lumxmsq_h1jet.o \
lumxmsq_w1jet.o \
lumxmsq_z1jet.o \
qqb_w1jet_vbis.o \
qqb_z1jet_vbis.o \
soft_ab_ggg.o \
soft_ab_qag.o \
soft_ab_qgq.o \
soft_nab_ggg.o \
soft_nab_qag.o \
soft_nab_qgq.o \
Wampqqbgsq.o \
Zampgggsq.o \
Zampqqbgsq.o

QTFILES = \
qtint.o \
coeffqt.o \
auxfuncqt.o

WBFFILES = \
getvbfpoint.o \
comparevbf.o \
qq_WZqqstrong.o \
qqWZggampf.o \
qqWZggamp.o \
jonezstrong.o \
jZexch2strong.o \
jww.o \
jWexch2.o \
jZexch2.o \
ampmidWZ.o \
jonewWZ.o \
qq_WZqq.o \
jonez.o \
qq_WWss.o \
qq_WWssstrong.o \
ampBdiags.o \
ampmidWpWp.o \
jtwodiagsWpWp.o \
jonewWpWp.o \
ampvbf.o \
jtwoWW.o \
jtwoWWstrong.o \
jZWZ.o \
jtwodiags.o \
jtwoWexch.o \
jonew.o \
jonewstrong.o \
ampmidWW.o \
qq_VVqq.o \
qq4lggampf.o \
qq4lggamp.o \
qq2l2nuggamp.o \
qq_WWqq.o \
qq_WWqqstrong.o \
ampmid.o \
nplotter_qqZZqq.o \
scaleset_s-hat.o \
setupmom.o \
setupzprops.o \
averageovertwoZ.o \
pgen.o \
pgen2.o \
qq_ZZqq.o \
qq_ZZqq_both.o \
qq_ZZqqstrong.o \
WWZZ.o \
ZZHZZamp.o \
spinorcurr.o \
jzero.o \
jone.o \
jtwo.o \
jtwo3456.o \
ZZSingleres.o \
ZZ_HZZ1.o

UTOOLSFILES=\
defpeta_bubs.o\
gen_del1del2_bub3m.o\
gen_kflat_spins.o\

MULTICHANFILES = \
singcheck.o \
dipoleconfig.o \
getperp.o \
multichan.o

FOURGAMFILES= \
qqb_fourgam.o \
qqb_fourgam_fragdips.o \
qqb_fourgam_frag.o \
qqb_fourgam_g.o \
qqb_fourgam_gs.o \
qqb_fourgam_v.o \
qqb_fourgam_z.o \
real_aaaaj.o \
real_aaajj_fill.o \
aaaa_fill.o \
aaaa_MHV_b_init.o \
aaaa_MHV_c_init.o \
aaaa_MHV_d_init.o \
aaaa_MHV_r_init.o \
aaaa_NMHV_b_init.o \
aaaa_NMHV_c_init.o \
aaaa_NMHV_d_init.o \
aaaa_NMHV_r_init.o \
Bubint_init.o \
Triint_init.o \
Boxint_init.o \

TRIGAMFILES=\
checkgvec.o \
gmgmjetn.o \
qqb_gmgmjt.o\
amp_2gam1g.o\
virt_gmgmjt.o\
qqb_gmgmjt_gs.o \
qqb_gmgmjt_v.o\
qqb_gmgmjt_z.o\
qqb_gmgmjt_gvec.o \
qqb_gmgmjt_frag_comb.o \
ampsq_2gam2q.o \
amp_2gam2q.o \
ampsq_2gam2g.o \
amp_2gam2g.o \
ampsq_3gam1g.o \
virt_trigam.o\
kinem_trigam.o\
amp_3gam1g.o \
qqb_dirgam_g_swap.o \
qqb_gmgmjt_g.o \
qqb_trigam.o \
qqb_trigam_frag.o\
qqb_trigam_fragdips.o \
qqb_trigam_g.o \
qqb_trigam_gs.o \
qqb_trigam_v.o \
qqb_trigam_z.o

PWGPLOTSFILES = \
boostrot.o \
pwhg_analysis_generic.o \
pwhg_analysis_KN.o \
pwhg_analysis_ST_sch_dk.o \
pwhg_analysis_ST_tch_dk.o \
pwhg_analysis_ST_wt_dk.o \
pwhg_bookhist-new.o \
pwhg_analysis_Z.o \
pwhgplotter.o \


TOPDKBSYFILES = \
A1Hggppmp.o \
A1Hggpppp.o \
A1Hqqppmp.o \
A1fggppmp.o \
A1fggpppp.o \
A1fqqppmp.o \
A1lcqqppmp.o \
A1slcqqppmp.o \
A41ggppmp.o \
A41ggpppp.o \
A41qqppmp.o \
A43ggppmp.o \
A43ggpppp.o \
ALggppmp.o \
ALggpppp.o \
ALslcggpppm.o \
ALslcggpppp.o \
ARggppmp.o \
ARggpppp.o \
BSYA0ggppmp.o \
BSYA0ggpppp.o \
BSYA0qedpppm.o \
BSYA0qedpppp.o \
BSYA0qqppmp.o \
BSYA1Hggppmp.o \
BSYA1Hggpppp.o \
BSYA1Hqqppmp.o \
BSYA1fggppmp.o \
BSYA1fggpppp.o \
BSYA1fqqppmp.o \
BSYA1lcqqppmp.o \
BSYA1slcqqppmp.o \
BSYALggppmp.o \
BSYALggppmphp.o \
BSYALggpppp.o \
BSYALggpppphp.o \
BSYALslggpppm.o \
BSYALslggpppp.o \
BSYARggppmp.o \
BSYARggppmphp.o \
BSYARggpppp.o \
BSYARggpppphp.o \
gluegen.o \
integralfill.o \
qqbgen.o \
spinorextend.o

TOPDECAYFILES = \
tdecayg.o \
adecayg.o \
tdecayW_v.o \
adecayW_v.o \
tdecayWg.o \
adecayWg.o \
coefsdkmass.o \
coefswdk.o \
lotopdecaywidth.o \
nloratiotopdecay.o \
Gamma0.o \
Gamma0int.o \
asGamma1.o \
asGamma1int.o \
wtransform_generic.o \
tdecay.o \
adecay.o \
tdecay_v.o \
adecay_v.o 

WBFROMCFILES = \
qqb_wbfromc.o \
qqb_wbfromc_v.o \
qqb_wbfromc_z.o \
qqb_wbfromc_g.o \
qqb_wbfromc_gvec.o \
qqb_wbfromc_gs.o 

FRAGFILES = \
GGdR_frag.o \
gdrg_NLO_frag.o\
NP_FragSetI.o \
NP_FragSetII.o \
Pert_Frag.o \
dipolesfrag.o \
dipolesfragx.o \
dopols.o \
fragdriver.o \
int_dips_ga.o \
locatemcfm.o \
transformfrag.o

WPWP2JFILES = \
qqqqampl.o \
qqb_wpwp_qqb.o\
qqb_wpwp_qqb_g.o \
qqqqgampl.o

F90FILES = \
consts_dp.o \
spinfns.o \
recurrenceA.o \
recurrenceB.o \
recurrenceC.o \
recurrence.o 

STOPBFILES = \
extend_trans_stopb.o \
qg_tbq.o \
qg_tbq_g.o \
qg_tbq_gs.o \
qg_tbq_gvec.o \
qg_tbq_v.o \
qg_tbq_z.o \
inter.o \
inter_gg.o \
inter_qq.o \
reals1.o \
reals2.o \
reals3.o \
reals_qq.o \
stopf0.o \
stopf1.o \
stop_def.o \
Ammm.o \
Ammp.o \
Ampm.o \
Ampp.o \
Apmm.o \
Apmp.o \
Appm.o \
Appp.o \
Bmmm.o \
Bmmp.o \
Bmpm.o \
Bmpp.o \
Bpmm.o \
Bpmp.o \
Bppm.o \
Bppp.o

STOPJETFILES = \
inters.o \
inters_qq.o \
qq_tbg.o \
qq_tbg_g.o \
qq_tbg_gs.o \
qq_tbg_gvec.o \
qq_tbg_v.o \
qq_tbg_z.o

EPEM3JFILES = \
epem3j.o \
epem3j_g.o \
epem3j_gs.o \
epem3j_gvec.o \
epem3j_v.o \
epem3j_z.o

HWWJETFILES = \
gg_hgagag.o \
gg_hgagagg.o \
gg_hgagag_gs.o \
gg_hgagag_gvec.o \
gg_hgagag_v.o \
gg_hZZg.o \
gg_hZZgg.o \
gg_hZZg_gs.o \
gg_hZZg_gvec.o \
gg_hZZg_v.o \
gg_hZZg_z.o \
gg_hWWg.o \
gg_hWWgg.o \
gg_hWWg_gs.o \
gg_hWWg_gvec.o \
gg_hWWg_v.o \
gg_hWWg_z.o \
gg_hWWgg_v.o \
gg_hWWgg_gvec.o \
gg_hWWggg.o \
gg_hWWgg_gs.o \
gg_hWWgg_z.o \
gg_hZZgg_v.o \
gg_hZZgg_gvec.o \
gg_hZZggg.o \
gg_hZZgg_gs.o \
gg_hZZgg_z.o \

BBHIGGSFILES = \
bbaqh.o \
bbbbh.o \
bbggh.o \
bbghvirt.o \
qqb_H_gvec.o \
qqb_Hg.o \
qqb_Hg_g.o \
qqb_Hg_gs.o \
qqb_Hg_v.o \
qqb_Hg_z.o

DIRGAMFILES = \
qqb_hflgam.o \
qqb_hflgam_g.o \
qqb_hflgam_gs.o \
qqb_hflgam_gvec.o \
qqb_hflgam_v.o \
qqb_hflgam_z.o \
qqb_2jnogg.o \
qqb_dirgam_swap.o \
qqb_dirgam.o \
qqb_twojet.o\
qqb_2jet.o\
qqb_2jet_swap.o\
qqb_dirgam_frag_comb.o\
qqb_dirgam_g.o \
qqb_dirgam_gs.o \
qqb_dirgam_gvec.o \
qqb_dirgam_v.o \
qqb_dirgam_z.o \
small.o \
Bigagam.o \
Bigbgam.o \
Bigcgam.o \
gagaAGTYassemble.o \
gagaAGTYassemble1loop.o \
gagaOneloopsq.o \
AGTYassemble_new.o \
AGTYtest.o \
AGTYAs.o \
AGTYAu.o \
AGTYBs.o \
AGTYBu.o \
AGTYCs.o \
AGTYCu.o \
AGTYD1s.o \
AGTYD1u.o \
AGTYD2s.o \
AGTYD2u.o \
AGTYE1s.o \
AGTYE1u.o \
AGTYE2s.o \
AGTYE2u.o \
AGTYE3s.o \
AGTYE3u.o \
AGTYF1s.o \
AGTYF1u.o \
AGTYG1s.o \
AGTYG1u.o \
AGTYG2s.o \
AGTYG2u.o \
AGTYG3s.o \
AGTYG3u.o \
AGTYX1s.o \
AGTYX1u.o \
AGTYX2s.o \
AGTYX2u.o \
AGTYX3s.o \
AGTYX3u.o \
Finite0x2qqgagas.o \
Finite0x2qqgagau.o \
Finite1x1qqgagas.o \
Finite1x1qqgagau.o \

DMFILES=\
gg_dm_monojet_top.o\
ehsv_dm.o\
qqb_dm_monojet_v_PSamps.o\
qqb_dm_monojet_nf_ax.o\
qqb_dm_monophot_g_Samps.o\
qqb_dm_monophot_g_PSamps.o\
qqb_dm_gg_PSamps.o\
qqb_dm_monophot_v_Samps.o\
qqb_dm_monophot_v_PSamps.o\
dmmonojn_scal.o\
dmmonojn_pscal.o\
qqb_dm_qqb_Samps.o\
qqb_dm_gg_Samps.o\
qqb_dm_monojet_Samps.o\
qqb_dm_monojet_slc_Samps.o\
qqb_dm_monojet_lc_Samps.o\
qqb_dm_monojet_slc_PSamps.o\
qqb_dm_monojet_lc_PSamps.o\
qqb_dm_monojet_v_Samps.o\
qqb_dm_monojet_PSamps.o\
scalar_dm.o\
qqb_dm_monojet.o\
qqb_dmprod.o\
qqb_dm_monojet_g.o\
qqb_dm_monophot_v.o\
qqb_dm_monophot_z.o\
qqb_dm_monophot_fragdips.o\
qqb_dm_monophot_frag.o\
qqb_dm_monojet_Vamps.o\
qqb_dm_qqb.o\
qqb_dm_gg_Vamps.o\
qqb_dm_monophot.o\
qqb_dm_monojet_v_Vamps.o\
qqb_dm_monophot_v_Vamps.o\
qqb_dm_monojet_v.o\
qqb_dm_monojet_z.o\
qqb_dm_monojet_lc_Vamps.o\
qqb_dm_monojet_gs.o\
qqb_dm_monojet_gvec.o\
qqb_dm_monojet_slc_Vamps.o\
qqb_dm_monophot_g.o\
qqb_dm_monophot_gs.o\
qqb_dm_monophot_g_Vamps.o\
gg_dm_monojet.o\
gg_dm_monojet_g.o\
gg_dm_monojet_gs.o\
gg_dm_monojet_gvec.o\
gg_dm_monojet_v.o\
gg_dm_monojet_z.o\
qqb_dm_monojet_AxAmp.o\
qqb_dm_monojet_slc_Axamps.o\
qqb_dm_monojet_lc_Axamps.o\
qqb_dm_monojet_v_Axamps.o\
qqb_dm_monophot_v_Axamps.o\
qqb_dm_monophot_g_Axamps.o\
qqb_dm_gg_Axamps.o\
qqb_dm_qq_Axamps.o\

GAMGAMFILES = \
gggaga_mt.o\
qqb_gamgam.o \
qqb_gamgam_frag.o \
qqb_gamgam_fragdips.o \
qqb_gamgam_g.o \
qqb_gamgam_gs.o \
qqb_gamgam_gvec.o \
qqb_gamgam_v.o \
qqb_gamgam_z.o \
gg_2gam.o \
gg_2gam_g.o \
gg_2gam_gs.o \
gg_2gam_gvec.o \
gg_2gam_v.o \
gg_2gam_z.o \
qqbgnfsq.o \
Cgamgam.o \
A51ppppp.o \
A51mpppp.o \
A51mmppp.o \
A51mpmpp.o \
Aboxfill.o \
ggtogagag.o \
Virtgggamgam.o \
twoloop.o \
oneloopep.o

GGHFILES = \
finitemtcorr.o \
gg_h.o \
gg_hg.o \
gg_h_gs.o \
gg_h_gvec.o \
gg_h_v.o \
gg_h_z.o \
hqqgg.o \
gg_hgamgam.o \
gg_hgamgamg.o \
gg_hgamgam_gs.o \
gg_hgamgam_gvec.o \
gg_hgamgam_v.o \
gg_hgamgam_z.o \
msqgamgam.o \
msqhgamgam.o \
Ftriangle.o

GGHZGAFILES = \
C0DDHK.o \
ffDDHK.o \
hzgamwidth.o \
HZgamMSQ.o \
gg_hzgam.o \
gg_hzgamg.o \
gg_hzgam_gs.o \
gg_hzgam_gvec.o \
gg_hzgam_v.o \
gg_hzgam_z.o 


GGHGFILES = \
amplonumer.o \
Hqarbsq.o \
Hqaqasq.o \
qqgghn.o \
gg_hgg.o \
gg_hgg_v.o \
gg_hgg_z.o \
gg_hgg_gs.o \
gg_hgg_gvec.o \
gg_hg_gs.o \
gg_hg_gvec.o \
gg_hg_v.o \
gg_hg_z.o \
H4plo.o \
HAQggnew.o \
h4gnew.o \
Amplo_AQgg.o \
Ampsq_AQaq.o \
h4q.o \
h4g.o \
hjetfill.o 

GGHGGREALFILES = \
Appppp.o \
nAmpppp.o \
nAmmppp.o \
na2q3g_mmmpp.o \
na2q3g_mmpmp.o \
na2q3g_mpmmp.o \
na2q3g_mpppp.o \
gg_hggg.o \
gg_ggg.o \
amp_h5g.o \
h5g.o \
iperm.o \
h2q3g.o \
h4qg.o \
q4ghppp1.o \
q4ghppp3.o \
q4ghpmp1.o \
q4ghpmp3.o \
extract.o \
gggghn_amp.o \
hqqggdfm.o \
ndveccur.o

GGHGGvirtFILES = \
checkegzres.o \
checkscheme.o \
H4prenorm.o \
Ampvirt_AQgg.o \
Ampvirt_gggg.o \
Ampvirtsq_AQaq.o \
Hggggvsqanal.o \
HAQggvsqanal.o \
Hqarbvsqanal.o \
Hqaqavsqanal.o \
Hqaggvsq.o \
Hggggvsq.o \
Hqarbvsq.o \
Hqaqavsq.o \
GZHggggvsqPoles.o \
GZHqaggvsqPoles.o

H4PCODEFILES = \
A0phiggggpppp.o \
A0phiggggpmmm.o \
A0phiggggmmpp.o \
A0phiggggmpmp.o \
A0phiggggmmmm.o \
A0Hggggpppp.o \
A0Hggggpmmm.o \
A0Hggggmmpp.o \
A0Hggggmpmp.o \
A1phiggggmmmm.o \
A1phiggggmmpp.o \
A1phiggggmpmp.o \
A1Hggggpppp.o \
A1Hggggpmmm.o \
A1Hggggmmpp.o \
A1Hggggmpmp.o \
A0phiAQgg.o \
A0phiAgQg.o \
A0HAQgg.o \
A1phiAQggmpmm.o \
A1phiAQggmpmp.o \
A1phiAQggmppm.o \
A1phiAQggmppp.o \
A1phiAgQgmmpp.o \
A1phiAgQgmppm.o \
A1phiAgQgmppp.o \
A41phiAQgg.o \
A41HAQgg.o \
A43phiAQggmpmm.o \
A43phiAQggmpmp.o \
A43phiAQggmppm.o \
A43HAQgg.o \
A0Hqarb.o \
A0phiqarb.o \
A41Hqarb.o \
A41phiqarb.o \
A42Hqarb.o \
A42phiqarb.o \
Afphiqarbmpmp.o \
Afphiqarbmppm.o \
Alcphiqarbmpmp.o \
Alcphiqarbmppm.o \
Aslcphiqarbmpmp.o \
Aslcphiqarbmppm.o \
BGRL.o \
F31m.o \
F33m.o \
F41m.o \
F41mF.o \
F42me.o \
F42meF.o \
F42mhF.o \
Lsm1DS.o

HHFILES = \
gg_HH.o \
HHamps.o \
PSZHHamps.o \
hgamgamdecay.o

HWWFILES = \
gg_WW_int.o \
qqb_hww.o \
dkqqb_hww_v.o \
dkqqb_hww_g.o \
dkqqb_hww_gs.o \
qqb_hww_g.o \
qqb_hww_gvec.o \
qqb_hww_gs.o \
qqb_hww_tb.o \
qqb_hww_z.o \
qqb_hww_v.o

HZZFILES = \
qqb_hzz.o \
qqb_hzz_g.o \
qqb_hzz_gvec.o \
qqb_hzz_gs.o \
qqb_hzz_z.o \
qqb_hzz_v.o

HTTBARFILES = \
qqb_higgs.o \
qqb_higgs_odd.o \
ehsv.o \
ehsv_odd.o

HJETFILES= \
specialpoint.o \
specialamp.o \
qqbggHamp.o \
qqbggHampsq.o \
DelDBfunctions.o \
DelDTfunctions.o \
real2.o\
integrals.o\
hjetmass_bubints_init.o \
hjetmass_boxints_init.o \
hjetmass_triints_init.o \
hjetmass_bubints_init_dp.o \
hjetmass_boxints_init_dp.o \
hjetmass_triints_init_dp.o \
hjetmass.o \
hjetmass_gvec.o \
dilog.o \
hjetmass_v.o \
hjetmass_highpt_v.o \
hjetmass_hel.o \
hjetmass_z.o \
hjetmass_r.o \
hjetmass_r_amp_qp.o \
hjetmass_r_amp_dp.o \
hjetmass_gs.o \
gggH_new.o \
gqqH_new.o \
spinoru_qp.o \
spinoru_dp.o \
hjetmass_qqbgg_boxints_init_dp.o \
hjetmass_qqbgg_bubints_init_dp.o \
hjetmass_qqbgg_triints_init_dp.o \
hjetmass_qqbgg_boxints_init.o \
hjetmass_qqbgg_bubints_init.o \
hjetmass_qqbgg_triints_init.o \
hjetmass_HqqbQQb_dp.o \
hjetmass_HqqbQQb.o \
hjetmass_qqbQQb_bubints_init_dp.o \
hjetmass_qqbQQb_triints_init_dp.o \
hjetmass_qqbQQb_bubints_init.o \
hjetmass_qqbQQb_triints_init.o \
hjetmass_box_onemass_pppp_dp.o \
hjetmass_box_onemass_pppp.o \
hjetmass_box_pmpm_0_0_0_s123_s12_s23_dp.o \
hjetmass_box_pmpm_0_0_0_s123_s12_s23.o \
hjetmass_box_pmpm_0_0_0_s124_s14_s12_dp.o \
hjetmass_box_pmpm_0_0_0_s124_s14_s12.o \
hjetmass_box_pmpm_0_0_0_s134_s34_s14_dp.o \
hjetmass_box_pmpm_0_0_0_s134_s34_s14.o \
hjetmass_box_pmpm_0_0_0_s234_s23_s34_dp.o \
hjetmass_box_pmpm_0_0_0_s234_s23_s34.o \
hjetmass_box_pmpm_0_0_s12_mhsq_s34_s123_dp.o \
hjetmass_box_pmpm_0_0_s12_mhsq_s34_s123.o \
hjetmass_box_pmpm_0_0_s12_mhsq_s34_s124_dp.o \
hjetmass_box_pmpm_0_0_s12_mhsq_s34_s124.o \
hjetmass_box_pmpm_0_0_s14_mhsq_s23_s124_dp.o \
hjetmass_box_pmpm_0_0_s14_mhsq_s23_s124.o \
hjetmass_box_pmpm_0_0_s14_mhsq_s23_s134_dp.o \
hjetmass_box_pmpm_0_0_s14_mhsq_s23_s134.o \
hjetmass_box_pmpm_0_0_s23_mhsq_s14_s123_dp.o \
hjetmass_box_pmpm_0_0_s23_mhsq_s14_s123.o \
hjetmass_box_pmpm_0_0_s23_mhsq_s14_s234_dp.o \
hjetmass_box_pmpm_0_0_s23_mhsq_s14_s234.o \
hjetmass_box_pmpm_0_0_s34_mhsq_s12_s134_dp.o \
hjetmass_box_pmpm_0_0_s34_mhsq_s12_s134.o \
hjetmass_box_pmpm_0_0_s34_mhsq_s12_s234_dp.o \
hjetmass_box_pmpm_0_0_s34_mhsq_s12_s234.o \
hjetmass_box_pmpm_0_0_s34_mhsq_s12_s234_qp.o \
hjetmass_box_pmpm_0_s12_0_mhsq_s123_s124_d.o \
hjetmass_box_pmpm_0_s12_0_mhsq_s123_s124.o \
hjetmass_box_pmpm_0_s14_0_mhsq_s124_s134_d.o \
hjetmass_box_pmpm_0_s14_0_mhsq_s124_s134.o \
hjetmass_box_pmpm_0_s23_0_mhsq_s123_s234_d.o \
hjetmass_box_pmpm_0_s23_0_mhsq_s123_s234.o \
hjetmass_box_pmpm_0_s34_0_mhsq_s134_s234_d.o \
hjetmass_box_pmpm_0_s34_0_mhsq_s134_s234.o \
hjetmass_box_ppmm_0_0_0_s123_s12_s23_dp.o \
hjetmass_box_ppmm_0_0_0_s123_s12_s23.o \
hjetmass_box_ppmm_0_0_0_s124_s14_s12_dp.o \
hjetmass_box_ppmm_0_0_0_s124_s14_s12.o \
hjetmass_box_ppmm_0_0_0_s134_s34_s14_dp.o \
hjetmass_box_ppmm_0_0_0_s134_s34_s14.o \
hjetmass_box_ppmm_0_0_0_s234_s23_s34_dp.o \
hjetmass_box_ppmm_0_0_0_s234_s23_s34.o \
hjetmass_box_ppmm_0_0_s12_mhsq_s34_s123_dp.o \
hjetmass_box_ppmm_0_0_s12_mhsq_s34_s123.o \
hjetmass_box_ppmm_0_0_s12_mhsq_s34_s124_dp.o \
hjetmass_box_ppmm_0_0_s12_mhsq_s34_s124.o \
hjetmass_box_ppmm_0_0_s14_mhsq_s23_s124_dp.o \
hjetmass_box_ppmm_0_0_s14_mhsq_s23_s124.o \
hjetmass_box_ppmm_0_0_s14_mhsq_s23_s134_dp.o \
hjetmass_box_ppmm_0_0_s14_mhsq_s23_s134.o \
hjetmass_box_ppmm_0_0_s23_mhsq_s14_s123_dp.o \
hjetmass_box_ppmm_0_0_s23_mhsq_s14_s123.o \
hjetmass_box_ppmm_0_0_s23_mhsq_s14_s234_dp.o \
hjetmass_box_ppmm_0_0_s23_mhsq_s14_s234.o \
hjetmass_box_ppmm_0_0_s34_mhsq_s12_s134_dp.o \
hjetmass_box_ppmm_0_0_s34_mhsq_s12_s134.o \
hjetmass_box_ppmm_0_0_s34_mhsq_s12_s234_dp.o \
hjetmass_box_ppmm_0_0_s34_mhsq_s12_s234.o \
hjetmass_box_ppmm_0_s12_0_mhsq_s123_s124_d.o \
hjetmass_box_ppmm_0_s12_0_mhsq_s123_s124.o \
hjetmass_box_ppmm_0_s14_0_mhsq_s124_s134_d.o \
hjetmass_box_ppmm_0_s14_0_mhsq_s124_s134.o \
hjetmass_box_ppmm_0_s23_0_mhsq_s123_s234_d.o \
hjetmass_box_ppmm_0_s23_0_mhsq_s123_s234.o \
hjetmass_box_ppmm_0_s34_0_mhsq_s134_s234_d.o \
hjetmass_box_ppmm_0_s34_0_mhsq_s134_s234.o \
hjetmass_box_pppm_0_0_0_s123_s12_s23_dp.o \
hjetmass_box_pppm_0_0_0_s123_s12_s23.o \
hjetmass_box_pppm_0_0_0_s124_s14_s12_dp.o \
hjetmass_box_pppm_0_0_0_s124_s14_s12.o \
hjetmass_box_pppm_0_0_0_s134_s34_s14_dp.o \
hjetmass_box_pppm_0_0_0_s134_s34_s14.o \
hjetmass_box_pppm_0_0_0_s234_s23_s34_dp.o \
hjetmass_box_pppm_0_0_0_s234_s23_s34.o \
hjetmass_box_pppm_0_0_s12_mhsq_s34_s123_dp.o \
hjetmass_box_pppm_0_0_s12_mhsq_s34_s123.o \
hjetmass_box_pppm_0_0_s12_mhsq_s34_s124_dp.o \
hjetmass_box_pppm_0_0_s12_mhsq_s34_s124.o \
hjetmass_box_pppm_0_0_s14_mhsq_s23_s124_dp.o \
hjetmass_box_pppm_0_0_s14_mhsq_s23_s124.o \
hjetmass_box_pppm_0_0_s14_mhsq_s23_s134_dp.o \
hjetmass_box_pppm_0_0_s14_mhsq_s23_s134.o \
hjetmass_box_pppm_0_0_s23_mhsq_s14_s123_dp.o \
hjetmass_box_pppm_0_0_s23_mhsq_s14_s123.o \
hjetmass_box_pppm_0_0_s23_mhsq_s14_s234_dp.o \
hjetmass_box_pppm_0_0_s23_mhsq_s14_s234.o \
hjetmass_box_pppm_0_0_s34_mhsq_s12_s134_dp.o \
hjetmass_box_pppm_0_0_s34_mhsq_s12_s134.o \
hjetmass_box_pppm_0_0_s34_mhsq_s12_s234_dp.o \
hjetmass_box_pppm_0_0_s34_mhsq_s12_s234.o \
hjetmass_box_pppm_0_s12_0_mhsq_s123_s124_d.o \
hjetmass_box_pppm_0_s12_0_mhsq_s123_s124.o \
hjetmass_box_pppm_0_s14_0_mhsq_s124_s134_d.o \
hjetmass_box_pppm_0_s14_0_mhsq_s124_s134.o \
hjetmass_box_pppm_0_s23_0_mhsq_s123_s234_d.o \
hjetmass_box_pppm_0_s23_0_mhsq_s123_s234.o \
hjetmass_box_pppm_0_s34_0_mhsq_s134_s234_d.o \
hjetmass_box_pppm_0_s34_0_mhsq_s134_s234.o \
hjetmass_box_pppp_0_0_0_s123_s12_s23_dp.o \
hjetmass_box_pppp_0_0_0_s123_s12_s23.o \
hjetmass_box_pppp_0_0_0_s124_s14_s12_dp.o \
hjetmass_box_pppp_0_0_0_s124_s14_s12.o \
hjetmass_box_pppp_0_0_0_s134_s34_s14_dp.o \
hjetmass_box_pppp_0_0_0_s134_s34_s14.o \
hjetmass_box_pppp_0_0_0_s234_s23_s34_dp.o \
hjetmass_box_pppp_0_0_0_s234_s23_s34.o \
hjetmass_box_pppp_0_0_s12_mhsq_s34_s123_dp.o \
hjetmass_box_pppp_0_0_s12_mhsq_s34_s123.o \
hjetmass_box_pppp_0_0_s12_mhsq_s34_s124_dp.o \
hjetmass_box_pppp_0_0_s12_mhsq_s34_s124.o \
hjetmass_box_pppp_0_0_s14_mhsq_s23_s124_dp.o \
hjetmass_box_pppp_0_0_s14_mhsq_s23_s124.o \
hjetmass_box_pppp_0_0_s14_mhsq_s23_s134_dp.o \
hjetmass_box_pppp_0_0_s14_mhsq_s23_s134.o \
hjetmass_box_pppp_0_0_s23_mhsq_s14_s123_dp.o \
hjetmass_box_pppp_0_0_s23_mhsq_s14_s123.o \
hjetmass_box_pppp_0_0_s23_mhsq_s14_s234_dp.o \
hjetmass_box_pppp_0_0_s23_mhsq_s14_s234.o \
hjetmass_box_pppp_0_0_s34_mhsq_s12_s134_dp.o \
hjetmass_box_pppp_0_0_s34_mhsq_s12_s134.o \
hjetmass_box_pppp_0_0_s34_mhsq_s12_s234_dp.o \
hjetmass_box_pppp_0_0_s34_mhsq_s12_s234.o \
hjetmass_box_pppp_0_s12_0_mhsq_s123_s124_d.o \
hjetmass_box_pppp_0_s12_0_mhsq_s123_s124.o \
hjetmass_box_pppp_0_s14_0_mhsq_s124_s134_d.o \
hjetmass_box_pppp_0_s14_0_mhsq_s124_s134.o \
hjetmass_box_pppp_0_s23_0_mhsq_s123_s234_d.o \
hjetmass_box_pppp_0_s23_0_mhsq_s123_s234.o \
hjetmass_box_pppp_0_s34_0_mhsq_s134_s234_d.o \
hjetmass_box_pppp_0_s34_0_mhsq_s134_s234.o \
hjetmass_box_twomass_easy_pppp_dp.o \
hjetmass_box_twomass_easy_pppp.o \
hjetmass_box_twomass_hard_pppp_dp.o \
hjetmass_box_twomass_hard_pppp.o \
hjetmass_bubble_pmpm_mhsq_dp.o \
hjetmass_bubble_pmpm_mhsq.o \
hjetmass_bubble_pmpm_s123_dp.o \
hjetmass_bubble_pmpm_s123.o \
hjetmass_bubble_pmpm_s124_dp.o \
hjetmass_bubble_pmpm_s124.o \
hjetmass_bubble_pmpm_s12_dp.o \
hjetmass_bubble_pmpm_s12.o \
hjetmass_bubble_pmpm_s134_dp.o \
hjetmass_bubble_pmpm_s134.o \
hjetmass_bubble_pmpm_s14_dp.o \
hjetmass_bubble_pmpm_s14.o \
hjetmass_bubble_pmpm_s234_dp.o \
hjetmass_bubble_pmpm_s234.o \
hjetmass_bubble_pmpm_s23_dp.o \
hjetmass_bubble_pmpm_s23.o \
hjetmass_bubble_pmpm_s34_dp.o \
hjetmass_bubble_pmpm_s34.o \
hjetmass_bubble_ppmm_mhsq_dp.o \
hjetmass_bubble_ppmm_mhsq.o \
hjetmass_bubble_ppmm_s123_dp.o \
hjetmass_bubble_ppmm_s123.o \
hjetmass_bubble_ppmm_s124_dp.o \
hjetmass_bubble_ppmm_s124.o \
hjetmass_bubble_ppmm_s134_dp.o \
hjetmass_bubble_ppmm_s134.o \
hjetmass_bubble_ppmm_s14_dp.o \
hjetmass_bubble_ppmm_s14.o \
hjetmass_bubble_ppmm_s234_dp.o \
hjetmass_bubble_ppmm_s234.o \
hjetmass_bubble_ppmm_s23_dp.o \
hjetmass_bubble_ppmm_s23.o \
hjetmass_bubble_pppm_mhsq_dp.o \
hjetmass_bubble_pppm_mhsq.o \
hjetmass_bubble_pppm_s124_dp.o \
hjetmass_bubble_pppm_s124.o \
hjetmass_bubble_pppm_s134_dp.o \
hjetmass_bubble_pppm_s134.o \
hjetmass_bubble_pppm_s14_dp.o \
hjetmass_bubble_pppm_s14.o \
hjetmass_bubble_pppm_s234_dp.o \
hjetmass_bubble_pppm_s234.o \
hjetmass_bubble_pppm_s34_dp.o \
hjetmass_bubble_pppm_s34.o \
hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s123.o \
hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s124.o \
hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s.o \
hjetmass_qqbgg_box_pmmp_0_s12_0_mhsq_s123_.o \
hjetmass_qqbgg_box_pmmp_0_s12_0_mhsq_s123_s124.o \
hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s123.o \
hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s124.o \
hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s.o \
hjetmass_qqbgg_box_pmpm_0_s12_0_mhsq_s123_.o \
hjetmass_qqbgg_box_pmpm_0_s12_0_mhsq_s123_s124.o \
hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s123.o \
hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s124.o \
hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s.o \
hjetmass_qqbgg_box_pmpp_0_s12_0_mhsq_s123_.o \
hjetmass_qqbgg_box_pmpp_0_s12_0_mhsq_s123_s124.o \
hjetmass_qqbgg_bubble_pmmp_mhsq_dp.o \
hjetmass_qqbgg_bubble_pmmp_mhsq.o \
hjetmass_qqbgg_bubble_pmmp_s123_dp.o \
hjetmass_qqbgg_bubble_pmmp_s123.o \
hjetmass_qqbgg_bubble_pmmp_s124_dp.o \
hjetmass_qqbgg_bubble_pmmp_s124.o \
hjetmass_qqbgg_bubble_pmmp_s12_dp.o \
hjetmass_qqbgg_bubble_pmmp_s12.o \
hjetmass_qqbgg_bubble_pmmp_s34_dp.o \
hjetmass_qqbgg_bubble_pmmp_s34.o \
hjetmass_qqbgg_bubble_pmpm_mhsq_dp.o \
hjetmass_qqbgg_bubble_pmpm_mhsq.o \
hjetmass_qqbgg_bubble_pmpm_s123_dp.o \
hjetmass_qqbgg_bubble_pmpm_s123.o \
hjetmass_qqbgg_bubble_pmpm_s124_dp.o \
hjetmass_qqbgg_bubble_pmpm_s124.o \
hjetmass_qqbgg_bubble_pmpm_s12_dp.o \
hjetmass_qqbgg_bubble_pmpm_s12.o \
hjetmass_qqbgg_bubble_pmpm_s34_dp.o \
hjetmass_qqbgg_bubble_pmpm_s34.o \
hjetmass_qqbgg_bubble_pmpp_mhsq_dp.o \
hjetmass_qqbgg_bubble_pmpp_mhsq.o \
hjetmass_qqbgg_bubble_pmpp_s123_dp.o \
hjetmass_qqbgg_bubble_pmpp_s123.o \
hjetmass_qqbgg_bubble_pmpp_s124_dp.o \
hjetmass_qqbgg_bubble_pmpp_s124.o \
hjetmass_qqbgg_bubble_pmpp_s12_dp.o \
hjetmass_qqbgg_bubble_pmpp_s12.o \
hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_d.o \
hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0.o \
hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_rat.o \
hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_r.o \
hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_d.o \
hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0.o \
hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat.o \
hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_r.o \
hjetmass_qqbgg_triangle_pmmp_s12_s123_0_dp.o \
hjetmass_qqbgg_triangle_pmmp_s12_s123_0.o \
hjetmass_qqbgg_triangle_pmmp_s12_s124_0_dp.o \
hjetmass_qqbgg_triangle_pmmp_s12_s124_0.o \
hjetmass_qqbgg_triangle_pmmp_s34_0_0_dp.o \
hjetmass_qqbgg_triangle_pmmp_s34_0_0.o \
hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_.o \
hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12.o \
hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_rat.o \
hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_d.o \
hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0.o \
hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_rat.o \
hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_r.o \
hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_d.o \
hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0.o \
hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_rat.o \
hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_r.o \
hjetmass_qqbgg_triangle_pmpm_s12_s123_0_dp.o \
hjetmass_qqbgg_triangle_pmpm_s12_s123_0.o \
hjetmass_qqbgg_triangle_pmpm_s12_s124_0_dp.o \
hjetmass_qqbgg_triangle_pmpm_s12_s124_0.o \
hjetmass_qqbgg_triangle_pmpm_s34_0_0_dp.o \
hjetmass_qqbgg_triangle_pmpm_s34_0_0.o \
hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_.o \
hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12.o \
hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_rat.o \
hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_d.o \
hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0.o \
hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_rat.o \
hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_r.o \
hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_d.o \
hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0.o \
hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_rat.o \
hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_r.o \
hjetmass_qqbgg_triangle_pmpp_s12_s123_0_dp.o \
hjetmass_qqbgg_triangle_pmpp_s12_s123_0.o \
hjetmass_qqbgg_triangle_pmpp_s12_s124_0_dp.o \
hjetmass_qqbgg_triangle_pmpp_s12_s124_0.o \
hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_.o \
hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12.o \
hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_rat.o \
hjetmass_triangle_pmpm_s12_0_0_dp.o \
hjetmass_triangle_pmpm_s12_0_0.o \
hjetmass_triangle_pmpm_s123_0_mhsq_dp.o \
hjetmass_triangle_pmpm_s123_0_mhsq.o \
hjetmass_triangle_pmpm_s123_0_mhsq_rat_dp.o \
hjetmass_triangle_pmpm_s123_0_mhsq_rat.o \
hjetmass_triangle_pmpm_s124_0_mhsq_dp.o \
hjetmass_triangle_pmpm_s124_0_mhsq.o \
hjetmass_triangle_pmpm_s124_0_mhsq_rat_dp.o \
hjetmass_triangle_pmpm_s124_0_mhsq_rat.o \
hjetmass_triangle_pmpm_s12_s123_0_dp.o \
hjetmass_triangle_pmpm_s12_s123_0.o \
hjetmass_triangle_pmpm_s12_s123_0_rat_dp.o \
hjetmass_triangle_pmpm_s12_s123_0_rat.o \
hjetmass_triangle_pmpm_s12_s124_0_dp.o \
hjetmass_triangle_pmpm_s12_s124_0.o \
hjetmass_triangle_pmpm_s12_s124_0_rat_dp.o \
hjetmass_triangle_pmpm_s12_s124_0_rat.o \
hjetmass_triangle_pmpm_s134_0_mhsq_dp.o \
hjetmass_triangle_pmpm_s134_0_mhsq.o \
hjetmass_triangle_pmpm_s134_0_mhsq_rat_dp.o \
hjetmass_triangle_pmpm_s134_0_mhsq_rat.o \
hjetmass_triangle_pmpm_s14_0_0_dp.o \
hjetmass_triangle_pmpm_s14_0_0.o \
hjetmass_triangle_pmpm_s14_s124_0_dp.o \
hjetmass_triangle_pmpm_s14_s124_0.o \
hjetmass_triangle_pmpm_s14_s124_0_rat_dp.o \
hjetmass_triangle_pmpm_s14_s124_0_rat.o \
hjetmass_triangle_pmpm_s14_s134_0_dp.o \
hjetmass_triangle_pmpm_s14_s134_0.o \
hjetmass_triangle_pmpm_s14_s134_0_rat_dp.o \
hjetmass_triangle_pmpm_s14_s134_0_rat.o \
hjetmass_triangle_pmpm_s23_0_0_dp.o \
hjetmass_triangle_pmpm_s23_0_0.o \
hjetmass_triangle_pmpm_s234_0_mhsq_dp.o \
hjetmass_triangle_pmpm_s234_0_mhsq.o \
hjetmass_triangle_pmpm_s234_0_mhsq_rat_dp.o \
hjetmass_triangle_pmpm_s234_0_mhsq_rat.o \
hjetmass_triangle_pmpm_s23_mhsq_s14_dp.o \
hjetmass_triangle_pmpm_s23_mhsq_s14.o \
hjetmass_triangle_pmpm_s23_mhsq_s14_rat_dp.o \
hjetmass_triangle_pmpm_s23_mhsq_s14_rat.o \
hjetmass_triangle_pmpm_s23_s123_0_dp.o \
hjetmass_triangle_pmpm_s23_s123_0.o \
hjetmass_triangle_pmpm_s23_s123_0_rat_dp.o \
hjetmass_triangle_pmpm_s23_s123_0_rat.o \
hjetmass_triangle_pmpm_s23_s234_0_dp.o \
hjetmass_triangle_pmpm_s23_s234_0.o \
hjetmass_triangle_pmpm_s23_s234_0_rat_dp.o \
hjetmass_triangle_pmpm_s23_s234_0_rat.o \
hjetmass_triangle_pmpm_s34_0_0_dp.o \
hjetmass_triangle_pmpm_s34_0_0.o \
hjetmass_triangle_pmpm_s34_mhsq_s12_dp.o \
hjetmass_triangle_pmpm_s34_mhsq_s12.o \
hjetmass_triangle_pmpm_s34_mhsq_s12_rat_dp.o \
hjetmass_triangle_pmpm_s34_mhsq_s12_rat.o \
hjetmass_triangle_pmpm_s34_s134_0_dp.o \
hjetmass_triangle_pmpm_s34_s134_0.o \
hjetmass_triangle_pmpm_s34_s134_0_rat_dp.o \
hjetmass_triangle_pmpm_s34_s134_0_rat.o \
hjetmass_triangle_pmpm_s34_s234_0_dp.o \
hjetmass_triangle_pmpm_s34_s234_0.o \
hjetmass_triangle_pmpm_s34_s234_0_rat_dp.o \
hjetmass_triangle_pmpm_s34_s234_0_rat.o \
hjetmass_triangle_ppmm_s123_0_mhsq_dp.o \
hjetmass_triangle_ppmm_s123_0_mhsq.o \
hjetmass_triangle_ppmm_s123_0_mhsq_rat_dp.o \
hjetmass_triangle_ppmm_s123_0_mhsq_rat.o \
hjetmass_triangle_ppmm_s124_0_mhsq_dp.o \
hjetmass_triangle_ppmm_s124_0_mhsq.o \
hjetmass_triangle_ppmm_s124_0_mhsq_rat_dp.o \
hjetmass_triangle_ppmm_s124_0_mhsq_rat.o \
hjetmass_triangle_ppmm_s134_0_mhsq_dp.o \
hjetmass_triangle_ppmm_s134_0_mhsq.o \
hjetmass_triangle_ppmm_s134_0_mhsq_rat_dp.o \
hjetmass_triangle_ppmm_s134_0_mhsq_rat.o \
hjetmass_triangle_ppmm_s14_0_0_dp.o \
hjetmass_triangle_ppmm_s14_0_0.o \
hjetmass_triangle_ppmm_s14_s124_0_dp.o \
hjetmass_triangle_ppmm_s14_s124_0.o \
hjetmass_triangle_ppmm_s14_s124_0_rat_dp.o \
hjetmass_triangle_ppmm_s14_s124_0_rat.o \
hjetmass_triangle_ppmm_s14_s134_0_dp.o \
hjetmass_triangle_ppmm_s14_s134_0.o \
hjetmass_triangle_ppmm_s14_s134_0_rat_dp.o \
hjetmass_triangle_ppmm_s14_s134_0_rat.o \
hjetmass_triangle_ppmm_s23_0_0_dp.o \
hjetmass_triangle_ppmm_s23_0_0.o \
hjetmass_triangle_ppmm_s234_0_mhsq_dp.o \
hjetmass_triangle_ppmm_s234_0_mhsq.o \
hjetmass_triangle_ppmm_s234_0_mhsq_rat_dp.o \
hjetmass_triangle_ppmm_s234_0_mhsq_rat.o \
hjetmass_triangle_ppmm_s23_mhsq_s14_dp.o \
hjetmass_triangle_ppmm_s23_mhsq_s14.o \
hjetmass_triangle_ppmm_s23_mhsq_s14_rat_dp.o \
hjetmass_triangle_ppmm_s23_mhsq_s14_rat.o \
hjetmass_triangle_ppmm_s23_s123_0_dp.o \
hjetmass_triangle_ppmm_s23_s123_0.o \
hjetmass_triangle_ppmm_s23_s123_0_rat_dp.o \
hjetmass_triangle_ppmm_s23_s123_0_rat.o \
hjetmass_triangle_ppmm_s23_s234_0_dp.o \
hjetmass_triangle_ppmm_s23_s234_0.o \
hjetmass_triangle_ppmm_s23_s234_0_rat_dp.o \
hjetmass_triangle_ppmm_s23_s234_0_rat.o \
hjetmass_triangle_pppm_s123_0_mhsq_dp.o \
hjetmass_triangle_pppm_s123_0_mhsq.o \
hjetmass_triangle_pppm_s123_0_mhsq_rat_dp.o \
hjetmass_triangle_pppm_s123_0_mhsq_rat.o \
hjetmass_triangle_pppm_s124_0_mhsq_dp.o \
hjetmass_triangle_pppm_s124_0_mhsq.o \
hjetmass_triangle_pppm_s124_0_mhsq_rat_dp.o \
hjetmass_triangle_pppm_s124_0_mhsq_rat.o \
hjetmass_triangle_pppm_s134_0_mhsq_dp.o \
hjetmass_triangle_pppm_s134_0_mhsq.o \
hjetmass_triangle_pppm_s134_0_mhsq_rat_dp.o \
hjetmass_triangle_pppm_s134_0_mhsq_rat.o \
hjetmass_triangle_pppm_s14_0_0_dp.o \
hjetmass_triangle_pppm_s14_0_0.o \
hjetmass_triangle_pppm_s14_s124_0_dp.o \
hjetmass_triangle_pppm_s14_s124_0.o \
hjetmass_triangle_pppm_s14_s124_0_rat_dp.o \
hjetmass_triangle_pppm_s14_s124_0_rat.o \
hjetmass_triangle_pppm_s14_s134_0_dp.o \
hjetmass_triangle_pppm_s14_s134_0.o \
hjetmass_triangle_pppm_s14_s134_0_rat_dp.o \
hjetmass_triangle_pppm_s14_s134_0_rat.o \
hjetmass_triangle_pppm_s234_0_mhsq_dp.o \
hjetmass_triangle_pppm_s234_0_mhsq.o \
hjetmass_triangle_pppm_s234_0_mhsq_rat_dp.o \
hjetmass_triangle_pppm_s234_0_mhsq_rat.o \
hjetmass_triangle_pppm_s23_mhsq_s14_dp.o \
hjetmass_triangle_pppm_s23_mhsq_s14.o \
hjetmass_triangle_pppm_s23_mhsq_s14_rat_dp.o \
hjetmass_triangle_pppm_s23_mhsq_s14_rat.o \
hjetmass_triangle_pppm_s34_0_0_dp.o \
hjetmass_triangle_pppm_s34_0_0.o \
hjetmass_triangle_pppm_s34_mhsq_s12_dp.o \
hjetmass_triangle_pppm_s34_mhsq_s12.o \
hjetmass_triangle_pppm_s34_mhsq_s12_rat_dp.o \
hjetmass_triangle_pppm_s34_mhsq_s12_rat.o \
hjetmass_triangle_pppm_s34_s134_0_dp.o \
hjetmass_triangle_pppm_s34_s134_0.o \
hjetmass_triangle_pppm_s34_s134_0_rat_dp.o \
hjetmass_triangle_pppm_s34_s134_0_rat.o \
hjetmass_triangle_pppm_s34_s234_0_dp.o \
hjetmass_triangle_pppm_s34_s234_0.o \
hjetmass_triangle_pppm_s34_s234_0_rat_dp.o \
hjetmass_triangle_pppm_s34_s234_0_rat.o \
hjetmass_triangle_pppp_s123_0_mhsq_dp.o \
hjetmass_triangle_pppp_s123_0_mhsq.o \
hjetmass_triangle_pppp_s123_0_mhsq_rat_dp.o \
hjetmass_triangle_pppp_s123_0_mhsq_rat.o \
hjetmass_triangle_pppp_s124_0_mhsq_dp.o \
hjetmass_triangle_pppp_s124_0_mhsq.o \
hjetmass_triangle_pppp_s124_0_mhsq_rat_dp.o \
hjetmass_triangle_pppp_s124_0_mhsq_rat.o \
hjetmass_triangle_pppp_s134_0_mhsq_dp.o \
hjetmass_triangle_pppp_s134_0_mhsq.o \
hjetmass_triangle_pppp_s134_0_mhsq_rat_dp.o \
hjetmass_triangle_pppp_s134_0_mhsq_rat.o \
hjetmass_triangle_pppp_s234_0_mhsq_dp.o \
hjetmass_triangle_pppp_s234_0_mhsq.o \
hjetmass_triangle_pppp_s234_0_mhsq_rat_dp.o \
hjetmass_triangle_pppp_s234_0_mhsq_rat.o \
hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s123_dp.o \
hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_dp.o \
hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s123_dp.o \
hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_dp.o \
hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s123_dp.o \
hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_dp.o


INTEGRATEFILES = \
dgauss.o \
ebook.o \
mbook.o \
mtxt.o \
cxx11_random.cpo \
cxx11_random_mod.o \
updatehisto.o

LIBFILES = \
cli2.o \
ddilog.o \
lenocc.o \
Li2.o \
Li3.o \
Li4.o \
WGPLG.o \
cln.o

SPINORFILES = \
fillgam.o \
UbKlSt.o \
VKlSt.o \
Ubkslash.o \
kslashU.o \
uspinor0.o \
ubarspinor0.o \
pol_massless.o \
pol_real.o \
cdot.o

NEEDFILES = \
usercode_f77.o \
olo_dummy.o \
lhcbcode.o \
set_anomcoup.o \
debugtools_m.o \
maketaucut.o \
maketaucut_bb.o \
mcfm_writelhe.o \
write_gg_lhe.o \
hbbdecay.o \
htautaudecay.o \
hwwdecay.o \
hbbdecay_g.o \
hbbdecay_v.o \
arraysort.o \
aveptjet.o \
banner.o \
boost.o \
boostx.o \
branch.o \
cdotpr.o \
checkndotp.o \
checkversion.o \
checkjets.o \
ckmfill.o \
computepdfuncertainty.o \
count_parts.o \
coupling.o \
coupling2.o \
couplz.o \
dclaus.o \
determinefilenames.o \
dipoles.o \
dipoles_fac.o \
dipoles_mass.o \
dipolesub.o \
dipolesubx.o \
dipolesubxx.o \
dipolesubx_new.o \
donothing_gvec.o \
dot.o \
dotem.o \
etmiss.o \
ff_alpha.o \
getbs.o \
getptilde.o \
getptildejet.o \
getptQ1.o \
higgsp.o \
higgsw.o \
histofin.o \
itransform.o \
includedipole.o \
interpolate_hto.o \
is_functions.o \
kcasestring.o \
kpartstring.o \
masscuts.o \
mcfmmain.o \
mcfmsub.o \
mcfm_exit.o \
mcfm_init.o \
mcfm_vegasnr.o \
mcfm_vegas_nnlo.o \
mcfm_vegas_adaptive.o \
vegasnr.o \
massfrun.o \
ptyrap.o \
r.o \
read_jetcuts.o \
read_dm_params.o\
readcoup.o \
reader_input.o \
reset_aem.o \
relativeewhisto.o \
Rgen.o \
scaleset.o \
scaleset_m34.o \
scaleset_m345.o \
scaleset_m3456.o \
scaleset_Msqpt34sq.o \
scaleset_Msqpt345sq.o \
scaleset_Msqpt5sq.o \
scaleset_Msqptj1sq.o \
scaleset_Msqsumptjsq.o \
scaleset_m34sqsumptjsq.o \
scaleset_ptj1.o \
scaleset_ptphoton.o \
scaleset_HT.o \
scaleset_ddis.o \
scaleset_mVpmH.o\
sethparams.o \
setmb_msbar.o \
setrunname.o \
setuptau.o \
setvdecay.o \
smalls.o \
smalltau.o \
smarthistos.o \
spinork.o \
spinoru.o \
spinorz.o \
storedip.o \
storedip_mass.o \
storeptilde.o \
store_zdip.o \
swapjet.o \
transform.o \
transform_mass.o \
writeinfo.o \
writeinput.o \
writeout.o \
writereference.o \
dips_mass.o \
mathfun.o

PARTONFILES = \
checkpath.o \
newton1.o

PHASEFILES = \
phase_realstop.o \
vetow_2gam.o \
breitw.o \
breitw_mod.o \
breitw1.o \
gen2.o \
gen2a.o \
gen2jet.o \
gen2m.o \
gen3.o \
gen3b.o \
gen3h.o \
gen_higgszgam.o \
gen_vgamj.o \
gen3jet.o \
gen3jetgaga.o \
gen3m.o \
gen3m_rap.o \
gen3mdk.o \
gen4.o \
gen4_intf.o \
gen4_3m.o \
gen4h.o \
gen4handc.o \
gen4m.o\
gen4mdk.o \
gen4mdkrad.o \
gen_HZgamj.o \
gen5.o \
gen5h.o \
phase5h.o \
phase5Vgam.o \
gen5mdk.o \
gen6.o \
gen6_rap.o \
gen7.o \
gen7m.o \
phase7m_alt.o \
gen7_rap.o \
gen8.o \
gen8_rap.o \
gen9_rap.o \
gen9dk_rap.o \
gen10.o \
genpt.o \
phase10.o \
genttvdk.o \
gen_photons_jets.o \
gen_njets.o \
gen_soft.o \
gen_Vphotons_jets.o \
gen_Vphotons_jets_dkrad.o \
gen_Vphotons_jets_dkrad2.o \
genVphoton.o \
gen_stop.o \
pstop.o \
phase2.o \
phase3.o \
phase3m.o \
phase4.o \
phase41.o \
phase4m.o \
phase4Vgam.o \
phase5.o \
phase5a.o \
phase51.o \
phase6.o \
phase6a.o \
phase6b.o \
phase7.o \
phase7a.o \
phase7b.o \
phase7dk.o \
phase7m.o \
phase8.o \
phi1_2.o \
phi1_2bis.o \
phi1_2bw.o \
phi1_2nobw.o \
phi1_2m.o \
phi1_2m_bw.o \
phi1_2m_nobw.o \
phi3m.o \
phi3m0.o \
wt2gen.o \
wt4gen.o \
wtgen.o

PROCDEPFILES = \
checkorder.o \
chooser.o \
fragint.o \
gen_lops.o \
gen_realps.o \
lowint.o \
realint.o \
virtint.o \
scetint.o 

QQHFILES = \
qq_Hqq_g.o \
VV_Hqq.o \
VV_Hqq_g.o \
VV_Hqq_gs.o \
VV_Hqq_v.o \
VV_Hqq_z.o \
WW_Hqq.o \
WW_Hqq_g.o \
WW_Hqq_gs.o \
ZZ_Hqq.o \
ZZ_Hqq_g.o \
ZZ_Hqq_gs.o

QQHGAGAHFILES = \
VV_Hgaga.o \
VV_Hgaga_g.o \
VV_Hgaga_gs.o \
VV_Hgaga_v.o \
VV_Hgaga_z.o \
WW_Hgaga.o \
WW_Hgaga_g.o \
WW_Hgaga_gs.o \
ZZ_Hgaga.o \
ZZ_Hgaga_g.o \
ZZ_Hgaga_gs.o

QQHWWFILES = \
VV_HWW.o \
WW_HWW.o \
ZZ_HWW.o \
VV_HWW_g.o \
WW_HWW_g.o \
ZZ_HWW_g.o \
VV_HWW_gs.o \
WW_HWW_gs.o \
ZZ_HWW_gs.o \
VV_HWW_v.o \
VV_HWW_z.o

QQHZZFILES = \
VV_HZZ.o \
WW_HZZ.o \
ZZ_HZZ.o \
VV_HZZ_g.o \
WW_HZZ_g.o \
ZZ_HZZ_g.o \
VV_HZZ_gs.o \
WW_HZZ_gs.o \
ZZ_HZZ_gs.o \
VV_HZZ_v.o \
VV_HZZ_z.o

QQZTTFILES = \
qqbZtt.o \
qqbZtt1.o

HELASFILES = \
oxxxxx.o \
ixxxxx.o \
jioxxx.o \
helas_ext.o

SINGLETOP2FILES = \
singletop2_m.o \
singletop2_realamps_m.o \
singletop2_ints_nonres_m.o \
singletop2_ints_new_m.o \
singletop2_virtamps_m.o \
singletop2_scale_m.o \
singletop2_realamps_nonres_m.o \
singletop2_decaywidth.o \
singletop2_ints_dd_m.o \
singletop2_ints_nonres_dd_m.o

SINGLETOPFILES = \
topwidth.o \
extend_trans.o \
qqb_tbb.o \
qqb_tbb_g.o \
qqb_tbb_v.o \
qqb_tbb_z.o \
qqb_tbb_gs.o \
qqb_tbb_gdk.o \
qqb_tbb_gsdk.o \
qqb_tbb_vdk.o \
coefsdk.o \
bq_tpq.o \
bq_tpq_v.o \
bq_tpq_z.o \
bq_tpq_gdk.o \
bq_tpq_gsdk.o \
bq_tpq_vdk.o \
dkqqb_tbbdk_v.o \
dkqqb_tbbdk_g.o \
dkqqb_tbbdk_gs.o \
qqb_tbbdk.o \
qqb_tbbdk_g.o \
qqb_tbbdk_gs.o \
qqb_tbbdk_v.o \
qqb_tbbdk_z.o \
schantoponshell.o \
schantoponshellv.o \
schanatoponshellv.o \
schantoponshellg.o \
schanatoponshellg.o \
schanatoponshell.o \
amp_dkqg_tbqdk_g.o \
dkqg_tbqdk_g.o \
dkqg_tbqdk_gs.o \
dkqg_tbqdk_v.o \
interdk.o \
interdk_gg.o \
interdk_qq.o \
qg_tbqdk.o \
qg_tbqdk_g.o \
qg_tbqdk_gs.o \
qg_tbqdk_gvec.o \
qg_tbqdk_v.o \
qg_tbqdk_z.o \
qg_tbqndk_amp.o \
qg_tbqndk_ampanti.o \
singleatoponshell.o \
singletoponshell.o

SINGLETOPHFILES = \
lowerh.o \
middleh.o \
qq_tchan_htq_dk.o \
qq_tchan_htq_dk_gs.o \
qq_tchan_htq_dk_v.o \
qq_tchan_htq_dk_z.o \
qq_tchan_htq.o \
qq_tchan_htqg_dk.o \
qq_tchan_htqg.o \
qq_tchan_htq_gs.o \
qq_tchan_htq_v.o \
qq_tchan_htq_z.o \
scalarh.o \
stringsh_dk.o \
stringsh.o \
ubhtdgsqdk.o \
ubhtdgsq.o \
ubhtdsqdk.o \
ubhtdsq.o \
ubthdamp_split.o

SINGLETOPZFILES = \
extend_trans_ztj.o \
extradk.o \
extra.o \
kininv.o \
lowerdk_partbox.o \
lowerdk_parttri.o \
lower_partbox.o \
lower_parttri.o \
middledk.o \
middle.o \
qq_tchan_ztq_dk.o \
qq_tchan_ztq_dk_gs.o \
qq_tchan_ztq_dk_v.o \
qq_tchan_ztq_dk_z.o \
qq_tchan_ztq.o \
qq_tchan_ztqg_dk.o \
qq_tchan_ztqg.o \
qq_tchan_ztq_gs.o \
qq_tchan_ztq_v.o \
qq_tchan_ztq_z.o \
scalardk.o \
scalar.o \
strings_dk.o \
strings.o \
ubtzdamp_dk.o \
ubtzdsq.o \
ubztdgsqdk.o \
ubztdgsq.o \
upperdk_partbox.o \
upperdk_parttri.o \
upper_partbox.o \
upper_parttri.o \
vertices_bt2.o \
vertices_tt1.o \
vertices_tt2.o

TAUTAUFILES = \
qqb_tautau.o \
tautauww.o \
dotks.o \
std.o

TDHPLFILES = \
fillcoeff2dhpl210.o \
fillcoeff2dhpl320e.o \
fillcoeff2dhpl321u.o \
fillcoeff2dhpl310.o \
fillcoeff2dhpl321.o \
fillcoeff2dhplaux.o \
fillcoeff2dhpl320.o \
fillcoeff2dhpl321e.o \
tdhpl.o

TOPFILES = \
qqb_QQb.o \
qqb_QQb_g.o \
qqb_QQb_gs.o \
qqb_QQb_gvec.o \
qqb_QQb_v.o \
qqb_QQb_z.o

TOPDKFILES = \
dk1qqb_QQb_g.o \
dk1qqb_QQb_gs.o \
dk2qqb_QQb_g.o \
dk2qqb_QQb_gs.o \
dkqqb_QQb_v.o \
dkW1qqb_QQb_g.o \
dkW1qqb_QQb_gs.o \
dkW2qqb_QQb_g.o \
dkW2qqb_QQb_gs.o \
dkWqqb_QQb_v.o \
ggttww1n.o \
toppaironshell.o \
qqb_QQbdk.o \
ttbqqbg_sq.o \
qqb_QQbdk_g.o \
qqb_QQbdk_gs.o \
qqb_QQbdk_gvec.o \
qqb_QQbdk_v.o \
qqb_QQbdk_z.o \
dkuqqb_QQb_v.o \
dk1uqqb_QQb_g.o \
dk2uqqb_QQb_g.o \
dk1uqqb_QQb_gs.o \
dk2uqqb_QQb_gs.o \
qqb_QQbdku.o \
qqb_QQbdku_g.o \
qqb_QQbdku_gs.o \
qqb_QQbdku_gvec.o \
qqb_QQbdku_v.o \
qqb_QQbdku_z.o \
extend_trans_ttb.o \
ttbgggppp.o \
ttbgggppm.o \
ttbgggpmm.o \
ttbgggmmm.o \
ttbgggmpp.o \
ttbgggmmp.o \
ttbgggmpm.o \
ttbgggpmp.o \
ttbqqbsqpp.o \
ttbqqbsqpm.o \
ttbqqbsqmp.o \
ttbqqbsqmm.o \
ttbqqbtqpp.o \
ttbqqbtqpm.o \
ttbqqbtqmp.o \
ttbqqbtqmm.o \
ttbqqbqqpp.o \
ttbqqbqqpm.o \
ttbqqbqqmp.o \
ttbqqbqqmm.o \
ttbqqbrqpp.o \
ttbqqbrqpm.o \
ttbqqbrqmp.o \
ttbqqbrqmm.o

TOPHFILES= \
ggtth.o \
qqbtth.o \
qqb_tth.o \
qqb_tottth.o \
ttggHamp.o \
ttggHdriver.o \
ttqqHampsq.o

TOPWFILES = \
qqb_ttw.o \
qqb_ttw_v.o \
dk1qqb_ttw_g.o \
dk2qqb_ttw_g.o \
dk1qqb_ttw_gs.o \
dk2qqb_ttw_gs.o \
dkqqb_ttw_v.o \
qqb_ttw_z.o \
qqb_ttw_g.o \
qqb_ttw_gs.o \
extend_trans_ttw.o \
ttWprod.o \
trodmsqm.o \
adecayrod.o \
tdecayrod.o \
tdecayrodg.o \
adecayrodg.o \
adecayrod_v.o \
tdecayrod_v.o 

TOPZFILES = \
qqb_ttz.o \
qqbttz.o \
ggttz.o

USERFILES = \
autoplots.o \
bookplot.o \
bookrelew.o \
cdfhwwcuts.o \
cms_higgsWW.o \
ATLAS_hww.o \
ATLAS_hww2013.o \
ATLAS_ss.o \
CMS_hzz.o \
CMS_hzz_vbf.o \
dm_cuts.o\
deltarj.o \
durhamalg.o \
etdoublebin.o \
METHTdoublebin.o \
fill_stdhep.o \
filterW_bjet.o \
filterWbbmas.o \
frix.o \
FBhisto.o \
genclust2.o \
genclust_kt.o \
genclust_cone.o \
genclustphotons.o \
gencuts.o \
gencuts_input.o \
gencuts_ATLAS_gaga2.o \
gencuts_Zt.o \
gencuts_WZjj.o \
gencuts_VHbb.o\
gencuts_VHWW.o\
genplots.o \
getet.o \
hwwcuts.o \
hwwjetplots.o \
integratehisto.o \
irregbins.o \
iso.o \
jetlabel_to_stdhep.o \
jetreorder.o \
mcfm_getunweighted.o \
mdata.o \
miscclust.o \
nplotter.o \
nplotter_4ftwdk.o \
nplotter_auto.o \
nplotter_dirgam.o \
nplotter_dm_monj.o\
nplotter_dm_mongam.o\
nplotter_gamgam.o \
nplotter_generic.o \
nplotter_gmgmjt.o \
nplotter_trigam.o \
nplotter_tbbar.o \
nplotter_ttbar.o \
nplotter_ttw.o \
nplotter_ttZ.o \
nplotter_tt_tot.o \
nplotter_twojet.o \
nplotter_Vgamma.o \
nplotter_VV.o \
nplotter_VHbbar.o \
nplotter_VHbbarHXSWG.o \
nplotter_W_only.o \
nplotter_Z_only.o \
nplotter_Wbbmas.o \
nplotter_wgamgam.o \
nplotter_Wjets.o \
nplotter_WPWP.o \
nplotter_WWjet.o \
nplotter_zgamgam.o \
nplotter_zgamjet.o \
nplotter_Zbbbar.o \
nplotter_Zjets.o \
nplotter_Ztj.o \
nplotter_Ztjdk.o \
nplotter_ZZlept.o \
nplotter_VHgaga.o\
nplotter_VHWW.o\
nplotter_ktopanom.o \
wrecons.o \
fasthist.o \
idjet.o \
photo_iso.o \
photo_iso_z.o \
photo_iso_phys.o \
photgenclust.o \
photoncuts.o \
setnotag.o \
singletopreconstruct.o \
stopcuts.o \
topreconstruct.o \
VBS.o \
wbfcuts.o \
wbfcuts_jeppe.o \
wconstruct.o

VOLFILES = \
qqb_vol.o \
vol.o \
vol3_mass.o \
vol_mass.o

VVFILES = \
gg_VV.o \
gg_VV_all.o \
gg_hvv_tb.o \
getggHWWamps.o \
getggWWamps.o \
gen4vv.o \
gen4handcvv.o \
gen4hvv.o

WFILES = \
qqb_w.o \
qqb_w_g.o \
ampqqbgll.o \
qqb_w_gs.o \
qqb_w_v.o \
qqb_w_z.o

W1JETFILES = \
A5NLO.o \
A51.o \
A52.o \
A53.o \
qqb_w1jet_gs.o \
qqb_w1jet_v.o \
qqb_w1jet_z.o \
qqb_w_gvec.o \
virt5.o 
 
WBJETFILES = \
addhel_wbj.o \
qqb_wbjet.o \
qqb_wbjet_g.o \
qqb_wbjet_gs.o \
qqb_wpbjet_gs.o \
qqb_wmbjet_gs.o \
qqb_wbjet_v.o \
qqb_wbjet_z.o
 
WCJETFILES = \
qqb_w_cjet.o \
qqb_w_cjet_g.o \
qqb_w_cjet_gs.o \
qqb_w_cjet_gvec.o \
qqb_w_cjet_v.o \
qqb_w_cjet_z.o \
qqb_w_cjet_massless.o \
qqb_w_cjet_massless_g.o \
wqq_sc.o \
w2jetsq_mass.o \
subqcdm.o

W2JETFILES = \
ampwqqqqg.o \
qqb_w2jet.o \
qqb_w2jet_g.o \
qqb_wp2jet_g.o \
qqb_wm2jet_g.o \
qqb_w2jet_gs.o \
qqb_w2jet_gs_new.o \
qqb_wp2jet_gs.o \
qqb_wp2jet_gs_new.o \
qqb_wm2jet_gs.o \
qqb_wm2jet_gs_new.o \
qqb_w2jet_gvec.o \
qqb_w2jet_gvecx.o \
qqb_w2jet_gvecx_new.o \
qqb_wp2jet_gvecx_new.o \
qqb_wm2jet_gvecx_new.o \
qqb_w2jet_v.o \
qqb_w2jet_z.o \
qqb_w2jetx.o \
qqb_wp2jetx_new.o \
qqb_wm2jetx_new.o \
qqbw2j_loop.o \
storecsv_px.o \
storecsv_qx.o \
subqcd.o \
subqcdn.o \
w2jetnx.o \
w2jetsq.o \
xwqqgg_v.o \
xwqqggg.o 


W2JETVIRTFILES = \
a6g.o \
a61g.o \
a6treeg.o \
fax.o \
faxsl.o \
fcc.o \
fcc_qpgmgpqm.o \
fcc_qpgpgmqm.o \
fcc_qpgpgpqm.o \
fcc_qpgpqmgm.o \
fcc_qpgpqmgp.o \
fcc_qpqmgmgp.o \
fcc_qpqmgpgm.o \
fcc_qpqmgpgp.o \
fsc.o \
fsc1.o \
fsc2.o \
fsc3.o \
fsc4.o \
fsc5.o \
fsc6.o \
fsc7.o \
fsc8.o \
fvf.o \
fvs.o \
vvg.o

WHBBARFILES = \
GammaHbb0.o \
GammaHbb1.o \
dkqqb_wh_g.o \
dkqqb_wh_gs.o \
dkqqb_wh_v.o \
dkqqb_wh_v_massless.o \
qqb_wh.o \
qqb_wh_g.o \
qqb_wh_gs.o \
qqb_wh_v.o \
qqb_wh_z.o\
qqb_wh_HtopEFT.o\

WHGAGAFILES = \
qqb_wh_gaga.o \
qqb_wh_gaga_g.o \
qqb_wh_gaga_gs.o \
qqb_wh_gaga_v.o \

WHWWFILES = \
qqb_wh_ww.o \
qqb_wh_ww_g.o \
qqb_wh_ww_gs.o \
qqb_wh_ww_v.o

WHZZFILES = \
qqb_wh_zz.o \
qqb_wh_zz_g.o \
qqb_wh_zz_gs.o \
qqb_wh_zz_v.o

WWFILES = \
dkWWamps.o \
dkqqb_ww_v.o \
dkqqb_ww_gs.o \
dkqqb_ww_g.o \
c7tree.o \
BigT.o \
L34_12.o \
a6loop.o \
a7tree.o \
amps_anom.o \
b7tree.o \
box1.o \
box3.o \
box5.o \
bub6.o \
comparenum.o \
d1six.o \
d2six.o \
d3six.o \
d4six.o \
fa.o \
gg_ww.o \
massivebox.o \
massivebox6.o \
massivebub.o \
massivetri.o \
massivetri6.o \
mbc.o \
qqb_ww.o \
qqb_ww_g.o \
qqb_ww_gs.o \
qqb_ww_unpol.o \
qqb_ww_v.o \
qqb_ww_z.o \
susana.o \
triangle11new.o \
triangle12.o \
triangle6.o \
triangle7new.o \
triangle9new.o \
vpole.o \
wwamps.o

DIBOSONFILES := $(DIBOSONFILES) \
zgamma_amps.o \
zgamma_amps_vdecay.o \
wgamma_amps.o \
wzgamma_scheme.o \
spinoru_s.o \
spinorz_s.o \
omega_wzgamma.o \
omega_zgamma_vdecay.o \
lumxmsq_wgamma.o \
lumxmsq_zgamma.o \
vvconfig_m.o

WZFILES = \
qqb_wz.o \
qqb_wz_g.o \
qqb_wz_gs.o \
qqb_wz_v.o \
qqb_wz_z.o \
wzamps.o \
dkqqb_wz_g.o \
dksrwz.o \
dksrwzc.o \
dkdrwz.o \

WBBFILES = \
a6.o \
a61.o \
a61LLL.o \
a61LRL.o \
a62.o \
a6routine.o \
a6tree.o \
aqqb_wbb.o \
atrLLL.o \
atrLRL.o \
atree.o \
fpm.o \
fpp.o \
fsl.o \
i3m.o \
lfunctions.o \
lnrat.o \
msqwbb.o \
nagyqqqqg.o \
qqb_wbb.o \
qqb_wbb_g.o \
wbbgamp.o \
qqb_wbb_gs.o \
qqb_wbb_v.o \
qqb_wbb_z.o \
t.o \
t1234.o \
vv.o

WBBMFILES = \
a61mass.o \
a6treemass.o \
Afh.o \
ALC_mm.o \
ALC_mp.o \
ALC_pm.o \
ALC_pp.o \
BDK1211mm.o \
BDK1211mp.o \
BDKfillmm.o \
BDKfillmp.o \
BLC_mm.o \
BLC_mp.o \
BLC_pm.o \
BLC_pp.o \
catani.o \
clearcoeffs.o \
computescalars.o \
qqb_wbbm.o \
qqb_wbbm_g.o \
qqb_wbbm_gs.o \
qqb_wbbm_v.o \
qqb_wbbm_z.o \
spstrng.o \
sumamp.o \
writetable.o \
Wbb.o \
zab.o \
zaba.o \
zabab.o \
zbab.o

WGAMFILES = \
qqb_wgam.o \
qqb_wgam_g.o \
qqb_wgam_gs.o \
qqb_wgam_v.o \
qqb_wgam_z.o \
qqb_wgam_fragdips.o \
qqb_wgam_frag.o \
fagamma.o \
fbgamma.o

WPMZBJFILES = \
qqb_WZjj.o \
WZddidmsq.o \
WZuuidmsq.o \
WZggmsq.o \
TWZggAB.o \
TWZggNAB.o \
TWZggsra.o \
TWZggsrl.o \
WZbbmsq.o \
WZccmsq.o \
TWZbbnr1.o \
TWZbbnr2.o \
TWZbbab.o \
TWZbbnab.o \
TWbbpZab.o \
TWbbmZab.o \
qqb_WZbb.o \
qqb_WZbj.o 

WTFILES = \
BBamps.o \
BBamps_nores.o \
C0fa2m.o \
C0fb2m.o \
I3me.o \
Lsm1_2m.o \
Lsm2_2m.o \
extend_trans_wt.o \
functions.o \
functions1.o \
gs_wc_dg.o \
gs_wt_prog.o \
gs_wt_prog_nores.o \
qb_wtq.o \
dkqqb_w_twdk_g.o \
dkqqb_w_twdk_gs.o \
dkqqb_w_twdk_v.o \
qqb_w_tndk.o \
qqb_w_tndk_g.o \
qqb_w_tndk_gs.o \
qqb_w_tndk_gvec.o \
qqb_w_tndk_v.o \
qqb_w_tndk_z.o \
qqb_w_twdk.o \
qqb_w_twdk_g.o \
qqb_w_twdk_gs.o \
qqb_w_twdk_gvec.o \
qqb_w_twdk_v.o \
qqb_w_twdk_z.o \
qqb_wtbndk.o \
qqb_wtbwdk.o \
Wtoponshell.o \
Watoponshell.o \
tree.o \
virt_mm.o \
virt_mp.o \
virt_pm.o \
virt_pp.o \
vol_wt.o \
wamp.o \
wampd.o \
wcjetn.o \
wtransform_wt.o

ZFILES = \
qqb_z.o \
qqb_z1jet.o \
qqb_z_gs.o \
qqb_z_v.o \
qqb_z_z.o

Z1JETFILES = \
F1anom.o \
qqb_z1jet_gs.o \
qqb_z1jet_v.o \
qqb_z1jet_z.o \
qqb_z_gvec.o

Z2JETFILES = \
A6axBDK.o \
BDKqqbggAxAmp.o \
qqbggAxbox3x12x4.o \
qqbggAxbox3x4x12.o \
qqbggAxslCoeffs.o \
qqbggAxtri123x4x56.o \
qqbggAxtri12x3x456.o \
qqbZgg_floop.o \
qqbZggtree.o \
Ltfunctions.o \
qqb_z2jetx_new.o \
qqb_z2jet_gvecx_new.o \
qqb_z2jet_gs_new.o \
fmtfull.o \
tloop.o \
Bdiff.o \
a61z.o \
a62z.o \
a63.o \
a63z.o \
a6ax.o \
atreez.o \
A6texact.o \
fmt.o \
fzip.o \
Ftexact.o \
makem.o \
makemb.o \
msq_ZqqQQg.o \
msq_z2jetx.o \
ps_check.o \
qqb_z2jet.o \
qqb_z2jet_g.o \
qqb_z2jet_gs.o \
qqb_z2jet_gvec.o \
qqb_z2jet_gvecx.o \
qqb_z2jet_v.o \
qqb_z2jet_z.o \
qqb_z2jetx.o \
storecsz.o \
xzqqgg_v_sym.o \
z2jetsq.o \
z2jetsqn.o

ZBBFILES = \
amp_qqggg.o \
ampqqb_qqb.o \
aqqb_zbb.o \
msq_qqQQg.o \
qqb_zbb.o \
qqb_zbb_g.o \
qqb_zbb_gs.o \
qqb_zbb_gvec.o \
qqb_zbb_v.o \
qqb_zbb_z.o \
xzqqgg.o \
xzqqgg_v.o \
xzqqggg.o

ZBBMFILES = \
gamps0.o \
gampsabc.o \
gampsdef.o \
gampsgh.o \
mamps.o \
qqb_zbbm.o \
qqb_zccm.o

ZGAMFILES = \
gg_zgam.o \
ggZgam_vec.o \
qqb_zgam.o \
qqb_zgam_g.o \
qqb_zgam_gs.o \
qqb_zgam_v.o \
qqb_zgam_z.o \
qqb_zgam_fragdips.o \
qqb_zgam_frag.o

ZGAMJETFILES = \
a6treega.o \
a6treeQL.o \
a6virtQL_slc.o \
qqb_z2jetx_swap.o \
qqb_zaj_frag.o \
qqb_zaj_fragdips.o \
amp_zqqagg_ql.o 

ZAJFILES = \
zaj_utils.o \
zaj_tree_mod.o \
zaj_tree.o \
zaj_tree_vdecay.o \
zaj_virt_mod.o \
zaj_virt_l_mod.o \
zaj_virt.o \
zajj_tree.o \
zajj_tree_vdecay.o \
zajj_tree_mod.o \
omega_ee3jet_m.o \
zajj_anomcoup_mod.o \
zajj_assembly_qqbQQb.o \
zajj_assembly_qqbQQb_vdecay.o \
zajj_assembly_qqbgg.o \
zajj_assembly_qqbgg_vdecay.o \
qqb_zaj_z.o \
qqb_zaj_gs.o \
qqb_zaj_gs_vdecay.o \
qqb_zaj_gvec.o \
qqb_zaj_gvec_vdecay.o

ZGAMGAMFILES = \
a6virtLL.o \
xzqqaa.o \
xzqqaa_v.o \
xzqqaag.o \
qqb_zaa.o \
qqb_zaa_v.o \
qqb_zaa_g.o \
qqb_zaa_gs.o \
qqb_zaa_z.o \
qqb_zaa_frag.o \
qqb_zaa_fragdips.o

WGAJETFILES = \
amp_wqqagg.o \
msq_wqqbagg.o \
qqb_waj_g.o \
xwqqagg.o

WGAMGAMFILES = \
msqWaa.o \
qqb_Waa.o \
spinoruorz.o

#WGAMGAMFILES += $(WGAMGAMFILESNLO)

WGAMGAMFILESNLO = \
onemassbox.o \
easy.o \
hard.o \
Waaintegralfill.o \
Waa_loamp.o \
make_badps.o\
new_coeffs.o\
msqWaa_g.o \
setparams.o \
msqWaa_frm.o \
qqb_Waa_frag_combo.o \
qqb_Waa_g.o \
qqb_Waa_gs.o \
qqb_Waa_v.o \
qqb_Waa_z.o \
Waa_vamps.o \
virt_Qusq.o \
virt_Qdsq.o \
virt_QuQd.o \
virt_QuQe.o \
virt_QdQe.o \
rat_Qusq_pp.o \
rat_Qusq_mp.o \
rat_Qdsq_pp.o \
rat_Qdsq_mp.o \
rat_QuQd_pp.o \
rat_QuQd_mp.o \
rat_QuQe.o \
rat_QdQe.o \
threemtri_coeffs_Wgg.o \
coeff_mp_QuQd_s25.o \
checkpiDpjk.o

ZHBBARFILES = \
ggHZ_mp_box.o\
ggHZ_pp_box.o\
gg_HZ_box.o\
gg_HZ_tri.o\
gg_zh.o\
dkqqb_zh_g.o \
dkqqb_zh_gs.o \
dkqqb_zh_v.o \
dkqqb_zh_v_massless.o\
qqb_ZH_VIItop.o\
qqb_ZH_VItop.o\
qqb_zh.o \
qqb_zh_g.o \
qqb_zh_gs.o \
qqb_zh_v.o \
qqb_zh_z.o

ZHGAGAFILES = \
qqb_zh_gaga.o \
qqb_zh_gaga_g.o \
qqb_zh_gaga_gs.o \
qqb_zh_gaga_v.o \

ZHWWFILES = \
qqb_zh_ww.o \
qqb_zh_ww_g.o \
qqb_zh_ww_gs.o \
qqb_zh_ww_v.o

ZHZZFILES = \
qqb_zh_zz.o \
qqb_zh_zz_g.o \
qqb_zh_zz_gs.o \
qqb_zh_zz_v.o

ZZFILES = \
Acalc.o \
getggHZZamps.o \
getggZZamps.o \
gg_hzz_tb.o \
gg_ZZ_all.o \
ggZZcapture.o \
gg_ZZ.o \
gg_ZZ_Hpi.o \
ggZZmassamp.o \
ggZZmassamp_new.o \
ggZZparsecheck.o \
ggZZwritetable.o \
higgsprop.o \
LRcalc.o \
qg_Cont_ZZj_amp.o \
qg_Hint_ZZ.o \
qg_HZZjet_amp.o \
qqb_zz.o \
qqb_zz_g.o \
qqb_zz_gs.o \
qqb_zz_v.o \
qqb_zz_z.o \
Tri3masscoeff.o \
ZZbox1LL.o \
ZZbox2LL.o \
ZZC012x34LLmp.o \
ZZC01x2LLmp.o \
ZZC01x34LLmp.o \
ZZC02x34LLmp.o \
ZZD02x1x34LLmp.o \
ZZD062x1x34LLmp.o \
zzgamp.o \
ZZintegraleval.o \
ZZmassivebox.o \
ZZmassiveboxtri.o \
ZZmassivebub.o \
ZZmassivetri.o \
ZZmbc.o \
ZZtri12_34LL.o \
ZZtri1_2LL.o \
ZZtri1_34LL.o \
gg_ZZ_int.o 

ZQFILES = \
gQ_zQ.o \
gQ_zQ_g.o \
gQ_zQ_gs.o \
gQ_zQ_v.o \
gQ_zQ_z.o

ZQJETFILES = \
genclust_hqrk.o \
msq_ZqqQQg_noid.o \
qqb_zbjet.o \
qqb_zbjet_g.o \
qqb_zbjet_gvec.o \
qqb_zbjet_gs.o \
qqb_zbjet_v.o \
qqb_zbjet_z.o

APPLGRID= \
gridwrap.o 
#gridwrap_dummy.o 
#mcfm_grid.o \
#fill_APPLgrid.o \

MODFILES= \
scaleset_m.o \
mod_qcdloop_c.o \
qsort_c.o \
types_mod.o \
constants_mod.o \
anomcoup_tbW.o \
jettagging.o \
bbfrac.o \
eftcouple.o \
ddmod.o \
ddext.o

LIBDIR=$(MCFMHOME)/qcdloop-2.0.2/local/lib
LIBFLAGS=-lqcdloop -lov -lpv  -lpvext -lsmallG -lsmallY -lsmallP -lsmallF -lstdc++ -lqdmod -lqd

QDDIR=$(MCFMHOME)/qd-2.3.22/local/lib


# Check NTUPLES flag
ifeq ($(NTUPLES),FROOT)
  USERFILES += mcfm_froot.o froot.cpo
  LIBFLAGS += $(ROOTLIBS)
  NTUPMSG='   ----> MCFM compiled with FROOT n-tuple output <----'
 else
  ifeq ($(NTUPLES),YES)
    ifeq ($(CERNLIB),)
      ERRORMSG=Please specify the path to CERNLIB to use n-tuples
      $(error $(ERRORMSG))
    endif
    USERFILES += dswhbook.o
    LIBDIR=$(CERNLIB)
# Note: some versions of the CERN libraries (e.g. Redhat v7.2)
#       require that "-lpacklib" be changed to "-lpacklib_noshift"
    LIBFLAGS += -lmathlib -lpacklib -lkernlib
    NTUPMSG='   ----> MCFM compiled with optional n-tuple output <----'
  else
    ifeq ($(NTUPLES),NO)
      USERFILES += dsw_dummy.o
    else
      ERRORMSG=Please set NTUPLES equal to NO/YES/FROOT
      $(error $(ERRORMSG))
    endif
    NTUPMSG = '   ----> MCFM compiled with histogram output only <----'
  endif
endif

ifeq ($(PDFROUTINES),LHAPDF)
   PARTONFILES += \
   alfamz_lhapdf.o \
   fdist_lhapdf.o \
   pdfwrap_lhapdf.o
   LIBFLAGS += -lLHAPDF
   LIBFLAGS += $(shell lhapdf-config --ldflags)
   ifeq ($(LHAPDFLIB),)
   else
     LIBDIR += -Wl,-rpath=$(LHAPDFLIB) -L$(LHAPDFLIB)
   endif
   PDFMSG='   ----> MCFM compiled with LHAPDF routines <----'
else
ifeq ($(PDFROUTINES),NATIVE)
   PARTONFILES += \
   alfamz.o \
   CT10Pdf.o \
   CT14Pdf.o \
   CT14Pdfqed.o \
   POLINT4F.o \
   Ctq4Fn.o \
   Ctq5Par.o \
   Ctq5Pdf.o \
   Cteq6Pdf-2008.o \
   cteq3.o \
   mrs96.o \
   mrs98.o \
   mrs98lo.o \
   mrs98ht.o \
   mrs99.o \
   mrsebh.o \
   mrsg.o \
   mrst2001lo.o \
   jeppelo.o \
   mrst2001.o \
   mrst2002.o \
   mrst2004.o \
   mrst2004f3.o \
   mrst2004f4.o \
   mrst2004qed.o \
   mstwpdf.o \
   MMHTmstwpdf.o \
   mt.o \
   NNPDFDriver.o \
   fdist_linux.o \
   pdfwrap_linux.o \
   eks98r.o
   PDFMSG='   ----> MCFM compiled with its own PDFs only <----'
else
   ERRORMSG=Please set PDFROUTINES equal to NATIVE/PDFLIB/LHAPDF
   $(error $(ERRORMSG))
endif
endif

OURCODE = $(LIBFILES) $(NEEDFILES)  $(PROCDEPFILES) $(SPINORFILES) \
          $(PHASEFILES) $(SINGLETOPFILES) $(SINGLETOP2FILES) $(HELASFILES) \
          $(TOPHFILES) $(TOPZFILES) $(TOPWFILES) $(TOPDKFILES) \
          $(USERFILES) $(VOLFILES) $(WFILES) $(W2JETFILES) \
          $(WCJETFILES) $(WBFROMCFILES) $(WBJETFILES) \
	  $(W2JETVIRTFILES) $(WHBBARFILES) $(WGAMFILES) $(ZGAMFILES) \
          $(WWFILES) $(WZFILES) $(ZFILES) $(ZHBBARFILES) \
		  $(DIBOSONFILES) \
		$(WGAJETFILES)\
          $(ZZFILES) $(ZGFILES) $(W1JETFILES) $(Z2JETFILES) \
	    $(Z1JETFILES) $(HWWFILES) $(HZZFILES) $(VVFILES) \
          $(TAUTAUFILES) $(HTTBARFILES) $(HJETFILES) \
          $(BBHIGGSFILES) $(WBBFILES) $(ZBBFILES) \
          $(QQHFILES) $(QQHWWFILES) $(QQHZZFILES) $(GGHFILES) $(GGHGFILES) \
          $(GGHGGREALFILES) $(GGHGGvirtFILES) $(H4PCODEFILES) \
	  $(TOPFILES) $(ZQFILES) $(ZQJETFILES) $(WTFILES) $(HWWJETFILES) \
	  $(WHWWFILES) $(ZHWWFILES) $(WHZZFILES) $(ZHZZFILES) \
	  $(WHGAGAFILES) $(ZHGAGAFILES) $(QQHGAGAHFILES) \
	  $(STOPBFILES) $(STOPJETFILES) $(EPEM3JFILES) $(QQZTTFILES) \
	  $(GAMGAMFILES) $(DIRGAMFILES) \
	  $(WBBMFILES) $(ZBBMFILES) $(FRAGFILES) \
	  $(TOPDKBSYFILES) $(TOPDECAYFILES) \
	  $(GGHZGAFILES) $(ZGAMJETFILES) $(ZAJFILES) $(ZGAMGAMFILES) \
	  $(WGAMGAMFILES) $(HHFILES) \
	  $(SINGLETOPHFILES) $(SINGLETOPZFILES) \
	  $(DMFILES) $(TRIGAMFILES) $(FOURGAMFILES) \
	  $(WPMZBJFILES) $(MULTICHANFILES) \
	  $(PWGPLOTSFILES) \
	  $(CHECKINGFILES) $(UTOOLSFILES) $(WBFFILES) \
          $(SCETFILES) $(SCET0jFILES) $(TDHPLFILES) $(WH1JETFILES) $(EWFILES) \
		  $(INTEGRATEFILES) $(MODFILES) $(APPLGRID)
          
OTHER = $(PARTONFILES) $(WPWP2JFILES) $(F90FILES) 

ALLMCFM = $(OURCODE) $(MAIN) $(OTHER)

ifeq ($(USEMPI),YES)
else
  OURCODE += mpi_stubs.o
endif

.DEFAULT_GOAL := Bin/mcfm_omp

#FTNCHEKPATH = /home/ellis/Fortran/Ftnchek/ftnchek-3.1.2
#FORCHKPATH = /home/ellis/bin/

# Specify the dependencies of the .o files and the rules to make them.

FOROPTS = -include=$(INCPATH) -nonovice -nopretty -quiet

#.SUFFIXES: .prj

# improved so that .prj files are moved out of src directory and
# into base, only for .f files that don't already exist there
#.f.prj:
#		$(FTNCHEKPATH)/ftnchek -project -noextern\
#            $(FOROPTS) $< ; \
#            if ! [ -e $(MCFMHOME)/$(notdir $<) ] ; then \
#            mv $(basename $<).prj $(MCFMHOME) ; fi

#PRJS =      $(OURCODE:.o=.prj)
#check:      $(PRJS)
#		$(FTNCHEKPATH)/ftnchek $(FOROPTS) $(PRJS)
#PRJSF =      $(OURCODE:.o=.f)
#checkf:
#		$(FORCHKPATH)/forchk -allc -I $(INCPATH) $(PRJSF)



Bin/mcfm_omp: $(ALLMCFM)
	@$(FC) $(FFLAGS) -L$(LIBDIR) -L$(QDDIR) -L$(PVDIR) -L$(RECURDIR) -L$(OVDIR) -o $@ \
	$(patsubst %,$(OBJNAME)/%,$(ALLMCFM)) $(LIBFLAGS) $(LDFLAGS)
	@echo $(PDFMSG)
	@echo $(NTUPMSG)
	@echo "   ----> Executable is mcfm_omp <----"

%.cpo: %.cpp
	$(CXX) -std=c++11 -c $(CXXFLAGS) $(INCPATH) -o $(OBJNAME)/$@ $<

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $(OBJNAME)/$@ $<

%.o: %.f
	$(FC) $(FFLAGS) -c -o $(OBJNAME)/$@ $<

clean:
	- rm -f obj/*.o obj/*.cpo obj/*.mod Bin/mcfm_omp

debug: Bin/mcfm_omp
	cd $(BIN) && OMP_NUM_THREADS=1 ./mcfm_omp

run: Bin/mcfm_omp
	cd $(BIN) && ./mcfm_omp

# -----------------------------------------------------------------------------

include depends.mk
