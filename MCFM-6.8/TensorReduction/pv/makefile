FORCHKPATH = /home/ellis/bin/
FFLAGS 	= -O2 -g
HERE = $(PWD)
INCDIR = $(HERE)/../Include/
FC = gfortran 
LIBRARY	      = libpv.a


# -----------------------------------------------------------------------------
# Specify the object files. 

OBJS = \
TRsetmaxindex.o \
pvCcheck.o \
pvDcheck.o \
pvEcheck.o \
checkaccuracy.o \
pvsetmudim.o \
pvsymfunctions.o \
pvdsvdcmp.o \
pvdpythag.o \
pvclearcache.o \
pvGramsing.o \
pvBackSubst.o \
pvIpiksrt.o \
pvGaussPivot.o \
pvXdet.o \
pvXInverse.o \
pvAtensor.o \
pvAcache.o \
pvAfill.o \
pvBtensor.o \
pvBcache.o \
pvBfill.o \
pvCtensor.o \
pvCcache.o \
pvCfill.o \
pvDtensor.o \
pvDcache.o \
pvDfill.o \
pvEtensor.o \
pvE0scalar.o \
pvvcalc.o \
pvdot.o \
pvfndd.o \
pvswitch1.o \
pvswitch2.o \
pvswitch3.o \
pvswitch4.o \
pvswitch5.o \
pvswitch6.o \
pvArraysetup.o \
pvYcalc.o \
pvswap.o


$(LIBRARY):	$(OBJS) 
		ar cru $(LIBRARY) $(OBJS) 
		ranlib $(LIBRARY)
 
.f.o:	
		$(FC) $(F77FLAGS) -I$(INCDIR) -c $*.f

# -----------------------------------------------------------------------------

all:		$(LIBRARY)

# -----------------------------------------------------------------------------

PRJS =      $(OBJS:.o=.f) 
PRJSF =      $(OBJS:.o=.f) 

check:      $(PRJS) 
		$(FTNCHEKPATH)/ftnchek -include=$(INCDIR) $(FOROPTS) $(PRJS)

checkf:      
		$(FORCHKPATH)/forchk -I $(INCDIR) $(PRJSF)


clean:
	- rm -f *.o *.s *.prj *.a *~ core

# -----------------------------------------------------------------------------

# DO NOT DELETE


