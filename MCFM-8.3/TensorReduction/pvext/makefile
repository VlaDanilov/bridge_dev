# Makefile routine.
VPATH		= $(PWD)/obj
INCPATH = $(PWD)/../Include/
QCDLOOP_DIR=$(PWD)/../../obj
OUTPUT_OPTION	= -o $(PWD)/obj/$@ 

ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif
FFLAGS 	= -O2 -Wall -fno-f2c -fopenmp -g -I$(INCPATH) -Iobj -Jobj -I$(QCDLOOP_DIR)

LIBRARY	      = libpvext.a

objs = \
pvextAfill.o \
pvextBfill.o \
pvextCfill.o \
pvextDfill.o \
pvextAcache.o \
pvextBcache.o \
pvextCcache.o \
pvextDcache.o \
pvextBtensor.o \
pvextCtensor.o \
pvextDtensor.o \
pvextclearcache.o
#pvextCcheck.o \
#pvextDcheck.o \
#TRsetmaxindex.o \
#pvEtensor.o \
#pvE0scalar.o \
#pvvcalc.o \
#pvdot.o \
#pvfndd.o \
#pvswitch1.o \
#pvswitch2.o \
#pvswitch3.o \
#pvswitch4.o \
#pvswitch5.o \
#pvswitch6.o \
#pvArraysetup.o \
#pvYcalc.o \
#pvswap.o
#pvCcheck.o \
#pvDcheck.o \
#pvEcheck.o \
#checkaccuracy.o \
#pvsetmudim.o \
#pvsymfunctions.o \
#pvdsvdcmp.o \
#pvdpythag.o \
#pvGramsing.o \
#pvBackSubst.o \
#pvIpiksrt.o \
#pvGaussPivot.o \
#pvXdet.o \
#pvXInverse.o \
#pvAtensor.o \

$(LIBRARY):	$(objs) 
		ar cru $(LIBRARY) $(patsubst %,obj/%,$(objs)) 
		ranlib $(LIBRARY)
 
all:		$(LIBRARY)

clean:
	- rm -f obj/*.o *.a

