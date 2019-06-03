VPATH		= $(PWD)/obj
INCPATH = $(PWD)/../Include/
QCDLOOP_DIR=$(PWD)/../../obj
OUTPUT_OPTION	= -o $(PWD)/obj/$@ 

ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif
FFLAGS 	= -O2 -fno-f2c -fopenmp -g -I$(INCPATH) -Iobj -Jobj -I$(QCDLOOP_DIR)

LIBRARY	      = libpv.a

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
		ar cru $(LIBRARY) $(patsubst %,obj/%,$(OBJS))
		ranlib $(LIBRARY)
 
all:		$(LIBRARY)

clean:
	- rm -f obj/*.o *.a

