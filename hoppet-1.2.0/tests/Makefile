# Makefile generated automatically with
# ../scripts/makef95makefile "lumi_tests" "-I../src" "-L../src" "-lhoppet_v1" "-lio"
# default program to compile
PROG =	lumi_tests

ALLPROG = 	bidirectional_mass_thresholds coupling determine_accuracies \
	intrinsic_charm lumi_tests msbar_tests split_mat_ops \
	streamlined_multisplit

# This will be used one day...
ALLPROGSRC =	bidirectional_mass_thresholds.f90 coupling.f90 \
	determine_accuracies.f90 intrinsic_charm.f90 lumi_tests.f90 \
	msbar_tests.f90 split_mat_ops.f90 streamlined_multisplit.f90

ALLPROGOBJ =	bidirectional_mass_thresholds.o coupling.o determine_accuracies.o \
	intrinsic_charm.o lumi_tests.o msbar_tests.o split_mat_ops.o \
	streamlined_multisplit.o

SRCS =	

POSTSRCS =	

OBJS =	

POSTOBJS =	
POSTLIB = 

LIBS = -L/Users/gsalam/utils/libf90-gfortran -liof90 -L../src -lhoppet_v1  

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = -O3 -fPIC -I../src -I/Users/gsalam/utils/libf90-gfortran 
F90 = gfortran
F90FLAGS = -O3 -fPIC -I../src -I/Users/gsalam/utils/libf90-gfortran 
LDFLAGS =  

all: $(ALLPROG)

ALL:  $(ALLPROG)

bidirectional_mass_thresholds: bidirectional_mass_thresholds.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o bidirectional_mass_thresholds bidirectional_mass_thresholds.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

coupling: coupling.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o coupling coupling.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

determine_accuracies: determine_accuracies.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o determine_accuracies determine_accuracies.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

intrinsic_charm: intrinsic_charm.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o intrinsic_charm intrinsic_charm.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

lumi_tests: lumi_tests.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o lumi_tests lumi_tests.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

msbar_tests: msbar_tests.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o msbar_tests msbar_tests.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

split_mat_ops: split_mat_ops.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o split_mat_ops split_mat_ops.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

streamlined_multisplit: streamlined_multisplit.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o streamlined_multisplit streamlined_multisplit.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

distclean: realclean
	rm -f Makefile
realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	../scripts/makef95makefile  "lumi_tests" "-I../src" "-L../src" "-lhoppet_v1" "-lio"

install:

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

bidirectional_mass_thresholds.o: 
coupling.o: 
determine_accuracies.o: 
intrinsic_charm.o: 
lumi_tests.o: 
msbar_tests.o: 
split_mat_ops.o: 
streamlined_multisplit.o: 
