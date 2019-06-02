# Makefile generated automatically with
# ../scripts/makef95makefile "prec_and_timing" "-L../src" "-lhoppet_v1" "-I../src" "-remake" "./mkmk"
# default program to compile
PROG =	prec_and_timing

ALLPROG = 	benchmark prec_and_timing

# This will be used one day...
ALLPROGSRC =	benchmark.f90 prec_and_timing.f90

ALLPROGOBJ =	benchmark.o prec_and_timing.o

SRCS =	io_utils.f90 lcl_dec.f90 NameSelect.f90

POSTSRCS =	

OBJS =	io_utils.o lcl_dec.o NameSelect.o

POSTOBJS =	
POSTLIB = 

LIBS = -L../src -lhoppet_v1  

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = -O3 -fPIC -I../src 
F90 = gfortran
F90FLAGS = -O3 -fPIC -I../src 
LDFLAGS =  

all: $(ALLPROG)

ALL:  $(ALLPROG)

benchmark: benchmark.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o benchmark benchmark.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

prec_and_timing: prec_and_timing.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o prec_and_timing prec_and_timing.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

libclean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) 

clean:
	rm -f  $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

distclean: realclean
	rm -f Makefile
realclean:
	rm -f  $(ALLPROG) $(ALLPROGOBJ) $(OBJS) $(POSTOBJS) *.mod *.d

make:
	./mkmk

install:

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $<

benchmark.o: NameSelect.o io_utils.o
io_utils.o: 
NameSelect.o: io_utils.o
prec_and_timing.o: io_utils.o
