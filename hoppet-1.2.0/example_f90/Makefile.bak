# Makefile generated automatically with
# ../scripts/makef95makefile "tabulation_example" "-L../src" "-lhoppet_v1" "-I../src" "-remake" "./mkmk"
# default program to compile
PROG =	tabulation_example

ALLPROG = 	sumrules tabulation_example tabulation_example_2 \
	tabulation_example_streamlined

# This will be used one day...
ALLPROGSRC =	sumrules.f90 tabulation_example.f90 tabulation_example_2.f90 \
	tabulation_example_streamlined.f90

ALLPROGOBJ =	sumrules.o tabulation_example.o tabulation_example_2.o \
	tabulation_example_streamlined.o

SRCS =	

POSTSRCS =	

OBJS =	

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

sumrules: sumrules.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o sumrules sumrules.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

tabulation_example: tabulation_example.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o tabulation_example tabulation_example.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

tabulation_example_2: tabulation_example_2.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o tabulation_example_2 tabulation_example_2.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

tabulation_example_streamlined: tabulation_example_streamlined.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o tabulation_example_streamlined tabulation_example_streamlined.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

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

sumrules.o: 
tabulation_example.o: 
tabulation_example_2.o: 
tabulation_example_streamlined.o: 
