# Makefile generated automatically with
# ../../scripts/makef95makefile "compare2files" "-L../../src" "-lhoppet_v1" "-I../../src" "-remake" "./mkmk"
# default program to compile
PROG =	compare2files

ALLPROG = 	compare2files compare2files_v2

# This will be used one day...
ALLPROGSRC =	compare2files.f90 compare2files_v2.f90

ALLPROGOBJ =	compare2files.o compare2files_v2.o

SRCS =	io_utils.f90 lcl_dec.f90

POSTSRCS =	

OBJS =	io_utils.o lcl_dec.o

POSTOBJS =	
POSTLIB = 

LIBS = -L../../src -lhoppet_v1  

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = -O3 -fPIC -I../../src 
F90 = gfortran
F90FLAGS = -O3 -fPIC -I../../src 
LDFLAGS =  

all: $(ALLPROG)

ALL:  $(ALLPROG)

compare2files: compare2files.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o compare2files compare2files.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

compare2files_v2: compare2files_v2.o $(OBJS) $(POSTOBJS) 
	$(F90) $(LDFLAGS) -o compare2files_v2 compare2files_v2.o $(OBJS) $(LIBS) $(POSTOBJS) $(POSTLIB)

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

compare2files.o: io_utils.o
compare2files_v2.o: io_utils.o
io_utils.o: 
