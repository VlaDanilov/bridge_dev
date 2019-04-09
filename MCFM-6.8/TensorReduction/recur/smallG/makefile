FORCHKPATH = /usr/local/bin/
FFLAGS 	= -O2 -g
HERE = $(PWD)
INCDIR = $(HERE)/../../Include/ $(HERE)/../Include
FC = gfortran 
LIBRARY	      = libsmallG.a


# -----------------------------------------------------------------------------
# Specify the object files. 

OBJS = \
determinant.o \
Array2dim.o \
Array3dim.o \
Carraysetup.o \
Darraysetup.o \
Cfill_recur.o \
Dfill_recur.o \
runC_0.o \
runC_00.o \
runC_0000.o \
runC_000000.o \
runC_0000i.o \
runC_0000ii.o \
runC_00i.o \
runC_00ii.o \
runC_00iii.o \
runC_00iiii.o \
runC_i.o \
runC_ii.o \
runC_iii.o \
runC_iiii.o \
runC_iiiii.o \
run_0.o \
run_00.o \
run_0000.o \
run_000000.o \
run_0000i.o \
run_0000ii.o \
run_00i.o \
run_00ii.o \
run_00iii.o \
run_00iiii.o \
run_i.o \
run_ii.o \
run_iii.o \
run_iiii.o \
run_iiiii.o


$(LIBRARY):	$(OBJS) 
		ar cru $(LIBRARY) $(OBJS) 
		ranlib $(LIBRARY)
		mv $(LIBRARY) ..
 

.f.o:	
		$(FC) $(F77FLAGS) $(addprefix -I,$(INCDIR)) -c $*.f


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


