FORCHKPATH = /usr/local/bin/
FFLAGS 	= -O2 -g
HERE = $(PWD)
INCDIR = $(HERE)/../../Include/ $(HERE)/../Include
FC = gfortran 
LIBRARY	      = libsmallF.a


# -----------------------------------------------------------------------------
# Specify the object files. 

OBJS = \
Cfill_recur4.o \
Dfill_recur4.o \
runCF_0.o \
runCF_00.o \
runCF_00i.o \
runCF_00ii.o \
runCF_00iii.o \
runCF_00iiii.o \
runCF_00iiiii.o \
runCF_i.o \
runCF_ii.o \
runCF_iii.o \
runCF_iiii.o \
runCF_iiiii.o \
runF_0.o \
runF_00.o \
runF_00i.o \
runF_00ii.o \
runF_00iii.o \
runF_00iiii.o \
runF_00iiiii.o \
runF_i.o \
runF_ii.o \
runF_iii.o \
runF_iiii.o \
runF_iiiii.o


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


