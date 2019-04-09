FORCHKPATH = /usr/local/bin/
FFLAGS 	= -O2 -g
HERE = $(PWD)
INCDIR = $(HERE)/../../Include/ $(HERE)/../Include
FC = gfortran 
LIBRARY	      = libsmallP.a


# -----------------------------------------------------------------------------
# Specify the object files. 

OBJS = \
Cfill_recur3.o \
Dfill_recur3.o \
runCP_0.o \
runCP_00.o \
runCP_0000.o \
runCP_0000i.o \
runCP_00i.o \
runCP_00ii.o \
runCP_00iii.o \
runCP_00iiii.o \
runCP_i.o \
runCP_ii.o \
runCP_iii.o \
runCP_iiii.o \
runCP_iiiii.o \
runP_0.o \
runP_00.o \
runP_0000.o \
runP_0000i.o \
runP_00i.o \
runP_00ii.o \
runP_00iii.o \
runP_00iiii.o \
runP_i.o \
runP_ii.o \
runP_iii.o \
runP_iiii.o \
runP_iiiii.o


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


