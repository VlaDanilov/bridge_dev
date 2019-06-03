FTNCHEKPATH = /home/ellis/Fortran/Ftnchek/ftnchek-3.1.2
FORCHKPATH = /usr/local/bin/
FFLAGS 	= -g -Wall 
HERE = $(PWD)
LIBNAME = QCDLoop
LIBDIR = $(HERE)/Lib
INCDIR = $(HERE)/pv/Include
PVDIR  = $(HERE)/pv
PVEXTDIR  = $(HERE)/pvext
OVDIR  = $(HERE)/ov
QLDIR  = $(HERE)/../QCDLoop/ql
FFDIR  = $(HERE)/../QCDLoop/ff
RRDIR  = $(HERE)/recur
ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif

libs:
	(cd $(RRDIR)/smallG; $(MAKE) -f makefile) 
	(cd $(RRDIR)/smallY; $(MAKE) -f makefile) 
	(cd $(RRDIR)/smallP; $(MAKE) -f makefile) 
	(cd $(RRDIR)/smallF; $(MAKE) -f makefile) 
	(cd $(PVDIR); $(MAKE) -f makefile) 
	(cd $(PVEXTDIR); $(MAKE) -f makefile) 
	(cd $(OVDIR); $(MAKE) -f makefile) 

clean:
	(cd $(RRDIR)/smallG; $(MAKE) -f makefile clean) 
	(cd $(RRDIR)/smallY; $(MAKE) -f makefile clean) 
	(cd $(RRDIR)/smallP; $(MAKE) -f makefile clean) 
	(cd $(RRDIR)/smallF; $(MAKE) -f makefile clean) 
	(cd $(PVDIR); $(MAKE) -f makefile clean) 
	(cd $(PVEXTDIR); $(MAKE) -f makefile clean) 
	(cd $(OVDIR); $(MAKE) -f makefile clean) 

