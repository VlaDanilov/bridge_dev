ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif
FFLAGS 	= -g

recurdir:
	(cd smallG; $(MAKE)) 
	(cd smallY; $(MAKE)) 
	(cd smallP; $(MAKE)) 
	(cd smallF; $(MAKE)) 

# -----------------------------------------------------------------------------

clean:
	- rm -f *.a */*.o

# -----------------------------------------------------------------------------

# DO NOT DELETE


