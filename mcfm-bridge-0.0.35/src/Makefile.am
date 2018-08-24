##########################################################################
#
#   File:         Makefile.am
#
#                 Copyright (C) M.Sutton (sutt@cern.ch) 
#
#   Description:  generate the Makefile for the amcatnlo - applgrid bridge
#                 
#
#   Created:    M.Sutton (sutt@cern.ch) Fri 10 July 2013 00:05:29 GMT
#   
#   $Id: makefile 20 2013-07-10 01:55:42Z sutt $                
#                   
#########################################################################


APPLCFLAGS = $(shell applgrid-config --cxxflags)
APPLLIBS   = $(shell applgrid-config --ldflags)

HOPPETCFLAGS = $(shell hoppet-config --cxxflags)
HOPPETLIB    = $(shell hoppet-config --libs)

LHAPDFCXXFLAGS = $(shell lhapdf-config --cppflags)
LHAPDFLIB      = $(shell lhapdf-config --ldflags)
LHAPDFPATH     = $(shell lhapdf-config --pdfsets-path)

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --nonew --libs)
ROOTGLIBS  = $(shell root-config --nonew --glibs)
ROOTARCH   = $(shell root-config --ldflags)

# NB For some *stupid* reason, although hoppet requests these two functions, 
#    the linker doesn't find them on lxplus - becuas ethe STUPID GFORTRAN
#    installation is broken 
FRTLIB   = $(shell gfortran -print-file-name=libgfortran.a)
# FRTLLIB  = -u _gfortran_runtime_error_at  -u _gfortran_os_error -L$(subst /libgfortran.a, ,$(FRTLIB) ) -lgfortran 
FRTLLIB  = -L$(subst /libgfortran.a, ,$(FRTLIB) ) -lgfortran 



bin_PROGRAMS        = standSimple
standSimple_SOURCES = standSimple.cxx 
standSimple_CXXFLAGS = -O2 -Wall -fPIC -DPDFPATH=\"$(LHAPDFPATH)\" -I. -I.. -I$(srcdir)  $(LHAPDFCXXFLAGS) $(ROOTCFLAGS) $(APPLCFLAGS) $(ROOTCFLAGS) $(HOPPETCFLAGS)
standSimple_LDFLAGS  = $(ROOTARCH) -O2  $(ROOTLIBS) $(APPLLIBS) $(LHAPDFLIB) $(HOPPETLIB) $(FRTLLIB) $(FRTLIB) 
# -L$(HOME)/lib -lAS

lib_LIBRARIES= libmcfmbridge.a 

libmcfmbridge_a_SOURCES = mcfm_interface.cxx mcfm_grid.cxx mcfm_grid.h LHAPDF.h

# library_includedir=$(includedir)/mcfmbridge
# library_include_HEADERS=



AM_CXXFLAGS = -O2 -Wall -fPIC  -I. -I.. -I$(srcdir)  $(ROOTCFLAGS) $(APPLCFLAGS) $(HOPPETCFLAGS) $(LHAPDFCXXFLAGS)
AM_CFLAGS   = -O2 -Wall -fPIC  -I. -I.. -I$(srcdir)  $(ROOTCFLAGS) $(APPLCFLAGS) $(HOPPETCFLAGS) $(LHAPDFCXXFLAGS)

AM_LDFLAGS = $(ROOTARCH) -O2  $(ROOTLIBS) $(FCLIBS) $(APPLLIBS) $(LHAPDFLIB) $(FRTLLIB) $(HOPPETLIB)

AM_FCFLAGS = -c
AM_FCLIBS  = $(FRTLLIB)
AM_SOFLAGS = -shared

clean-local:
	rm -rf *.o *.lo 

.cxx.o : 
	$(CXX) $(AM_CXXFLAGS) -c $<

.c.o :
	$(CC) $(AM_CFLAGS) -c $<

