// emacs: this is -*- c++ -*-
//
//   @file    mcfm_grid.h        
//
//                   
//  
//   Copyright (C) 2011 M.Sutton (sutt@cern.ch)    
//
//   $Id: mcfm_grid.h, v0.0   Fri 12 Aug 2011 07:50:44 CEST sutt $


#ifndef __MCFM_GRID_H
#define __MCFM_GRID_H

#include <iostream>

#include "appl_grid/appl_grid.h"





/// fortran structures and common blocks...

static const int __nf__   = 5;
static const int __nf2__  = 11;
static const int __maxd__ = 41;


typedef struct {
  double weightfactor;
  double weightb [ __nf2__ ][ __nf2__ ];
  double weightv [ __nf2__ ][ __nf2__ ];
  double weightv1[ __nf2__ ][ __nf2__ ];
  double weightv2[ __nf2__ ][ __nf2__ ];
  double weightr [ __nf2__ ][ __nf2__ ][ __maxd__ ];
} __gridweight__;


typedef struct {
  double Vud, Vus, Vub, Vcd, Vcs, Vcb; 
} __cabib__;


typedef struct {
  double vsq[ __nf2__ ][ __nf2__ ], vsum[ __nf2__ ];
} __ckm__;


typedef struct  {
  double ag_xx1,ag_xx2,ag_x1z,ag_x2z,ag_scale,refwt,refwt2;
  int    contrib, dipole;
} __gridevent__;


typedef struct  { 
  int nproc;
} __nproc__;


typedef struct {
  int nflav;
} __nflav__;


typedef struct  {
    double sqrts;
} __energy__;



typedef struct  {
    int itmx1, ncall1, itmx2, ncall2;
} __iterat__;


 

extern "C" __ckm__   ckm_;
extern "C" __cabib__ cabib_;
extern "C" __gridevent__  gridevent_;
extern "C" __gridweight__ gridweight_;
extern "C" __nproc__  nproc_;
extern "C" __nflav__  nflav_;
extern "C" __energy__ energy_;
extern "C" __iterat__ iterat_;



namespace appl { 

class mcfm_grid : public  grid {

public:

  //  mcfm_grid() { } 

  mcfm_grid(int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0,  int Q2order=5,  
	    int Nx=50,   double xmin=1e-5,     double xmax=0.9,          int xorder=5,
	    int Nobs=20, double obsmin=100.0,  double obsmax=7000.0, 
	    std::string genpdf="mcfm_pdf", 
	    int leading_order=0, int nloops=1, 
	    std::string transform="f2")  
    : grid( NQ2,  Q2min,  Q2max,  Q2order, Nx, xmin, xmax, xorder,
	    Nobs, obsmin, obsmax, genpdf,  leading_order,  nloops, 
	    transform )  {  }

  mcfm_grid( int Nobs, const double* obsbins,
	     int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0, int Q2order=5,
	     int Nx=50,   double xmin=1e-5,     double xmax=0.9,         int xorder=5, 
	     std::string genpdf="mcfm_pdf",
	     int leading_order=0, int nloops=1, 
	     std::string transform="f2" )
    : grid( Nobs, obsbins, NQ2,  Q2min,  Q2max,  Q2order, Nx, xmin, xmax, xorder,
	    genpdf,  leading_order,  nloops, 
	    transform )  {  }
  
  mcfm_grid( const std::vector<double> obs,
	     int NQ2=50,  double Q2min=10000.0, double Q2max=25000000.0,   int Q2order=5, 
	     int Nx=50,   double xmin=1e-5,     double xmax=0.9,           int xorder=5, 
	     std::string genpdf="mcfm_pdf", 
	     int leading_order=0, int nloops=1, 
	     std::string transform="f2" )
    : grid( obs, NQ2,  Q2min,  Q2max,  Q2order, Nx, xmin, xmax, xorder,
	    genpdf,  leading_order,  nloops, 
	    transform )  {  }
  
  // build a grid but don't build the internal igrids - these can be added later
  mcfm_grid( const std::vector<double> obs,
	     std::string genpdf="nlojet_pdf", 
	     int leading_order=0, int nloops=1, 
	     std::string transform="f2" )
    : grid( obs,  genpdf,  leading_order,  nloops, transform )  {  }

  
  mcfm_grid(const grid& g) : grid(g) {  }
  
  mcfm_grid(const std::string& filename="./grid.root", const std::string& dirname="grid") : grid(filename, dirname) {  }
  

  // destructor
  virtual ~mcfm_grid() { } 

  
  //  fill weight from MCFM common block  
  void fillMCFM( double obs);
  void collectWeight   ( const int, const int, double* wt, int =0 );
  void decideSubProcess( const int, const int, int&, double&, int =0 );
  

};


};



#endif  // __MCFM_GRID_H 










