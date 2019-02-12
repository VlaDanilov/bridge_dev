//
//   mcfm_procmap.cxx        
//
//                   
// 
//   2017 V.Danilov (vladyslav.danilov@cern.ch)    
//
//   $Id: mcfm_procdev.cxx, v   Fri  11 Mai 2017 10:45


#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <map>
#include <sys/stat.h>

#include <cstdlib> 
#include <sys/time.h> 


#include "TFile.h"
#include "TH1D.h"
 #include "TMatrixT.h"
#include "TVectorT.h"
#include "TString.h"

#include "mcfm_procmap.h"
#include "mcfm_grid.h"

info whichProcess(int  proc)
{
  std::map<int, info> whPr;
  
  double x_low  = 1.0e-9;
  double x_up   = 1.0   ; 
  double q2_Low = 0.0   ;
  double q2_Up  = 0.0   ;
  unsigned int nq2_Bins    = 0 ;
  unsigned int q_order     = 0 ;
  unsigned int LowestOrder = 0 ;  
  unsigned int n_XBins     = 0 ;
  unsigned int x_Order     = 0 ;  
  unsigned int n_loops     = 0 ;

  n_XBins = 40;
  x_Order = 6 ; 
  n_loops = 1 ;
 
  // W production
  q2_Low = 6299.99;
  q2_Up  = 6800.01;
  nq2_Bins = 3;
  q_order  = 1;
  LowestOrder = 0;  
  whPr[1]   = {	"W+ production", "mcfmwp.config", "-Wplus", 		  q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[6]   = {	"W- production", "mcfmwm.config", "-Wminus",  		  q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  
  q2_Low = 1.0;
  q2_Up  = 4000*4000;
  nq2_Bins = 15;
  q_order  = 3;
  LowestOrder = 0;
  whPr[11] = {	" W+ + jet production", "mcfm-wpjet", "-WplusJet", 	  q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder,  n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[13] = {         " W+ + Cbar production", "mcfm-wpc", "-WplusCbar", q2_Low, q2_Up, nq2_Bins, q_order, 1,		  n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[16] = {	" W- + jet production", "mcfm-wmjet", "-WminusJet", 	  q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder,  n_XBins, x_Order, x_low, x_up, n_loops}; 
  whPr[18] = {	"W- + C production", "mcfm-wmc", "-WminusC", 	 	  q2_Low, q2_Up, nq2_Bins, q_order, 1, 		  n_XBins, x_Order, x_low, x_up, n_loops};
  
  //
  //Z production  
  q2_Low = 8280.99;
  q2_Up  = 8281.01;
  nq2_Bins = 3;
  q_order  = 1;
  LowestOrder = 0;  
  whPr[31] = {	" Z production", "mcfm-z", "-Z0",                         q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  q2_Low = 1.0;
  q2_Up  = 4000*4000;
  nq2_Bins = 15;
  q_order  = 3;
  whPr[41] = {	" Z-jet production", "mcfm-zjet",  "-Zjet_41",            q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[42] = {	" Z-jet production", "mcfm-zjet",  "-Zjet_42", 	          q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[43] = {	" Z-jet production", "mcfm-zjet",  "-Zjet_43", 		  q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops}; 
  
  //
  //Higgs
  q2_Low = 1.0;
  q2_Up  = 4000*4000;
  nq2_Bins = 15;
  q_order  = 3;
  //LowestOrder = 2;  
  LowestOrder = 1;  
  //LowestOrder = 0;  
  n_XBins = 40;
  x_Order = 6;
  //whPr[111] = {	" Higgs production into bbar decay", "mcfm-H", "-Higgs-111", 						q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  //whPr[112] = {	" Higgs production into tau-taubar decay", "mcfm-H", "-Higgs-112", 					q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  
  //
  // TTbar, BBbar & CCbar production
  q2_Low = 1.0;
  q2_Up = 4000*4000;
  nq2_Bins = 15;
  q_order = 3;
  LowestOrder = 2;  
  n_XBins = 40;
  x_Order= 6;
  whPr[141] = {	" TTbar production with 2 semi-leptonic decays", "mcfm-TT", "-TTbar-141", 							q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[142] = {	" TTbar production with 2 semi-leptonic decays, corrections only in decays", "mcfm-TT", "-TTbar-142", 				q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[144] = {	" TTbar production with 2 semi-leptonic decays, no correlations", "mcfm-TT", "-TTbar-144", 					q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[145] = {	" TTbar production with 2 semi-leptonic decays, no spin correlations in top decays", "mcfm-TT", "-TTbar-145", 			q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[146] = {	" TTbar production with Tbar hadronic decay, radiative corrections in production and decay", "mcfm-TT", "-TTbar-146",     	q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[147] = {	" TTbar production with Tbar hadronic decay, radiative corrections in Tbar decay", "mcfm-TT", "-TTbar-147", 			q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[148] = {	" TTbar production with Tbar hadronic decay, radiative corrections in W decay", "mcfm-TT", "-TTbar-148", 			q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[149] = {	" TTbar production with T hadronic decay, radiative corrections in production and decay", "mcfm-TT", "-TTbar-149", 		q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[150] = {	" TTbar production with T hadronic decay, radiative corrections in T decay", "mcfm-TT", "-TTbar-nq2_Bins0", 			q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[151] = {	" TTbar production with T hadronic decay, radiative correstions in W decay", "mcfm-TT", "-TTbar-nq2_Bins1", 			q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[157] = {	" TTbar production", "mcfm-TT", "-TTbar", 											q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[158] = {	" BBbar production", "mcfm-BB", "-BBbar", 											q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[159] = {	" CCbar production", "mcfm-CC", "-CCbar", 											q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops}; 

  //
  //tau^-+tau^+
  /*
  q2_Low = 1.0;
  q2_Up = 4000*4000;
  nq2_Bins = 15;
  q_order = 3;
  LowestOrder = 0;  
  n_XBins = 40;
  x_Order= 6;
  whPr[221] = {	" +Tau-T	au production into leptons decay", "//mcfm-tautau", "-Tau+Tau-221", 						q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
*/

//
//  Photon production
  q2_Low = 1.0;
  q2_Up = 4000*4000;
  nq2_Bins = 15;
  q_order = 3;
  LowestOrder = 0;
  whPr[280] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_280",q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[281] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_281",q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[282] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_282",q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[283] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_283",q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[284] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_284",q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[285] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_285",q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  whPr[286] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_286",q2_Low, q2_Up, nq2_Bins, q_order, LowestOrder, n_XBins, x_Order, x_low, x_up, n_loops};
  
  
  return whPr[proc];
}

