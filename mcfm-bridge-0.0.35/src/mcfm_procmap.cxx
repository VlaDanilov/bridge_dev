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

//#include "mcfm_procmap.h"
#include "mcfm_grid.h"

info whichProcess(int  proc)
{
  std::map<int, info> whPr;
  
  whPr[1]   = {	"W+ production", "mcfmwp.config", "-Wplus", 			6299.99, 6800.01, 3, 1, 0};
  whPr[6]   = {	"W- production", "mcfmwm.config", "-Wminus",  			6399.99, 6400.01, 3, 1, 0};
  whPr[11] = {	" W+ + jet production", "mcfm-wpjet", "-WplusJet", 		 1.0, 4000*4000, 15, 3, 0};
  whPr[16] = {	" W- + jet production", "mcfm-wmjet", "-WminusJet", 		 1.0, 4000*4000, 15, 3, 0}; 
  whPr[31] = {	" Z production", "mcfm-z", "-Z0", 						8280.99, 8281.01, 3, 1, 0};
  whPr[41] = {	" Z-jet production", "mcfm-zjet",  "-Zjet_41", 			8280.99, 8281.01, 3, 1, 0};
  whPr[42] = {	" Z-jet production", "mcfm-zjet",  "-Zjet_42", 			8280.99, 8281.01, 3, 1, 0};
  whPr[43] = {	" Z-jet production", "mcfm-zjet",  "-Zjet_43", 			8280.99, 8281.01, 3, 1, 0};  
  whPr[13] = {         " W+ + Cbar production", "mcfm-wpc", "-WplusCbar", 	1.0, 4000*4000, 15, 3, 1 };
  whPr[18] = {	"W- + C production", "mcfm-wmc", "-WminusC", 		1.0, 4000*4000, 15, 3, 1 };
  
  whPr[141] = {	" TTbar production with 2 semi-leptonic decays", "mcfm-TT", "-TTbar-141", 										1.0, 4000*4000, 15, 3, 2};
  whPr[142] = {	" TTbar production with 2 semi-leptonic decays, corrections only in decays", "mcfm-TT", "-TTbar-142", 					1.0, 4000*4000, 15, 3, 2};
  whPr[144] = {	" TTbar production with 2 semi-leptonic decays, no correlations", "mcfm-TT", "-TTbar-144", 							1.0, 4000*4000, 15, 3, 2};
  whPr[145] = {	" TTbar production with 2 semi-leptonic decays, no spin correlations in top decays", "mcfm-TT", "-TTbar-145", 			1.0, 4000*4000, 15, 3, 2};
  whPr[146] = {	" TTbar production with Tbar hadronic decay, radiative corrections in production and decay", "mcfm-TT", "-TTbar-146", 	1.0, 4000*4000, 15, 3, 2};
  whPr[147] = {	" TTbar production with Tbar hadronic decay, radiative corrections in Tbar decay", "mcfm-TT", "-TTbar-147", 			1.0, 4000*4000, 15, 3, 2};
  whPr[148] = {	" TTbar production with Tbar hadronic decay, radiative corrections in W decay", "mcfm-TT", "-TTbar-148", 				1.0, 4000*4000, 15, 3, 2};
  whPr[149] = {	" TTbar production with T hadronic decay, radiative corrections in production and decay", "mcfm-TT", "-TTbar-149", 		1.0, 4000*4000, 15, 3, 2};
  whPr[150] = {	" TTbar production with T hadronic decay, radiative corrections in T decay", "mcfm-TT", "-TTbar-150", 					1.0, 4000*4000, 15, 3, 2};
  whPr[151] = {	" TTbar production with T hadronic decay, radiative correstions in W decay", "mcfm-TT", "-TTbar-151", 					1.0, 4000*4000, 15, 3, 2};
  whPr[157] = {	" TTbar production", "mcfm-TT", "-TTbar", 																	1.0, 4000*4000, 15, 3, 2};
  whPr[158] = {	" BBbar production", "mcfm-BB", "-BBbar", 																	1.0, 4000*4000, 15, 3, 2};
  whPr[159] = {	" CCbar production", "mcfm-CC", "-CCbar", 																	1.0, 4000*4000, 15, 3, 2}; 
  return whPr[proc];
}
  
  
info2 whichProcess2(int proc)
{
  std::map<int, info2> whPr2;
  whPr2[280] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_280", 1.0, 4000*4000, 15, 3, 0, 30, 6, 3, {12, 13, 10},  { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 }};
  whPr2[281] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_281", 1.0, 4000*4000, 15, 3, 0, 30, 6, 3, {12, 13, 10},  { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 }};
  whPr2[282] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_282", 1.0, 4000*4000, 15, 3, 0, 30, 6, 3, {12, 13, 10},  { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 }};
  whPr2[283] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_283", 1.0, 4000*4000, 15, 3, 0, 30, 6, 3, {12, 13, 10},  { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 }};
  whPr2[284] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_284", 1.0, 4000*4000, 15, 3, 0, 30, 6, 3, {12, 13, 10},  { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 }};
  whPr2[285] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_285", 1.0, 4000*4000, 15, 3, 0, 30, 6, 3, {12, 13, 10},  { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 }};
  whPr2[286] = {	" Photon production", "photonLO.config:photonNLO.config", "-GammaProd_286", 1.0, 4000*4000, 15, 3, 0, 30, 6, 3, {12, 13, 10},  { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 },  { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 }};
  return whPr2[proc];
}

