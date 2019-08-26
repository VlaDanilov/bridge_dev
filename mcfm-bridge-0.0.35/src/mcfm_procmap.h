#ifndef MCFM_PROCMAP_H
#define MCFM_PROCMAP_H

#include <iostream>
#include <string>
#include <map>


struct info{
      std::string chan;
      std::string pdf_fun;
      std::string glab;
      double q2low;
      double q2up;
      int nq2bins;
      unsigned int qOrder;
      unsigned int LowestOrder;      
      int nxbins;
      unsigned int xOrder;
      double xlow;
      double xup;
      unsigned int nloops;
};
 
std::map<int,info> procmap =
{
  // W production
  {1, { " W+ production", "mcfmwp.config", "-Wplus",  6299.99, 6800.01, 3, 1, 0, 40, 6, 1.0e-9, 1.0, 1}},
  {6, { " W- production", "mcfmwm.config", "-Wminus", 6299.99, 6800.01, 3, 1, 0, 40, 6, 1.0e-9, 1.0, 1}},
  
  // W+jet
  { 11, { " W+ + jet production", "mcfm-wpjet", "-WplusJet",  1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 13, { " W+ + Cbar production", "mcfm-wpc", "-WplusCbar",  1.0, 1.6e7, 15, 3, 1, 40, 6, 1.0e-9, 1.0, 1}},
  { 16, { " W- + jet production", "mcfm-wmjet", "-WminusJet", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}}, 
  { 18, { " W- + C production", "mcfm-wmc", "-WminusC",       1.0, 1.6e7, 15, 3, 1, 40, 6, 1.0e-9, 1.0, 1}},
  
  // Z production
  { 31, { " Z production", "mcfm-z", "-Z0",                   8280.99, 8281.01, 3, 1, 0, 40, 6, 1.0e-9, 1.0, 1}},
  
  // Z+jet
  { 41, { " Z-jet production", "mcfm-zjet",  "-Zjet_41", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 42, { " Z-jet production", "mcfm-zjet",  "-Zjet_42", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 43, { " Z-jet production", "mcfm-zjet",  "-Zjet_43", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}}, 
  
  // Higgs
  { 111, { " Higgs production into bbar decay", "mcfm-H", "-Higgs-111",       1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 112, { " Higgs production into tau-taubar decay", "mcfm-H", "-Higgs-112", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
           
  // TTbar     
  { 141, { " TTbar production with 2 semi-leptonic decays", "mcfm-TT", "-TTbar-141", 					        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 142, { " TTbar production with 2 semi-leptonic decays, corrections only in decays", "mcfm-TT", "-TTbar-142", 	        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 144, { " TTbar production with 2 semi-leptonic decays, no correlations", "mcfm-TT", "-TTbar-144", 			        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 145, { " TTbar production with 2 semi-leptonic decays, no spin correlations in top decays", "mcfm-TT", "-TTbar-145",        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 146, { " TTbar production with Tbar hadronic decay, radiative corrections in production and decay", "mcfm-TT", "-TTbar-146",1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 147, { " TTbar production with Tbar hadronic decay, radiative corrections in Tbar decay", "mcfm-TT", "-TTbar-147", 	        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 148, { " TTbar production with Tbar hadronic decay, radiative corrections in W decay", "mcfm-TT", "-TTbar-148", 	        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 149, { " TTbar production with T hadronic decay, radiative corrections in production and decay", "mcfm-TT", "-TTbar-149",   1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 150, { " TTbar production with T hadronic decay, radiative corrections in T decay", "mcfm-TT", "-TTbar-nq2_Bins0", 	        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 151, { " TTbar production with T hadronic decay, radiative correstions in W decay", "mcfm-TT", "-TTbar-nq2_Bins1", 	        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 157, { " TTbar production", "mcfm-TT", "-TTbar", 									        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 158, { " BBbar production", "mcfm-BB", "-BBbar", 									        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}},
  { 159, { " CCbar production", "mcfm-CC", "-CCbar", 									        1.0, 1.6e7, 15, 3, 2, 40, 6, 1.0e-9, 1.0, 1}}, 

  // Tau-tau
  { 221, { " +Tau-Tau production into leptons decay", "mcfm-tautau", "-Tau+Tau-221", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
 
  // Photon production
  { 280, { " Photon production", "photonLO.config:photonNLO.config", "-GammaProd_280", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 281, { " Photon production", "photonLO.config:photonNLO.config", "-GammaProd_281", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 282, { " Photon production", "photonLO.config:photonNLO.config", "-GammaProd_282", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 283, { " Photon production", "photonLO.config:photonNLO.config", "-GammaProd_283", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 284, { " Photon production", "photonLO.config:photonNLO.config", "-GammaProd_284", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 285, { " Photon production", "photonLO.config:photonNLO.config", "-GammaProd_285", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}},
  { 286, { " Photon production", "photonLO.config:photonNLO.config", "-GammaProd_286", 1.0, 1.6e7, 15, 3, 0, 40, 6, 1.0e-9, 1.0, 1}}
};  


#endif // MCFM_PROCMAP_H
