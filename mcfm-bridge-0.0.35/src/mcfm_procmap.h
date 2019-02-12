
#ifndef MCFM_PROCMAP_H
#define MCFM_PROCMAP_H

#include <iostream>
#include <string>

#include "appl_grid/appl_grid.h"


 extern struct info{
      std::string chan;
      std::string pdf_fun;
      std::string glab;
      double q2low;
      double q2up;
      unsigned int nq2bins;
      unsigned int qOrder;
      unsigned int LowestOrder;      
      unsigned int nxbins;
      unsigned int xOrder;
      double xlow;
      double xup;
      unsigned int nloops;
} Proc1;
 
info whichProcess(int proc);


#endif // MCFM_PROCMAP_H
