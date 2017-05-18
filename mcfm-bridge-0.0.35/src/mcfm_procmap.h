
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
      int nq2bins;
      int qOrder;
      int Lowest_Order;
} Proc1;

extern struct info2{
      std::string chan;
      std::string pdf_fun;
      std::string glab;
      double q2low;
      double q2up;
      int nq2bins;
      int qOrder;
      int Lowest_Order;
      
       int nxbins;
       int xOrder;
       int ngrids;
       
       int NobsBins[3];
       double Eta[13];
       double  Pt[14];
       double Ptf[11];
} Proc2;

 info whichProcess(int proc);
 info2 whichProcess2(int proc);


#endif // MCFM_PROCMAP_H
