// emacs: this is -*- c++ -*-
//
//   mcfmzjet_pdf.h        
//
//   pdf transform functions                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: mcfmz_pdf.h, v1.0   Mon Dec 10 01:36:04 GMT 2007 sutt $


#ifndef __MCFMZJET_PDF_H
#define __MCFMZJET_PDF_H

#include <cmath>

#include "appl_grid/appl_pdf.h" 

#include "TFile.h"
#include "TVectorT.h"
#include "TMatrixT.h"

//
// MCFM Z-jet production
//
class mcfmzjet_pdf : public appl::appl_pdf { 

public:
	mcfmzjet_pdf() : appl_pdf("mcfm-zjet") { m_Nproc = 33; } 
	
	~mcfmzjet_pdf() { } 
	
	virtual void evaluate(const double* fA, const double* fB, double* H) const;
	
};


// actual funtion to evaluate the pdf combinations 

inline  void mcfmzjet_pdf::evaluate(const double* fA, const double* fB, double* H) const {
	//  const int nQuark = 6;
	const int iQuark = 5; 
	
	// offset psd ptrs so can use [-6..6] indexing rather than [nQuark+..] indexing
	fA += 6;
	fB += 6;
	
	double GA=fA[0];
	double GB=fB[0];
	
	double UpA=0; double UpB=0; double DnA=0; double DnB=0;
	double UpbarA=0; double UpbarB=0; double DnbarA=0; double DnbarB=0;
	
	// zero H first
	for(int i = 0; i < m_Nproc; i++)  H[i] = 0;
	
	for (int i = 1; i <= iQuark; i++) {
	  if ((i % 2) == 0) {
	    UpA += fA[i];
	    UpB += fB[i];
	  }else{
	    DnA += fA[i];
	    DnB += fB[i];
	  }
	}
	
	for(int i = -iQuark; i < 0; i++) {
	  if (((int)(std::abs((double)i)) % 2) == 0) {
	    UpbarA += fA[i];
	    UpbarB += fB[i];
	  } else {
	    DnbarA += fA[i];
	    DnbarB += fB[i];
	  }
	}
	
	// zero H first
	for(int i = 0; i<m_Nproc;++i) H[i]=0;
	
	static int  _choice[13]  = { 2, 3, 2, 3, 2, 3,    0,    1, 0, 1, 0, 1, 0 } ;
	static int*  choice      = _choice+6;
	// second diagonal (QbQb and QQ) for real contribution
	static int  _choiceR[13] = { 14, 13, 14, 13, 14, 13,    0,    15, 16, 15, 16, 15, 16 } ;
	static int*  choiceR     = _choiceR+6;
	
	for(int i = -iQuark; i <= iQuark; i++) {
	  if (i == 0) continue;
	  //    int Choice = ( i>0 ? i%2 : abs(i%2)+2 );
	  int Choice = choice[i];
	  H[Choice] += fA[i] * fB[-i];
	  int ChoiceR = choiceR[i];
	  H[ChoiceR] += fA[i] * fB[i];
	}
	// H0 UUbar  H1 DDbar H2 UbarU H3 DbarD
	// H14 UbarUbar  H13 DbarDbar H16 UU H5 DD
	H[4]  = GA * UpB;
	H[5]  = GA * UpbarB;
	H[6]  = GA * DnB;
	H[7]  = GA * DnbarB;
	H[8]  =    UpA *GB;
	H[9]  = UpbarA *GB;
	H[10] =    DnA *GB;
	H[11] = DnbarA *GB;
	
	H[12] = GA*GB;
	// V1 contribution 
//	H[13] = (UpA + UpbarA + DnA + DnbarA)*( UpbarB ) - H[0]  ;
//	H[14] = (UpA + UpbarA + DnA + DnbarA)*( DnbarB ) - H[1]  ;
//	H[15] = (UpA + UpbarA + DnA + DnbarA)*( UpB )    - H[2]  ;
//	H[16] = (UpA + UpbarA + DnA + DnbarA)*( DnB )    - H[3]  ;
        // V2 contribution
//	H[17] = (UpB + UpbarB + DnB + DnbarB)*( UpbarA ) - H[2]  ;
//	H[18] = (UpB + UpbarB + DnB + DnbarB)*( DnbarA ) - H[3]  ;
//	H[19] = (UpB + UpbarB + DnB + DnbarB)*( UpA )    - H[0]  ;
//	H[20] = (UpB + UpbarB + DnB + DnbarB)*( DnA )    - H[1]  ;

	// Real contributions
	//ND1
//	H[25] = H[14] - H[21];
//	H[26] = H[13] - H[22];
//	H[27] = (UpA + DnA)*DnB - H[23];
//	H[28] = (UpA + DnA)*UpB - H[24];
	//ND 4
//	H[29] = H[20] - H[23];
//	H[30] = H[19] - H[24];
//	H[31] = (UpbarB+DnbarB)*DnbarA - H[21];
//	H[32] = (UpbarB+DnbarB)*UpbarA - H[22];
	// ND 3 
//	H[33] = (UpbarA + DnbarA)*(DnB) - H[3];
//	H[34] = (UpbarA + DnbarA)*(UpB) - H[2];
	// ND 2
//	H[35] = (UpB + DnB)*(DnbarA) - H[3];
//	H[36] = (UpB + DnB)*(UpbarA) - H[2];

	// ND0
	//
	H[17] = DnA*DnB - H[15];
	H[18] = UpA*DnB;
	H[19] = DnA*UpB;
	H[20] = UpA*UpB - H[16];
	//
	H[21] = DnbarA*DnB - H[3];
	H[22] = UpbarA*DnB;
	H[23] = DnbarA*UpB;
	H[24] = UpbarA*UpB - H[2];
	//
	H[25] = DnA*DnbarB - H[1];
	H[26] = UpA*DnbarB;
	H[27] = DnA*UpbarB;
	H[28] = UpA*UpbarB - H[0];
	//
	H[29] = DnbarA*DnbarB - H[13];
	H[30] = UpbarA*DnbarB;
	H[31] = DnbarA*UpbarB;
	H[32] = UpbarA*UpbarB - H[14];

	return;
}


extern "C" void fmcfmzjet_pdf__(const double* fA, const double* fB, double* H);

#endif  // __MCFMZJET_PDF_H

