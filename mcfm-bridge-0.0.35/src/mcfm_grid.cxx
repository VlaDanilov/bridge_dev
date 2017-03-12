//
//   @file    mcfm_grid.cxx         
//             
//
//   @author M.Sutton
// 
//   Copyright (C) 2011 M.Sutton (sutt@cern.ch)    
//
//   $Id: mcfm_grid.cxx, v0.0   Fri 12 Aug 2011 08:03:37 CEST sutt $

#include <cstdlib>


#include "mcfm_grid.h"
#include "appl_grid/lumi_pdf.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

static const bool _debug_ = false;

/// special mcfm specific methods

void appl::mcfm_grid::fillMCFM(double obs)
{
  //  bool _debug_ = true;

  if (_debug_) std::cout << __PRETTY_FUNCTION__ << "   --------------------------------------------- " << std::endl;
  
  //  double* weight = new double[ subProcesses() ];

  std::vector<double> _weightLO( subProcesses(0),0 );
  double* weightLO = &_weightLO[0];

  std::vector<double> _weightNLO( subProcesses(1),0 );
  double* weightNLO = &_weightNLO[0];


  
  double  scale2 =  gridevent_.ag_scale * gridevent_.ag_scale;
  double _x1 =  gridevent_.ag_xx1, _x2 =  gridevent_.ag_xx2;
  int flag = 0;
  
  if (_debug_) 
     std::cout <<" x1 = "<< _x1
	       <<" x2 = "<< _x2
	       <<" Q  = "<< gridevent_.ag_scale
	       <<" CON = "<< gridevent_.contrib
	       <<std::endl;
  

  if      ( gridevent_.contrib == 100 ) {
    // BORN
    flag = 0;
    collectWeight( gridevent_.contrib , flag, weightLO, 0 );
    fill_grid( _x1, _x2, scale2, obs, weightLO, 0 );
    
  }
  else if ( gridevent_.contrib == 200 ) { 
    //REAL
    collectWeight( gridevent_.contrib, gridevent_.dipole, weightNLO, 1 );
    fill_grid( _x1, _x2, scale2, obs, weightNLO, 1 );
  }
  else if ( gridevent_.contrib == 300 ) { 
    //VIRTUAL   
    // VIRT BORN
    flag = 0;
    collectWeight( gridevent_.contrib, flag, weightLO, 0 );
    fill_grid( _x1, _x2, scale2, obs, weightLO, 0 );

    //VIRT X1X2
    flag = -1;
    collectWeight( gridevent_.contrib, flag, weightNLO, 1 );
    fill_grid( _x1, _x2, scale2, obs, weightNLO, 1 );

    
    //VIRT X1onZ
    flag = -2;
    _x1 =  gridevent_.ag_x1z; 
    _x2 =  gridevent_.ag_xx2;
    collectWeight( gridevent_.contrib, flag, weightNLO, 1 );
    fill_grid( _x1, _x2, scale2, obs, weightNLO, 1 );
    
    //VIRT X2onZ
    flag = -3;
    _x1 =  gridevent_.ag_xx1; 
    _x2 =  gridevent_.ag_x2z;
    collectWeight( gridevent_.contrib, flag, weightNLO, 1 );
    fill_grid( _x1, _x2, scale2, obs, weightNLO, 1 );

  }
  else  { 
    // UNKNOWN
    std::cout <<__PRETTY_FUNCTION__<<" Unknown contribution!!!!!!!!!"<< std::endl; 
    std::cout <<__PRETTY_FUNCTION__<<" Unknown contribution!!!!!!!!!"<< std::endl; 
    std::cout <<__PRETTY_FUNCTION__<<" Unknown contribution!!!!!!!!!"<< std::endl; 
  }

  //  double binWidth = deltaobs(obsbin(obs));
  getReference()->Fill( obs, gridevent_.refwt );

  if (false)
    std::cout << "   --------------------------------------------- " << std::endl;
  

  //  std::cout << "delete" << std::endl;
  //  delete [] weight;

  return;
}


void appl::mcfm_grid::collectWeight(const int contribution, const int id, double* wt, int iorder )
{
  double factor = 1.0;
  int iproc = -1;

  ///  use the order dependent subProcess(iorder) multiplicity
  ///  rather than  just the number of sub processes from the lowest 
  ///  order, which may not be the same
  for ( int i=0 ; i<subProcesses(iorder); i++ ) wt[i] =  0.0;

  //   std::cout<<" \n( ";
  //   for (int jj =  0 ; jj < subProcesses();jj++) std::cout << wt[jj] <<" , ";
  //   std::cout<<")\n";
  
  double _count = 0, _sumfactor = 0;

  for ( int iflav = -nflav_.nflav; iflav <= nflav_.nflav; iflav++)
    for ( int jflav = -nflav_.nflav; jflav <= nflav_.nflav; jflav++)
      {

	decideSubProcess( iflav, jflav, iproc, factor, iorder );

	//	std::cout << "." << "\n";

	double _born_ = gridweight_.weightb [ __nf__ + jflav ][ __nf__ + iflav ];
	double _virt_ = gridweight_.weightv [ __nf__ + jflav ][ __nf__ + iflav ];
	double _vir1_ = gridweight_.weightv1[ __nf__ + jflav ][ __nf__ + iflav ];
	double _vir2_ = gridweight_.weightv2[ __nf__ + jflav ][ __nf__ + iflav ];
	if(false)
	  {
	    if ( _born_!=0 || _virt_!=0 || _vir1_!=0 || _vir2_!=0 )
	      {
		std::cout<<"( "<<iflav<<" , "<<jflav<<" ) Weights :  born = "<<_born_
			 << " _v_ = "<<_virt_<<" _v1_ = "<<_vir1_<<" _v2_ "<<_vir2_ 
			 <<" \n "
			 <<std::endl;
	      }
	  }

	if ( iproc < 0 ) continue;

        _count+=1; _sumfactor += 1;
	
	if ( contribution == 100 ) 
	  wt[ iproc ] +=  factor * gridweight_.weightb[ __nf__ + jflav ][ __nf__ + iflav ];
	else if ( contribution == 200 ) 
	  wt[ iproc ] +=  factor * gridweight_.weightr[ __nf__ + jflav ][ __nf__ + iflav ][ id ];
	else if ( contribution == 300 )
	  {
	    if      ( id ==  0 )
	      wt[ iproc ] +=  factor * gridweight_.weightb [ __nf__ + jflav ][ __nf__ + iflav ];
	    else if ( id == -1 )
	      wt[ iproc ] +=  factor * gridweight_.weightv [ __nf__ + jflav ][ __nf__ + iflav ];
	    else if ( id == -2 )
	      wt[ iproc ] +=  factor * gridweight_.weightv1[ __nf__ + jflav ][ __nf__ + iflav ];
	    else if ( id == -3 )
	      wt[ iproc ] +=  factor * gridweight_.weightv2[ __nf__ + jflav ][ __nf__ + iflav ];
	    else
	      {
		std::cout <<"\t\t : "
			  <<__PRETTY_FUNCTION__<< " : \t Unknown id for virtual process "
			  << std::endl;
		exit(-1);	
	      }
	  }
	else
	  {
	    std::cout <<"\t\t : "
		      <<__PRETTY_FUNCTION__<< " : \t Unknown contribution "
		      << std::endl;
	    exit(-1);
	  }
	
	if ( false )
	  {
	    //if (0 != gridweight_.weightb[ __nf__ + jflav ][ __nf__ + iflav ]) 
	    std::cout <<" ( i= "<<iflav<<" j= "<<jflav<<" : "
		      << gridweight_.weightb[ __nf__ + jflav ][ __nf__ + iflav ]<<" ) proc = "<<iproc;
	    
	    //std::cout<<" \n( ";
	    //for (int jj =  0 ; jj < subProcesses(order);jj++) std::cout << wt[jj] <<" , ";
	    //std::cout<<")\n";
	    
	  }
      }

  for (int jj=0 ; jj<subProcesses(iorder) ; jj++ )  wt[jj] *= gridweight_.weightfactor ;

  //std::cout<<" \n\t\t factor = "<< gridweight_.weightfactor<<std::endl;
  if (_debug_)  
    { 
      std::cout << "\t\t collectWeight :  PROC = " << nproc_.nproc <<" W = ( dip = "<<id<<" , ";
      for (int jj =  0 ; jj < subProcesses(iorder);jj++) std::cout << wt[jj] <<" , ";
      std::cout <<" ) ifact = "<<_sumfactor<<"/"<<_count<<" psw = "<< gridweight_.weightfactor <<" refwt = "<<gridevent_.refwt  << std::endl;
    }
  //   std::cout <<" me(2,-1,0)"<< gridweight_.weightr[ __nf__ - 1][ __nf__ + 2][0] <<std::endl;
  //   std::cout <<" me(2,-1,1)"<< gridweight_.weightr[ __nf__ - 1][ __nf__ + 2][1] <<std::endl;
  
  return ;
}


void appl::mcfm_grid::decideSubProcess(const int iflav1, const int iflav2, int& iProcess , double& factor, int iorder )
{

  iProcess = -1;
  factor   =  0.;
  
  if ( getGenpdf().find("basic.config") != std::string::npos ) { 
    
    lumi_pdf* pdf = dynamic_cast<lumi_pdf*>(genpdf(iorder));
    
    if ( pdf==0 ) { 
      std::cerr << __PRETTY_FUNCTION__ << " cannot retrieve lumi_pdf" << std::endl; 
    }

    //    static std::string duff = pdf->name();
    //    if ( pdf->name()!=duff ) 
    //    std::cout << "iorder " << iorder << "\tpdf " << pdf << " " << pdf->name() << "  " << pdf->Nproc() << std::endl; 

    iProcess = pdf->decideSubProcess( iflav1, iflav2 );
    if ( iProcess != -1 ) { 
      factor = (*pdf)[iProcess].size(); 
      if ( factor>0 ) factor = 1.0/factor; // hmmm
    }

    //    std::cout << "iorder " << iorder << "\tiproc " << iProcess << "\tfactor " << factor << std::endl; 
    //    duff = pdf->name();

  }
  else if( getGenpdf().find("basic") != std::string::npos )  {
    
    if ( std::fabs(iflav1)<6 && std::fabs(iflav2)<6 ) {  
      iProcess = 11*(iflav1+5)+(iflav2+5);
      factor = 1.0;
    }

  }
  else if ( (nproc_.nproc == 1)  ||  (nproc_.nproc == 11 ))  {
  
    if( (iflav1 == 0) && (iflav2 ==  2) )  {
      factor = 1.0/ckm_.vsum[ __nf__ + iflav2];
      iProcess = 5;
    }        
    else if( (iflav1 == 0) && (iflav2 == -1) )	{
      factor = 1.0/ckm_.vsum[ __nf__ + iflav2];
      iProcess = 4;
    }
    else if( (iflav2 == 0) && (iflav1 == 2) )  {
      factor = 1.0/ckm_.vsum[ __nf__ + iflav1];
      iProcess = 3;
    }
    else if( (iflav2 == 0) && (iflav1 == -1) )  {
      factor = 1.0/ckm_.vsum[ __nf__ + iflav1];
      iProcess = 2;
    }
    else if( (iflav1 == 2) && (iflav2 == -1) )  {
      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ + iflav1];
      iProcess = 1;
    }
    else if( (iflav1 == -1) && (iflav2 == 2) )  {
      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ + iflav1];
      iProcess = 0;
    }
    else if ( (nproc_.nproc == 11) && ((iflav1 == 0) && (iflav2 == 0)) )
      iProcess = 6;

#if 0
    if ( iProcess!=-1 ) { 
      std::cout << "\tiflav1 " << iflav1 << "\tiflav2 " << iflav2 
		<< "\tip "  << iProcess 
		<< "\tfac " << factor 
		<< "\tckmsum " <<  1/ckm_.vsum[ __nf__ + iflav2]
		<< "\tckmsum " <<  1/ckm_.vsum[ __nf__ + iflav2]
		<< "\tckm2 "   <<  1/ckm_.vsq[ __nf__ + iflav2][ __nf__ + iflav1] 
		<<  std::endl;
    }
#endif

  }
  else if ( (nproc_.nproc == 6) ||  (nproc_.nproc == 16) ) {
    
    if((iflav1 == 0) && (iflav2 == -2))  {
      factor = 1.0/ckm_.vsum[__nf__ + iflav2];
      iProcess = 5;
    }
    else if((iflav1 == 0) && (iflav2 == 1))  {
      factor = 1.0/ckm_.vsum[__nf__ + iflav2];
      iProcess = 4;
    }
    else if((iflav2 == 0) && (iflav1 == -2)) { 
      factor = 1.0/ckm_.vsum[__nf__ + iflav1];
      iProcess = 3;
    }
    else if((iflav2 == 0) && (iflav1 == 1)) { 
      factor = 1.0/ckm_.vsum[__nf__ + iflav1];
      iProcess = 2;
    }
    else if ((iflav1 == -2) && (iflav2 == 1)) { 
      factor = 1.0/ckm_.vsq[__nf__ + iflav2][__nf__ + iflav1];
      iProcess = 1;
    }
    else if ((iflav1 == 1) && (iflav2 == -2)) { 
      factor = 1.0/ckm_.vsq[__nf__ + iflav1][ __nf__ + iflav2];
      iProcess = 0;
    }
    else if ( (nproc_.nproc == 16) && ((iflav1 == 0) && (iflav2 == 0)) )
      iProcess = 6;
  }
  else if ( (nproc_.nproc == 31) || ( (nproc_.nproc >= 41) && (nproc_.nproc <= 43) ))    {
    if ( iflav2 == 0 )	{
      if      (iflav1 == -1) {iProcess = 11;}
      else if (iflav1 ==  1) {iProcess = 10;}
      else if (iflav1 == -2) {iProcess = 9;}
      else if (iflav1 ==  2) {iProcess = 8;}
    }
    else if ( iflav1 == 0 )	{
      if      (iflav2 == -1) iProcess = 7;
      else if (iflav2 ==  1) iProcess = 6;
      else if (iflav2 == -2) iProcess = 5;
      else if (iflav2 ==  2) iProcess = 4;
    }
    else if ( (iflav1 != 0 ) && ( iflav2 != 0 ) )  {
      if      (iflav1 == -1) iProcess = 3;
      else if (iflav1 == -2) iProcess = 2;
      else if (iflav1 ==  1) iProcess = 1;
      else if (iflav1 ==  2) iProcess = 0;
    }
    else if ( (nproc_.nproc != 31) && ( (iflav1 == 0 ) && ( iflav2 == 0 ) ) )
      iProcess = 12;

    factor = 1.0;
  } 
  else if ( 
	   (nproc_.nproc == 157) || (nproc_.nproc == 158) || (nproc_.nproc==159) || 
	   ((nproc_.nproc>=141)&&(nproc_.nproc<=151)&&(nproc_.nproc!=143)&&(nproc_.nproc!=145)) )
    {
      //factor = std::pow(0.5*std::pow(FourPi, 2), 2)/std::pow(qcdcouple_.gsq, 2); 
      factor = 1.0;
      if      ( (iflav1 ==  0 ) && (iflav2 == 0)  ) iProcess = 0;
      else if ( (iflav1 ==  1 ) && (iflav2 ==  0) ) iProcess = 1;
      else if ( (iflav1 ==  0 ) && (iflav2 ==  1) ) iProcess = 2;
      else if ( (iflav1 == -1 ) && (iflav2 ==  0) ) iProcess = 3;
      else if ( (iflav1 ==  0 ) && (iflav2 == -1) ) iProcess = 4;
      else if ( (iflav1 ==  1 ) && (iflav2 == -1) ) iProcess = 5;
      else if ( (iflav1 == -1 ) && (iflav2 ==  1) ) iProcess = 6;
      else {factor = 0.0;iProcess=-1;};

      //      std::cout << "\t\t *** \t"<<iflav1<<" <> "<<iflav2<<" iproc = "<<iProcess<<" factor = "<<factor<< std::endl; 
    }
  else if ( (nproc_.nproc == 13) )
    {
      if ( iflav1*iflav2 == 0 )
	{
	  if ( (iflav1 ==  0 ) && ( iflav2 == -1))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ + 4];
	      iProcess = 0;
	    }
	  else if ( (iflav1 ==  -1 ) && ( iflav2 == 0) )
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + 4][ __nf__ + iflav1];
	      iProcess = 1;
	    }
	  else if ( (iflav1 ==  0 ) && ( iflav2 == 0) )
	    {
	      factor = 1.0;
	      iProcess = 2;
	    }
	}
      else if ( iflav1*iflav2 < 0 )
	{
	  if ( (iflav1 ==  -1) && (iflav2 == 3) )
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + 4][ __nf__ + iflav1];
	      iProcess = 3;
	    }
	  else if ( (iflav1 ==  3 ) && ( iflav2 == -1))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ + 4];
	      iProcess = 4;
	    }
	} 
      else if (iflav1*iflav2 > 0 )
	{
	  if ( (iflav1 ==  -2 ) && ( iflav2 == -1))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ + 4];
	      iProcess = 5;
	    }
	  else if ( (iflav1 ==  -1 ) && ( iflav2 == -2))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + 4][ __nf__ + iflav1];
	      iProcess = 6;
	    }
	  else if ( (iflav1 ==  -1 ) && ( iflav2 == -1))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ + 4];
	      iProcess = 7;
	    }
	  else if ( (iflav1 ==  -1 ) && ( iflav2 == -3))
	    {
	      factor = 1.0;
	      iProcess = 8;
	    }
	  else if ( (iflav1 ==  -3 ) && ( iflav2 == -1))
	    {
	      factor = 1.0;
	      iProcess = 9;
	    }
	}
       
    }
  else if ( (nproc_.nproc == 18) )
    {
      if ( iflav1*iflav2 == 0 )
	{
	  if ( (iflav1 ==  0 ) && ( iflav2 == 1))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ - 4];
	      iProcess = 0;
	    }
	  else if ( (iflav1 ==  1 ) && ( iflav2 == 0) )
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ - 4][ __nf__ + iflav1];
	      iProcess = 1;
	    }
	  else if ( (iflav1 ==  0 ) && ( iflav2 == 0))
	    {
	      factor = 1.0;
	      iProcess = 2;
	    }
	}
      else if ( iflav1*iflav2 < 0 )
	{
	  if ( (iflav1 ==  1) && (iflav2 == -3) )
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ - 4][ __nf__ + iflav1];
	      iProcess = 3;
	    }
	  else if ( (iflav1 ==  -3 ) && ( iflav2 == 1))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ - 4];
	      iProcess = 4;
	    }
	}
      else if (iflav1*iflav2 > 0 )
	{
	  if ( (iflav1 ==  2 ) && ( iflav2 == 1))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ - 4];
	      iProcess = 5;
	    }
	  else if ( (iflav1 ==  1 ) && ( iflav2 == 2))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ - 4][ __nf__ + iflav1];
	      iProcess = 6;
	    }
	  else if ( (iflav1 ==  1 ) && ( iflav2 == 1))
	    {
	      factor = 1.0/ckm_.vsq[ __nf__ + iflav2][ __nf__ - 4];
	      iProcess = 7;
	    }
	  else if ( (iflav1 ==  1 ) && ( iflav2 == 3))
	    {
	      factor = 1.0;//ckm_.vsq[ __nf__ + iflav2][ __nf__ - 4];
	      iProcess = 8;
	    }
	  else if ( (iflav1 ==  3 ) && ( iflav2 == 1))
	    {
	      factor = 1.0;//ckm_.vsq[ __nf__ + iflav2][ __nf__ - 4];
	      iProcess = 9;
	    }
	}
      // cout << " fl1 = "<<iflav1<<" fl2 = "<<iflav2 
      //  	   << "   factor = "<<factor<<" iProc = "<<iProcess
      // 	   <<endl;

       
    }
  // else
  //   {
  //     iProcess = -1;
  //     factor   = 0.;
  //   }
  //  std::cout << "\t\t *** \t"<<iflav1<<" <> "<<iflav2<<" iproc = "<<iProcess<< std::endl; 
  
  return;
}

