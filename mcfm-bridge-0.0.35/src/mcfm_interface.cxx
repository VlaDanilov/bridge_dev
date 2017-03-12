//
//   mcfm_interface.cxx        
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)    
//
//   $Id: mcfm_interfce.cxx, v   Fri  8 Nov 2013 09:07:01 CET sutt


#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <sys/stat.h>

#include <cstdlib> 
#include <sys/time.h> 


#include "TFile.h"
#include "TH1D.h"
// #include "TMatrixT.h"
#include "TVectorT.h"
#include "TString.h"

#include "mcfm_grid.h"

  
// extern "C" struct {
//  bool creategrid;
//  int nSubProcess;
// } grid_;


// bool file_exists(const std::string& s) {   
//
//  if ( FILE* testfile=fopen(s.c_str(),"r") ) { 
//    fclose(testfile);
//    return true;
//  }
//  else return false;
// }




bool file_exists( const std::string& filename ) { 
  struct stat sb;
  if ( stat( filename.c_str(), &sb)==0 ) return true; // && S_ISREG(sb.st_mode ))
  else return false;
}


/// mcfm changed the the mxpart size between 6.7 and 6.8
// extern "C" int getmxpart_(void);

static const int mxpart = 14;    // mcfm parameter : max number of partons in event record. defined in Inc/constants.f

/// was 12 for mcfm 6.7
/// static const int mxpart = 12;    // mcfm parameter : max number of partons in event record. defined in Inc/constants.f

static const int _Ngrids = 6;
static       int  Ngrids = 6;
appl::mcfm_grid* mygrid[_Ngrids];

static const char* gridFiles[_Ngrids] = {
    "_eta3.root",
    "_pt3-etalt137.root",
    //    "_eta4.root",
    "_pt3-eta152-237.root",
    "_pt4.root",
    "_eta5.root",
    "_pt5.root"
};


static double Observable[_Ngrids] = {  0,  0,  0,  0,  0,  0 };   // observable array
int             nObsBins[_Ngrids] = { 40, 45, 40, 45, 40, 45 }; // eta4, pt4 cental eta-bin, pt4 forward eta-bin

static const double eta[41] =  {
  -4.0, -3.8, -3.6, -3.4, -3.2, 
  -3.0, -2.8, -2.6, -2.4, -2.2, 
  -2.0, -1.8, -1.6, -1.4, -1.2, 
  -1.0, -0.8, -0.6, -0.4, -0.2,  0.0, 
   0.2,  0.4,  0.6,  0.8,  1.0,  
   1.2,  1.4,  1.6,  1.8,  2.0,  
   2.2,  2.4,  2.6,  2.8,  3.0,  
   3.2,  3.4,  3.6,  3.8,  4.0  
};

// int nObsBins[_Ngrids] = {13, 45, 40, 45, 40, 45};      // eta4, pt4 cental eta-bin, pt4 forward eta-bin
// static const double eta[] =   {
//     -4.0,  -3.4, -2.8,  -2.2,  -1.6, -1.0, -0.4,  0.2,   0.8,  
//      1.4,  2.0,  2.6,  3.2,   3.8
// };


static const double pt[46] =  { 
      0.0,   5.0,   10.0,   15.0,   20.0, 
     25.0,  30.0,   35.0,   40.0,   45.0, 
     50.0,  55.0,   60.0,   65.0,   70.0, 
     75.0,  80.0,   85.0,   90.0,   95.0, 
    100.0, 105.0,  110.0,  115.0,  120.0, 
    125.0, 130.0,  135.0,  140.0,  145.0,
    150.0, 175.0,  200.0,  225.0,  250.0,
    275.0, 300.0,  350.0,  400.0,  450.0,
    500.0, 750.0, 1000.0, 1500.0, 2000.0,
    2500.0
};


long unsigned int runs  =  0;
bool isBooked           =  false;
std::string glabel      =  "";

void getObservable( const double evt[][mxpart] );

int  cuts(int);


std::string date() { 
  time_t _t;
  time(&_t);
  return ctime(&_t);
}

void book_grid()  // inital grid booking
{
  if (isBooked) return;
  
  //  time_t _t;
  //  time(&_t);
  
  std::cout<<" ***********************************************"<<std::endl;
  std::cout<<" booking the grids " << date() << std::endl;
  
  // binning information for the grid constructor
  double xLow    = 1.0e-9, xUp = 1.0;
  int    nXbins  = 40;
  int    xorder  = 6;
  double q2Low   = 1.0*1.0;
  double q2Up    = 4000*4000;
  int    nQ2bins = 15;
  int    qorder  = 3;
  // set transform2 value
  double apramval=5.;
  appl::grid::transformvar(apramval);

  // lowest order in alphas	
  int lowest_order = 0;
  // how many loops
  int nloops = 1;
  
  // number of observables and binning for observables  
  const double *obsBins[_Ngrids] = { eta, pt, eta, pt, eta, pt };
    
  std::string pdf_function;

  glabel = "grid-40-6-15-3";

  const char* basename = std::getenv("appl_basename");
  if ( basename && std::string(basename)!="" ) glabel = basename;

  const char* q2upper = std::getenv("appl_q2up");
  if ( q2upper && std::string(q2upper)!="" ) q2Up = std::atof(q2upper);

  const char* q2lower = std::getenv("appl_q2low");
  if ( q2lower && std::string(q2lower)!="" ) q2Low = std::atof(q2lower);

  const char* q2order = std::getenv("appl_q2order");
  if ( q2order && std::string(q2order)!="" ) qorder = std::atoi(q2order);
  

  std::cout << "q2low " << q2lower << "\tq2up " << q2upper << std::endl;

  std::cout << "Process : " << nproc_.nproc << std::endl;

  if      ( nproc_.nproc == 1 )  
    {
      std::cout << " W+ production"; 
      pdf_function = "mcfmwp.config"; 
      glabel+="-Wplus";
      q2Low   = 6399.99, q2Up = 6400.01;
      nQ2bins = 3;
      qorder  = 1;
    }  
  else if ( nproc_.nproc == 6 )  
    {
      std::cout << " W- production"; 
      pdf_function = "mcfmwm.config"; 
      glabel+="-Wminus";
      q2Low   = 6399.99, q2Up = 6400.01;
      nQ2bins = 3;
      qorder  = 1;
    }  
  else if      ( nproc_.nproc == 11 )  
    {
      std::cout << " W+ + jet production"; 
      pdf_function = "mcfm-wpjet"; 
      glabel+="-WplusJet";
    }
  else if ( nproc_.nproc == 16 )  
    {
      std::cout << " W- + jet production"; 
      pdf_function = "mcfm-wmjet"; 
      glabel+="-WminusJet";
    }  
  else if ( nproc_.nproc == 31 ) 
    {
      std::cout << " Z production"; 
      pdf_function = "mcfm-z"; 
      glabel+="-Z0";
      q2Low = 8280.99, q2Up = 8281.01;
      nQ2bins = 3;
      qorder  = 1;
    }  
  else if ( (nproc_.nproc >= 41) && (nproc_.nproc <= 43)) 
    {
      std::cout << " Z-jet production"; 
      pdf_function = "mcfm-zjet"; 
      glabel+=TString::Format("-Zjet_%d",nproc_.nproc).Data();
      //      q2Low = 8280.99, q2Up = 8281.01;
    }  
  else if ( (nproc_.nproc >= 280) && (nproc_.nproc <= 286)) 
    {
      nXbins  = 30;
      xorder  = 6;

      //      q2Low = 8280.99, q2Up = 8281.01;
      //      q2Up    = 7000*7000;
      //      q2Up    = 1.01e8;
      //      nQ2bins = 15;
      //      qorder  = 4;

      Ngrids  = 3;

      nObsBins[0] = 12;
      nObsBins[1] = 13;
      nObsBins[2] = 10;
      
      double _eta[13] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 };
      double  _pt[14] = { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 };
      double _ptf[11] = { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 };

      obsBins[0] = _eta;
      obsBins[1] = _pt;
      obsBins[2] = _ptf;

      std::cout << " Photon production"; 
      pdf_function = "photonLO.config:photonNLO.config";
      glabel += TString::Format("-GammaProd_%d",nproc_.nproc).Data();
    }  
  else if ( nproc_.nproc == 13 )  
    {
      std::cout << " W+ + Cbar production"; 
      pdf_function = "mcfm-wpc"; 
      glabel+="-WplusCbar";
      lowest_order = 1;
    }  
  else if ( nproc_.nproc == 18 )  
    {
      std::cout << " W- + C production"; 
      pdf_function = "mcfm-wmc"; 
      glabel+="-WminusC";
      
      lowest_order = 1;
    }
  else if ( (nproc_.nproc == 141) || (nproc_.nproc == 142) || (nproc_.nproc==144) ||
	    (nproc_.nproc == 145) || (nproc_.nproc == 146) || (nproc_.nproc==147) ||
	    (nproc_.nproc == 148) || (nproc_.nproc == 149) || (nproc_.nproc==150) ||(nproc_.nproc == 151) ||
	    (nproc_.nproc == 157) || (nproc_.nproc == 158) || (nproc_.nproc==159) )
    {
      lowest_order = 2;
      
      if (nproc_.nproc == 141)
	{
	  std::cout << " TTbar production with 2 semi-leptonic decays"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-141";
	}
      else if (nproc_.nproc == 142)
	{
	  std::cout << " TTbar production with 2 semi-leptonic decays, corrections only in decays"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-142";
	}
      else if (nproc_.nproc == 144)
	{
	  std::cout << " TTbar production with 2 semi-leptonic decays, no correlations"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-144";
	}
      else if (nproc_.nproc == 145)
	{
	  std::cout << " TTbar production with 2 semi-leptonic decays, no spin correlations in top decays"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-145";
	}
      else if (nproc_.nproc == 146)
	{
	  std::cout << " TTbar production with Tbar hadronic decay, radiative corrections in production and decay"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-146";
	}
      else if (nproc_.nproc == 147)
	{
	  std::cout << " TTbar production with Tbar hadronic decay, radiative corrections in Tbar decay"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-147";
	}
      else if (nproc_.nproc == 148)
	{
	  std::cout << " TTbar production with Tbar hadronic decay, radiative corrections in W decay"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-148";
	}
      else if (nproc_.nproc == 149)
	{
	  std::cout << " TTbar production with T hadronic decay, radiative corrections in production and decay"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-149";
	}
      else if (nproc_.nproc == 150)
	{
	  std::cout << " TTbar production with T hadronic decay, radiative corrections in T decay"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-150";
	}
      else if (nproc_.nproc == 151)
	{
	  std::cout << " TTbar production with T hadronic decay, radiative correstions in W decay"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar-151";
	}
      else if (nproc_.nproc == 157)
	{
	  std::cout << " TTbar production"; 
	  pdf_function = "mcfm-TT"; 
	  glabel+="-TTbar";
	}
      else if(nproc_.nproc == 158)
	{
	  std::cout << " BBbar production"; 
	  pdf_function = "mcfm-BB"; 
	  glabel+="-BBbar";
	}
      else if (nproc_.nproc == 159)
	{
	  std::cout << " CCbar production"; 
	  pdf_function = "mcfm-CC"; 
	  glabel+="-CCbar";
	}
    }  
  else                           
    { 
      std::cerr << "don't know which process" << std::endl; 
      std::exit(-1); 
    } 
  std::cout << std::endl;

  /// Read the ckm matrix from mcfm to store in the grid automatically
  /// NB: we store 13 x 13 ckm matrix - mcfm only stores 11 x 11 so we 
  ///     must add 1 to each index to keep them aligned

  std::vector< std::vector<double> > ckm_vsq( 13, std::vector<double>( 13, 0 ) );
  
  for ( int ic=0 ; ic<__nf2__ ; ic++ ) { 
    for ( int ic1=0 ; ic1<__nf2__ ; ic1++ ) ckm_vsq[ic+1][ic1+1] = ckm_.vsq[ic][ic1];
  } 

  std::vector<std::vector<double> > __ckm( 3, std::vector<double>(3, 0) );
  __ckm[0][0] = cabib_.Vud;
  __ckm[0][1] = cabib_.Vus;
  __ckm[0][2] = cabib_.Vub;
  __ckm[1][0] = cabib_.Vcd;
  __ckm[1][1] = cabib_.Vcs;
  __ckm[1][2] = cabib_.Vcb;


  for(int igrid=0; igrid < Ngrids; igrid++) 
    {
    
      bool create_new = false;

      // if the file does not exist, create a new grid...
      if ( !file_exists(glabel+gridFiles[igrid]) )  create_new = true;

      // or if it does exists but root file is a zombie...
      if ( !create_new ) {  
	TFile testFile( (glabel+gridFiles[igrid]).c_str() );
	if ( testFile.IsZombie() ) create_new = true;
	testFile.Close();
      }

      if ( create_new ) 
	{ 
	  std::cout << "Creating NEW grid... " << std::endl;
	  
	  std::cout << "grid interpolation: " 
		    << "\tQ2 " << nQ2bins << " " <<  q2Low << " " <<  q2Up << " " <<  qorder   
		    << "\tx "  <<  nXbins << " " <<   xLow << " " <<   xUp << " " <<  xorder
		    << std::endl; 
	    


	  mygrid[igrid] = new appl::mcfm_grid( nObsBins[igrid], obsBins[igrid],      // obs bins
					       nQ2bins, q2Low, q2Up, qorder,         // Q2 bins and interpolation order
					       nXbins,   xLow,  xUp, xorder,         // x bins and interpolation order
					       pdf_function, lowest_order, nloops ); 
	  /// try reweighting for a bit
	  mygrid[igrid]->reweight(true);
	  mygrid[igrid]->setCMSScale( energy_.sqrts );


	  /// store the ckm matrix
	  //	  mygrid[igrid]->setckm2( ckm_vsq );

	  mygrid[igrid]->setckm( __ckm );

	  //	  grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  
	  std::cout << "reference histo name = " 
	       << mygrid[igrid]->getReference()->GetName() << std::endl;
	  
	  std::cout<<*mygrid[igrid]<<std::endl;  
	}
      else 
	{
	  std::cout << "Using existing grid file " << (glabel+gridFiles[igrid]) << std::endl;
	  
	  mygrid[igrid] = new appl::mcfm_grid(glabel+gridFiles[igrid]); //optimise grid x,Q2 bins
	  //       grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  mygrid[igrid]->getReference()->Reset();
	  mygrid[igrid]->optimise(nQ2bins, nXbins);
	  
	  std::cout<<*(mygrid[igrid])<<std::endl;  
	}
      // CTEQ like reweighting
      //      mygrid[igrid]->reweight( false );
    }

  runs = 0;
  isBooked = true;
  std::cout<<" ***********************************************"<<std::endl;
}


void fill_grid( const double evt[][mxpart] )
{

  static unsigned evtcounter = 1;
  if ( evtcounter%50000==0 ) std::cout << "fill_grid() filled " << evtcounter << " weights " << date(); 
  evtcounter++;

  if (!isBooked) 
    {    
      book_grid();
      return;
    }

  getObservable( evt );
  
  for(int igrid = 0; igrid < Ngrids; igrid++)
    if(cuts(igrid)){
 
      mygrid[igrid]->fillMCFM( Observable[igrid] );
    }

  runs++; // counter of number of events (shouldn't this be after cuts)? or is it the number of runs?"
}


//
// just normalise to bin width
//
void Normalise(TH1D* h) 
{ 
  for ( int ibin=1 ; ibin<=h->GetNbinsX() ; ibin++ ) 
    { 
      double width = h->GetBinLowEdge(ibin+1) - h->GetBinLowEdge(ibin);
      h->SetBinContent( ibin, h->GetBinContent(ibin)/width );
    }
  return;
}



void write_grid(double& xstotal)   // writes out grid after some events
{
  std::cout<<"Write out grids ..."<<std::endl;

  for(int igrid = 0; igrid < Ngrids; igrid++)
    {
      std::cout << "saving grid N=" << igrid+1 << "\tof " << Ngrids << "\t";

      std::system("sleep 1");
      

      mygrid[igrid]->setNormalised( false );
      mygrid[igrid]->run() = (iterat_.ncall2)*(iterat_.itmx2);
      
      mygrid[igrid]->untrim();
      int untrim_size = mygrid[igrid]->size();

      mygrid[igrid]->trim();
      int trim_size = mygrid[igrid]->size();

      /// scale up by number of weights
      (*mygrid[igrid]) *= mygrid[igrid]->run();
      
      // normalise the reference histogram by bin width
      Normalise( mygrid[igrid]->getReference() );

      /// now scale *down* the reference histogram because we've just 
      /// scaled it up ...
      //      mygrid[igrid]->getReference()->Scale( 1/mygrid[igrid]->run() );

      std::string filename = glabel+gridFiles[igrid];

#if 0

      std::string newpdfname = "";

      /// automatically optimise subprocesses - doesn't quite work yet
      /// when it does it may be moved into the grid itself
      if ( mygrid[igrid]->getGenpdf().find("basic")!=std::string::npos ) { 
	
       	std::stringstream ss;
	ss << "proc" << nproc_.nproc;
	newpdfname = ss.str();

	std::cout << "appl::grid::Write() " << newpdfname << std::endl;

	mygrid[igrid]->Write( filename, "grid", newpdfname );
      }
      else { 
	mygrid[igrid]->Write( glabel+gridFiles[igrid] );
      }

#else

      mygrid[igrid]->Write( filename );

#endif


      std::cout << "size(untrimmed)=" << untrim_size 
		<< "\tsize(trimmed)=" << trim_size 
		<< "\tfraction="      << 100.*trim_size/untrim_size << " %" << std::endl;

      //      int nsub = mygrid[igrid]->subProcesses();

      delete mygrid[igrid];
      
    }
  
  time_t _t;
  time(&_t);
  
  std::cout<<" ***********************************************"<<std::endl;
  std::cout<<" saved grids " << ctime(&_t);
  std::cout<<" ***********************************************"<<std::endl;
}
 


//
// ----------------------------------------------
//    analysis
// ----------------------------------------------
//

void getObservable(const double evt[][mxpart])
{
  // evt[momentum][particle number-1]
  // momentum[0,1,2,3] = (x,y,z,E)
  //

  // calculate observables
  for(int igrid = 0; igrid < Ngrids; igrid++)Observable[igrid] = 0.0; // initialize
  
  double p3[4] = {evt[3][2],evt[0][2],evt[1][2],evt[2][2]}; // (E,x,y,z)
  double p4[4] = {evt[3][3],evt[0][3],evt[1][3],evt[2][3]};
  double p5[4] = {evt[3][4],evt[0][4],evt[1][4],evt[2][4]};
  double p6[4] = {evt[3][5],evt[0][5],evt[1][5],evt[2][5]};
  double p7[4] = {evt[3][6],evt[0][6],evt[1][6],evt[2][6]};
  double p8[4] = {evt[3][7],evt[0][7],evt[1][7],evt[2][7]};

  double rapidity3 = 0.0;
  rapidity3 = (p3[0] + p3[3])/(p3[0] - p3[3]);
  (rapidity3 < 1e-13) ? rapidity3 = 100.0 : rapidity3 = 0.5*std::log(rapidity3);
  
  double rapidity4 = 0.0;
  rapidity4 = (p4[0] + p4[3])/(p4[0] - p4[3]);
  (rapidity4 < 1e-13) ? rapidity4 = 100.0 : rapidity4 = 0.5*std::log(rapidity4);

  double rapidity5 = 0.0;
  rapidity5 = (p5[0] + p5[3])/(p5[0] - p5[3]);
  (rapidity5 < 1e-13) ? rapidity3 = 100.0 : rapidity5 = 0.5*std::log(rapidity5);

  double rapidity6 = 0.0;
  rapidity6 = (p6[0] + p6[3])/(p6[0] - p6[3]);
  (rapidity6 < 1e-13) ? rapidity6 = 100.0 : rapidity6 = 0.5*std::log(rapidity6);

  double rapidity7 = 0.0;
  rapidity7 = (p7[0] + p7[3])/(p7[0] - p7[3]);
  (rapidity7 < 1e-13) ? rapidity7 = 100.0 : rapidity7 = 0.5*std::log(rapidity7);
  
  double rapidity8 = 0.0;
  rapidity8 = (p8[0] + p8[3])/(p8[0] - p8[3]);
  (rapidity8 < 1e-13) ? rapidity8 = 100.0 : rapidity8 = 0.5*std::log(rapidity8);
  
  double rapidity34 = 0.0;                      // rapidity of particle (3+4) in event record
  rapidity34  = (p3[0] + p4[0]) + (p3[3] + p4[3]);
  rapidity34 /= (p3[0] + p4[0]) - (p3[3] + p4[3]);  
  (rapidity34 < 1e-13) ? rapidity34 = 100.0 : rapidity34 = 0.5*std::log(rapidity34);

  double rapidity56 = 0.0;                      // rapidity of particle (5+6) in event record
  rapidity56  = (p5[0] + p6[0]) + (p5[3] + p6[3]);
  rapidity56 /= (p5[0] + p6[0]) - (p5[3] + p6[3]);
  (rapidity56 < 1e-13) ? rapidity56 = 100.0 : rapidity56 = 0.5*std::log(rapidity56);

  double rapidity345 = 0.0;                      // rapidity of particle (3+4+5) in event record
  rapidity345  = (p3[0] + p4[0] + p5[0]) + (p3[3] + p4[3] + p5[3]);
  rapidity345 /= (p3[0] + p4[0] + p5[0]) - (p3[3] + p4[3] + p5[3]);  
  (rapidity345 < 1e-13) ? rapidity345 = 100.0 : rapidity345 = 0.5*std::log(rapidity345);

  double rapidity678 = 0.0;                      // rapidity of particle (7+6+8) in event record
  rapidity678  = (p6[0] + p7[0] + p8[0]) + (p6[3] + p7[3] + p8[3]);
  rapidity678 /= (p6[0] + p7[0] + p8[0]) - (p6[3] + p7[3] + p8[3]);  
  (rapidity678 < 1e-13) ? rapidity678 = 100.0 : rapidity678 = 0.5*std::log(rapidity678);

  
  double pt3 = 0;
  pt3 = std::sqrt( p3[1]*p3[1] + p3[2]*p3[2] );
  
  double pt4 = 0;
  pt4 = std::sqrt( p4[1]*p4[1] + p4[2]*p4[2] );

  double pt5 = 0;
  pt5 = std::sqrt( p5[1]*p5[1] + p5[2]*p5[2] );

  double pt6 = 0;
  pt6 = std::sqrt( p6[1]*p6[1] + p6[2]*p6[2] );
  
  double pt7 = 0;
  pt7 = std::sqrt( p7[1]*p7[1] + p7[2]*p7[2] );

  double pt8 = 0;
  pt8 = std::sqrt( p8[1]*p8[1] + p8[2]*p8[2] );
  
  double pt34 = 0;
  pt34 = std::sqrt( std::pow(p3[1] + p4[1],2) + std::pow(p3[2] + p4[2],2) );

  double pt56 = 0;
  pt56 = std::sqrt( std::pow(p5[1] + p6[1],2) + std::pow(p5[2] + p6[2],2) );
  
  double pt345 = 0;
  pt345 = std::sqrt( std::pow(p3[1] + p4[1] + p5[1],2) + std::pow(p3[2] + p4[2] + p5[2],2) );

  double pt678 = 0;
  pt678 = std::sqrt( std::pow(p6[1] + p7[1] + p8[1],2) + std::pow(p6[2] + p7[2] + p8[2],2) );

  if ( (nproc_.nproc >= 280) && (nproc_.nproc <= 286)) {
    Observable[ 0 ] = std::fabs(rapidity3);
    Observable[ 1 ] = pt3;
    Observable[ 2 ] = pt3;
  }
  else {
    Observable[ 0 ] = rapidity3;
    Observable[ 1 ] = pt3;
    Observable[ 2 ] = rapidity4;
  }
  
  Observable[ 3 ] = pt4;
  Observable[ 4 ] = rapidity5;
  Observable[ 5 ] = pt5;

}

int cuts(int igrid)
{
  int fill = 0;
  switch(igrid)
    {
    case(0):
      if ( (nproc_.nproc >= 280) && (nproc_.nproc <= 286) ) { 
	if ( Observable[1]>100 ) fill = 1;
      }
      else fill = 1;
      break;
    case(1):
      if ( (nproc_.nproc >= 280) && (nproc_.nproc <= 286) ) { 
	if ( std::fabs(Observable[0])<1.37 ) fill = 1;
      }
      else fill = 1;
      break;
    case(2):
      if ( (nproc_.nproc >= 280) && (nproc_.nproc <= 286) ) { 
	if ( std::fabs(Observable[0])>=1.52 && std::fabs(Observable[0])<2.37 ) fill = 1;
      }
      else fill = 1;
      break;
    case(3):
      break;
    case(4):
      fill = 1;
      break;
    case(5):
      fill = 1;
      break;
    default: 
      std::cerr<<" In gridwrap.cpp::cuts(int). No such process : "<<igrid<<std::endl;
      std::exit(-1);
    }
  return fill;
}



// namespace mcfm_bridge;

/// function pointer hooks - set to 0 when no functions defined and applgrid not linked
extern void (*book_gridptr)();                         
extern void (*fill_gridptr)(const double evt[][mxpart] );
extern void (*write_gridptr)(double& );   


extern "C" bool setup_mcfmbridge() { 
  std::cout << "setup_mcfmbridge()" << std::endl;
  book_gridptr  = book_grid;
  fill_gridptr  = fill_grid;
  write_gridptr = write_grid;
  return true;
}

extern "C" bool setup_mcfmbridge_() { 
  std::cout << "setup_mcfmbridge()" << std::endl;
  return setup_mcfmbridge();
}
 

bool mcfm_bridge_status = setup_mcfmbridge();


