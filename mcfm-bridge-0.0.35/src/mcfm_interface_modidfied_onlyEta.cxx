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
#include <map>

#include <cstdlib> 
#include <sys/time.h> 


#include "TFile.h"
#include "TH1D.h"
// #include "TMatrixT.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TCanvas.h"	

//#include "mcfm_procmap.h"
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

 double E;
 double Mass_ee;
 double delta_R[5];
 double Phi[7];
 double pseudorapidity[7];
 double pT[8];
 static const int mxpart = 14;    // mcfm parameter : max number of partons in event record. defined in Inc/constants.f

/// was 12 for mcfm 6.7
/// static const int mxpart = 12;    // mcfm parameter : max number of partons in event record. defined in Inc/constants.f

static const int _Ngrids =1; 
const int Ngrids = 1;

/// mcfm changed the the mxpart size between 6.7 and 6.8
// extern "C" int getmxpart_(void);
 TH1D LeptEta[Ngrids];

appl::mcfm_grid* mygrid[_Ngrids];

static const char* gridFiles[_Ngrids] = {
    "_LeptEta.root"};


static double Observable[_Ngrids] = { 0};//,  0,  0,  0,  0,  0,  0, 0};   // observable array    //  <= Change here for var.
// int             nObsBins[_Ngrids] = { 40}, 45, 40, 45, 40, 45 }; // eta4, pt4 cental eta-bin, pt4 forward eta-bin
int  nObsBins[_Ngrids] = { 11};//, 9,  9, 12, 10, 8,10, 9} ;//  <= Change here for var.



static const double Eta_l[12] =  { 
      0.0, 0.2, 0.4,
      0.6,  0.8,  1.0,
     1.2, 1.4,  1.6,
      1.85, 2.1, 2.4
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
  std::map<int, info> myMap;
  
  std::cout<<" ***********************************************"<<std::endl;
  std::cout<<" booking the grids " << date() << std::endl;
  
  // binning information for the grid constructor
  double xLow    = 1.0e-9, xUp = 1.0;
  
  // set transform2 value
  double apramval=5.;
  appl::grid::transformvar(apramval);

  // lowest order in alphas	
  gridorder_.lowest_order=0;
  
  // how many loops
  int nloops = 1;
  // number of observables and binning for observables  
  const double *obsBins[_Ngrids] = {Eta_l};

  glabel = "grid-40-6-15-3"; 

  std::cout << "Process : " << nproc_.nproc << std::endl;
  
  myMap[nproc_.nproc] = whichProcess(nproc_.nproc);
  gridorder_.lowest_order=myMap[nproc_.nproc].LowestOrder;     
  
 //<<--------------------------   Here You may change the default parameters  -------------------->>//
 // myMap[nproc_.nproc].chan="channel";
 // myMap[nproc_.nproc].pdf_fun="basic";
 // myMap[nproc_.nproc].glab="-smth";
 // myMap[nproc_.nproc].q2up=6500.0*6500.0*5;
 // myMap[nproc_.nproc].q2low=1.0;
 // myMap[nproc_.nproc].nq2bins=20;
 // myMap[nproc_.nproc].qOrder=5;
 // myMap[nproc_.nproc].nxbins=40;
 // myMap[nproc_.nproc].xOrder=0;
 
 
  std::cout << "myMap[nproc_.nproc].chan              = " << myMap[nproc_.nproc].chan << std::endl;
  std::cout << "myMap[nproc_.nproc].pdf_fun          = " << myMap[nproc_.nproc].pdf_fun << std::endl;
  std::cout << "myMap[nproc_.nproc].glab               = " << myMap[nproc_.nproc].glab  << std::endl;
 
  std::cout << "myMap[nproc_.nproc].q2low            = " << myMap[nproc_.nproc].q2low << std::endl;
  std::cout << "myMap[nproc_.nproc].q2up              = " << myMap[nproc_.nproc].q2up << std::endl;
  std::cout << "myMap[nproc_.nproc].nq2bins         = " << myMap[nproc_.nproc].nq2bins << std::endl;
  std::cout << "myMap[nproc_.nproc].qOrder           = " << myMap[nproc_.nproc].qOrder << std::endl;
  std::cout << "myMap[nproc_.nproc].nxbins           = " << myMap[nproc_.nproc].nxbins << std::endl;
  std::cout << "myMap[nproc_.nproc].xOrder           = " << myMap[nproc_.nproc].xOrder << std::endl;
  std::cout << "myMap[nproc_.nproc].LowestOrder = " << myMap[nproc_.nproc].LowestOrder << std::endl;

  const char* basename = std::getenv("appl_basename");
  if ( basename && std::string(basename)!="" ) glabel = basename;

  const char* q2upper = std::getenv("appl_q2up");
  if ( q2upper && std::string(q2upper)!="" ) myMap[nproc_.nproc].q2up = std::atof(q2upper);

  const char* q2lower = std::getenv("appl_q2low");
  if ( q2lower && std::string(q2lower)!="" ) myMap[nproc_.nproc].q2low = std::atof(q2lower);

  const char* q2order = std::getenv("appl_q2order");
  if ( q2order && std::string(q2order)!="" ) myMap[nproc_.nproc].qOrder = std::atoi(q2order);  

 // std::cout << "q2low " << q2lower << "\tq2up " << q2upper << std::endl;
  std::cout<<std::endl;
	 
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
		    << "\tQ2 " << myMap[nproc_.nproc].nq2bins << " " <<  myMap[nproc_.nproc].q2low << " " <<  myMap[nproc_.nproc].q2up << " " <<  myMap[nproc_.nproc].qOrder   
		    << "\tx "  <<  myMap[nproc_.nproc].nxbins << " " <<   xLow << " " <<   xUp << " " <<  myMap[nproc_.nproc].xOrder
		    << std::endl; 
	    


	  mygrid[igrid] = new appl::mcfm_grid( nObsBins[igrid], obsBins[igrid],      // obs bins
					       myMap[nproc_.nproc].nq2bins, myMap[nproc_.nproc].q2low, myMap[nproc_.nproc].q2up, myMap[nproc_.nproc].qOrder,         // Q2 bins and interpolation order
					       myMap[nproc_.nproc].nxbins,   xLow,  xUp, myMap[nproc_.nproc].xOrder,         // x bins and interpolation order
					       myMap[nproc_.nproc].pdf_fun, gridorder_.lowest_order, nloops ); 
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
	  mygrid[igrid]->optimise(myMap[nproc_.nproc].nq2bins, myMap[nproc_.nproc].nxbins);
	  
	  std::cout<<*(mygrid[igrid])<<std::endl;  
	}
       //*******Richards-type histograms
       
	LeptEta[igrid] = TH1D("LeptEta","LeptEta",11, Eta_l);	
	LeptEta[igrid].SetDirectory(0);
	
	//***************************
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
        mygrid[igrid]->fillMCFM(  std::fabs(Observable[0])  );
	LeptEta[igrid].Fill(  std::fabs(  Observable[0]  ),gridevent_.refwt); //,gridevent_.refwt
    }
  runs++; // !!!
}


//
// just normalise to bin width
//
void Normalise(TH1D* h) 
{ 
  for ( int ibin=1 ; ibin<=h->GetNbinsX() ; ibin++ ) 
    { 
      double width = h->GetBinLowEdge(ibin+1) - h->GetBinLowEdge(ibin);    /// UNCOMMENT	
      h->SetBinContent( ibin, h->GetBinContent(ibin)/(width) );
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
      TFile *Kinematics = new TFile(filename.c_str(), "UPDATE");
       
       //********Richard-type plots*********************
	
	Normalise(&LeptEta[igrid]);
	
	LeptEta[igrid].Write("",TObject::kOverwrite);	
	
       //**************************************
	//****
       Kinematics->Close();

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

  // calculate observables
  for(int igrid = 0; igrid < Ngrids; igrid++)Observable[igrid] = 0.0; // initialize
  
 
 TLorentzVector p3(evt[0][2],evt[1][2],evt[2][2],evt[3][2]); // neutrino
 TLorentzVector p4(evt[0][3],evt[1][3],evt[2][3],evt[3][3]); // positron
 
//**********
 
 double pseudorapidity3 = 0.0; 
 double pseudorapidity4 = 0.0; 
 
 pseudorapidity3 = p3.PseudoRapidity();
 pseudorapidity4 = p4.PseudoRapidity(); 
 
 //*********
 double pT3 = 0.0;
 double pT4 = 0.0;
 
 pT3 = p3.Pt();
 pT4 = p4.Pt();
  
 //****Assignment to general Array   
  
 pseudorapidity[0] = pseudorapidity3;
 pseudorapidity[1] = pseudorapidity4;
 
 Observable[0] = pseudorapidity[1];
 
 pT[0] = pT3;
 pT[1] = pT4;
  
}

//###############################

int cuts(int igrid)
{
  int fill = 1;
    if(std::fabs(pseudorapidity[1]) >= 2.4 ){return 0;}
    if( pT[1] < 25){return 0;}  
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


