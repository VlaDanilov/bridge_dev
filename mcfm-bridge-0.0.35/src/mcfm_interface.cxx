//
//   mcfm_interface.cxx        
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)    
//
//   $Id: mcfm_interfce.cxx, v   Fri  8 Nov 2013 09:07:01 CET sutt


//
//   modified by Vladyslav Danilov(vladyslav.danilov@cern.ch)
//   Monday 19 August 2019
//

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
#include "TVectorT.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "mcfm_grid.h"
#include "grid_input.h"   // Steering file reader class
#include "mcfm_procmap.h" // Map with process' default parameters   


bool file_exists( const std::string& filename ) { 
  struct stat sb;
  if ( stat( filename.c_str(), &sb)==0 ) return true; // && S_ISREG(sb.st_mode ))
  else return false;
}

// Check if the process number exists in process map
// Maybe it's better to make an option for user to define his own parameters for a certain(undefined in the map) process but with basic pdf decomposition?
bool proc_exists(int &num) {
  if(procmap.count(num) > 0) return false;
  else return true;
}

 static const int mxpart = 14;    // mcfm parameter : max number of partons in event record. defined in Inc/constants.f

 // Initialize object of a reader class
 grid_input *steeringFile   = new grid_input();
 unsigned short int  ngrids = std::stoi(steeringFile->gridNumber());

 // Grids
 appl::mcfm_grid* onegrid;
 std::vector<appl::mcfm_grid*> mygrid(ngrids);
 std::vector<std::string> gridFiles(ngrids);

 // Particles to make a kinematic cuts on
 std::vector <unsigned short int> k_zoo;
 // Numbers of particles to fill a grid(according to process description MCFM)
 std::vector <unsigned short int> g_zoo;
 
 // Observable
 std::vector <double> Observable(ngrids);
 
 // Containers for values of kinematic cuts
 std::vector<double> Kinematics_pt(std::stoi(steeringFile->nKinPar()));
 std::vector<double> Kinematics_eta(std::stoi(steeringFile->nKinPar()));

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

  // Filling containers 
  k_zoo.clear(); g_zoo.clear();
  for(int nk = 0 ; nk < std::stoi(steeringFile->nKinPar()); nk++) k_zoo.push_back(std::stoi(steeringFile->customKinematics(nk,0)));
  for(int ng = 0 ; ng < ngrids ; ng++)                            g_zoo.push_back(std::stoi(steeringFile->customGrid(ng,1)));

  std::cout << "-------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << " Date: " << date() << std::endl;
  std::cout << " Booking " << ngrids << " grids" << std::endl;     
  
  // set transform2 value
  double apramval=5.;
  appl::grid::transformvar(apramval);
   
  std::cout << "  MCFM process number : " << nproc_.nproc << std::endl;
  std::cout << "  Number of particles with kinematic cuts = " << steeringFile->nKinPar() << std::endl;
 
  // filling structure 'info' with a certain process parameters from the process' map 
  info myprocess = procmap[nproc_.nproc];
  
  // Checking if the process exists in the map
  if(proc_exists(nproc_.nproc)){
    std::cerr << " Error! Process " << nproc_.nproc << " doesn't exist ! " << std::endl;
    exit(0);
  }
 
  //Lowest order in alpha_s
  gridorder_.lowest_order=myprocess.LowestOrder; 

  // If grid_setup.txt has values written as 'default' - those values will be taken from process map.
  if(steeringFile->pdfFun()      != "default") myprocess.pdf_fun     =           steeringFile->pdfFun()      ;
  if(steeringFile->gLab()        != "default") myprocess.glab        =           steeringFile->gLab()        ;
  if(steeringFile->q2Up()        != "default") myprocess.q2up        = std::stod(steeringFile->q2Up())       ;
  if(steeringFile->q2Low()       != "default") myprocess.q2low       = std::stod(steeringFile->q2Low())      ;
  if(steeringFile->nQ2Bins()     != "default") myprocess.nq2bins     = std::stoi(steeringFile->nQ2Bins())    ;
  if(steeringFile->qOrder()      != "default") myprocess.qOrder      = std::stoi(steeringFile->qOrder())     ;
  if(steeringFile->nXBins()      != "default") myprocess.nxbins      = std::stoi(steeringFile->nXBins())     ;
  if(steeringFile->xOrder()      != "default") myprocess.xOrder      = std::stoi(steeringFile->xOrder())     ;
  if(steeringFile->lowestOrder() != "default") myprocess.LowestOrder = std::stoi(steeringFile->lowestOrder());
   
  std::cout << " Channel            = " << myprocess.chan        << std::endl;
  std::cout << " Pdf fun            = " << myprocess.pdf_fun     << std::endl;
  std::cout << " glab               = " << myprocess.glab        << std::endl;
  std::cout << " Q^2 Low            = " << myprocess.q2low       << std::endl;
  std::cout << " Q^2 Up             = " << myprocess.q2up        << std::endl;
  std::cout << " Number of Q^2 bins = " << myprocess.nq2bins     << std::endl;
  std::cout << " Interpolation order of Q = " << myprocess.qOrder      << std::endl;
  std::cout << " Number of x bins   = " << myprocess.nxbins      << std::endl;
  std::cout << " Interpolation order of x = " << myprocess.xOrder      << std::endl;
  std::cout << " Lowest order       = " << myprocess.LowestOrder << std::endl;
  
  std::cout << " ------------------------------------------------------------------------------------------------" << std::endl;
 
  glabel = "grid-"+std::to_string(myprocess.nxbins)+"-"+std::to_string(myprocess.xOrder)+"-"+std::to_string(myprocess.nq2bins)+"-"+std::to_string(myprocess.qOrder); 
 
  const char* basename = std::getenv("appl_basename");
  if ( basename && std::string(basename)!="" ) glabel = basename;
   
  const char* q2upper = std::getenv("appl_q2up");
  if ( q2upper && std::string(q2upper)!="" ) myprocess.q2up = std::atof(q2upper);
   
  const char* q2lower = std::getenv("appl_q2low");
  if ( q2lower && std::string(q2lower)!="" ) myprocess.q2low = std::atof(q2lower);
   
  const char* q2order = std::getenv("appl_q2order");
  if ( q2order && std::string(q2order)!="" ) myprocess.qOrder = std::atoi(q2order); 

  // std::cout << "q2low " << q2lower << "\tq2up " << q2upper << std::endl;
  

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

  std::fill(mygrid.begin(),mygrid.end(),onegrid);

  for(int igrid=0; igrid < ngrids; igrid++) 
    {
      std::cout << " ************************************************************************************************ "<<std::endl;
      gridFiles[igrid] = "_p"+steeringFile->customGrid(igrid,1)+"_"+steeringFile->customGrid(igrid,0)+".root";
      std::cout <<" Grid file suffix : " << gridFiles[igrid] << std::endl;
      std::cout <<" Number of observable bins = " << std::stoi(steeringFile->customGrid(igrid,3))  << std::endl; 

      bool create_new = false;

      // if the file does not exist, create a new grid...
      if ( !file_exists(glabel+gridFiles[igrid]) )  create_new = true;

      // or if it does exists but root file is a zombie...
      if ( !create_new ) {  
	TFile testFile( (glabel+gridFiles[igrid]).c_str() );
	if ( testFile.IsZombie() ) create_new = true;
	testFile.Close();
      }

      // Define a custom size for observable
      double *obs = new double[std::stoi(steeringFile->customGrid(igrid,3))+1];

      // Read, fill and print observable binning
      std::cout << " < Binning for observable " << steeringFile->customGrid(igrid,0) << " of particle " << steeringFile->customGrid(igrid,1) << " in process " <<  myprocess.chan << " > :"  << std::endl;
      for(int b=0; b<(std::stoi(steeringFile->customGrid(igrid,3))+1); b++){
           obs[b] =  std::stod(steeringFile->customGrid(igrid,b+4));
           if(b==0) continue;
           std::cout << " Bin " << b << " = [" <<  std::stod(steeringFile->customGrid(igrid,b+3)) << " - "<< std::stod(steeringFile->customGrid(igrid,b+4))<< "];" << std::endl;
      }

      if ( create_new ) 
	{ 
	  std::cout << " Creating NEW grid... " << std::endl;
	  
          std::cout << " grid interpolation: " 
		    << "\n Q^2: " << myprocess.nq2bins << " bins in range[" <<  myprocess.q2low   << "-" <<  myprocess.q2up   << "], interpolation order " <<  myprocess.qOrder   
		    << "\n   x: " << myprocess.nxbins  << " bins in range[" <<  std::stod(steeringFile->xLow()) << "-" <<  std::stod(steeringFile->xUp()) << "], interpolation order " <<  myprocess.xOrder
		    << std::endl;	  

	  mygrid[igrid] = new appl::mcfm_grid( std::stoi(steeringFile->customGrid(igrid,3)), obs,      // obs bins
					       myprocess.nq2bins, myprocess.q2low,   myprocess.q2up,   myprocess.qOrder,         // Q2 bins and interpolation order
					       myprocess.nxbins,  std::stod(steeringFile->xLow()), std::stod(steeringFile->xUp()), myprocess.xOrder,         // x bins and interpolation order
					       myprocess.pdf_fun, gridorder_.lowest_order,     std::stoi(steeringFile->nLoops()) ); 
	  
          /// try reweighting for a bit
          mygrid[igrid]->reweight(true);
	  mygrid[igrid]->setCMSScale( energy_.sqrts );

	  /// store the ckm matrix
	  // mygrid[igrid]->setckm2( ckm_vsq );

	  mygrid[igrid]->setckm( __ckm );

	  // grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  
	  std::cout << " reference histo name : " << mygrid[igrid]->getReference()->GetName() << std::endl;

	  //std::cout<<*mygrid[igrid]<<std::endl;  
	}
      else 
	{
	  std::cout << " Using existing grid file : " << (glabel+gridFiles[igrid]) << std::endl;
	  
	  mygrid[igrid] = new appl::mcfm_grid(glabel+gridFiles[igrid]); //optimise grid x,Q2 bins
	  // grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  mygrid[igrid]->getReference()->Reset();
	  mygrid[igrid]->optimise(myprocess.nq2bins, myprocess.nxbins);
	  
	  // std::cout<<*(mygrid[igrid])<<std::endl;  
	}
      // CTEQ like reweighting
      // mygrid[igrid]->reweight( false );
        
      delete [] obs;
    }

  runs = 0;
  isBooked = true;
  std::cout << " ************************************************************************************************ "<<std::endl;
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
  
  for(int igrid = 0; igrid < ngrids; igrid++)
    if(cuts(igrid)) mygrid[igrid]->fillMCFM( Observable[igrid] );
    
  runs++;
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
  
  for(int igrid = 0; igrid < ngrids; igrid++)
    {
      std::cout << "saving grid N = " << igrid+1 << "\tof " << ngrids << "\t" << std::endl;

      std::system("sleep 1");

      mygrid[igrid]->setNormalised( false );
 
      // Define number of iterations times number of events and 
      // make scaling condition to be true.
      // In Wirte() method this condition can be used to scale before writing the grid    
      mygrid[igrid]->run() = (iterat_.ncall2)*(iterat_.itmx2);
      
      mygrid[igrid]->untrim();
      int untrim_size = mygrid[igrid]->size();

      mygrid[igrid]->trim();
      int trim_size = mygrid[igrid]->size();

      /// scale up by number of weights
      (*mygrid[igrid]) *= mygrid[igrid]->run();
      
      // normalise the reference histogram by bin width
      Normalise( mygrid[igrid]->getReference() );

      // now scale *down* the reference histogram because we've just 
      // scaled it up ...
      /// We had to scale it before this procedure was moved into the Write() method, now it's a rudiment piece of code. 
      //  mygrid[igrid]->getReference()->Scale( 1/mygrid[igrid]->run() );

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
      
      delete mygrid[igrid];
    }
  
  delete steeringFile;
 
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
 // momentum[0,1,2,3] = (px,py,pz,E)

  // Select particles for kinematic cuts
  for(unsigned short int ik=0; ik<std::stoi(steeringFile->nKinPar()); ik++){
    TLorentzVector kp(evt[0][k_zoo[ik]-1],evt[1][k_zoo[ik]-1],evt[2][k_zoo[ik]-1],evt[3][k_zoo[ik]-1]);
    Kinematics_pt[ik]  = kp.Pt();   
    Kinematics_eta[ik] = kp.Eta();
  }

  // Select particles as observables
  for(unsigned short int ig=0; ig<ngrids; ig++){
    Observable[ig] = 0.0;
    TLorentzVector gp(evt[0][g_zoo[ig]-1],evt[1][g_zoo[ig]-1],evt[2][g_zoo[ig]-1],evt[3][g_zoo[ig]-1]); //particle to fill grid  (number-1)

    // Matching and defining* requested observable 
    // *Here for observables calculation ROOT methods are used 
         if( (steeringFile->customGrid(ig,0) == "Eta")      || ((steeringFile->customGrid(ig,0)) == "eta") ) Observable[ig] = gp.Eta(); 
    else if( (steeringFile->customGrid(ig,0) == "Pt")       || ((steeringFile->customGrid(ig,0)) == "pt")  ) Observable[ig] = gp.Pt(); 
    else if( (steeringFile->customGrid(ig,0) == "Et")                                                      ) Observable[ig] = gp.Et(); 
    else if( (steeringFile->customGrid(ig,0) == "M")        || ((steeringFile->customGrid(ig,0)) == "m")   ) Observable[ig] = gp.M(); 
    else if( (steeringFile->customGrid(ig,0) == "Energy")   || ((steeringFile->customGrid(ig,0)) == "E")   ) Observable[ig] = gp.E(); 
    else if( (steeringFile->customGrid(ig,0) == "Rapidity") || ((steeringFile->customGrid(ig,0)) == "y")   ) Observable[ig] = gp.Rapidity();
    else {
          std::cerr << "\n"; 
          std::cerr << " =====================> ERROR <=================== \n";
          std::cerr << "  \"The Observable name '" << (steeringFile->customGrid(ig,0))  << "' is UNKNOWN\"\n ";
          std::cerr << "  Available names:\n";
          std::cerr << "   rapidity: 'Rapidity', 'y'\n" ;
          std::cerr << "   energy:   'Energy', 'E'\n" ;
          std::cerr << "   pseudorapidity: 'Eta', 'eta'\n" ;
          std::cerr << "   transverse momentum: 'Pt',  'pt'\n" ;
          std::cerr << "   mass : 'M',   'm'\n" ;
          std::cerr << "   transverse energy: 'Et'\n" ;
          std::cerr << " \"Please check the steering file\"\n "  ;
          std::cerr << " Terminating...\n" ;
          std::cerr << " ================================================= \n";
          std::cerr << " " << std::endl; 
          exit(0);
         }
   
    if(std::stoi(steeringFile->customGrid(ig,2)) == 1) Observable[ig] = std::fabs(Observable[ig]);

  }
}

int cuts(int igrid)
{
  int fill = 1;
  // Loop over kinematic particles to perform pt and eta custom cuts defined in a steering file(grid_setup.txt)
  for(unsigned short int kp=0; kp<std::stoi(steeringFile->nKinPar()); kp++){
    if(          Kinematics_pt[kp]   <  std::stod(steeringFile->customKinematics(kp,1))){
      fill = 0;
    }
    // Here, probably, is better to remove abs value of eta as default option and add a boolian condition, 
    // then kinematic parameters in steering file would look like [numberOfParticle ptCut -etaCut +etaCut doAbsEta(1 or 0)] 
    if( Kinematics_eta[kp] <= std::stod(steeringFile->customKinematics(kp,2)) || Kinematics_eta[kp] >= std::stod(steeringFile->customKinematics(kp,3)) ){
      fill = 0;
    }   
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


