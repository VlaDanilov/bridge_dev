//
//   mcfm_interface.cxx        
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)    
//
//   $Id: mcfm_interfce.cxx, v   Fri  8 Nov 2013 09:07:01 CET sutt


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <map>
#include <algorithm>
#include <vector>

#include <cstdlib> 
#include <sys/time.h> 


#include "TFile.h"
#include "TH1D.h"
// #include "TMatrixT.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TCanvas.h"	

//#include "mcfm_procmap.h"
#include "mcfm_grid.h"
#include "grid_struct.h"
  
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
 //int ngrids=0;
 //int qfgrid(int ngrids);

// mcfm changed the the mxpart size between 6.7 and 6.8
// extern "C" int getmxpart_(void);i
// ***************************************************************
 GridStruct *myStruct = new GridStruct();
 unsigned short int  ngrids          = myStruct->gridNumber();
 appl::mcfm_grid* onegrid;
 std::vector<appl::mcfm_grid*> mygrid(ngrids);
 std::vector<TH1D> gridObs(ngrids);
 std::vector<std::string> gridFiles(ngrids);
 std::vector<double> observable(ngrids);
 std::vector<unsigned short int>    nobsBins(ngrids);
 std::vector<std::vector<double>> obsBinning;
 double *ObsBins = new double; 

 std::vector <unsigned short int> k_zoo; 
 std::vector <unsigned short int> g_zoo; 
 std::vector <double> Observable(ngrids);

 std::vector<double> Kinematics_pt(std::stoi(myStruct->nKinPar()));
 std::vector<double> Kinematics_eta(std::stoi(myStruct->nKinPar()));
 
 unsigned short int maxbins=1;
 unsigned short int maxpar =1;

 unsigned short int Nk_zoo = k_zoo.size();
 unsigned short int Ng_zoo = g_zoo.size(); 
 // double *obsbinning;
 // ************************************************************* 

 long unsigned int runs  =  0;
 bool isBooked           =  false;
 std::string glabel      =  "";

 void getObservable( const double evt[][mxpart] );

 int cuts(int);

 std::string date() { 
  time_t _t;
  time(&_t);
  return ctime(&_t);
}


void book_grid()  // inital grid booking
{ 
 /**************************************************************
 *******************  Grids Initialisation *********************/
  std::fill(mygrid.begin(),mygrid.end(),onegrid);
  for(int w=0; w<ngrids; w++){ 
     if(maxbins<std::stoi(myStruct->customGrid(w,2))){
        maxbins=std::stoi(myStruct->customGrid(w,2));
     }
  }
  obsBinning.resize(ngrids);
  std::cout<< "nKinPar = " << myStruct->nKinPar() << std::endl; 
  for(int g=0; g<ngrids; g++){
     
     gridFiles[g] = "_"+myStruct->customGrid(g,0)+"_"+myStruct->customGrid(g,1)+".root";
     std::cout << " gridFiles[g] = " << gridFiles[g] << std::endl;
     observable[g]= 0.0; 
     nobsBins[g]  = std::stoi(myStruct->customGrid(g,2));   
     std::cout <<" nobsBins[g]  = " << nobsBins[g]  << std::endl; 
  }
  /**************************************************************/
  k_zoo.clear(); g_zoo.clear();
  for(int nk = 0 ; nk < (std::stoi(myStruct->nKinPar())); nk++) k_zoo.push_back(std::stoi(myStruct->customKinematics(nk,0)));
  for(int ng = 0 ; ng < ngrids ; ng++)                          g_zoo.push_back(std::stoi(myStruct->customGrid(ng,1)));
 
  if (isBooked) return;
  
  //GridStruct *myStruct = new GridStruct();
  double xLow          = std::stod(myStruct->xLow());
  double xUp           = std::stod(myStruct->xUp());
  int  nloops          = std::stoi(myStruct->nLoops());
  // const int  ngrids    = myStruct->gridNumber();
  //std::string line     = myStruct->customGrid(1);  
 
  //  time_t _t;
  //  time(&_t);

  std::map<int, info> myMap;
  
  std::cout << " ***********************************************"<<std::endl;
  std::cout << " booking the grids " << date() << std::endl;
  std::cout << ngrids << std::endl;
  // set transform2 value
  double apramval=5.;
  appl::grid::transformvar(apramval);

  // lowest order in alphas	
  gridorder_.lowest_order=0;
  
  // how many loops
  // number of observables and binning for observables  

  std::cout << "Process : " << nproc_.nproc << std::endl;
  
  myMap[nproc_.nproc] = whichProcess(nproc_.nproc);
  gridorder_.lowest_order=myMap[nproc_.nproc].LowestOrder;     
  
  //<<--------------------------   Here You may change the default parameters  -------------------->>//
  if(myStruct->pdf_Fun()      !="default") myMap[nproc_.nproc].pdf_fun      =           myStruct->pdf_Fun()     ;
  if(myStruct->gLab()         !="default") myMap[nproc_.nproc].glab         =           myStruct->gLab()        ;
  if(myStruct->q2Up()         !="default") myMap[nproc_.nproc].q2up         = std::stod(myStruct->q2Up())       ;
  if(myStruct->q2Low()        !="default") myMap[nproc_.nproc].q2low        = std::stod(myStruct->q2Low())      ;
  if(myStruct->nQ2Bins()      !="default") myMap[nproc_.nproc].nq2bins      = std::stoi(myStruct->nQ2Bins())    ;
  if(myStruct->qOrder()       !="default") myMap[nproc_.nproc].qOrder       = std::stoi(myStruct->qOrder())     ;
  if(myStruct->nXBins()       !="default") myMap[nproc_.nproc].nxbins       = std::stoi(myStruct->nXBins())     ;
  if(myStruct->xOrder()       !="default") myMap[nproc_.nproc].xOrder       = std::stoi(myStruct->xOrder())     ;
  if(myStruct->LowestOrder()  !="default") myMap[nproc_.nproc].LowestOrder  = std::stoi(myStruct->LowestOrder());

  std::cout << " Channel            = " << myMap[nproc_.nproc].chan        << std::endl;
  std::cout << " Pdf fun            = " << myMap[nproc_.nproc].pdf_fun     << std::endl;
  std::cout << " glab               = " << myMap[nproc_.nproc].glab        << std::endl;
  std::cout << " Q^2 Low            = " << myMap[nproc_.nproc].q2low       << std::endl;
  std::cout << " Q^2 Up             = " << myMap[nproc_.nproc].q2up        << std::endl;
  std::cout << " Number of Q^2 bins = " << myMap[nproc_.nproc].nq2bins     << std::endl;
  std::cout << " q Order            = " << myMap[nproc_.nproc].qOrder      << std::endl;
  std::cout << " Number of x bins   = " << myMap[nproc_.nproc].nxbins      << std::endl;
  std::cout << " x Order            = " << myMap[nproc_.nproc].xOrder      << std::endl;
  std::cout << " Lowest order       = " << myMap[nproc_.nproc].LowestOrder << std::endl;

  glabel = "grid-"+myStruct->nXBins()+"-"+myStruct->xOrder()+"-"+myStruct->nQ2Bins()+"-"+myStruct->qOrder(); 
  
  const char* basename = std::getenv("appl_basename");
  if ( basename && std::string(basename)!="" ) glabel = basename;

  const char* q2upper = std::getenv("appl_q2up");
  if ( q2upper && std::string(q2upper)!="" ) myMap[nproc_.nproc].q2up = std::atof(q2upper);

  const char* q2lower = std::getenv("appl_q2low");
  if ( q2lower && std::string(q2lower)!="" ) myMap[nproc_.nproc].q2low = std::atof(q2lower);

  const char* q2order = std::getenv("appl_q2order");
  if ( q2order && std::string(q2order)!="" ) myMap[nproc_.nproc].qOrder = std::atoi(q2order);  

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


  for(int igrid=0; igrid < ngrids; igrid++) 
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

      double *obs = new double[std::stoi(myStruct->customGrid(igrid,2))+1];

        for(int b=0; b<(std::stoi(myStruct->customGrid(igrid,2))+1); b++){
           obs[b] =  std::stod(myStruct->customGrid(igrid,b+3));
           std::cout << "here it is = " << std::stod(myStruct->customGrid(igrid,b+3)) << std::endl;
      }

      if ( create_new ) 
	{ 
	  std::cout << "Creating NEW grid... " << std::endl;
	  
	  std::cout << "grid interpolation: " 
		    << "\tQ2 " << myMap[nproc_.nproc].nq2bins << " " <<  myMap[nproc_.nproc].q2low << " " <<  myMap[nproc_.nproc].q2up << " " <<  myMap[nproc_.nproc].qOrder   
		    << "\tx "  <<  myMap[nproc_.nproc].nxbins << " " <<   xLow << " " <<   xUp << " " <<  myMap[nproc_.nproc].xOrder
		    << std::endl; 

	  mygrid[igrid] = new appl::mcfm_grid( std::stoi(myStruct->customGrid(igrid,2)), obs,      // obs bins
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
	  
	  std::cout <<*mygrid[igrid]<<std::endl;  
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

        gridObs[igrid] = TH1D(myStruct->customGrid(igrid,0).c_str(), myStruct->customGrid(igrid,0).c_str(), std::stoi(myStruct->customGrid(igrid,2)), obs);
        gridObs[igrid].SetDirectory(0);       	
        
        delete [] obs;
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
  for(int igrid = 0; igrid < ngrids; igrid++)
    if(cuts(igrid)){  
        mygrid[igrid]->fillMCFM( std::fabs(Observable[0]) );
        gridObs[igrid].Fill(     std::fabs(Observable[0]),gridevent_.refwt); 
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

  for(int igrid = 0; igrid < ngrids; igrid++)
    {
      std::cout << "saving grid N=" << igrid+1 << "\tof " << ngrids << "\t";

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
       
       //****************Plots*****************
        Normalise(&gridObs[igrid]);	
        gridObs[igrid].Write("",TObject::kOverwrite);	
       //**************************************
       Kinematics->Close();

#endif


      std::cout << "size(untrimmed)=" << untrim_size 
		<< "\tsize(trimmed)=" << trim_size 
		<< "\tfraction="      << 100.*trim_size/untrim_size << " %" << std::endl;
     
      //      int nsub = mygrid[igrid]->subProcesses();

    }
  
  delete[] ObsBins;
  
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

void getObservable( const double evt[][mxpart])
{
 // evt[momentum][particle number-1]
 // momentum[0,1,2,3] = (x,y,z,E)

 // calculate observables
 //************My obeservables************
 
 for(int igrid = 0; igrid < ngrids; igrid++)Observable[igrid] = 0.0; // initialize
 
 for(unsigned short int ik=0; ik<std::stoi(myStruct->nKinPar()); ik++){
   TLorentzVector kp(evt[0][k_zoo[ik]-1],evt[1][k_zoo[ik]-1],evt[2][k_zoo[ik]-1],evt[3][k_zoo[ik]-1]);
   Kinematics_pt[ik]  = kp.Pt();   
   Kinematics_eta[ik] = kp.Eta();
 }

 for(unsigned short int ig=0; ig<ngrids; ig++){
   Observable[ig] = 0.0;
   TLorentzVector gp(evt[0][g_zoo[ig]-1],evt[1][g_zoo[ig]-1],evt[2][g_zoo[ig]-1],evt[3][g_zoo[ig]-1]); //particle to fill grid  (number-1)

        if( (myStruct->customGrid(ig,0) == "Eta")    || ((myStruct->customGrid(ig,0)) == "PseudoRapidity") ) Observable[ig] = gp.Eta(); 
   else if( (myStruct->customGrid(ig,0) == "Pt")     || ((myStruct->customGrid(ig,0)) == "pT") )             Observable[ig] = gp.Pt(); 
   else if( (myStruct->customGrid(ig,0) == "Mass")   || ((myStruct->customGrid(ig,0)) == "M") )              Observable[ig] = gp.M(); 
   else if( (myStruct->customGrid(ig,0) == "Energy") || ((myStruct->customGrid(ig,0)) == "E") )              Observable[ig] = gp.E(); 
   else if( (myStruct->customGrid(ig,0) == "Et")                                              )              Observable[ig] = gp.Et(); 
   else if( (myStruct->customGrid(ig,0) == "Rapidity") )                                                     Observable[ig] = gp.Rapidity();
   else {
         std::cout << "=====================> ERROR <=================== " << std::endl;
         std::cout << " \"The variable name '" << (myStruct->customGrid(ig,0))  << "' is UNKNOWN\" " << std::endl;
         std::cout << " \"Please check the steering file\" " << std::endl << std::endl; 
         std::cout << "Terminating..." << std::endl;
         std::cout << "================================================= " << std::endl; 
         std::terminate();
        } 
 }
}


int cuts(int igrid)
{
  int fill = 1;
  for(unsigned short int kp=0; kp<std::stoi(myStruct->nKinPar()); kp++){
      if(          Kinematics_pt[kp]   <  std::stod(myStruct->customKinematics(kp,1))){return 0;}   
      if(std::fabs(Kinematics_eta[kp]) >= std::stod(myStruct->customKinematics(kp,2))){return 0;}   
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


