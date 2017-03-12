#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>

#include <stdio.h>

#include "amconfig.h"


#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"


#include "appl_grid/appl_timer.h"



#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TPad.h"
#include "TStyle.h"

#include "hoppet_v1.h"

// #include "DrawLabel.h"
// #include "Atlas/AtlasStyle.h"

// lhapdf routines
#include "LHAPDF/LHAPDF.h"
extern "C" void evolvepdf_(const double& x, const double& Q, double* xf); 
extern "C" double alphaspdf_(const double& Q);

// double fnalphas_(const double& Q) { return alphaspdf_(Q); }
// void   fnpdf_(const double& x, const double& Q, double *f) { return evolvepdf_(x, Q, f); }


#define DBG false


static const int nFlavours = 5;


std::string fname = "";

void smooth(TH1D* h) { 
  for ( int i=2 ; i<h->GetNbinsX() ; i++ ) {
    double y0 = h->GetBinContent(i-1);
    double y  = h->GetBinContent(i);
    double y1 = h->GetBinContent(i+1);

    if ( std::fabs(y-y0)>0.2 ) { 
      y = 0.5*(y0+y1);
      h->SetBinContent( i,  y );
    }
  }
}


void binwidth(TH1D* h) { 
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) {
    double dx = h->GetBinLowEdge(i+1) - h->GetBinLowEdge(i);
    h->SetBinContent( i,  h->GetBinContent(i)/dx );
    h->SetBinError( i,  h->GetBinError(i)/dx );
  }
}


double getRealMinimum( TH1D* h ) { 
  double _min = 0;
  bool _set = false;
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) {
    double y = h->GetBinContent(i);
    if ( y!=0 ) { 
      if ( _set ) { 
	if ( y<_min ) _min = y;
      }
      else { 
	_min = y;
	_set = true;
      }
    }
  }

  return _min;
}


double getRealMaximum( TH1D* h ) { 
  double _max = 0;
  bool _set = false;
  for ( int i=1 ; i<=h->GetNbinsX() ; i++ ) {
    double y = h->GetBinContent(i);
    if ( y!=0 ) { 
      if ( _set ) { 
	if ( y>_max ) _max = y;
      }
      else { 
	_max = y;
	_set = true; 
      }
    }
  }

  return _max;
}


double scale_factor = 1;

TH1D* divide( const TH1D* h1, const TH1D* h2 ) {

  if ( h1==NULL || h2==NULL ) return NULL;
 
  TH1D* h = (TH1D*)h1->Clone();

  if ( DBG ) std::cout << "histograms " << h1->GetTitle() << " " << h2->GetTitle() << std::endl;

  
  for ( int i=1 ; i<=h1->GetNbinsX() ; i++ ) { 
    double b  = h2->GetBinContent(i);
    //    double be = h2->GetBinError(i);
    double t  = h1->GetBinContent(i);
    //  double te = h1->GetBinError(i);

    double r  = ( b!=0 ? t/b : 0 );
    //    double re = ( b!=0 ? sqrt((r+1)*r/b) : 0 );
    double re = 0;
    
    h->SetBinContent( i, r );
    h->SetBinError( i, re ) ;

    //    if ( DBG ) std::cout << "\tx=" << h->GetBinCenter(i) << "\tratio=" << r << std::endl;
  } 

  double hmin = h->GetBinContent(1);
  double hmax = h->GetBinContent(1);
  
  for ( int i=2 ; i<=h->GetNbinsX() ; i++ ) { 
    double d = h->GetBinContent(i);
    if ( hmin>d ) hmin=d;
    if ( hmax<d ) hmax=d; 
  }

  if ( DBG ) std::cout << "\tmin ratio = " << hmin << "\tmax ratio = " << hmax << std::endl;
  
  double _max = getRealMaximum( h );
  double _min = getRealMinimum( h );

  if ( DBG ) std::cout << "\tmin ratio = " << _min << "\tmax ratio = " << _max << std::endl;
  
  double mindelta = 0;

  for ( double delta=100 ; delta>1e-7 ; delta*=0.1 ) {

    if ( _min>(1-delta) && _max<(1+delta) ) { 
      h->SetMaximum(1+delta);
      h->SetMinimum(1-delta);

      /// aha!! sutt debug code !!!
      h->DrawCopy();
      char duff[64];
      sprintf(duff,"ratio-%f.pdf", delta );
      gPad->Print(duff);      

      mindelta = delta;
    }

  }

  scale_factor = 1;

  if ( mindelta>2e-3 ) { 
    std::cerr << "WARNING!! poor convolution agreement " << fname << " " << mindelta << std::endl;  

    if ( mindelta>0.05 ) {
      if ( std::fabs(_min-1)>std::fabs(_max-1) ) scale_factor = _min;
      else                                       scale_factor = _max;
    } 

  }

  return h;
}


void GetPdf(const double& x, const double& Q, double* xf) { 
  evolvepdf_( x, Q, xf);    
  //hoppeteval_( x, Q, xf);    
  return; 
}



void processphoton( std::map<std::string,std::vector<int> >&  mLO,
		    std::map<std::string,std::vector<int> >& mNLO );


void photon_subprocesses( appl::grid& g, const std::string& basename, double rscale=1, double fscale=1 );


/// return the head of a full path
std::string head( std::string s ) { 
  std::string h;
  while( s.find("/")!=std::string::npos ) { 
    h+=s.substr(0,s.find("/")+1);
    s.erase(0,s.find("/")+1);
  }
  if ( h.size() ) return h;
  else            return s;
}


/// return the tail of a full path
std::string tail( std::string s ) {
  while( s.find("/")!=std::string::npos ) s.erase(0,s.find("/")+1);
  return s;
}


int usage( int i=0, std::ostream& s=std::cout ){ 
  s << "Usage: standSimple [OPTIONS] -o  input_grid.root [input_grid1.root ... input_gridN.root]\n\n";
  s << "  standSimple performs the convolution for the grids entered on the command line and creates many output plots and a root file\n\n";
  s << "Options: \n";
  s << "    -rs, --rscale value\t  use a renormalisation scale factor of value, \n";
  s << "    -fs, --fscale value\t  use a factorisation scale factor of value, \n";
  s << "    -p,  --pdf value\t  use the the particular set from the specified file,\n";
  s << "    -s  value\t\t  use the particular pdf value from the full set,\n";
  s << "    -d  value\t\t  for a pt or scale dependent x axis variable, use a dynamic scale for the convolution in each bin,\n";
  s << "    -sp value\t\t  plot of the subprocess contributions,\n";
  s << "    -c           \t  apply the multiplicative corrections from the grid,\n";
  s << "    -v, --verbose \t  display verbose output,\n";
  s << "    -h, --help    \t  display this help\n";
  s << "\nSee " << PACKAGE_URL << " for more details\n"; 
  s << "\nReport bugs to <" << PACKAGE_BUGREPORT << ">";
  s << std::endl;
  return i;
}


int main(int argc, char** argv) { 

  //  SetAtlasStyle();

  //  if ( argc<2 ) return -1;

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);


  // setup pdf

  std::string pdfpath = std::string(PDFPATH) + "/";  
  std::string _pdfname = pdfpath + "CT10.LHgrid";  
  // std::string _pdfname = pdfpath + "CT10nlo.LHgrid";  
  // std::string _pdfname = pdfpath + "cteq66.LHgrid";  
  int Npdf = 0;
 

  std::vector<std::string> grids;

  bool prepare_subproc = false;

  double dynamicScale = 0;

  bool corrections = false;

  bool verbose = false;

  double fscale = 1;
  double rscale = 1;


  for ( int ia=1 ; ia<argc ; ia++ ) {  
    if      ( std::string(argv[ia])=="-s" && ia<(argc-1) ) Npdf = std::atoi(argv[++ia]);
    else if ( (std::string(argv[ia])=="-p" || std::string(argv[ia])=="--pdf") && ia<(argc-1) ) _pdfname = pdfpath + argv[++ia];  
    else if ( std::string(argv[ia])=="-d" && ia<(argc-1) ) dynamicScale = std::atof(argv[++ia]);  
    else if ( std::string(argv[ia])=="-sp" ) prepare_subproc = true;  
    else if ( std::string(argv[ia])=="-c"  ) corrections = true;  
    else if ( std::string(argv[ia])=="-v" ||  std::string(argv[ia])=="--verbose" ) verbose = true;
    else if ( (std::string(argv[ia])=="-rs" || std::string(argv[ia])=="--rscale") && ia<(argc-1) ) rscale = std::atof(argv[++ia]);  
    else if ( (std::string(argv[ia])=="-fs" || std::string(argv[ia])=="--fscale") && ia<(argc-1) ) fscale = std::atof(argv[++ia]);  
    else if ( std::string(argv[ia])=="-h" || std::string(argv[ia])=="--help"  ) return usage(0);
    else grids.push_back(argv[ia]);
  }
  
  std::cout << "trying to setup pdf set >" << _pdfname << "<" << " " << Npdf << std::endl; 
  
 
  LHAPDF::initPDFSet(_pdfname.c_str(), Npdf );


  std::cout << "setup pdf" << std::endl; 


#if 0
  { 
    double xf[13];

    for ( double Q2=100 ; Q2<=1e10 ; Q2*=10 ) { 
      for ( double y=-5 ; y<-0.01 ; y+=0.001 ) { 
	double x  = std::pow(10,y);
	evolvepdf_( x, Q2, xf );
	
	std::cout << x << " " << Q2 << " : ";
	for ( int i=0 ; i<13 ; i++ ) std::cout << " " << xf[i];
	std::cout << std::endl; 
      }
    } 
    
  }
#endif

  for ( unsigned ia=0 ; ia<grids.size() ; ia++ ) { 

    std::string gridName(grids[ia]);
    std::string basename  = tail(gridName).substr(0, gridName.find(".root") );
    std::string directory = head(gridName);

    fname = basename;

    try { 
      
    appl::grid g( gridName.c_str() );
    if ( dynamicScale ) g.setDynamicScale( dynamicScale );
    
    bool istrimmed = g.isTrimmed();
    
    g.untrim();
    int gusize = g.size();
    
    g.trim(); 
    
    int gtsize = g.size();
    
    if ( gtsize != gusize ) std::cout << "gridsize : " 
				      << "untrimmed " << gusize 
				      << "\ttrimmed " << gtsize 
				      << "  (" << (gtsize*100./gusize) << "%)\tistrimmed " << istrimmed << std::endl;
    
    if ( verbose ) std::cout << g.getDocumentation() << std::endl;

    
    // get all the reference histograms
    
    TFile* f;
    f = new TFile(grids[ia].c_str());
    
    std::string outFileName = gridName.substr(0,gridName.find(".root"));
    outFileName += std::string("_standSimple") + ( corrections ? "_corr.root" : ".root" );
 
   
    TFile* fout = new TFile(outFileName.c_str(), "recreate");
    
    
    Directory ref("reference");
    ref.push();
    
    TH1D* reference = (TH1D*)f->Get("grid/reference")->Clone("reference");
    
    //  binwidth( reference );
    
    reference->Write();
    ref.pop();


  
    // now calculate all the cross sections
  


    Directory xsDir("xSection");
    xsDir.push();
    
    const int nLoops = g.nloops();
    
    std::cout << "performing convolution ..." << std::endl; 
    

    struct timeval mytimer = appl_timer_start();
    
    TH1D* xsec = g.convolute(evolvepdf_, alphaspdf_, nLoops, rscale, fscale ); 
    
    double mytime = appl_timer_stop(mytimer);

    std::cout << "done (" << mytime << " ms)" << std::endl; 

    //    std::cout << __LINE__ << std::endl; 

    //    std::cout << xsec << std::endl; 

    xsec->SetName("xsec");
    xsec->SetTitle(reference->GetTitle());
    
    //    std::cout << __LINE__ << std::endl; 

    xsec->SetLineColor(kBlue);

    xsec->GetXaxis()->SetMoreLogLabels(true);

    //    std::cout << __LINE__ << std::endl; 

    xsec->Write();  

    //    std::cout << __LINE__ << std::endl; 

    std::vector<TH1D*> xsec_corr; 

    if ( corrections ) { 
      std::cout << "applying corrections" << std::endl;
      xsec_corr.resize( g.corrections().size() );
      for ( unsigned j=0 ; j<g.corrections().size() ; j++ ) { 
	g.setApplyCorrection( j, true );
 
	struct timeval mytimer = appl_timer_start();
    
	xsec_corr[j] = g.convolute(evolvepdf_, alphaspdf_, nLoops, rscale, fscale); 
	
	double mytime = appl_timer_stop(mytimer);

	char hname[64];
	sprintf( hname, "xsec-corr%02d-%s", j, g.correctionLabels()[j].c_str() );

	xsec_corr[j]->SetName( hname );
	xsec_corr[j]->GetXaxis()->SetMoreLogLabels(true);

	xsec_corr[j]->Write();

	std::cout << "done (" << mytime << " ms)" << std::endl; 
       
      }
    }


    if ( verbose ) { 
      std::cout << "xsection has " << xsec->GetNbinsX() << " bins" << std::endl;
      for ( int ih=0 ; ih<xsec->GetNbinsX() ; ih++ ) { 
	std::cout << "\txsec[" << ih << "]";
	if ( ih<10  && xsec->GetNbinsX()>9 )  std::cout << " "; 
	if ( ih<100 && xsec->GetNbinsX()>99 ) std::cout << " "; 
	std::cout << " = " << std::fixed << xsec->GetBinContent(ih+1) << std::endl;
      }
    }    

    xsec->DrawCopy();
    reference->SetLineColor(kBlack);
    reference->SetLineStyle(2);
    reference->GetXaxis()->SetMoreLogLabels(true);

    reference->DrawCopy("samee1");
    
    //    if ( basename.find("pt")!=std::string::npos ) { 
    //      gPad->SetLogy(true);
    //      gPad->SetLogx(true);
    //    } 

    //  gPad->Print("xsec.pdf");
    gPad->Print(("xsec"+basename+".pdf").c_str());
    gPad->Print((directory+"xsec"+basename+".pdf").c_str());
    
    std::cout << "head " << (directory+"xsec"+basename+".pdf") << std::endl;

    gPad->SetLogy(false);
    gPad->SetLogx(false);
    
    // now take all the ratios etc
    
    Directory ratiodir("ratio");
    ratiodir.push();
    
    TH1D* ratio = divide( xsec, reference );

    if ( ratio ) {
      ratio->SetName("ratio");
      ratio->Write();
    }

    ratio->DrawCopy();
    gPad->Print( ("ratio-"+basename+".pdf").c_str() );
    gPad->Print( ("ratio-"+basename+".png").c_str() );
    gPad->Print( (directory+"ratio-"+basename+".pdf").c_str() );

    
    if ( corrections ) { 
      for ( unsigned j=0 ; j<xsec_corr.size() ; j++ ) { 
	TH1D* ratio = divide( xsec_corr[j], reference );
	ratio->SetName( ("ratio-" + g.correctionLabels()[j]).c_str() );
	ratio->Write();

	ratio->DrawCopy();
	gPad->Print( ("ratio-"+basename+g.correctionLabels()[j]+".pdf").c_str() );
	gPad->Print( ("ratio-"+basename+g.correctionLabels()[j]+".png").c_str() );
	gPad->Print( (directory+"ratio-"+basename+g.correctionLabels()[j]+".pdf").c_str() );
      }
    }
     
    
    ratiodir.pop();
    
    std::cout << "individual components..." << std::endl;
    
    TH1D* hLOxsec   = g.convolute(evolvepdf_, alphaspdf_, 0, rscale, fscale); 
    TH1D* hNLOxsec  = 0;
    TH1D* hNLOOxsec = 0;
    
    if ( g.nloops()>0 ) { 
      hNLOxsec  = g.convolute(evolvepdf_, alphaspdf_,  1, rscale, fscale); 
      hNLOOxsec = g.convolute(evolvepdf_, alphaspdf_, -1, rscale, fscale); 
    }

    hLOxsec->SetTitle("");
    hNLOxsec->SetTitle("");
    hNLOOxsec->SetTitle("");

    
    hLOxsec->SetName("xsecLO");
    hNLOxsec->SetName("xsecNLO");
    hNLOOxsec->SetName("xsecNLOO");
    
    hLOxsec->Write();
    hNLOxsec->Write();
    hNLOOxsec->Write();
    
    hLOxsec->GetXaxis()->SetMoreLogLabels(true);
    hNLOxsec->GetXaxis()->SetMoreLogLabels(true);
    hNLOOxsec->GetXaxis()->SetMoreLogLabels(true);
    
    hLOxsec->DrawCopy();
    
    if ( basename.find("pt")!=std::string::npos ) { 
      gPad->SetLogy(true);
      gPad->SetLogx(true);
    }    


    // gPad->Print("LO.pdf");
    gPad->Print(("LO"+basename+".pdf").c_str());
    gPad->Print((directory+"LO"+basename+".pdf").c_str());
    
    if ( hNLOxsec ) {

      reference->GetXaxis()->SetMoreLogLabels(true);

      hNLOxsec->DrawCopy();

      reference->SetLineColor(kBlack);
      reference->SetLineStyle(2);

      reference->DrawCopy("samehist");
      reference->DrawCopy("samee1");
      
      //    gPad->Print("NLO.pdf");
      gPad->Print(("NLO"+basename+".pdf").c_str());
      gPad->Print((directory+"NLO"+basename+".pdf").c_str());

    }    
    
    if ( hNLOOxsec ) {
      hNLOOxsec->DrawCopy();
      // gPad->Print("NLO-only.pdf");
      gPad->Print(("NLO-only"+basename+".pdf").c_str());
      gPad->Print((directory+"NLO-only"+basename+".pdf").c_str());
    }
    
    gPad->SetLogy(false);
    gPad->SetLogx(false);
    
    xsDir.pop();
    
    
    /// write out invididual subprocesses
    
    if ( prepare_subproc && g.getGenpdf().find("photon")!=std::string::npos )   photon_subprocesses( g, basename, rscale, fscale ); 
    
    
    fout->Write();
    fout->Close();
    
    std::cout << "write out ..." << std::endl;


    }    
    catch ( appl::grid::exception e ) { 
      std::cout << "caught exception " << e.what() << " for grid " << gridName << std::endl;
      continue;
    }
  }
  
  return 0;
}






void processphoton( std::map<std::string,std::vector<int> >&  mLO,
		    std::map<std::string,std::vector<int> >& mNLO )
{ 
  
  std::vector<int> v(2);

  v[0] = 2;  v[1] = 5;
  mLO.insert( std::map<std::string,std::vector<int> >::value_type( "ug", v ) ); 
  v[0] = 4;  v[1] = 0;
  mLO.insert( std::map<std::string,std::vector<int> >::value_type( "dg", v ) ); 
  
  v.resize(1);
  v[0] = 1;
  mLO.insert( std::map<std::string,std::vector<int> >::value_type( "dd", v ) ); 
  v[0] = 3;
  mLO.insert( std::map<std::string,std::vector<int> >::value_type( "uu", v ) ); 

  v.resize(4);
  v[0] = 10;  v[1] = 15;  v[2] = 18;  v[3] = 29;
  mNLO.insert( std::map<std::string,std::vector<int> >::value_type( "ug", v ) ); 
  v[0] = 3;  v[1] = 14;  v[2] = 17;  v[3] = 22;
  mNLO.insert( std::map<std::string,std::vector<int> >::value_type( "dg", v ) ); 

  v.resize(3);
  v[0] = 0;  v[1] = 6;  v[2] = 23;
  mNLO.insert( std::map<std::string,std::vector<int> >::value_type( "dd", v ) ); 
  v[0] = 8;  v[1] = 13;  v[2] = 31;
  mNLO.insert( std::map<std::string,std::vector<int> >::value_type( "uu", v ) ); 

  v.resize(5);
  v[0] = 2;  v[1] = 4;  v[2] = 19;  v[3] = 21;  v[4] = 25;
  mNLO.insert( std::map<std::string,std::vector<int> >::value_type( "dd'", v ) ); 
  v[0] = 9;  v[1] = 12;  v[2] = 27;  v[3] = 28;  v[4] = 32;
  mNLO.insert( std::map<std::string,std::vector<int> >::value_type( "uu'", v ) ); 

  v.resize(8);
  v[0] = 1;  v[1] = 5;  v[2] = 7;  v[3] = 11;  v[4] = 20;  v[5] = 24;  v[6] = 26;  v[7] = 30;
  mNLO.insert( std::map<std::string,std::vector<int> >::value_type( "du", v ) ); 

  v.resize(1);
  v[0] = 16;
  mNLO.insert( std::map<std::string,std::vector<int> >::value_type( "gg", v ) ); 
 
}




void photon_subprocesses( appl::grid& g, const std::string& basename, double rscale, double fscale ) { 
  
  std::cout << "photon subprocesses ..." << std::endl;

  // if ( iorder==0 ) nsub = 4;
  // if ( iorder==1 ) nsub = 8;
  
  TH1D* hh = g.convolute( evolvepdf_, alphaspdf_, 1, rscale, fscale);
  //  hh->Rebin(2);


  TH1D* hhf[8];

  std::map<std::string,std::vector<int> >  mLO;
  std::map<std::string,std::vector<int> > mNLO;
  
  processphoton( mLO, mNLO );

  std::string labs[8] = { "ug", "dg", "dd", "uu", "du", "dd'", "uu'", "gg" };

  for ( int ip=0 ; ip<8 ; ip++ ) { 

      std::cout << labs[ip] << " ";
      if ( labs[ip].size()<3 ) std::cout << " ";
    
      TH1D* hlo= 0;
      TH1D* hnlo= 0;

      std::map<std::string,std::vector<int> >::const_iterator  itr = mLO.find(labs[ip]);
      if ( itr!=mLO.end() ) { 
	std::vector<int> vec = itr->second;
	std::cout << "LO ";
	for ( unsigned is=0 ; is<vec.size() ; is++ ) {
	  std::cout << " " << vec[is];
	  if ( hlo==0 ) hlo = g.convolute_subproc( vec[is], evolvepdf_, alphaspdf_, 0, rscale, fscale);
	  else          hlo->Add( g.convolute_subproc( vec[is], evolvepdf_, alphaspdf_, 0, rscale, fscale) );
	} 
	std::cout << "\t";
      }

      itr = mNLO.find(labs[ip]);
      if ( itr!=mNLO.end() ) { 
	std::vector<int> vec = itr->second;
	std::cout << "NLO  ";
	for ( unsigned is=0 ; is<vec.size() ; is++ ) {
	  std::cout << " " << vec[is];
	  if ( hnlo==0 ) hnlo = g.convolute_subproc( vec[is], evolvepdf_, alphaspdf_, -1, rscale, fscale);
	  else           hnlo->Add( g.convolute_subproc( vec[is], evolvepdf_, alphaspdf_, -1, rscale, fscale) );
	} 
	std::cout << std::endl; 
      }
    
      if ( hlo!=0 && hnlo!=0 ) hnlo->Add( hlo );

      //      hnlo->Rebin(2);
      hnlo->SetName( labs[ip].c_str() );
      hnlo->Write();
      
      hnlo->Divide( hh );
      hnlo->SetName( ("f"+labs[ip]).c_str() );

      //      smooth( hnlo );

      hnlo->Write();

      hhf[ip] = hnlo;
  }

  //  hhf[0]->GetXaxis()->SetRangeUser(20, 450);
  hhf[0]->SetMaximum(1);
  hhf[0]->SetMinimum(-0.11);

  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  
  if ( hhf[0]->GetBinCenter(1)>100 ) hhf[0]->SetTitle(";E_{T}^{#gamma} [GeV];parton subprocess fraction");
  else                               hhf[0]->SetTitle(";|#eta^{#gamma}|;parton subprocess fraction");
  
  for ( int ip=0 ; ip<8 ; ip++ ) { 
    hhf[ip]->GetXaxis()->SetMoreLogLabels(true);
    hhf[ip]->SetLineColor( ip<4 ? ip+1 : ip+2 );
    hhf[ip]->SetLineWidth( ip<4 ? 2 : 4 );
    if ( ip==0 ) hhf[ip]->DrawCopy("lhist");
    else         hhf[ip]->DrawCopy("samelhist");

    // DrawLabel(0.8, 0.85-ip*0.05, labs[ip],  ip<4 ? ip+1 : ip+2, 0.038 );
  }

  if ( basename.find("pt")!=std::string::npos ) { 
    //  gPad->SetLogy(true);
    gPad->SetLogx(true);
  }    
  

  // gPad->Print("fractions.pdf"); 
  gPad->Print(("fractions"+basename+".pdf").c_str());
  // gPad->Print((directory+"fractions"+basename+".pdf").c_str());

  gPad->SetLogy(false);
  gPad->SetLogx(false);


}
