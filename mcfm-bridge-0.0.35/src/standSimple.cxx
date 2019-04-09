#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>

#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"


#include "appl_grid/appl_timer.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TPad.h>

#include "LHAPDF/LHAPDF.h"
#include "hoppet_v1.h"

#define DBG false


#include "scales.h"


extern "C" void evolvepdf_(const double& x, const double& Q, double *xf);
extern "C" double alphaspdf_(const double& Q);

static const int nLoops    = 1;
static const int nFlavours = 5;


void GetPdf(const double& x, const double& Q, double* xf) { 
   double _x_(x);
   if (_x_<1e-7) _x_=1e-7;
  evolvepdf_( x, Q, xf);
  xf[0]=0.0; xf[12]=0.0;
  for(int i=0;i<13;i++)if(xf[i]>1e10)xf[i]=0;
  //hoppeteval_( x, Q, xf);    
  return; 
}


TH1D* divide( const TH1D* h1, const TH1D* h2 ) {

  if ( h1==NULL || h2==NULL ) return NULL;
 
  TH1D* h = (TH1D*)h1->Clone();

  if ( DBG ) std::cout << "histograms " << h1->GetTitle() << " " << h2->GetTitle() << std::endl;

  

  for ( int i=1 ; i<=h1->GetNbinsX() ; i++ ) { 
    double b  = h2->GetBinContent(i);
    double be = h2->GetBinError(i);
    double t  = h1->GetBinContent(i);
    double te = h1->GetBinError(i);

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
  
  if ( h->GetMaximum()<1.01 ) h->SetMaximum(1.01);
  if ( h->GetMinimum()>0.99 ) h->SetMinimum(0.99);

  return h;
}

// extern double Escale;

int main(int argc, char** argv) { 

  if ( argc<2 ) return -1;

  TString gridName(argv[1]);
 

  // <---GRID---> 
  appl::grid g(gridName.Data());
  g.trim();
  //


  std::cout << g << std::endl;

  // get all the reference histograms

  TFile* f;
  f = new TFile(argv[1]);

  TString outFileName = gridName(0,gridName.Index(".root"));
  outFileName += "_standSimple.root";

  std::cout <<" *** \t \t *** : Input file = "<<gridName<<" output file = "<<outFileName<<std::endl;
  TFile* fout = new TFile(outFileName, "recreate");


  Directory ref("reference");
  ref.push();

  TH1D* reference = (TH1D*)f->Get("grid/reference");
//  TH1D* referenceScale[Nscales];

//  for ( int i=0 ; i<Nscales ; i++  ) {
//    char hname[64];
//    sprintf( hname, "addReference/NreferenceScale_%d", i); 
//    referenceScale[i] = (TH1D*)f->Get(hname);
//    referenceScale[i]->Write();
//  }

  reference->Write();
  ref.pop();


  
  // now calculate all the cross sections

    const std::string _pdfname = "CT14nnlo_as_0118";
//    const std::string _pdfname = "CT10nlo";
//    const std::string _pdfname = "cteq66";
  int Npdf = 0;
  // setup gavins code
  LHAPDF::initPDFSet(_pdfname.c_str(), Npdf);
//  initPDF(Npdf);

  Directory xsDir("xSection");
  xsDir.push();
  
//  TH1D* xsec = g.convolute(evolvepdf_, alphaspdf_, nLoops); 
  TH1D* xsec = g.convolute(GetPdf, alphaspdf_, nLoops); 
  xsec->SetName("xsec");
  xsec->SetTitle(reference->GetTitle());

/*  TH1D* xsec_scale[Nscales];
  for ( int i=0; i < Nscales ; i++ ) 
    { 
      char hname[64];
      sprintf( hname, "xsec_scale_%d", i); 
      
      xsec_scale[i] = g.convolute(GetPdf, alphasPDF , nLoops, mur[i], muf[i]);
      xsec_scale[i]->SetName(hname);
      xsec_scale[i]->SetTitle(referenceScale[i]->GetTitle());
      xsec_scale[i]->Write();
    }
    */
  xsec->Write();  
  xsDir.pop();

  // now take all the ratios etc

  Directory ratiodir("ratio");
  ratiodir.push();

  TH1D* ratio = divide( xsec, reference ); 
  if ( ratio ) {
    ratio->SetName("ratio");
    ratio->Write();
  }
  
/*  TH1D* ratio_scale[Nscales];

  for ( int i=0; i < Nscales ; i++ ) 
    { 
      
      char hname[64];
      sprintf( hname, "ratio_scale_%d", i); 
      
      ratio_scale[i] = divide( xsec_scale[i], referenceScale[i] );
      if ( ratio_scale[i] ) 
	{
	  ratio_scale[i]->SetName(hname);
	  ratio_scale[i]->Write();
	}
    }
  */
  ratiodir.pop();

  fout->Write();
  fout->Close();

  return 0;
}
