//
//   @file    lumi_pdf.cxx         
//   
//
//   @author M.Sutton
// 
//   Copyright (C) 2013 M.Sutton (sutt@cern.ch)    
//
//   $Id: lumi_pdf.cxx, v0.0   Tue  9 Jul 2013 08:14:47 CEST sutt $


#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>



std::string str_replace( std::string s ) { 
  std::string a = s;
  for ( int i=a.size() ; i-- ; ) if ( a[i]=='_' ) a[i]='-';
  return a;
}

#include "appl_grid/lumi_pdf.h"

void latex( const lumi_pdf& p, const std::string& d);

bool lumi_pdf::m_runlatex = false;

lumi_pdf::lumi_pdf(const std::string& s, const std::vector<int>& combinations ) : // , int Wcharge ) :  //, bool amcflag ) : 
  appl_pdf(s), m_filename(s), 
  m_lookup( std::vector<std::vector<std::vector<int> > >(0) ) 
  //,  m_amcflag(amcflag)
{
  
  /// need to decode the input std::vector

  if ( combinations.size() ) { 

    /// std::vector initialised from serialised std::vector
    unsigned iv = 0;

    unsigned nproc = combinations[iv++];

    for ( unsigned i=0 ; i<nproc && iv<combinations.size() ; i++ ) {
      int index  = combinations[iv++];
      int npairs = combinations[iv++];
      std::vector<int> v(npairs*2+2);
      v[0] = index;
      v[1] = npairs;
      
      for ( int j=0 ; j<npairs ; j++ ) { 
	v[j*2+2] = combinations[iv++];
	v[j*2+3] = combinations[iv++];
      }

      combination c(v);
      if ( c.size() ) add(c); 
    } 

    /// extra value on the end for the W+- charge (if required) - flags that 
    /// ckm matrix is to be used    
    if ( iv<combinations.size() ) m_ckmcharge = combinations[iv];

  }
  else if ( m_filename!="" ) {  
    /// else read from file ...
  
    std::ifstream& infile = openpdf( m_filename  );

    /// will never fail, appl_pdf::open() would have thrown already
    ///  if ( infile.fail() ) throw exception( std::cerr << "lumi_pdf::lumi_pdf() cannot open file " << m_filename << std::endl ); 
 
    std::string   line;

    infile >> m_ckmcharge;
    
    ///    std::cout << "ckmcharge " << m_ckmcharge << std::endl;  

    while (std::getline(infile, line)) {
      //    std::cout << "line: " << line << std::endl;
      combination c( line );
      if ( c.size() ) add(c); 
    }
    
    infile.close();
  }

  if ( m_ckmcharge>0 ) { 
    std::cout << "lumi_pdf::lumi_pdf() setting W+ cmk matrix" << std::endl;
    make_ckm(true);
  }
  else if ( m_ckmcharge<0  ) { 
    std::cout << "lumi_pdf::lumi_pdf() setting W- cmk matrix" << std::endl;
    make_ckm(false);
  }

  
  // some checking
  
  //  for ( int i=0 ; i<m_combinations.size() ; i++ ) { 
  //    if ( m_combinations[i].
  //  }

  m_Nproc = m_combinations.size();

  std::cout << "lumi_pdf::lumi_pdf()         " << s << "\tcombinations " << size() << std::endl;  


  // create the reverse lookup 

  create_lookup();

  //  std::cout << "decideSuprocess " << decideSubProcess( 0, 0 ) << std::endl;
  //  std::cout << "lumi_pdf::lumi_pdf() " << s << "\tv size " << m_combinations.size() << " lookup size " << m_lookup.size() << std::endl; 
  //  std::cout << *this << std::endl;

  //  lumi_pdf* _pdf = dynamic_cast<lumi_pdf*>(appl::appl_pdf::getpdf(name()));
  //  std::cout << "done " << _pdf << _pdf->decideSubProcess( 0, 0 ) << std::endl;

  //  std::cout << *this << std::endl;

}





lumi_pdf::lumi_pdf(const std::string& s, const std::vector<combination>& combinations, int ckmcharge ) : 
  appl_pdf(s), m_filename(s), m_combinations(combinations),
  m_lookup( std::vector<std::vector<std::vector<int> > >(0) ) 
{
  
  //  std::cout << "lumi_pdf::lumi_pdf() " << s << "\tv size " << combinations.size() << " lookup size " << m_lookup.size() << " " << this << std::endl; 

  /// no need need to decode the input std::vector

  /// the W+- charge (if required) - flags that 
  /// ckm matrix is to be used    
  m_ckmcharge = ckmcharge;

  if ( m_ckmcharge>0 ) { 
    std::cout << "lumi_pdf::lumi_pdf() setting W+ cmk matrix" << std::endl;
    make_ckm(true);
  }
  else if ( m_ckmcharge<0  ) { 
    std::cout << "lumi_pdf::lumi_pdf() setting W- cmk matrix" << std::endl;
    make_ckm(false);
  }

  m_Nproc = m_combinations.size();

  create_lookup();

  //  std::cout << *this << std::endl; 
  
}



void lumi_pdf::create_lookup() { 
  if ( m_lookup.size()==0 ) { 
    /// create a 14 x 14 lookup table (including a photon contribution! 
    m_lookup = std::vector<std::vector<std::vector<int> > >(14, std::vector<std::vector<int> >(14) ); 
    for ( unsigned i=size() ; i-- ; ) { 
      const combination& c = m_combinations[i];
      for ( unsigned j=c.size() ; j-- ; ) m_lookup[ c[j].first+6 ][ c[j].second+6 ].push_back(i);
    } 
  }

  m_proclookup.clear();

  for ( unsigned i=size() ; i-- ; ) { 
    const combination& c = m_combinations[i];
    for ( unsigned j=c.index().size() ; j-- ; ) { 
      std::map<int,int>::iterator itr = m_proclookup.find(c.index()[j]);
      if ( itr==m_proclookup.end() ) m_proclookup.insert( std::map<int,int>::value_type( c.index()[j], i ) );
    }
  } 
  
}



int lumi_pdf::decideSubProcess(const int iproc ) const {
  std::map<int,int>::const_iterator itr = m_proclookup.find(iproc);
  if ( itr==m_proclookup.end() ) return -1;
  else return itr->second;
}


void lumi_pdf::evaluate(const double* xfA, const double* xfB, double* H) const { 
  /// if need to include the ckm matrix ...
  if ( m_ckmcharge==0 )  {
    for ( unsigned i=size() ; i-- ; ) { 
      H[i] = m_combinations[i].evaluate( xfA, xfB ); 
    }
  }
  else { 
   for ( unsigned i=size() ; i-- ; ) { 
     H[i] = m_combinations[i].evaluate( xfA, xfB, m_ckmsum, m_ckm2 ); 
    }
  }
}


int  lumi_pdf::decideSubProcess(const int iflav1, const int iflav2) const { 
  //  std::cout << "lumi_pdf::decideSubProcess() " << name() << " " << m_lookup.size() << std::endl;
  if ( m_lookup[iflav1+6][iflav2+6].size()==1 ) return m_lookup[iflav1+6][iflav2+6][0];
  else                                          return -1;
}


size_t  lumi_pdf::nSubProcesses(const int iflav1, const int iflav2) const { 
  //  std::cout << "lumi_pdf::decideSubProcess() " << name() << " " << m_lookup.size() << std::endl;
  return m_lookup[iflav1+6][iflav2+6].size();
}


std::vector<int> lumi_pdf::decideSubProcesses(const int iflav1, const int iflav2) const { 
  //  std::cout << "lumi_pdf::decideSubProcess() " << name() << " " << m_lookup.size() << std::endl;
  return m_lookup[iflav1+6][iflav2+6];
}


std::vector<int> lumi_pdf::serialise() const  { 

  std::vector<int> v;

  v.push_back( Nproc() );

  for ( int i=0 ; i<Nproc() ; i++ ) { 

    const combination& c = m_combinations[i];

    v.push_back( c.index()[0] );

    v.push_back( c.size() );
    for ( unsigned j=0 ; j<c.size() ; j++ ) { 
      v.push_back( c[j].first );
      v.push_back( c[j].second );
    }
  }  
  
  if      ( m_ckmcharge>0 ) v.push_back(1);
  else if ( m_ckmcharge<0 ) v.push_back(-1);
  else                      v.push_back(0);
			 
  return v;
}



std::vector<std::vector<int> > lumi_pdf::vectorise() const  { 

  std::vector<std::vector<int> > v;

  for ( int i=0 ; i<Nproc() ; i++ ) { 

    std::vector<int> v0;
    v0.push_back( i );

    const combination& c = m_combinations[i];
    for ( unsigned j=0 ; j<c.size() ; j++ ) { 
      v0.push_back( c[j].first );
      v0.push_back( c[j].second );
    }
    v.push_back( v0 );
  }  
  
  return v;
}





void lumi_pdf::write(std::ostream& s) const { 
  
  s << m_ckmcharge << "\n";

  for ( unsigned i=0 ; i<m_combinations.size() ; i++ ) { 
    s << m_combinations[i].index() << " ";
    s << m_combinations[i].size()  << " ";

    for ( unsigned j=0 ; j<m_combinations[i].size() ; j++ ) { 
      s << "  " << m_combinations[i][j].first << " " << m_combinations[i][j].second;
    }

    s << "\n";

  }

}


void lumi_pdf::write(const std::string& filename) const {  
  std::ofstream s(filename.c_str());
  write(s);
}


// std::string lumi_pdf::summary(std::ostream& s=std::cout) const { 
std::string lumi_pdf::summary() const { 
  std::stringstream s_;
  s_ << "lumi_pdf::lumi_pdf()\t" << name() << "\tcombinations " << m_combinations.size() << "\tlookup size " << m_lookup.size() << "\taddr: " << this; 
  return s_.str();
}



void lumi_pdf::removeDuplicates() { 

  std::cout << "lumi_pdf::removeDuplicates() " << name() << "\tsize " << size() << " -> ";

  std::vector<combination> combinations;

  for ( unsigned i=0 ; i<size() ; i++ ) {
    bool unique = true;
    for ( unsigned j=0 ; j<combinations.size() ; j++ ) {
      if ( at(i) == combinations[j] ) {
	unique = false;
	// std::cout << "duplicate:\n" << at(i) << "\n" << at(j) << std::endl;
	combinations[j].add_index( at(i).index()[0] );
      }
    }
    if ( unique ) combinations.push_back( at(i) );
  }

  m_combinations = combinations;

  m_Nproc = m_combinations.size();

  create_lookup();

  if ( m_runlatex ) latex( *this, ".pdf" );

  std::cout << size() << std::endl; 

}
   




void lumi_pdf::restoreDuplicates() { 

  std::vector<combination> combinations;

  for ( unsigned i=0 ; i<size() ; i++ ) {

    std::vector<int> indices = at(i).index();

    for ( unsigned j=0 ; j<indices.size() ; j++ ) {

      combination c = at(i);

      c.index().clear();
      c.index().push_back(indices[j]);

      combinations.push_back( c );
    }
  }

  std::sort( combinations.begin(), combinations.end() );

  m_combinations = combinations;

  m_Nproc = m_combinations.size();

  create_lookup();
  
}


   



#include <fstream>
#include <cstdlib>

/// decode as a latex table, and create pdf file

void latex( const lumi_pdf& p, const std::string& d) {

  std::ofstream zj( (d+p.name()+".tex").c_str() );
 
  std::string _f[13] = { 
    "\\bar{t}", 
    "\\bar{b}", 
    "\\bar{c}", 
    "\\bar{s}", 
    "\\bar{u}", 
    "\\bar{d}",
    "g",
    "d",
    "u", 
    "s", 
    "c", 
    "b", 
    "t" };

  const std::string* f = _f+6; 


  //  std::cout << "\\bigskip\n";

  std::cout << "Contribution:   " << p.name() << "\t processes " << p.Nproc() << "\n";

  zj << "\\documentclass[7pt,a4paper,landscape]{article}\n\n";

  zj << "\\usepackage{makecell}\n\n";
  zj << "\\usepackage[a4paper,landscape]{geometry}\n\n";

  zj << "\\topmargin=-3cm\n";
  zj << "\\textheight=18cm\n";
  zj << "\\textwidth=28cm\n";
  zj << "\\oddsidemargin=-2.2cm\n";
  

  zj << "\\begin{document}\n";

  zj << "\\ \\\\ \\ \\\\ \\ \\\\\n";

  zj << "{\\footnotesize\n";


  size_t maxproc = 0;

  size_t totalproc = 0;

  for ( int i=0 ; i<p.Nproc() ; i++ ) { 
    const combination& c = p[i];

    const std::vector<int>& ind = c.index();

    if ( ind.size()>maxproc ) maxproc = ind.size();

    totalproc += ind.size();
  }  

  size_t maxpairs = 20-maxproc;


  zj << "\\hspace{-12cm}";
  zj << "\\begin{minipage}[t]{18cm}\n";
  zj << "pdf : " << str_replace(p.name()) << "\tnprocesses: " << totalproc  << "\\\\" << std::endl; 
  zj << "\\begin{tabular}{cll}\\hline\\\\\n";


  for ( int i=0 ; i<p.Nproc() ; i++ ) { 

    const combination& c = p[i];

    const std::vector<int>& ind = c.index();

    int nrows = (ind.size()+27)/28;

    zj << "\\makecell{ " << i;
    for ( int ig=1 ; ig<nrows ; ig++ ) zj << " \\\\ \\ ";
    zj << "}\t&\t";


    zj << "\\makecell[l]{";
    // zj << "\\thead[l]{";

    for ( unsigned j=0 ; j<ind.size() ; j++ ) {
      zj << ind[j] << " ";
      if ( (j+1)%20==0 ) zj << "\\\\";
    }


    zj << "}\t&\t";


    zj << "\\makecell[l]{ ";
 
    for ( unsigned j=0 ; j<c.size() ; j++ ) {

      int p0 = c[j].first;
      int p1 = c[j].second;
      
      if ( j>0 ) zj << " + ";

      if ( (j+1)%maxpairs==0 ) zj << " \\\\ \\ ";

      zj << "($" << f[p0] << "$, $" << f[p1] << "$)\t";  
      
    }
    for ( int ig=1 ; ig<nrows ; ig++ ) zj << " \\\\ \\ ";
    zj << "}";

    zj << "\\\\\n";
    
  }

  zj << "\\\\\n";
  zj << "\\hline\n";
  zj << "\\end{tabular}\n";
  zj << "\\end{minipage}\n";

  zj << "}\n";
  zj << "\\end{document}\n";

  zj.close();

  std::system( (std::string("pdflatex ")+ d+p.name()+".tex").c_str() );

}




