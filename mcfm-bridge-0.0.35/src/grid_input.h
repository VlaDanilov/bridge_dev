
#ifndef GRID_INPUT_H
#define GRID_INPUT_H

#include <string>

class grid_input
{
  bool fakeLine;
  int  line=0, inf=0;
  std::string s, str, obsNames;

public:

  // Not good, array's are constant!
  std::string gridSetup[42]; 
  std::string customGrids[42][184];
  std::string customKinematic[42][184];
 
  grid_input();

  ~grid_input(){}
 
  std::string xLow();
  std::string xUp();
  std::string nLoops();
  std::string pdfFun();
  std::string gLab();
  std::string q2Up();
  std::string q2Low();
  std::string nQ2Bins();
  std::string qOrder();
  std::string nXBins();
  std::string xOrder();
  std::string lowestOrder();
  std::string gridNumber();

  // returns quantity of particles with kinematic(pt, eta) cuts
  std::string nKinPar(); 

  // customKinematics returns parameters for kinematic cuts: number of particle in process(according to MCFM numbering), value of Pt cut, value of eta cut
  std::string customKinematics(int p, int c);

  // customGrid returns observable number(MCFM numbering) and binning
  std::string customGrid(int k, int t);
  
};

#endif


