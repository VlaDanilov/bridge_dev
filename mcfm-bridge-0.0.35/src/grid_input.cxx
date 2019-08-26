

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <algorithm>

#include "grid_input.h"


 grid_input::grid_input()
 {
   ////////////////////////////////////////////////////////////////
   // execute mcfmbridge-config --version to get mcfmbridge version

   FILE* f = popen("mcfmbridge-config --version","r");
     if ( f == 0 ) {
         fprintf( stderr, "Could not execute 'mcfmbridge-config --version'\n" );
         exit(0);
     }
     // Buffer size for version number (like 0.0.35)
     unsigned short int BUFSIZE = 8;
     char buf[BUFSIZE];
     std::string version;
     while( fgets( buf, BUFSIZE,  f ) )  version = buf;
     pclose( f );
   
   /////////////////////////////////////////////////////////////////
   // Read grid_setup.DAT file

   version.erase(std::remove(version.begin(),version.end(),'\n'),version.end());
   std::string bridge_version = "mcfm-bridge-"+version;
   // The path to mcfm-bridge here is hardcoded from the assumption mcfm-bridge and MCFM are in the same dir
   std::string bridgepath = "../../"+bridge_version+"/src/grid_setup.DAT";
   std::ifstream input_file(bridgepath);
   // check if file exists
   if(!input_file){
   std::cerr << " << ERROR! File 'grid_setup.DAT' does not exists! >> " << std::endl;
   exit(-1);
   }

   /////////////////////////////////////////////////////////////////
   // Loop over lines of the file

   while(std::getline(input_file,str)){
     fakeLine=false;
     if(str.length()==0){fakeLine=true; continue;}
     for(char& c : str) {
       if(c == '\'')  {fakeLine=true; break;   }
       if(c != ' ' && line<14 )   s=s+c;
       else if(line<14){gridSetup[line] = s; s.clear(); break;}
       else if((line>13) && (line<=(13+std::stoi(gridSetup[13])))){
         if(c=='[') continue;
         if(c==']'){customKinematic[line-14][inf]=s; s.clear(); line++; inf=0; break;}
         if(c!=' ') s=s+c;
         else {     customKinematic[line-14][inf]=s; s.clear(); inf++; }
       }
       else {
         if(c=='[') continue;
         if(c==']'){customGrids[line-(14+std::stoi(gridSetup[13]))][inf]=s; s.clear(); line++; inf=0; break;}
         if(c!=' ') s=s+c;
         else { customGrids[line-(14+std::stoi(gridSetup[13]))][inf]=s; s.clear(); inf++; }
       }
     }
     if(line<14 && !fakeLine)line++;
   }
   input_file.close(); 
 }

 std::string grid_input::xLow()
 {
   return gridSetup[0];
 }

 std::string grid_input::xUp()
 {
   return gridSetup[1];
 }

 std::string grid_input::nLoops()
 {
   return gridSetup[2];
 } 

 std::string grid_input::pdfFun()
 {
   return gridSetup[3];
 }

 std::string grid_input::gLab()
 {
   return gridSetup[4];
 }

 std::string grid_input::q2Up()
 { 
   return gridSetup[5];
 }

 std::string grid_input::q2Low()
 {
   return gridSetup[6];
 }

 std::string grid_input::nQ2Bins()
 {
   return gridSetup[7];
 }

 std::string grid_input::qOrder()
 {
   return gridSetup[8];
 }

 std::string grid_input::nXBins()
 {
   return gridSetup[9];
 }

 std::string grid_input::xOrder()
 {
   return gridSetup[10];
 }

 std::string grid_input::lowestOrder()
 {
   return gridSetup[11];
 }

 std::string grid_input::gridNumber()
 {
   return gridSetup[12];
 }

 std::string grid_input::nKinPar()
 {
   return gridSetup[13];
 }

 // customKinematics returns parameters for kinematic cuts: number of particle in process(according to MCFM numbering), value of Pt cut, value of eta cut
 std::string grid_input::customKinematics(int p, int c)
 {
   return customKinematic[p][c];
 }
 
 // customGrid returns observable number(MCFM numbering) and binning
 std::string grid_input::customGrid(int k, int t) // custom grid
 {
   return customGrids[k][t];
 }


