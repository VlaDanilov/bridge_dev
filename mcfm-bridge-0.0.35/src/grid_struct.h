#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

class GridStruct
{
  private:
        bool fakeLine;
        int  line=0, inf=0;
        std::string s, str, obsNames; 
  public:

        std::string gridSetup[42]; //the answer to The Ultimate Question of Life, the Universe, and Everything
        std::string customGrids[42][184];
        std::string customKinematic[42][184];
        GridStruct()
        {
          bool fakeLine;
          int  line=0, inf=0;
          std::string s, str, obsNames;
	  std::ifstream input_file("/nfs/dust/cms/user/volinar/tools/Test_AREA/mcfm-bridge-0.0.35/src/grid_setup.txt");
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
                       //std::cout << "line number = " << line << std::endl;
                       //std::cout <<"<<<< "<<c<<" >>>>"<<std::endl; 
                       if(c=='[') continue;
                       if(c==']'){customGrids[line-(14+std::stoi(gridSetup[13]))][inf]=s; s.clear(); line++; inf=0; break;}
                       if(c!=' ') s=s+c;
                       else { customGrids[line-(14+std::stoi(gridSetup[13]))][inf]=s; s.clear(); inf++; }
                       //std::cout <<"<"<<customGrids[line-15][0]<<">"<<std::endl; 
                       //std::cout <<"<"<<customGrids[line-15][1]<<">"<<std::endl; 
                       //std::cout <<"<"<<customGrids[line-15][2]<<">"<<std::endl; 
                       //std::cout <<"<"<<customGrids[line-15][3]<<">"<<std::endl; 
                       //std::cout <<"<"<<customGrids[line-15][4]<<">"<<std::endl; 
                       ///if(c!=char(0))  s=s+c ;
                       ///else{ customGrids[line-13+inf]=s; s.clear(); inf++; }
                     }
            }
            if(line<14 && !fakeLine)line++;
          }
          input_file.close();
        }
         ~GridStruct(){}

 std::string xLow()
 { 
   return gridSetup[0];
 }
 std::string xUp()
 {
   return gridSetup[1];
 }
 std::string nLoops()
 {
   return gridSetup[2];
 }
 std::string pdf_Fun()
 {
   return gridSetup[3];
 }
 std::string gLab()
 {
   return gridSetup[4];
 }
 std::string q2Up()
 {
   return gridSetup[5];
 }
 std::string q2Low()
 {
   return gridSetup[6];
 }
 std::string nQ2Bins()
 {
   return gridSetup[7];
 }
 std::string qOrder()
 {
   return gridSetup[8];
 }
 std::string nXBins()
 {
   return gridSetup[9];
 }
 std::string xOrder()
 {
   return gridSetup[10];
 }
 std::string LowestOrder()
 {
   return gridSetup[11];
 }
 int gridNumber()
 {
   return std::stoi(gridSetup[12]);
 }
 std::string nKinPar()
 {
   return gridSetup[13];
 }
 std::string customKinematics(int p, int c) // custom grid
 {    
   return customKinematic[p][c];
 }
 std::string customGrid(int k, int t) // custom grid
 {
   return customGrids[k][t];
 } 








};
