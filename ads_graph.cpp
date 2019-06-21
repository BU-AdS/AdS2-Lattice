#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <cstring>
#include <random>
#include <fstream>

using namespace std;
mt19937_64 rng(137);   //pseudo-random 64-bit number generator 
uniform_real_distribution<double> unif;
#define Float long double

#include <util.h>
#include <graph.h>
#include <cg.h>
#include <cg_multishift.h>
#include <eigen.h>

//void latticeScaling(vector<Vertex> &NodeList, Param& p);


// Begin Main Program
//====================================================
int main(int argc, char **argv) 
{
  Param p;    
  if(argc > 1) p.init(argc, argv);   

  p.S1 = endNode(p.Levels,p) + 1 - (endNode(p.Levels-1,p) + 1);
  p.SurfaceVol = p.S1;
  p.AdSVol = (endNode(p.Levels,p) + 1);
  p.latVol = p.AdSVol;
  
  //Print parameters
  p.print();

  //If the specified source position is < 0, place the point source on the outer circumference.
  //if(p.src_pos < 0) p.src_pos = endNode(p.Levels-1,p) + (endNode(p.Levels,p) - endNode(p.Levels-1,p) )/2;  

  int bdry_node_num = endNode(p.Levels,p)-endNode(p.Levels-1,p);
  cout<<"bdry_node_num is: "<<endNode(p.Levels,p)<<"-"<<endNode(p.Levels-1,p)<<"="<<bdry_node_num<<"\n";
  if(p.src_pos < 0) p.src_pos = endNode(p.Levels-1,p)+1;

  //Put source at 6th level for certain calculations
  //p.src_pos = endNode(5,p)+1;
  cout<<"source position is: "<<p.src_pos<<"\n";


  //Setting up initial source locations for 4-point contact term; only the third term will vary
  int* src_array = new int[4];
  src_array[0] = p.src_pos;
  src_array[1] = p.src_pos + (bdry_node_num)/4;
  src_array[2] = p.src_pos + (bdry_node_num)/4 + 1;
  src_array[3] = p.src_pos + 3*(bdry_node_num)/4;

  int outer_configs = src_array[3]-src_array[1];
 
  int src_pos1 = src_array[0];
  int src_pos2 = src_array[1];
  int src_pos3 = src_array[2];
  int src_pos4 = src_array[3];

  //Print graph endnode info
  //for(int i=1; i<20; i++) cout<<"Endnode("<<i<<") = "<<endNode(i,p)<<endl;
  int TotNumber = (endNode(p.Levels,p) + 1);
  
  vector<Vertex> NodeList(TotNumber);
  vector<Vertex> NodeList1(TotNumber);
  vector<Vertex> NodeList2(TotNumber);
  vector<Vertex> NodeList3(TotNumber);
  vector<Vertex> NodeList4(TotNumber);

  //Initialise. -1 in NodeList indicates that node is not yet populated.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < p.q; mu++) { 
      NodeList[n].nn[mu] = -1;
      NodeList1[n].nn[mu] = -1;
      NodeList2[n].nn[mu] = -1;
      NodeList3[n].nn[mu] = -1;
      NodeList4[n].nn[mu] = -1;
    }

  //Construct neighbour table and z-coords
  BuildGraph(NodeList, p);           //populates NodeList[].nn[]
  GetComplexPositions(NodeList, p);  //populates NodeList[].z
  BuildGraph(NodeList1, p);          
  GetComplexPositions(NodeList1, p);  
  BuildGraph(NodeList2, p);           
  GetComplexPositions(NodeList2, p);  
  BuildGraph(NodeList3, p);           
  GetComplexPositions(NodeList3, p);  
  BuildGraph(NodeList4, p);           
  GetComplexPositions(NodeList4, p); 
  
  //for(int i=endNode(p.Levels-1,p)+1; i<TotNumber; i++)
  //NodeList[i].z = NodeList[i].z/abs(NodeList[i].z);

  //Debug tools
  // FIXME: temporal stuff removed
  //ConnectivityCheck(NodeList, p);
  //CheckEdgeLength(NodeList, p);
  //CheckArea(NodeList, p);  
  if(p.verbosity) 
    {
      PrintNodeTables(NodeList, p);  
      PrintComplexPositions(NodeList, p);
    }
  
  if(p.src_pos < 0 || p.src_pos > TotNumber) 
    {
      cout<<"ERROR: Source Position must be g.e. 0 and l.e. "<<TotNumber-1;
      cout<<endl;
      exit(0);
    }
  
  // Lattice Scaling
  latticeScaling(NodeList, p);   
  
  return 0;
}
