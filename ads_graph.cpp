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

  //Print graph endnode info
  //for(int i=1; i<20; i++) cout<<"Endnode("<<i<<") = "<<endNode(i,p)<<endl;
  int TotNumber = (endNode(p.Levels,p) + 1);
  
  vector<Vertex> NodeList(TotNumber);

  //Initialise. -1 in NodeList indicates that node is not yet populated.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < p.q; mu++) { 
      NodeList[n].nn[mu] = -1;
    }

  //Construct neighbour table and z-coords
  BuildGraph(NodeList, p);           //populates NodeList[].nn[]
  GetComplexPositions(NodeList, p);  //populates NodeList[].z
  
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
  //latticeScaling(NodeList, p);   

  
  //-------------//
  // CG routines //
  //-------------//

 
  /*
  Float* phi0 = new Float[TotNumber];   
  Float* b0 = new Float[TotNumber];
  for(int i=0; i<TotNumber; i++)
    {
      phi0[i] = 0.0;
      b0[i] = 0.0;
    }

  b0[p.src_pos] = 1;
  Float truesq1 = 0.0;
  truesq1 = Minv_phi(phi0, b0, NodeList, p);
  Bulk2Bdry(NodeList, phi0, p);
  */
  
  Mphi_ev(NodeList, p);

  /*
  Float** phi1 = new Float*[TotNumber];
  Float* b1 = new Float[TotNumber];
  for(int i=0; i<TotNumber; i++)
    {
      b1[i] = 0.0;
      phi1[i] = new Float[TotNumber];
      for(int j=0; j<TotNumber; j++)
	{
	  phi1[i][j] = 0.0;
	}
    }

  for(int i=0; i<TotNumber; i++)
    {
      b1[i] = 1.0;
      Minv_phi(phi1[i], b1, NodeList, p);
      b1[i] = 0.0;
    }

  */


  ofstream distfile;
  distfile.open("distance_m"+to_string(p.msqr)+"_L"+to_string(p.Levels)+"_q"+to_string(p.q)+".dat");
  for(int i=0; i<TotNumber; i++)
    {
      for(int j=0; j<TotNumber; j++)
	{
	  distfile << d12(NodeList[i].z, NodeList[j].z) << " ";
	}
      distfile << "\n";
    }
  distfile.close();

  return 0;
}
