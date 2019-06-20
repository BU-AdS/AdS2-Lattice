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
  ConnectivityCheck(NodeList, p);
  CheckEdgeLength(NodeList, p);
  CheckArea(NodeList, p);  
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
  
#if 1


  
  // Lattice Scaling
  latticeScaling(NodeList, p);   




  /******************************/

  //   for(int iteration=0; iteration < 1; iteration++){
  //     cout<<"Iteration "<<iteration<<"\n";
  //   //----------------//
  //   // (H)MC routines //
  //   //----------------//

  //   //Fake out
  //   //If we want to compute the lattice action, we need a momentum field.  
  //   Float mom[p.AdSVol];
  //   zeroField(mom, p);

  //   int accept = 0;
  //   int count1 = 0;
  
  //   //Initialize Lattice
  //   for(int i = 0; i<p.latVol; i++) {
  //     NodeList[i].phi = 2.0*drand48() - 1.0;
  //     NodeList1[i].phi = 2.0*drand48() - 1.0;
  //     NodeList2[i].phi = 2.0*drand48() - 1.0;
  //     NodeList3[i].phi = 2.0*drand48() - 1.0;
  //     NodeList4[i].phi = 2.0*drand48() - 1.0;
  //   }
  //   //Loop over warmup iterations
  //   for(int iter=0; iter<p.n_therm; iter++) {
  //     //accept = hmc(Lattice, p, iter);    
  //     heatbath(NodeList, p, iter);
  //     heatbath(NodeList1, p, iter);
  //     heatbath(NodeList2, p, iter);
  //     heatbath(NodeList3, p, iter);
  //     heatbath(NodeList4, p, iter);
  //   }
  //   //Loop over warmup iterations
  //   for(int iter=p.n_therm; iter<p.n_therm + p.n_meas; iter++) {

  //     count1++;
  //     double Hold = measH(NodeList, mom, p);
  //     heatbath(NodeList, p, iter);
  //     heatbath(NodeList1, p, iter);
  //     heatbath(NodeList2, p, iter);
  //     heatbath(NodeList3, p, iter);
  //     heatbath(NodeList4, p, iter);
  //     double H = measH(NodeList, mom, p);
  //     //accepted += hmc(Lattice, p, iter);
  //     //cout << iter << " " << setprecision(16) << (1.0*accepted_metropolis)/((iter+1)*p.AdSVol) << " " << H << " " << (H-Hold) << endl; 
  //   }


  //   //-------------//
  //   // CG routines //
  //   //-------------//

  //   //cout<<"In the CG routines"<<"\n";
  //   double phi_prop = 0.0;
  //   double dist = 0.0;
  //   double **phi_prop_array = new double*[TotNumber]; 
  //   double **dist_array = new double*[TotNumber];
  //   for(int i=0; i<TotNumber; i++)
  //     {
  //       phi_prop_array[i] = new double[TotNumber];
  //       dist_array[i] = new double[TotNumber];
  //     }
  //   for(int i=0; i<TotNumber; i++)
  //     {
  //       for(int j=0; j<TotNumber; j++)
  // 	{
  // 	  phi_prop_array[i][j] = 0.0;
  // 	  dist_array[i][j] = 0.0;
  // 	}
  //     }
  
  //   for(int i=0; i<TotNumber; i++)
  //     {
  //       for(int j=0; j<i; j++)
  // 	{
  // 	  //cout<<"i and j are: "<<i<<" "<<j<<"\n";
  // 	  phi_prop_array[i][j] = NodeList[i].phi*NodeList[j].phi;
  // 	  dist_array[i][j] = abs(NodeList[i].z - NodeList[j].z);  
  // 	  //cout<<"NodesList[i].phi is: "<<phi_prop_array[i][j]<<"\n";
  // 	}
  //     }
  
  
  //   cout<<"POOP: "<<phi_prop_array[1][1]<<" "<<dist_array[1][1]<<"\n";
	
  //   ofstream myfile10;
  //   myfile10.open("BrutePhi"+std::to_string(iteration)+".dat");
  //   for(int i=0; i<TotNumber; i++)
  //   {
  //     for(int j=0; j<i; j++)
  //       {
  // 	myfile10<<phi_prop_array[i][j]<<" "<<dist_array[i][j]<<"\n";
  //        }
  //   }
  //   myfile10.close();


  //   //Bulk-boundary propagator
  //   Float* phi0 = new Float[TotNumber];   
  //   Float* b0 = new Float[TotNumber];
  //   for(int i=0; i<TotNumber; i++)
  //     {
  //       phi0[i] = NodeList[i].phi;
  //       b0[i] = 0.0;
  //     }

  //   b0[p.src_pos] = 1;
  //   Float truesq1 = 0.0;
  //   truesq1 = Minv_phi(phi0, b0, NodeList, p);
  //   Bulk2Bdry(NodeList, phi0, p);
  


  //cout<<"After bulk-boundary, before boundary-boundary"<<"\n";
  

  /*  
  //Boundary-boundary propagator
  Float* phibdry = new Float[TotNumber];
  Float* bbdry = new Float[TotNumber];
  for(int i=0; i<TotNumber; i++)
  {
  //phibdry[i] = 0.0;
  phibdry[i] = NodeList[i].phi;
  bbdry[i] = 0.0;
  }
  bbdry[p.src_pos] = 100000000.0;
  Float truesq2 = 0.0;
  truesq2 = Minv_phi(phibdry, bbdry, NodeList, p);
  Bulk2Bdry(NodeList, phibdry, p);
  
  cout<<"After boundary-boundary before radial boundary"<<"\n";
  */
  

//   // ******Laplacian ***********************
//   Float* philap = new Float[TotNumber];
//   Float* blap = new Float[TotNumber];
//   Float** M = new Float*[TotNumber];
//   Float** geo = new Float*[TotNumber];
//   for(int n=0; n<TotNumber; n++)
//     {
//       M[n] = new Float[TotNumber];
//       geo[n] = new Float[TotNumber];
//     }
//   for(int i=0; i<TotNumber; i++)
//     {
//       philap[i] = 0.0;
//       blap[i] = 0.0;
//     }
//   for(int i=0; i<TotNumber; i++)
//     {
//       blap[i] = 1.0;
//       Minv_phi(philap, blap, NodeList, p);
//       for(int j=0; j<TotNumber; j++)
// 	{
// 	  //if(i != j)
// 	  //{
// 	  //Minv_phi(philap, blap, NodeList, p);
// 	      geo[i][j] = d12(NodeList[i].z, NodeList[j].z);
// 	      M[j][i] = philap[j];
// 	      //}
// 	}
//       blap[i] = 0.0;
//     }


//   ofstream wavefiles;
//   wavefiles.open("lap_matrix.dat");
//   for(int i=0; i<p.latVol; i++)
//     {
//       for(int j=0; j<p.latVol; j++)
// 	{
// 	 wavefiles<<M[i][j]<<" ";
// 	}
//       wavefiles<<"\n";
//     }
//   wavefiles.close();


//   ofstream geofiles;
//   geofiles.open("geo_matrix.dat");
//   for(int i=0; i<p.latVol; i++)
//     {
//       for(int j=0; j<p.latVol; j++)
// 	{
// 	  geofiles<<geo[i][j]<<" ";
// 	}
//       geofiles<<"\n";
//     }
//   geofiles.close();
    
  
//   //**********End Laplacian ***********************



  //Radial Bulk-Bdry: characterize as a function of theta
  
  Float* phirad = new Float[TotNumber];
  Float* brad = new Float[TotNumber];
  for(int i=0; i<TotNumber; i++)
    {
      phirad[i] = 0.0;
      //phirad[i] = NodeList[i].phi;
      brad[i] = 0.0;
    }
  brad[0] = 1.0;
  Float truesq = 0.0;
  truesq = Minv_phi(phirad, brad, NodeList, p);
  //Bulk2Bdry(NodeList, phirad, p);
  
  int outer_num = endNode(p.Levels,p)-endNode(p.Levels-1,p)+1;
  Float* phirad1_array = new Float[outer_num];
  /*
    for(int j=0; j<outer_num; j++)
    {
    phirad1_array[j] = phirad[j+endNode(p.Levels-1,p)+1];
 
    Float* phirad1 = new Float[TotNumber];
    Float* brad1 = new Float[TotNumber];
    for(int i=0; i<TotNumber; i++)
    {
    phirad1[i] = 0.0;
    brad1[i] = 0.0;
    }
    brad1[j+endNode(p.Levels-1,p)+1] = 1;
    Float truesq2 = 0.0;
    truesq2 = Minv_phi(phirad1, brad1, NodeList, p);
    //Bulk2Bdry(NodeList, phirad1, p);
    //cout<<"here "<<j<<"\n";
    phirad1_array[j] = phirad1[0];
    delete[] phirad1;
    delete[] brad1;
      
    //cout<<"At number "<<j<<" out of "<<outer_num<<"\n";
    }
  */
  
  
  ofstream myfile1;
  myfile1.open ("phirad.dat");
  complex<long double> one(1,0);
  //for(int i=0; i<TotNumber; i++)
    for(int i=endNode(p.Levels-1,p)+1; i<endNode(p.Levels,p)+1; i++) 
    {
      complex<long double> ratio = NodeList[i].z/NodeList[0].z;
      long double theta = atan2( ratio.imag(), ratio.real());
      //cout<<"i is: "<<i<<"and "<<i-endNode(p.Levels-1,p)-1<<"\n";
      myfile1<<d12(NodeList[0].z, NodeList[i].z)<<" "
	     <<abs((one-(NodeList[i].z)*conj(NodeList[i].z))
		   /((one+NodeList[i].z)*(one+conj(NodeList[i].z))))<<" "
	     <<phirad[i]<<" "
	     <<phirad1_array[i-endNode(p.Levels-1,p)-1]<<" "
	     <<i<<" "
	     <<abs(NodeList[i].z)<<" "
	     <<theta<<" "
	     <<(abs(one)-(long double)abs(NodeList[i].z)*abs(NodeList[i].z))
	/(abs(one)+(long double)(abs(NodeList[i].z)*abs(NodeList[i].z))-2*abs(NodeList[i].z)*cos(3.141592))<<
	" "<<abs(one)/(abs(one)+abs(NodeList[i].z))<<" "<<r(s(NodeList[i].z))<<"\n";
    }   
  myfile1.close();
    
  //1 geo dist b/t center and bdry point
  //2 |1-z*zb|/|(1+z)*(1+zb)|
  //3 phirad[i]
  //4 phirad1_array[i-endNode(p.Levels-1,p)-1]
  //5 i
  //6 |z|
  //7 theta 
  //8 |1-r^2|/(1+r^2-2*r*Cos[theta])
  //9 1/(1+r)
  //10 geodesic dist from origin
 


  //cout<<"After radial bulk-bound"<<"\n";
  

  //cout<<"LATICE SCALING"<<endl;
  //latticeScaling(NodeList, p);



//   //****WAVE EQUATION****
  
//   Float *wave = new Float[p.AdSVol];
//   wave = wave_eqn(NodeList, phirad, p);

//   ofstream wavefile;
//   wavefile.open("wave_eqn.dat");
//   for(int i=0; i<p.latVol; i++)
//     wavefile<<wave[i]<<" "<<phirad[i]<<"\n";
//   wavefile.close();

//   cout<<"After wave"<<"\n";

//   int latVol = p.latVol;
//   int lower = endNode(p.Levels-1,p)+1;

//   Float* phi_ave = new Float[latVol];
//   Float* phi = new Float[latVol];
//   Float* b   = new Float[latVol];
      
//   //Take the average data from all sources (points)
//   int sources = latVol;
//   //int sources =  (endNode(p.Levels,p) - endNode(p.Levels-1,p));

//   //initialise, invert, average
//   for(int i=0; i<latVol; i++) 
//     phi_ave[i] = 0.0;
//   for(int s=0; s<sources; s++) 
//     {
//       for(int i=0; i<latVol; i++) 
//         {
//           b[i] = 0.0;
//           phi[i] = 0.0;
//         }
//       b[s] = 1.0;
//       //b[lower + s] = 1.0;
//       Minv_phi(phi, b, NodeList, p);
//       for(int i=0; i<latVol; i++) 
//         phi_ave[i] += phi[i];
//     }
//   for(int i=0; i<latVol; i++)
//     phi_ave[i] = 0.5 * phi_ave[i] / (latVol*sources);    // double-counted in for loops

//   // phi_ave[n] is thus the average of props from all points to point n
//   // also the same as the average value of phi moving the source around
//   // all points


//   // now that we have this configuration, ask what the wave eq is
//   // to see what the lowest eigenvalue is
//   Float *wave1 = new Float[p.AdSVol];
//   wave1 = wave_eqn(NodeList, phi_ave, p);


//   ofstream wavefile1;
//   wavefile1.open("wave_eqn_1.dat");
//   for(int i=0; i<p.latVol; i++)
//     wavefile1<<wave1[i]<<" "<<phi_ave[i]<<"\n";
//   wavefile1.close();



//   // ***** END WAVE EQUATION ****




  //   //Four-point contact term
  //   //Strategy: loop over configurations. For each configuration, loop over
  //   //which of the four sources CG is called on. 
  //   int temp = src_array[2];
  //   int count = 0;
  //   //Float Amp4 = 0.0;
  //   Float prop = 1.0;
  //   Float* Amp_ana = new Float[outer_configs];
  //   Float* Amp_array = new Float[outer_configs]; 
  //   Float* Amp_pos_array = new Float[outer_configs];
  //   Float* u_array = new Float[outer_configs];
  //   Float* src_pos3_array = new Float[outer_configs];
  //   Float* fourpt_info = new Float[outer_configs];
  //   Float* fourpt_num = new Float[outer_configs];
  //   for(int src_pos3 = temp; src_pos3 < src_array[3]; src_pos3++)
  //     {
  //       //cout<<"src_pos3 is: "<<src_pos3<<" out of "<<src_array[3]<<"\n";
  //       Float xr12 = NodeList[src_pos1].z.real()-NodeList[src_pos2].z.real();
  //       Float xr34 = NodeList[src_pos3].z.real()-NodeList[src_pos4].z.real();
  //       Float xr13 = NodeList[src_pos1].z.real()-NodeList[src_pos3].z.real();
  //       Float xr24 = NodeList[src_pos2].z.real()-NodeList[src_pos4].z.real();
  //       Float xi12 = NodeList[src_pos1].z.imag()-NodeList[src_pos2].z.imag();
  //       Float xi34 = NodeList[src_pos3].z.imag()-NodeList[src_pos4].z.imag();
  //       Float xi13 = NodeList[src_pos1].z.imag()-NodeList[src_pos3].z.imag();
  //       Float xi24 = NodeList[srco_pos2].z.imag()-NodeList[src_pos4].z.imag();
  //       Float u_cr = ((xr12*xr12+xi12*xi12)*(xr34*xr34+xi34*xi34))/((xr13*xr13+xi13*xi13)*(xr24*xr24+xi24*xi24));
  //       Float Dfunction =    (sqrt(M_PI)/(2*(xr13*xr13+xi13*xi13)*(xr24*xr24+xi24*xi24)))*((-2*log(1- sqrt(u_cr)))/(sqrt(u_cr)) + log(u_cr)/(sqrt(u_cr)-1));

  //       Float* prop = new Float[TotNumber-4];
  //       Float* prop1 = new Float[TotNumber-4];    
  //       Float* prop2 = new Float[TotNumber-4];    
  //       Float* prop3 = new Float[TotNumber-4];    
  //       Float* phi = new Float[TotNumber];
  //       Float* phi1 = new Float[TotNumber];
  //       Float* phi2 = new Float[TotNumber];
  //       Float* phi3 = new Float[TotNumber];
  //       Float* b   = new Float[TotNumber];
  //       Float* b1   = new Float[TotNumber]; 
  //       Float* b2   = new Float[TotNumber]; 
  //       Float* b3   = new Float[TotNumber]; 

  //       for(int i=0; i<TotNumber; i++) 
  // 	{
	  
  // 	  phi[i] = 0.0;
  // 	  phi1[i] = 0.0;
  // 	  phi2[i] = 0.0;
  // 	  phi3[i] = 0.0;
  // 	  b[i] = 0.0;
  // 	  b1[i] = 0.0;
  // 	  b2[i] = 0.0;
  // 	  b3[i] = 0.0;
	  
  // 	  if(i<TotNumber-4) 
  // 	    {
  // 	      prop[i] = 0.0;
  // 	      prop1[i] = 0.0;
  // 	      prop2[i] = 0.0;
  // 	      prop3[i] = 0.0;
  // 	    }
  // 	}
      
  //       for(int src_num=0; src_num<4; src_num++)
  // 	{
  // 	  if(src_num == 0)
  // 	    {
  // 	      b[src_pos1] = 100000000;
  // 	      //cout<<"src_num is: "<<src_num<<" and src_array[src_num] is: "<<src_array[src_num]<<"\n";
  // 	      Float truesq = 0.0;
  // 	      truesq = Minv_phi(phi, b, NodeList1, p);   //Residual to this CG solution
  // 	      //cout<<"Tolerance = "<<p.tol<<" True Residual = "<<sqrt(truesq)<<endl;
	      
  // 	      //FourPoint(NodeList, phi, p, src_num, src_array, src_pos3);

  // 	      //cout<<phi[0]<<phi[1]<<"\n";
  // 	      //cout<<(Float)p.N_latt*(Float)phi[0]<<"\n";
  // 	      for(int i=0;i<TotNumber-4;i++)  
  // 		{
  // 		  if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))
  // 		    {
  // 		      prop[i] = (Float)p.N_latt*(Float)phi[i];
  // 		    }     
  // 		}
  // 	      //cout<<prop[0]<<"\n";
  // 	    }
	  
  // 	  if(src_num == 1)
  // 	    {
  // 	      //cout<<"src_num is: "<<src_num<<" and src_array[src_num] is: "<<src_array[src_num]<<"\n";   
  // 	      b1[src_pos2] = 100000000;
  // 	      Float truesq = 0.0;
  // 	      truesq = Minv_phi(phi1, b1, NodeList2, p);
  // 	      //FourPoint(NodeList, phi1, p, src_num, src_array, src_pos3);
  // 	      for(int i=0;i<TotNumber-4;i++) 
  // 		{
  // 		  if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))
  // 		    {  
  // 		      prop1[i] = (Float)p.N_latt*(Float)phi1[i];
  // 		    }
  // 		}
  // 	    }
	  
  // 	  if(src_num == 2)
  // 	    {
  // 	      //cout<<"src_num is: "<<src_num<<" and src_array[src_num] is: "<<src_array[src_num]<<"\n";
  // 	      //cout<<"src_pos2 is: "<<src_pos2<<"\n";
  // 	      b2[src_pos3] = 100000000;
  // 	      Float truesq = 0.0;
  // 	      truesq = Minv_phi(phi2, b2, NodeList3, p);
  // 	      //FourPoint(NodeList, phi2, p, src_num, src_array, src_pos3);
  // 	      for(int i=0;i<TotNumber-4;i++)  
  // 		{
  // 		  if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))  
  // 		    {
  // 		      prop2[i] = (Float)p.N_latt*(Float)phi2[i];
  // 		    }
  // 		}
  // 	    }
  // 	  if(src_num == 3)
  // 	    {
  // 	      //cout<<"src_num is: "<<src_num<<" and src_array[src_num] is: "<<src_array[src_num]<<"\n";   
  // 	      b3[src_pos4] = 100000000;
  // 	      Float truesq = 0.0;
  // 	      truesq = Minv_phi(phi3, b3, NodeList4, p);
  // 	      //FourPoint(NodeList, phi3, p, src_num, src_array, src_pos3);
  // 	      for(int i=0;i<TotNumber-4;i++) 
  // 		{
  // 		  if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))  
  // 		    {
  // 		      prop3[i] = (Float)p.N_latt*(Float)phi3[i];
  // 		    }
  // 		}
  // 	    }
  // 	  //cout<<"AT END "<<"\n";
  // 	  //cout<<prop[0]<<", "<<prop1[0]<<", "<<prop2[0]<<", "<<prop3[0]<<"\n";
  // 	}
  //       Float Amp4 = 0.0;
  //       for(int i=0; i<TotNumber; i++)
  // 	{
  // 	  //cout<<"prop["<<i<<"] is: "<<prop[i]<<"\n";
  // 	  //cout<<"prop1["<<i<<"] is: "<<prop1[i]<<"\n";    
  // 	  //cout<<"prop2["<<i<<"] is: "<<prop2[i]<<"\n";    
  // 	  //cout<<"prop3["<<i<<"] is: "<<prop3[i]<<"\n";    
  // 	  Amp4 += prop[i]*prop1[i]*prop2[i]*prop3[i]*area3q(p.q);
  // 	  //Amp4 += (Float)p.N_latt*(Float)phi[i]*(Float)p.N_latt*(Float)phi1[i]*(Float)p.N_latt*(Float)phi2[i]*(Float)p.N_latt*(Float)phi3[i];
  // 	} 
  //       //cout<<"Amp4 is: "<<Amp4<<"\n";
  //       Amp_array[count] = Amp4;
  //       Amp_pos_array[count] = abs(NodeList[src_pos3].z);
  //       //cout<<"count is: "<<count<<"\n";
  //       //cout<<"u_cr is: "<<u_cr<<"\n";
  //       u_array[count] = u_cr;
  //       Amp_ana[count] = Dfunction;
  //       src_pos3_array[count] = src_pos3;
  
  

  //       /*
  //       for(int i=0; i<TotNumber-4; i++)
  // 	{
  // 	  cout<<"in prop loop"<<"\n";
  // 	  cout<<(Float)p.N_latt*(Float)phi[i]<<"\n";
  // 	  prop += (Float)p.N_latt*(Float)phi[i]*(Float)p.N_latt*(Float)phi1[i]*(Float)p.N_latt*(Float)phi2[i]*(Float)p.N_latt*(Float)phi3[i];  
  // 	}
  //       */
  
  
  //       //fourpt_info[count] = Amp4;
  //       fourpt_num[count] = src_pos3;
  //       count++;

  //       delete[] prop;
  //       delete[] prop1;
  //       delete[] prop2;
  //       delete[] prop3;
  //       delete[] phi;
  //       delete[] phi1; 
  //       delete[] phi2; 
  //       delete[] phi3;
  //       delete[] b;
  //       delete[] b1; 
  //       delete[] b2; 
  //       delete[] b3; 
  //     }
  
  //   //End of 4-pt contact term

  //   cout<<"LOOK HERE: "<<d12(NodeList[0].z, NodeList[1].z)<<"\n";

  //   /*
  //   int n_shift = p.n_shift;
  //   cout<<"n_shift is: "<<n_shift<<"\n";
  //   Float** phi_ms = new Float*[n_shift];
  //   for(int i=0; i<n_shift; i++) {
  //     phi_ms[i] = new long double[TotNumber];
  //     for(int j=0; j<TotNumber; j++) phi_ms[i][j] = 0.0;
  //   }
  //   */
  //   cout<<"outer_configs is: "<<outer_configs<<"\n";

  
  //   ofstream myfile;
  //   myfile.open ("Amp4_HB_"+std::to_string(iteration)+".dat");
  //   for(int i=0; i<outer_configs; i++) 
  //     {
  //       myfile<<u_array[i]<<" "<<Amp_array[i]<<" "<<Amp_ana[i]<<" "<<src_pos3_array[i]<<" "<<Amp_pos_array[i]<<"\n";
  //     }
  //   myfile.close();
 



  //Solve lhs = A^(-1) rhs using multishift CG 
  /*
    Float minbdryradius = GetMinBoundaryRadius(NodeList,p);
    Float maxbdryradius = GetMaxBoundaryRadius(NodeList,p);
    cout<<"Minimum boundary radius is: "<<minbdryradius<<"\n";
    cout<<"Maximum boundary radius is: "<<maxbdryradius<<"\n";
  */

  /*
    ofstream myfile1;
    myfile1.open ("phirad_heatbath"+std::to_string(iteration)+".dat");
    complex<long double> one(1,0);
    for(int i=endNode(5,p)+1; i<endNode(6,p)+1; i++) 
    {
    complex<long double> ratio = NodeList[i].z/NodeList[0].z;
    long double theta = atan2( ratio.imag(), ratio.real());
    //cout<<"i is: "<<i<<"and "<<i-endNode(p.Levels-1,p)-1<<"\n";
    myfile1<<d12(NodeList[0].z, NodeList[i].z)<<" "
    <<abs((one-(NodeList[i].z)*conj(NodeList[i].z))
    /((one+NodeList[i].z)*(one+conj(NodeList[i].z))))<<" "
    <<phirad[i]<<" "
    <<phirad1_array[i-endNode(p.Levels-1,p)-1]<<" "
    <<i<<" "
    <<abs(NodeList[i].z)<<" "
    <<theta<<" "
    <<(abs(one)-(long double)abs(NodeList[i].z)*abs(NodeList[i].z))
    /(abs(one)+(long double)(abs(NodeList[i].z)*abs(NodeList[i].z))-2*abs(NodeList[i].z)*cos(3.141592))<<
    " "<<abs(one)/(abs(one)+abs(NodeList[i].z))<<"\n";
    }   
    myfile1.close();
    
    //1 geo dist b/t center and bdry point
    //2 |1-z*zb|/|(1+z)*(1+zb)|
    //3 phirad[i]
    //4 phirad1_array[i-endNode(p.Levels-1,p)-1]
    //5 i
    //6 |z|
    //7 theta 
    //8 |1-r^2|/(1+r^2-2*r*Cos[theta])
    //9 1/(1+r)
    //10

    delete[] phirad1_array;

  */

  /*
    ofstream myfile2;
    myfile2.open("phibdry_heatbath"+std::to_string(iteration)+".dat");
    for(int i=endNode(p.Levels-1,p)+1; i<endNode(p.Levels,p)+1; i++)
    {
    complex<long double> ratio = NodeList[i].z/NodeList[p.src_pos].z;
    long double theta = atan2(ratio.imag(), ratio.real());
    myfile2<<d12(NodeList[p.src_pos].z, NodeList[i].z)<<" "
    <<theta<<" "
    <<phibdry[i]<<" "
    <<abs(NodeList[i].z)<<" "
    <<NodeList[i].z.real()<<" "
    <<NodeList[i].z.imag()<<" "
    <<NodeList[p.src_pos].z.real()<<" "
    <<NodeList[p.src_pos].z.imag()<<"\n";
    }
    myfile2.close();
  
    delete[] phibdry;
    //1 geodesic distance
    //2 theta
    //3 lattice bound-bound prop
    //4 |z|
    //5 Re[z]
    //6 Im[z]

    */

  //}


  //Write file of Phi[X] vs d12(Phi[X], Source[src_pos]), X on boundary
  /*
  
    ofstream myfile1;
    myfile1.open ("phirad_Lvl6_Lvl9.dat");
    complex<long double> one(1,0);
    for(int i=endNode(5,p)+1; i<endNode(6,p)+1; i++) 
    {
    complex<long double> ratio = NodeList[i].z/NodeList[0].z;
    long double theta = atan2( ratio.imag(), ratio.real());
    //cout<<"i is: "<<i<<"and "<<i-endNode(p.Levels-1,p)-1<<"\n";
    myfile1<<d12(NodeList[0].z, NodeList[i].z)<<" "
    <<abs((one-(NodeList[i].z)*conj(NodeList[i].z))
    /((one+NodeList[i].z)*(one+conj(NodeList[i].z))))<<" "
    <<phirad[i]<<" "
    <<phirad1_array[i-endNode(p.Levels-1,p)-1]<<" "
    <<i<<" "
    <<abs(NodeList[i].z)<<" "
    <<theta<<" "
    <<(abs(one)-(long double)abs(NodeList[i].z)*abs(NodeList[i].z))
    /(abs(one)+(long double)(abs(NodeList[i].z)*abs(NodeList[i].z))-2*abs(NodeList[i].z)*cos(3.141592))<<
    " "<<abs(one)/(abs(one)+abs(NodeList[i].z))<<"\n";
    }   
    myfile1.close();
    
    //1 geo dist b/t center and bdry point
    //2 |1-z*zb|/|(1+z)*(1+zb)|
    //3 phirad[i]
    //4 phirad1_array[i-endNode(p.Levels-1,p)-1]
    //5 i
    //6 |z|
    //7 theta 
    //8 |1-r^2|/(1+r^2-2*r*Cos[theta])
    //9 1/(1+r)
    //10

  
    ofstream myfile2;
    myfile2.open ("phirad_Lvl5_Lvl9.dat");
    //complex<long double> one(1,0);
    for(int i=endNode(5,p)+1; i<endNode(6,p)+1; i++) 
    {
    complex<long double> ratio = NodeList[i].z/NodeList[0].z;
    long double theta = atan2( ratio.imag(), ratio.real());
    //cout<<"i is: "<<i<<"and "<<i-endNode(p.Levels-1,p)-1<<"\n";
    myfile2<<d12(NodeList[0].z, NodeList[i].z)<<" "
    <<abs((one-(NodeList[i].z)*conj(NodeList[i].z))
    /((one+NodeList[i].z)*(one+conj(NodeList[i].z))))<<" "
    <<phirad[i]<<" "
    <<phirad1_array[i-endNode(p.Levels-1,p)-1]<<" "
    <<i<<" "
    <<abs(NodeList[i].z)<<" "
    <<theta<<" "
    <<(abs(one)-(long double)abs(NodeList[i].z)*abs(NodeList[i].z))
    /(abs(one)+(long double)(abs(NodeList[i].z)*abs(NodeList[i].z))-2*abs(NodeList[i].z)*cos(3.141592))<<
    " "<<abs(one)/(abs(one)+abs(NodeList[i].z))<<"\n";
    }   
    myfile2.close();
  
  
    ofstream myfile3;
    myfile3.open ("phirad_Lvl6_Lvl9.dat");
    //complex<long double> one(1,0);
    for(int i=endNode(6,p)+1; i<endNode(7,p)+1; i++) 
    {
    complex<long double> ratio = NodeList[i].z/NodeList[0].z;
    long double theta = atan2( ratio.imag(), ratio.real());
    //cout<<"i is: "<<i<<"and "<<i-endNode(p.Levels-1,p)-1<<"\n";
    myfile3<<d12(NodeList[0].z, NodeList[i].z)<<" "
    <<abs((one-(NodeList[i].z)*conj(NodeList[i].z))
    /((one+NodeList[i].z)*(one+conj(NodeList[i].z))))<<" "
    <<phirad[i]<<" "
    <<phirad1_array[i-endNode(p.Levels-1,p)-1]<<" "
    <<i<<" "
    <<abs(NodeList[i].z)<<" "
    <<theta<<" "
    <<(abs(one)-(long double)abs(NodeList[i].z)*abs(NodeList[i].z))
    /(abs(one)+(long double)(abs(NodeList[i].z)*abs(NodeList[i].z))-2*abs(NodeList[i].z)*cos(3.141592))<<
    " "<<abs(one)/(abs(one)+abs(NodeList[i].z))<<"\n";
    }   
    myfile3.close();


    ofstream myfile4;
    myfile4.open ("phirad_Lvl7_Lvl9.dat");
    //complex<long double> one(1,0);
    for(int i=endNode(7,p)+1; i<endNode(8,p)+1; i++) 
    {
    complex<long double> ratio = NodeList[i].z/NodeList[0].z;
    long double theta = atan2( ratio.imag(), ratio.real());
    //cout<<"i is: "<<i<<"and "<<i-endNode(p.Levels-1,p)-1<<"\n";
    myfile4<<d12(NodeList[0].z, NodeList[i].z)<<" "
    <<abs((one-(NodeList[i].z)*conj(NodeList[i].z))
    /((one+NodeList[i].z)*(one+conj(NodeList[i].z))))<<" "
    <<phirad[i]<<" "
    <<phirad1_array[i-endNode(p.Levels-1,p)-1]<<" "
    <<i<<" "
    <<abs(NodeList[i].z)<<" "
    <<theta<<" "
    <<(abs(one)-(long double)abs(NodeList[i].z)*abs(NodeList[i].z))
    /(abs(one)+(long double)(abs(NodeList[i].z)*abs(NodeList[i].z))-2*abs(NodeList[i].z)*cos(3.141592))<<
    " "<<abs(one)/(abs(one)+abs(NodeList[i].z))<<"\n";
    }   
    myfile4.close();
  

    ofstream myfile5;
    myfile5.open ("phirad_Lvl8_Lvl9.dat");
    //complex<long double> one(1,0);
    for(int i=endNode(8,p)+1; i<TotNumber; i++) 
    {
    complex<long double> ratio = NodeList[i].z/NodeList[0].z;
    long double theta = atan2( ratio.imag(), ratio.real());
    //cout<<"i is: "<<i<<"and "<<i-endNode(p.Levels-1,p)-1<<"\n";
    myfile5<<d12(NodeList[0].z, NodeList[i].z)<<" "
    <<abs((one-(NodeList[i].z)*conj(NodeList[i].z))
    /((one+NodeList[i].z)*(one+conj(NodeList[i].z))))<<" "
    <<phirad[i]<<" "
    <<phirad1_array[i-endNode(p.Levels-1,p)-1]<<" "
    <<i<<" "
    <<abs(NodeList[i].z)<<" "
    <<theta<<" "
    <<(abs(one)-(long double)abs(NodeList[i].z)*abs(NodeList[i].z))
    /(abs(one)+(long double)(abs(NodeList[i].z)*abs(NodeList[i].z))-2*abs(NodeList[i].z)*cos(3.141592))<<
    " "<<abs(one)/(abs(one)+abs(NodeList[i].z))<<"\n";
    }   
    myfile5.close();
  

    ofstream myfile6;
    myfile6.open ("phirad_Lvl9_Lvl9.dat");
    //complex<long double> one(1,0);
    for(int i=endNode(9,p)+1; i<TotNumber; i++) 
    {
    complex<long double> ratio = NodeList[i].z/NodeList[0].z;
    long double theta = atan2( ratio.imag(), ratio.real());
    //cout<<"i is: "<<i<<"and "<<i-endNode(p.Levels-1,p)-1<<"\n";
    myfile5<<d12(NodeList[0].z, NodeList[i].z)<<" "
    <<abs((one-(NodeList[i].z)*conj(NodeList[i].z))
    /((one+NodeList[i].z)*(one+conj(NodeList[i].z))))<<" "
    <<phirad[i]<<" "
    <<phirad1_array[i-endNode(p.Levels-1,p)-1]<<" "
    <<i<<" "
    <<abs(NodeList[i].z)<<" "
    <<theta<<" "
    <<(abs(one)-(long double)abs(NodeList[i].z)*abs(NodeList[i].z))
    /(abs(one)+(long double)(abs(NodeList[i].z)*abs(NodeList[i].z))-2*abs(NodeList[i].z)*cos(3.141592))<<
    " "<<abs(one)/(abs(one)+abs(NodeList[i].z))<<"\n";
    }   
    myfile6.close();
  */
  

  /*

    ofstream myfile2;
    myfile2.open("phibdry_Lvl6_Lvl9.dat");
    for(int i=endNode(5,p)+1; i<endNode(6,p)+1; i++)
    {
    complex<long double> ratio = NodeList[i].z/NodeList[p.src_pos].z;
    long double theta = atan2(ratio.imag(), ratio.real());
    myfile2<<d12(NodeList[p.src_pos].z, NodeList[i].z)<<" "
    <<theta<<" "
    <<phibdry[i]<<" "
    <<abs(NodeList[i].z)<<" "
    <<NodeList[i].z.real()<<" "
    <<NodeList[i].z.imag()<<" "
    <<NodeList[p.src_pos].z.real()<<" "
    <<NodeList[p.src_pos].z.imag()<<"\n";
    }
    myfile2.close();
  */
  //1 geodesic distance
  //2 theta
  //3 lattice bound-bound prop
  //4 |z|
  //5 Re[z]
  //6 Im[z]
  
  /*
    cout<<"before phibdry_Lvl9.dat"<<"\n";
    ofstream myfile9;
    myfile9.open("phibdry_Lvl9.dat");
    for(int i=endNode(8,p)+1; i<TotNumber; i++)
    {
    complex<long double> ratio = NodeList[i].z/NodeList[p.src_pos].z;
    long double theta = atan2(ratio.imag(), ratio.real());
    myfile9<<d12(NodeList[p.src_pos].z, NodeList[i].z)<<" "
    <<theta<<" "
    <<phibdry[i]<<" "
    <<abs(NodeList[i].z)<<" "
    <<NodeList[i].z.real()<<" "
    <<NodeList[i].z.imag()<<" "
    <<NodeList[p.src_pos].z.real()<<" "
    <<NodeList[p.src_pos].z.imag()<<"\n";
    }
    myfile9.close();
  */
  //1 geodesic distance
  //2 theta
  //3 lattice bound-bound prop
  //4 |z|
  //5 Re[z]
  //6 Im[z]

 
  //delete[] phi0;
  //delete[] b0;
  //delete[] phi_ms;
  //delete[] phirad;
  //delete[] brad;

#endif

  
  
  return 0;
}



// void latticeScaling(vector<Vertex> &NodeList, Param& p)
// {
//   int latVol  = p.latVol;
//   int lower = endNode(p.Levels-1,p)+1;
//   long unsigned int TotNumber = (endNode(p.Levels,p) + 1) * p.t;
//   int outer_cirum = endNode(p.Levels,p)-endNode(p.Levels-1,p);

//   int t1, t2, delta_t;
//   double delta = 1.0 + sqrt(1 + p.msqr);
//   double theta, r, r_p;
//   complex<double> ratio;
//   complex<double> snk;
//   double* analytic_prop = new double[outer_cirum*p.Lt/2];
//   double* xi_invariant  = new double[outer_cirum*p.Lt/2];

 
//   // Construct the xi invariant and the analytic propagator
  
//   int j = lower;
//   complex<Float> src = NodeList[j].z;  
//   //Loop over timeslices
//   for(int t=0; t<p.Lt/2; t++) {
//     int T_offset = (endNode(p.Levels,p) + 1) * t;

//     //Loop over outer circle of the Poincare disk
//     for(long unsigned int k = 0; k < outer_cirum; k++) {

//       //Construct i (full lattice index)
//       int i = T_offset + lower + k;
//       //Construct s (surface lattice index)
//       int s = t*outer_cirum + k;
      
//       ratio = NodeList[i].z/NodeList[j].z;
//       theta = atan2( ratio.imag() , ratio.real() );
//       //index divided by disk size, using the int floor feature/bug,
//       //gives the timeslice for each index.
//       int t1 = j / (TotNumber/p.t);
//       int t2 = i / (TotNumber/p.t);

//       //Assume PBC.
//       int delta_t = (t2-t1);
//       complex<Float> snk = NodeList[i].z;
//       Float r   = abs(NodeList[i].z);
//       Float r_p = abs(NodeList[j].z);
      
//       Float xi_invariant  = log( ((1-r)*(1-r_p))/(cosh(delta_t)*(1+r)*(1+r_p)
// 					       - 4*r*r_p*cos(theta)) );
      
//       Float analytic_prop = log( exp(-delta*sigma(src,snk,delta_t)) /
// 			      (1 - exp(-2*sigma(src,snk,delta_t))));
      
//       //if(s<10) cout<<"s="<<s<<" xi="<<xi_invariant[s]<<" ap="<<analytic_prop[s]<<endl;
//     }
//   }

//   double* c = new double[2];
//   double* cov_ssq = new double[4];
//   gsl_fit_linear(xi_invariant, 1, analytic_prop, 1, outer_cirum*p.Lt/2, &c[0], &c[1],
// 		 &cov_ssq[0], &cov_ssq[1], &cov_ssq[2], &cov_ssq[3]);

//   cout<<"Target data"<<endl;
//   cout<<"GSL data: C="<<c[0]<<" M="<<c[1]<<endl;
//   cout<<"          covar00 = "<<cov_ssq[0]<<endl;
//   cout<<"          covar01 = "<<cov_ssq[1]<<endl;
//   cout<<"          covar11 = "<<cov_ssq[2]<<endl;
//   cout<<"          sum_sq  = "<<cov_ssq[3]<<endl;

//   double grad = c[1];
//   double d_grad = 100;

//   double grad_tol = 1e-4;

//   //Search in wisdom file for a shorcut
//   bool preTuned = false;
//   bool tuneWin = false;
//   ifstream fileInput;
//   string line;
//   char params[256];
//   sprintf(params, "%d %d %d %.4f", p.q, p.Levels, p.Lt, p.msqr);
//   char* search = params;
//   unsigned int curLine = 0;
//   // open file to search
//   fileInput.open("ads_wisdom");
//   if(fileInput.is_open()) {
//     while(getline(fileInput, line)) { 
//       curLine++;
//       if (line.find(search, 0) != string::npos) {
// 	cout<<"Found data for your problem! Lucky you..."<<endl;
// 	preTuned = true;
// 	getline(fileInput, line);
// 	p.N_latt = stof(line);
// 	getline(fileInput, line);
// 	p.C_msqr = stof(line);
// 	cout<<"Using N_latt="<<p.N_latt<<endl;
// 	cout<<"Using C_msqr="<<p.C_msqr<<endl;
// 	fileInput.close();
//       }
//     }
//     if(!preTuned)
//       cout<<endl<<"AdS wisdom data not found. Strap in for some tuning..."<<endl<<endl;
//   }
//   else cout<<endl<<"AdS wisdom file not found. Strap in for some tuning..."<<endl<<endl;

//   //Begin search for correct scaling factors.
//   int iter = 0;
//   //while(abs(d_grad) > grad_tol || abs(d_inter) > inter_tol) {
//   while(abs(d_grad) > grad_tol) {
    
//     double* phi_ave = new double[latVol];  
//     double* phi = new double[latVol];
//     double* b   = new double[latVol];

//     //Take the average data from sources in the the qth sector
//     //of the outer level. 
//     int sources = (endNode(p.Levels,p) - endNode(p.Levels-1,p)) / p.q ;

//     //initialise, invert, average.
//     for(int i=0; i<latVol; i++) phi_ave[i] = 0.0;    
//     for(int s=0; s<sources; s++) {
//       for(int i=0; i<latVol; i++) {
// 	b[i] = 0.0;
// 	phi[i] = 0.0;
//       }
//       b[lower + s] = 1.0;    
//       Minv_phi(phi, b, NodeList, p);
//       for(int i=0; i<latVol; i++) phi_ave[i] += phi[i]/sources;
//     }
    
//     //Use current lattice normalisation.
//     for(int i=0; i<latVol; i++) {
//       phi_ave[i] = log(p.N_latt*phi_ave[i]);
//     }
    
//     //phi_ave now contains an averaged solution vector. We now
//     //perform linear regression on this vector and the analytic
//     //prop (w.r.t the xi invariant) to compute the relevant
//     //scaling and normalsation.
    
//     double* latt_prop  = new double[outer_cirum*p.Lt/2];
//     //Loop over timeslices
//     for(int t=0; t<p.Lt/2; t++) {
//       int T_offset = (endNode(p.Levels,p) + 1) * t;
//       //Loop over H2 disk
//       for(long unsigned int k = 0; k < outer_cirum; k++) {
	
// 	//Construct i (full lattice AdS2p1 index)
// 	int i = T_offset + lower + k;
// 	//Construct s (surface lattice 2D index)
// 	int s = t*outer_cirum + k;
	
// 	latt_prop[s] = phi_ave[i];
//       }
//     }
  
//     //Extract linear fit data from log-log plot.
//     gsl_fit_linear(xi_invariant, 1, latt_prop, 1, outer_cirum*p.Lt/2, &c[0], &c[1],
// 		   &cov_ssq[0], &cov_ssq[1], &cov_ssq[2], &cov_ssq[3]);
    
//     cout<<"At iteration "<<iter<<endl;
//     cout<<"GSL data: C="<<c[0]<<" M="<<c[1]<<endl;
//     cout<<"          covar00 = "<<cov_ssq[0]<<endl;
//     cout<<"          covar01 = "<<cov_ssq[1]<<endl;
//     cout<<"          covar11 = "<<cov_ssq[2]<<endl;
//     cout<<"          sum_sq  = "<<cov_ssq[3]<<endl;
    
//     //Adjust the parameters and start over.
//     d_grad = (c[1] - grad)/abs(grad);
//     cout<<"D grad = "<<d_grad<<endl;
//     cout<<"C_msqr = "<<p.C_msqr<<endl<<endl;
    
//     double fac1 = abs(p.C_msqr);
    
//     if( abs(d_grad) > grad_tol ) {
//       if(c[1] > grad) p.msqr < 0 ? p.C_msqr += (fac1*abs(d_grad) + 1*grad_tol) : p.C_msqr -= (fac1*abs(d_grad) + 1*grad_tol);
//       if(c[1] < grad) p.msqr < 0 ? p.C_msqr -= (fac1*abs(d_grad) + 1*grad_tol) : p.C_msqr += (fac1*abs(d_grad) + 1*grad_tol);
//     }
    
//     if(cov_ssq[0] != cov_ssq[0]) {
//       cout<<"GSL failed!"<<endl;
//       exit(0);
//     }
    
//     //Did it tune properly?
//     if(d_grad < grad_tol) tuneWin = true;
    
//     delete[] phi_ave;
//     delete[] phi;
//     delete[] b;
//     delete[] latt_prop;

//     iter++;
    
//   }
  
//   //If this is a new problem, and it tuned properly, save the data.
//   if(!preTuned && tuneWin) {
//     FILE *file = fopen("ads_wisdom", "a");
//     fprintf(file,"%d %d %d %.4f\n%f\n%f\n",
// 	    p.q, p.Levels, p.Lt, p.msqr,
// 	    p.N_latt,
// 	    p.C_msqr);
//     fclose(file);
//   }

//   delete[] c;
//   delete[] cov_ssq;

//   delete[] analytic_prop;
//   delete[] xi_invariant;
// }



