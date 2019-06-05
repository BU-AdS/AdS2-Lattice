#ifndef UTIL_H
#define UTIL_H
#include <complex>
#include <cstring>
#include <math.h>


//#include <util.h>
//#include <graph.h>
#include "cg.h"
//#include <cg_multishift.h>
//#include <eigen.h>
//#include <hmc_util.h>


//#include "gsl_fit.h"

using namespace std;

#define I complex<Float>(0.0,1.0)



class Param{

 public:

  int q = 7;
  
  bool bc = true;       //if true, use Dirichlet. If false, use Neumann
  bool Vcentre = true;  //if true, place vertex at centre. If false, use circumcentre.
  bool verbosity = false;  //if true, print all data. If false, print summary.
  int MaxIter = 100000;
  int n_shift = 100;
  Float tol = pow(10,-6);
  int t = 1;
  Float msqr = 0.1;
  Float lambda = 1.0;
  Float C_msqr = 1.0;
  Float N_latt = 1.0;
  int Levels = 3;
  int src_pos = -1;
  Float DiskScale = 1.0;
  char fname[256];

  int Lt = 32;
  int S1 = 32;
  int SurfaceVol = 0;
  int AdSVol = 0;
  int latVol = 0;
  // double lambda = 1.0;
  //double musqr = 1.0;
  int *cluster ;    // Swendsen Wang Data Struture
  int *stack ;     // Wolf Data Struture
  int NumClusters ;

  //HMC 
  int n_metro_cool = 0;
  int n_therm = 1000;
  int n_meas = 1000;
  int n_write = 100;
  int n_skip = 100;
  int n_cluster = 8;
  int n_jkblock = 10;
  int n_step = 10;
  double tau = 1.0;
  double dt = 0.1;
  double delta_phi = 1.5;

  
  
  void print(){
    cout<<"Parameter status:"<<endl;
    cout<<"Triangulation = "<<q<<endl;
    cout<<"B.C. = "<< (bc ? ("Dirichlet") : ("Neumann") ) << endl;
    cout<<"Centre = "<< (Vcentre ? ("Vertex") : ("Circum") ) << endl;
    cout<<"MaxIter = "<<MaxIter<<endl;
    cout<<"Number of Shifts = "<<n_shift<<endl;
    cout<<"Tol = "<<tol<<endl;
    cout<<"TimeSlices = "<<t<<endl;   
    cout<<"Mass squared = "<<msqr<<endl;
    cout<<"lambda = "<<lambda<<endl;
    cout<<"Levels = "<<Levels<<endl;
    cout<<"Source Position = "<<src_pos<<endl;
    cout<<"Mass squared Correction = "<<C_msqr<<endl;
    cout<<"Lattice normalisation = "<<N_latt<<endl;
    cout<<"DiskScale = "<<DiskScale<<endl;
  }
  
  void init(int argc, char **argv) {
    
    std::string BC(argv[1]);
    if (BC == "D" || BC == "d") {
      bc = true;
    } else if (BC == "N" || BC == "n") {
      bc = false;
    } else {
      cout<<"Invalid boundary condition given. Use D/d for Dirichlet or N/n for Neumann."<<endl;
      exit(0);
    }
    
    std::string Centre(argv[2]);
    if (Centre == "V" || Centre == "v") {
      Vcentre = true;
    } else if (Centre == "C" || Centre == "c") {
      Vcentre = false;
    } else {
      cout<<"Invalid centre condition given. Use V/v for Vertexcentred or C/c for Circumcentred."<<endl;
      exit(0);
    }

    std::string verbose(argv[3]);
    if (verbose == "V" || verbose == "v") {
      verbosity = true;
    } else if(verbose == "Q" || verbose == "q") {
      verbosity = false;
    } else {
      cout<<"Invalid Verbosity conditions given. Use verbose/quiet"<<endl;
      exit(0);
    }

    MaxIter = atoi(argv[4]);
    tol     = atof(argv[5]);
    t       = atoi(argv[6]);    
    msqr    = atof(argv[7]);
    lambda  = atof(argv[8]);
    Levels  = atoi(argv[9]);
    src_pos = atoi(argv[10]);
    
    //if(atof(argv[11]) == 0) C_msqr = -0.0126762/msqr + 0.0689398*msqr + 2.02509;
    if(atof(argv[11]) == 0) C_msqr = -0.0126762/msqr + 0.0689398*msqr + 2.02509;
    else C_msqr = atof(argv[11]);
    
    //if(atof(argv[12]) == 0) N_latt = 0.294452/(msqr + 0.766901) + 0.0788137;
    if(atof(argv[12]) == 0) N_latt = 0.294452/(msqr + 0.766901) + 0.0788137;
    else N_latt = atof(argv[12]);
    
    q = atoi(argv[13]);
    n_shift = atoi(argv[14]);
    
  }
};

class Vertex{
 public:
  //int nn[11] = {0,0,0,0,0,0,0,0,0,0,0};
  //int fwdLinks;
  complex<Float> z;
  double temporal_weight = 1.0;

  //If the pos value is -1, it is not connected
  //to the graph.
  int pos = -1;
  
  //Nearest neighbours for up to q=9 and 2 temporal directions.
  int nn[11] = {0,0,0,0,0,0,0,0,0,0,0};

  //How many forward links (important in the
  //buildGraph() function.
  int fwdLinks;

  //Positon on the Poincare disk.
  //std::complex<double> z;
  
  //Phi field value at this vertex.
  double phi = 0.0;
  
  //Old phi field value at this vertex (HMC).
  double phi_old = 0.0;
  
  //Ising field value at this vertex.
  int ising = 0.0;



};


typedef vector<Vertex> Graph;

complex<Float> T(complex<Float> z,  complex<Float> w);
complex<Float> R(complex<Float> z, complex<Float> omega);
complex<Float> flip(complex<Float> z, complex<Float> z1, complex<Float> z2);
Float s(complex<Float> z);
Float r(Float s );
Float d12(complex<Float> z1, complex<Float> z2);
Float sigma(complex<Float> z1, complex<Float> z2, int t);
Float s3p(int q);
Float area3q(int q);
Float areaGeneral(Param P, Float A, Float B, Float C);
Float centralRad(Float s);
complex<Float> DisktoUHP(complex<Float> z);
complex<Float> UHPtoDisk(complex<Float> u);
complex<Float> inversion(complex<Float> z0, Float r);
complex<Float> squareInversion(complex<Float>z0, Float r1, Float r2 );
Float greens2D(complex<Float> z, complex<Float> w);
Float greensM2D(complex<Float> z, complex<Float> w, Param p);
complex<Float> newVertex(complex<Float> z,complex<Float> z0,int k, int q);

void PrintNodeTables(const vector<Vertex> NodeList, Param P);

//- Edge length from center z = 0
Float edgeLength(int q) {
  return sqrt( 1 - 4*sin(M_PI/q)*sin(M_PI/q) );
}

//Using the formula c(n) = (q-4)*c(n-1) - c(n-2) where c is the
//number of nodes on circumference at level n, we can construct
//the address of the end node on a given level for triangulation q:
//EN(lev,q) = SUM c(n) n=0..level
long unsigned int endNode(int lev, Param P) { 
  
  int q = P.q;
  
  //This is the case where a vertex is at the centre of the disk
  if(P.Vcentre == true) {
    //Explicit results for level <= 2.
    if(lev==0) return 0;
    if(lev==1) return q;
    if(lev==2) return q*(q-4) + q;
    
    //level >= 3
    int p=q-4; //A convenience  
    long unsigned int c[20]; //Array to hold the circumference info
    
    //initialise array
    for(int i=0; i<20; i++) c[i] = 0;
    //manually set 0,1,2 entries:
    c[0] = 0;
    c[1] = q;
    c[2] = p*q;
    
    long unsigned int EN = 0; //address of node at end of level;  
    //Loop up to lev
    for(int n=3; n<lev+1; n++) {
      c[n] = p*c[n-1] - c[n-2];
      EN += c[n];
    }
    EN += (c[1] + c[2]);
    return EN;
  }

  //This is the case where a circumcentre is at the centre of the disk.
  else {
    //Explicit results for level <= 1.
    if(lev==0) return 2;
    if(lev==1) return (q-3)*3 + 2;
    
    //level >= 2
    int p=q-4; //A convenience  
    long unsigned int c[20]; //Array to hold the circumference info
    
    //initialise array
    for(int i=0; i<20; i++) c[i] = 0;
    //manually set 0,1,2 entries:
    //NB!!! There are three nodes on level 0, but the END NODE is addressed
    //as 2. Therefore, when correcting the number of nodes on a
    //circumference, we use 3 (there are three nodes)
    //but when giving the endNode count, we use 2 (we count from 0)
    c[0] = 3;       
    c[1] = (q-3)*3;
    
    long unsigned int EN = 0; //address of node at end of level;  
    //Loop up to lev
    for(int n=2; n<lev+1; n++) {
      c[n] = p*c[n-1] - c[n-2];
      EN += c[n];
    }
    EN += (2 + c[1]); //0,1,2 nodes to add from circumference 0
    return EN;
  }
}

//- Get the z coordinates of every node on the Poincare disk 
void GetComplexPositions(Graph &NodeList, Param& P){

  int q = P.q;
  int Levels = P.Levels;
  int T_offset = endNode(P.Levels,P)+1;
  
  if(P.Vcentre == true) {
    //Assume for now that the origin (level 0) is a vertex
    NodeList[0].z = 0.0;
    //Assert that node 1 is on the real axis
    complex<Float> init(edgeLength(q),0.0);
    NodeList[1].z = init;
    //Rotate to create level level 1
    for(int k=1; k<q+1; k++) {
      NodeList[k].z = newVertex(init, 0.0, k-1, q);
    }
    //For every node on level >=1, the zeroth
    //nn is the (n-1)th link on the same level. If
    //n=1, the node address is the endnode value.
    for(int l=1; l<Levels+1; l++) {
      for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
	for(int k=0; k<q; k++) {
	  if(NodeList[n].nn[k] != -1) {
	    NodeList[NodeList[n].nn[k]].z = newVertex(NodeList[NodeList[n].nn[0]].z, NodeList[n].z, k, q);
	    //NodeList[NodeList[n].nn[k]].temporal_weight = (1+pow(abs(NodeList[NodeList[n].nn[k]].z),2))/(1-pow(abs(NodeList[NodeList[n].nn[k]].z),2));
	    //NodeList[NodeList[n].nn[k]].temporal_weight = 1.0 / pow(cos(M_PI*abs(NodeList[NodeList[n].nn[k]].z)/2),1);
	  }
	}
      }
    }    
  }
  else {
    
    Float numer = sqrt(cos(M_PI*(q+6)/(6*q)) - sin(M_PI/q));
    Float denom = sqrt(sin(M_PI/q) + sin(M_PI*(q+3)/(3*q)));    
    Float init_mod = sqrt(norm(numer/denom));
    
    //Assume that node 0 lies on +ve real axis
    complex<Float> init_0(init_mod,0.0);
    NodeList[0].z = init_0;
    
    //Assert that node 1 is node 0 rotated by 2*PI/3
    complex<Float>init_1(init_mod*cos(2.0*M_PI/3.0),
			  init_mod*sin(2.0*M_PI/3.0));
    
    NodeList[1].z = init_1;
    //Rotate node 1 about node 0 to create level 0 (the equilateral triangle)
    NodeList[2].z = newVertex(init_1, init_0, 1, q);

    //For every node on level >=1, the zeroth
    //nn is the (n-1)th link on the same level. If
    //n=1, the node address is the endnode value.
    for(long unsigned int n=0; n<endNode(1,P)+1; n++) {
      for(int k=0; k<q; k++) {
	if(NodeList[n].nn[k] != -1) {
	  NodeList[NodeList[n].nn[k]].z = newVertex(NodeList[NodeList[n].nn[0]].z, NodeList[n].z, k, q);
	}
      }
    }
    for(int l=1; l<Levels+1; l++) {
      for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
	for(int k=0; k<q; k++) {
	  if(NodeList[n].nn[k] != -1) {
	    NodeList[NodeList[n].nn[k]].z = newVertex(NodeList[NodeList[n].nn[0]].z, NodeList[n].z, k, q);
	    //NodeList[NodeList[n].nn[k]].temporal_weight = (1+pow(abs(NodeList[NodeList[n].nn[k]].z),2))/(1-pow(abs(NodeList[NodeList[n].nn[k]].z),2));
	  }
	}
      }
    }
  }

  if(P.t > 1) {
    //Copy all 2D complex positions along the cylinder
    for(long unsigned int n=0; n<endNode(P.Levels,P)+1; n++) 
      for(int t=1; t<P.t; t++) NodeList[n + T_offset*t].z = NodeList[n].z;
  }
}

//- For each node n, with a link to another node,
//  it checks that the neighbour table on the linked
//  node contains the original node n as a neighbour.
void ConnectivityCheck(Graph &NodeList, Param P){

  int q = P.q;
  int Levels = P.Levels;
  int T = P.t;
  int TotNumber = T*(endNode(Levels,P)+1);
  int t_offset  = 0;
  T == 1 ? t_offset = 0 : t_offset = 2;
  
  //Object to hold boolean values of graph connectivity.
  vector<Vertex> AuxNodeList(TotNumber);
  //Initialise to 0.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < q+t_offset; mu++) {
      AuxNodeList[n].nn[mu] = 0;
    }
  
  for(long unsigned int n=0; n<TotNumber; n++) {
    for(int m=0; m<q+t_offset; m++) {
      //Check that the link is valid
      if(NodeList[n].nn[m] != -1) {
	for(int p=0; p<q+t_offset; p++) {
	  //Loop over all links on the linked node,
	  //check if original node exists in neighbour
	  //table.
	  if( n == NodeList[ NodeList[n].nn[m] ].nn[p] ) {
	    AuxNodeList[n].nn[m] = 1;
	  }
	}
      }
    }
  }

  if(P.verbosity) PrintNodeTables(AuxNodeList, P);
}

void PrintNodeTables(const vector<Vertex> NodeList, Param P) {

  int q = P.q;
  int Levels = P.Levels;  
  int T = P.t;
  int t_offset  = 0;
  T == 1 ? t_offset = 0 : t_offset = 2;
  
  for(int t=0; t<T; t++) {

    int offset = t*( endNode(Levels,P) + 1 );
    
    cout << endl << "lev = " << 0 << "  T = " << t << endl;
    
    if(P.Vcentre) {
      cout << endl<< " Node number = " << 0 + offset << " : ";
      for(int i = 0; i < q+t_offset; i++) cout << NodeList[offset + 0].nn[i] << "  ";
    }      
    else {
      for(long unsigned int n = 0; n < endNode(0,P)+1; n++) {
	cout << endl<< " Node number = " << n + offset << " FL="<<NodeList[n].fwdLinks<<" : ";
	for(int i = 0; i < q+t_offset; i++) cout << NodeList[offset + n].nn[i] << "  ";
      } 
    }
    for(int lev = 1; lev < Levels+1; lev++)  {
      cout << endl << "lev = " << lev << "  T = " << t << endl;
      for(long unsigned int n = endNode(lev-1,P)+1; n < endNode(lev,P)+1; n++) {
	cout << endl<< " Node number = " << n + offset << " FL="<<NodeList[n].fwdLinks<<" : ";
	for(int i = 0; i < q+t_offset; i++) cout << NodeList[offset + n].nn[i] << "  ";
      }
    }      
  }  
  cout<<endl;
}

void PrintComplexPositions(const vector<Vertex> NodeList, Param P) {

  int Levels = P.Levels;
  
  if(P.verbosity) cout<<endl<<"#Printing for Level 0"<<endl;
  for(long unsigned int n=0; n<endNode(0,P)+1; n++) {
    if(P.verbosity) {
      cout<<"n="<<n<<" z="<<NodeList[n].z.real()<<","<<NodeList[n].z.imag();
      cout<<" |z|="<<abs(NodeList[n].z)<<" phi="<<arg(NodeList[n].z);
    }
  }
  for(int l=1; l<Levels+1; l++) {
    if(P.verbosity) cout<<endl<<"Printing for Level "<<l<<endl;
    for(long unsigned int n=endNode(l-1,P)+1; n<endNode(l,P)+1; n++) {
      if(P.verbosity) {
	cout<<"n="<<n<<" z="<<NodeList[n].z.real()<<","<<NodeList[n].z.imag();
	cout<<" |z|="<<abs(NodeList[n].z)<<" phi="<<arg(NodeList[n].z);
	cout<<endl;
      }
    }
  }
}

void CheckArea(const vector<Vertex> NodeList, Param P) {

  Float length_01 = 0.0;
  Float length_02 = 0.0;
  Float length_12 = 0.0;
  Float equi_area = area3q(P.q);
  Float ave       = 0.0;

  Float sig1 = 0.0;
  Float sig2 = 0.0;
  int count = 0;

  if(P.verbosity) cout<<endl<<"Checking boundary areas"<<endl;
  for(long unsigned int n=endNode(P.Levels-2,P) + 1; n<endNode(P.Levels-1,P)+1; n++) {
    for(int k=0; k<NodeList[n].fwdLinks; k++) {
      length_01 = d12(NodeList[n].z, NodeList[NodeList[n].nn[k+1]].z);
      length_02 = d12(NodeList[n].z, NodeList[NodeList[n].nn[k+2]].z);
      length_12 = d12(NodeList[NodeList[n].nn[k+1]].z, NodeList[NodeList[n].nn[k+2]].z);      
      ave += areaGeneral(P, length_01, length_02, length_12);
      count++;
    }
  }
  
  ave /= count;
  
  for(long unsigned int n=endNode(P.Levels-2,P) + 1; n<endNode(P.Levels-1,P)+1; n++) {
    for(int k=0; k<NodeList[n].fwdLinks; k++) {
      length_01 = d12(NodeList[n].z, NodeList[NodeList[n].nn[k+1]].z);
      length_02 = d12(NodeList[n].z, NodeList[NodeList[n].nn[k+2]].z);
      length_12 = d12(NodeList[NodeList[n].nn[k+1]].z, NodeList[NodeList[n].nn[k+2]].z);      
      if(P.verbosity) cout<<"n="<<n<<" area "<<k+1<<" = "<<areaGeneral(P, length_01, length_02, length_12)<<endl;
      sig1 += pow(equi_area - areaGeneral(P, length_01, length_02, length_12),2);
      sig2 += pow(ave       - areaGeneral(P, length_01, length_02, length_12),2);
    }
  }

  sig1 /= count - 1;
  sig2 /= count - 1;
  sig1 = sqrt(sig1);
  sig2 = sqrt(sig2);

  cout<<"Boundary areas"<<endl;
  cout<<"AREA EQUI = "<<equi_area<<endl;
  cout<<"AREA STD DEV W.R.T. EQUI = "<<sig1<<endl;
  cout<<"AREA AVE = "<<ave<<endl;  
  cout<<"AREA STD DEV W.R.T AVE = "<<sig2<<endl;  

}

void CheckEdgeLength(const vector<Vertex> NodeList, Param P) {
  
  int q = P.q;
  int Levels = P.Levels;
  Float length = 0.0;
  Float sig = 0.0;
  int  nn_node;
  bool Vcentre = P.Vcentre;
  Float length_0 = d12(NodeList[0].z, NodeList[1].z);
  Float tol = 1e-2;

  //Level 0 is specific to how the graph is centred.
  if(Vcentre) {
    if(P.verbosity) cout<<" lev =  " << 0 << endl;
    if(P.verbosity) cout<<endl<<" Node number = "<<0<<" : "<<endl;
    for(int i = 0; i < q; i++){
      nn_node = NodeList[0].nn[i];
      length = d12(NodeList[0].z, NodeList[nn_node].z);
      if(P.verbosity) cout<<" "<<NodeList[0].nn[i]<<" > "<<length<<"  ";
    }
  }
  else {
    if(P.verbosity) cout<<" lev = "<<0<<endl;
    if(P.verbosity) cout<<endl<<" Node number = "<<0<<" : "<<endl;
    for(int i = 0; i < 2; i++){
      nn_node = NodeList[0].nn[i];
      length = d12(NodeList[0].z, NodeList[nn_node].z);
      if(P.verbosity) cout << NodeList[0].nn[i] << " >  " << length<< "  ";
    }
  }
  
  for(int lev = 1; lev < Levels+1; lev++)  {
    if(P.verbosity) cout<<endl<<endl<<" lev = "<<lev<<endl;      
    for(long unsigned int n = endNode(lev-1,P) + 1;n < endNode(lev,P) + 1 ;n++) {
      if(P.verbosity) cout<<endl<<" Node number = "<<n<<":"<<endl;
      sig += pow( length_0 - d12(NodeList[n].z, NodeList[NodeList[n].nn[q-1]].z), 2);
      
      for(int i = 0; i <q; i++){
	nn_node = NodeList[n].nn[i];
	if(NodeList[n].nn[i] != -1 ) {
	  length = d12(NodeList[n].z, NodeList[nn_node].z);
	  if(P.verbosity) {
	    cout<<" to "<<NodeList[n].nn[i]<<" = "<<length<<" ";
	    if(abs(length - length_0)/length_0 > tol) cout<<"<-! "<<endl;
	    else cout<<"    "<<endl;
	  }
	}
      }
    }
  }
  sig /= endNode(Levels,P);
  sig = sqrt(sig);
  cout<<endl<<"LENGTH STD DEV = "<<sig<<endl;
  if(sig>tol) {
    cout<<"WARNING: Hypergeometric length STD_DEV has diverged over "<<tol<<endl;
    //exit(0);
  }
}


//Data file for lattice/analytical propagator data,
void Bulk2Bdry(vector<Vertex> NodeList, Float *phi, Param p) 
{
  long unsigned int TotNumber = (endNode(p.Levels,p) + 1) * p.t;
  Float norm = 0.0;
  for(long unsigned int i = 0;i < TotNumber; i++) norm += phi[i]*phi[i];
  for(long unsigned int i = 0;i < TotNumber; i++) phi[i] /= sqrt(norm); 
  long unsigned int j = p.src_pos;

  int T = p.t;
  int T_offset = 0;
  Float theta = 0.0;
  Float delta = p.t < 2 ? 0.5 + sqrt(0.25 + p.msqr) : 1.0 + sqrt(1 + p.msqr);
  complex<Float> ratio;
  complex<Float> src = NodeList[j].z;

  //Loop over timeslices
  for(int t=0; t<T; t++) 
    {
      T_offset = (endNode(p.Levels,p) + 1) * t;

      //Loop over circumference levels
      for(int lev=0; lev<p.Levels; lev++) 
	{

	  //Loop over H2 disk
	  for(long unsigned int k = endNode(lev,p)+1; k < endNode(lev+1,p)+1; k++) 
	    {
	      
	      //Construct i
	      int i = k + T_offset;
	      complex<Float> snk = NodeList[i].z;
	      ratio = NodeList[i].z/NodeList[j].z;
	      theta = atan2( ratio.imag() , ratio.real() );
	
	      //index divided by disk size, using the int floor feature/bug,
	      //gives the timeslice for each index.
	      int t1 = j / (TotNumber/p.t);
	      int t2 = i / (TotNumber/p.t);
	      //Assume PBC.
	      int delta_t = (t2-t1) > p.t/2 ? (t2-t1) - p.t : (t2-t1);
	      
	      Float geo_dist = d12(src,snk);
	      Float r = abs(NodeList[i].z);
	      Float r_p = abs(NodeList[j].z);
	      Float xi = (cosh(delta_t)*(1+r)*(1+r_p) - 4*r*r_p*cos(theta)) / ((1-r)*(1-r_p));
	      Float analytic_prop = log( exp(-delta*sigma(src,snk,delta_t)) /
					 (1 - exp(-2*sigma(src,snk,delta_t))));	   
	    }
	}
    }
}


void FourPoint(vector<Vertex> NodeList, Float *phi, Param p, int *src_array)    
{
  //Four-point contact term
  //Strategy: loop over configurations. For each configuration, loop over
  //which of the four sources the CG inverter is called on.
  int TotNumber = (endNode(p.Levels, p) + 1) * p.t;

  
  //vector<Vertex> NodeList(TotNumber);
  vector<Vertex> NodeList1(TotNumber);
  vector<Vertex> NodeList2(TotNumber);
  vector<Vertex> NodeList3(TotNumber);
  vector<Vertex> NodeList4(TotNumber);

  //Initialise. -1 in NodeList indicates that node is not yet populated.
  for(int n = 0; n <TotNumber;n++)
    for(int mu = 0; mu < p.q+2; mu++) 
      { 
	NodeList[n].nn[mu] = -1;
	NodeList1[n].nn[mu] = -1;
	NodeList2[n].nn[mu] = -1;
	NodeList3[n].nn[mu] = -1;
	NodeList4[n].nn[mu] = -1;
      }

  int src_pos1 = src_array[0];
  int src_pos2 = src_array[1];
  int src_pos3 = src_array[2];
  int src_pos4 = src_array[3];

  int src_configs = src_array[3] - src_array[1];

  int temp = src_array[2];
  int count = 0;
  //Float Amp4 = 0.0;
  Float prop = 1.0;
  Float* Amp_ana = new Float[src_configs];
  Float* Amp_array = new Float[src_configs]; 
  Float* Amp_pos_array = new Float[src_configs];
  Float* u_array = new Float[src_configs];
  Float* src_pos3_array = new Float[src_configs];
  Float* fourpt_info = new Float[src_configs];
  Float* fourpt_num = new Float[src_configs];
  for(int src_pos3 = temp; src_pos3 < src_array[3]; src_pos3++)
    {
      //cout<<"src_pos3 is: "<<src_pos3<<" out of "<<src_array[3]<<"\n";
      Float xr12 = NodeList[src_pos1].z.real()-NodeList[src_pos2].z.real();
      Float xr34 = NodeList[src_pos3].z.real()-NodeList[src_pos4].z.real();
      Float xr13 = NodeList[src_pos1].z.real()-NodeList[src_pos3].z.real();
      Float xr24 = NodeList[src_pos2].z.real()-NodeList[src_pos4].z.real();
      Float xi12 = NodeList[src_pos1].z.imag()-NodeList[src_pos2].z.imag();
      Float xi34 = NodeList[src_pos3].z.imag()-NodeList[src_pos4].z.imag();
      Float xi13 = NodeList[src_pos1].z.imag()-NodeList[src_pos3].z.imag();
      Float xi24 = NodeList[src_pos2].z.imag()-NodeList[src_pos4].z.imag();
      Float u_cr = ((xr12*xr12+xi12*xi12)*(xr34*xr34+xi34*xi34))/((xr13*xr13+xi13*xi13)*(xr24*xr24+xi24*xi24));
      Float Dfunction =    (sqrt(M_PI)/(2*(xr13*xr13+xi13*xi13)*(xr24*xr24+xi24*xi24)))*((-2*log(1- sqrt(u_cr)))/(sqrt(u_cr)) + log(u_cr)/(sqrt(u_cr)-1));

      Float* prop = new Float[TotNumber-4];
      Float* prop1 = new Float[TotNumber-4];    
      Float* prop2 = new Float[TotNumber-4];    
      Float* prop3 = new Float[TotNumber-4];    
      Float* phi = new Float[TotNumber];
      Float* phi1 = new Float[TotNumber];
      Float* phi2 = new Float[TotNumber];
      Float* phi3 = new Float[TotNumber];
      Float* b   = new Float[TotNumber];
      Float* b1   = new Float[TotNumber]; 
      Float* b2   = new Float[TotNumber]; 
      Float* b3   = new Float[TotNumber]; 

      for(int i=0; i<TotNumber; i++) 
	{
	  
	  phi[i] = 0.0;
	  phi1[i] = 0.0;
	  phi2[i] = 0.0;
	  phi3[i] = 0.0;
	  b[i] = 0.0;
	  b1[i] = 0.0;
	  b2[i] = 0.0;
	  b3[i] = 0.0;
	  
	  /*
	  phi[i] = NodeList1[i].phi;
	  phi1[i] = NodeList2[i].phi;
	  phi2[i] = NodeList3[i].phi;
	  phi3[i] = NodeList4[i].phi;
	  */
	  if(i<TotNumber-4) {
	    prop[i] = 0.0;
	    prop1[i] = 0.0;
	    prop2[i] = 0.0;
	    prop3[i] = 0.0;
	  }
	}
      
      for(int src_num=0; src_num<4; src_num++)
	{
	  if(src_num == 0)
	    {
	      b[src_pos1] = 100000000;
	      //cout<<"src_num is: "<<src_num<<" and src_array[src_num] is: "<<src_array[src_num]<<"\n";
	      Float truesq = 0.0;
	      truesq = Minv_phi(phi, b, NodeList1, p);   //Residual to this CG solution
	      //cout<<"Tolerance = "<<p.tol<<" True Residual = "<<sqrt(truesq)<<endl;
	      
	      //FourPoint(NodeList, phi, p, src_num, src_array, src_pos3);

	      //cout<<phi[0]<<phi[1]<<"\n";
	      //cout<<(Float)p.N_latt*(Float)phi[0]<<"\n";
	      for(int i=0;i<TotNumber-4;i++)  
		{
		  if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))
		    {
		      prop[i] = (Float)p.N_latt*(Float)phi[i];
		    }     
		}
	      //cout<<prop[0]<<"\n";
	    }
	  
	  if(src_num == 1)
	    {
	      //cout<<"src_num is: "<<src_num<<" and src_array[src_num] is: "<<src_array[src_num]<<"\n";   
	      b1[src_pos2] = 100000000;
	      Float truesq = 0.0;
	      truesq = Minv_phi(phi1, b1, NodeList2, p);
	      //FourPoint(NodeList, phi1, p, src_num, src_array, src_pos3);
	      for(int i=0;i<TotNumber-4;i++) 
		{
		  if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))
		    {  
		      prop1[i] = (Float)p.N_latt*(Float)phi1[i];
		    }
		}
	    }
	  
	  if(src_num == 2)
	    {
	      //cout<<"src_num is: "<<src_num<<" and src_array[src_num] is: "<<src_array[src_num]<<"\n";
	      //cout<<"src_pos2 is: "<<src_pos2<<"\n";
	      b2[src_pos3] = 100000000;
	      Float truesq = 0.0;
	      truesq = Minv_phi(phi2, b2, NodeList3, p);
	      //FourPoint(NodeList, phi2, p, src_num, src_array, src_pos3);
	      for(int i=0;i<TotNumber-4;i++)  
		{
		  if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))  
		    {
		      prop2[i] = (Float)p.N_latt*(Float)phi2[i];
		    }
		}
	    }
	  if(src_num == 3)
	    {
	      //cout<<"src_num is: "<<src_num<<" and src_array[src_num] is: "<<src_array[src_num]<<"\n";   
	      b3[src_pos4] = 100000000;
	      Float truesq = 0.0;
	      truesq = Minv_phi(phi3, b3, NodeList4, p);
	      //FourPoint(NodeList, phi3, p, src_num, src_array, src_pos3);
	      for(int i=0;i<TotNumber-4;i++) 
		{
		  if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))  
		    {
		      prop3[i] = (Float)p.N_latt*(Float)phi3[i];
		    }
		}
	    }
	  //cout<<"AT END "<<"\n";
	  //cout<<prop[0]<<", "<<prop1[0]<<", "<<prop2[0]<<", "<<prop3[0]<<"\n";
	}
      Float Amp4 = 0.0;
      for(int i=0; i<TotNumber; i++)
	{
	  //cout<<"prop["<<i<<"] is: "<<prop[i]<<"\n";
	  //cout<<"prop1["<<i<<"] is: "<<prop1[i]<<"\n";    
	  //cout<<"prop2["<<i<<"] is: "<<prop2[i]<<"\n";    
	  //cout<<"prop3["<<i<<"] is: "<<prop3[i]<<"\n";    
	  Amp4 += prop[i]*prop1[i]*prop2[i]*prop3[i]*area3q(p.q);
	  //Amp4 += (Float)p.N_latt*(Float)phi[i]*(Float)p.N_latt*(Float)phi1[i]*(Float)p.N_latt*(Float)phi2[i]*(Float)p.N_latt*(Float)phi3[i];
	} 
      //cout<<"Amp4 is: "<<Amp4<<"\n";
      Amp_array[count] = Amp4;
      Amp_pos_array[count] = abs(NodeList[src_pos3].z);
      //cout<<"count is: "<<count<<"\n";
      //cout<<"u_cr is: "<<u_cr<<"\n";
      u_array[count] = u_cr;
      Amp_ana[count] = Dfunction;
      src_pos3_array[count] = src_pos3;
  
  

      /*
      for(int i=0; i<TotNumber-4; i++)
	{
	  cout<<"in prop loop"<<"\n";
	  cout<<(Float)p.N_latt*(Float)phi[i]<<"\n";
	  prop += (Float)p.N_latt*(Float)phi[i]*(Float)p.N_latt*(Float)phi1[i]*(Float)p.N_latt*(Float)phi2[i]*(Float)p.N_latt*(Float)phi3[i];  
	}
      */
  
  
      //fourpt_info[count] = Amp4;
      fourpt_num[count] = src_pos3;
      count++;

      delete[] prop;
      delete[] prop1;
      delete[] prop2;
      delete[] prop3;
      delete[] phi;
      delete[] phi1; 
      delete[] phi2; 
      delete[] phi3;
      delete[] b;
      delete[] b1; 
      delete[] b2; 
      delete[] b3; 
    }
  
  //End of 4-pt contact term

  cout<<"LOOK HERE: "<<d12(NodeList[0].z, NodeList[1].z)<<"\n";

  /*
  int n_shift = p.n_shift;
  cout<<"n_shift is: "<<n_shift<<"\n";
  Float** phi_ms = new Float*[n_shift];
  for(int i=0; i<n_shift; i++) {
    phi_ms[i] = new long double[TotNumber];
    for(int j=0; j<TotNumber; j++) phi_ms[i][j] = 0.0;
  }
  */
  cout<<"src_configs is: "<<src_configs<<"\n";

    
  ofstream myfile;
  myfile.open ("Amp4_HB_"+std::to_string(iteration)+".dat");
  for(int i=0; i<src_configs; i++) 
    {
      myfile<<u_array[i]<<" "<<Amp_array[i]<<" "<<Amp_ana[i]<<" "<<src_pos3_array[i]<<" "<<Amp_pos_array[i]<<"\n";
    }
  myfile.close();
  
}




// /********************************************
// Basic Hyperbolic Algebra. 

// Moebius Transforms 

// Hyperbolic reflection for geodesic in UHP

// Given Line: (x-a)^2 + y^2 = r^2 
// Determined by x = a \pm r at y = 0
//               x = a, y = r


//  RL(a,r, u) = a + r^2/(conj(u) - a)

// Preserves the line and swap a  to infty.

// Map to Disc: z = D(u) =  (u -I)/(1 -I * u) with u = x + i y
// Map to UHP   u = U(z) = (z + I)/(1 + I * z);

// Find UHP circle that hits  +/- theta_0 on  Disc

// |z - A|^2 = R^2  
// |z|^2 - |z| |A| cos(theta) + |A|^2 = R^2
// boundary 1 - |A| cos (theta) + |A|^2 = R^2
// pick A = real. when a = 0 with map

// Need 3 point to define the mobius. Circle to Circle. 

// ****************************************************************/

// // Translate w to 0 
// complex<Float> T(complex<Float> z,  complex<Float> w)
// { //translate w to 0
//   return (z - w)/(z*conj(w) + (Float)1.0);
// }

// // Rotate at z = 0
// complex<Float> R(complex<Float> z, complex<Float> omega)
// {
//   // rotate by omega = exp [i theta] about z= 0
//   return omega*z;
//  }

// //Reflection z accross the z1 -- z2 line
// complex<Float> flip(complex<Float> z, complex<Float> z1, complex<Float> z2)
// {
//   // reflection (or flip)  z across (z1,z2)
//   return T( R( conj( T(z,z1) ), z2/conj(z2) ), -z1 );
// }

// //Geodesic from z = 0 to z
// Float  s(complex<Float> z)
// {
//    return log(((Float)1.0+abs(z))/((Float)1.0-abs(z)));
// }

// //Geodesic distance s from origin
// Float r(Float s)
// {
//   return tanh(s/2);
// }

// //Geodesic distance from z1 to z2
// Float d12(complex<Float> z, complex<Float> w)
// {
//   return log ( (abs((Float)1.0-conj(z)*w) + abs(z-w))/(abs((Float)1.0-conj(z)*w) - abs(z-w)));
// }

// //Geodesic distance from z1,t1 to z2,t2
// Float sigma(complex<Float> z, complex<Float> w, int delta_t) {

//   Float theta = atan2( (w/z).imag() , (w/z).real() );
//   Float r = abs(z);
//   Float r_p = abs(w);  
//   Float xi = (cosh(delta_t)*(1+r)*(1+r_p) - 4*r*r_p*cos(theta)) / ((1-r)*(1-r_p)); 
  
//   return acosh(xi);
    
// }

// // length of arc q fold triangle to origin.
// Float s3p(int q)
// {  //vertex centered Arc lengeth
//   return (Float)2.0*acosh((Float)1.0/sin(M_PI/(Float)q));
// }

// // Area equilateral triangle with angles 2 pi/q
// Float area3q(int q)
// {
//   //pi - (3 * hyp_angle) = defect
//   return M_PI - (Float)3.0*(2.0*M_PI/(Float)q);
// }

// // Area non-equilateral triangle with side length a,b,c
// Float areaGeneral(Param P, Float a, Float b, Float c) {
//   //pi - (A+B+C) = defect
  
//   // use general cosine law:
//   // cos(C) = (cosh(c) - cosh(a)cosh(b)) / sinh(a)sinh(b)
//   Float C = acos( -(cosh(c) - cosh(a)*cosh(b)) / (sinh(a)*sinh(b)) );
  
//   // use general sine law:
//   // sin(A)/sinh(a) = sin(B)/sinh(b) = ...
//   Float B = asin( sinh(b)*sin(C)/sinh(c) );
//   Float A = asin( sinh(a)*sin(C)/sinh(c) );

//   return M_PI - (A+B+C);
// }

// //
// Float centralRad(Float s)
// {
//   return (sqrt( cosh(s/2.0) - (Float)1.0/4.0) - 0.5*sqrt(3.0))/sqrt(cosh(s/2.0) -(Float)1.0);
// }

// //
// complex<Float> DisktoUHP(complex<Float> z)
// {
//   // z = -1 -i, 1, i maps to u =  -1, 0 1, infty
//   return (z + I)/((Float)1.0 + I * z);
// }
// //
// complex<Float> UHPtoDisk(complex<Float> u)
// {
//   // u = 0, 1, infty  maps to -1 -i , 1, i  
//   return (u - I)/((Float)1.0 - I*u); 
// }

// //- Rotate z about z0 by 2*k*pi/q 
// complex<Float> newVertex(complex<Float> z,complex<Float> z0, int k, int q) 
// {
//   complex<Float> w( 0.0, 2.0 * sin(k * M_PI/q) );
//   complex<Float> a( cos(k*M_PI/q)*((Float)1.0 - norm(z0)), sin(k*M_PI/q)*((Float)1.0 + norm(z0)) ); 
//   w = w*z0;
  
//   //cout<<"New z = "<<-(a*z - w)/(conj(w)*z - conj(a))<<endl;
//   return - (a*z - w)/(conj(w)*z - conj(a)); 
// }


// complex<Float> inversion(complex<Float> z0, Float r)
// {
//   // z_image conj(z0) = r^2
//   return r*2/conj(z0);
// }

// complex<Float> squareInversion(complex<Float>z0,Float r1,Float r2 )
// {
//   return inversion(inversion(z0, r1),r2);
// }

// Float greens2D(complex<Float> z, complex<Float> w)
// {
//   return -log( tanh ( log ( (abs((Float)1.0-conj(z)*w) + abs(z-w))/(abs((Float)1.0-conj(z)*w) - abs(z-w)) )/2 ) );    
// }

// Float greensM2D(complex<Float> z, complex<Float> w, Param p)
// {

//   //First compute 2F1  

//   Float delta = 0.5 + sqrt(0.25 + p.msqr);
//   Float h = 1;
//   Float result = 0.0;
//   Float result_0 = 0.0;
//   Float geo = exp(-2*d12(z,w));
//   Float a,b,c;
//   Float tol = 1e-10;
//   int n=0;
//   bool conv = false;

//   while( !conv && n < 10000 ) {    
//     result_0 = result;
//     a = tgamma(delta + n)/tgamma(delta);
//     b = tgamma(h + n)/tgamma(h);
//     c = tgamma(delta+1-h + n)/tgamma(delta+1-h);
//     result += ( (a*b) / (c*tgamma(n+1)) ) * pow(geo,n);
//     if( abs(result_0 - result)/result_0 < tol ) conv = true;
//     n++;
//     if(n%10000 == 0) cout<<n<<" 2F1 iters "<<geo<<" "<<abs(result_0 - result)/result_0<<endl; 
//   }
  
//   //Then compute coefficient. 
//   result *= pow(geo,delta/2) * tgamma(delta) / (2*pow(M_PI,h)*tgamma(delta+1-h));

//   return result;
// }


// /*
// //D-function
// Float DFunc(complex<Float> z, complex<Float> w)
// {
//   return ((z * conj(w))/(z-w)) * (2*PolyLog[2,z] - 2*PolyLog[2,w]+Log[z*w]*Log[(1-z)/(1-w)])
// }
// */


// /* (3,p)  packing lenght side = 2 ch^{-1} (1/2\sin (\pi/p)) 

// cosh(side/2) = 1/(2 sin(\pi/p)) 

// r_0=[ sqrt{ 4 cosh^2(side/2) -1} -(sqrt 3)/2]/ sqrt{cosh^2(side/2) -1}

// */


// //Find the minimum value of an array of floats
// float GetMinBoundaryRadius(vector<Vertex> NodeList, Param p)
// {
//   int Levels = p.Levels;
//   int size = endNode(Levels,p)-endNode(Levels-1,p);   //number of nodes in outer level
//   float outerCoords[endNode(Levels,p)-endNode(Levels-1,p)];
//   for(long unsigned int i=endNode(Levels-1,p)+1; i<endNode(Levels,p)+1; i++)
//     {
//       outerCoords[i-endNode(Levels-1,p)-1] = abs(NodeList[i].z);
//     }
//   float minValue = 1;
//   for (int i = 0; i < size; i++) 
//     {
//       //Find minValue
//       if(outerCoords[i] != 0)
// 	{
// 	  if (outerCoords[i] < minValue)
// 	    {
// 	      minValue = outerCoords[i];
// 	    }
// 	}
//     }
//   return minValue;
// }

// //Find the maximum value of an array of floats
// float GetMaxBoundaryRadius(vector<Vertex> NodeList, Param p)
// {
//   int Levels = p.Levels;
//   int size = endNode(Levels,p)-endNode(Levels-1,p);   //number of nodes in outer level
//   float outerCoords[endNode(Levels,p)-endNode(Levels-1,p)];
//   for(long unsigned int i=endNode(Levels-1,p)+1; i<endNode(Levels,p)+1; i++)
//     {
//       outerCoords[i-endNode(Levels-1,p)-1] = abs(NodeList[i].z);
//     }
//   float maxValue = outerCoords[0];
//   for (int i = 0; i < size; i++) 
//     {
//       //Find maxValue 
//       if (outerCoords[i] > maxValue) 
// 	{
// 	  maxValue = outerCoords[i];
// 	}
//     }
//   return maxValue;
// }


//Data file for lattice/analytical propagator data,
void FourPoint(vector<Vertex> NodeList, Float *phi, Param p, int src_number, int *src_array, int src_position3) 
{
  int src_pos3 = src_position3;
  int src_num = src_number;
  long unsigned int TotNumber = (endNode(p.Levels,p) + 1) * p.t;
  int outer_num = endNode(p.Levels,p)-endNode(p.Levels-1,p);

  int src_pos4 = src_array[3];
  int src_pos2 = src_array[1];
  int src_pos1 = src_array[0];
  //int src_pos2 = src_array[1];

  //Normalize the phi
  Float norm = 0.0;
  for(long unsigned int i = 0;i < TotNumber; i++) norm += phi[i]*phi[i];
  for(long unsigned int i = 0;i < TotNumber; i++) phi[i] /= sqrt(norm); 
  
  //long unsigned int j = p.src_pos;
  long unsigned int j = src_pos3;

  int T = p.t;
  int T_offset = 0;
  Float theta = 0.0;
  Float delta = p.t < 2 ? 0.5 + sqrt(0.25 + p.msqr) : 1.0 + sqrt(1 + p.msqr);
  complex<Float> ratio;
  complex<Float> src = NodeList[j].z;
 
  //Loop over timeslices
  //for(int t=0; t<T; t++) 
  int t=0;
  // {
      T_offset = (endNode(p.Levels,p) + 1) * t;

      //Loop over circumference levels
      //for(int lev=0; lev<p.Levels; lev++)
      //for(int src_pos2 = src_pos3 + 1; src_pos2 < src_pos1; src_pos2++)
      //	{
	  sprintf(p.fname, "q%d_Lev%d_T%d_msqr%.3Le_srct0_srcpos%d_srcnum%d_sinkt%dLev%d_%s_%s.dat",
		  p.q,
		  p.Levels,
		  p.t,
		  (Float)p.msqr,
		  src_pos3,
		  src_num,
		  t,
		  p.Levels+1,
		  p.bc == true ? "Dirichlet" : "Neumann",
		  p.Vcentre == true ? "Vertex" : "Circum");
	  FILE *fp1;
	  fp1=fopen(p.fname, "w");
      
	  Float prop = 1.0;

	  //Loop over H2 disk
	  for(long unsigned int k = endNode(0,p)+1; k < endNode(p.Levels,p)+1; k++) 
	    {
	      //Construct i - sink index
	      int i = k + T_offset;
	      
	      ratio = NodeList[i].z/NodeList[j].z;
	      theta = atan2( ratio.imag() , ratio.real() );
	      complex<Float> snk = NodeList[i].z;   //location of current sink (endpoint)
	      
	      //index divided by disk size, using the int floor feature/bug,
	      //gives the timeslice for each index.
	      int t1 = j / (TotNumber/p.t);
	      int t2 = i / (TotNumber/p.t);
	      //Assume PBC.
	      int delta_t = (t2-t1) > p.t/2 ? (t2-t1) - p.t : (t2-t1);
	      
	      NodeList[j].z;

	      Float r = abs(NodeList[i].z);   //Radius of source
	      Float r_p = abs(NodeList[j].z); //Radius of sink
	      Float xi = (cosh(delta_t)*(1+r)*(1+r_p) - 4*r*r_p*cos(theta)) / ((1-r)*(1-r_p));
	  
	      Float analytic_prop = log( exp(-delta*sigma(src,snk,delta_t)) /
				      (1 - exp(-2*sigma(src,snk,delta_t))));
    
	      //cout<<"src_pos1 is: "<<src_pos1<<" src_pos 2 is: "<<src_pos2<<" src_pos3 is: "<<src_pos3<<" src_pos4 is: "<<src_pos4<<"\n";

	      Float xr12 = NodeList[src_pos1].z.real()-NodeList[src_pos2].z.real();
	      Float xr34 = NodeList[src_pos3].z.real()-NodeList[src_pos4].z.real();
              Float xr13 = NodeList[src_pos1].z.real()-NodeList[src_pos3].z.real();
              Float xr24 = NodeList[src_pos2].z.real()-NodeList[src_pos4].z.real();
	      Float xi12 = NodeList[src_pos1].z.imag()-NodeList[src_pos2].z.imag();
	      Float xi34 = NodeList[src_pos3].z.imag()-NodeList[src_pos4].z.imag();
	      Float xi13 = NodeList[src_pos1].z.imag()-NodeList[src_pos3].z.imag();
	      Float xi24 = NodeList[src_pos2].z.imag()-NodeList[src_pos4].z.imag();

	      Float u_cr = ((xr12*xr12+xi12*xi12)*(xr34*xr34+xi34*xi34))/((xr13*xr13+xi13*xi13)*(xr24*xr24+xi24*xi24));
	      Float Dfunction = (-2*log(1- sqrt(u_cr)))/(sqrt(u_cr)) + log(u_cr)/(sqrt(u_cr)-1);
	      
	      
	      
	      
	      //cout<<"Before IF. src_num is: "<<src_num<<". src_array[src_num] is: "<<src_array[src_num]<<"\n";
		if( (i != src_pos1) && (i != src_pos3) && (i != src_pos4) && (i != src_pos2))  
		{
		  fprintf(fp1, "%d %d %.8Le %.8Le %.8Le %.8Le %.8Le %.8Le %.8Le %.8Le %.8Le %.8Le %.8Le %d %.8Le \n",
			  //1 Timeslice, 2 H2 pos
			  t, i,
			  
			  //3 source/sink angle
			  (Float)theta,
			  
			  //4 lattice prop
			  (Float)p.N_latt*(Float)phi[i],
			  
			  //5 invariant
			  1.0/xi,
			  //(Float)( ((Float)1.0-abs(snk))*((Float)1.0-abs(src))/(pow(abs(snk - src),2))),
			  
			  //6 AdS2p1 formula
			  (Float)(exp(-delta*sigma(src,snk,delta_t)) / (1 - exp(-2*sigma(src,snk,delta_t)))),
			  
			  //7 geodesic
			  (Float)sigma(src,snk,delta_t),

			  //8 
			  NodeList[i].z.real(),

			  //9 
			  NodeList[i].z.imag(),

			  //10
			  NodeList[j].z.real(),

			  //11
			  NodeList[j].z.imag(),

			  //12 u cross-ratio 
			  u_cr,

			  //13 D-function 
			  Dfunction,
			  
			  //14 j
			  src_array[src_num],

			  //15 Analytic Propagator
			  analytic_prop
			  );
		}
	    }
	  //prop = prop*(Float)p.N_latt*(Float)phi[i];
	  fclose(fp1);
	  //	}
}

/*
void latticeScaling(vector<Vertex> &NodeList, Param& p)
{
  int latVol  = p.latVol;
  int lower = endNode(p.Levels-1,p)+1;
  long unsigned int TotNumber = (endNode(p.Levels,p) + 1) * p.t;
  int outer_cirum = endNode(p.Levels,p)-endNode(p.Levels-1,p);

  int t1, t2, delta_t;
  double delta = 1.0 + sqrt(1 + p.msqr);
  double theta, r, r_p;
  complex<double> ratio;
  complex<double> snk;
  double* analytic_prop = new double[outer_cirum*p.Lt/2];
  double* xi_invariant  = new double[outer_cirum*p.Lt/2];

 
  // Construct the xi invariant and the analytic propagator
  
  int j = lower;
  complex<Float> src = NodeList[j].z;  
  //Loop over timeslices
  for(int t=0; t<p.Lt/2; t++) {
    int T_offset = (endNode(p.Levels,p) + 1) * t;

    //Loop over outer circle of the Poincare disk
    for(long unsigned int k = 0; k < outer_cirum; k++) {

      //Construct i (full lattice index)
      int i = T_offset + lower + k;
      //Construct s (surface lattice index)
      int s = t*outer_cirum + k;
      
      ratio = NodeList[i].z/NodeList[j].z;
      theta = atan2( ratio.imag() , ratio.real() );
      //index divided by disk size, using the int floor feature/bug,
      //gives the timeslice for each index.
      int t1 = j / (TotNumber/p.t);
      int t2 = i / (TotNumber/p.t);

      //Assume PBC.
      int delta_t = (t2-t1);
      complex<Float> snk = NodeList[i].z;
      Float r   = abs(NodeList[i].z);
      Float r_p = abs(NodeList[j].z);
      
      Float xi_invariant  = log( ((1-r)*(1-r_p))/(cosh(delta_t)*(1+r)*(1+r_p)
					       - 4*r*r_p*cos(theta)) );
      
      Float analytic_prop = log( exp(-delta*sigma(src,snk,delta_t)) /
			      (1 - exp(-2*sigma(src,snk,delta_t))));
      
      //if(s<10) cout<<"s="<<s<<" xi="<<xi_invariant[s]<<" ap="<<analytic_prop[s]<<endl;
    }
  }

  double* c = new double[2];
  double* cov_ssq = new double[4];
  gsl_fit_linear(xi_invariant, 1, analytic_prop, 1, outer_cirum*p.Lt/2, &c[0], &c[1],
		 &cov_ssq[0], &cov_ssq[1], &cov_ssq[2], &cov_ssq[3]);

  cout<<"Target data"<<endl;
  cout<<"GSL data: C="<<c[0]<<" M="<<c[1]<<endl;
  cout<<"          covar00 = "<<cov_ssq[0]<<endl;
  cout<<"          covar01 = "<<cov_ssq[1]<<endl;
  cout<<"          covar11 = "<<cov_ssq[2]<<endl;
  cout<<"          sum_sq  = "<<cov_ssq[3]<<endl;

  double grad = c[1];
  double d_grad = 100;

  double grad_tol = 1e-4;

  //Search in wisdom file for a shorcut
  bool preTuned = false;
  bool tuneWin = false;
  ifstream fileInput;
  string line;
  char params[256];
  sprintf(params, "%d %d %d %.4f", p.q, p.Levels, p.Lt, p.msqr);
  char* search = params;
  unsigned int curLine = 0;
  // open file to search
  fileInput.open("ads_wisdom");
  if(fileInput.is_open()) {
    while(getline(fileInput, line)) { 
      curLine++;
      if (line.find(search, 0) != string::npos) {
	cout<<"Found data for your problem! Lucky you..."<<endl;
	preTuned = true;
	getline(fileInput, line);
	p.N_latt = stof(line);
	getline(fileInput, line);
	p.C_msqr = stof(line);
	cout<<"Using N_latt="<<p.N_latt<<endl;
	cout<<"Using C_msqr="<<p.C_msqr<<endl;
	fileInput.close();
      }
    }
    if(!preTuned)
      cout<<endl<<"AdS wisdom data not found. Strap in for some tuning..."<<endl<<endl;
  }
  else cout<<endl<<"AdS wisdom file not found. Strap in for some tuning..."<<endl<<endl;

  //Begin search for correct scaling factors.
  int iter = 0;
  //while(abs(d_grad) > grad_tol || abs(d_inter) > inter_tol) {
  while(abs(d_grad) > grad_tol) {
    
    double* phi_ave = new double[latVol];  
    double* phi = new double[latVol];
    double* b   = new double[latVol];

    //Take the average data from sources in the the qth sector
    //of the outer level. 
    int sources = (endNode(p.Levels,p) - endNode(p.Levels-1,p)) / p.q ;

    //initialise, invert, average.
    for(int i=0; i<latVol; i++) phi_ave[i] = 0.0;    
    for(int s=0; s<sources; s++) {
      for(int i=0; i<latVol; i++) {
	b[i] = 0.0;
	phi[i] = 0.0;
      }
      b[lower + s] = 1.0;    
      Minv_phi(phi, b, NodeList, p);
      for(int i=0; i<latVol; i++) phi_ave[i] += phi[i]/sources;
    }
    
    //Use current lattice normalisation.
    for(int i=0; i<latVol; i++) {
      phi_ave[i] = log(p.N_latt*phi_ave[i]);
    }
    
    //phi_ave now contains an averaged solution vector. We now
    //perform linear regression on this vector and the analytic
    //prop (w.r.t the xi invariant) to compute the relevant
    //scaling and normalsation.
    
    double* latt_prop  = new double[outer_cirum*p.Lt/2];
    //Loop over timeslices
    for(int t=0; t<p.Lt/2; t++) {
      int T_offset = (endNode(p.Levels,p) + 1) * t;
      //Loop over H2 disk
      for(long unsigned int k = 0; k < outer_cirum; k++) {
	
	//Construct i (full lattice AdS2p1 index)
	int i = T_offset + lower + k;
	//Construct s (surface lattice 2D index)
	int s = t*outer_cirum + k;
	
	latt_prop[s] = phi_ave[i];
      }
    }
  
    //Extract linear fit data from log-log plot.
    gsl_fit_linear(xi_invariant, 1, latt_prop, 1, outer_cirum*p.Lt/2, &c[0], &c[1],
		   &cov_ssq[0], &cov_ssq[1], &cov_ssq[2], &cov_ssq[3]);
    
    cout<<"At iteration "<<iter<<endl;
    cout<<"GSL data: C="<<c[0]<<" M="<<c[1]<<endl;
    cout<<"          covar00 = "<<cov_ssq[0]<<endl;
    cout<<"          covar01 = "<<cov_ssq[1]<<endl;
    cout<<"          covar11 = "<<cov_ssq[2]<<endl;
    cout<<"          sum_sq  = "<<cov_ssq[3]<<endl;
    
    //Adjust the parameters and start over.
    d_grad = (c[1] - grad)/abs(grad);
    cout<<"D grad = "<<d_grad<<endl;
    cout<<"C_msqr = "<<p.C_msqr<<endl<<endl;
    
    double fac1 = abs(p.C_msqr);
    
    if( abs(d_grad) > grad_tol ) {
      if(c[1] > grad) p.msqr < 0 ? p.C_msqr += (fac1*abs(d_grad) + 1*grad_tol) : p.C_msqr -= (fac1*abs(d_grad) + 1*grad_tol);
      if(c[1] < grad) p.msqr < 0 ? p.C_msqr -= (fac1*abs(d_grad) + 1*grad_tol) : p.C_msqr += (fac1*abs(d_grad) + 1*grad_tol);
    }
    
    if(cov_ssq[0] != cov_ssq[0]) {
      cout<<"GSL failed!"<<endl;
      exit(0);
    }
    
    //Did it tune properly?
    if(d_grad < grad_tol) tuneWin = true;
    
    delete[] phi_ave;
    delete[] phi;
    delete[] b;
    delete[] latt_prop;

    iter++;
    
  }
  
  //If this is a new problem, and it tuned properly, save the data.
  if(!preTuned && tuneWin) {
    FILE *file = fopen("ads_wisdom", "a");
    fprintf(file,"%d %d %d %.4f\n%f\n%f\n",
	    p.q, p.Levels, p.Lt, p.msqr,
	    p.N_latt,
	    p.C_msqr);
    fclose(file);
  }

  delete[] c;
  delete[] cov_ssq;

  delete[] analytic_prop;
  delete[] xi_invariant;
}
*/
#endif
