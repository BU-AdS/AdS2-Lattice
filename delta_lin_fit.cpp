#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

//#include "delta_lin_fit.h"

using namespace std;

int main(int argc, char *argv[])
{
  
  /***** Initialize and read in files *****/

  string dir = "/Users/hattrick/Desktop/Dean/adsSpectrum/AdS2-Lattice/";

  if(argc != 4)
    {
      cout<<"Correct format is Mass, Level, Triangulation."<<endl;
      cout<<"Please try again."<<endl;
      return -1;
    }
  
  double M = atof(argv[1]);
  int L = atof(argv[2]);
  int Q = atof(argv[3]);
  
  //stringstream stream;
  //stream<<fixed<<setprecision(5)<
 
 

   stringstream mstream, lstream, qstream;
	     
   mstream<<fixed<<setprecision(6)<<M;
   lstream<<fixed<<setprecision(6)<<L;
   qstream<<fixed<<setprecision(6)<<Q;
	     
   string ham_loc=dir+"ham_m"+mstream.str()+"_L"+lstream.str()+"_q"+qstream.str()+".dat";
   string dist_loc=dir+"distance_m"+mstream.str()+"_L"+lstream.str()+"_q"+qstream.str()+".dat";
   string green_loc=dir+"green_m"+mstream.str()+"_L"+lstream.str()+"_q"+qstream.str()+".dat";
		      
   ifstream hamin(ham_loc);
   ifstream distin(dist_loc);
   ifstream greenin(green_loc);
   double num1;
   double num2;
   double num3;
   vector<double> hamvec;
   vector<double> distvec;
   vector<double> greenvec;

   while(hamin >> num1)
     {
       hamvec.push_back(num1);
     }
   while(distin >> num2)
     {
       distvec.push_back(num2);
     }
   while(greenin >> num3)
     {
       greenvec.push_back(num3);
     }

   int TotNumber = sqrt(hamvec.size());
   if(sqrt(hamvec.size()) != sqrt(distvec.size()))
     cout << "ERROR: hamiltonian and distance matrices not equal!" << endl;

   float** ham = new float*[TotNumber];
   float** dist = new float*[TotNumber];
   float** datalin = new float*[TotNumber];
   for(int i=0; i<TotNumber; i++)
     {
       ham[i] = new float[TotNumber];
       dist[i] = new float[TotNumber];
       datalin[i] = new float[TotNumber];
     }

   for(int i=0; i<TotNumber; i++)
     {
       for(int j=0; j<TotNumber; j++)
	 {
	   ham[i][j] = hamvec[TotNumber*i+j];
	   dist[i][j] = distvec[TotNumber*i+j];
	   datalin[i][j] = 0;
	 }
     }


   /*  
   cout << "count = " << vec.size() << endl;
   cout << "first = " << vec[0] << endl;
   cout << "last = " << vec[vec.size()-1] << endl;
   */
	     
  /***** Linearize *****/
        
   for(int i=0; i<TotNumber; i++)
     {
       for(int j=0; j<TotNumber; j++)
	 {
	   datalin[i][j] = -dist[i][j] - log(1 - exp(-2*dist[i][j]));
	 }
     }

   ofstream outfile;
   ofstream hamfile;
   ofstream distfile;
   outfile.open("datalin.dat");
   hamfile.open("hamfile.dat");
   distfile.open("distfile.dat");
   for(int i=0; i<TotNumber; i++)
     {
       for(int j=0; j<TotNumber; j++)
	 {
	   outfile<<datalin[i][j]<<" ";
	   hamfile<<ham[i][j]<<" ";
	   distfile<<dist[i][j]<<" ";
	 }
       outfile<<"\n";
       hamfile<<"\n";
       distfile<<"\n";
     }
   outfile.close();
   hamfile.close();
   distfile.close();

  return 0;
}
