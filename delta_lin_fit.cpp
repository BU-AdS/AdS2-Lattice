#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_fit.h>

using namespace std;

int main(int argc, char *argv[])
{
  /***** Initialize *****/

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
 
  // Search delta.dat for a shortcut
  bool preTuned = false;
  ifstream fileInput;
  string line;
  char params[256];
  sprintf(params, "%2.4f %d %d", M, L, Q);
  cout<<"The parameters (M,L,Q) you have input are: ";
  printf(params, "%2.4f %d %d", M, L, Q);
  cout<<endl;
  char* search = params;
  unsigned int curLine=0;

  // Open file to search
  fileInput.open("delta.dat");
  if(fileInput.is_open())
    {
      while(getline(fileInput, line))
	{
	  curLine++;
	  //if(line.find(search, 0) != string::npos)
	  cout<<"line is: "<<line<<endl;
	  cout<<"search is: "<<search<<endl;
	  if(line.compare(0,9,search,0,9) == 0)
	    {
	      cout<<"Found data!"<<endl;
	      preTuned = true;
	    }
	}
      if(!preTuned)
	cout<<endl<<"Data not found. Proceeding with generating data."<<endl;
    }
  else cout<<endl<<"Data not found. Proceeding with generating data.."<<endl;
  
  if(preTuned == true)
    cout<<"Done."<<endl;
  
  if(preTuned == false)
    {
      /***** Data reading *****/

      stringstream mstream, lstream, qstream;
      
      mstream<<fixed<<setprecision(6)<<M;   //for naming format
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
     
      vector<vector<float>> dist(TotNumber);
      vector<vector<float>> ham(TotNumber);
      vector<vector<float>> datalin(TotNumber);
      vector<vector<float>> green(TotNumber);
      
      for(int i=0; i<TotNumber; i++)
	{
	  dist[i].resize(TotNumber);
	  ham[i].resize(TotNumber);
	  datalin[i].resize(TotNumber);
	  green[i].resize(TotNumber);
	}
          
      for(int i=0; i<TotNumber; i++)
	{
	  for(int j=0; j<TotNumber; j++)
	    {
	      ham[i][j] = hamvec[TotNumber*i+j];
	      dist[i][j] = distvec[TotNumber*i+j];
	      datalin[i][j] = 0;
	      green[i][j] = greenvec[TotNumber*i+j];
	    }
	}
      
      // Exit if there is no data file 
      // Maybe there is a better implementation using try/catch?
      /*
      if(hamvec.size() == distvec.size() == greenvec.size() == 0)
	{
	  cout<<"Data has not be generated yet. Generate first."<<endl;
	  return -1;
	}
      */

      /***** Linearize *****/
      
      for(int i=0; i<TotNumber; i++)
	{
	  for(int j=0; j<TotNumber; j++)
	    {
	      datalin[i][j] = -dist[i][j] - log(1 - exp(-2*dist[i][j]));
	    }
	}
     
      // Delete infinites      
      for(int i=0; i<TotNumber; i++)
	{     
	  datalin[i].erase( datalin[i].begin() + i );
	  green[i].erase( green[i].begin() + i );
	}
      
      double** datalinpoint = new double*[TotNumber-1];
      double** greenpoint = new double*[TotNumber-1];
      for(int i=0; i<TotNumber-1; i++)
	{
	  datalinpoint[i] = new double[TotNumber-1];
	  greenpoint[i] = new double[TotNumber-1];
	}
      
      for(int i=0; i<TotNumber-1; i++)
	{
	  for(int j=0; j<TotNumber-1; j++)
	    {
	      datalinpoint[i][j] = datalin[i][j];
	      greenpoint[i][j] = log(green[i][j]);
	    }
	}


      /***** Linear fit *****/

      double** c = new double*[2];
      double** cov_ssq = new double*[4];
      
      for(int i=0; i<TotNumber-1; i++)
	{
	  c[i] = new double[2];
	  cov_ssq[i] = new double[4];
	}
 
      double c0_avg = 0.0;
      double c1_avg = 0.0;
      double temp0 = 0.0;
      double temp1 = 0.0;
      
      for(int i=0; i<TotNumber-1; i++)
	{
	  gsl_fit_linear(datalinpoint[i], 1, greenpoint[i], 1, TotNumber-1, &c[0][i], &c[1][i],
			 &cov_ssq[0][i], &cov_ssq[1][i], &cov_ssq[2][i], &cov_ssq[3][i]);
	  temp0 += c[0][i];
	  temp1 += c[1][i];
   	}

      c0_avg = temp0/(TotNumber-1);
      c1_avg = temp1/(TotNumber-1);
      
      cout<<"Average intercept is: "<<c0_avg<<endl;
      cout<<"Average slope (Delta) is: "<<c1_avg<<endl;
      

      /***** Debugging files *****/
      
      /*
      ofstream outfile;
      ofstream hamfile;
      ofstream distfile;
      ofstream greenfile;
      outfile.open("datalin.dat");
      hamfile.open("hamfile.dat");
      distfile.open("distfile.dat");
      greenfile.open("greenfile.dat");
      for(int i=0; i<TotNumber; i++)
	{
	  for(int j=0; j<TotNumber; j++)
	    {
	      outfile<<datalin[i][j]<<" ";
	      hamfile<<ham[i][j]<<" ";
	      distfile<<dist[i][j]<<" ";
	      greenfile<<greenpoint[i][j]<<" ";
	    }
	  outfile<<"\n";
	  hamfile<<"\n";
	  distfile<<"\n";
	  greenfile<<"\n";
	}
      outfile.close();
      hamfile.close();
      distfile.close();
      greenfile.close();
      */
  

      /***** Write results *****/
      
      ofstream file;
      file.open("delta.dat", ios::app);
      file<<fixed<<setprecision(4)<<M<<" "<<L<<" "<<Q<<" "<<c0_avg<<" "<<c1_avg<<endl;
      file.close();
      
      
    }

  return 0;
}
