#ifndef ARP_INT_ADS_H
#define ARP_INT_ADS_H


//====================================================================================
//
// Thu Aug 23 Dean Howarth
//  
// ARPACK interafce for AdS2-Lattice.
// 
//====================================================================================

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <complex>
#include <cg.h>

using namespace std;

void arpackErrorHelpSAUPD(int *iparam_);
void arpackErrorHelpSEUPD(int *iparam_);

#define ARPACK(s) s ## _

#ifdef __cplusplus
extern "C" {
#endif

//  Interface functions to the external ARPACK library. These functions utilize 
//  ARPACK's implemntation of the Implicitly Restarted Arnoldi Method to compute a 
//  number of eigenvectors/eigenvalues with user specified features, such as those 
//  with small real part, small magnitude etc. Parallel (OMP/MPI) versions
//  are also supported.

  
  //Serial, single prec real eigenvectors
  extern int ARPACK(ssaupd) (int *ido, char *bmat, int *n, char *which, int *nev,
			     float *tol, float *resid, int *ncv,
			     float *v, int *ldv, int *iparam, int *ipntr,
			     float *workd, float *workl,
			     int *lworkl, int *info);

  //Serial, double prec real eigenvectors
  extern int ARPACK(dsaupd)(int *ido, char *bmat, int *n, char *which, int *nev,
			    double *tol, double *resid, int *ncv,
			    double *v, int *ldv, int *iparam, int *ipntr,
			    double *workd, double *workl, 
			    int *lworkl, int *info);
  
  //Serial, single prec real eigenvalues
  extern int ARPACK(sseupd) (int *comp_evecs, char *howmany, int *select,
			     float *evals, float *v,
			     int *ldv, float *sigma,
			     float *workev, char *bmat, int *n,
			     char *which, int *nev, float *tol,
			     float *resid, int *ncv,
			     float *v1, int *ldv1, int *iparam,
			     int *ipntr, float *workd,
			     float *workl, int *lworkl,
			     int *info);			
  
  //Serial, double prec real eigenvalues
  extern int ARPACK(dseupd) (int *comp_evecs, char *howmany, int *select,
			     double *evals, double *v,
			     int *ldv, double *sigma,
			     char *bmat, int *n,
			     char *which, int *nev, double *tol,
			     double *resid, int *ncv,
			     double *v1, int *ldv1, int *iparam,
			     int *ipntr, double *workd,
			     double *workl, int *lworkl,
			     int *info);
  
  
  extern int ARPACK(mcinitdebug)(int*,int*,int*,int*,int*,int*,int*,int*);
    
  //ARPACK initlog and finilog routines for printing the ARPACK log  
  extern int ARPACK(initlog) (int*, char*, int);
  extern int ARPACK(finilog) (int*);
  
#ifdef __cplusplus
}
#endif


void arpackErrorHelpSAUPD(int *iparam_) {
  printf("\nError help NAUPD\n\n");
  printf("INFO Integer.  (INPUT/OUTPUT)\n");
  printf("     If INFO .EQ. 0, a randomly initial residual vector is used.\n");
  printf("     If INFO .NE. 0, RESID contains the initial residual vector,\n");
  printf("                        possibly from a previous run.\n");
  printf("     Error flag on output.\n");
  printf("     =  0: Normal exit.\n");
  printf("     =  1: Maximum number of iterations taken.\n");
  printf("        All possible eigenvalues of OP has been found. IPARAM(5)\n");
  printf("        returns the number of wanted converged Ritz values.\n");
  printf("     =  2: No longer an informational error. Deprecated starting\n");
  printf("        with release 2 of ARPACK.\n");
  printf("     =  3: No shifts could be applied during a cycle of the\n");
  printf("        Implicitly restarted Arnoldi iteration. One possibility\n");
  printf("        is to increase the size of NCV relative to NEV.\n");
  printf("        See remark 4 below.\n");
  printf("     = -1: N must be positive.\n");
  printf("     = -2: NEV must be positive.\n");
  printf("     = -3: NCV-NEV >= 1 and less than or equal to N.\n");
  printf("     = -4: The maximum number of Arnoldi update iteration\n");
  printf("        must be greater than zero.\n");
  printf("     = -5: WHICH must be 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
  printf("     = -6: BMAT must be one of 'I' or 'G'.\n");
  printf("     = -7: Length of private work array is not sufficient.\n");
  printf("     = -8: Error return from LAPACK eigenvalue calculation;\n");
  printf("     = -9: Starting vector is zero.\n");
  printf("     = -10: IPARAM(7) must be 1,2,3.\n");
  printf("     = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
  printf("     = -12: IPARAM(1) must be equal to 0 or 1.\n");
  printf("     = -9999: Could not build an Arnoldi factorization.\n");
  printf("        User input error highly likely.  Please\n");
  printf("        check actual array dimensions and layout.\n");
  printf("        IPARAM(5) returns the size of the current Arnoldi\n");
  printf("        factorization.\n");
  printf("        iparam_[5] = %d\n", iparam_[4]);
}

void arpackErrorHelpSEUPD(int *iparam_) {
  printf("\nError help NEUPD\n\n");
  printf("INFO Integer.  (OUTPUT)\n");
  printf("     Error flag on output.\n");
  printf("     =  0: Normal exit.\n");
  printf("     =  1: The Schur form computed by LAPACK routine csheqr\n");
  printf("        could not be reordered by LAPACK routine ztrsen.\n");
  printf("        Re-enter subroutine zneupd with IPARAM(5)=NCV and\n");
  printf("        increase the size of the array D to have\n");
  printf("        dimension at least dimension NCV and allocate at\n");
  printf("        least NCV\n");
  printf("        columns for Z. NOTE: Not necessary if Z and V share\n");
  printf("        the same space. Please notify the authors if this\n");
  printf("        error occurs.\n");
  printf("     = -1: N must be positive.\n");
  printf("     = -2: NEV must be positive.\n");
  printf("     = -3: NCV-NEV >= 1 and less than or equal to N.\n");
  printf("     = -5: WHICH must be 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
  printf("     = -6: BMAT must be one of 'I' or 'G'.\n");
  printf("     = -7: Length of private work WORKL array is inufficient.\n");
  printf("     = -8: Error return from LAPACK eigenvalue calculation.\n");
  printf("        This should never happened.\n");
  printf("     = -9: Error return from calculation of eigenvectors.\n");
  printf("        Informational error from LAPACK routine ztrevc.\n");
  printf("     = -10: IPARAM(7) must be 1,2,3\n");
  printf("     = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
  printf("     = -12: HOWMNY = 'S' not yet implemented\n");
  printf("     = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.\n");
  printf("     = -14: ZNAUPD did not find any eigenvalues to sufficient\n");
  printf("        accuracy.\n");
  printf("     = -15: ZNEUPD got a different count of the number of\n");
  printf("        converged Ritz values than ZNAUPD got. This\n");
  printf("        indicates the user probably made an error in\n");
  printf("        passing data from ZNAUPD to ZNEUPD or that the\n");
  printf("        data was modified before entering ZNEUPD\n");
  printf("        iparam_[5] = %d\n", iparam_[4]);
}

void arpack_solve(vector<Vertex> NodeList, Param p) {

  //Construct parameters and memory allocation
  //------------------------------------------
  
  // all FORTRAN communication uses underscored 
  int ido_;
  int info_;
  int *ipntr_ = (int*)malloc(11*sizeof(int));
  int *iparam_ = (int*)malloc(11*sizeof(int));
  int n_    = endNode(p.Levels,p) + 1,
    nev_    = n_-100,
    nkv_    = n_,
    ldv_    = n_,
    lworkl_ = (nkv_ * (nkv_ + 8)),
    rvec_   = 1;
  int max_iter = p.arpackMaxiter;

  double tol_ = p.arpackTol;

  //ARPACK workspace
  double sigma_ = 0.0;
  double *resid_ = (double *) malloc(ldv_*sizeof(double));
  double *w_workd_ = (double *) malloc(3*ldv_*sizeof(double));
  double *w_workl_ = (double *) malloc(lworkl_*sizeof(double)); 
  double *w_workev_= (double *) malloc(3*nkv_*sizeof(double));
  double *w_rwork_ = (double *) malloc(nkv_*sizeof(double));    
  int *select_ = (int*)malloc(nkv_*sizeof(int));
  
  double *evecs = (double *) malloc(nkv_*n_*sizeof(double));
  double *evals = (double *) malloc(nkv_   *sizeof(double));
  
  for(int n=0; n<nkv_; n++) {
    evals[n] = 0;
    for(int i=0; i<n_; i++) {
      evecs[n*n_ + i] = 0;
    }
  }
  
  //Alias pointers
  double *evecs_ = nullptr;
  evecs_ = (double*)(evecs);    
  double *evals_ = nullptr;
  evals_ = (double*)(evals);
  
  //Memory checks
  if((iparam_ == nullptr) ||
     (ipntr_ == nullptr) || 
     (resid_ == nullptr) ||  
     (w_workd_ == nullptr) || 
     (w_workl_ == nullptr) ||
     (w_workev_ == nullptr) ||
     (w_rwork_ == nullptr) || 
     (select_ == nullptr) ) {
    printf("eigenSolver: not enough memory for ARPACK workspace.\n");
    exit(0);
  }    

  //Assign values to ARPACK params 
  ido_        = 0;
  info_       = 0;
  iparam_[0]  = 1;
  iparam_[2]  = max_iter;
  iparam_[3]  = 1;
  iparam_[6]  = 1;
  
  //ARPACK problem type to be solved
  char howmany = 'P';
  char bmat = 'I';
  char spectrum[2] = {'S','A'};
  int iter_cnt= 0;

  //Start ARPACK routines
  //---------------------------------------------------------------------------------
 
  double *psi1;
  double *psi2;

  double *psi1_cpy = (double*)malloc(n_*sizeof(double));
  double *psi2_cpy = (double*)malloc(n_*sizeof(double));
  
  for(int i=0; i<n_; i++) {
    psi1_cpy[i] = 1.0;
    psi2_cpy[i] = 1.0;
  }
  
  psi1 = w_workd_;
  psi2 = w_workd_ + n_;
  
  double t1;
  double time = 0.0;;
  do {
    
    t1 = -((double)clock());
    
    //Interface to arpack routines
    //----------------------------
    
    ARPACK(dsaupd)(&ido_, &bmat, &n_, spectrum, &nev_, &tol_, resid_, &nkv_,
		   evecs_, &ldv_, iparam_, ipntr_, w_workd_, w_workl_, &lworkl_,
		   &info_);

    if (info_ != 0) {
      printf("\nError in dsaupd info = %d. Exiting...\n",info_);
      arpackErrorHelpSAUPD(iparam_);
      exit(0);
    }
    
    if (ido_ == 99 || info_ == 1)
      break;
    
    if (ido_ == -1 || ido_ == 1) {

      //Copy from Arpack workspace
      for(int i=0; i<n_; i++) {
	psi1_cpy[i] = *(psi1 + i);
      }
      
      //Apply matrix-vector operation
      Mphi(psi2_cpy, psi1_cpy, NodeList, p);
      
      //Copy to Arpack workspace
      for(int i=0; i<n_; i++) {
	*(psi2 + i) = psi2_cpy[i];
      }
    }
    
    t1 += clock();
    time += t1;
    printf("Arpack Iteration: %d (%e secs)\n", iter_cnt, time/(CLOCKS_PER_SEC));
    iter_cnt++;
    
  } while (99 != ido_ && iter_cnt < max_iter);
  
  //Subspace calulated sucessfully. Compute nEv eigenvectors and values
  printf("ARPACK Finished in %e secs: iter=%04d  info=%d  ido=%d\n", time/(CLOCKS_PER_SEC), iter_cnt, info_, ido_);

  printf("ARPACK Computing Eigenvlaues\n");
  ARPACK(dseupd)(&rvec_, &howmany, select_, evals_, evecs_, &n_, &sigma_,
		 &bmat, &n_, spectrum, &nev_, &tol_,
		 resid_, &nkv_, evecs_, &n_, iparam_, ipntr_, w_workd_,
		 w_workl_, &lworkl_, &info_);
  if (info_ == -15) {
    printf("\nError in dseupd info = %d. You likely need to\n"
	   "increase the maximum ARPACK iterations. Exiting...\n", info_);
    arpackErrorHelpSEUPD(iparam_);
    exit(0);
  } else if (info_ != 0) {
    printf("\nError in dseupd info = %d. Exiting...\n", info_);
    arpackErrorHelpSEUPD(iparam_);
  }

  // Print additional convergence information.
  if(info_ == 1){
    printf("Maximum number of iterations reached.\n");
  }
  else{
    if(info_ == 3){
      printf("Error: No shifts could be applied during implicit\n");
      printf("Error: Arnoldi update, try increasing NkV.\n");
    }
  }
  
  //Print Evalues  
  for(int i=0; i<nev_ ;i++){
    printf("%04d %+.16e\n", i, evals_[i]);    
  }
  
  t1 += clock();
  printf("\n*************************************************\n");
  printf("%d Eigenvalues of hamiltonian computed in: %f sec\n", nev_, t1/(CLOCKS_PER_SEC));
  printf("Total time spent in ARPACK: %f sec\n", (time+t1)/(CLOCKS_PER_SEC));
  printf("*************************************************\n");
  
  // cleanup 
  free(ipntr_);
  free(iparam_);
  free(resid_);
  free(w_workd_);
  free(w_workl_);
  free(w_workev_);
  free(w_rwork_);
  free(select_);
  
  return;
}

#endif
