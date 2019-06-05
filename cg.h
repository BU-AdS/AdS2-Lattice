#ifndef CG_H
#define CG_H
#include <complex>
#include <cstring>
//#include <util.h>
#include <gsl/gsl_fit.h>

using namespace std;

int Mphi(Float *phi, const Float *phi0,
	 const vector<Vertex> NodeList, Param P) {

  int Levels = P.Levels;
  int q = P.q;
  int T = P.t;
  Float msqr = P.msqr;
  Float C_msqr = P.C_msqr;
  bool bc = P.bc;
  int InternalNodes = endNode(Levels-1,P)+1;
  int TotNumber = endNode(Levels,P)+1;
  int offset = TotNumber;
  int T_offset = 0;
  T == 1 ? T_offset = 0 : T_offset = 2;
  
  for(int i=0; i<T*TotNumber; i++) {    
    //cout<<"Initial Mphi pass phi at "<<i<<" is "<<phi[i]<<endl;
    //cout<<"Initial Mphi pass phi0 at "<<i<<" is "<<phi0[i]<<endl;
  }
  
  
  for(int t = 0; t<T; t++) {

    //loop over interior nodes on all disks    
    for(int i = t*offset; i < t*offset + InternalNodes; i++) {
      
      //mass term
      phi[i] = C_msqr*msqr * phi0[i];    
      
      //Spatial links
      for(int mu = 0; mu < q; mu++) {
	//cout<<"i="<<i<<" mu="<<mu<<" phi0="<<phi0[i]<<" phi0["<<NodeList[i].nn[mu]<<"]="<<phi0[NodeList[i].nn[mu]]<<endl;
	phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
      }
      
      //Temporal links
      if(T>1) {
	for(int mu = q; mu < q+T_offset; mu++) {
	  //cout<<"i="<<i<<" mu="<<mu<<" phi0="<<phi0[i]<<" phi0["<<NodeList[i].nn[mu]<<"]="<<phi0[NodeList[i].nn[mu]]<<endl;
	  phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
	}
      }
    }
    
    //Dirichlet or Neuman at Exterior Nodes.  
    for(int i = t*offset + InternalNodes; i < t*offset + TotNumber; i++){
      
      //cout<<"Exterior i="<<i<<" t="<<t<<endl;
      //mass term
      phi[i] = C_msqr*msqr * phi0[i];
      
      //Spatial links
      //the zeroth link is always on the same level.
      phi[i] += phi0[i] - phi0[NodeList[i].nn[0]];
      
      //The q-1 (and q-2) link(s) go back one level,
      //the q-3 (or q-2) link is on the same level. These
      //links are always computed in the same way.
      for(int mu = q-1; mu > (NodeList[i].fwdLinks); mu--) {
	//cout<<"i="<<i<<" mu="<<mu<<" phi0="<<phi0[i]<<" phi0["<<NodeList[i].nn[mu]<<"]="<<phi0[NodeList[i].nn[mu]]<<endl;
	phi[i] += phi0[i] - phi0[NodeList[i].nn[mu]];
      }      
      
      //We use the member data fwdLinks to apply the boundary
      //condition. For Dirichlet, the field value is 0. For
      //Neumann, the derivative is zero.
      for(int mu = NodeList[i].fwdLinks; mu > 0; mu--) {
	if(bc == true) {
	  //Apply Dirchlet BC.
	  phi[i] += phi0[i];
	} else {
	  //Apply Neumann
	  phi[i] += 0.0;
	}
      }
      
      //Temporal links at exterior
      if(T>1) {
	for(int mu = q; mu < q+T_offset; mu++) {
	  //cout<<"i="<<i<<" mu="<<mu<<" phi0="<<phi0[i]<<" phi0["<<NodeList[i].nn[mu]<<"]="<<phi0[NodeList[i].nn[mu]]<<endl;
	  phi[i] += (phi0[i] - phi0[NodeList[i].nn[mu]]);
	}
      }    
    }
  }
  return 0;  
}


Float Minv_phi(Float *phi, Float *b,
	       const vector<Vertex> NodeList, Param P) {
  // CG solutions to Mphi = b 
  //  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  int Levels = P.Levels;
  int diskN = endNode(Levels,P) + 1;
  int TotNumber = P.t*diskN;
  
  Float *res, *resNew, *pvec, *Mpvec, *pvec_tmp;  
  res      = new Float[TotNumber];
  resNew   = new Float[TotNumber];
  pvec     = new Float[TotNumber];
  Mpvec    = new Float[TotNumber];
  pvec_tmp = new Float[TotNumber];

  for(int i=0; i<TotNumber; i++) {
    res[i]      = 0.0;
    resNew[i]   = 0.0;
    pvec[i]     = 0.0;  
    Mpvec[i]    = 0.0;
    pvec_tmp[i] = 0.0;
    //cout<<"phi "<<phi[i]<<" b "<<b[i]<<endl;
  }
  
  Float alpha = 0.0, beta = 0.0, denom = 0.0;
  Float rsq = 0.0, rsqNew = 0.0, bsqrt = 0.0, truersq = 0.0;
  int  i;
  
  for(i = 0; i<TotNumber; i++){
    res[i] = b[i];
    pvec[i] = res[i];
    bsqrt += b[i]*b[i];
  }
  bsqrt = sqrt(bsqrt);
  
  int maxIter = P.MaxIter;
  Float resEpsilon = P.tol;
  // iterate till convergence
  rsqNew = 100.0;
  int k = 0;
  
  while( (k<maxIter)&&(sqrt(rsqNew) > resEpsilon*bsqrt) ){
    
    k++;
    rsq = 0;
    for (int i = 0; i < TotNumber; i++) rsq += res[i]*res[i];
    
    //cout << endl << "In CG at iteration = "<< k <<" Residue Squared  = " << rsq << endl;

    //Mat-Vec operation
    Mphi(Mpvec, pvec, NodeList, P);  
    
    denom = 0;
    for(i=0; i<TotNumber; i++) {
      denom += pvec[i]*Mpvec[i];
      //if(!pvec[i])  cout<<"it's pvec at "<<i<<endl;
      //if(!Mpvec[i]) cout<<"it's Mpvec at "<<i<<endl;
    }
    //cout<<"Denom "<<k<<" = "<<denom<<endl;
    //exit(0);
    
    alpha = rsq/denom;
    
    for(i=0; i < TotNumber; i++) phi[i] +=  alpha * pvec[i];
    for(i=0; i < TotNumber; i++) resNew[i] = res[i]- alpha*Mpvec[i];
    
    // Exit if new residual is small enough
    rsqNew = 0;
    for (i = 0; i < TotNumber; i++) rsqNew += resNew[i]*resNew[i];
    
    // Update vec using new residual
    beta = rsqNew / rsq;
    for (i = 0; i < TotNumber; i++) {
      pvec[i] = resNew[i] + beta * pvec[i];
      res[i] = resNew[i];
    }
  }
  
  if(k == maxIter) {
    printf("CG: Failed to converge iter = %d, rsq = %Le\n", k, (Float)rsq); 
    //  Failed convergence 
  }
  
  Mphi(Mpvec, phi, NodeList, P);  
  for(int i=0; i < TotNumber ; i++) truersq += (Mpvec[i] - b[i])*(Mpvec[i] - b[i]);
  
  //printf("# CG: Converged iter = %d, rsq = %Le, truersq = %Le\n",k,(Float)rsq,(Float)truersq);

  return truersq; // Convergence 
}


void latticeScaling(const vector<Vertex> NodeList, Param p)
{
  cout<<"In!"<<"\n";
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
    cout<<"t is: "<<t<<"\n";
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
      cout<<"i is: "<<i<<"\n";
      //Assume PBC.
      int delta_t = (t2-t1);
      complex<Float> snk = NodeList[i].z;
      Float r   = abs(NodeList[i].z);
      Float r_p = abs(NodeList[j].z);
      
      xi_invariant[s]  = log( ((1-r)*(1-r_p))/(cosh(delta_t)*(1+r)*(1+r_p)
					       - 4*r*r_p*cos(theta)) );
      
      analytic_prop[s] = log( exp(-delta*sigma(src,snk,delta_t)) /
			      (1 - exp(-2*sigma(src,snk,delta_t))));
      
      //if(s<10) cout<<"s="<<s<<" xi="<<xi_invariant[s]<<" ap="<<analytic_prop[s]<<endl;
    }
  }
  
  cout<<"xi invariant is "<<xi_invariant[10]<<"\n";
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
    
    Float* phi_ave = new Float[latVol];
    Float* phi = new Float[latVol];
    Float* b   = new Float[latVol];

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



#endif
