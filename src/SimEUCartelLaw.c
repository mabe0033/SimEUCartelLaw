#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <R_ext/Callbacks.h>

#define MyEPS sqrt(DOUBLE_EPS)

double rtunifur(double lb1, double rb1, double lb2, double rb2, int truncate,   // simulate (truncated) uniform distribution based on unifrand
                double unifrand) {
  if (truncate) {
    return(fmax2(lb1,lb2)+(fmin2(rb1,rb2)-fmax2(lb1,lb2))*unifrand);
  } else {
    return(lb1+(rb1-lb1)*unifrand);
  }  
}

void update(double p1, double p2, double c, int gg, double rrho, double rphi,   // update result array
            double rksi, double rchi, double rM, double rG, double rA, 
            double* res, int kk, int mm) {
  res[0*mm+kk] = p1;   res[1*mm+kk]  = p2;   res[2*mm+kk] = c;
  res[3*mm+kk] = gg;   res[4*mm+kk]  = rrho; res[5*mm+kk] = rphi;
  res[6*mm+kk] = rksi; res[7*mm+kk]  = rchi; res[8*mm+kk] = rM;
  res[9*mm+kk] = rG;   res[10*mm+kk] = rA;
}

void LGame(double rho, double phi, double ksi, double chi, double B, double G,  // solve legal exemption game
           double A, double *c, double *p1, double *p2, int *gg) {
  if (A <= ksi*chi*B+MyEPS) {                                                   // compare in consideration of machine accuracy
    c[0] = 1.; p1[0] = 0.; p2[0] = 0.; gg[0] = 1;                               // EQ-1
  } else if (A <= (chi+(1.-chi)*(rho-phi))*ksi*B+MyEPS) {                       // compare in consideration of machine accuracy
    double gamma2 = (1.-phi)/(2.-rho-phi);                                      // EQ-3.1
    double beta1 = 1. + (ksi*chi*B-A)/((1.-chi)*(rho-phi)*ksi*B);
    c[0] = gamma2;
    p1[0] = gamma2*ksi*(1.-chi)*(1.-rho)*(1.-beta1);
    p2[0] = (1.-gamma2)*(1.-ksi*(1.-(1.-chi)*(phi+beta1*(1.-phi))));
    gg[0] = 2;
  } else {
    c[0] = 0.; p1[0] = 0.; p2[0] = 1.-ksi; gg[0] = 3;                           // EQ-5
  }
}

SEXP Simulate(SEXP Rho, SEXP Phi, SEXP Ksi, SEXP Chi, SEXP MB, SEXP GB, 
              SEXP AB, SEXP NRuns, SEXP Korr, SEXP QMC) {
  double *rho, *phi, *ksi, *chi, *M, *G, *A, *korr, *res;                       // get parameters
  rho = REAL(Rho); phi = REAL(Phi); ksi = REAL(Ksi); chi = REAL(Chi);
  M   = REAL(MB);  G   = REAL(GB);  A   = REAL(AB);  korr  = REAL(Korr); 
  int nruns = INTEGER(AS_INTEGER(NRuns))[0];
  int qmc   = INTEGER(AS_INTEGER(QMC))[0];
  int nsteps, i, j, k, k2, gg;                                                  // define local variables
  double inc, offset, c, p1, p2, rrho, rphi, rksi, rchi, rM, rG, rA;
  double rand1[7], rand2[7];
  if (qmc) {                                                                    // adjust nruns for QMC and calculate offset and increment
    nsteps = nruns;
    inc = 1.0/nsteps; offset = inc/2.0;
    nruns  = pow(nsteps,7);
  }
  SEXP Res;
  PROTECT (Res = allocVector(REALSXP, 11*nruns));                               // allocate memory for result array
  res = REAL(Res);  

  GetRNGstate();
  for (k = 0; k < nruns; k++) {                                                 // (Q)MC loop
    R_CheckUserInterrupt();
    k2   = k;

    for (i=0;i<7;i++) {                                                         // simulate (correlated) *uniform* distributions based on normal random variates
      if (qmc) {                                                                // for QMC, use grid of 'uniform' numbers and transform to normal random variates
        rand1[i] = qnorm(offset + inc * (k2 % nsteps), 0.0, 1.0, 1, 0);         // this is alg. 22 on p. 169 in Gilli et al. (2011)
        k2 /= nsteps;
      } else rand1[i] = norm_rand();
    }
    for (i=0;i<7;i++) {
      rand2[i] = 0.;
      for (j=0;j<7;j++) rand2[i] += rand1[j] * korr[i*7+j];
      rand2[i] = pnorm(rand2[i],0.,1.,1,0);
    }

    rG   = rtunifur(G[0],G[1],0.,0.,0,rand2[0]);                                // simulate G, A, phi, rho, ksi, chi and M based on (correlated) uniform numbers
    rA   = rtunifur(A[0],A[1],0.,0.,0,rand2[1]);
    rphi = rtunifur(phi[0],phi[1],0.,0.,0,rand2[2]);
    rrho = rtunifur(rho[0],rho[1],rphi,rho[1],1,rand2[3]);                      // nota bene: if rho[0]>=phi[1], there is no truncation
    rksi = rtunifur(ksi[0],ksi[1],0.,0.,0,rand2[4]);
    rchi = rtunifur(chi[0],chi[1],0.,0.,0,rand2[5]);
    rM   = rtunifur(M[0],M[1],0.,0.,0,rand2[6]);

    LGame(rrho,rphi,rksi,rchi,rM,rG,rA,&c,&p1,&p2,&gg);                         // solve the legal exemption game
    update(p1,p2,c,gg,rrho,rphi,rksi,rchi,rM,rG,rA,res,k,nruns);                // update result array
  }
  PutRNGstate();
  UNPROTECT(1);
  return(Res);
}

static R_CallMethodDef callMethods[] = {                                        // Prepare registering native C routines
  {"Simulate", (DL_FUNC) &Simulate, 10},
  {NULL, NULL, 0}
};

void R_init_SimEUCartelLaw(DllInfo *info) {                                     // Register native C routines
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
