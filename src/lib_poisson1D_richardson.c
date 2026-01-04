/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include<string.h>

void eig_poisson1D(double* eigval, int *la){
    int n = *la;
    double h=1.0/(n+1);
    for(int i=0; i<n; i++){
        eigval[i]=2.0*(1.0 - cos((i+1)*M_PI*h));//i=k-1
    }




}

double eigmax_poisson1D(int *la){
    int n = *la;
    double h=1.0/(n+1); 
    double max_eig = 2.0*(1.0 - cos(n*M_PI/(n+1)));


    return max_eig;


  return 0;
}

double eigmin_poisson1D(int *la){

  int n = *la;
  double h=1.0/(n+1); 
  double min_eig = 2.0*(1.0 - cos(M_PI/(n+1)));
  
  return min_eig;

   

  return 0;
}

double richardson_alpha_opt(int *la){
  
   double eigmax = eigmax_poisson1D(la);
   double eigmin = eigmin_poisson1D(la);
   double alpha_opt = 2.0/(eigmax + eigmin);
   
   return alpha_opt;



  return 0;
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
      int n = *la;
      double alpha = *alpha_rich;
      double *AX = (double *)malloc(n * sizeof(double));
      double *R = (double *)malloc(n * sizeof(double));
     
      double normb = cblas_dnrm2(n,RHS,1);

      int i;
      for(i=0;i<*maxit;i++){
         cblas_dgbmv(CblasColMajor,CblasNoTrans,n,n,*kl,*ku,1.0,AB,*lab,X,1,0.0,AX,1); // AX = A*X
         // R = RHS - A*X
         for(int j=0;j<n;j++){
             R[j] = RHS[j] - AX[j];
         }
         
         double res = cblas_dnrm2(n,R,1)/ normb; // ||R||_2 / ||b||_2
         resvec[i] = res;


          if(res < *tol){ 
            i++;
            break;
          }

         
          cblas_daxpy(n,alpha,R,1,X,1); // X = X + alpha*R

      }      

      *nbite = i;
      free(AX);
      free(R);


}  

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    
  int n = *la;
  int ldab = *lab;


  memset(MB, 0, (size_t)ldab * (size_t)n * sizeof(double));

   int diag_src = *ku;
   int diag_dst = *kv + 1;

  for (int j = 0; j < n; j++) {
    MB[diag_dst + j * ldab] = AB[diag_src + j * ldab];
  }
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  int n = *la;
  int ldab = *lab;
  
  for(int j=0; j<n; j++){
      for(int i=0; i<*lab; i++){
          MB[i + j*ldab] = 0.0;
      }
  }

  for(int j=0; j<n; j++){
      // diagonal
      MB[(*ku) + j*ldab] = AB[(*ku) + j*ldab];
      // sub-diagonal
      if(j>0){
          MB[(*ku)+1 + j*ldab] = AB[(*ku)+1 + j*ldab];
      }
  }


}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  
    int n = *la;
    int ldab = *lab;
    int ku_val = *ku;
    int kl_val = *kl;

    double *Ax   = (double *)malloc((size_t)n * sizeof(double));
    double *r    = (double *)malloc((size_t)n * sizeof(double));
    double *z    = (double *)malloc((size_t)n * sizeof(double));
    double *diag = (double *)malloc((size_t)n * sizeof(double));
    double *sub  = (double *)malloc((size_t)(n > 1 ? n - 1 : 1) * sizeof(double));

  
    int diag_row = 0;
     while (diag_row < ldab && MB[diag_row] == 0.0) diag_row++;
    if (diag_row >= ldab) diag_row = ku_val; 
    int sub_row = (diag_row + 1 < ldab) ? diag_row + 1 : -1;

    for (int j = 0; j < n; j++) diag[j] = MB[diag_row + j * ldab];

    int has_sub = 0;
      if (sub_row != -1) {
       for (int j = 1; j < n; j++) {
         sub[j - 1] = MB[sub_row + j * ldab];
         if (sub[j - 1] != 0.0) has_sub = 1; 
       }
     }

    double normb = cblas_dnrm2(n, RHS, 1);
    if (normb == 0.0) normb = 1.0;

    int k = 0;
    for (k = 0; k < *maxit; k++) {
    /* Ax = A * x */
    cblas_dgbmv(CblasColMajor, CblasNoTrans, n, n,kl_val, ku_val, 1.0, AB, ldab, X, 1, 0.0, Ax, 1);

    /* r = b - Ax */
    for (int i = 0; i < n; i++) r[i] = RHS[i] - Ax[i];

    const double res = cblas_dnrm2(n, r, 1) / normb;
    resvec[k] = res;
    if (res < *tol) { k++; break; }

 
    if (has_sub) {
      
      z[0] = r[0] / diag[0];
      for (int i = 1; i < n; i++) {
        z[i] = (r[i] - sub[i - 1] * z[i - 1]) / diag[i];
      }
    } 
    else {
    
      for (int i = 0; i < n; i++) z[i] = r[i] / diag[i];
    }

    /* x = x + z */
    cblas_daxpy(n, 1.0, z, 1, X, 1);
   }

  *nbite = k;

  free(Ax);
  free(r);
  free(z);
  free(diag);
  free(sub);
}

void richardson_alpha_csr(CSRMatrix *A, double *b, double *x,double *alpha, double *tol, int *maxit,double *resvec, int *nbite){
    int n = A->n;
    double *Ax = malloc(n * sizeof(double));
    double *r  = malloc(n * sizeof(double));

    double normb = cblas_dnrm2(n, b, 1);

    int k;
    for (k = 0; k < *maxit; k++) {
        dcsrmv(A, x, Ax);

        for (int i = 0; i < n; i++)
            r[i] = b[i] - Ax[i];

        double res = cblas_dnrm2(n, r, 1) / normb;
        resvec[k] = res;

        if (res < *tol) break;

        cblas_daxpy(n, *alpha, r, 1, x, 1);
    }
    *nbite = k;

    free(Ax);
    free(r);
}


void richardson_alpha_csc(CSCMatrix *A, double *b, double *x,double *alpha, double *tol, int *maxit,double *resvec, int *nbite){
    int n = A->n;
    double *Ax = malloc(n * sizeof(double));
    double *r  = malloc(n * sizeof(double));

    double normb = cblas_dnrm2(n, b, 1);

    int k;
    for (k = 0; k < *maxit; k++) {
        dcscmv(A, x, Ax);

        for (int i = 0; i < n; i++)
            r[i] = b[i] - Ax[i];

        double res = cblas_dnrm2(n, r, 1) / normb;
        resvec[k] = res;

        if (res < *tol) break;

        cblas_daxpy(n, *alpha, r, 1, x, 1);
    }
    *nbite = k;

    free(Ax);
    free(r);
}
  



