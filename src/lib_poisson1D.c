/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
      int n=*la;
      int l=*lab;
      int k= *kv;

      int diag=*kv+1;
      int super=diag+1;
      int sub=diag-1;
     
      double h = 1.0 / (n + 1);
      double coef = 1.0 / (h * h);
      
     for (int j = 0; j < n; ++j)
        for (int i = 0; i < l; ++i)
            AB[indexABCol(i, j, lab)] = 0.0;

    
     for (int j = 0; j < n; ++j)
        AB[indexABCol(diag, j, lab)] = 2.0 * coef;

  
    for (int j = 0; j < n; ++j) {
        if (j < n-1) {
            AB[indexABCol(super, j, lab)] = -1.0 * coef;
        }
    }


    for (int j = 0; j < n; ++j) {
        if (j > 0) {
            AB[indexABCol(sub, j, lab)] = -1.0 * coef;
        }
    }
    


}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
      
      int n=*la;
      int l=*lab;
      int k= *kv;

      int diag=*kv+1;
      int super=diag+1;
      int sub=diag-1;
     
      double h = 1.0 / (n + 1);
      double coef = 1.0 / (h * h);
      
     for (int j = 0; j < n; ++j)
        for (int i = 0; i < l; ++i)
            AB[indexABCol(i, j, lab)] = 0.0;

    
     for (int j = 0; j < n; ++j)
        AB[indexABCol(diag, j, lab)] = 1.0 ;

  
    for (int j = 0; j < n; ++j) {
        if (j < n-1) {
            AB[indexABCol(super, j, lab)] = 0.0 ;
        }
    }


    for (int j = 0; j < n; ++j) {
        if (j > 0) {
            AB[indexABCol(sub, j, lab)] = 0.0 ;
        }
    }

}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
   int n = *la;
    double h = 1.0 / (n + 1);
    double coef = 1.0 / (h * h);  

    for (int i = 0; i < n; i++) {
        RHS[i] = 0.0;
    }

  
    RHS[0]     = (*BC0) * coef;   
    RHS[n - 1] = (*BC1) * coef;    
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
    int n = *la;
    double a = *BC0;
    double b = *BC1;

    for (int i = 0; i < n; i++) {
        EX_SOL[i] = a + (b - a) * X[i]; 
    }
}  

void set_grid_points_1D(double* x, int* la){
   int n = *la;
    double h = 1.0 / (n + 1);   // h = 1/(nbpoints-1) with nbpoints = n+2

    for (int i = 0; i < n; i++) {
        x[i] = (i + 1) * h;     // interior points: h, 2h, ..., n*h = 1-h
    }
}

double relative_forward_error(double* x, double* y, int* la){
   int n = *la;
    double num = 0.0;

    for (int i = 0; i < n; ++i) {
        double point_err = x[i] - y[i];
        num += point_err * point_err;
    }

    num = sqrt(num);

    double den = cblas_dnrm2(n, y, 1);

    if (den == 0.0) {
        return -1.0;
    }

    return num / den;
}

int indexABCol(int i, int j, int *lab){
  
  return (i + j*(*lab));
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){

    int N = *la;
    int ldab = *lab; 
  
   
    int diag = *kl+*ku;     // =1+1=2
    int super = diag - 1;     // =1
    int sub   = diag + 1;     // =3

    *info = 0;
   for (int i = 0; i < N; ++i)
        ipiv[i] = i + 1;

    // Thomas algorithm
    for (int i = 1; i < N; ++i)
    {
        double a_i   = AB[indexABCol(sub,  i-1, lab)];   //  a_i
        double u_prev = AB[indexABCol(diag, i-1, lab)]; // u_{i-1}

        if (u_prev == 0.0) {
           *info = i;
           return *info;
        }

       double l_i = a_i / u_prev;
       AB[indexABCol(sub, i-1, lab)] = l_i;

       double c_prev = AB[indexABCol(super, i, lab)];   // c_{i-1}
       double b_i    = AB[indexABCol(diag,  i, lab)];   // b_i

       AB[indexABCol(diag, i, lab)] = b_i - l_i * c_prev;
    }
  
    return *info;
}
