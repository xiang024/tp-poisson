/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	      AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
  double temp, relres;
  temp = cblas_ddot(*la, x, 1, x,1);
  temp = sqrt(temp);
  cblas_daxpy(*la, -1.0, x, 1, y, 1);
  relres = cblas_ddot(*la, y, 1, y,1);
  relres = sqrt(relres);
  relres = relres / temp;
  return relres;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
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

void set_CSR_operator_poisson1D(CSRMatrix *A, int *la) {
    int n = *la;
    int nnz = 3*n - 2;

    A->n = n;
    A->nnz = nnz;
    A->values = malloc(nnz * sizeof(double));
    A->col_ind = malloc(nnz * sizeof(int));
    A->row_ptr = malloc((n+1) * sizeof(int));

    int k = 0;
    A->row_ptr[0] = 0;

    for (int i = 0; i < n; i++) {
        if (i > 0) {
            A->values[k] = -1.0;
            A->col_ind[k++] = i - 1;
        }
        A->values[k] = 2.0;
        A->col_ind[k++] = i;

        if (i < n - 1) {
            A->values[k] = -1.0;
            A->col_ind[k++] = i + 1;
        }
        A->row_ptr[i + 1] = k;
    }
}
void set_CSC_operator_poisson1D(CSCMatrix *A, int *la) {
    int n = *la;
    int nnz = 3*n - 2;

    A->n = n;
    A->nnz = nnz;
    A->values = malloc(nnz * sizeof(double));
    A->row_ind = malloc(nnz * sizeof(int));
    A->col_ptr = malloc((n+1) * sizeof(int));

    int k = 0;
    A->col_ptr[0] = 0;

    for (int j = 0; j < n; j++) {
        if (j > 0) {
            A->values[k] = -1.0;
            A->row_ind[k++] = j - 1;
        }
        A->values[k] = 2.0;
        A->row_ind[k++] = j;

        if (j < n - 1) {
            A->values[k] = -1.0;
            A->row_ind[k++] = j + 1;
        }
        A->col_ptr[j + 1] = k;
    }
}

void dcsrmv(const CSRMatrix *A, const double *x, double *y) {
    for (int i = 0; i < A->n; i++) {
        y[i] = 0.0;
        for (int k = A->row_ptr[i]; k < A->row_ptr[i+1]; k++) {
            y[i] += A->values[k] * x[A->col_ind[k]];
        }
    }
}


void dcscmv(const CSCMatrix *A, const double *x, double *y) {
    for (int i = 0; i < A->n; i++) y[i] = 0.0;

    for (int j = 0; j < A->n; j++) {
        for (int k = A->col_ptr[j]; k < A->col_ptr[j+1]; k++) {
            y[A->row_ind[k]] += A->values[k] * x[j];
        }
    }
}

