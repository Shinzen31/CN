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
  return 0;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}

// Function to compute the error compared to the analytical solution
double compute_error(double *X, double *ANALYTICAL_SOLUTION, int *la) {
    int n = *la;
    double error = 0.0;
    for (int i = 0; i < n; ++i) {
        error += pow(X[i] - ANALYTICAL_SOLUTION[i], 2);
    }
    error = sqrt(error);
    return error;
}

// Function to plot the convergence history using gnuplot
void plot_convergence(double *resvec, int *nbite) {
    FILE *gnuplot = popen("gnuplot -persist", "w");
    if (gnuplot != NULL) {
        fprintf(gnuplot, "set title 'Convergence History'\n");
        fprintf(gnuplot, "set xlabel 'Iteration'\n");
        fprintf(gnuplot, "set ylabel 'Residual Norm'\n");
        fprintf(gnuplot, "plot '-' with lines title 'Residual Norm'\n");
        for (int i = 0; i < *nbite; ++i) {
            fprintf(gnuplot, "%d %lf\n", i + 1, resvec[i]);
        }
        fprintf(gnuplot, "e\n");
        fflush(gnuplot);
        getchar();
        pclose(gnuplot);
    } else {
        printf("Error: Could not open gnuplot.\n");
    }
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
    // Extract matrix for Jacobi iteration from GB format matrix
    for (int i = 0; i < *la; ++i) {
        MB[i] = AB[i * (*lab) + *kv]; // Assuming main diagonal is stored at *kv
    }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
    // Extract matrix for Gauss-Seidel iteration from GB format matrix
    for (int i = 0; i < *la; ++i) {
        MB[i] = 0.0;
        for (int j = fmax(0, i - *kl); j <= fmin(*la - 1, i + *ku); ++j) {
            if (i != j) {
                MB[i] += AB[*ku + *kl + 1 + i - j + j * (*lab)];
            }
        }
    }
}
