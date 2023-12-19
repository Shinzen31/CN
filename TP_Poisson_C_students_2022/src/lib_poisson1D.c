/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, int* kv) {
    int i;
    for (i = 0; i < *la; ++i) {
        // Diagonale sup��rieure
        if (i != 0)
            AB[i * (*lab) + (*kv) - 1] = -1.0; // -1 sur la diagonale sup��rieure
        // Diagonale principale
        AB[i * (*lab) + (*kv)] = 2.0; // 2 sur la diagonale principale
        // Diagonale inf��rieure
        if (i != (*la) - 1)
            AB[i * (*lab) + (*kv) + 1] = -1.0; // -1 sur la diagonale inf��rieure
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int* la, int* kv) {
    int i;
    for (i = 0; i < *la; ++i) {
        AB[i * (*lab) + (*kv)] = 1.0; // 1 sur la diagonale principale
    }
}


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1) {
    // Appliquer les conditions aux limites de Dirichlet
    RHS[0] = *BC0;
    RHS[*la - 1] = *BC1;
}



void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1) {
    double T0 = *BC0;
    double T1 = *BC1;

    for (int i = 0; i < *la; ++i) {
        double x_i = X[i];
        EX_SOL[i] = T0 + (T1 - T0) * x_i;
    }
}


void set_grid_points_1D(double* x, int* la) {
    int i;
    double h = 1.0 / (*la + 1); // Exemple de calcul de pas de grille
    for (i = 0; i < *la; ++i) {
        x[i] = h * (i + 1); // D��finir les points de la grille
    }
}


void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

int indexABCol(int i, int j, int *lab){
  return 0;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename) {
    FILE* file;
    int i, j;

    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    for (j = 0; j < *la; j++) {
        for (i = 0; i < *lab; i++) {
            fprintf(file, "%lf ", AB[j * (*lab) + i]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

double compute_relative_error(double* numerical, double* analytical, int size) {
    double error = 0.0;
    double norm_numerical = 0.0;
    double norm_analytical = 0.0;

    for (int i = 0; i < size; i++) {
        error += (numerical[i] - analytical[i]) * (numerical[i] - analytical[i]);
        norm_numerical += numerical[i] * numerical[i];
        norm_analytical += analytical[i] * analytical[i];
    }
    error = sqrt(error);
    norm_numerical = sqrt(norm_numerical);
    norm_analytical = sqrt(norm_analytical);

    if (norm_analytical != 0) {
        return error / norm_analytical;
    }
    else {
        return norm_numerical;
    }
}
