/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){ // lambda_max
    double h = 1.0 / ((*la) +1.0);
    return 4*sin(((*la) * M_PI * h)/2)*sin(((*la) * M_PI * h)/2);
}

double eigmin_poisson1D(int *la){ // lambda_min
    double h = 1.0 / ((*la) +1.0);
    return 4*sin((M_PI * h)/2)*sin(( M_PI * h)/2);
}

double richardson_alpha_opt(int *la){
    // return 0.5 | Démonstration dans le rapport
    return 2/(eigmax_poisson1D(la)+eigmin_poisson1D(la)); 
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    double *B = malloc(*la* sizeof(double));
    double norme_B = cblas_dnrm2(*la,RHS,1);

    for((*nbite) = 0; (*nbite) < *maxit ;(*nbite)++){
        //Copie de RHS dans B
        cblas_dcopy(*la,  RHS, 1, B,1);
        //b = b - Ax
        cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,B,1);
        //Calcul du residu
        resvec[*nbite] = cblas_dnrm2(*la,B,1) / norme_B;
        //x = x + alpha * b | x = x + alpha*(b-Ax)
        cblas_daxpy(*la,*alpha_rich,B,1,X,1);
        if(resvec[*nbite]<=*tol) break;
    }
    free(B);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    //M = D
    for (int i = 0; i < *la; i++)
        MB[i*(*lab)+1] = AB[i*(*lab)+1];
}
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
    //M = D - E 
    for(int i = 0; i < *la; i++){
        MB[*lab*i+1] = AB[*lab*i+1]; // M = D
        MB[*lab*i+2] = AB[*lab*i+2]; // M = M - E
    }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    double *B = malloc(*la * sizeof(double));
    double norme_B = cblas_dnrm2(*la,RHS,1);
    int * ipiv = malloc(*la * sizeof(int));
    int info = 0;
    int NRHS = 1;
    int ku_minus = *ku-1; // nécessaire pour MB

    // LU factorization of MB
    dgbtrf_(la, la, kl, &ku_minus, MB, lab, ipiv, &info);
    for((*nbite) = 0; (*nbite) < *maxit ;(*nbite)++){
        //Copie de RHS pour garder le même RHS à chaque iteration
        cblas_dcopy(*la,  RHS, 1, B,1);
        //b = b - Ax
        cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,-1.0,AB,*lab,X,1,1.0,B,1);
        //Calcul du residu
        resvec[*nbite] = cblas_dnrm2(*la,B,1) / norme_B;
        //b = b/M | b = (b - Ax)/M
        dgbtrs_("N", la, kl, &ku_minus, &NRHS, MB, lab, ipiv, B, la, &info);
        // x = x + b | x = x + (b - Ax)/M
        cblas_daxpy(*la,1,B,1,X,1);
        if(resvec[*nbite]<=*tol) break;
    }
    free(B);
    free(ipiv);
}
