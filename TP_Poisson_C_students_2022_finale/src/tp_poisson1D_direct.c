/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include "../include/atlas_headers.h"

int main(int argc,char *argv[])
    /* ** argc: Nombre d'arguments */
    /* ** argv: Valeur des arguments */
{
    int ierr;
    int jj;
    int nbpoints, la;
    int ku, kl, kv, lab;
    int *ipiv;
    int info;
    int NRHS;
    double T0, T1;
    double *RHS, *EX_SOL, *X, *RHS_TMP;
    double **AAB;
    double *AB;

    double temp, relres;

    double norm_exsol;
    double norm_sol;

    struct timespec start, end;
    double cpu_time_used;

    NRHS=1;
    if(argc==2){
        nbpoints = atoi(argv[1]);
        if(nbpoints<3)nbpoints = 10;
    }
    else nbpoints=10;
    la=nbpoints-2;
    T0=-5.0;
    T1=5.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS=(double *) malloc(sizeof(double)*la);
    RHS_TMP = malloc(la*sizeof(double));
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

    // TODO : you have to implement those functions
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    write_vec(RHS, &la, "RHS.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(X, &la, "X_grid.dat");

    kv=1;
    ku=1;
    kl=1;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
    

    printf("Solution with LAPACK\n");
    printf("Function\tTime\t\tRelres\n");
    //----------------------TEST ON DGBMV----------------------------
    //PROBLEME AVEC DGBMV TO RESOLVE
    cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1.0,AB+1,lab,EX_SOL,1,0.0,RHS_TMP,1);
    write_vec(RHS_TMP, &la, "RHS.dat");
    //On compare le résultat, on devrait retomber sur RHS
    cblas_daxpy(la,-1,RHS,1,RHS_TMP,1);
    for(int i = 0; i < la; i++)
        if (RHS_TMP[i] >1e-10 ) printf("DGBMV Failed to find the good result\n");

    //---------------------TEST with DGBTRF---------------------------
    /* LU Factorization */
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    info=0;
    ipiv = (int *) calloc(la, sizeof(int));
    
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    /* Solution (Triangular) */
    if (info==0){
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
        printf("\n INFO = %d\n",info);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    cpu_time_used = ( (double) end.tv_sec + (double) end.tv_nsec/1e9) - ((double) start.tv_sec + (double)start.tv_nsec/1e9);
    cblas_dcopy(la,RHS,1,RHS_TMP,1);//Copy the result in RHS_TMP

    //Calcul de l'erreur avant
    // norm_exsol = ||x-x⁰||
    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    cblas_dscal(la,-1.0,RHS,1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    //norm_sol = ||x||
    norm_sol = cblas_dnrm2(la, RHS, 1);
    //relres = ||x-x⁰||/||x||
    relres = norm_sol / norm_exsol;
    printf("DGBTRF&DGBTRS\t%f\t%e\n",cpu_time_used,relres);

    //----------------DGBTRFTRIDIAG---------------------
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    info=0;
    /* memset(ipiv,0,sizeof(int)); */

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    /* Solution (Triangular) */
    if (info==0){
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
        printf("\n INFO = %d\n",info);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    cpu_time_used = ( (double) end.tv_sec + (double) end.tv_nsec/1e9) - ((double) start.tv_sec + (double)start.tv_nsec/1e9);
    //comparing DGBTRIDIAG and DGBTRF results
    cblas_daxpy(la,-1,RHS,1,RHS_TMP,1);
    for(int i = 0; i < la; i++)
        if (RHS_TMP[i] != 0) printf("The solution provided by dgbtrftridiag is different from the one produced by dgbtrf\n");
    

    //Calcul de l'erreur avant
    // norm_exsol = ||x-x⁰||
    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    cblas_dscal(la,-1.0,RHS,1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    //norm_sol = ||x||
    norm_sol = cblas_dnrm2(la, RHS, 1);
    //relres = ||x-x⁰||/||x||
    relres = norm_sol / norm_exsol;
    printf("Tridiag&DGBTRS\t%f\t%e\n",cpu_time_used,relres);
    
    //-----------------------------------TEST with DGBSV------------------
    /* It can also be solved with dgbsv */
    // TODO : use dgbsv
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    cpu_time_used = ( (double) end.tv_sec + (double) end.tv_nsec/1e9) - ((double) start.tv_sec + (double)start.tv_nsec/1e9);
    write_xy(RHS, X, &la, "SOL.dat");

    if (info!=0){printf("\n INFO DGBSV = %d\n",info);}


    /* Relative forward error */
    // TODO : Compute relative norm of the residual

    //Calcul de l'erreur avant
    // norm_exsol = ||x-x⁰||
    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    cblas_dscal(la,-1.0,RHS,1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    //norm_sol = ||x||
    norm_sol = cblas_dnrm2(la, RHS, 1);
    //relres = ||x-x⁰||/||x||
    relres = norm_sol / norm_exsol;
    printf("DGBSV   \t%f\t%e\n",cpu_time_used,relres);

    /* printf("\nThe relative forward error is relres = %e\n",relres); */

    free(RHS);
    free(RHS_TMP);
    free(EX_SOL);
    free(X);
    free(AB);
    free(ipiv);
    printf("\n\n--------- End -----------\n");
}
