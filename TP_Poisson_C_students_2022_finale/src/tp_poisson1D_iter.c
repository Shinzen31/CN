/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include "../include/atlas_headers.h"

int main(int argc,char *argv[])
    /* ** argc: Number of arguments */
    /* ** argv: Values of arguments */
{
    int ierr;
    int jj;
    int nbpoints, la;
    int ku, kl, lab, kv;
    int *ipiv;
    int info;
    int NRHS;
    double T0, T1;
    double *RHS, *SOL, *EX_SOL, *X;
    double *AB;
    double *MB;

    double norm_exsol;
    double norm_sol;

    double temp, relres;

    double opt_alpha;

    struct timespec start, end;
    double elapsed;

    /* Size of the problem */
    NRHS=1;
    if(argc==2)nbpoints=atoi(argv[1]);
    else nbpoints=12;
    la=nbpoints-2;

    /* Dirichlet Boundary conditions */
    T0=5.0;
    T1=20.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS=(double *) malloc(sizeof(double)*la);
    SOL=(double *) calloc(la, sizeof(double)); 
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

    /* Setup the Poisson 1D problem */
    /* General Band Storage */
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    write_vec(RHS, &la, "../data/RHS.dat");
    write_vec(EX_SOL, &la, "../data/EX_SOL.dat");
    write_vec(X, &la, "../data/X_grid.dat");

    kv=0;
    ku=1;
    kl=1;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    /* uncomment the following to check matrix A */
    // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

    /********************************************/
    /* Solution (Richardson with optimal alpha) */

    /* Computation of optimum alpha */
    opt_alpha = richardson_alpha_opt(&la);
    printf("Optimal alpha for simple Richardson iteration is : %lf\n",opt_alpha); 
    printf("Function\t\tErreur\t\tnbite\tTime\n");
    /* Solve */
    double tol=1e-3;
    int maxit=1000;
    double *resvec;
    int nbite=0;

    resvec=(double *) calloc(maxit, sizeof(double));

    /*-----------------Solve with Richardson alpha----------------------*/
    clock_gettime(CLOCK_MONOTONIC_RAW,&start);
    richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    clock_gettime(CLOCK_MONOTONIC_RAW,&end);
    elapsed = ((double) end.tv_sec + (double) end.tv_nsec/1e9) - ((double)start.tv_sec + (double) start.tv_nsec/1e9);

    //------CALCUL D'ERREUR-----------
    // norm_exsol = ||x-x⁰||
    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    cblas_dscal(la,-1.0,RHS,1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    //norm_sol = ||x||
    norm_sol = cblas_dnrm2(la, RHS, 1);
    //relres = ||x-x⁰||/||x||
    relres = norm_sol / norm_exsol;
    printf("Richardson_alpha\t%e\t%d\t%e\n",relres,nbite,elapsed);

    /* Write solution */
    write_vec(SOL, &la, "../data/SOLalpha.dat");

    /* Write convergence history */
    write_vec(resvec, &nbite, "../data/RESVECalpha.dat");

    /*-----------------Richardson General Tridiag-----------------------*/

    /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
    kv = 1;
    ku = 1;
    kl = 1;
    MB = (double *) malloc(sizeof(double)*(lab)*la);

    /*-----------------          Jacobi          -----------------------*/
    memset(SOL,0,la*sizeof(double));
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    nbite=0;
    extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    /* write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "EXTRACT.dat"); */

    /* Solve with General Richardson */
    clock_gettime(CLOCK_MONOTONIC_RAW,&start);
    richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    clock_gettime(CLOCK_MONOTONIC_RAW,&end);
    elapsed = ((double) end.tv_sec + (double) end.tv_nsec/1e9) - ((double)start.tv_sec + (double) start.tv_nsec/1e9);
    //------CALCUL D'ERREUR-----------
    // norm_exsol = ||x-x⁰||
    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    cblas_dscal(la,-1.0,RHS,1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    //norm_sol = ||x||
    norm_sol = cblas_dnrm2(la, RHS, 1);
    //relres = ||x-x⁰||/||x||
    relres = norm_sol / norm_exsol;
    printf("JACOBI\t\t\t%e\t%d\t%e\n",relres,nbite,elapsed);

    /* Write solution */
    write_vec(SOL, &la, "../data/SOLJac.dat");

    /* Write convergence history */
    write_vec(resvec, &nbite, "../data/RESVECJac.dat");

    /*-----------------       Gauss-Seidel       -----------------------*/
    memset(SOL,0,la*sizeof(double));
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    nbite=0;
    extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    /* write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "EXTRACT.dat"); */

    /* Solve with General Richardson */
    clock_gettime(CLOCK_MONOTONIC_RAW,&start);
    richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    clock_gettime(CLOCK_MONOTONIC_RAW,&end);
    elapsed = ((double) end.tv_sec + (double) end.tv_nsec/1e9) - ((double)start.tv_sec + (double) start.tv_nsec/1e9);
    //------CALCUL D'ERREUR-----------
    // norm_exsol = ||x-x⁰||
    norm_exsol = cblas_dnrm2(la, EX_SOL, 1);
    cblas_dscal(la,-1.0,RHS,1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    //norm_sol = ||x||
    norm_sol = cblas_dnrm2(la, RHS, 1);
    //relres = ||x-x⁰||/||x||
    relres = norm_sol / norm_exsol;
    printf("GAUSS-SEIDEL\t\t%e\t%d\t%e\n",relres,nbite,elapsed);

    /* Write solution */
    write_vec(SOL, &la, "../data/SOLGS.dat");

    /* Write convergence history */
    write_vec(resvec, &nbite, "../data/RESVECGS.dat");


    free(resvec);
    free(RHS);
    free(SOL);
    free(EX_SOL);
    free(X);
    free(AB);
    free(MB);
    printf("\n\n--------- End -----------\n");
}
