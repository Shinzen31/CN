#include "lib_poisson1D.h"
#include "atlas_headers.h"

int main(int argc, char *argv[]) {
    int ierr;
    int jj;
    int nbpoints, la;
    int ku, kl, kv, lab;
    int *ipiv;
    int info;
    int NRHS;
    double T0, T1;
    double *RHS, *EX_SOL, *X;
    double *AB;
    double **AAB;
    double *y;
    double temp, relres;

    NRHS = 1;
    nbpoints = 10;
    la = nbpoints - 2;
    T0 = -5.0;
    T1 = 5.0;

    printf("--------- Poisson 1D ---------\n\n");

    RHS = (double *)malloc(sizeof(double) * la);
    EX_SOL = (double *)malloc(sizeof(double) * la);
    X = (double *)malloc(sizeof(double) * la);

    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    write_vec(RHS, &la, "RHS.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(X, &la, "X_grid.dat");

    kv = 1;
    ku = 1;
    kl = 1;
    lab = kv + kl + ku + 1;

    AB = (double *)malloc(sizeof(double) * lab * la);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

    int m = la;
    int n = la;
    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1;
    int incy = 1;
    int lda = lab;
    y = (double *)malloc(la * sizeof(double));

    cblas_dgbmv(CblasColMajor, CblasNoTrans, m, n, kl, ku, alpha, AB, lda, RHS, incx, beta, y, incy);

    relres = compute_relative_error(y, EX_SOL, la);
    printf("Relative error: %e\n", relres);

    printf("Solution with LAPACK\n");
    info = 0;
    ipiv = (int *)malloc(sizeof(int) * la); // Allocate ipiv array

    // Perform LU decomposition
    char trans = 'N';
    LAPACK_dgbtrf(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    if (info == 0) {
        // Solve the linear system using dgbtrs_
        LAPACK_dgbtrs(&trans, &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    }

    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");

    if (info == 0) {
        // Validate LU decomposition by comparing L * U to the original matrix AB
        double *LU = (double *)malloc(sizeof(double) * lab * la);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, la, la, la, 1.0, AB, lda, AB, lda, 0.0, LU, lda);

        // Check if LU matches the original matrix AB
        int validation_passed = 1;
        for (int i = 0; i < lab * la; i++) {
            if (fabs(LU[i] - AB[i]) > 1e-10) {
                validation_passed = 0;
                break;
            }
        }

        if (validation_passed) {
            printf("LU decomposition validation passed.\n");
        } else {
            printf("LU decomposition validation failed.\n");
        }

        free(LU);

        // Solve the linear system again using dgbtrs_
        LAPACK_dgbtrs(&trans, &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        if (info != 0) {
            printf("\n INFO LAPACK_dgbtrs = %d\n", info);
        }
    } else {
        printf("\n INFO = %d\n", info);
    }

    // You can use dgbsv_ if you want to use a different solver, but it's optional
    // dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    // if (info != 0) {
    //     printf("\n ERROR dgbsv = %d\n", info);
    // }

    write_xy(RHS, X, &la, "SOL.dat");

    relres = compute_relative_error(RHS, EX_SOL, la);
    printf("\nThe relative forward error is relres = %e\n", relres);

    free(y);
    free(RHS);
    free(EX_SOL);
    free(X);
    free(AB);
    free(ipiv);

    printf("\n\n--------- End -----------\n");
}
