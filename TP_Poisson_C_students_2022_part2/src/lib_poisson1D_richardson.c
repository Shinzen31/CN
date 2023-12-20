/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Function to compute the Richardson iteration
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite){
    int n = *la;
    int max_iterations = *maxit;
    double tolerance = *tol;

    // Initialize temporary vectors
    double *temp = (double *)malloc(n * sizeof(double));
    double *temp2 = (double *)malloc(n * sizeof(double));

    // Perform Richardson iteration
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute A * X
        for (int i = 0; i < n; ++i) {
            temp[i] = 0.0;
            for (int j = fmax(0, i - *kl); j <= fmin(n - 1, i + *ku); ++j) {
                temp[i] += AB[*ku + *kl + 1 + i - j + j * (*lab)] * X[j];
            }
        }

        // Compute the residual: RHS - A * X
        for (int i = 0; i < n; ++i) {
            temp2[i] = RHS[i] - temp[i];
        }

        // Update X with the Richardson formula
        for (int i = 0; i < n; ++i) {
            X[i] += *alpha_rich * temp2[i];
        }

        // Compute the 2-norm of the residual
        double residual_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            residual_norm += temp2[i] * temp2[i];
        }
        residual_norm = sqrt(residual_norm);

        // Store the residual norm in resvec
        resvec[iter] = residual_norm;

        // Check convergence
        if (residual_norm < tolerance) {
            *nbite = iter + 1;  // Number of iterations
            break;
        }
    }

    // Free temporary memory
    free(temp);
    free(temp2);
}

// Function to compute the Jacobi iteration
void jacobi_iteration(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite){
    int n = *la;
    int max_iterations = *maxit;
    double tolerance = *tol;

    // Initialize temporary vectors
    double *temp = (double *)malloc(n * sizeof(double));
    double *temp2 = (double *)malloc(n * sizeof(double));

    // Perform Jacobi iteration
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute A * X
        for (int i = 0; i < n; ++i) {
            temp[i] = 0.0;
            for (int j = fmax(0, i - *kl); j <= fmin(n - 1, i + *ku); ++j) {
                if (i != j) {
                    temp[i] += AB[*ku + *kl + 1 + i - j + j * (*lab)] * X[j];
                }
            }
        }

        // Compute the residual: RHS - A * X
        for (int i = 0; i < n; ++i) {
            temp2[i] = RHS[i] - temp[i];
        }

        // Update X with the Jacobi formula
        for (int i = 0; i < n; ++i) {
            X[i] = temp2[i] / AB[*ku + *kl + 1 + i];
        }

        // Compute the 2-norm of the residual
        double residual_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            residual_norm += temp2[i] * temp2[i];
        }
        residual_norm = sqrt(residual_norm);

        // Store the residual norm in resvec
        resvec[iter] = residual_norm;

        // Check convergence
        if (residual_norm < tolerance) {
            *nbite = iter + 1;  // Number of iterations
            break;
        }
    }

    // Free temporary memory
    free(temp);
    free(temp2);
}

// Function to compute the Gauss-Seidel iteration
void gauss_seidel_iteration(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite){
    int n = *la;
    int max_iterations = *maxit;
    double tolerance = *tol;

    // Initialize temporary vectors
    double *temp = (double *)malloc(n * sizeof(double));
    double *temp2 = (double *)malloc(n * sizeof(double));

    // Perform Gauss-Seidel iteration
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute A * X
        for (int i = 0; i < n; ++i) {
            temp[i] = 0.0;
            for (int j = fmax(0, i - *kl); j <= fmin(n - 1, i + *ku); ++j) {
                if (i != j) {
                    temp[i] += AB[*ku + *kl + 1 + i - j + j * (*lab)] * X[j];
                }
            }
        }

        // Compute the residual: RHS - A * X
        for (int i = 0; i < n; ++i) {
            temp2[i] = RHS[i] - temp[i];
        }

        // Update X with the Gauss-Seidel formula
        for (int i = 0; i < n; ++i) {
            X[i] = temp2[i] / AB[*ku + *kl + 1 + i];
        }

        // Compute the 2-norm of the residual
        double residual_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            residual_norm += temp2[i] * temp2[i];
        }
        residual_norm = sqrt(residual_norm);

        // Store the residual norm in resvec
        resvec[iter] = residual_norm;

        // Check convergence
        if (residual_norm < tolerance) {
            *nbite = iter + 1;  // Number of iterations
            break;
        }
    }

    // Free temporary memory
    free(temp);
    free(temp2);
}


double richardson_alpha_opt(int *la) {
    // Compute the optimal alpha value for Richardson iteration
    // This is a placeholder implementation, adjust as per matrix properties
    double alpha = 0.5; // Example value, adjust based on your problem specifics
    return alpha;
}

void richardson_MB(double *AB, double *RHS, double *SOL, double *MB, int *lab, int *la, int *ku, int *kl, double *alpha, int *maxit, double *resvec, int *nbite) {
    // Richardson iteration method implementation
    // Initialize solution and residual vectors
    memset(SOL, 0, (*la) * sizeof(double));
    double *residual = (double *)malloc((*la) * sizeof(double));
    int iter;
    for (iter = 0; iter < *maxit; ++iter) {
        // Update solution SOL
        for (int i = 0; i < *la; ++i) {
            SOL[i] = SOL[i] + (*alpha) * (RHS[i] - MB[i] * SOL[i]);
        }
        // Compute residual
        double norm = 0.0;
        for (int i = 0; i < *la; ++i) {
            residual[i] = RHS[i] - MB[i] * SOL[i];
            norm += residual[i] * residual[i];
        }
        norm = sqrt(norm);
        resvec[iter] = norm;
        if (norm < 1e-6) { // Convergence criteria, adjust as necessary
            break;
        }
    }
    *nbite = iter;
    free(residual);
}
