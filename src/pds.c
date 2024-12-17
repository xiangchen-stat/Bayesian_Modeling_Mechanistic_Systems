#include <stdio.h>
#include <stdlib.h>

/* check_pds takes in the flattened array of a square matrix C,
 * its one-sided dimension p, and a tolerance (ideally, say, 10^-9).
 */

double* check_sym(double* C, int p, double tol) {
    double F_C = 0;
    int i,j;
    for (i=0; i < p; i++) {
        for (j=0; j < p; j++) {
            // Get the Froebenius norm of C - t(C).
            F_C += (C[i*p + j] - C[j * p + i]) * (C[i*p + j] - C[j * p + i]);
	}
    }
    if (F_C < tol) {
        // Allocate C_new.
	double* C_new = malloc(sizeof(C));
        
	for (i=0; i<p; i++) for (j=0; j<p; j++) C_new[i*p + j] = (C[i*p + j] + C[j*p + i])/2;
	return C_new;
    } else {
        printf("%dx%d input matrix is not symmetric within tolerance %.12f.\n", p, p, tol);
	exit(1);
    }
}
