/* test-matrix.c */

/**
 * Test program for matrix utilities
 *                                                           
 * Author:                                                                  
 */

/* includes of system headers */
#include <stdio.h>
#include <gmp.h> /* For mpz-matrix */

/* includes of project headers */

/* includes of local headers */
#include "mpz-matrix.h"

/* Types and constant definitions */

/* Global variables */

/* Functions */

int main(int argc, char ** argv)
{
  mpz_matrix* mat;
  mpz_t a_mpz_int;
  mat = mpz_matrix_alloc(10,10);
  if (mat == NULL) return(1);
  mpz_init_set_si(a_mpz_int, 5);
  mpz_matrix_set_at(mat, 5, 5, a_mpz_int);
  mpz_clear(a_mpz_int);
  mpz_init(a_mpz_int);
  mpz_set(a_mpz_int, *mpz_matrix_get_at(mat, 5, 5));
  mpz_out_str(stdout, 10, a_mpz_int);
  printf("\n");
  mpz_matrix_print(stdout, 10, mat);
  mpz_matrix_free(mat);
  return(0);
}
