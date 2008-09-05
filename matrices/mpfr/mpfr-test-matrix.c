/* mpfr-test-matrix.c */

/**
 * Test program for matrix utilities
 *                                                           
 * Author:                                                                  
 */

/* includes of system headers */
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h> /* For mpfr-matrix */
#include <mpfr.h>

/* includes of project headers */

/* includes of local headers */
#include "mpfr-matrix.h"

/* Types and constant definitions */

#define VECTOR_SIZE  5
#define DEFAULT_PRECISION 200
#define NUM_ROWS 10
#define NUM_COLUMNS 10

/* Global variables */

/* Functions */

int main(int argc, char ** argv)
{
  mpfr_matrix* mat;
  mpfr_t a_mpfr_int;
  mpfr_t mpfr_array[NUM_COLUMNS];
  mpfr_vector* v1;
  mpfr_vector* v2;
  FILE* the_file;
  unsigned int i;
  /* Test the set/get accessors for a "default" matrix. */ 
  mat = mpfr_matrix_alloc(NUM_ROWS, NUM_COLUMNS, DEFAULT_PRECISION);
  if (mat == NULL) return(1);
  mpfr_init_set_si(a_mpfr_int, 5, GMP_RNDN);
  mpfr_matrix_set_at(mat, 5, 5, a_mpfr_int, GMP_RNDN);
  mpfr_clear(a_mpfr_int);
  mpfr_init(a_mpfr_int);
  mpfr_set(a_mpfr_int, *mpfr_matrix_get_at(mat, 5, 5), GMP_RNDN);
  mpfr_out_str(stdout, 10, 0, a_mpfr_int, GMP_RNDN);
  printf("\n");
  mpfr_matrix_print(stdout, 10, 0, GMP_RNDN, mat);
  mpfr_matrix_free(mat);
  /* Try to read a matrix from a file.*/
  the_file = fopen("matrix.txt", "ro");
  if (the_file == NULL)
  {
    fprintf(stderr, "\n%s: can't open matrix.txt. Aborting the program!\n\n", argv[0]);
    exit(1);
  }
  mat = mpfr_matrix_read_from_file(10, 
                                    GMP_RNDN, 
                                    DEFAULT_PRECISION, the_file);
  fclose(the_file);
  printf("\n");
  mpfr_matrix_print_no_dims(stdout, 10, 0, GMP_RNDN, mat);
  mpfr_matrix_free(mat);
  /* Check the vector <-> matrix communicaton */
  mat = mpfr_matrix_alloc(NUM_ROWS, NUM_COLUMNS, DEFAULT_PRECISION);
  if (mat == NULL) return(1);
  the_file = fopen("matrix.txt", "ro");
  if (the_file == NULL)
  {
    fprintf(stderr, "\n%s: can't open matrix.txt. Aborting the program!\n\n", argv[0]);
    exit(1);
  }
  mat = mpfr_matrix_read_from_file(10, GMP_RNDN, DEFAULT_PRECISION, the_file);
  for (i = 0 ; i < NUM_COLUMNS ; i++)
  {
    mpfr_init_set_ui(mpfr_array[i], 1000 + i, GMP_RNDN);
  }
  v1 = mpfr_vector_from_array(mpfr_array, NUM_COLUMNS, DEFAULT_PRECISION);
  mpfr_matrix_set_row_from_vector(mat, NUM_ROWS - 1, v1);
  printf("\n");
  mpfr_matrix_print_no_dims(stdout, 10, 0, GMP_RNDN, mat);
  v2 = mpfr_matrix_get_vector_from_row(mat, NUM_ROWS - 1, DEFAULT_PRECISION);
  printf("\n");
  mpfr_vector_print(stdout, 10, 10, GMP_RNDN, v2);
  return(0);
}
