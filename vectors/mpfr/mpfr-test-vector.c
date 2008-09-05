/* mpfr-test-vector.c */

/**
 * Test program for vector utilities
 *                                                           
 * Author:                                                                  
 */

/* includes of system headers */
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h> 
#include <mpfr.h>


/* includes of project headers */

/* includes of local headers */
#include "mpfr-vector.h"

/* Types and constant definitions */

#define VECTOR_SIZE  5
#define DEFAULT_PRECISION 200

/* Global variables */

/* Functions */

int main(int argc, char ** argv)
{
  mpfr_t array1[VECTOR_SIZE], array2[VECTOR_SIZE];
  mpfr_vector* v1 = NULL;
  mpfr_vector* v2 = NULL;
  mpfr_t vector_prod;
  mpfr_t vector_elems_sum;
  int i;
  
  mpfr_set_default_prec(DEFAULT_PRECISION);
  
  for (i = 0 ; i < VECTOR_SIZE ; i++)
  {
    mpfr_init_set_ui(array1[i], 1, GMP_RNDN);  
    mpfr_init_set_ui(array2[i], 2, GMP_RNDN);  
  }
  v1 = mpfr_vector_from_array(array1, VECTOR_SIZE, DEFAULT_PRECISION);
  v2 = mpfr_vector_from_array(array2, VECTOR_SIZE, DEFAULT_PRECISION);
  if ((v1 == NULL) || (v2 == NULL))
  {
    fprintf(stderr, 
            "%s: can not allocate vectors. Aborting the program!\n\n", 
            argv[0]);
    exit(1);
  }
  mpfr_vector_print_no_dims(stdout,
                            10,
                             5,
                            GMP_RNDN,
                            v1);
  mpfr_init(vector_prod);
  mpfr_init(vector_elems_sum);
  mpfr_vector_scalar_product(v1, v2, &vector_prod, DEFAULT_PRECISION);
  fprintf(stdout, "Scalar product -> should be 10: ");
  mpfr_out_str(stdout, 10, 0, vector_prod, GMP_RNDN);
  printf("\n\n");

  mpfr_vector_elems_sum(v1, &vector_elems_sum, GMP_RNDN);
  fprintf(stdout, "Elements sum -> should be 5: ");
  mpfr_out_str(stdout, 10, 0, vector_elems_sum, GMP_RNDN);
  printf("\n\n");
  
  mpfr_vector_set_to_max(v1, v2);
  fprintf(stdout, "All elements should be 2.\n");
  mpfr_vector_print(stdout, 10, 5, GMP_RNDN, v1);
  printf("\n");
  mpfr_vector_free(v1);
  mpfr_vector_free(v2);
  mpfr_clear(vector_prod);
  mpfr_clear(vector_elems_sum);
  
  return(0);
}
