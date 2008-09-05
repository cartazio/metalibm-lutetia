/* test-matrix.c */

/**
 * Test program for matrix utilities
 *                                                           
 * Author:                                                                  
 */

/* includes of system headers */
#include <stdio.h>

/* includes of project headers */

/* includes of local headers */
#include "si-vector.h"

/* Types and constant definitions */

#define VECTOR_SIZE  5

/* Global variables */

/* Functions */

int main(int argc, char ** argv)
{
  int array_1[VECTOR_SIZE] = { 7, -2, 1, 12, -10};
  si_vector* vector_1;
  si_vector* abs_vector_1;
  /* Test vector elements sort */
  vector_1 = si_vector_from_array(array_1, VECTOR_SIZE);
  si_vector_print_no_dims(stderr, vector_1);
  /*si_vector_shift(vector_1, 4, SI_VECTOR_DOWN);*/
  si_vector_insertion_sort(vector_1, SI_VECTOR_DOWN);
  si_vector_print_no_dims(stderr, vector_1);
  /* Test compute abs of vector elements */
  abs_vector_1 = si_vector_get_abs(vector_1);
  si_vector_print_no_dims(stderr, abs_vector_1);
  si_vector_insertion_sort(abs_vector_1, SI_VECTOR_UP);
  si_vector_print_no_dims(stderr, abs_vector_1);
  /* Test si_vector_are_equal */
  if (si_vector_are_equal(vector_1, abs_vector_1))
  {
    si_vector_print_no_dims(stderr, vector_1);
    si_vector_print_no_dims(stderr, abs_vector_1);
    fprintf(stderr, "Come on, these vectors can't be equal!\n");
  }
  if (! si_vector_are_equal(vector_1, vector_1))
  {
    si_vector_print_no_dims(stderr, vector_1);
    si_vector_print_no_dims(stderr, vector_1);
    fprintf(stderr, "Come on, these vectors can't be different!\n");
  }
  return(0);
}
