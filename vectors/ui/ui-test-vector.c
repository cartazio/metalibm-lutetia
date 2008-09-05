/* test-matrix.c */

/**
 * Test program for matrix utilities
 *                                                           
 * Author:                                                                  
 */

/* includes of system headers */
#include <stdlib.h>
#include <stdio.h>

/* includes of project headers */

/* includes of local headers */
#include "ui-vector.h"

/* Types and constant definitions */

#define VECTOR_SIZE  5

/* Global variables */

/* Functions */

int main(int argc, char ** argv)
{
  unsigned int array1[VECTOR_SIZE], array2[VECTOR_SIZE];
  ui_vector* v1 = NULL;
  ui_vector* v2 = NULL;
  unsigned int vector_prod;
  int i;
  
  
  for (i = 0 ; i < VECTOR_SIZE ; i++)
  {
    array1[i] = 1;  
    array2[i] = 2;  
  }
  v1 = ui_vector_from_array(array1, VECTOR_SIZE);
  v2 = ui_vector_from_array(array2, VECTOR_SIZE);
  if ((v1 == NULL) || (v2 == NULL))
  {
    fprintf(stderr, "%s: can not allocate vectors. Aborting the program!\n\n", argv[0]);
    exit(1);
  }
  ui_vector_scalar_product(v1, v2, &vector_prod);
  fprintf(stdout, "Should be 10: ");
  fprintf(stdout, "%u", vector_prod);
  fprintf(stdout, "\n");
  ui_vector_free(v1);
  ui_vector_free(v2);
  return(0);
}
