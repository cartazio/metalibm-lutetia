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
#include "../../vectors/ui/ui-vector.h"
#include "ui-matrix.h"

/* Types and constant definitions */

#define VECTOR_SIZE  5
#define NUM_ROWS 10
#define NUM_COLUMNS 10

/* Global variables */

/* Functions */

int main(int argc, char ** argv)
{
  ui_matrix* mat;
  unsigned int an_int, i;
  unsigned int ui_array[NUM_COLUMNS];
  ui_vector* v1;
  ui_vector* v2;
  FILE* the_file;

  mat = ui_matrix_alloc(NUM_ROWS, NUM_COLUMNS);
  if (mat == NULL) return(1);
  an_int = 5;
  ui_matrix_set_at(mat, 5, 5, an_int);
  an_int = *(ui_matrix_get_at(mat, 5, 5));
  fprintf(stdout, "%u", an_int);
  printf("\n");
  ui_matrix_print(stdout, 10, mat);
  ui_matrix_free(mat);
  /* Try to read a matrix from a file.*/
  the_file = fopen("matrix.txt", "ro");
  if (the_file == NULL)
  {
    fprintf(stderr, "\n%s: can't open matrix.txt. Aborting the program!\n\n", argv[0]);
    exit(1);
  }
  mat = ui_matrix_read_from_file(10, the_file);
  fclose(the_file);
  printf("\n");
  ui_matrix_print_no_dims(stdout, 10, mat);
  ui_matrix_free(mat);
  /* Check the vector <-> matrix communicaton */
  mat = ui_matrix_alloc(NUM_ROWS, NUM_COLUMNS);
  if (mat == NULL) return(1);
  the_file = fopen("matrix.txt", "ro");
  if (the_file == NULL)
  {
    fprintf(stderr, "\n%s: can't open matrix.txt. Aborting the program!\n\n", argv[0]);
    exit(1);
  }
  mat = ui_matrix_read_from_file(10, the_file);
  for (i = 0 ; i < NUM_COLUMNS ; i++)
  {
    ui_array[i] = 1000 + i;
  }
  v1 = ui_vector_from_array(ui_array, NUM_COLUMNS);
  ui_matrix_set_row_from_vector(mat, NUM_ROWS - 1, v1);
  printf("\n");
  ui_matrix_print_no_dims(stdout, 10, mat);
  v2 = ui_matrix_get_vector_from_row(mat, NUM_ROWS - 1);
  printf("\n");
  ui_vector_print(stdout, v2);
  ui_matrix_free(mat);
  ui_vector_free(v1);
  ui_vector_free(v2);
  return(0);
}
