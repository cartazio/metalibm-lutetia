/* test-matrix.c */

/**
 * Test program for matrix utilities
 *                                                           
 * Author:                                                                  
 */

/* includes of system headers */
#include <stdio.h>
#include <gmp.h> /* For mpq-matrix */


/* includes of project headers */
#include "../../vectors/mpq/mpq-vector.h"
#include "../../misc-utils/std-exit-errors.h"

/* includes of local headers */
#include "mpq-matrix.h"

/* Types and constant definitions */

#define VECTOR_SIZE  5
#define RADIX 10

/* Global variables */

/* Functions */

int main(int argc, char ** argv)
{
  mpq_matrix* mat1  = NULL;
  mpq_matrix* mat2  = NULL;
  mpq_matrix* mat3  = NULL;
  mpq_matrix* mat4  = NULL;
  /*FILE*       fmat1 = NULL;*/
  FILE*       fmat2 = NULL;
  char*       vdmf  = "vanderMatrix.txt";
  /* Test matrix allocation and reading from stdin */
  mat1 = mpq_matrix_read_from_stdin(10);
  if (mat1 == NULL)
  {
    fprintf(stderr, "\n%s: Could not read a matrix from stdin!\n", argv[0]);
    return(EX_SOFTWARE);
  }
  mpq_matrix_print(stderr, RADIX, mat1);
  /* Test the swap funtion. */
  mpq_matrix_swap_elements(mat1, 0, 0, 4, 4);
  fprintf(stderr,
          "\n%s: after swapping elements (0,0) and (4,4)\n",
          argv[0]);   
  mpq_matrix_print(stderr, 10, mat1);
  
  mpq_matrix_free(mat1);
  /* Test reading from a file */
  fmat2 = fopen(vdmf, "r");
  if (fmat2 == NULL)
  {
    fprintf(stderr, 
            "\n%s: could not open %s. Aborting the program!\n\n",
            argv[0],
            vdmf);
    return(EX_SOFTWARE);
  }
  mat2 = mpq_matrix_read_from_file(RADIX, fmat2);
  if (mat2 == NULL)
  {
    fprintf(stderr, 
            "\n%s: could not read the matrix in %s. Aborting the program!\n\n",
            argv[0],
            vdmf);
    return(EX_SOFTWARE);
  }
  mat3 = mpq_matrix_clone(mat2);
  if (mat3 == NULL)
  {
    fprintf(stderr, 
            "\n%s: could not clone the matrix read from %s. Aborting the program!\n\n",
            argv[0],
            vdmf);
    return(EX_SOFTWARE);
  }
  if (mpq_matrix_inv_gauss_jordan(mat2))
  {
    fprintf(stderr, 
            "\n%s: could not invert the matrix read from %s. Aborting the program!\n\n",
            argv[0],
            vdmf);
    return(EX_SOFTWARE);
  }
  mat4 = mpq_matrix_naive_prod(mat2, mat3);
  if (mat4 == NULL)
  {
    fprintf(stderr, 
            "\n%s: could make the product between a matrix and it's inverse.\n",
            argv[0]);
    fprintf(stderr,
            "Aborting the program!\n\n");
    return(EX_SOFTWARE);
  }
  mpq_matrix_print(stderr, RADIX, mat4);
  return(0);
}
