/* mpfr-matrix-types.h */
/******************************************************************************/

/** @file
 * Types used for mpfr-matrix that do not form part of the module interface
 * (not present in mpfr-matrix.h).
 *
 * @author ST 
 * @date 2005-07-21
 *
 */
 
 #ifndef MPFR_MATRIX_TYPES_h
 #define MPFR_MATRIX_TYPES_h
 
 /* user must include: gmp.h  
                       mpfr.h
   before this file.
  */

/**
 * The mpfr_matrix type holds MPFR float elements.
 * See MPFR documentation for their behaviour.
 *
 * The p member allows for row permutations without actually moving the elements
 * inside the matrix array. As a consequence, element access must always be done
 * through the p member, selecting the right row and, only then, the right column.
 *
 */
typedef struct fmatrix {
  unsigned NbRows;    /** The number of rows*/
  unsigned NbColumns; /** The number of columns*/
  mpfr_t **p;         /** A pointer on an array of pointers, each one pointing to the first
                          of a row of elements. */
  mpfr_t *p_Init;     /** A pointer to the array holding all the matrix elements. */
  int p_Init_size;    /** needed to free the memory allocated by mpfr_init */
} mpfr_matrix;

 #endif /* MPFR_MATRIX_TYPES_h */
 
