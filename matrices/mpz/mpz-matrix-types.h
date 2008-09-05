/* mpz-matrix-types.h */

/**
 * Types used for mpz-matrix that do not form part of the module interface
 * (not present in mpz-matrix.h.
 */
 
 #ifndef MPZ_MATRIX_TYPES_h
 #define MPZ_MATRIX_TYPES_h
 
 /* user must include: gmp.h  
                      mpfr.h
   before this file.
*/

/**
 * The mpz_matix type holds GMP integer elements.
 * See GMP documentation for their behaviour.
 *
 * The p member allows for row permutations without actually moving the elements
 * inside the matrix array. As a consequence, element access must always be done
 * through the p member, selecting the right row and, only then, the right column.
 *
 */
typedef struct zmatrix {
  unsigned NbRows;    /** The number of rows*/
  unsigned NbColumns; /** The number of columns*/
  mpz_t **p;          /** A pointer on an array of pointers, each one pointing to the first
                          of a row of elements. */
  mpz_t *p_Init;      /** A pointer to the array holding all the matrix elements. */
  int p_Init_size;    /** needed to free the memory allocated by mpz_init */
} mpz_matrix;

 #endif /* MPZ_MATRIX_TYPES_h */
