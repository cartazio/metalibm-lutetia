/* mpz-matrix-types.h */

/**
 * Types used for mpz-matrix that do not form part of the module interface
 * (not present in mpz-matrix.h.
 */
 
 #ifndef MPQ_MATRIX_TYPES_h
 #define MPQ_MATRIX_TYPES_h
 
 /* user must include: gmp.h  
                      mpfr.h
   before this file.
*/

/**
 * The mpq_matrix type holds GMP rational elements.
 * See GMP documentation for their behaviour.
 *
 * The p member allows for row permutations without actually moving the elements
 * inside the matrix array. As a consequence, element access must always be done
 * through the p member, selecting the right row and, only then, the right column.
 *
 */
typedef struct qmatrix {
  unsigned NbRows;    /** The number of rows*/
  unsigned NbColumns; /** The number of columns*/
  mpq_t **p;          /** A pointer on an array of pointers, each one pointing to the first
                          of a row of elements. */
  mpq_t *p_Init;      /** A pointer to the array holding all the matrix elements. */
  int p_Init_size;    /** needed to free the memory allocated by mpz_init */
} mpq_matrix;

 #endif /* MQZ_MATRIX_TYPES_h */
