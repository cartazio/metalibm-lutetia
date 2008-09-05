/* mpz-vector-types.h */

/**
 * Types used for mpz-vector that do not form part of the module interface
 * (not present in mpz-vector.h.
 */
 
 #ifndef MPZ_VECTOR_TYPES_h
 #define MPZ_VECTOR_TYPES_h
 
 /* user must include: gmp.h  
   before this file.
*/

/**
 * The mpz_vector type holds GMP integer elements.
 * See GMP documentation for their behaviour.
 */
typedef struct zvector {
  unsigned NbElems;   /** The number of elements*/
  mpz_t *p_Init;      /** A pointer to the array holding all the vector elements. */
  int p_Init_size;    /** needed to free the memory allocated by mpfr_init */
} mpz_vector;

 #endif /* MPZ_VECTOR_TYPES_h */
 
