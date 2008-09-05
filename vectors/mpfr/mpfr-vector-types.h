/* mpz-vector-types.h */

/**
 * Types used for mpfr-vector that do not belong to the module interface
 * (not present in mpfr-vector.h).
 */
 
 #ifndef MPFR_VECTOR_TYPES_h
 #define MPFR_VECTOR_TYPES_h
 
 /* user must include: gmp.h  
                      mpfr.h
   before this file.
*/

/**
 * The mpfr_vector type holds MPFR float elements.
 * See MPFR documentation for their behaviour.
 */
typedef struct fvector {
  unsigned NbElems;   /** The number of elements*/
  mpfr_t *p_Init;     /** A pointer to the array holding all the vector elements. */
} mpfr_vector;

 #endif /* MPFR_VECTOR_TYPES_h */
 
