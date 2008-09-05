/* mpq-vector-types.h */

/**
 * Types used for mpq-vector that do not form part of the module interface
 * (not present in mpq-vector.h.
 */
 
 #ifndef MPQ_VECTOR_TYPES_h
 #define MPQ_VECTOR_TYPES_h
 
 /* user must include: gmp.h  
   before this file.
*/

/**
 * The mpfr_vector type holds GMP rational elements.
 * See GMP documentation for their behaviour.
 */
typedef struct qvector {
  unsigned NbElems;   /** The number of elemens*/
  mpq_t *p_Init;      /** A pointer to the array holding all the vector elements. */
  int p_Init_size;    /** needed to free the memory allocated by mpfr_init */
} mpq_vector;

 #endif /* MPQ_VECTOR_TYPES_h */
 
