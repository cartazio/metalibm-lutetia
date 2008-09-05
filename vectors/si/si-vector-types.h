/* si-vector-types.h */

/**
 * Types used for si-vector that do not form part of the module interface
 * (not present in si-vector.h.
 */
 
 #ifndef SI_VECTOR_TYPES_h
 #define SI_VECTOR_TYPES_h
 

/**
 * The mpfr_vector type holds signed int elements.
 */
typedef struct  sivector {
  unsigned    NbElems;      /* The number of elements*/
  signed int  *p_Init;      /* A pointer to the array holding all 
                               the vector elements. */
  int         p_Init_size;  /* not really needed here but used 
                               for consistence with other "vectors"  */
} si_vector;

 #endif /* SI_VECTOR_TYPES_h */
 
