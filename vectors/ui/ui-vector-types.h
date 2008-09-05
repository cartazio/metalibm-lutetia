/* ui-vector-types.h */

/**
 * Types used for ui-vector that do not form part of the module interface
 * (not present in ui-vector.h.
 */
 
 #ifndef UI_VECTOR_TYPES_h
 #define UI_VECTOR_TYPES_h
 

/**
 * The ui_vector type holds unsigned int elements.
 */
typedef struct  uivector {
  unsigned int  NbElems;      /* The number of elements*/
  unsigned int  *p_Init;      /* A pointer to the array holding all 
                                 the vector elements. */
} ui_vector;

 #endif /* UI_VECTOR_TYPES_h */
 
