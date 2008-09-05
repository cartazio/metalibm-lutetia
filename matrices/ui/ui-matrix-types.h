/* ui-matrix-types.h */

/**
 * Types used for ui_matrix that do not form part of the module interface
 * (not present in ui-matrix.h.
 */
 
 #ifndef UI_MATRIX_TYPES_h
 #define UI_MATRIX_TYPES_h
 
 /* user must include: (empty list)
                      
   before this file.
*/

/**
 * The ui_matix type holds unsigned integer elements.
 *
 * The p member allows for row permutations without actually moving the elements
 * inside the matrix array. As a consequence, element access must always be done
 * through the p member, selecting the right row and, only then, the right column.
 *
 */
typedef struct uimatrix {
  unsigned NbRows;    /** The number of rows*/
  unsigned NbColumns; /** The number of columns*/
  unsigned int **p;   /** A pointer on an array of pointers, each one pointing
                          to the first element of a row. */
  unsigned *p_Init;   /** A pointer to the array holding all the matrix elements. */
  int p_Init_size;    /** needed to free the memory allocated by mpz_init */
} ui_matrix;

 #endif /* UI_MATRIX_TYPES_h */
