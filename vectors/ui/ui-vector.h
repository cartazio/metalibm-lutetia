/* ui_vector.h */
/******************************************************************************/

#ifndef UI_VECTOR_h
#define UI_VECTOR_h

#define UI_VECTOR_UP 0
#define UI_VECTOR_DOWN 1

#include "ui-vector-types.h"

/**
 * Create a new ui vector, allocating memory, initializing
 * all the elements and setting them to 0.
 * @param NbElems - the number of elements;
 * @return a pointer on the vector if all parameters have acceptable values 
 * and if all allocation and initialization operations succeed, 
 * returns NULL otherwise.
 */
ui_vector*
ui_vector_alloc(unsigned NbElems);

ui_vector*
ui_vector_from_array(unsigned int* si_array, unsigned int NbElems);

void
ui_vector_free(ui_vector *Vect);

void
ui_vector_print(FILE *Dst, 
                ui_vector *Vect);

void 
ui_vector_print_no_dims(FILE *Dst, 
                        ui_vector *Vect);

void 
ui_vector_print_matlab(FILE *Dst, 
                        ui_vector *Vect);


ui_vector* 
ui_vector_read_from_stdin(int radix);

void 
ui_vector_read_input_from_file(ui_vector *Vect,
                                int radix, 
                                FILE *Scr);

ui_vector*
ui_vector_read_from_file(int radix, 
                          FILE *Scr);

unsigned int*
ui_vector_get_at(ui_vector* vect, 
                  unsigned int pos);

void 
ui_vector_set_at(ui_vector* vect, 
                  unsigned int pos, 
                  unsigned int value);

/**
 * Check whether there are duplicate values in a vector
 * 
 * @param vect - the vector to check;
 * @return 1 if there are duplicates, 0 otherwise
 */
int
ui_vector_has_duplicates(ui_vector* vect);

unsigned int 
ui_vector_get_size(ui_vector* vect);

/**
 * Shift an element of the vector of one position, overwriting
 * the neigbouring element.
 *
 * No return value.
 *
 * @param vect      - the vector;
 * @param index     - the index of the element to shift;
 * @param direction - the direction of the shift: UI_VECTOR_UP -> to an higher index
 *                                                UI_VECTOR_DOWN -> to a lower index.
 */
void
ui_vector_shift(ui_vector* vect, 
      unsigned int index, 
      int direction);

/**
 * Order the elements of the vector.
 *
 * Sorting is made in-place (the original vector is altered) using 
 * the insertion sort alogorithm and there is no return value.
 * 
 *
 * @param - vect, the vector to order
 * @param - direction, the sorting direction: UI_VECTOR_UP ->  increasing
 *                                            UI_VECTOR_DOWN -> decreasing.
 */
void
ui_vector_insertion_sort(ui_vector* vi, int direction);

/**
 * Compute the product of two vectors.
 * Both vector must have the same size. No check is done
 * @param v1     - the first vector (an array of si-t);
 * @param v2     - the second vector (an array of si-t);
 * @param result - is a pointer on a signed int number allocated by the caller
 * (and initialized by the callee).
 * @return 0 if successful.
 */
int ui_vector_scalar_product(ui_vector* v1, ui_vector* v2, unsigned int* result); 

#endif /* ui_vector_h */
