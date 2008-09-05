/* si_vector.h */

#ifndef SI_VECTOR_h
#define SI_VECTOR_h

#define SI_VECTOR_UP 0
#define SI_VECTOR_DOWN 1

#include "si-vector-types.h"

/**
 * Create a new si vector, allocating memory, initializing all the
 * elements and setting them to 0.
 * @param NbElems, the number of elements;
 * @return a pointer on the vector if all parameters have acceptable values 
 * and if all allocation and initialization operations succeed, 
 * returns NULL otherwise.
 */
si_vector* si_vector_alloc(unsigned NbElems);

si_vector* si_vector_from_array(signed int* si_array, unsigned int NbElems);

void si_vector_free(si_vector *Vect);

void si_vector_print(FILE *Dst, 
                        si_vector *Vect);

void si_vector_print_no_dims(FILE *Dst, 
                        si_vector *Vect);

void si_vector_print_matlab(FILE *Dst, 
                        si_vector *Vect);


si_vector *si_vector_read_from_stdin(int radix);

void si_vector_read_input_from_file(si_vector *Vect, 
                                      int radix, 
                                      FILE *Scr);

si_vector* si_vector_read_from_file(int radix, 
                                    FILE *Scr);

signed int* si_vector_get_at(si_vector* vect, 
                              unsigned int pos);

void si_vector_set_at(si_vector* vect, 
                        unsigned int pos, 
                        signed int value);

unsigned int si_vector_get_NbElems(si_vector* vect);

unsigned int si_vector_get_size(si_vector* vect);

/**
 * Check whether there are duplicate values in a vector
 * 
 * @param vect - the vector to check;
 * @return 1 if there are duplicates, 0 otherwise
 */
int
si_vector_has_duplicates(si_vector* vect);

/**
 * Compute a new vector whose components are all positive
 * (as the abolute values of the input vector).
 */
si_vector*
si_vector_get_abs(si_vector* vect);

/**
 * Shift an element of the vector of one position, overwriting
 * the neigbouring element.
 *
 * No return value.
 *
 * @param vect      - the vector;
 * @param index     - the index of the element to shift;
 * @param direction - the direction of the shift: SI_VECTOR_UP -> to an higher index
 *                                                SI_VECTOR_DOWN -> to a lower index.
 */
void
si_vector_shift(si_vector* vect, 
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
 * @param - direction, the sorting direction: SI_VECTORE_UP ->  increasing
 *                                            SI_VECTOR_DOWN -> decreasing.
 */
void
si_vector_insertion_sort(si_vector* vi, int direction);

/**
 * Check wether two vectors are equal.
 *
 * Equal vectors have:
 * - the same size (or be both NULL pointers);
 * - the same values at the same indices.
 * @param vect_1 - the first vector pointer;
 * @param vect_2 - the second vector pointer;
 * @return 1 if vect_1 and vect_2 are equal;
 *         0 if they differ.
 */
int 
si_vector_are_equal(si_vector* vect_1, si_vector* vect_2);

/**
 * Clone a vector.
 *
 *
 */
si_vector*
si_vector_clone(si_vector* orig_vect);

/**
 * Compute the product of two vectors.
 * Both vector must have the same size. No check is done
 * @v1 is the first vector (an array of si-t);
 * @v2 is the second vector (an array of si-t);
 * @result is a pointer on a signed int number allocated by the caller
 * (and initialized by the callee).
 * @size is the number of elements of both vectors
 */
int si_vector_vector_product(si_vector* v1, si_vector* v2, signed int* result); 

#endif /* SI_VECTOR_h */
