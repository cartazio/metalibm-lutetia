/* mpq_vector.h */

#ifndef MPQ_VECTOR_h
#define MPQ_VECTOR_h

/* user must include: gmp.h before this file. */

#include "mpq-vector-types.h"

/**
 * Create a new float vector, allocating memory, initializing all the
 * elements and setting them to 0.
 * @param NbElems, the number of elements;
 * @return a pointer on the vector if all parameters have acceptable values 
 * and if all allocation and initialization operations succeed, 
 * returns NULL otherwise.
 */
mpq_vector* mpq_vector_alloc(unsigned NbElems);

mpq_vector* mpq_vector_from_array(mpq_t* mpq_array, unsigned int NbElems);

void mpq_vector_free(mpq_vector *Vect);

void mpq_vector_print(FILE *Dst, 
                        int radix, 
                        mpq_vector *Vect);

void mpq_vector_print_no_dims(FILE *Dst, 
                        int radix, 
                        mpq_vector *Vect);

void mpq_vector_print_matlab(FILE *Dst, 
                        int radix, 
                        mpq_vector *Vect);


mpq_vector* mpq_vector_read_from_stdin(int radix);

void mpq_vector_read_input_from_file(mpq_vector *Vect, 
                                        int radix, 
                                        FILE *Scr);

mpq_vector* mpq_vector_read_from_file(int radix, 
                                          FILE *Scr);

mpq_t* mpq_vector_get_at(mpq_vector* vect, 
                          unsigned int pos);

void mpq_vector_set_at(mpq_vector* vect, 
                        unsigned int pos, 
                        mpq_t value);

unsigned int 
mpq_vector_get_NbElems(mpq_vector* vect);

unsigned int
mpq_vector_get_size(mpq_vector* vect);

/**
 * Check whether the same value occurs twice (or more) in a vector.
 *
 * @param vect - a pointer on the vector;
 * @return 0 if there are no duplicates
 *         1 otherwise.
 */
int
mpq_vector_has_duplicates(mpq_vector* vect);

/**
 * Compute a vector whose components are all positive
 * (as the abolute values of the input vector).
 */
mpq_vector*
mpq_vector_get_abs(mpq_vector* vect);

/**
 * Order the elements of the vector.
 *
 * Sorting is made in-place (the original vector is altered) and
 * there is no return value.
 * using the insertion sort alogorithm.
 *
 * @param - vect, the vector to order
 * @param - direction, the sorting direction: if direction == 0, increasing
 *                                            else, decreasing.
 */
void
mpq_vector_insertion_sort(mpq_vector* vect, int direction);

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
mpq_vector_are_equal(mpq_vector* vect_1, mpq_vector* vect_2);

/**
 * Compute the product of two vectors.
 * Both vector must have the same size. No check is done
 * @v1 is the first vector (an array of mpq_t);
 * @v2 is the second vector (an array of mpq_t);
 * @result is a pointer on a mpq_t number allocated by the caller
 * (and initialized by the callee).
 * @size is the number of elements of both vectors
 */
int mpq_vector_vector_product(mpq_vector* v1, 
                                mpq_vector* v2, 
                                mpq_t* result);


#endif /* MPQ_VECTOR_h */
