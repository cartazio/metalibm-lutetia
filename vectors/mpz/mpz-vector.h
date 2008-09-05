/* mpz_vector.h */

#ifndef MPZ_VECTOR_h
#define MPZ_VECTOR_h

/* user must include: gmp.h before this file. */

#include "mpz-vector-types.h"

/**
 * Create a new float vector, allocating memory, initializing all the
 * elements and setting them to 0.
 * @param NbElems - the number of elements;
 * @return a pointer on the vector if all parameters have acceptable values 
 * and if all allocation and initialization operations succeed, 
 * returns NULL otherwise.
 */
mpz_vector* mpz_vector_alloc(unsigned NbElems);

mpz_vector* mpz_vector_from_array(mpz_t* mpz_array, unsigned int NbElems);

void mpz_vector_free(mpz_vector *Vect);

void mpz_vector_print(FILE *Dst, 
                        int radix, 
                        mpz_vector *Vect);

void mpz_vector_print_no_dims(FILE *Dst, 
                        int radix, 
                        mpz_vector *Vect);

void mpz_vector_print_matlab(FILE *Dst, 
                        int radix, 
                        mpz_vector *Vect);


mpz_vector* mpz_vector_read_from_stdin(int radix);

void mpz_vector_read_input_from_file(mpz_vector *Vect, 
                                        int radix, 
                                        FILE *Scr);

mpz_vector* mpz_vector_read_from_file(int radix, 
                                          FILE *Scr);

mpz_t* mpz_vector_get_at(mpz_vector* vect, 
                          unsigned int pos);

void mpz_vector_set_at(mpz_vector* vect, 
                        unsigned int pos, 
                        mpz_t value);

/**
 * Compute the product of two vectors.
 * Both vector must have the same size. No check is done
 * @param v1     - is the first vector (an array of mpz_t);
 * @param v2     - is the second vector (an array of mpz_t);
 * @param result - a pointer on a mpz_t number allocated by the caller
 * (and initialized by the callee).
 */
int mpz_vector_vector_product(mpz_vector* v1, 
                                mpz_vector* v2, 
                                mpz_t* result);

unsigned int mpz_vector_get_NbElems(mpz_vector* vect);

unsigned int mpz_vector_get_size(mpz_vector* vect);

#endif /* MPZ_VECTOR_h */
