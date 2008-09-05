/* mpfr_vector.h */
/******************************************************************************/

/* This version is much improved over those used in previous projects 
   ST, 2006-07-31 */

#ifndef MPFR_VECTOR_h
#define MPFR_VECTOR_h

/* user must include: mpfr.h before this file. */

#include "mpfr-vector-types.h"

/**
 * Create a new mpfr_vector, allocating memory, initializing all the
 * elements and setting them to 0.
 * @param NbElems   - the number of elements;
 * @param precision - the precision of the created MPFR elements (in
 *                    number of bits);
 * @return a pointer on the vector if all parameters have acceptable values 
 * and if all allocation and initialization operations succeed, 
 * returns NULL otherwise.
 */
mpfr_vector* 
mpfr_vector_alloc(unsigned NbElems,
                  unsigned precision);

/**
 * Create a new mpfr_vector from an array of mpfr_t elements.
 * The contents the vector is a copy of the array.
 * @param mpfr_array - pointer on the array;
 * @param NbElems    - the number of elements in the array;
 * @param precision  - the precision, in bits, of the vector elements.
 * @return a pointer on the newly allocated and initialized array.
 */ 
mpfr_vector* 
mpfr_vector_from_array(mpfr_t* mpfr_array, 
                        unsigned int NbElems,
                        unsigned int precision);

void 
mpfr_vector_free(mpfr_vector *Vect);

void 
mpfr_vector_print(FILE *Dst, 
                  int radix, 
                  size_t digits, 
                  mp_rnd_t rounding_mode,
                  mpfr_vector *Vect);

void 
mpfr_vector_print_no_dims(FILE *Dst, 
                          int radix, 
                          size_t digits, 
                          mp_rnd_t rounding_mode,
                          mpfr_vector *Vect);

void 
mpfr_vector_print_matlab(FILE *Dst, 
                          int radix, 
                          size_t digits, 
                          mp_rnd_t rounding_mode,
                          mpfr_vector *Vect);


mpfr_vector*
mpfr_vector_read_from_stdin(int radix, 
                            mp_rnd_t rounding_mode, 
                            unsigned int precision);

void 
mpfr_vector_read_input_from_file(mpfr_vector *Vect, 
                                  int radix, 
                                  mp_rnd_t rounding_mode,
                                  unsigned int precision,
                                  FILE *Scr);

mpfr_vector*
mpfr_vector_read_from_file(int radix, 
                            mp_rnd_t rounding_mode,
                            unsigned int precision,
                            FILE *Scr);

mpfr_t* 
mpfr_vector_get_at(mpfr_vector* vect, 
                    unsigned int pos);

void 
mpfr_vector_set_at(mpfr_vector* vect, 
                    unsigned int pos, 
                    mpfr_t value,
                    mp_rnd_t rounding_mode);

/**
 * Check whether there are duplicate values in a vector
 * 
 * @param vect - the vector to check;
 * @return 1 if there are duplicates, 0 otherwise
 */
int
mpfr_vector_has_duplicates(mpfr_vector* vect);

unsigned int
mpfr_vector_get_size(mpfr_vector* vect);

/**
 * Compute a new vector whose components are all positive
 * (as the abolute values of the input vector).
 */
mpfr_vector*
mpfr_vector_get_abs(mpfr_vector* vect,
                    unsigned int precision);


/**
 * Compute the scalar product of two vectors.
 * Both vector must have the same size. No check is done
 * @param v1        - the first vector;
 * @param v2        - the second vector;
 * @param result    - a pointer on a mpfr_t number allocated by the caller
 * (and assigned by the callee); 
 * @param precision - the precision for computations in MPFR.
 * @return 0 if OK a different value if something goes wrong.
 */
int 
mpfr_vector_scalar_product(mpfr_vector* v1, 
                           mpfr_vector* v2, 
                           mpfr_t* result,
                           unsigned int precision);

/**
 * Compute the (scalar) sum of all the elements of a vector.
 * @param vect  - the vector;
 * @param sum - a pointer on a mpfr_t number allocated by the caller
 * (and assigned by the callee); 
 * @param precision - the precision for computations in MPFR;
 * @return 0 if OK a different value if something goes wrong.
 */
int 
mpfr_vector_elems_sum(mpfr_vector* vect, 
                           mpfr_t* sum,
                           unsigned int precision);

/**
 * Set each element of a vector to the maximum value of
 * the corresponding elements from two vectors.
 *
 * Here the elements of the first vector are set to the
 * maximum value of the corresponding elements from
 * itself and another vector.
 * 
 * Ex: v1 = (1, 3, 5, 7) v2 = (2, 2, 6, 6) ->
 * v1 = (2, 3, 6, 7).
 *
 * @param in_out_vect - the first vector to pick
 *                      the maximum from and the
 *                      vector that will hold it.
 * @param other_vect  - the other vector to pick
 *                      the maximum from.
 * @return in_out_vector (if successful) or NULL 
 * (if something goes wrong).
 *  
 */
mpfr_vector*
mpfr_vector_set_to_min (mpfr_vector* in_out_vect, 
                        mpfr_vector* other_vect);
/**
 * Set each element of a vector to the minimum value of
 * the corresponding elements from two vectors.
 *
 * Here the elements of the first vector are set to the
 * minimum value of the corresponding elements from
 * itself and another vector.
 * 
 * Ex: v1 = (1, 3, 5, 7) v2 = (2, 2, 6, 6) ->
 * v1 = (1, 2, 5, 6).
 *
 * @param in_out_vect - the first vector to pick
 *                      the minimum from and the
 *                      vector that will hold it.
 * @param other_vect  - the other vector to pick
 *                      the minimum from.
 * @return in_out_vector (if successful) or NULL 
 * (if something goes wrong).
 *  
 */
mpfr_vector*
mpfr_vector_set_to_max (mpfr_vector* in_out_vect, 
                        mpfr_vector* other_vect);
#endif /* MPFR_VECTOR_h */
