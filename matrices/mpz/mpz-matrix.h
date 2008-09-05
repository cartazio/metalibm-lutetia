/* mpz_matrix.h */

#ifndef MPZ_MATRIX_h
#define MPZ_MATRIX_h

/* user must include: gmp.h before this file. */

#include "mpz-matrix-types.h"

/**
 * Create a new integer matrix, allocating memory and initializing all the
 * elements.
 * @param NbRows    - the number of rows;
 * @param NbColumns - the number of columns;
 * @return a pointer on the matrix if all parameters have acceptable values and all
 * allocation and initialization operations succeed, NULL otherwise.
 */
mpz_matrix* mpz_matrix_alloc(unsigned NbRows,unsigned NbColumns);

void mpz_matrix_free(mpz_matrix *Mat);

void mpz_matrix_print(FILE *Dst, int radix, mpz_matrix *Mat);

void mpz_matrix_print_no_dims(FILE *Dst, int radix, mpz_matrix *Mat);

void mpz_matrix_print_matlab(FILE *Dst, int radix, mpz_matrix *Mat);

void mpz_matrix_read_input(mpz_matrix *Mat, int radix);

mpz_matrix *mpz_matrix_read(int radix);

mpz_t* mpz_matrix_get_at(mpz_matrix* mat, 
                          unsigned int row, 
                          unsigned int column);

void mpz_matrix_set_at(mpz_matrix* mat, 
                        unsigned int row, 
                        unsigned int column, 
                        mpz_t value);

#endif /* MPZ_MATRIX_h */
