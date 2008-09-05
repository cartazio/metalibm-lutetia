/* mqz_matrix.h */

#ifndef MQZ_MATRIX_h
#define MQZ_MATRIX_h

/* user must include: gmp.h before this file. */

#include "mpq-matrix-types.h"
#include "../../vectors/mpq/mpq-vector.h"

/**
 * Create a new rational matrix, allocating memory and initializing all the
 * elements.
 * @param NbRows, the number of rows;
 * @param NbColumns, the number of columns;
 * @return a pointer on the matrix if all parameters have acceptable values and all
 * allocation and initialization operations succeed, NULL otherwise.
 */
mpq_matrix* mpq_matrix_alloc(unsigned NbRows,unsigned NbColumns);

void mpq_matrix_free(mpq_matrix *Mat);

void mpq_matrix_print(FILE *Dst, int radix, mpq_matrix *Mat);

void mpq_matrix_print_no_dims(FILE *Dst, int radix, mpq_matrix *Mat);

void 
mpq_matrix_print_matlab(FILE *Dst, int radix, mpq_matrix *Mat);

/* 
 * Load the contents of the 'Src' FILE into
 * the 'Mat' mpq_matrix.
 * 
 * The 'Mat' mpq_matrix has already been allocated.
 */
void 
mpq_matrix_read_input_from_file(mpq_matrix *Mat, 
                                int radix, 
                                FILE* Src);
/* 
 * Create a mpq_matrix with the data read from stdin. 
 *
 */
mpq_matrix*
mpq_matrix_read_from_stdin(int radix);

/* 
 * Create a mpq_matrix with the date read from a file. 
 *
 */
mpq_matrix*
mpq_matrix_read_from_file(int radix, 
                          FILE* Src);

mpq_t* 
mpq_matrix_get_at(mpq_matrix* mat, 
                  unsigned int row, 
                  unsigned int column);

void 
mpq_matrix_set_at(mpq_matrix* mat, 
                  unsigned int row, 
                  unsigned int column, 
                  mpq_t value);

unsigned int
mpq_matrix_get_num_columns(mpq_matrix* mat);

unsigned int
mpq_matrix_get_num_rows(mpq_matrix* mat);

/**
 * Set a row of a matrix with the elements of a vector.
 *
 * The row and the vector must hold the same number of elements.
 *
 */
int
mpq_matrix_set_row_from_vector(mpq_matrix* mat, 
                                unsigned int row, 
                                mpq_vector* vect); 

/**
 * Get (allocate and initialize) a vector from a matrix row.
 *
 */
mpq_vector*
mpq_matrix_get_vector_from_row(mpq_matrix* mat, 
                                unsigned int row);

/**
 * Inverse a square matrix with the Gauss-Jordan elimination with full
 * pivoting method.
 *
 * The routine is an adaptation to mpq_matrix of the "gaussj" function
 * described in "Numerical recipes in C 2nd edition" pages 39-40. 
 *
 * This is a simplified version that only inverses the matrix. 
 * It does not compute the solutions.
 *
 * @param a - the square matrix to inverse; the inversion is done "in
 *            place": the original content is lost.
 * @return - 0 if everything goes well, 1 if something wrong happens;
 *           the most obvious reason for a failure is if the matrix is
 *           is singular; other causes are invalid parameters and internal
 *           software errors. 
 */
int mpq_matrix_inv_gauss_jordan(mpq_matrix* a);

/**
 * Swap to matrix elements. 
 *
 * @param mat   - the matrix;
 * @param row_1 - the row index of the first element
 * @param col_1 - the column index of the first element
 * @param row_2 - the row index of the second element
 * @param col_2 - the column index of the second element
 */
void mpq_matrix_swap_elements(mpq_matrix* mat,
                              unsigned int row_1,
                              unsigned int col_1,
                              unsigned int row_2,
                              unsigned int col_2);
                      

mpq_matrix*
mpq_matrix_naive_prod(mpq_matrix* mat1,
                        mpq_matrix* mat2);

mpq_matrix*
mpq_matrix_clone(mpq_matrix* orig_matrix);

#endif /* MPQ_MATRIX_h */
