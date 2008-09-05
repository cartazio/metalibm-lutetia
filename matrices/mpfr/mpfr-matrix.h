/* mpfr_matrix.h */

#ifndef MPFR_MATRIX_h
#define MPFR_MATRIX_h

/* user must include: gmp.h 
                      mpfr.h
   before this file. */

#include "mpfr-matrix-types.h"
#include "../../vectors/mpfr/mpfr-vector.h"

/**
 * Create a new float matrix, allocating memory, initializing all the
 * elements and setting them to 0.
 * @param NbRows    - the number of rows;
 * @param NbColumns - the number of columns;
 * @param prec      - the precision (in bits);
 * @return a pointer on the matrix if all parameters have acceptable values and all
 * allocation and initialization operations succeed, NULL otherwise.
 */
mpfr_matrix* 
mpfr_matrix_alloc(unsigned NbRows,
                  unsigned NbColumns,
                  mp_prec_t prec);

void 
mpfr_matrix_free(mpfr_matrix *Mat);

void 
mpfr_matrix_print(FILE *Dst, 
                  int radix, 
                  size_t digits, 
                  mp_rnd_t rounding_mode,
                  mpfr_matrix *Mat);

void
mpfr_matrix_print_no_dims(FILE *Dst, 
                          int radix, 
                          size_t digits, 
                          mp_rnd_t rounding_mode,
                          mpfr_matrix *Mat);

void 
mpfr_matrix_print_matlab(FILE *Dst, 
                          int radix, 
                          size_t digits, 
                          mp_rnd_t rounding_mode,
                          mpfr_matrix *Mat);


mpfr_matrix*
mpfr_matrix_read_from_stdin(int radix, 
                              mp_rnd_t rounding_mode, 
                              unsigned int precision);

/* 
 * Load the contents of the FILE Src into the existing matrix Mat. 
 *
 * Matrix elements are expressed as strings convertible into mpfr_t
 * or into a mpq_t. In the latter case, the rational numbers are
 * converted into floats wiht the given precision.
 *
 * The file is a text file where each line holds a row.
 *
 * @param Mat - a pointer to the matrix to fill (the matrix must have
 *              been created and allocate beforehand);
 * @param radix - the radix used to write the number in the file;
 * @param rounding_mode - the rounding mode for approximate conversion;
 * @param precision - the precision (in bits) for approximate conversion;
 * @param Src  the source file that holds the matix elements.
 * @return - 0 if all data where correctly read;
 *         - 1 should a problem arise.
 *
 * TODO: switch all the read_input_from _file to this model.
 */
int
mpfr_matrix_read_input_from_file(mpfr_matrix *Mat, 
                                        int radix, 
                                        mp_rnd_t rounding_mode,
                                        unsigned int precision,
                                        FILE *Src);

mpfr_matrix* 
mpfr_matrix_read_from_file(int radix, 
                           mp_rnd_t rounding_mode,
                           unsigned int precision,
                           FILE *Scr);

mpfr_t* 
mpfr_matrix_get_at(mpfr_matrix* mat, 
                    unsigned int row, 
                    unsigned int column);

void 
mpfr_matrix_set_at(mpfr_matrix* mat, 
                    unsigned int row, 
                    unsigned int column, 
                    mpfr_t value,
                    mp_rnd_t rounding_mode);

unsigned int
mpfr_matrix_get_num_columns(mpfr_matrix* mat);

unsigned int
mpfr_matrix_get_num_rows(mpfr_matrix* mat);

int
mpfr_matrix_set_row_from_vector(mpfr_matrix* mat, 
                                unsigned int row, 
                                mpfr_vector* vect); 

mpfr_vector*
mpfr_matrix_get_vector_from_row(mpfr_matrix* mat, 
                                unsigned int row,
                                unsigned int precison);
/**
 * Is a matrix an identity matrix for some givent precision.
 *
 * If all the elements of the matrix except the elements of
 * the diagonal are equal to 0 and the elements of the diagonal
 * are equal to 1, then the matrix is an identity matrix.
 * Equality to 0 or to 1 is evaluated for a given precision
 * @param mat - the matrix ;
 * @param precision - the precision (in number of bits)
 * @return - 1 if the matrix is an identity matrix, 0 otherwise
 *         (even in the case of an error).
 */     
int
mpfr_matrix_is_id(mpfr_matrix* mat,
                  unsigned int precision);

/**
 * The maximum absolute value found in a matrix diagonal excepted.
 * 
 * Will find the maximum absolute value in the matrix without taking
 * into account the elements of the diagonal.
 * 
 * @param mat - the matrix
 * @param max_value - the mpfr_t value holding the result (must be
 *                    initialized by the caller)
 * @return - 1 if there is an error.
 */
 int
 mpfr_matrix_max_abs_no_diag(mpfr_matrix* mat,
                              mpfr_t max_value);
 
/**
 * Inverse a square matrix with the Gauss-Jordan elimination with full
 * pivoting method.
 *
 * The routine is an adaptation to mpfr_matrix of the "gaussj" function
 * described in "Numerical recipes in C 2nd edition" pages 39-40. 
 *
 * This is a simplified version that only inverses the matrix. 
 * It does not compute the solutions.
 *
 * @param mat - the square matrix to inverse; the inversion is done "in
 *              place": the original content is lost;
 * @param zero_threshold - due to rounding errors, testing singularity
 *                         against a zero value is not very robust;
 *                         instead, we use a user defined "problem dependent
 *                         convenient small value";
 * @return - 0 if everything goes well, 1 if something wrong happens;
 *           the most obvious reason for a failure is if the matrix is
 *           is singular; other causes are invalid parameters and internal
 *           software errors. 
 */
int mpfr_matrix_inv_gauss_jordan(mpfr_matrix* mat, mpfr_t zero_threshold);

/**
 * Swap to matrix elements. 
 *
 * @param mat   - the matrix;
 * @param row_1 - the row index of the first element
 * @param col_1 - the column index of the first element
 * @param row_2 - the row index of the second element
 * @param col_2 - the column index of the second element
 */
void mpfr_matrix_swap_elements(mpfr_matrix* mat,
                                unsigned int row_1,
                                unsigned int col_1,
                                unsigned int row_2,
                                unsigned int col_2);
                      
/**
 * A naive implementation of the matrix product.
 *
 * @param mat1 - the first matrix (m,n);
 * @param mat2 - the second matrix (n,m);
 * @param prec - the precision, in bits, of the
 *               elements of the product matrix.
 * @return a pointer to the product matrix (or
 * NULL if something went wrong).
 */

mpfr_matrix*
mpfr_matrix_naive_prod(mpfr_matrix* mat1, 
                        mpfr_matrix* mat2,
                        mp_prec_t prec);

mpfr_matrix*
mpfr_matrix_clone(mpfr_matrix* orig_matrix, 
                  mp_prec_t prec);
#endif /* MPFR_MATRIX_h */
