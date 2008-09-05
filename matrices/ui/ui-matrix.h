/* ui_matrix.h */

#ifndef UI_MATRIX_h
#define UI_MATRIX_h

#include "ui-matrix-types.h"

/**
 * Create a new integer matrix, allocating memory and initializing all the
 * elements.
 * @param NbRows    - the number of rows;
 * @param NbColumns - the number of columns;
 * @return a pointer on the matrix if all parameters have acceptable values and all
 * allocation and initialization operations succeed, NULL otherwise.
 */
ui_matrix* 
ui_matrix_alloc(unsigned NbRows,unsigned NbColumns);

void 
ui_matrix_free(ui_matrix *Mat);

void 
ui_matrix_print(FILE *Dst, int radix, ui_matrix *Mat);

void 
ui_matrix_print_no_dims(FILE *Dst, int radix, ui_matrix *Mat);

void 
ui_matrix_print_matlab(FILE *Dst, int radix, ui_matrix *Mat);

void 
ui_matrix_read_input_from_file(ui_matrix *Mat, int radix, FILE *Src);

ui_matrix*
ui_matrix_read_from_stdin(int radix, FILE *Src);

ui_matrix*
ui_matrix_read_from_file(int radix, FILE *Src);

unsigned int* 
ui_matrix_get_at(ui_matrix* mat, 
                  unsigned int row, 
                  unsigned int column);

void
ui_matrix_set_at(ui_matrix* mat, 
                  unsigned int row, 
                  unsigned int column, 
                  unsigned int value);

unsigned int
ui_matrix_get_num_colums(ui_matrix* mat);

unsigned int
ui_matrix_get_num_rows(ui_matrix* mat);

int
ui_matrix_set_row_from_vector(ui_matrix* mat, 
                              unsigned int row, 
                              ui_vector* vect); 

ui_vector*
ui_matrix_get_vector_from_row(ui_matrix* mat, unsigned int row);

#endif /* UI_MATRIX_h */
