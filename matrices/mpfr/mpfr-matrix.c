/* mpfr_matrix.c */
/******************************************************************************/

/** @file
 * Implementation of mpfr_matrix behaviour.
 *
 * @author ST 
 * @date 2005-07-21
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <gmp.h>
#include <mpfr.h>

#include "mpfr-matrix.h"
#include "../../misc-utils/std-exit-errors.h"

void
mpfr_matrix_errormsg(char *f , char *msgname, char *msg);

/* 
 * Allocate space for matrix dimensioned by 'NbRows X NbColumns'.
 */
mpfr_matrix* mpfr_matrix_alloc(unsigned NbRows,
                                unsigned NbColumns,
                                mp_prec_t prec) 
{
  mpfr_matrix *Mat;
  mpfr_t *p, **q;
  int i,j;

  mpfr_set_default_prec(prec);
  /* Deal with pathological parameter */
  if (NbRows < 0)
  {
    mpfr_matrix_errormsg("mpfr_matrix_alloc", "invalid parameter", "NbRows must be >= 0");
    return 0;
  }
  if (NbRows < 0)
  {
    mpfr_matrix_errormsg("mpfr_matrix_alloc", "invalid parameter", "NbColumns must be >= 0");
    return 0;
  }
  /* Create the matrix struct */
  Mat=(mpfr_matrix *)malloc(sizeof(mpfr_matrix));
  if(!Mat) {	
    mpfr_matrix_errormsg("mpfr_matrix_Alloc", "outofmem", "out of memory space");
    return 0;
  }
  Mat->NbRows=NbRows;
  Mat->NbColumns=NbColumns;
  /* NbRows == 0: created an "empty" matrix. */
  if(NbRows==0) 
  {
    Mat->p = (mpfr_t **)0;
    Mat->p_Init= (mpfr_t *)0;
  }  
  else 
  {
    /* NbColumns == 0: created an "empty" matrix. */
    if(NbColumns==0) 
    {
      Mat->p = (mpfr_t **)0;
      Mat->p_Init= (mpfr_t *)0;
    }
    else 
    {
      /* Create the elements array that holds the matrix elements */
      q = (mpfr_t **)malloc(NbRows * sizeof(*q));
      if(!q) 
      {
        free(Mat);
        mpfr_matrix_errormsg("mpfr_matrix_Alloc", "outofmem", "out of memory space");
        return 0;
      }
      p = (mpfr_t *)malloc(NbRows * NbColumns * sizeof(mpfr_t));
      if(!p) 
      {
        free(q);
        free(Mat);
        mpfr_matrix_errormsg("mpfr_matrix_Alloc", "outofmem", "out of memory space");
        return 0;
      }
      Mat->p = q;
      Mat->p_Init = p;
      for (i=0;i<NbRows;i++) 
      {
        *q++ = p;
        for (j=0;j<NbColumns;j++)
          mpfr_init_set_ui(*(p+j),0,GMP_RNDN);
        p += NbColumns;
      }
    }
  }
  p = NULL;
  q = NULL;
  Mat->p_Init_size = NbColumns*NbRows;
  /*fprintf(stderr, 
          "mpfr_matrix_alloc: current precision = %u\n", 
          mpfr_get_default_prec()); */
  return Mat;
} /* mpfr_matrix_alloc */

/* 
 * Free the memory space occupied by mpfr_matrix 'Mat' 
 */
void mpfr_matrix_free(mpfr_matrix *Mat) {
  
  int i;
  mpfr_t *p;

  if(Mat == NULL) return;

  if( Mat->p )
  {
    p = *(Mat->p);
    for(i=0 ; i < Mat->p_Init_size ; i++) {
      mpfr_clear(*p++);
    }
  }

  if (Mat->p_Init)
    free(Mat->p_Init);
  if (Mat->p)
    free(Mat->p);
  free(Mat);

} /* mpfr_matrix_Free */

/* 
 * Print the contents of the mpfr_matrix 'Mat'.
 */
void mpfr_matrix_print(FILE *Dst, 
                        int radix, 
                        size_t digits, 
                        mp_rnd_t rounding_mode,
                        mpfr_matrix *Mat) 
{
  
  mpfr_t *p;
  int i, j;
  unsigned NbRows, NbColumns;
  
  fprintf(Dst,"%d %d\n", NbRows=Mat->NbRows, NbColumns=Mat->NbColumns);
  for (i=0 ; i < NbRows ; i++) 
  {
    p=*(Mat->p+i);
    for (j=0;j<NbColumns;j++) 
    {
      mpfr_out_str(Dst, radix, digits, *p++, rounding_mode);
      fprintf(Dst, " ");
    }
    fprintf(Dst, "\n");
  }
} /* mpfr_matrix_print */

/* 
 * Print the contents of the mpfr_matrix 'Mat' without printing the
 * the dimensions.
 */
void mpfr_matrix_print_no_dims(FILE *Dst, 
                                int radix, 
                                size_t digits, 
                                mp_rnd_t rounding_mode,
                                mpfr_matrix *Mat) 
{
  
  mpfr_t *p;
  int i, j;
  unsigned NbRows, NbColumns;
  
  NbRows    = Mat->NbRows; 
  NbColumns = Mat->NbColumns;
  for (i=0 ; i<NbRows ; i++) 
  {
    p=*(Mat->p+i);
    for (j=0;j<NbColumns;j++) 
    {
      mpfr_out_str(Dst, radix, digits, *p++, rounding_mode);
      fprintf(Dst, " ");
    }
    fprintf(Dst, "\n");
  }
} /* mpfr_matrix_print_no_dims */

/* 
 * Print the contents of the mpfr_matrix 'Mat' in the Matlab
 * input format.
 */
void mpfr_matrix_print_matlab(FILE *Dst, 
                                int radix, 
                                size_t digits, 
                                mp_rnd_t rounding_mode,
                                mpfr_matrix *Mat) 
{
  
  mpfr_t *p;
  int i, j;
  unsigned NbRows, NbColumns;
  
  NbRows    = Mat->NbRows; 
  NbColumns = Mat->NbColumns;
  fprintf(Dst, "[ ");
  for (i=0 ; i < NbRows ; i++) 
  {
    p=*(Mat->p+i);
    for (j=0 ; j < NbColumns ; j++) 
    {
      mpfr_out_str(Dst, radix, digits, *p++, rounding_mode);
      if (j < (NbColumns - 1))
      {
        fprintf(Dst, ", ");
      }
      else
      {
        if (i < (NbRows - 1))
        {
          fprintf(Dst, "; ");
        }
        else
        {
          fprintf(Dst, " ]\n");
        }
      }
    }
  }
} /* mpfr_matrix_print_matlab */

/*
 * @see mpfr-matrix.h#mpfr_matrix_read_input_from_file
 */
int
mpfr_matrix_read_input_from_file(mpfr_matrix *Mat, 
                                  int radix, 
                                  mp_rnd_t rounding_mode,
                                  unsigned int precision,
                                  FILE* Src) 
{
  
  mpfr_t *p;
  mpq_t   rat_temp;
  int i,j,n;
  char *c, s[1024],str[1024];

  mpfr_set_default_prec(precision);

  p = Mat->p_Init;
  for (i=0 ; i < Mat->NbRows ; i++) 
  {
    do 
    {
      c = fgets(s, 1024, Src);
      while(isspace(*c) && *c != '\n')
	           ++c;
    } 
    while(c && (*c == '#' || *c == '\n'));
    
    if (!c) 
    {
      mpfr_matrix_errormsg( "mpfr_matrix_readinput_from_file", 
                            "baddim", 
                            "not enough rows" );
      return(1);
    }
    for (j=0 ; j < Mat->NbColumns ; j++) 
    {
      if(!c || *c=='\n' || *c=='#') 
      {
        mpfr_matrix_errormsg("mpfr_matrix_read_input_from_file", 
                              "baddim", 
                              "not enough columns");
        return(1);
      }
      if (sscanf(c,"%s%n",str,&n) == 0) 
      {
        mpfr_matrix_errormsg( "mpfr_matrix_Read", 
                              "baddim", 
                              "not enough columns" );
        return(1);
      }
      /* Try to read the float directly: mpfr_set_str returns 0 if OK. */
      if (mpfr_set_str(*p, str, radix, rounding_mode))
      {
        /* Try to read it as a rational: mpq_set_str returns 0 if OK. */
        mpq_init(rat_temp);
        if (mpq_set_str(rat_temp, str, radix))
        {
          mpfr_matrix_errormsg("mpfr_vector_read_input_from_file", 
                               "badformat", 
                               "number in a bad format");
          mpq_clear(rat_temp);
          return(1);
        }
        mpfr_set_q(*p, rat_temp, GMP_RNDN);
        mpq_clear(rat_temp);
      }
      p++;
      c += n;
    } /* End for j */
  } /* End for i */
  return(0);
} /* mpfr_matrix_read_input_from_file */

/* 
 * Create a new matrix filled with the data read from standard input. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the matrix contents, one row per line.
 */
mpfr_matrix *mpfr_matrix_read_from_stdin(int radix, 
                                          mp_rnd_t rounding_mode, 
                                          unsigned int precision) 
{
  
  mpfr_matrix *Mat;
  unsigned NbRows, NbColumns;
  char s[1024];

  mpfr_set_default_prec(precision);
  
  while(fgets(s, 1024, stdin)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d %d", &NbRows, &NbColumns)<2))
    fgets(s, 1024, stdin);
  Mat = mpfr_matrix_alloc(NbRows, NbColumns, precision);
  if(Mat == NULL) 
  {
    mpfr_matrix_errormsg("mpfr_matrix_read", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  if (mpfr_matrix_read_input_from_file(Mat, 
                                       radix,
                                       rounding_mode,
                                       precision, 
                                       stdin))
  {
    mpfr_matrix_free(Mat);
    return(NULL);
  }
  return Mat;
} /* mpfr_matrix_read_from_stdin */

/* 
 * Read the contents of the matrix 'Mat' from FILE pointed to by Src. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the matrix contents, one row per line.
 * Matrix contents is made of strings directly convertible into mpfr_t.
 */
mpfr_matrix*
mpfr_matrix_read_from_file(int radix, 
                            mp_rnd_t rounding_mode,
                            unsigned int precision, 
                            FILE* Src) 
{
  
  mpfr_matrix *Mat;
  unsigned NbRows, NbColumns;
  char s[1024];
  
  mpfr_set_default_prec(precision);
  
  while(fgets(s, 1024, Src)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d %d", &NbRows, &NbColumns)<2))
    fgets(s, 1024, Src);
  Mat = mpfr_matrix_alloc(NbRows, NbColumns, precision);
  if(Mat == NULL) 
  {
    mpfr_matrix_errormsg("mpfr_matrix_read_from_file", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  if (mpfr_matrix_read_input_from_file(Mat, 
                                        radix,
                                        rounding_mode,
                                        precision, 
                                        Src))
  {
    mpfr_matrix_free(Mat);
    return(NULL);
  }
  return Mat;
} /* mpfr_matrix_read_from_file */

mpfr_t*
mpfr_matrix_get_at(mpfr_matrix* mat, unsigned int row, unsigned int column)
{
  mpfr_t** matRow;

  mpfr_set_default_prec(200);
  
  if (mat == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_at", "null pointer parameter", "mat parameter is NULL");
    return(NULL);
  }
  if (mat->p_Init_size == 0)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_at", "empty matrix", "can't return an element from an empty matrix");
    return(NULL);
  }
  if (row >= mat->NbRows)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_at", "index out of range", "row index beyond bounds");
    return(NULL);
  }
  if (column >= mat->NbColumns)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_at", "index out of range", "column index beyond bounds");
    return(NULL);
  }
  matRow = mat->p + row; 
  return(*matRow + column);
} /* End mpfr_matrix_get_at */

/*
 */
void 
mpfr_matrix_set_at(mpfr_matrix* mat, 
                    unsigned int row, 
                    unsigned int column, 
                    mpfr_t value,
                    mp_rnd_t rounding_mode)
{
  mpfr_t** matRow;

  mpfr_set_default_prec(200);

  if (mat == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_set_at", "nullPointerParameter", "mat parameter is NULL");
  }
  if (mat->p_Init_size == 0)
  {
    mpfr_matrix_errormsg("mpfr_matrix_set_at", "empty matrix", "can't set an element from an empty matrix");
  }
  if (row >= mat->NbRows)
  {
    mpfr_matrix_errormsg("mpfr_matrix_set_at", "index out of range", "row index beyond bounds");
  }
  if (column >= mat->NbColumns)
  {
    mpfr_matrix_errormsg("mpfr_matrix_set_at", "index out of range", "column index beyond bounds");
  }
  matRow = mat->p + row; 
  mpfr_set(*(*matRow + column), value, rounding_mode);
} /* End mpfr_matrix_set_at */

unsigned int
mpfr_matrix_get_num_columns(mpfr_matrix* mat)
{
  if (mat == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_num_columns", 
                          "nullPointerParameter", 
                          "mat parameter is NULL");
    return(0);
  }
  return(mat->NbColumns);
} /* End mpfr_matrix_get_num_columns. */

unsigned int
mpfr_matrix_get_num_rows(mpfr_matrix* mat)
{
  if (mat == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_num_rows", 
                          "nullPointerParameter", 
                          "mat parameter is NULL");
    return(0);
  }
  return(mat->NbRows);
} /* End mpfr_matrix_get_num_rows */

/*
 */
int
mpfr_matrix_set_row_from_vector(mpfr_matrix* mat, unsigned int row, mpfr_vector* vect)
{
  unsigned int i;
  if (mat == NULL || vect == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_set_row_from_vector", 
                          "nullPointerParameter", 
                          "one of mat or vect is a NULL pointer");
    return(1);
  }
  if (mpfr_vector_get_size(vect) != mat->NbColumns)
  {
    mpfr_matrix_errormsg("mpfr_matrix_set_row_from_vector", 
                          "incompSize", 
                          "incompatible vector and matrix row sizes");
    return(1);
  }
  if (mat->NbRows <= row)
  {
    mpfr_matrix_errormsg("mpfr_matrix_set_row_from_vector", 
                          "index out of range",
                          "row index beyond bounds");
    return(1);
  }
  for (i = 0 ; i < mat->NbColumns ; i ++)
  {
    mpfr_matrix_set_at(mat, row, i, *(mpfr_vector_get_at(vect, i)), GMP_RNDN);
  }
  return(0);
} /* End mpfr_matrix_set_row_from_vector */

/*
 */
mpfr_vector*
mpfr_matrix_get_vector_from_row(mpfr_matrix* mat, 
                                unsigned int row,
                                unsigned int precision)
{
  mpfr_vector* vect = NULL;
  unsigned int i;
  if (mat == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_vector_from_row", 
                          "nullPointerParameter", 
                          "mat is a NULL pointer");
    return(NULL);
  }
  if (mat->NbRows <= row)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_vector_from_row", 
                          "index out of range",
                          "row index beyond bounds");
    return(NULL);
  }
  vect = mpfr_vector_alloc(mat->NbColumns, precision);
  if (vect == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_get_vector_from_row", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  for (i = 0 ; i < mat->NbColumns ; i++)
  {
    mpfr_vector_set_at(vect, i, *(mpfr_matrix_get_at(mat, row, i)), GMP_RNDN);
  }
  return(vect);
} /* End mpfr_matrix_get_vector_from_row */

/**
  * see mpfr-matrix.c#mpfr_matrix_is_id
  */
int
mpfr_matrix_is_id(mpfr_matrix* mat,
                  unsigned int precision)
{
  unsigned int matrix_order = 0;
  unsigned int i, j;
  mpfr_t tmp;
  mpfr_t current_element;
  mpfr_t error_margin;
  /* Deal with pathological parameters */
  /* A NULL matrix can not be an identity matrix. */
  if (mat == NULL)
  {
    return(0);
  }
  /* An empty matrix can not be an indentity matrix. */
  if ((matrix_order = mpfr_matrix_get_num_columns(mat)) == 0)
  {
    return(0);
  }
  /* A non square matrix can not be an identity matrix. */
  if (matrix_order != mpfr_matrix_get_num_rows(mat))
  {
    return(0);
  }
  /* If precision < 2 then return false, this is a bit arbitrary. */
  if (precision < 2)
  {
    return(0);
  }
  mpfr_init(tmp);
  mpfr_init(current_element);
  mpfr_init(error_margin);
  /* Compute the error margin from the precision. */
  mpfr_set_ui(error_margin, 1, GMP_RNDN);
  mpfr_div_2ui(error_margin, error_margin, precision, GMP_RNDN);
  for (i = 0 ; i < matrix_order ; i++)
  {
    for (j = 0 ; j < matrix_order ; j++)
    {
      mpfr_set(current_element,
                *(mpfr_matrix_get_at(mat, i, j)),
                GMP_RNDN);
      mpfr_abs(current_element, current_element, GMP_RNDN);
      if (i != j) /* current_element must == 0 */
      {
        if (mpfr_cmp(current_element, error_margin) > 0)
        {
          mpfr_clear(tmp);
          mpfr_clear(current_element);
          mpfr_clear(error_margin);
          return(0);
        }
      }
      else /* (i == j) */
      {
        mpfr_sub_ui(tmp, current_element, 1, GMP_RNDN);
        mpfr_abs(tmp, tmp, GMP_RNDN);
        if (mpfr_cmp(tmp, error_margin) > 0)
        {
          mpfr_clear(tmp);
          mpfr_clear(current_element);
          mpfr_clear(error_margin);
          return(0);
        }
      }
    } /* End for j. */
  } /* End for i. */
  mpfr_clear(tmp);
  mpfr_clear(current_element);
  mpfr_clear(error_margin);
  return(1);
} /* End mpfr_matrix_is_id */

/**
 * @see mpfr_matrix.h#mpfr_matrix_max_abs_no_diag
 */
 int
 mpfr_matrix_max_abs_no_diag(mpfr_matrix* mat,
                              mpfr_t max_value)
{
  unsigned int n, i, j;
  mpfr_t abs_value;
  /* Deal with pathological parameters */
  if (mat == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_max_abs_no_diag",
                          "null pointer parameter",
                          "\"mat\" parameter is NULL");
    return(1);
  }
  
  if ((mat->NbRows == 0) || (mat->NbColumns == 0))
  {
    mpfr_matrix_errormsg("mpfr_matrix_max_abs_no_diag",
                          "empty matrix",
                          "\"mat\"a parameter points to an empty matrix");
    return(1);
  }
  if (mat->NbRows != mat->NbColumns)
  {
    mpfr_matrix_errormsg("mpfr_matrix_max_abs_no_diag",
                          "non square matrix",
                          "\"mat\" parameter points to a non square matrix");
    return(1);
  }
  /* Let's do some work */
  n = mat->NbRows;
  mpfr_init(abs_value);
  mpfr_set_ui(max_value, 0, GMP_RNDN);
  for (i = 0 ; i < n ; i++)
  {
    for (j = 0 ; j < n ; j++)
    {
      if (i != j)
      {
        mpfr_abs(abs_value,
                  *(mpfr_matrix_get_at(mat, i, j)),
                  GMP_RNDN);
        if (mpfr_cmp(abs_value, max_value) > 0)
        {
          mpfr_set(max_value, abs_value, GMP_RNDN);
        }
      } /* End i != j */
    } /* End for j. */
  } /* End for i. **/
  mpfr_clear(abs_value);
  return(0);
} /* End mpfr_matrix_max_abs_no_diag */
 
/**
 * @see mpfr_matrix.h#mpfr_matrix_inv_gauss_jordan
 *
 */
int 
mpfr_matrix_inv_gauss_jordan(mpfr_matrix* a,
                              mpfr_t zero_threshold)
{
  int *indxc, *indxr, *ipiv;
  int i, icol, irow, j, k , l, ll, n;
  mpfr_t big, dum, pivinv, temp, temp_abs, elem_abs, mpfr_zero, mpfr_one;

  /* Check parameters */
  if (a == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_inv_gaus_jordan",
                          "null pointer parameter",
                          "a parameter is NULL");
    return(1);
  }
  
  if ((a->NbRows == 0) || (a->NbColumns == 0))
  {
    mpfr_matrix_errormsg("mpfr_matrix_inv_gaus_jordan",
                          "empty matrix",
                          "a parameter points to an empty matrix");
    return(1);
  }
  if (a->NbRows != a->NbColumns)
  {
    mpfr_matrix_errormsg("mpfr_matrix_inv_gaus_jordan",
                          "non square matrix",
                          "\"a\" parameter points to a non square matrix");
    return(1);
  }
  n = a->NbRows;
  indxc = (int*) malloc(n * sizeof(int));
  indxr = (int*) malloc(n * sizeof(int));
  ipiv  = (int*) malloc(n * sizeof(int));
  if ((indxc == NULL) ||
      (indxr == NULL) ||
      (ipiv == NULL))
  {
    mpfr_matrix_errormsg("mpfr_matrix_inv_gaus_jordan",
                          "outofmem",
                          "out of memory space");
    if (indxc != NULL) free(indxc);
    if (indxr != NULL) free(indxr);
    if (ipiv  != NULL) free(ipiv);
    return(1);
  }
  mpfr_init(big);
  mpfr_init(dum);
  mpfr_init(pivinv);
  mpfr_init(temp);
  mpfr_init(temp_abs);
  mpfr_init(elem_abs);
  mpfr_init(mpfr_zero);
  mpfr_set_ui(mpfr_zero, 0, GMP_RNDN);
  mpfr_init(mpfr_one);
  mpfr_set_ui(mpfr_one, 1, GMP_RNDN);
  
  for (j = 0 ; j < n ; j++) ipiv[j] = 0;
  
  for (i = 0 ; i < n ; i++)
  {
    mpfr_set_ui(big, 0, GMP_RNDN);
    for (j = 0 ; j < n ; j++)
    {
      if (ipiv[j] != 1)
      {
        for (k = 0 ; k < n ; k++)
        {
          if (ipiv[k] == 0)
          {
            mpfr_abs(temp_abs, *(mpfr_matrix_get_at(a, j, k)), GMP_RNDN);
            if (mpfr_cmp(temp_abs, big) >= 0)
            {
              mpfr_set(big, temp_abs, GMP_RNDN);
              irow = j;
              icol = k;
            }
          } /* End if */
        } /* End for k */
      } /* End if */
    } /* End for j */
  
    ++(ipiv[icol]);
    
    if (irow != icol)
    {
      for (l = 0 ; l < n ; l++)
      {
        mpfr_matrix_swap_elements(a, irow, l, icol, l);
      }
    } /* End if. */
    indxr[i] = irow;
    indxc[i] = icol;
    /* if a[icol][icol] == 0, singular matrix. */
    mpfr_abs(elem_abs, *(mpfr_matrix_get_at(a, icol, icol)), GMP_RNDN);
    if (mpfr_cmp(elem_abs, zero_threshold ) <= 0)
    {
      fprintf(stderr, 
              "\nmpfr_matrix_inv_gauss_jordan: singular matrix! Giving up!\n");
      free(indxc);
      free(indxr);
      free(ipiv);
      mpfr_clear(big);
      mpfr_clear(dum);
      mpfr_clear(pivinv);
      mpfr_clear(temp);
      mpfr_clear(temp_abs);
      mpfr_clear(elem_abs);
      mpfr_clear(mpfr_zero);
      mpfr_clear(mpfr_one);
      return(1);
    }

    mpfr_div(pivinv, mpfr_one, *(mpfr_matrix_get_at(a, icol, icol)), GMP_RNDN);
    mpfr_set(*(mpfr_matrix_get_at(a, icol, icol)), mpfr_one,GMP_RNDN);

    for (l = 0 ; l < n ; l++)
    {
      mpfr_mul(*(mpfr_matrix_get_at(a, icol, l)), 
              *(mpfr_matrix_get_at(a, icol, l)),
              pivinv,
              GMP_RNDN);
    }
    for (ll = 0 ; ll < n ; ll++)
    {
      if ( ll != icol)
      {
        mpfr_set(dum, *(mpfr_matrix_get_at(a, ll, icol)), GMP_RNDN);
        mpfr_matrix_set_at(a, ll, icol, mpfr_zero, GMP_RNDN);
        for (l = 0 ; l < n ; l++)
        {
          mpfr_mul(temp,
                  *(mpfr_matrix_get_at(a, icol, l)),
                  dum,
                  GMP_RNDN);
          mpfr_sub(*(mpfr_matrix_get_at(a, ll, l)),
                  *(mpfr_matrix_get_at(a, ll, l)),
                  temp,
                  GMP_RNDN);
        }
      } /* End if. */
    } /* End for ll. */
  
  } /* End for i */
  for (l = n - 1 ; l >= 0 ; l--)
  {
    if (indxr[l] != indxc[l])
    {
      for (k = 0 ; k < n ; k++)
      {
        mpfr_matrix_swap_elements(a, k, indxr[l], k, indxc[l]);
      }
    }
  }
  free(indxc);
  free(indxr);
  free(ipiv);
  mpfr_clear(big);
  mpfr_clear(dum);
  mpfr_clear(pivinv);
  mpfr_clear(temp);
  mpfr_clear(temp_abs);
  mpfr_clear(elem_abs);
  mpfr_clear(mpfr_zero);
  mpfr_clear(mpfr_one);
  return(0);
} /* mpfr_matrix_inv_gauss_jordan */

/**
 * @see mpfr_matrix.h#mpfr_matrix_swap
 *
 */
void mpfr_matrix_swap_elements(mpfr_matrix* mat,
                                unsigned int row_1,
                                unsigned int col_1,
                                unsigned int row_2,
                                unsigned int col_2)
{
  mpfr_t *point_1, *point_2, temp;
  
  mpfr_init(temp);

  /* Only attempt the swap if both elements exist */
  point_1 = mpfr_matrix_get_at(mat, row_1, col_1);
  point_2 = mpfr_matrix_get_at(mat, row_2, col_2);
  if ((point_1 != NULL) && (point_2 != NULL))
  {
    mpfr_set(temp, *(point_1), GMP_RNDN);
    mpfr_matrix_set_at(mat, row_1, col_1, *point_2, GMP_RNDN);
    mpfr_matrix_set_at(mat, row_2, col_2, temp, GMP_RNDN);
  }
  mpfr_clear(temp);
} /* End mpfr_matrix_swap */

mpfr_matrix*
mpfr_matrix_naive_prod(mpfr_matrix* mat1,
                        mpfr_matrix* mat2,
                        mp_prec_t prec)
{
  unsigned int mat1_num_rows;
  unsigned int mat2_num_rows;
  unsigned int mat1_num_columns;
  unsigned int mat2_num_columns;
  unsigned int i, j, k;
  mpfr_matrix* prod_mat;
  mpfr_t tmp, mpfr_zero;
  
  /* Deal with pathological parameters */
  if ((mat1 == NULL) ||
      (mat2 == NULL))
  {
    mpfr_matrix_errormsg("mpfr_matrix_naive_prod", 
                          "nullPointerParameter", 
                          "one of mat1 or mat2 is a NULL pointer");
    return(NULL);
  }
  
  mat1_num_rows     = mpfr_matrix_get_num_rows(mat1);
  mat2_num_rows     = mpfr_matrix_get_num_rows(mat2);
  mat1_num_columns  = mpfr_matrix_get_num_rows(mat1);
  mat2_num_columns  = mpfr_matrix_get_num_rows(mat2);
  /* Check the dimensions compatibility */
  if (mat1_num_columns != mat2_num_rows)
  {
    mpfr_matrix_errormsg("mpfr_matrix_naive_prod", 
                          "incompatDimensions", 
                          "mat1 and mat2 have incompatible dimensions");
    return(NULL);
  }
  /* Create the product matrix */
  prod_mat = mpfr_matrix_alloc(mat1_num_rows, mat2_num_columns, prec);
  if (prod_mat == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_naive_prod", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  mpfr_init(tmp);
  mpfr_init_set_ui(mpfr_zero, 0, GMP_RNDN);
  /* Compute the product */
  for (i = 0 ; i < mat1_num_rows ; i++)
  {
    for (j = 0 ; j < mat2_num_columns ; j++)
    {
      mpfr_matrix_set_at(prod_mat, i, j, mpfr_zero, GMP_RNDN);
      for (k = 0 ; k < mat1_num_columns ; k++)
      {
        mpfr_mul(tmp, 
                  *(mpfr_matrix_get_at(mat1, i, k)),
                  *(mpfr_matrix_get_at(mat2, k, j)),
                  GMP_RNDN);
        mpfr_add(*(mpfr_matrix_get_at(prod_mat, i, j)),
                  *(mpfr_matrix_get_at(prod_mat, i, j)),
                  tmp,
                  GMP_RNDN);
      } /* En for k. */
    }
  } /* End ofr i */
  
  mpfr_clear(tmp);
  mpfr_clear(mpfr_zero);
  return(prod_mat);
}

mpfr_matrix*
mpfr_matrix_clone(mpfr_matrix* orig_matrix,
                  mp_prec_t prec)
{
  mpfr_matrix*  clone_matrix;
  unsigned int  num_rows;
  unsigned int  num_cols;
  unsigned int  i, j;
  /* Deal with pathological parameters */
  if (orig_matrix == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_clone",
                          "null pointer parameter",
                          "a parameter is NULL");
    return(NULL);
  }
  num_rows = mpfr_matrix_get_num_rows(orig_matrix);
  num_cols = mpfr_matrix_get_num_columns(orig_matrix);
  if ((num_rows == 0) || (num_cols == 0))
  {    
    mpfr_matrix_errormsg("mpfr_matrix_clone",
                          "empty matrix",
                          "a parameter points to an empty matrix");
    return(NULL);
  }
  clone_matrix = mpfr_matrix_alloc(num_rows, num_cols, prec);
  if (clone_matrix == NULL)
  {
    mpfr_matrix_errormsg("mpfr_matrix_clone",
                          "outofmem",
                          "out of memory space");
    return(NULL);
  }
  for (i = 0 ; i < num_rows ; i++)
  {
    for (j = 0 ; j < num_cols ; j++)
    {
      mpfr_matrix_set_at(clone_matrix,
                          i,
                          j,
                          *(mpfr_matrix_get_at(orig_matrix,
                                            i,
                                            j)),
                          GMP_RNDN);
    } /* End for j. */ 
  
  } /* End for i. */
  return(clone_matrix);
} /* End mpfr_matrix_clone.*/

void 
mpfr_matrix_errormsg(char *f , char *msgname, char *msg) 
{
  fprintf(stderr, "?%s: %s\n", f, msg);
}


#if 0
/* 
 * Basic hermite engine 
 */
static int hermite(mpfr_matrix *H,mpfr_matrix *U,mpfr_matrix *Q) {
  
  int nc, nr, i, j, k, rank, reduced, pivotrow;
  mpfr_t pivot,x,aux;
  mpfr_t *temp1, *temp2;
  
  /*                     T                     -1   T */
  /* Computes form: A = Q H  and U A = H  and U  = Q  */
  
  if (!H) { 
    mpfr_matrix_errormsg("Domlib", "nullH", "hermite: ? Null H");
    return -1;
  }
  nc = H->NbColumns;
  nr = H->NbRows;
  temp1 = (mpfr_t *) malloc(nc * sizeof(mpfr_t));
  temp2 = (mpfr_t *) malloc(nr * sizeof(mpfr_t));
  if (!temp1 ||!temp2) {
    mpfr_matrix_errormsg("Domlib", "outofmem", "out of memory space");
    return -1;
  }
  
  /* Initialize all the 'mpfr_t' variables */
  value_init(pivot); value_init(x); 
  value_init(aux);   
  for(i=0;i<nc;i++)
    value_init(temp1[i]);
  for(i=0;i<nr;i++)
    value_init(temp2[i]);
  
#ifdef DEBUG
  fprintf(stderr,"Start  -----------\n");
  mpfr_matrix_Print(stderr,0,H);
#endif
  for (k=0, rank=0; k<nc && rank<nr; k=k+1) {
    reduced = 1;	/* go through loop the first time */
#ifdef DEBUG
    fprintf(stderr, "Working on col %d.  Rank=%d ----------\n", k+1, rank+1);
#endif
    while (reduced) {
      reduced=0;
      
      /* 1. find pivot row */
      value_absolute(pivot,H->p[rank][k]);
      
      /* the kth-diagonal element */
      pivotrow = rank;
      
      /* find the row i>rank with smallest nonzero element in col k */
      for (i=rank+1; i<nr; i++) {
	value_absolute(x,H->p[i][k]);
	if (value_notzero_p(x) &&
	    (value_lt(x,pivot) || value_zero_p(pivot))) {
	  value_assign(pivot,x);
	  pivotrow = i;
	}
      }
      
      /* 2. Bring pivot to diagonal (exchange rows pivotrow and rank) */
      if (pivotrow != rank) {
	Vector_Exchange(H->p[pivotrow],H->p[rank],nc);
	if (U)
	  Vector_Exchange(U->p[pivotrow],U->p[rank],nr);
	if (Q)
	  Vector_Exchange(Q->p[pivotrow],Q->p[rank],nr);

#ifdef DEBUG
	fprintf(stderr,"Exchange rows %d and %d  -----------\n", rank+1, pivotrow+1);
	mpfr_matrix_Print(stderr,0,H);
#endif
      }
      value_assign(pivot,H->p[rank][k]);	/* actual ( no abs() ) pivot */
      
      /* 3. Invert the row 'rank' if pivot is negative */
      if (value_neg_p(pivot)) {
	value_oppose(pivot,pivot); /* pivot = -pivot */
	for (j=0; j<nc; j++)
	  value_oppose(H->p[rank][j],H->p[rank][j]);
	
	/* H->p[rank][j] = -(H->p[rank][j]); */
	if (U)
	  for (j=0; j<nr; j++)
	    value_oppose(U->p[rank][j],U->p[rank][j]);
	
	/* U->p[rank][j] = -(U->p[rank][j]); */
	if (Q)
	  for (j=0; j<nr; j++)
	    value_oppose(Q->p[rank][j],Q->p[rank][j]);
	
	/* Q->p[rank][j] = -(Q->p[rank][j]); */
#ifdef DEBUG
	fprintf(stderr,"Negate row %d  -----------\n", rank+1);
	mpfr_matrix_Print(stderr,0,H);
#endif

      }      
      if (value_notzero_p(pivot)) {
	
	/* 4. Reduce the column modulo the pivot */
	/*    This eventually zeros out everything below the */
	/*    diagonal and produces an upper triangular matrix */
	
	for (i=rank+1;i<nr;i++) {
	  value_assign(x,H->p[i][k]);
	  if (value_notzero_p(x)) {	    
	    value_modulus(aux,x,pivot);
	    
	    /* floor[integer division] (corrected for neg x) */
	    if (value_neg_p(x) && value_notzero_p(aux)) {
	      
	      /* x=(x/pivot)-1; */
	      value_division(x,x,pivot);
	      value_decrement(x,x);
	    }	
	    else 
	      value_division(x,x,pivot);
	    for (j=0; j<nc; j++) {
	      value_multiply(aux,x,H->p[rank][j]);
	      value_substract(H->p[i][j],H->p[i][j],aux);
	    }
	    
	    /* U->p[i][j] -= (x * U->p[rank][j]); */
	    if (U)
	      for (j=0; j<nr; j++) {
		value_multiply(aux,x,U->p[rank][j]);
		value_substract(U->p[i][j],U->p[i][j],aux);
	      }
	    
	    /* Q->p[rank][j] += (x * Q->p[i][j]); */
	    if (Q)
	      for(j=0;j<nr;j++) {
		value_multiply(aux,x,Q->p[i][j]);
		value_addto(Q->p[rank][j],Q->p[rank][j],aux);
	      }
	    reduced = 1;

#ifdef DEBUG
	    fprintf(stderr,
		    "row %d = row %d - %d row %d -----------\n", i+1, i+1, x, rank+1);
	    mpfr_matrix_Print(stderr,0,H);
#endif
	
	  } /* if (x) */
	} /* for (i) */
      } /* if (pivot != 0) */
    } /* while (reduced) */
    
    /* Last finish up this column */
    /* 5. Make pivot column positive (above pivot row) */
    /*    x should be zero for i>k */
    
    if (value_notzero_p(pivot)) {
      for (i=0; i<rank; i++) {
	value_assign(x,H->p[i][k]);
	if (value_notzero_p(x)) { 	  
	  value_modulus(aux,x,pivot);
	  
	  /* floor[integer division] (corrected for neg x) */
	  if (value_neg_p(x) && value_notzero_p(aux)) {
	    value_division(x,x,pivot);
	    value_decrement(x,x);
	    
	    /* x=(x/pivot)-1; */
	  }
	  else
	    value_division(x,x,pivot);
	  
	  /* H->p[i][j] -= x * H->p[rank][j]; */
	  for (j=0; j<nc; j++) {
	    value_multiply(aux,x,H->p[rank][j]);
	    value_substract(H->p[i][j],H->p[i][j],aux);
	  }
	  
	  /* U->p[i][j] -= x * U->p[rank][j]; */
	  if (U)
	    for (j=0; j<nr; j++) {
	      value_multiply(aux,x,U->p[rank][j]);
	      value_substract(U->p[i][j],U->p[i][j],aux);
	    }
	  
	  /* Q->p[rank][j] += x * Q->p[i][j]; */
	  if (Q)
	    for (j=0; j<nr; j++) {
	      value_multiply(aux,x,Q->p[i][j]);
	      value_addto(Q->p[rank][j],Q->p[rank][j],aux);
	    }  
#ifdef DEBUG
	  fprintf(stderr,
		  "row %d = row %d - %d row %d -----------\n", i+1, i+1, x, rank+1);
	  mpfr_matrix_Print(stderr,0,H);
#endif
	} /* if (x) */
      } /* for (i) */
      rank++;
    } /* if (pivot!=0) */
  } /* for (k) */
  
  /* Clear all the 'mpfr_t' variables */
  value_clear(pivot); value_clear(x); 
  value_clear(aux); 
  for(i=0;i<nc;i++)
    value_clear(temp1[i]);
  for(i=0;i<nr;i++)
    value_clear(temp2[i]);
  free(temp2);
  free(temp1);
  return rank;
} /* Hermite */ 

void right_hermite(mpfr_matrix *A,mpfr_matrix **Hp,mpfr_matrix **Up,mpfr_matrix **Qp) {
  
  mpfr_matrix *H, *Q, *U;
  int i, j, nr, nc, rank;
  mpfr_t tmp;
  
  /* Computes form: A = QH , UA = H */  
  nc = A->NbColumns;
  nr = A->NbRows;
  
  /* H = A */
  *Hp = H = mpfr_matrix_Alloc(nr,nc);
  if (!H) { 
    mpfr_matrix_errormsg("DomRightHermite", "outofmem", "out of memory space");
    return;
  }
  
  /* Initialize all the 'mpfr_t' variables */
  value_init(tmp);
  
  Vector_Copy(A->p_Init,H->p_Init,nr*nc);
  
  /* U = I */
  if (Up) {
    *Up = U = mpfr_matrix_Alloc(nr, nr);
    if (!U) {
      mpfr_matrix_errormsg("DomRightHermite", "outofmem", "out of memory space");
      value_clear(tmp);
      return;
    }
    Vector_Set(U->p_Init,0,nr*nr);             /* zero's */
    for(i=0;i<nr;i++)                          /* with diagonal of 1's */
      value_set_si(U->p[i][i],1);
  }
  else
    U = (mpfr_matrix *)0;
  
  /* Q = I */
  /* Actually I compute Q transpose... its easier */
  if (Qp) {
    *Qp = Q = mpfr_matrix_Alloc(nr,nr);
    if (!Q) {
      mpfr_matrix_errormsg("DomRightHermite", "outofmem", "out of memory space");
      value_clear(tmp);
      return;
    }
    Vector_Set(Q->p_Init,0,nr*nr);            /* zero's */
    for (i=0;i<nr;i++)                      /* with diagonal of 1's */
      value_set_si(Q->p[i][i],1);
  }
  else
    Q = (mpfr_matrix *)0;
  
  rank = hermite(H,U,Q);
  
  /* Q is returned transposed */ 
  /* Transpose Q */
  if (Q) {
    for (i=0; i<nr; i++) {
      for (j=i+1; j<nr; j++) {
	value_assign(tmp,Q->p[i][j]);
	value_assign(Q->p[i][j],Q->p[j][i] );
	value_assign(Q->p[j][i],tmp);
      }
    }
  }
  value_clear(tmp);
  return;
} /* right_hermite */

void left_hermite(mpfr_matrix *A,mpfr_matrix **Hp,mpfr_matrix **Qp,mpfr_matrix **Up) {
  
  mpfr_matrix *H, *HT, *Q, *U;
  int i, j, nc, nr, rank;
  mpfr_t tmp;
  
  /* Computes left form: A = HQ , AU = H , 
     T    T T    T T   T
     using right form  A  = Q H  , U A = H */
  
  nr = A->NbRows;
  nc = A->NbColumns;
  
  /* HT = A transpose */
  HT = mpfr_matrix_Alloc(nc, nr);
  if (!HT) {
    mpfr_matrix_errormsg("DomLeftHermite", "outofmem", "out of memory space");
    return;
  }
  value_init(tmp);
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      value_assign(HT->p[j][i],A->p[i][j]);
  
  /* U = I */
  if (Up) {
    *Up = U = mpfr_matrix_Alloc(nc,nc);
    if (!U) {
      mpfr_matrix_errormsg("DomLeftHermite", "outofmem", "out of memory space");
      value_clear(tmp);
      return;
    }
    Vector_Set(U->p_Init,0,nc*nc);            /* zero's */
    for (i=0;i<nc;i++)                        /* with diagonal of 1's */
      value_set_si(U->p[i][i],1);
  }
  else U=(mpfr_matrix *)0;
  
  /* Q = I */
  if (Qp) {
    *Qp = Q = mpfr_matrix_Alloc(nc, nc);
    if (!Q) {
      mpfr_matrix_errormsg("DomLeftHermite", "outofmem", "out of memory space");
      value_clear(tmp);
      return;
    }
    Vector_Set(Q->p_Init,0,nc*nc);            /* zero's */
    for (i=0;i<nc;i++)                        /* with diagonal of 1's */
      value_set_si(Q->p[i][i],1);
  }
  else Q=(mpfr_matrix *)0;
  rank = hermite(HT,U,Q);
  
  /* H = HT transpose */
  *Hp = H = mpfr_matrix_Alloc(nr,nc);
  if (!H) {
    mpfr_matrix_errormsg("DomLeftHermite", "outofmem", "out of memory space");
    value_clear(tmp);
    return;
  }
  for (i=0; i<nr; i++)
    for (j=0;j<nc;j++)
      value_assign(H->p[i][j],HT->p[j][i]);
  mpfr_matrix_Free(HT);
  
  /* Transpose U */
  if (U) {
    for (i=0; i<nc; i++) {
      for (j=i+1; j<nc; j++) {
	value_assign(tmp,U->p[i][j]);
	value_assign(U->p[i][j],U->p[j][i] );
	value_assign(U->p[j][i],tmp);
      }
    }
  }
} /* left_hermite */

/*
 * Given a integer matrix 'Mat'(k x k), compute its inverse rational matrix 
 * 'MatInv' k x (k+1). The last column of each row in matrix MatInv is used 
 * to store the common denominator of the entries in a row. The output is 1,
 * if 'Mat' is non-singular (invertible), otherwise the output is 0. Note:: 
 * (1) mpfr_matrix 'Mat' is modified during the inverse operation.
 * (2) mpfr_matrix 'MatInv' must be preallocated before passing into this function.
 */
int MatInverse(mpfr_matrix *Mat,mpfr_matrix *MatInv ) {
  
  int i, k, j, c;
  mpfr_t x, gcd, piv;
  mpfr_t m1,m2;
  
  if(Mat->NbRows != Mat->NbColumns) {
   fprintf(stderr,"Trying to invert a non-square matrix !\n");
    return 0;
  }
  
  /* Initialize all the 'mpfr_t' variables */
  value_init(x);  value_init(gcd); value_init(piv);
  value_init(m1); value_init(m2);

  k = Mat->NbRows; 

  /* Initialise MatInv */
  Vector_Set(MatInv->p[0],0,k*(k+1));

  /* Initialize 'MatInv' to Identity matrix form. Each diagonal entry is set*/
  /* to 1. Last column of each row (denominator of each entry in a row) is  */
  /* also set to 1.                                                         */ 
  for(i=0;i<k;++i) {
    value_set_si(MatInv->p[i][i],1);	
    value_set_si(MatInv->p[i][k],1);	/* denum */
  }  
  /* Apply Gauss-Jordan elimination method on the two matrices 'Mat' and  */
  /* 'MatInv' in parallel.                                                */
  for(i=0;i<k;++i) {
    
    /* Check if the diagonal entry (new pivot) is non-zero or not */
    if(value_zero_p(Mat->p[i][i])) {   	
      
      /* Search for a non-zero pivot down the column(i) */
      for(j=i;j<k;++j)      
	if(value_notzero_p(Mat->p[j][i]))
	  break;
      
      /* If no non-zero pivot is found, the matrix 'Mat' is non-invertible */
      /* Return 0.                                                         */
      if(j==k) {
	
	/* Clear all the 'mpfr_t' variables */
	value_clear(x);  value_clear(gcd); value_clear(piv);
	value_clear(m1); value_clear(m2);
	return 0;
      }	
      
      /* Exchange the rows, row(i) and row(j) so that the diagonal element */
      /* Mat->p[i][i] (pivot) is non-zero. Repeat the same operations on    */
      /* matrix 'MatInv'.                                                   */
      for(c=0;c<k;++c) {

	/* Interchange rows, row(i) and row(j) of matrix 'Mat'    */
	value_assign(x,Mat->p[j][c]);
	value_assign(Mat->p[j][c],Mat->p[i][c]);
	value_assign(Mat->p[i][c],x);
	
	/* Interchange rows, row(i) and row(j) of matrix 'MatInv' */
	value_assign(x,MatInv->p[j][c]);
	value_assign(MatInv->p[j][c],MatInv->p[i][c]);
	value_assign(MatInv->p[i][c],x);
      }
    }
    
    /* Make all the entries in column(i) of matrix 'Mat' zero except the */
    /* diagonal entry. Repeat the same sequence of operations on matrix  */
    /* 'MatInv'.                                                         */
    for(j=0;j<k;++j) {
      if (j==i) continue;	         /* Skip the pivot */
      value_assign(x,Mat->p[j][i]);
      if(value_notzero_p(x)) {
	value_assign(piv,Mat->p[i][i]);
	Gcd(x,piv,&gcd);
	if (value_notone_p(gcd) ) {
	  value_division(x,x,gcd);
	  value_division(piv,piv,gcd);
	}
	for(c=((j>i)?i:0);c<k;++c) {
	  value_multiply(m1,piv,Mat->p[j][c]);
	  value_multiply(m2,x,Mat->p[i][c]);
	  value_substract(Mat->p[j][c],m1,m2); 
	}
	for(c=0;c<k;++c) {
	  value_multiply(m1,piv,MatInv->p[j][c]);
	  value_multiply(m2,x,MatInv->p[i][c]);
	  value_substract(MatInv->p[j][c],m1,m2);
	}
	      
	/* Simplify row(j) of the two matrices 'Mat' and 'MatInv' by */
	/* dividing the rows with the common GCD.                     */
	Vector_Gcd(&MatInv->p[j][0],k,&m1);
	Vector_Gcd(&Mat->p[j][0],k,&m2);
	Gcd(m1,m2,&gcd);
	if(value_notone_p(gcd)) {
	  for(c=0;c<k;++c) {
	    value_division(Mat->p[j][c],Mat->p[j][c],gcd);
	    value_division(MatInv->p[j][c],MatInv->p[j][c],gcd);
	  }
	}
      }
    }
  }
  
  /* Simplify every row so that 'Mat' reduces to Identity matrix. Perform  */
  /* the same sequence of operations on the matrix 'MatInv'.               */
  for(j=0;j<k;++j) {
    value_assign(MatInv->p[j][k],Mat->p[j][j]);
    
    /* Make the last column (denominator of each entry) of every row greater */
    /* than zero.                                                            */
    Vector_Normalize_Positive(&MatInv->p[j][0],k+1,k); 
  }
  
  /* Clear all the 'mpfr_t' variables */
  value_clear(x);  value_clear(gcd); value_clear(piv);
  value_clear(m1); value_clear(m2);

  return 1;
} /* Mat_Inverse */

/*
 * Given (m x n) integer matrix 'X' and n x (k+1) rational matrix 'P', compute
 * the rational m x (k+1) rational matrix  'S'. The last column in each row of
 * the rational matrices is used to store the common denominator of elements
 * in a row.                              
 */
void rat_prodmat(mpfr_matrix *S,mpfr_matrix *X,mpfr_matrix *P) {
  
  int i,j,k;
  int last_column_index = P->NbColumns - 1;
  mpfr_t lcm, old_lcm,gcd,last_column_entry,s1,s2;
  mpfr_t m1,m2;
  
  /* Initialize all the 'mpfr_t' variables */
  value_init(lcm); value_init(old_lcm); value_init(gcd);
  value_init(last_column_entry); value_init(s1); value_init(s2);
  value_init(m1); value_init(m2);

  /* Compute the LCM of last column entries (denominators) of rows */
  value_assign(lcm,P->p[0][last_column_index]);	
  for(k=1;k<P->NbRows;++k) {
    value_assign(old_lcm,lcm);
    value_assign(last_column_entry,P->p[k][last_column_index]);
    Gcd(lcm,last_column_entry,&gcd);
    value_division(m1,last_column_entry,gcd);
    value_multiply(lcm,lcm,m1);
  }
  
  /* S[i][j] = Sum(X[i][k] * P[k][j] where Sum is extended over k = 1..nbrows*/
  for(i=0;i<X->NbRows;++i)
    for(j=0;j<P->NbColumns-1;++j) {
      
      /* Initialize s1 to zero. */
      value_set_si(s1,0);
      for(k=0;k<P->NbRows;++k) {
	
	/* If the LCM of last column entries is one, simply add the products */
	if(value_one_p(lcm)) {
	  value_set_si(s2,0);
	  value_multiply(s2,X->p[i][k],P->p[k][j]);
          value_addto(s1,s1,s2);
	}  
	
	/* Numerator (num) and denominator (denom) of S[i][j] is given by :- */
	/* numerator  = Sum(X[i][k]*P[k][j]*lcm/P[k][last_column_index]) and */
	/* denominator= lcm where Sum is extended over k = 1..nbrows.        */
	else {
	  value_multiply(m1,X->p[i][k],P->p[k][j]);
	  value_division(m2,lcm,P->p[k][last_column_index]);
	  value_multiply(s2,m1,m2);
	  value_addto(s1,s1,s2);
	}
      }	
      value_assign(S->p[i][j],s1);
    }
  
  for(i=0;i<S->NbRows;++i) {
    value_assign(S->p[i][last_column_index],lcm);

    /* Normalize the rows so that last element >=0 */
    Vector_Normalize_Positive(&S->p[i][0],S->NbColumns,S->NbColumns-1);
  }
  
  /* Clear all the 'mpfr_t' variables */
  value_clear(lcm); value_clear(old_lcm); value_clear(gcd);
  value_clear(last_column_entry); value_clear(s1); value_clear(s2);
  value_clear(m1); value_clear(m2);
 
  return;
} /* rat_prodmat */

/*
 * Given a matrix 'Mat' and vector 'p1', compute the matrix-vector product 
 * and store the result in vector 'p2'. 
 */
void mpfr_matrix_Vector_Product(mpfr_matrix *Mat,mpfr_t *p1,mpfr_t *p2) {

  int NbRows, NbColumns, i, j;
  mpfr_t **cm, *q, *cp1, *cp2;
  mpfr_t s;
  
  value_init(s);
  NbRows=Mat->NbRows;
  NbColumns=Mat->NbColumns;
  
  cm = Mat->p;
  cp2 = p2;
  for(i=0;i<NbRows;i++) {
    q = *cm++;
    cp1 = p1;
    value_multiply(*cp2,*q,*cp1);
    q++;
    cp1++;
    
    /* *cp2 = *q++ * *cp1++ */
    for(j=1;j<NbColumns;j++) {
      
      value_set_si(s,0);
      value_multiply(s,*q, *cp1);
      value_addto(*cp2,*cp2,s);
      q++;
      cp1++;
    }
    cp2++;
  }
  value_clear(s);
  return;
} /* mpfr_matrix_Vector_Product */

/*
 * Given a vector 'p1' and a matrix 'Mat', compute the vector-matrix product 
 * and store the result in vector 'p2'
 */
void Vector_mpfr_matrix_Product(mpfr_t *p1,mpfr_matrix *Mat,mpfr_t *p2) {
  
  int NbRows, NbColumns, i, j;
  mpfr_t **cm, *cp1, *cp2;
  mpfr_t s;
  
  value_init(s);
  NbRows=Mat->NbRows;
  NbColumns=Mat->NbColumns;
  cp2 = p2;
  cm  = Mat->p;
  for (j=0;j<NbColumns;j++) {
    cp1 = p1;
    value_multiply(*cp2,*(*cm+j),*cp1);
    cp1++;
    
    /* *cp2= *(*cm+j) * *cp1++; */
    for (i=1;i<NbRows;i++) {
      
      value_set_si(s,0);
      value_multiply(s,*(*(cm+i)+j),*cp1);
      value_addto(*cp2,*cp2,s);
      cp1++;
    }
    cp2++;
  }
  value_clear(s);
  return;
} /* Vector_mpfr_matrix_Product */

/* 
 * Given matrices 'Mat1' and 'Mat2', compute the matrix product and store in 
 * matrix 'Mat3' 
 */
void mpfr_matrix_Product(mpfr_matrix *Mat1,mpfr_matrix *Mat2,mpfr_matrix *Mat3) {
  
  int Size, i, j, k;
  unsigned NbRows, NbColumns;
  mpfr_t **q1, **q2, *p1, *p3,sum,s;
  
  NbRows    = Mat1->NbRows;
  NbColumns = Mat2->NbColumns;
  
  Size      = Mat1->NbColumns;
  if(Mat2->NbRows!=Size||Mat3->NbRows!=NbRows||Mat3->NbColumns!=NbColumns) {
    fprintf(stderr, "? mpfr_matrix_Product : incompatable matrix dimension\n");
    return;
  }     
  value_init(sum); value_init(s);
  p3 = Mat3->p_Init;
  q1 = Mat1->p;
  q2 = Mat2->p;
  
  /* Mat3[i][j] = Sum(Mat1[i][k]*Mat2[k][j] where sum is over k = 1..nbrows */
  for (i=0;i<NbRows;i++) {
    for (j=0;j<NbColumns;j++) {
      p1 = *(q1+i);
      value_set_si(sum,0);
      for (k=0;k<Size;k++) {
	
	/* sum+=*p1++ * *(*(q2+k)+j); */
	value_set_si(s,0);
	value_multiply(s,*p1, *(*(q2+k)+j));
	value_addto(sum,sum,s);
	p1++;
      }
      value_assign(*p3,sum);
      p3++;
    }
  }
  value_clear(sum); value_clear(s);
  return;
} /* mpfr_matrix_Product */
  
/*
 * Given a rational matrix 'Mat'(k x k), compute its inverse rational matrix 
 * 'MatInv' k x k. The last column of each row in matrix MatInv is used 
 * to store the common denominator of the entries in a row. The output is 1,
 * if 'Mat' is non-singular (invertible), otherwise the output is 0. Note:: 
 * (1) mpfr_matrix 'Mat' is modified during the inverse operation.
 * (2) mpfr_matrix 'MatInv' must be preallocated before passing into this function.
 */
int mpfr_matrix_Inverse(mpfr_matrix *Mat,mpfr_matrix *MatInv ) {
  
  int i, k, j, c;
  mpfr_t x, gcd, piv;
  mpfr_t m1,m2;
  mpfr_t *den;
  
  if(Mat->NbRows != Mat->NbColumns) {
   fprintf(stderr,"Trying to invert a non-square matrix !\n");
    return 0;
  }
  
  /* Initialize all the 'mpfr_t' variables */
  value_init(x);  value_init(gcd); value_init(piv);
  value_init(m1); value_init(m2);

  k = Mat->NbRows; 

  /* Initialise MatInv */
  Vector_Set(MatInv->p[0],0,k*k);

  /* Initialize 'MatInv' to Identity matrix form. Each diagonal entry is set*/
  /* to 1. Last column of each row (denominator of each entry in a row) is  */
  /* also set to 1.                                                         */ 
  for(i=0;i<k;++i) {
    value_set_si(MatInv->p[i][i],1);	
    // value_set_si(MatInv->p[i][k],1);	/* denum */
  }  
  /* Apply Gauss-Jordan elimination method on the two matrices 'Mat' and  */
  /* 'MatInv' in parallel.                                                */
  for(i=0;i<k;++i) {
    
    /* Check if the diagonal entry (new pivot) is non-zero or not */
    if(value_zero_p(Mat->p[i][i])) {   	
      
      /* Search for a non-zero pivot down the column(i) */
      for(j=i;j<k;++j)      
	if(value_notzero_p(Mat->p[j][i]))
	  break;
      
      /* If no non-zero pivot is found, the matrix 'Mat' is non-invertible */
      /* Return 0.                                                         */
      if(j==k) {
	
	/* Clear all the 'mpfr_t' variables */
	value_clear(x);  value_clear(gcd); value_clear(piv);
	value_clear(m1); value_clear(m2);
	return 0;
      }	
      
      /* Exchange the rows, row(i) and row(j) so that the diagonal element */
      /* Mat->p[i][i] (pivot) is non-zero. Repeat the same operations on    */
      /* matrix 'MatInv'.                                                   */
      for(c=0;c<k;++c) {

	/* Interchange rows, row(i) and row(j) of matrix 'Mat'    */
	value_assign(x,Mat->p[j][c]);
	value_assign(Mat->p[j][c],Mat->p[i][c]);
	value_assign(Mat->p[i][c],x);
	
	/* Interchange rows, row(i) and row(j) of matrix 'MatInv' */
	value_assign(x,MatInv->p[j][c]);
	value_assign(MatInv->p[j][c],MatInv->p[i][c]);
	value_assign(MatInv->p[i][c],x);
      }
    }
    
    /* Make all the entries in column(i) of matrix 'Mat' zero except the */
    /* diagonal entry. Repeat the same sequence of operations on matrix  */
    /* 'MatInv'.                                                         */
    for(j=0;j<k;++j) {
      if (j==i) continue;	         /* Skip the pivot */
      value_assign(x,Mat->p[j][i]);
      if(value_notzero_p(x)) {
	value_assign(piv,Mat->p[i][i]);
	Gcd(x,piv,&gcd);
	if (value_notone_p(gcd) ) {
	  value_division(x,x,gcd);
	  value_division(piv,piv,gcd);
	}
	for(c=((j>i)?i:0);c<k;++c) {
	  value_multiply(m1,piv,Mat->p[j][c]);
	  value_multiply(m2,x,Mat->p[i][c]);
	  value_substract(Mat->p[j][c],m1,m2); 
	}
	for(c=0;c<k;++c) {
	  value_multiply(m1,piv,MatInv->p[j][c]);
	  value_multiply(m2,x,MatInv->p[i][c]);
	  value_substract(MatInv->p[j][c],m1,m2);
	}
	      
	/* Simplify row(j) of the two matrices 'Mat' and 'MatInv' by */
	/* dividing the rows with the common GCD.                     */
	Vector_Gcd(&MatInv->p[j][0],k,&m1);
	Vector_Gcd(&Mat->p[j][0],k,&m2);
	Gcd(m1,m2,&gcd);
	if(value_notone_p(gcd)) {
	  for(c=0;c<k;++c) {
	    value_division(Mat->p[j][c],Mat->p[j][c],gcd);
	    value_division(MatInv->p[j][c],MatInv->p[j][c],gcd);
	  }
	}
      }
    }
  }
  
  /* Find common denom for each row */ 
   den = (mpfr_t *)malloc(k*sizeof(mpfr_t));
   value_set_si(x,1);
   for(j=0 ; j<k ; ++j) {
     value_init(den[j]);
     value_assign(den[j],Mat->p[j][j]);
     
     /* gcd is always positive */
     Vector_Gcd(&MatInv->p[j][0],k,&gcd);
     Gcd(gcd,den[j],&gcd);
     if (value_neg_p(den[j])) 
       value_oppose(gcd,gcd); /* make denominator positive */
     if (value_notone_p(gcd)) {
       for (c=0; c<k; c++) 
	 value_division(MatInv->p[j][c],MatInv->p[j][c],gcd); /* normalize */
       value_division(den[j],den[j],gcd);
     }  
     Gcd(x,den[j],&gcd);
     value_division(m1,den[j],gcd);
     value_multiply(x,x,m1);
   }
   if (value_notone_p(x)) 
     for(j=0 ; j<k ; ++j) {       
       for (c=0; c<k; c++) {
	 value_division(m1,x,den[j]);
	 value_multiply(MatInv->p[j][c],MatInv->p[j][c],m1);  /* normalize */
       }
     }

   /* Clear all the 'mpfr_t' variables */
   for(j=0 ; j<k ; ++j) {
     value_clear(den[j]);
   }  
   value_clear(x);  value_clear(gcd); value_clear(piv);
   value_clear(m1); value_clear(m2);
   free(den);
   
   return 1;
} /* mpfr_matrix_Inverse */

#endif /* #if 0 */






