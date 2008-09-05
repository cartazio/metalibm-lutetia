/*
 * Copyright 2008 by 
 * 
 * Laboratoire de l'Informatique du Parallélisme, 
 * UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668
 *
 * Utilities for the vandercoeff program.
 *
 * Contributor: Serge Torres (ENS Lyon) -- serge.torres@ens-lyon.fr
 *
 * This file is an integrated part of the metalibm library developed by the 
 * Arenaire project at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

/* utils.c */
/******************************************************************************/

/**
 * Name & purpose                                                           
 * Author:                                                                  
 *
 */

/* includes of system headers */
#include <stdio.h>
#include <mpfr.h>
/* includes of project headers */

#include "./vectors/mpfr/mpfr-vector.h"
#include "./vectors/mpz/mpz-vector.h"
#include "./vectors/ui/ui-vector.h"
#include "./matrices/mpfr/mpfr-matrix.h"
#include "./matrices/mpz/mpz-matrix.h"
#include "./matrices/ui/ui-matrix.h"
/* includes of local headers */

#include "utils.h"

/* Types, constants and macros definitions */

/* Global variables */

/* Local function prototypes */

/**
 *
 *
 */
void
vander_coeffs_errormsg(char *f , char *msgname, char *msg);

/* Function definitions */

/**
 *@see utils.h#create_vandermonde_style_matrix
 *
 */

mpfr_matrix*  
create_vandermonde_style_matrix(unsigned int degree, 
                                mpfr_vector* vect, 
                                unsigned int precision)
{
  mpfr_matrix* mat          = NULL;
  mpfr_vector* powers_vect  = NULL;
  unsigned int nb_rows      = 0;
  unsigned int nb_columns   = 0;
  unsigned int i;
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    vander_coeffs_errormsg("create_vandermonde_style_matrix",
                            "nullPointerParameter",
                            "the generator vector is a NULL pointer");
    return(NULL);
  }
  nb_rows = mpfr_vector_get_size(vect);
  if (nb_rows == 0)
  {
    vander_coeffs_errormsg("create_vandermonde_style_matrix",
                            "emptyGenVector",
                            "the generator vector is empty");
    return(NULL);
  }
  nb_columns = degree + 1;
  if (nb_columns < 1 )
  {
    vander_coeffs_errormsg("create_vandermonde_style_matrix",
                            "invalidDegree",
                            "degree is invalid (should be >= 1)");
    return(NULL);
  }
  mat = mpfr_matrix_alloc(nb_rows, nb_columns, precision);
  if (mat == NULL)
  {
    vander_coeffs_errormsg("create_vandermonde_style_matrix", 
                            "outofmem", 
                            "out of memory space");
  }
  /* */
  for(i = 0 ; i < nb_rows ; i++)
  {
    powers_vect = mpfr_first_powers_of(*(mpfr_vector_get_at(vect, i)), 
                                        degree, 
                                        precision);
    mpfr_matrix_set_row_from_vector(mat, i, powers_vect);
    mpfr_vector_free(powers_vect);
  }
  return(mat);
} /* End create_vandermonde_style_matrix */

/**
 *@see utils.h#create_vandermonde_style_matrix
 *
 */

mpfr_matrix*  
create_sparse_vandermonde_style_matrix(mpfr_vector* points,
                                        ui_vector* exponents, 
                                        unsigned int precision)
{
  mpfr_matrix* mat          = NULL;
  mpfr_vector* powers_vect  = NULL;
  unsigned int nb_rows      = 0;
  unsigned int nb_columns   = 0;
  unsigned int i;
  /* Deal with pathological parameters */
  if ((points == NULL) || (exponents == NULL))
  {
    vander_coeffs_errormsg("create_sparse_vandermonde_style_matrix",
                            "nullPointerParameter",
                            "the generator vector or the exponent vector is a NULL pointer");
    return(NULL);
  }
  nb_rows = mpfr_vector_get_size(points);
  if (nb_rows == 0)
  {
    vander_coeffs_errormsg("create_sparse_vandermonde_style_matrix",
                            "emptyGenVector",
                            "the generator vector is empty");
    return(NULL);
  }
  nb_columns = ui_vector_get_size(exponents);
  if (nb_columns < 1 )
  {
    vander_coeffs_errormsg("create_sparse_vandermonde_style_matrix",
                            "invalidExponentsNum",
                            "the numver so exponent is invalid (should be >= 1)");
    return(NULL);
  }
  mat = mpfr_matrix_alloc(nb_rows, nb_columns, precision);
  if (mat == NULL)
  {
    vander_coeffs_errormsg("create_sparse_vandermonde_style_matrix", 
                            "outofmem", 
                            "out of memory space");
  }
  /* */
  for(i = 0 ; i < nb_rows ; i++)
  {
    powers_vect = mpfr_arbitrary_int_powers_of(*(mpfr_vector_get_at(points, i)), 
                                              exponents, 
                                              precision);
    mpfr_matrix_set_row_from_vector(mat, i, powers_vect);
    mpfr_vector_free(powers_vect);
  }
  return(mat);
} /* End create_sparse_vandermonde_style_matrix */

/**
 * @see utils.h#first_powers_of
 *
 */
/**
 * @see utils.h#mpfr_first_powers_of
 *
 */

mpfr_vector*
mpfr_first_powers_of(mpfr_t x, unsigned int max_power, unsigned int precision)
{
  mpfr_vector* new_vector = NULL;
  mpfr_t current_power;
  unsigned int i;
  /* Deal with pathological parameters */
  if (max_power < 0)
  {
    vander_coeffs_errormsg("first_powers_of", 
                            "invalidMaxPow", 
                            "invalid maximum power (must be >= 0)");
    return(NULL);
  }

  mpfr_set_default_prec(precision);
  new_vector = mpfr_vector_alloc(max_power + 1, precision);

  if(new_vector == NULL) 
  {
    vander_coeffs_errormsg("first_powers_of", 
                            "outofmem", 
                            "out of memory space");
    return NULL;
  }
  mpfr_init_set_ui(current_power, 1, GMP_RNDN);
  // x^0 = 1
  mpfr_vector_set_at(new_vector,
                      0,
                      current_power,
                      GMP_RNDN);
  mpfr_init_set(current_power, x, GMP_RNDN);
  for (i = 1 ; i <= max_power ; i++)
  {
    mpfr_vector_set_at(new_vector,
                        i,
                        current_power,
                        GMP_RNDN);
    mpfr_mul(current_power, current_power, x, GMP_RNDN);
  }
  return(new_vector);
} /* End mpfr_first_powers_of */

/**
 * @see utils.h#mpfr_arbitrary_int_powers_of
 *
 */
mpfr_vector*
mpfr_arbitrary_int_powers_of(mpfr_t x,
                          ui_vector* powers_exponents,
                          unsigned int precision)
{
  mpfr_vector* powers_values = NULL;
  int i;
  unsigned int nb_elems;
  mpfr_t temp;
  /* Parameters check */
  if (powers_exponents == NULL)
  {
    vander_coeffs_errormsg("mpfr_arbitrary_int_powers_of", 
                            "nullPointerParam", 
                            "the powers_exponents vector pointer is NULL");
    return(NULL);
  }
  if (powers_exponents->NbElems == 0)
  {
    return(NULL);
  }
  else
  {
    nb_elems = powers_exponents->NbElems;
  }
  if (precision == 0)
  {
    precision = mpfr_get_default_prec();
  }
  /* Start working */
  powers_values = mpfr_vector_alloc(nb_elems,
                                    precision);
  if (powers_values == NULL)
  {
    vander_coeffs_errormsg("mpfr_arbitrary_int_powers_of", 
                            "outofmem", 
                            "the powers_values vector could not be allocated.");
    return(NULL);
  }
  mpfr_init(temp);
  for (i = 0 ; i < nb_elems ; i++)
  {
    mpfr_pow_ui(temp,
                x,
                *(ui_vector_get_at(powers_exponents, i)),
                GMP_RNDN);
    mpfr_vector_set_at(powers_values, i, temp, GMP_RNDN);
    /*mpfr_pow_si(*(mpfr_vector_get_at(powers_values, i)),
                x,
                *(si_vector_get_at(powers_exponents, i)),
                GMP_RNDN);
    mpfr_vector_set_at(powers_values, i, temp, GMP_RNDN);*/
  }
  mpfr_clear(temp);
  return(powers_values);
} /* End mpfr_arbitrary_int_powers_of */
 
 /**
 * 
 *
 *
 */

ui_vector*
successor(ui_vector* vect, unsigned int n, unsigned int k)
{
  unsigned int i, j;
  unsigned int offset;
  unsigned int current_value;
  /* We scan vector elements from the last to the first */
  for (i = k - 1 ; i >= 0 ; i--)
  {
    current_value = *(ui_vector_get_at(vect, i));
    /* According to their position into the vector, elements
       can not go beyond a given value. The offset variable
       makes it easy to compute it from the current index i. */
    offset = k - i;
    /* If an element is already beyond it's largest possible
       value, it is impossible to compute a valid successor
       for it. */
    if (current_value > (n - offset))
    {
      vander_coeffs_errormsg("successor", 
                              "invalidParam", 
                              "one element is already too large");
    }
    /* If an element has reached it's largest possible value,
       move to the previous one (continue) unless we already 
       are at the first element (return (NULL)). */
    if (current_value == (n - offset))
    {
      if (i == 0)
      {
        return(NULL);
      }
      else
      {
        continue;
      }
    }
    /* If the element has not reached it's maximum value, increment it.
       and set all the following elements to their "correct value" 
       (previous_element + 1). */
    else
    {
      current_value ++;
      ui_vector_set_at(vect, i, current_value);
      for (j = i + 1 ; j < k ; j++)
      {
        current_value ++;
        ui_vector_set_at(vect, j, current_value);
      } /* End for j */
      return(vect);
    }
  } /* End for i */
} /* End successor */

/**
 *
 *
 */
ui_matrix*
generate_combinations(unsigned int n, unsigned int k)
{
  ui_matrix* combinations = NULL;
  unsigned int i, ui_num_combinations;
  ui_vector* initial_vector = NULL;
  ui_vector* current_combination_vector = NULL;
  mpz_t mpz_num_combinations;
  
  /* Deal with pathological parameters */
  if (n <= 0)
  {
    vander_coeffs_errormsg("generate_combinations", 
                            "invalidParam", 
                            "n must be > 0");
    return(NULL);
  }
  if (k <= 0)
  {
    vander_coeffs_errormsg("generate_combinations", 
                            "invalidParam", 
                            "k must be > 0");
    return(NULL);
  }
  if (k > n)
  {
    vander_coeffs_errormsg("generate_combinations", 
                            "invalidParam", 
                            "k must be <= n");
    return(NULL);
  }
  /* Create the initial vector. */
  initial_vector = ui_vector_alloc(k);
  if (initial_vector == NULL)
  {
    vander_coeffs_errormsg("generate_combinations", 
                            "outofmem", 
                            "out of memory space");
    return(NULL);
  }
  for (i = 0 ; i < k ; i++ )
  {
    ui_vector_set_at(initial_vector, i, i);
  }
  /* Compute the number of combinations */
  /* We use GMP :
     - because it implements an efficient algorithm to compute
       the number of combinations;
     - to avoid an overflow in intermediate computations.
     We finaly convert the result into an unsigned int since,
     if, at this point, an overflow happens, the rest of the
     algorithm is pointless (an impratical computation time). */
  mpz_init(mpz_num_combinations);
  mpz_bin_uiui(mpz_num_combinations, n, k);
  ui_num_combinations = mpz_get_ui(mpz_num_combinations);
  mpz_clear(mpz_num_combinations);
  /* Create the combinations matrix */
  combinations = ui_matrix_alloc(ui_num_combinations, k);
  if (combinations == NULL)
  {
    vander_coeffs_errormsg("generate_combinations", 
                            "outofmem", 
                            "out of memory space");
    return(NULL);
  }
  current_combination_vector = initial_vector;
  i = 0;
  while(current_combination_vector != NULL)
  {
    if (ui_matrix_set_row_from_vector(combinations, 
                                      i, 
                                      current_combination_vector))
    {
      vander_coeffs_errormsg("generate_combinations", 
                              "combigenfailure", 
                              "could not generated the combinations matrix");
      ui_matrix_free(combinations);
      ui_vector_free(initial_vector);
      return(NULL);
    }
    current_combination_vector = successor(initial_vector, n, k);
    i++;
  }
  ui_vector_free(initial_vector);
  return(combinations);
} /* End generate_combinations */

/**
 *
 *
 */

mpfr_matrix*
mpfr_sub_matrix(mpfr_matrix* mat, 
                ui_vector* vect,
                unsigned precision)
{
  unsigned int i, vect_size, current_index;
  mpfr_matrix* sub_matrix     = NULL;
  mpfr_vector* current_vector = NULL;
  
  /* Deal with pathological cases */
  if (mat == NULL)
  {
    vander_coeffs_errormsg("mpfr_sub_matrix", 
                            "invalid param", 
                            "mat should not be a NULL pointer");
    return(NULL);
  }
    
  if (vect == NULL)
  {
    vander_coeffs_errormsg("mpfr_sub_matrix", 
                            "invalid param", 
                            "vect should not be a NULL pointer");
    return(NULL);
  }

  vect_size = ui_vector_get_size(vect);
  if (vect_size == 0)
  {
    return(NULL);
  }
  
  if ((mpfr_matrix_get_num_rows(mat) == 0) ||
        (mpfr_matrix_get_num_columns(mat) == 0))
  {
    return(NULL);
  }
  /* Usefull stuff now. */
  sub_matrix = mpfr_matrix_alloc(vect_size,
                                  mpfr_matrix_get_num_columns(mat),
                                  precision);
  if (sub_matrix == NULL)
  {
    vander_coeffs_errormsg("mpfr_sub_matrix", 
                            "outofmem", 
                            "out of memory space");
    return(NULL);
  }
  for (i = 0 ; i < vect_size ; i++)
  {
    current_index = *(ui_vector_get_at(vect, i));
    current_vector = 
        mpfr_matrix_get_vector_from_row(mat,
                                        *(ui_vector_get_at(vect, i)),
                                        precision);
    mpfr_matrix_set_row_from_vector(sub_matrix, i, current_vector);
    mpfr_vector_free(current_vector);
  }
  return(sub_matrix);
}


/**
 *
 *
 */
  
mpfr_vector*
mpfr_sub_vector(mpfr_vector* mpfr_vect, 
                ui_vector* ui_vect,
                unsigned precision)
{
  unsigned int i, vect_size, current_index;
  mpfr_vector* sub_vector     = NULL;
  
  /* Deal with pathological cases */
  if (mpfr_vect == NULL)
  {
    vander_coeffs_errormsg("mpfr_sub_vector", 
                            "invalid param", 
                            "mpfr_vect should not be a NULL pointer");
    return(NULL);
  }
    
  if (ui_vect == NULL)
  {
    vander_coeffs_errormsg("mpfr_sub_vector", 
                            "invalid param", 
                            "ui_vect should not be a NULL pointer");
    return(NULL);
  }

  vect_size = ui_vector_get_size(ui_vect);
  if (vect_size == 0)
  {
    return(NULL);
  }
  
  if (mpfr_vector_get_size(mpfr_vect) == 0)
  {
    return(NULL);
  }
  
  /* Some usefull stuff now. */
  sub_vector = mpfr_vector_alloc(vect_size,
                                  precision);
  if (sub_vector == NULL)
  {
    vander_coeffs_errormsg("mpfr_sub_vector", 
                            "outofmem", 
                            "out of memory space");
    return(NULL);
  }
  for (i = 0 ; i < vect_size ; i++)
  {
    current_index = *(ui_vector_get_at(ui_vect, i));
    if ((current_index >= mpfr_vector_get_size(mpfr_vect)) ||
          (current_index < 0))
    {
      vander_coeffs_errormsg("mpfr_sub_vector", 
                              "accessBeyondBounds", 
                              "trying to access mpfr_vect beyond bounds");
    }
    mpfr_vector_set_at(sub_vector,
                        i,
                        *(mpfr_vector_get_at(mpfr_vect, current_index)),
                        GMP_RNDN);
  } /* End for i */
  return(sub_vector);
}

/**
 *
 * The algorithm used here follows the outlines given in
 * the N.J.\ Higham's  Accuracy and Stability of Numerical 
 * Algorithms, pages 426-427.
 *
 * To compute the master polynomial according to the algorithm
 * you have to go one step further than it's actual or usefull 
 * size. Here, we introduce an extra test to avoid going "a step to far".
 * @param mat - the matrix;
 * @param precision - the precision, in bits, of the computations.
 * @return a vector (or NULL if something went wrong).
*/

mpfr_vector*
mpfr_matrix_build_master_polynomial(mpfr_matrix* mat, unsigned precision)
{
  unsigned int num_rows;
  unsigned int num_columns;
  unsigned int k, j, last_index;
  mpfr_t mpfr_one;
  mpfr_vector* master_polynomial;
  /* Deal with pathological parameters */
  if (mat == NULL)
  {
    vander_coeffs_errormsg("mpfr_matrix_build_master_polynomial", 
                            "invalid param", 
                            "mat should not be a NULL pointer");
    return(NULL);
  }
  num_rows    = mpfr_matrix_get_num_rows(mat);
  num_columns = mpfr_matrix_get_num_columns(mat);
  if ((num_rows <= 0) || 
      (num_columns <= 0))
  {
    vander_coeffs_errormsg("mpfr_matrix_build_master_polynomial", 
                            "invalid param", 
                            "mat should have a strictly positive number of rows and columns");
    return(NULL);
  }
  master_polynomial = mpfr_vector_alloc(num_rows, precision);
  if (master_polynomial == NULL)
  {
    vander_coeffs_errormsg("mpfr_matrix_build_master_polynomial", 
                            "outofmem", 
                            "out of memory space for master polynomial");
  }
  mpfr_init_set_ui(mpfr_one, 1, GMP_RNDN);
  /* The matrix generator is the 1 column index element of each row
     of the Vandermonde matrix. */
  // a0 = -alpha0
  mpfr_vector_set_at(master_polynomial,
                      0,
                      *(mpfr_matrix_get_at(mat, 0, 1)), 
                      GMP_RNDN);
  mpfr_neg(*(mpfr_vector_get_at(master_polynomial,
                               0)),
           *(mpfr_vector_get_at(master_polynomial,
                               0)),
           GMP_RNDN);
     
  // a1 = 1
  mpfr_vector_set_at(master_polynomial,
                      1, 
                      mpfr_one,
                      GMP_RNDN);
  
  last_index = num_columns - 1;
  /* The matrix generator is the 1 column index element of each row
     of the Vandermonde matrix. */
  for (k = 1 ; k < num_columns ; k++)
  {
    // a[k + 1] = 1. The test prevents an access beyond vector bounds.
    if (k < last_index)
    {
      mpfr_vector_set_at(master_polynomial,
                          k + 1, 
                          mpfr_one,
                          GMP_RNDN);
    }
    for (j = k ; j >= 1 ; j--)
    {
      // alpha[k] * a[j]
      mpfr_mul(*(mpfr_vector_get_at(master_polynomial, j)),
                *(mpfr_matrix_get_at(mat, k, 1)),
                *(mpfr_vector_get_at(master_polynomial, j)),
                GMP_RNDN);
      // a[j] = a[j - 1] - alpha[k] * a[j]
      mpfr_sub(*(mpfr_vector_get_at(master_polynomial, j)),
                *(mpfr_vector_get_at(master_polynomial, j - 1)),
                *(mpfr_vector_get_at(master_polynomial, j)),
                GMP_RNDN);
    } /* End for j */
    // a0 = - (alpha[k] * a[0])
    mpfr_mul(*(mpfr_vector_get_at(master_polynomial, 0)),
              *(mpfr_matrix_get_at(mat, k, 1)),
              *(mpfr_vector_get_at(master_polynomial, 0)),
              GMP_RNDN);
    mpfr_neg(*(mpfr_vector_get_at(master_polynomial, 0)),
              *(mpfr_vector_get_at(master_polynomial, 0)),
              GMP_RNDN);
  } /* End for k */
  mpfr_clear(mpfr_one);
  return(master_polynomial);
} /* End mpfr_matrix_build_master_polynomial. */

/**
 *
 * The algorithm used here follows the outlines given in
 * the N.J.\ Higham's  Accuracy and Stability of Numerical 
 * Algorithms, pages 426-427.
 *
 * Beware that the definition of the Vandermonde matrix used
 * in this reference is the transposed matrix of that we
 * use here. Hence all row/column indexes of the inv_mat
 * are switched.
 *
 * @param mat       - a pointer to the original matrix;
 * @param precision - the precision, in bits, of the inverse
 *                    matrix elements.
 * @return a pointer to the inverse matrix (or NULL if
 * something went wrong).

 */
mpfr_matrix* 
mpfr_matrix_inv_vander(mpfr_matrix* mat, unsigned precision)
{
  mpfr_vector* master_polynomial = NULL;
  mpfr_matrix* inv_mat = NULL;
  mpfr_t s, tmp;
  int i, j;
  unsigned int num_rows;
  unsigned int num_columns;
  /* Deal with pathological parameters */
  if (mat == NULL)
  {
    vander_coeffs_errormsg("mpfr_matrix_inv_vander", 
                            "invalid param", 
                            "mat should not be a NULL pointer");
    return(NULL);
  }
  num_rows    = mpfr_matrix_get_num_rows(mat);
  num_columns = mpfr_matrix_get_num_columns(mat);
  if ((num_rows <= 0) || 
      (num_columns <= 0))
  {
    vander_coeffs_errormsg("mpfr_matrix_inv_vander", 
                            "invalid param", 
                            "mat should have a strictly positive number of rows and columns");
    return(NULL);
  }
  /* Not a square matrix ! */
  if (num_rows != num_columns)
  {
    vander_coeffs_errormsg("mpfr_matrix_inv_vander", 
                            "invalid param", 
                            "mat should be a square matrix");
    return(NULL);
  }
  master_polynomial = 
      mpfr_matrix_build_master_polynomial(mat, precision);
  if (master_polynomial == NULL)
  {
    vander_coeffs_errormsg("mpfr_matrix_inv_vander", 
                            "outofmem", 
                            "out of memory space for master polynomial");
    return(NULL);
  }
  /* 
  fprintf(stderr, "Master polynomial:\n");
  mpfr_vector_print(stderr, 10, 5, GMP_RNDN, master_polynomial); */
  inv_mat = mpfr_matrix_alloc(num_rows, num_columns, precision);
  if (inv_mat == NULL)
  {  
    vander_coeffs_errormsg("mpfr_matrix_inv_vander", 
                            "outofmem", 
                            "out of memory space for inverse matrix");
    mpfr_vector_free(master_polynomial);
    return(NULL);
  }
  mpfr_init(s);
  mpfr_init(tmp);
  for (i = 0 ; i < num_columns ; i++)
  {
    // s = 1
    mpfr_set_ui(s, 1, GMP_RNDN);
    //w[i,n] = 1
    mpfr_matrix_set_at(inv_mat, num_columns - 1, i, s, GMP_RNDN);
    for (j = (num_columns - 2) ; j >= 0 ; j--)
    {
      // alpha[i] * w[i,j + 1]
      mpfr_mul(tmp, 
                *(mpfr_matrix_get_at(mat, i, 1)),         // alpha[i]
                *(mpfr_matrix_get_at(inv_mat, j + 1, i)), // w[i, j + 1]
                GMP_RNDN);
      // a[j + 1] + alpha[i] * w[i,j + 1]
      mpfr_add(tmp, 
                *(mpfr_vector_get_at(master_polynomial, j + 1)),
                tmp,
                GMP_RNDN);
      // w[i,j] = a[j + 1] + alpha[i] * w[i,j + 1]
      mpfr_matrix_set_at(inv_mat, j, i, tmp, GMP_RNDN);
      // s = alpha[i] * s
      mpfr_mul(s, *(mpfr_matrix_get_at(mat, i, 1)), s, GMP_RNDN); 
      // s = alpha[i] * s + w[i,j]
      mpfr_add(s, s, *(mpfr_matrix_get_at(inv_mat, j, i)), GMP_RNDN); 
    } /* End for j */ 
    for (j = 0 ; j < num_columns ; j++)
    {
      mpfr_div(*(mpfr_matrix_get_at(inv_mat, j, i)),
                *(mpfr_matrix_get_at(inv_mat, j, i)),
                s,
                GMP_RNDN);
    } /* End for j. */
  } /* End for i. */
  mpfr_vector_free(master_polynomial);
  mpfr_clear(s); 
  mpfr_clear(tmp); 
  return(inv_mat);
  
} /* End mpfr_matrix_inv_vander */

/**
 *
 *
 */
void 
vander_coeffs_errormsg(char *f , char *msgname, char *msg) 
{
  fprintf(stderr, "?%s: %s\n", f, msg);
} /* End vander_coeffs_errormsg */



