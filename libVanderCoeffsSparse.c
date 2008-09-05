/*
 * Copyright 2008 by 
 * 
 * Laboratoire de l'Informatique du Parallélisme, 
 * UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668
 *
 * The vandercoeff program.
 *
 * This software is based on scientific research made by
 * Serge Torres and Nicolas Brisebarre.
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

/* libVanderCoeffsSparce.c */
/******************************************************************************/
/**
 * Name & purpose                                                           
 * @author ST
 * @date 2008-02-14
 *
 */

/* includes of system headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

/* includes of project headers */
#include "./vectors/mpfr/mpfr-vector.h"
#include "./matrices/mpfr/mpfr-matrix.h"
#include "./vectors/ui/ui-vector.h"
#include "./matrices/ui/ui-matrix.h"
#include "./misc-utils/std-exit-errors.h"
#include "./utils.h"


/* includes of local headers */

#include "./libVanderCoeffsSparse.h"

/* Types, constants and macros definitions */
#define RADIX 10

/* Global variables */

/* Functions */

/**
 * @see libVanderCoeffsSparse.h#vanderCoeffsSparseAbs
 */

int 
vanderCoeffsSparseAbs(unsigned int  points_count,
                      mpfr_t*       points_array, 
                      mpfr_t*       func_evals_array,
                      unsigned int  degrees_count,
                      unsigned int* degrees_array,
                      mpfr_t        max_err,
                      mpfr_t*       lower_bounds_array,
                      mpfr_t*       upper_bounds_array)
{
  mpfr_matrix*  vander_style_matrix     = NULL;
  mpfr_matrix*  calibration_matrix_1    = NULL;
  mpfr_matrix*  calibration_matrix_2    = NULL;
  mpfr_matrix*  calibration_matrix_3    = NULL;
  mpfr_t        calibration_number;
  mpfr_vector*  mpfr_points_vector      = NULL;
  mpfr_vector*  func_evaluation         = NULL;
  mpfr_vector*  abs_func_evaluation     = NULL;
  mpfr_vector*  sub_func_evaluation     = NULL;
  mpfr_vector*  eval_upper_bounds       = NULL;
  mpfr_vector*  eval_lower_bounds       = NULL;
  mpfr_vector*  coeff_upper_bounds      = NULL;
  mpfr_vector*  coeff_lower_bounds      = NULL;
  mpfr_matrix*  current_sub_matrix      = NULL;
  mpfr_matrix*  inv_current_sub_matrix  = NULL;
  mpfr_vector*  current_inv_mat_line    = NULL;
  mpfr_vector*  current_abs_inv_mat_line= NULL;
  mpfr_vector*  current_upper_bounds    = NULL;
  mpfr_vector*  current_lower_bounds    = NULL;
  ui_matrix*    combinations            = NULL;
  ui_vector*    ui_degrees_vector       = NULL;
  ui_vector*    current_combination     = NULL;
  mpfr_t        scal_prod;
  mpfr_t        correction_term;
  mp_prec_t     precision;
  unsigned int  combinations_count      = 0;
  int           i, j;
  /* Check the parameters */
  /* Use the current default precision whenever a precision is needed. */
  precision = mpfr_get_default_prec();
  /* The number of points must be >= degree+1. */
  if (points_count < degrees_count)
  {
    /* Nothing to deallocate. */
    return(1); 
  }
  /* The input arrays must not be NULL. */
  if ((points_array       == NULL)  || 
      (func_evals_array   == NULL)  ||
      (degrees_array      == NULL)  ||
      (upper_bounds_array == NULL)  ||
      (lower_bounds_array == NULL))
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: one of the parameter pointers is NULL.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    /* Nothing to deallocate. */
    return(1); 
  }
  /* The max_err must be > 0 */
  if (mpfr_sgn(max_err) <= 0)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: the error is not > 0.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    /* Nothing to deallocate. */
    return(1); 
  }
  /* The number of degrees must be > 0 */
  if (degrees_count <= 0)
  {
    #ifndef NDEBUG
      fprintf(stderr,
            "%s: the degrees count is not > 0.\n",
            "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    /* Nothing to deallocate. */
    return(1); 
  }
  /* End input data check. */
  /* Create the points vector */
  mpfr_points_vector = 
        mpfr_vector_from_array(points_array, points_count, precision);
  if (mpfr_points_vector == NULL)
  {
    /* Nothing to deallocate, so far. */
    return(2); 
  }
  /* Check that the points vector has no duplicates. This is still input data
     check but it can not conveniently be done before. */
  if (mpfr_vector_has_duplicates(mpfr_points_vector))
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: the points array has duplicates.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    return(1); 
  }
  /* Create the func_evaluation vector. */
  func_evaluation = mpfr_vector_from_array(func_evals_array, 
                                            points_count, 
                                            precision);
  if (func_evaluation == NULL)
  {
    mpfr_vector_free(mpfr_points_vector);
    return(2); 
  }
  #ifndef NDEBUG
    mpfr_vector_print(stderr, 
                      RADIX,
                      5,
                      GMP_RNDN,
                      func_evaluation);
  #endif
  /* Compute the degrees vector */
  ui_degrees_vector = ui_vector_from_array(degrees_array,
                                            degrees_count);
  if (ui_degrees_vector == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not create the ui_degrees_vector.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    return(2);
  }
  /* Check that the degrees vector has no duplicates. This is still input data
     check but it can not conveniently be done before. */
  if (ui_vector_has_duplicates(ui_degrees_vector))
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: the degrees array has duplicates.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    return(1); 
  }
  /* Compute the absolute value of the function evaluation vector. */
  abs_func_evaluation = mpfr_vector_get_abs(func_evaluation, precision);
  if (abs_func_evaluation == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not create the absolute value of the function evaluation vector.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    return(2); 
  }

  /*
   * Calibration stuff
   */
  /* Create a Vandermonde matrix for calibration */
  calibration_matrix_1 = 
      create_vandermonde_style_matrix(mpfr_vector_get_size(mpfr_points_vector) - 1,
                                      mpfr_points_vector,
                                      precision);
  if (calibration_matrix_1 == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not allocate callibration_matrix_1.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    return(2);
  }
  /* Inverse the matrix */
  calibration_matrix_2 = mpfr_matrix_inv_vander(calibration_matrix_1,
                                                  precision);
  if (calibration_matrix_2 == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not compute the inverse matrix callibration_matrix_2.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_matrix_free(calibration_matrix_1);
    return(2);
  }
  /* Compute the product of a matrix and it's inverse */
  calibration_matrix_3 = mpfr_matrix_naive_prod(calibration_matrix_1,
                                                calibration_matrix_2,
                                                precision);
  if (calibration_matrix_3 == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not compute the identity matrix callibration_matrix_3.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_matrix_free(calibration_matrix_1);
    mpfr_matrix_free(calibration_matrix_2);
    return(2);
  }
  /* Get the maximum absolute value out of the diagonal.
   * This should be "the maximum value of 0", kinf of speak. */
  mpfr_init(calibration_number);
  if (mpfr_matrix_max_abs_no_diag(calibration_matrix_3, calibration_number))
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not compute the calibration number.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_matrix_free(calibration_matrix_1);
    mpfr_matrix_free(calibration_matrix_2);
    mpfr_matrix_free(calibration_matrix_3);
    mpfr_clear(calibration_number);
    return(2);
  }
  /* Get rid of the matrices used for calibration. */
  mpfr_matrix_free(calibration_matrix_1);
  mpfr_matrix_free(calibration_matrix_2);
  mpfr_matrix_free(calibration_matrix_3);
  #ifndef NDEBUG
    fprintf(stderr, "calibration zero = ");
    mpfr_out_str(stderr, 10, 0, calibration_number, GMP_RNDN);
    fprintf(stderr,"\n");
  #endif
  /* Create a "pseudo" "sparse" Vandermonde style matrix */
  vander_style_matrix = 
            create_sparse_vandermonde_style_matrix(mpfr_points_vector,
                                                   ui_degrees_vector,
                                                   precision);
  if (vander_style_matrix == NULL)
  {
    #ifndef NDEBUG
    fprintf(stderr,
            "%s: could not compute the sparse Vandermonde matrix.\n",
            "vanderCoeffsSparseAbs");
    fprintf(stderr,
            "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    return(2); 
  }
  #ifndef NDEBUG
  mpfr_matrix_print(stderr,
                    10,
                    5,
                    GMP_RNDN,
                    vander_style_matrix);
    
  #endif
  /**
   * Compute the combinations suff.
   */
  combinations = generate_combinations(points_count, degrees_count);
  if (combinations == NULL)
  {
    #ifndef NDEBUG
    fprintf(stderr,
            "%s: could not compute the combinations.\n",
            "vanderCoeffsSparseAbs");
    fprintf(stderr,
            "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    mpfr_matrix_free(vander_style_matrix);
    return(2);    
  }
  #ifndef NDEBUG
    fprintf(stderr, "Combinations list:\n");
    ui_matrix_print(stderr, 10, combinations);
  #endif
  combinations_count = ui_matrix_get_num_rows(combinations);
  
  /* Compute the bounds */
  /* Create the eval upper bounds vector */
  eval_upper_bounds = mpfr_vector_alloc(degrees_count,
                                        precision);
  if (eval_upper_bounds == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
      "%s: could not allocate the upper bounds vector.\n",
      "vanderCoeffsSparseAbs");
      fprintf(stderr,
        "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    mpfr_matrix_free(vander_style_matrix);
    return(2); 
  }
  /* Compute the eval upper_bounds */
  for (i = 0 ; i < degrees_count ; i++)
  {
    mpfr_add(*(mpfr_vector_get_at(eval_upper_bounds, i)),
              *(mpfr_vector_get_at(func_evaluation,i)),
              max_err,
              GMP_RNDU);
  } /* End for i. */
  #ifndef NDEBUG
    fprintf(stderr, "Eval upper bounds:\n");
    mpfr_vector_print(stderr, 
                      10,
                      5,
                      GMP_RNDN,
                      eval_upper_bounds);
  #endif

  /* Create the eval lower bounds vector */
  eval_lower_bounds = mpfr_vector_alloc(degrees_count,
                                        precision);
  if (eval_lower_bounds == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not allocate the lower bounds vector.\n",
              "vanderCoeffsSparseAbs");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    mpfr_vector_free(eval_upper_bounds);
    return(2); 
  }
  /* Compute the eval lower bounds vector */
  for (i = 0 ; i < degrees_count ; i++)
  {
    mpfr_sub(*(mpfr_vector_get_at(eval_lower_bounds, i)),
              *(mpfr_vector_get_at(func_evaluation,i)),
              max_err,
              GMP_RNDU);
  } /* End for i. */
  #ifndef NDEBUG
    fprintf(stderr, "Eval lower bounds:\n");
    mpfr_vector_print(stderr, 
                      10,
                      5,
                      GMP_RNDN,
                      eval_lower_bounds);
  #endif
  /* Create the coeff_upper_bounds and coeff_lower_bounds vectors */
  coeff_upper_bounds = mpfr_vector_alloc(degrees_count, precision);
  if (coeff_upper_bounds == NULL)
  {
    #ifndef NDEBUG
    fprintf(stderr,
            "%s: could not allocate the coeff upper bounds vector.\n",
            "vanderCoeffsSparseAbs");
    fprintf(stderr,
            "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    mpfr_vector_free(eval_upper_bounds);
    mpfr_vector_free(eval_lower_bounds);
    return(2);
  }
  coeff_lower_bounds = mpfr_vector_alloc(degrees_count, precision);
  if (coeff_lower_bounds == NULL)
  {
    #ifndef NDEBUG
    fprintf(stderr,
            "%s: could not allocate the coeff lower bounds vector.\n",
            "vanderCoeffsSparseAbs");
    fprintf(stderr,
            "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    mpfr_vector_free(eval_upper_bounds);
    mpfr_vector_free(eval_lower_bounds);
    mpfr_vector_free(coeff_upper_bounds);
    return(2);
  }
  /* Set the coeff_upper_bounds elements to +infinity as a starting
     value and set coeff_lower_bounds elements to -infinity. */
  for (i = 0 ; i < degrees_count ; i++)
  {
    mpfr_set_inf(*(mpfr_vector_get_at(coeff_lower_bounds, i)),
                    -1);
    mpfr_set_inf(*(mpfr_vector_get_at(coeff_upper_bounds, i)),
                     1);
  } /* End for i */ 

  mpfr_init(scal_prod);
  mpfr_init(correction_term);
  
  for (i = 0 ; i < combinations_count ; i++)
  {
    #ifndef NDEBUG 
      fprintf(stderr, "Current combination #: %u\n", i);
    #endif
    /* Extract a pseudo Vandermonde submatrix */
    current_combination = ui_matrix_get_vector_from_row(combinations, i);
    if (current_combination == NULL)
    {
      #ifndef NDEBUG 
        fprintf(stderr,
                "%s: could not extract the current combination.\n",
                "vanderCoeffsSparseAbs");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
      #endif
      mpfr_vector_free(mpfr_points_vector);
      mpfr_vector_free(func_evaluation);
      ui_vector_free(ui_degrees_vector);
      mpfr_vector_free(abs_func_evaluation);
      mpfr_clear(calibration_number);
      mpfr_vector_free(eval_upper_bounds);
      mpfr_vector_free(eval_lower_bounds);
      mpfr_vector_free(coeff_upper_bounds);
      mpfr_vector_free(coeff_lower_bounds);
      mpfr_clear(scal_prod);
      mpfr_clear(correction_term);
      return(2);
    }                                      
    current_sub_matrix = mpfr_sub_matrix(vander_style_matrix,
                                          current_combination,
                                          precision);
    if (current_sub_matrix == NULL)
    {
      #ifndef NDEBUG 
        fprintf(stderr,
                "%s: could not extract a sub matrix from the main matrix.\n",
                "vanderCoeffsSparseAbs");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
      #endif
      mpfr_vector_free(mpfr_points_vector);
      mpfr_vector_free(func_evaluation);
      ui_vector_free(ui_degrees_vector);
      mpfr_vector_free(abs_func_evaluation);
      mpfr_clear(calibration_number);
      mpfr_vector_free(eval_upper_bounds);
      mpfr_vector_free(eval_lower_bounds);
      mpfr_vector_free(coeff_upper_bounds);
      mpfr_vector_free(coeff_lower_bounds);
      mpfr_clear(scal_prod);
      mpfr_clear(correction_term);
      ui_vector_free(current_combination);
      return(2);
    }                                      
    #ifndef NDEBUG 
      fprintf(stderr, "current_sub_matrix:\n");
      mpfr_matrix_print(stderr,
                        RADIX,
                        8,
                        GMP_RNDN,
                        current_sub_matrix);
    #endif
    
    /* Inverse the matrix */
      //fprintf(stderr,"Blib\n");
    inv_current_sub_matrix = mpfr_matrix_clone(current_sub_matrix,
                                                precision);       
    //fprintf(stderr,"Blab\n");                                         
    if (inv_current_sub_matrix == NULL)
    {
      #ifndef NDEBUG 
        fprintf(stderr,
                "%s: could not make a copy of the current sub matrix.\n",
                "vanderCoeffsSparseAbs");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
      #endif
      mpfr_vector_free(mpfr_points_vector);
      mpfr_vector_free(func_evaluation);
      ui_vector_free(ui_degrees_vector);
      mpfr_vector_free(abs_func_evaluation);
      mpfr_clear(calibration_number);
      mpfr_vector_free(eval_upper_bounds);
      mpfr_vector_free(eval_lower_bounds);
      mpfr_vector_free(coeff_upper_bounds);
      mpfr_vector_free(coeff_lower_bounds);
      mpfr_clear(scal_prod);
      mpfr_clear(correction_term);
      ui_vector_free(current_combination);
      mpfr_matrix_free(current_sub_matrix);
      return(2);
    }
    /* If we could not inverse the current sub matrix, just switch to the 
     * next one. */
    if (mpfr_matrix_inv_gauss_jordan(inv_current_sub_matrix, calibration_number))
    {
      /* Simple version without change */
      #ifndef NDEBUG
        fprintf(stderr,
                "%s: could not inverse the sub matrix:\n",
                "vanderCoeffsSparseAbs");
        mpfr_matrix_print(stderr,
                          RADIX,
                          8,
                          GMP_RNDN,
                          inv_current_sub_matrix);
        fprintf(stderr, "calibration zero = ");
        mpfr_out_str(stderr, 10, 0, calibration_number, GMP_RNDN);
        fprintf(stderr,"\n");
      #endif
      ui_vector_free(current_combination);
      mpfr_matrix_free(current_sub_matrix);
      mpfr_matrix_free(inv_current_sub_matrix);
    }
    else
    {
      /* Extract the current function evaluation sub vector. */
      sub_func_evaluation = mpfr_sub_vector(func_evaluation,
                                            current_combination,
                                            precision);
      if (sub_func_evaluation == NULL)
      {
        
        #ifndef NDEBUG      
        fprintf(stderr,
                "%s: could not extract a sub vector from the function ",
                "vanderCoeffsSparseAbs");
        fprintf(stderr,
                "evaluation vector.\n ");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
        #endif
        mpfr_vector_free(mpfr_points_vector);
        mpfr_vector_free(func_evaluation);
        ui_vector_free(ui_degrees_vector);
        mpfr_vector_free(abs_func_evaluation);
        mpfr_clear(calibration_number);
        mpfr_vector_free(eval_upper_bounds);
        mpfr_vector_free(eval_lower_bounds);
        mpfr_vector_free(coeff_upper_bounds);
        mpfr_vector_free(coeff_lower_bounds);
        mpfr_clear(scal_prod);
        mpfr_clear(correction_term);
        ui_vector_free(current_combination);
        mpfr_matrix_free(current_sub_matrix);
        mpfr_matrix_free(inv_current_sub_matrix);
        return(2);
      }
      /* Compute the min and max vector candidates */
      current_upper_bounds = mpfr_vector_alloc(degrees_count,
                                                precision);
      if (current_upper_bounds == NULL)
      {
        
        #ifndef NDEBUG
          fprintf(stderr,
                  "%s: could not allocate the current upper bounds vector.\n",
                  "vanderCoeffsSparseAbs");
          fprintf(stderr,
                  "Returning to caller with an error!\n\n");
        #endif
        mpfr_vector_free(mpfr_points_vector);
        mpfr_vector_free(func_evaluation);
        ui_vector_free(ui_degrees_vector);
        mpfr_vector_free(abs_func_evaluation);
        mpfr_clear(calibration_number);
        mpfr_vector_free(eval_upper_bounds);
        mpfr_vector_free(eval_lower_bounds);
        mpfr_vector_free(coeff_upper_bounds);
        mpfr_vector_free(coeff_lower_bounds);
        mpfr_clear(scal_prod);
        mpfr_clear(correction_term);
        ui_vector_free(current_combination);
        mpfr_matrix_free(current_sub_matrix);
        mpfr_matrix_free(inv_current_sub_matrix);
        mpfr_vector_free(sub_func_evaluation);
        return(2);
      }
      current_lower_bounds = mpfr_vector_alloc(degrees_count,
                                                precision);
      if (current_lower_bounds == NULL)
      {
        
        #ifndef NDEBUG
          fprintf(stderr,
                  "%s: could not allocate the current lower bounds vector.\n",
                  "vanderCoeffsSparseAbs");
          fprintf(stderr,
                  "Returning to caller with an error!\n\n");
        #endif
        mpfr_vector_free(mpfr_points_vector);
        mpfr_vector_free(func_evaluation);
        ui_vector_free(ui_degrees_vector);
        mpfr_vector_free(abs_func_evaluation);
        mpfr_clear(calibration_number);
        mpfr_vector_free(eval_upper_bounds);
        mpfr_vector_free(eval_lower_bounds);
        mpfr_vector_free(coeff_upper_bounds);
        mpfr_vector_free(coeff_lower_bounds);
        mpfr_clear(scal_prod);
        mpfr_clear(correction_term);
        ui_vector_free(current_combination);
        mpfr_matrix_free(current_sub_matrix);
        mpfr_matrix_free(inv_current_sub_matrix);
        mpfr_vector_free(sub_func_evaluation);
        mpfr_vector_free(current_upper_bounds);
        return(2);
      }
      for (j = 0 ; j < degrees_count ; j++)
      {
        /* Extract a line from the inverse matrix. */
        current_inv_mat_line = 
            mpfr_matrix_get_vector_from_row(inv_current_sub_matrix,
                                            j,
                                            precision);
        if (current_inv_mat_line == NULL)
        {
          #ifndef NDEBUG
          fprintf(stderr,
                  "%s: could not extract a line from the current inv matrix.\n",
                  "vanderCoeffsSparseAbs");
          fprintf(stderr,
                  "Returning to caller with an error!\n\n");
          #endif
          mpfr_vector_free(mpfr_points_vector);
          mpfr_vector_free(func_evaluation);
          ui_vector_free(ui_degrees_vector);
          mpfr_vector_free(abs_func_evaluation);
          mpfr_clear(calibration_number);
          mpfr_vector_free(eval_upper_bounds);
          mpfr_vector_free(eval_lower_bounds);
          mpfr_vector_free(coeff_upper_bounds);
          mpfr_vector_free(coeff_lower_bounds);
          mpfr_clear(scal_prod);
          mpfr_clear(correction_term);
          ui_vector_free(current_combination);
          mpfr_matrix_free(current_sub_matrix);
          mpfr_matrix_free(inv_current_sub_matrix);
          mpfr_vector_free(sub_func_evaluation);
          mpfr_vector_free(current_upper_bounds);
          mpfr_vector_free(current_lower_bounds);
          return(2);
        }
        /* Compute another vector only made of the absolute values
           of the preceeding one. */
        current_abs_inv_mat_line = mpfr_vector_get_abs(current_inv_mat_line,
                                                      precision);
        if (current_abs_inv_mat_line == NULL)
        {
          #ifndef NDEBUG
            fprintf(stderr,
                    "%s: could not compute an absolute values vector.\n",
                    "vanderCoeffsSparseAbs");
            fprintf(stderr,
                    "Returning to caller with an error!\n\n");
          #endif
          mpfr_vector_free(mpfr_points_vector);
          mpfr_vector_free(func_evaluation);
          ui_vector_free(ui_degrees_vector);
          mpfr_vector_free(abs_func_evaluation);
          mpfr_clear(calibration_number);
          mpfr_vector_free(eval_upper_bounds);
          mpfr_vector_free(eval_lower_bounds);
          mpfr_vector_free(coeff_upper_bounds);
          mpfr_vector_free(coeff_lower_bounds);
          mpfr_clear(scal_prod);
          mpfr_clear(correction_term);
          ui_vector_free(current_combination);
          mpfr_matrix_free(current_sub_matrix);
          mpfr_matrix_free(inv_current_sub_matrix);
          mpfr_vector_free(sub_func_evaluation);
          mpfr_vector_free(current_upper_bounds);
          mpfr_vector_free(current_lower_bounds);
          mpfr_vector_free(current_inv_mat_line);
          return(2);
        }
        /* Compute the lambda[j,k] * f(x[k]) sum over k. In fact, it
           is the scalar product of the the sub_func_evalution and
           current_inv_mat_line vectors. */
        mpfr_vector_scalar_product(sub_func_evaluation,
                                    current_inv_mat_line,
                                    &scal_prod,
                                    precision);
        /* Compute the |lambda[j,k]| sum over k. */
        mpfr_vector_elems_sum(current_abs_inv_mat_line,
                              &correction_term,
                              precision);
        /* Compute the min and max candidates */
        mpfr_mul(*(mpfr_vector_get_at(current_upper_bounds, j)),
                  correction_term,
                  max_err,
                  GMP_RNDN);
      
        mpfr_mul(*(mpfr_vector_get_at(current_lower_bounds, j)),
                  correction_term,
                  max_err,
                  GMP_RNDN);

        mpfr_neg(*(mpfr_vector_get_at(current_lower_bounds, j)),
                  *(mpfr_vector_get_at(current_lower_bounds, j)),
                  GMP_RNDN);

        mpfr_add(*(mpfr_vector_get_at(current_upper_bounds, j)),
                  *(mpfr_vector_get_at(current_upper_bounds, j)),
                  scal_prod,
                  GMP_RNDN);
      
        mpfr_add(*(mpfr_vector_get_at(current_lower_bounds, j)),
                  *(mpfr_vector_get_at(current_lower_bounds, j)),
                  scal_prod,
                  GMP_RNDN);
        /* Free the vectors computed inside the loops */
        mpfr_vector_free(current_inv_mat_line);
        mpfr_vector_free(current_abs_inv_mat_line);
      } /* End for j. */
      /* Update the min and max vectors */
      mpfr_vector_set_to_min(coeff_upper_bounds,
                              current_upper_bounds);
      mpfr_vector_set_to_max(coeff_lower_bounds,
                              current_lower_bounds);
      #ifndef NDEBUG
        fprintf(stderr, "Coeffs lower bounds:\n");
        mpfr_vector_print_no_dims(stderr,
                                  10,
                                  10,
                                  GMP_RNDN,
                                  coeff_lower_bounds);
        mpfr_vector_print_no_dims(stderr,
                                  10,
                                  10,
                                  GMP_RNDN,
                                  current_lower_bounds);
        fprintf(stderr, "Coeffs upper bounds:\n");
        mpfr_vector_print_no_dims(stderr,
                                  10,
                                  10,
                                  GMP_RNDN,
                                  coeff_upper_bounds);
        mpfr_vector_print_no_dims(stderr,
                                  10,
                                  10,
                                  GMP_RNDN,
                                  current_upper_bounds);
      #endif
      ui_vector_free(current_combination);
      mpfr_matrix_free(current_sub_matrix);
      mpfr_matrix_free(inv_current_sub_matrix);
      mpfr_vector_free(sub_func_evaluation);
      mpfr_vector_free(current_upper_bounds);
      mpfr_vector_free(current_lower_bounds);
    } /* End if (mpfr_matrix_inv_gauss_jordan */
  } /* End for i. */
  /* Free uneeded data. */
  mpfr_vector_free(mpfr_points_vector);
  mpfr_vector_free(func_evaluation);
  ui_vector_free(ui_degrees_vector);
  mpfr_vector_free(abs_func_evaluation);
  mpfr_clear(calibration_number);
  mpfr_matrix_free(vander_style_matrix);
  mpfr_vector_free(eval_upper_bounds);
  mpfr_vector_free(eval_lower_bounds);
  mpfr_clear(scal_prod);
  mpfr_clear(correction_term);
  ui_matrix_free(combinations);
  /* Write down into the bounds vectors */
  for (i = 0 ; i < degrees_count ; i ++)
  {
    mpfr_init(upper_bounds_array[i]);
    mpfr_set(upper_bounds_array[i],
              *(mpfr_vector_get_at(coeff_upper_bounds, i)),
              GMP_RNDN);
    mpfr_init(lower_bounds_array[i]);
    mpfr_set(lower_bounds_array[i],
              *(mpfr_vector_get_at(coeff_lower_bounds, i)),
              GMP_RNDN);
  } /* End for i. */
  mpfr_vector_free(coeff_upper_bounds);
  mpfr_vector_free(coeff_lower_bounds);
  return(0);
} /* End vanderCoeffsSparseAbs */

/**
 * @see libVanderCoeffsSparse.h#vanderCoeffsSparseRel
 */

int vanderCoeffsSparseRel(unsigned int  points_count,
                          mpfr_t*       points_array, 
                          mpfr_t*       func_evals_array,
                          unsigned int  degrees_count,
                          unsigned int* degrees_array,
                          mpfr_t        max_err,
                          mpfr_t*       lower_bounds_array,
                          mpfr_t*       upper_bounds_array)
{
  mpfr_matrix*  vander_style_matrix     = NULL;
  mpfr_matrix*  calibration_matrix_1    = NULL;
  mpfr_matrix*  calibration_matrix_2    = NULL;
  mpfr_matrix*  calibration_matrix_3    = NULL;
  mpfr_t        calibration_number;
  mpfr_vector*  mpfr_points_vector      = NULL;
  mpfr_vector*  func_evaluation         = NULL;
  mpfr_vector*  abs_func_evaluation     = NULL;
  mpfr_vector*  sub_func_evaluation     = NULL;
  mpfr_vector*  sub_abs_func_evaluation = NULL;
  mpfr_vector*  eval_upper_bounds       = NULL;
  mpfr_vector*  eval_lower_bounds       = NULL;
  ui_matrix*    combinations            = NULL;
  ui_vector*    ui_degrees_vector       = NULL;
  ui_vector*    current_combination     = NULL;
  mpfr_vector*  coeff_upper_bounds      = NULL;
  mpfr_vector*  coeff_lower_bounds      = NULL;
  mpfr_matrix*  current_sub_matrix      = NULL;
  mpfr_matrix*  inv_current_sub_matrix  = NULL;
  mpfr_vector*  current_inv_mat_line    = NULL;
  mpfr_vector*  current_abs_inv_mat_line= NULL;
  mpfr_vector*  current_upper_bounds    = NULL;
  mpfr_vector*  current_lower_bounds    = NULL;
  mpfr_t        scal_prod;
  mpfr_t        correction_term;
  mp_prec_t     precision;
  unsigned int  combinations_count      = 0;
  int           i, j;
  /* Check the parameters */
  /* Use the current default precision whenever a precision is needed. */
  precision = mpfr_get_default_prec();
  /* The number of points must be >= degree+1. */
  if (points_count < degrees_count)
  {
    /* Nothing to deallocate. */
    return(1); 
  }
  /* The input arrays must not be NULL. */
  if ((points_array == NULL)        || 
      (func_evals_array == NULL)    ||
      (degrees_array      == NULL)  ||
      (upper_bounds_array == NULL)  ||
      (lower_bounds_array == NULL))
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: one of the parameter pointers is NULL.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    /* Nothing to deallocate. */
    return(1); 
  }
  /* The max_err must be > 0 */
  if (mpfr_sgn(max_err) <= 0)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: the error is not > 0.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    /* Nothing to deallocate. */
    return(1); 
  }
  /* The number of degrees must be > 0 */
  if (degrees_count <= 0)
  {
    #ifndef NDEBUG
      fprintf(stderr,
            "%s: the degrees count is not > 0.\n",
            "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    /* Nothing to deallocate. */
    return(1); 
  }
  /* End input data check. */
  /* Create the points vector */
  mpfr_points_vector = mpfr_vector_from_array(points_array, 
                                              points_count, 
                                              precision);
  if (mpfr_points_vector == NULL)
  {
    /* Nothing to deallocate, so far. */
    return(2); 
  }
  /* Check that the points vector has no duplicates. This is still input data
     check but it can not conveniently be done before. */
  if (mpfr_vector_has_duplicates(mpfr_points_vector))
  {
    #ifndef NDEBUG
    fprintf(stderr,
            "%s: the points array has duplicates.\n",
            "vanderCoeffsSparseRel");
    fprintf(stderr,
            "Returning to caller with an error!\n\n");
  #endif
  mpfr_vector_free(mpfr_points_vector);
    return(1); 
  }
  /* Create the func_evaluation vector. */
  func_evaluation = mpfr_vector_from_array(func_evals_array, 
                                            points_count, 
                                            precision);
  if (func_evaluation == NULL)
  {
    mpfr_vector_free(mpfr_points_vector);
    return(2); 
  }
  #ifndef NDEBUG
    mpfr_vector_print(stderr, 
                      RADIX,
                      5,
                      GMP_RNDN,
                      func_evaluation);
  #endif
  /* Compute the degrees vector */
  ui_degrees_vector = ui_vector_from_array(degrees_array,
                                            degrees_count);
  if (ui_degrees_vector == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not create the ui_degrees_vector.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    return(2);
  }
  /* Check that the degrees vector has no duplicates. This is still input data
     check but it can not conveniently be done before. */
  if (ui_vector_has_duplicates(ui_degrees_vector))
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: the degrees array has duplicates.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    return(1); 
  }
  /* Compute the absolute value of the function evaluation vector. */
  abs_func_evaluation = mpfr_vector_get_abs(func_evaluation, precision);
  if (abs_func_evaluation == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not create the absolute value of the function evaluation vector.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    return(2); 
  }
  /*
   * Calibration stuff
   */
  /* Create a Vandermonde matrix for calibration */
  calibration_matrix_1 = 
      create_vandermonde_style_matrix(mpfr_vector_get_size(mpfr_points_vector) - 1,
                                      mpfr_points_vector,
                                      precision);
  if (calibration_matrix_1 == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not allocate callibration_matrix_1.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    return(2);
  }
  /* Inverse the matrix */
  calibration_matrix_2 = mpfr_matrix_inv_vander(calibration_matrix_1,
                                                  precision);
  if (calibration_matrix_2 == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not compute the inverse matrix callibration_matrix_2.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_matrix_free(calibration_matrix_1);
    return(2);
  }
  /* Compute the product of a matrix and it's inverse */
  calibration_matrix_3 = mpfr_matrix_naive_prod(calibration_matrix_1,
                                                calibration_matrix_2,
                                                precision);
  if (calibration_matrix_3 == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not compute the identity matrix callibration_matrix_3.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_matrix_free(calibration_matrix_1);
    mpfr_matrix_free(calibration_matrix_2);
    return(2);
  }
  /* Get the maximum absolute value out of the diagonal.
   * This should be "the maximum value of 0", kinf of speak. */
  mpfr_init(calibration_number);
  if (mpfr_matrix_max_abs_no_diag(calibration_matrix_3, calibration_number))
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not compute the calibration number.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_matrix_free(calibration_matrix_1);
    mpfr_matrix_free(calibration_matrix_2);
    mpfr_matrix_free(calibration_matrix_3);
    mpfr_clear(calibration_number);
    return(2);
  }
  /* Get rid of the matrices used for calibration. */
  mpfr_matrix_free(calibration_matrix_1);
  mpfr_matrix_free(calibration_matrix_2);
  mpfr_matrix_free(calibration_matrix_3);
  #ifndef NDEBUG
    fprintf(stderr, "calibration zero = ");
    mpfr_out_str(stderr, 10, 0, calibration_number, GMP_RNDN);
    fprintf(stderr,"\n");
  #endif

  /* Create a "pseudo" "sparse" Vandermonde style matrix */
  vander_style_matrix = 
            create_sparse_vandermonde_style_matrix(mpfr_points_vector,
                                                   ui_degrees_vector,
                                                   precision);
  if (vander_style_matrix == NULL)
  {
    #ifndef NDEBUG
    fprintf(stderr,
            "%s: could not compute the sparse Vandermonde matrix.\n",
            "vanderCoeffsSparseRel");
    fprintf(stderr,
            "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    return(2); 
  }
  #ifndef NDEBUG
  mpfr_matrix_print(stderr,
                    10,
                    5,
                    GMP_RNDN,
                    vander_style_matrix);
      
  #endif
  /**
   * Compute the combinations suff.
   */
  combinations = generate_combinations(points_count, degrees_count);
  if (combinations == NULL)
  {
    #ifndef NDEBUG
    fprintf(stderr,
            "%s: could not compute the combinations.\n",
            "vanderCoeffsSparseRel");
    fprintf(stderr,
            "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    mpfr_matrix_free(vander_style_matrix);
    return(2);    
  }
  #ifndef NDEBUG
    fprintf(stderr, "Combinations list:\n");
    ui_matrix_print(stderr, 10, combinations);
  #endif
  combinations_count = ui_matrix_get_num_rows(combinations);
  /* Compute the bounds */
  /* Create the upper bounds vector */
  eval_upper_bounds = mpfr_vector_alloc(points_count,
                                        precision);
  if (eval_upper_bounds == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
      "%s: could not allocate the upper bounds vector.\n",
      "vanderCoeffsSparseRel");
      fprintf(stderr,
        "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    mpfr_matrix_free(vander_style_matrix);
    return(2); 
  }
  /* Compute the upper bounds */
  for (i = 0 ; i < points_count ; i++)
  {
    /* K * |f(x)| */
    mpfr_mul(*(mpfr_vector_get_at(eval_upper_bounds, i)),
              max_err,
              *(mpfr_vector_get_at(abs_func_evaluation, i)),
              GMP_RNDN);
    mpfr_add(*(mpfr_vector_get_at(eval_upper_bounds, i)),
              *(mpfr_vector_get_at(func_evaluation,i)),
              *(mpfr_vector_get_at(eval_upper_bounds, i)),
              GMP_RNDN);
  } /* End for i. */
  #ifndef NDEBUG
    fprintf(stderr, "Eval upper bounds:\n");
    mpfr_vector_print(stderr, 
                      10,
                      5,
                      GMP_RNDN,
                      eval_upper_bounds);
  #endif
  /* Create the lower bounds vector */
  eval_lower_bounds = mpfr_vector_alloc(points_count,
                                        precision);
  if (eval_lower_bounds == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not allocate the lower bounds vector.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    ui_vector_free(ui_degrees_vector);
    mpfr_vector_free(abs_func_evaluation);
    mpfr_clear(calibration_number);
    mpfr_matrix_free(vander_style_matrix);
    mpfr_vector_free(eval_upper_bounds);
    return(2); 
  }
  /* Compute the lower bounds. */
  for (i = 0 ; i < points_count ; i++)
  {
    /* K * |f(x)| */
    mpfr_mul(*(mpfr_vector_get_at(eval_lower_bounds, i)),
              max_err,
              *(mpfr_vector_get_at(abs_func_evaluation, i)),
              GMP_RNDN);
    mpfr_sub(*(mpfr_vector_get_at(eval_lower_bounds, i)),
              *(mpfr_vector_get_at(func_evaluation, i)),
              *(mpfr_vector_get_at(eval_lower_bounds, i)),
              GMP_RNDN);
  } /* End for i. */
  #ifndef NDEBUG
    mpfr_vector_print(stderr, 
                      RADIX,
                      5,
                      GMP_RNDN,
                      eval_lower_bounds);
    mpfr_vector_print(stderr, 
                      RADIX,
                      5,
                      GMP_RNDN,
                      eval_lower_bounds);
  #endif
  /* Create the coeff bounds vectors */  
  coeff_upper_bounds = mpfr_vector_alloc(degrees_count,
                                          precision);
  if (coeff_upper_bounds == NULL)
  {
    #ifndef NDEBUG
    fprintf(stderr,
            "%s: could not allocate the coeff upper bounds vector.\n",
            "vanderCoeffsSparseRel");
    fprintf(stderr,
            "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    mpfr_matrix_free(vander_style_matrix);
    mpfr_vector_free(eval_upper_bounds);
    mpfr_vector_free(eval_lower_bounds);
    ui_matrix_free(combinations);
    return(2); 
  }
  coeff_lower_bounds = mpfr_vector_alloc(degrees_count,
                                          precision);
  if (coeff_lower_bounds == NULL)
  {
    #ifndef NDEBUG
      fprintf(stderr,
              "%s: could not allocate the coeff lower bounds vector.\n",
              "vanderCoeffsSparseRel");
      fprintf(stderr,
              "Returning to caller with an error!\n\n");
    #endif
    mpfr_vector_free(mpfr_points_vector);
    mpfr_vector_free(func_evaluation);
    mpfr_matrix_free(vander_style_matrix);
    mpfr_vector_free(eval_upper_bounds);
    mpfr_vector_free(eval_lower_bounds);
    ui_matrix_free(combinations);
    mpfr_vector_free(coeff_upper_bounds);
    return(2); 
  }
  /* Initialize the lower (resp. upper) bounds to -infinity 
     (resp. +infinity). */
  for (i = 0 ; i < degrees_count ; i++)
  {
    mpfr_set_inf(*(mpfr_vector_get_at(coeff_lower_bounds, i)), -1);
    mpfr_set_inf(*(mpfr_vector_get_at(coeff_upper_bounds, i)),  1);
  } /* End for i */ 

  mpfr_init(scal_prod);
  mpfr_init(correction_term);
  
  for (i = 0 ; i < combinations_count ; i++)
  {
    #ifndef NDEBUG 
      fprintf(stderr, "Current combination #: %u\n", i);
    #endif
    /* Extract a pseudo Vandermonde submatrix */
    current_combination = ui_matrix_get_vector_from_row(combinations, i);
    if (current_combination == NULL)
    {
      #ifndef NDEBUG 
        fprintf(stderr,
                "%s: could not extract the current combination.\n",
                "vanderCoeffsSparseRel");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
      #endif
      mpfr_vector_free(mpfr_points_vector);
      mpfr_vector_free(func_evaluation);
      ui_vector_free(ui_degrees_vector);
      mpfr_vector_free(abs_func_evaluation);
      mpfr_clear(calibration_number);
      mpfr_vector_free(eval_upper_bounds);
      mpfr_vector_free(eval_lower_bounds);
      mpfr_vector_free(coeff_upper_bounds);
      mpfr_vector_free(coeff_lower_bounds);
      mpfr_clear(scal_prod);
      mpfr_clear(correction_term);
      return(2);
    }                                      
    current_sub_matrix = mpfr_sub_matrix(vander_style_matrix,
                                          current_combination,
                                          precision);
    if (current_sub_matrix == NULL)
    {
      #ifndef NDEBUG 
        fprintf(stderr,
                "%s: could not extract a sub matrix from the main matrix.\n",
                "vanderCoeffsSparseRel");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
      #endif
      mpfr_vector_free(mpfr_points_vector);
      mpfr_vector_free(func_evaluation);
      ui_vector_free(ui_degrees_vector);
      mpfr_vector_free(abs_func_evaluation);
      mpfr_clear(calibration_number);
      mpfr_vector_free(eval_upper_bounds);
      mpfr_vector_free(eval_lower_bounds);
      mpfr_vector_free(coeff_upper_bounds);
      mpfr_vector_free(coeff_lower_bounds);
      mpfr_clear(scal_prod);
      mpfr_clear(correction_term);
      ui_vector_free(current_combination);
      return(2);
    }                                      
    #ifndef NDEBUG 
      fprintf(stderr, "current_sub_matrix:\n");
      mpfr_matrix_print(stderr,
                        RADIX,
                        8,
                        GMP_RNDN,
                        current_sub_matrix);
    #endif
    
    /* Inverse the matrix */

      //fprintf(stderr,"Blib\n");

    inv_current_sub_matrix = mpfr_matrix_clone(current_sub_matrix,
                                                precision);                                                

    //fprintf(stderr,"Blab\n");
    if (inv_current_sub_matrix == NULL)
    {
      #ifndef NDEBUG 
        fprintf(stderr,
                "%s: could not make a copy of the current sub matrix.\n",
                "vanderCoeffsSparseRel");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
      #endif
      mpfr_vector_free(mpfr_points_vector);
      mpfr_vector_free(func_evaluation);
      ui_vector_free(ui_degrees_vector);
      mpfr_vector_free(abs_func_evaluation);
      mpfr_clear(calibration_number);
      mpfr_vector_free(eval_upper_bounds);
      mpfr_vector_free(eval_lower_bounds);
      mpfr_vector_free(coeff_upper_bounds);
      mpfr_vector_free(coeff_lower_bounds);
      mpfr_clear(scal_prod);
      mpfr_clear(correction_term);
      ui_vector_free(current_combination);
      mpfr_matrix_free(current_sub_matrix);
      return(2);
    }
    /* If we could not inverse the current sub matrix, just switch to the 
     * next one. */
    if (mpfr_matrix_inv_gauss_jordan(inv_current_sub_matrix, calibration_number))
    {
      //fprintf(stderr,"Blub\n");
      /* Simple version without change */
      #ifndef NDEBUG
        fprintf(stderr,
                "%s: could not inverse the sub matrix:\n",
                "vanderCoeffsSparseRel");
        mpfr_matrix_print(stderr,
                          RADIX,
                          8,
                          GMP_RNDN,
                          inv_current_sub_matrix);
        fprintf(stderr, "calibration zero = ");
        mpfr_out_str(stderr, 10, 0, calibration_number, GMP_RNDN);
        fprintf(stderr,"\n");
      #endif
      ui_vector_free(current_combination);
      mpfr_matrix_free(current_sub_matrix);
      mpfr_matrix_free(inv_current_sub_matrix);

      //fprintf(stderr,"Blieb\n");
      
    }
    else
    {
      //fprintf(stderr,"Bleb\n");
      /* Extract the current function evaluation sub vector. */
      sub_func_evaluation = mpfr_sub_vector(func_evaluation,
                                            current_combination,
                                            precision);
      if (sub_func_evaluation == NULL)
      {
        
        #ifndef NDEBUG      
        fprintf(stderr,
                "%s: could not extract a sub vector from the function ",
                "vanderCoeffsSparseRel");
        fprintf(stderr,
                "evaluation vector.\n ");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
        #endif
        mpfr_vector_free(mpfr_points_vector);
        mpfr_vector_free(func_evaluation);
        ui_vector_free(ui_degrees_vector);
        mpfr_vector_free(abs_func_evaluation);
        mpfr_clear(calibration_number);
        mpfr_vector_free(eval_upper_bounds);
        mpfr_vector_free(eval_lower_bounds);
        mpfr_vector_free(coeff_upper_bounds);
        mpfr_vector_free(coeff_lower_bounds);
        mpfr_clear(scal_prod);
        mpfr_clear(correction_term);
        ui_vector_free(current_combination);
        mpfr_matrix_free(current_sub_matrix);
        mpfr_matrix_free(inv_current_sub_matrix);
        return(2);
      }
      /* Extract the current function evaluation sub vector. */
      sub_abs_func_evaluation = mpfr_sub_vector(abs_func_evaluation,
                                                current_combination,
                                                precision);
      if (sub_abs_func_evaluation == NULL)
      {
        
        #ifndef NDEBUG      
        fprintf(stderr,
                "%s: could not extract a sub vector from the function ",
                "vanderCoeffsSparseRel");
        fprintf(stderr,
                "evaluation absolute values vector.\n ");
        fprintf(stderr,
                "Returning to caller with an error!\n\n");
        #endif
        mpfr_vector_free(mpfr_points_vector);
        mpfr_vector_free(func_evaluation);
        ui_vector_free(ui_degrees_vector);
        mpfr_vector_free(abs_func_evaluation);
        mpfr_clear(calibration_number);
        mpfr_vector_free(eval_upper_bounds);
        mpfr_vector_free(eval_lower_bounds);
        mpfr_vector_free(coeff_upper_bounds);
        mpfr_vector_free(coeff_lower_bounds);
        mpfr_clear(scal_prod);
        mpfr_clear(correction_term);
        ui_vector_free(current_combination);
        mpfr_matrix_free(current_sub_matrix);
        mpfr_matrix_free(inv_current_sub_matrix);
        mpfr_vector_free(abs_func_evaluation);
        return(2);
      }
      /* Compute the min and max vector candidates */
      current_upper_bounds = mpfr_vector_alloc(degrees_count,
                                                precision);
      if (current_upper_bounds == NULL)
      {
        
        #ifndef NDEBUG
          fprintf(stderr,
                  "%s: could not allocate the current upper bounds vector.\n",
                  "vanderCoeffsSparseRel");
          fprintf(stderr,
                  "Returning to caller with an error!\n\n");
        #endif
        mpfr_vector_free(mpfr_points_vector);
        mpfr_vector_free(func_evaluation);
        ui_vector_free(ui_degrees_vector);
        mpfr_vector_free(abs_func_evaluation);
        mpfr_clear(calibration_number);
        mpfr_vector_free(eval_upper_bounds);
        mpfr_vector_free(eval_lower_bounds);
        mpfr_vector_free(coeff_upper_bounds);
        mpfr_vector_free(coeff_lower_bounds);
        mpfr_clear(scal_prod);
        mpfr_clear(correction_term);
        ui_vector_free(current_combination);
        mpfr_matrix_free(current_sub_matrix);
        mpfr_matrix_free(inv_current_sub_matrix);
        mpfr_vector_free(sub_func_evaluation);
        return(2);
      }
      current_lower_bounds = mpfr_vector_alloc(degrees_count,
                                                precision);
      if (current_lower_bounds == NULL)
      {
        
        #ifndef NDEBUG
          fprintf(stderr,
                  "%s: could not allocate the current lower bounds vector.\n",
                  "vanderCoeffsSparseRel");
          fprintf(stderr,
                  "Returning to caller with an error!\n\n");
        #endif
        mpfr_vector_free(mpfr_points_vector);
        mpfr_vector_free(func_evaluation);
        ui_vector_free(ui_degrees_vector);
        mpfr_vector_free(abs_func_evaluation);
        mpfr_clear(calibration_number);
        mpfr_vector_free(eval_upper_bounds);
        mpfr_vector_free(eval_lower_bounds);
        mpfr_vector_free(coeff_upper_bounds);
        mpfr_vector_free(coeff_lower_bounds);
        mpfr_clear(scal_prod);
        mpfr_clear(correction_term);
        ui_vector_free(current_combination);
        mpfr_matrix_free(current_sub_matrix);
        mpfr_matrix_free(inv_current_sub_matrix);
        mpfr_vector_free(sub_func_evaluation);
        mpfr_vector_free(current_upper_bounds);
        return(2);
      }
      for (j = 0 ; j < degrees_count ; j++)
      {
        /* Extract a line from the inverse matrix. */
        current_inv_mat_line = 
            mpfr_matrix_get_vector_from_row(inv_current_sub_matrix,
                                            j,
                                            precision);
        if (current_inv_mat_line == NULL)
        {
          #ifndef NDEBUG
          fprintf(stderr,
                  "%s: could not extract a line from the current inv matrix.\n",
                  "vanderCoeffsSparseRel");
          fprintf(stderr,
                  "Returning to caller with an error!\n\n");
          #endif
          mpfr_vector_free(mpfr_points_vector);
          mpfr_vector_free(func_evaluation);
          ui_vector_free(ui_degrees_vector);
          mpfr_vector_free(abs_func_evaluation);
          mpfr_clear(calibration_number);
          mpfr_vector_free(eval_upper_bounds);
          mpfr_vector_free(eval_lower_bounds);
          mpfr_vector_free(coeff_upper_bounds);
          mpfr_vector_free(coeff_lower_bounds);
          mpfr_clear(scal_prod);
          mpfr_clear(correction_term);
          ui_vector_free(current_combination);
          mpfr_matrix_free(current_sub_matrix);
          mpfr_matrix_free(inv_current_sub_matrix);
          mpfr_vector_free(sub_func_evaluation);
          mpfr_vector_free(current_upper_bounds);
          mpfr_vector_free(current_lower_bounds);
          return(2);
        }
        /* Compute another vector only made of the absolute values
           of the preceeding one. */
        current_abs_inv_mat_line = mpfr_vector_get_abs(current_inv_mat_line,
                                                      precision);
        if (current_abs_inv_mat_line == NULL)
        {
          #ifndef NDEBUG
            fprintf(stderr,
                    "%s: could not compute an absolute values vector.\n",
                    "vanderCoeffsSparseRel");
            fprintf(stderr,
                    "Returning to caller with an error!\n\n");
          #endif
          mpfr_vector_free(mpfr_points_vector);
          mpfr_vector_free(func_evaluation);
          ui_vector_free(ui_degrees_vector);
          mpfr_vector_free(abs_func_evaluation);
          mpfr_clear(calibration_number);
          mpfr_vector_free(eval_upper_bounds);
          mpfr_vector_free(eval_lower_bounds);
          mpfr_vector_free(coeff_upper_bounds);
          mpfr_vector_free(coeff_lower_bounds);
          mpfr_clear(scal_prod);
          mpfr_clear(correction_term);
          ui_vector_free(current_combination);
          mpfr_matrix_free(current_sub_matrix);
          mpfr_matrix_free(inv_current_sub_matrix);
          mpfr_vector_free(sub_func_evaluation);
          mpfr_vector_free(current_upper_bounds);
          mpfr_vector_free(current_lower_bounds);
          mpfr_vector_free(current_inv_mat_line);
          return(2);
        }
        /* Compute the min and max candidates */
        /* Compute the lambda[j,k] * f(x[k]) sum over k. In fact, it
           is the scalar product of the the sub_func_evalution and
           current_inv_mat_line vectors. */
        mpfr_vector_scalar_product(sub_func_evaluation,
                                    current_inv_mat_line,
                                    &scal_prod,
                                    precision);
        /* Compute the |lambda[j,k]| * |f(x[k])| sum over k. In fact, is
           is the scalar product of the sub_abs_func_evaluation and
           the current current_abs_inv_mat_line vectors.*/
        mpfr_vector_scalar_product(sub_abs_func_evaluation,
                                    current_abs_inv_mat_line,
                                    &correction_term,
                                    precision);
        /* Compute the min candidate */
        mpfr_mul(*(mpfr_vector_get_at(current_lower_bounds, j)),
                  correction_term,
                  max_err,
                  GMP_RNDN);

        mpfr_mul(*(mpfr_vector_get_at(current_upper_bounds, j)),
                  correction_term,
                  max_err,
                  GMP_RNDN);
  
        /* For the min candidate, take the opposite of the correction term. */          
        mpfr_neg(*(mpfr_vector_get_at(current_lower_bounds, j)),
                  *(mpfr_vector_get_at(current_lower_bounds, j)),
                  GMP_RNDN);

        mpfr_add(*(mpfr_vector_get_at(current_upper_bounds, j)),
                  *(mpfr_vector_get_at(current_upper_bounds, j)),
                  scal_prod,
                  GMP_RNDN);        
  
        mpfr_add(*(mpfr_vector_get_at(current_lower_bounds, j)),
                  *(mpfr_vector_get_at(current_lower_bounds, j)),
                  scal_prod,
                  GMP_RNDN);
        mpfr_vector_free(current_inv_mat_line);
        mpfr_vector_free(current_abs_inv_mat_line);
      } /* End for j. */
      /* Update the min and max vectors */
      mpfr_vector_set_to_min(coeff_upper_bounds,
                              current_upper_bounds);
      mpfr_vector_set_to_max(coeff_lower_bounds,
                              current_lower_bounds);
      #ifndef NDEBUG
        fprintf(stderr, "Coeffs lower bounds:\n");
        mpfr_vector_print_no_dims(stderr,
                                  10,
                                  10,
                                  GMP_RNDN,
                                  coeff_lower_bounds);
        mpfr_vector_print_no_dims(stderr,
                                  10,
                                  10,
                                  GMP_RNDN,
                                  current_lower_bounds);
        fprintf(stderr, "Coeffs upper bounds:\n");
        mpfr_vector_print_no_dims(stderr,
                                  10,
                                  10,
                                  GMP_RNDN,
                                  coeff_upper_bounds);
        mpfr_vector_print_no_dims(stderr,
                                  10,
                                  10,
                                  GMP_RNDN,
                                  current_upper_bounds);
      #endif
      ui_vector_free(current_combination);
      mpfr_matrix_free(current_sub_matrix);
      mpfr_matrix_free(inv_current_sub_matrix);
      mpfr_vector_free(sub_func_evaluation);
      mpfr_vector_free(current_upper_bounds);
      mpfr_vector_free(current_lower_bounds);
    } /* End if (mpfr_matrix_inv_gauss_jordan */

    //fprintf(stderr,"Blueb\n");

  } /* End for i. */

  //fprintf(stderr,"Bleib\n");

  /* Free uneeded data. */
  mpfr_vector_free(mpfr_points_vector);
  mpfr_vector_free(func_evaluation);
  ui_vector_free(ui_degrees_vector);
  mpfr_vector_free(abs_func_evaluation);
  mpfr_clear(calibration_number);

  //fprintf(stderr,"Blueng\n");

  mpfr_matrix_free(vander_style_matrix);
  mpfr_vector_free(eval_upper_bounds);
  mpfr_vector_free(eval_lower_bounds);
  mpfr_clear(scal_prod);
  mpfr_clear(correction_term);
  ui_matrix_free(combinations);
/* Write down into the bounds vectors */
  for (i = 0 ; i < degrees_count ; i ++)
  {
    
    //fprintf(stderr,"Bloing\n");

    mpfr_init(upper_bounds_array[i]);
    mpfr_set(upper_bounds_array[i],
              *(mpfr_vector_get_at(coeff_upper_bounds, i)),
              GMP_RNDN);
    mpfr_init(lower_bounds_array[i]);
    mpfr_set(lower_bounds_array[i],
              *(mpfr_vector_get_at(coeff_lower_bounds, i)),
              GMP_RNDN);
  } /* End for i. */
  mpfr_vector_free(coeff_upper_bounds);
  mpfr_vector_free(coeff_lower_bounds);

  //fprintf(stderr,"Bluerb\n");

  return(0);
} /* End vanderCoeffsSparseRel */

