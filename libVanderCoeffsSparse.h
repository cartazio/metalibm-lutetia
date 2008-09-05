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

/* libVanderCoeffsSparse.h */
/******************************************************************************/
/**
 * Name & purpose
 * Author:
 */
 
/*
 * Prerequisites for use:
 *   The user must include
 */

#ifndef LIB_VANDER_COEFFS_SPARSE_LIB_h
#define LIB_VANDER_COEFFS_SPARSE_LIB_h

/**
 * Computes the lower bounds and upper bounds for the coefficients of a
 * polynomial approximation for an absolute error.
 * The computed values are unchecked for consistency. If any of the upper
 * bounds is smaller than the corresponding lower bound, the results is
 * inconsistent. It is the responsability of the caller to make this check
 * and to take the most sensible action (quite often, calling the function
 * again with a larger error helps).
 * The precision used to carry out MPFR internal computation is the current 
 * precision in the caller at call time.
 * @param points_cout ->  the size of the array holding the "interpolation" points;
 *                        set by the caller ;
 *                        read-only by the callee;
 * @param points_array -> the array holding the "interpolation" points ; 
 *                        must hold at least points_count elements ;
 *                        must not hold duplicates (singular matrix failure); 
 *                        allocated and assigned by the caller;
 *                        read-only by the callee;
 * @param func_evals_array -> the array holding the values of the function to
 *                            to approximate evaluated at the points given
 *                            in points_array;
 *                            must hold at least points_count elements ;
 *                            indices must match in points_array and func_evals_array;
 *                            allocated and assigned by the caller; 
 *                            read-only by the callee;
 * @param degrees_count ->  the size of the array holding the degrees of the monomials
 *                          the approximation polynomial is made of;
 *                          set by the caller;
 *                          read-only by the callee;
 * @param degrees_array ->  the array holding the degrees of the monomials the 
 *                          approximation polynomial is madeof;
 *                          must hold at least degrees_count elements;
 *                          must not hold duplicates (singular matrix failure); 
 *                          allocated and assigned by the caller;
 *                          read-only by the callee;
 * @param max_err ->  the maximum (absolute) error of the approximation polynomial;
 *                    set by the caller;
 *                    read-only by the callee;
 * @param lower_bounds -> the array holding the lower bounds for the coefficients
 *                        of the monomials the approximation polynomial is made of;
 *                        must hold at least degrees_count elements;
 *                        it's indices match those of degrees_array;
 *                        allocated (but not MPFR_inited) by the caller;
 *                        MPFR_inited and assigned by the callee if return
 *                        value is 0;
 *                        untouched by the callee if the return value != 0;
 * @param upper_bounds -> the array holding the upper bounds for the coefficients
 *                        of the monomials the approximation polynomial is made of;
 *                        must hold at least degrees_count elements;
 *                        it's indices match those of degrees_array;
 *                        allocated (but not MPFR_inited) by the caller;
 *                        MPFR_inited and assigned by the callee if return
 *                        value is 0;
 *                        untouched by the callee if the return value != 0;
 * @return -> 0 of everything is OK;
 *            != 0 is something goes wrong; in this case lower_bounds and 
 *            upper_bounds are untouched (specifically: not MPFR_inited).
 */ 
int 
vanderCoeffsSparseAbs(unsigned int  points_count,
                      mpfr_t*       points_array, 
                      mpfr_t*       func_evals_array,
                      unsigned int  degrees_count,
                      unsigned int* degrees_array,
                      mpfr_t        max_err,
                      mpfr_t*       lower_bounds,
                      mpfr_t*       upper_bounds);
                    
/**
 * Computes the lower bounds and upper bounds for the coefficients of a
 * polynomial approximation for a relative error.
 * The computed values are unchecked for consistency. If any of the upper
 * bounds is smaller than the corresponding lower bound, the results is
 * inconsistent. It is the responsability of the caller to make this check
 * and to take the most sensible action (quite often, calling the function
 * again with a larger error helps).
 * The precision used to carry out MPFR internal computation is the current 
 * precision in the caller at call time.
 * 
 * @param points_cout ->  the size of the array holding the "interpolation" points;
 *                        set by the caller ;
 *                        read-only by the callee;
 * @param points_array -> the array holding the "interpolation" points ; 
 *                        must hold at least points_count elements ; 
 *                        must not hold duplicates (singular matrix failure); 
 *                        allocated and assigned by the caller;
 *                        read-only by the callee;
 * @param func_evals_array -> the array holding the values of the function to
 *                            to approximate evaluated at the points given
 *                            in points_array;
 *                            must hold at least points_count elements ;
 *                            indices must match in points_array and func_evals_array;
 *                            allocated and assigned by the caller; 
 *                            read-only by the callee;
 * @param degrees_count ->  the size of the array holding the degrees of the monomials
 *                          the approximation polynomial is made of;
 *                          set by the caller;
 *                          read-only by the callee;
 * @param degrees_array ->  the array holding the degrees of the monomials the 
 *                          approximation polynomial is madeof;
 *                          must hold at least degrees_count elements;
 *                          must not hold duplicates (singular matrix failure); 
 *                          allocated and assigned by the caller;
 *                          read-only by the callee;
 * @param max_err ->  the (relative) maximum error of the approximation polynomial;
 *                    set by the caller;
 *                    read-only by the callee;
 * @param lower_bounds -> the array holding the lower bounds for the coefficients
 *                        of the monomials the approximation polynomial is made of;
 *                        must hold at least degrees_count elements;
 *                        it's indices match those of degrees_array;
 *                        allocated (but not MPFR_inited) by the caller;
 *                        MPFR_inited and assigned by the callee if return
 *                        value is 0;
 *                        untouched by the callee if the return value != 0;
 * @param upper_bounds -> the array holding the upper bounds for the coefficients
 *                        of the monomials the approximation polynomial is made of;
 *                        must hold at least degrees_count elements;
 *                        it's indices match those of degrees_array;
 *                        allocated (but not MPFR_inited) by the caller;
 *                        MPFR_inited and assigned by the callee if return
 *                        value is 0;
 *                        untouched by the callee if the return value != 0;
 * @return -> 0 of everything is OK;
 *            != 0 is something goes wrong; in this case lower_bounds and 
 *            upper_bounds are untouched (specifically: not MPFR_inited).
 */ 
int 
vanderCoeffsSparseRel(unsigned int  points_count,
                      mpfr_t*       points_array, 
                      mpfr_t*       func_evals_array,
                      unsigned int  degrees_count,
                      unsigned int* degrees_array,
                      mpfr_t        max_err,
                      mpfr_t*       lower_bounds,
                      mpfr_t*       upper_bounds);


int 
vanderCoeffsSparseAbs(unsigned int  points_count,
                      mpfr_t*       points_array, 
                      mpfr_t*       func_evals_array,
                      unsigned int  degrees_count,
                      unsigned int* degrees_array,
                      mpfr_t        max_err,
                      mpfr_t*       lower_bounds,
                      mpfr_t*       upper_bounds);

                    

#endif
  
