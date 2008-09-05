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

/* utils.h */
/******************************************************************************/

/** @file
 * Utilities for the programme computing an interval for the coefficients of 
 * an interpolation polynomial.
 * @author ST 
 * @date 2006-07-21
 */
 
/*
 * Prerequisites for use:
 *   The user must include:
 *      mpfr-vector.h
 *      mpfr-matrix.h
 *      ui-vector.h
 *      ui-matrix.h
 *      gmp.h
 */

#ifndef VANDER_COEFFS_UTILS_h
#define VANDER_COEFFS_UTILS_h  

/**
 * Create a "Vandermonde like" matrix from a vector of MPFR "reals" and
 * the given degree.
 * Each line of the matrix is made of the powers of the corresponding element
 * in the input vector.
 * It is not a Vandermonde matrix since (degree+1) may differ from the number
 * of elements of the vector. Hence the matrix is not square.
 *
 * @param degree    - the maximum power in each line;
 * @param vect      - the "generator" vector;
 * @param precision - the precision used in MPFR to compute the powers.
 * @return a pointer to the created matrix or NULL if the creation failed for some reason. 
 */

mpfr_matrix*  
create_vandermonde_style_matrix(unsigned int degree, 
                                mpfr_vector* vect,
                                unsigned int precision);

/**
 * Create a "Vandermonde like" matrix from a vector of MPFR floats and
 * a vector of exponents.
 * Each line of the matrix is made of the powers, according of the expoents
 * vector, of the corresponding element in the points vector.
 * This is not a Vandermonde matrix at all: 
 * - the expoents list differs for 0, 1, 2...;
 * - the number of points may differ (must be greater or equal to) from the 
 *   number of powers, hence the matrix is not square.
 *
 * @param degree    - the maximum power in each line;
 * @param vect      - the "generator" vector;
 * @param precision - the precision used in MPFR to compute the powers.
 * @return a pointer to the created matrix or NULL if the creation failed for some reason. 
 */

mpfr_matrix*  
create_sparse_vandermonde_style_matrix( mpfr_vector* points,
                                        ui_vector*   exponents,
                                        unsigned int precision);

/**
 * Build the master polynomial of a Vandermonde matrix.
 *
 * @param mat - the matrix;
 * @param precision - the precision, in bits, of the computations.
 * @return a vector (or NULL if something went wrong).
 */

mpfr_vector*
mpfr_matrix_build_master_polynomial(mpfr_matrix* mat, unsigned precision);

/**
 * Create a vector whose elements are the 0, 1,.\ .\ .\ , <i>max_power</i> 
 * powers
 * of <i>x</i>.
 *
 * @param x         - the number to compute the powers of;
 * @param max_power - the maximum power computed;
 * @param precision - the precision, in bits, of the vector elements.
 * @return a vector holding the first powers of <i>x</i>. The index
 * in the vector matches the power (e.g. element at index 5 is x^5).
 */
mpfr_vector*
mpfr_first_powers_of(mpfr_t x, unsigned int max_power, unsigned int precision);

/**
 * Generate a matrix of unsigned ints where each row is an element of the 
 * combinations of k elements from the {0, 1,.\ .\ .\ ,n-1} set.
 *
 * @param n - the number of element of the base set: {0, 1,... , n-1}
 *            (e.g. {0, 1, 2, 3, 4, 5}).
 * @param k - the number of elements of a combination (e.g. if k == 4,
 *            {2, 3, 4, 5} is a valid combination).
 * @return a matrix holding all the possible combinations (or NULL
 * if something went wrong).
 */
ui_matrix*
generate_combinations(unsigned int n, unsigned int k);

/**
 * Compute the "lexicographic" successor a given vector 
 * in the list of the combinations k elements taken
 * from the {0, 1,.\ .\ .\ , n-1) set}.
 * Combinations are computed as ui_vectors. The successor is
 * computed "in place", overwriting the predecessor.
 * @param vect - the in-out vector holding the predecessor
 *               (upon calling) and the successor (when
 *               the function returns successfully);
 * @param n    - the size and bound of the initial "collection"
 *               (e.g. if n == 6 then the set is {0, 1, 2, 3, 4, 5});
 * @param k    - the number of element of the combination
 *               (e.g. if k == 4, {2, 3, 4, 5} is a valid combination).
 * @return the vect parameter (or NULL if a successor can not be computed).
 *
 */
ui_vector*
successor(ui_vector* vect, unsigned int n, unsigned int k);

/**
 * Generate a matrix from another (parameter mat) picking up
 * rows whose index in given by the elements of a ui_vector
 * (parameter vect).
 *
 * @param mat       - the "source" matrix;
 * @param vect      - the vector holding the indexes of the rows to pickup;
 * @param precision - the precision, in bits, of the elements of the returned
 *                    matrix.
 * @return a matrix made of the rows of mat whose index is given in
 * vect (or NULL if something went wrong).
 */

mpfr_matrix*
mpfr_sub_matrix(mpfr_matrix* mat, 
                ui_vector* vect,
                unsigned precision);

/**
 * Generate a vector from another (parameter mpfr_vect) picking up
 * the elements whose index in given by the elements of a ui_vector
 * (parameter ui_vect).
 *
 * @param mpfr_vect - the "source" vector;
 * @param ui_vect   - the vector holding the indexes of the rows to pickup;
 * @param precision - the precision, in bits, of the computed elements.
 * @return a vector made of the elements of <i>mpfr_vect</i> whose index 
 * is given in <i>ui_vect</i> (or NULL if something went wrong).
 */
mpfr_vector*
mpfr_sub_vector(mpfr_vector* mpfr_vect, 
                ui_vector* ui_vect,
                unsigned precision);

mpfr_vector*
mpfr_first_powers_of(mpfr_t x, unsigned int max_power, unsigned int precision);

/**
 * Create a vector whose elements are arbitrary integer powers of <i>x</i>.
 * 
 * The exponents list is given in the <i>powers_exponents</i> vector.
 *
 * @param x                - the number to compute the powers of;
 * @param powers_exponents - the list of the exponents;
 * @param precision        - the precision, in bits, of the vector elements.
 * @return a vector holding the arbitrary powers of <i>x</i>. 
 */
mpfr_vector*
mpfr_arbitrary_int_powers_of(mpfr_t x,
                              ui_vector* powers_exponents, 
                              unsigned int precision);


/**
 * Generate the inverse matrix of a Vandermonde matrix.
 *
 * @param mat       - a pointer to the original matrix;
 * @param precision - the precision, in bits, of the inverse
 *                    matrix elements.
 * @return a pointer to the inverse matrix (or NULL if
 * something went wrong).
 */
mpfr_matrix*
mpfr_matrix_inv_vander(mpfr_matrix* mat, unsigned precision);

#endif /* VANDER_COEFFS_UTILS_h */

