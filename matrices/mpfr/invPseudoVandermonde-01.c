/* invPseudoVandermonde.c */
/******************************************************************************/
/**
 * Name & purpose                                                           
 * Author:                                                                  
 *
 */

/* includes of system headers */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>
#include <mpfr.h>

/* includes of project headers */
#include "../vectors/mpfr/mpfr-vector.h"
#include "../vectors/mpq/mpq-vector.h"
#include "../vectors/si/si-vector.h"
#include "../vectors/ui/ui-vector.h"
#include "../matrices/mpfr/mpfr-matrix.h"
#include "../matrices/mpq/mpq-matrix.h"
#include "../matrices/ui/ui-matrix.h"
#include "../misc-utils/std-exit-errors.h"
#include "../misc-utils/utilities.h"
#include "../proj-utils/utils.h"

/* includes of local headers */

/* Types, constants and macros definitions */

#define ALL_NUMBERS_BASE 10
#define INDEX_FIRST_EXPONENT 5
#define DRIFT_SLICE 100
/* Global variables */

/* Main function */
int main(int argc, char** argv)
{
  mpq_t         rat_lower_bound;
  mpq_t         rat_upper_bound;
  mpq_t         rat_drift_increment;
  mpq_t         rat_points_min_spacing;
  mpq_t         rat_temp_1;
  mpq_t         rat_temp_2;
  si_vector*    powers_exponents;
  mpq_matrix*   pseudo_vander_matrix;
  mpq_matrix*   inverted_pseudo_vander_matrix;
  mpq_matrix*   identity_matrix;
  mpq_vector*   powers_values;
  mpq_vector*   points;
  int           num_exponents;
  unsigned int  i;
  int           exponent;
  char*         endptr[1];
  si_vector*    drift_directions;
  FILE*         inv_matrix_file;
  FILE*         points_file;

  if (argc < 7)
  {
    fprintf(stderr,
            "\n\nUsage: %s invMatrixFile pointsFile lowerBound upperBound powersExponentsList\n\n",
            argv[0]);
    fprintf(stderr,
            "       -invMatrixFile: the name of the file holding the matrix;\n");
    fprintf(stderr,
            "       -pointsFile: the name of the file holding the points");
    fprintf(stderr,
            "        used to compute the matrix;\n");
    fprintf(stderr,
            "       -lowerBound: the lower bound of the interval as a string\n");
    fprintf(stderr,
            "        convertible into an \"mpq_t rational\";\n");
    fprintf (stderr,
            "       -upperBound: the upper bound of the interval as a string\n");
    fprintf(stderr,
            "        convertible into an \"mpq_t rational\";\n");
    fprintf(stderr,
            "       -powersExponentsList: the list of the exponents of the monomials\n");
    fprintf(stderr,
            "        of the approximation polynomial as string convertible into\n");
    fprintf(stderr,
            "        integers (at least two are needed)\n\n");
    fprintf(stderr,
            "The base for all numbers is 10\n\n");
    return(EX_USAGE);
  } /* End if argc... */
  /* Read the command line parameters */
  mpq_init(rat_lower_bound);
  if (mpq_set_str(rat_lower_bound, argv[3], ALL_NUMBERS_BASE))
  {
    fprintf(stderr,
            "%s: can not convert the lower bound \"%s\" into an mpq_t\n",
            argv[0],
            argv[3]);
    fprintf(stderr,
            "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    return(EX_USAGE);
  } /* End read rat_lower_bound. */
  mpq_init(rat_upper_bound);
  if (mpq_set_str(rat_upper_bound, argv[4], ALL_NUMBERS_BASE))
  {
    fprintf(stderr,
            "%s: can not convert the upper bound \"%s\" into an mpq_t\n",
            argv[0],
            argv[4]);
    fprintf(stderr,
            "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    return(EX_USAGE);
  } /* End read rat_upper_bound. */
  /* Work out the exponents stuff */
  num_exponents = argc - INDEX_FIRST_EXPONENT;

  powers_exponents = si_vector_alloc(num_exponents);
  if (powers_exponents == NULL)
  {
    fprintf(stderr,
            "%s: can not allocate the powers exponents vector.\n",
            argv[0]);
    fprintf(stderr,
            "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    return(EX_SOFTWARE);
  }
  /* Read the exponents */
  for (i = 0 ; i < num_exponents ; i++)
  {
    exponent = strtol(argv[INDEX_FIRST_EXPONENT + i], endptr, 10);
    /* For a completely safe conversion *endptr must point to a 0 char value
      (not the '0' char but '\0' char value!). */
    if (**endptr != '\0')
    {
      fprintf(stderr,
              "\n\n%s: Could not safely convert the \"%s\" value into ",
              argv[0],
              argv[INDEX_FIRST_EXPONENT + i]);
      fprintf(stderr,
              "an int.\n");
      fprintf(stderr, "Aborting the program!\n\n");
      mpq_clear(rat_lower_bound);
      mpq_clear(rat_upper_bound);
      si_vector_free(powers_exponents);
      return(EX_USAGE);
    } /* End individual exponent reading. */

    si_vector_set_at(powers_exponents, i, exponent);
  }
  /* End of command line parameters reading. */
  /* Handle and check the command line parameters */
  if ((inv_matrix_file = fopen(argv[1], "w")) == NULL)
  {
    fprintf(stderr,
            "\n\n%s: Could not open the \"%s\" file\n.", argv[0], argv[1]);
    fprintf(stderr, "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    si_vector_free(powers_exponents);
    return(EX_USAGE);
  }
  if ((points_file = fopen(argv[2], "w")) == NULL)
  {
    fprintf(stderr,
            "\n\n%s: Could not open the \"%s\" file\n.", argv[0], argv[2]);
    fprintf(stderr, "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    si_vector_free(powers_exponents);
    fclose(inv_matrix_file);
    unlink(argv[1]);
    return(EX_USAGE);
  }
  /* Canonicalize the bounds and check that the upper_bound 
     is greater than the lower_bound. */
  mpq_canonicalize(rat_upper_bound);
  mpq_canonicalize(rat_lower_bound);
  if (mpq_cmp(rat_lower_bound, rat_upper_bound) >= 0)
  {
    fprintf(stderr,
            "\n\n%s: The lower bound (",
            argv[0]);
    mpq_out_str(stderr, 10, rat_lower_bound);
    fprintf(stderr,
            ") is greater or equal to the upper bound (");
    mpq_out_str(stderr, 10, rat_upper_bound);
    fprintf(stderr,
            ").\n");
    fprintf(stderr, "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    si_vector_free(powers_exponents);
    fclose(inv_matrix_file);
    unlink(argv[1]);
    fclose(points_file);
    unlink(argv[2]);
    return(EX_DATAERR);
  }
  /* Check that there are no duplicates in the powers: that would
     automaticaly lead to a singular matrix */
  if (si_vector_has_duplicates(powers_exponents))
  {
    fprintf(stderr,
            "\n\n%s: There is, at least, a duplicate in your powers exponents list: \n" ,
            argv[0]);
    si_vector_print_no_dims(stderr, powers_exponents);        
    fprintf(stderr,
            "This will produce a singular matrix. ");
    fprintf(stderr, "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    si_vector_free(powers_exponents);
    fclose(inv_matrix_file);
    unlink(argv[1]);
    fclose(points_file);
    unlink(argv[2]);
    return(EX_DATAERR);
  }
  /* Sorted the powers exponents */
  si_vector_insertion_sort(powers_exponents, SI_VECTOR_UP);
  /* Allocate the pseudo Vandermonde matrix. */
  pseudo_vander_matrix = mpq_matrix_alloc(num_exponents, num_exponents);
  if (pseudo_vander_matrix == NULL)
  {
    fprintf(stderr,
            "\n\n%s: Could not allocate the pseudo_vander_matrix.\n",
            argv[0]);
    fprintf(stderr, "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    si_vector_free(powers_exponents);
    fclose(inv_matrix_file);
    unlink(argv[1]);
    fclose(points_file);
    unlink(argv[2]);
    return(EX_SOFTWARE);
  }
  /* Get (compute) the points 
     There must be as many point as num_exponents (the number of monomials). 
  */
  points = compute_rat_points(rat_lower_bound,
                              rat_upper_bound,
                              num_exponents);
  if (points == NULL)
  {
    fprintf(stderr,
            "\n\n%s: Can not compute the points list.\n",
            argv[0]);
    fprintf(stderr,
            "Aborting the program!\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    si_vector_free(powers_exponents);
    fclose(inv_matrix_file);
    unlink(argv[1]);
    fclose(points_file);
    unlink(argv[2]);
    mpq_matrix_free(pseudo_vander_matrix);
    return(EX_SOFTWARE);
  }
  /* Compute the drift related stuff. */
  /* The drift increment is the value we add to (resp. substract from)
     the points to create a new matrix, expecting that this time
     it is no singular. */
  mpq_init(rat_drift_increment);
  /* The minimum spacing is the minimum value between to points as they
     drift one towards each other. If they drift to closely, we abort the
     the program. */
  mpq_init(rat_points_min_spacing);
  mpq_init(rat_temp_1);
  mpq_set_si(rat_temp_1, num_exponents - 1, 1);
  mpq_sub(rat_drift_increment, rat_upper_bound, rat_lower_bound);
  mpq_div(rat_drift_increment, rat_drift_increment, rat_temp_1);
  mpq_set_si(rat_temp_1, DRIFT_SLICE, 1);
  mpq_div(rat_drift_increment, rat_drift_increment, rat_temp_1);
  mpq_canonicalize(rat_drift_increment);
  drift_directions = si_vector_alloc(num_exponents);
  if (drift_directions == NULL)
  {
    fprintf(stderr,
            "\n\n%s: Could not allocate the drift_direction_vector.\n",
            argv[0]);
    fprintf(stderr, "Aborting the program!\n\n");
    mpq_clear(rat_lower_bound);
    mpq_clear(rat_upper_bound);
    si_vector_free(powers_exponents);
    fclose(inv_matrix_file);
    unlink(argv[1]);
    fclose(points_file);
    unlink(argv[2]);
    mpq_matrix_free(pseudo_vander_matrix);
    return(EX_SOFTWARE);
  }
  /* The drift direction is selected at random for each point. Bounds
     do not drift. */
  srandom(63 & time(NULL));
  mpq_sub(rat_points_min_spacing, rat_upper_bound, rat_lower_bound);
  mpq_set_si(rat_temp_1, num_exponents, 1);
  mpq_div(rat_points_min_spacing, rat_points_min_spacing, rat_temp_1);
  mpq_set_si(rat_temp_1, 2, 1);
  mpq_div(rat_points_min_spacing, rat_points_min_spacing, rat_temp_1);
  /* The drift direction is selected at random for each point. Bounds
     do not drift (see below). */
  for (i = 0 ; i < num_exponents ; i++)
  {
    if (random() & 1)
    {
      si_vector_set_at(drift_directions, i, 1);
    }
    else
    {
      si_vector_set_at(drift_directions, i, 0);
    }
  }
  /* For each point, compute the power values and add them to the 
     pseudo Vandermonde matrix. */
  for (i = 0 ; i < num_exponents  ; i++)
  {
    powers_values = mpq_arbitrary_int_powers_of(*(mpq_vector_get_at(points, i)),
                                            powers_exponents);
    if (powers_values == NULL)
    {
      fprintf(stderr,
              "\n\n%s: Can not compute the powers for value ",
              argv[0]);
      mpq_out_str(stderr, 
                  ALL_NUMBERS_BASE,
                  *(mpq_vector_get_at(points, i)));
      fprintf(stderr,
              "Aborting the program!\n");
      mpq_clear(rat_lower_bound);
      mpq_clear(rat_upper_bound);
      si_vector_free(powers_exponents);
      fclose(inv_matrix_file);
      unlink(argv[1]);
      fclose(points_file);
      unlink(argv[2]);
      mpq_matrix_free(pseudo_vander_matrix);
      mpq_clear(rat_temp_1);
      si_vector_free(drift_directions);
      return(EX_SOFTWARE);
    }
    mpq_matrix_set_row_from_vector(pseudo_vander_matrix, 
                                    i,
                                    powers_values);
    mpq_vector_free(powers_values);
  } /* End for i */
  inverted_pseudo_vander_matrix = mpq_matrix_clone(pseudo_vander_matrix);
  mpq_matrix_print(stderr, 10, inverted_pseudo_vander_matrix);
  mpq_init(rat_temp_2);
  while (mpq_matrix_inv_gauss_jordan(inverted_pseudo_vander_matrix))
  {
    mpq_vector_print_no_dims(stderr, 10, points);
    /* If the matrix is not invertible, try a new build by drifting
       the points. The points can not drift to much :
       - obviously, not beyond the bounds;
       - but not too close to each other as well. */
    if (mpq_vector_get_size(points) < 3)
    {
      fprintf(stderr,
              "\n\n%s: The matrix is singular and there are ",
              argv[0]);
      fprintf(stderr, "only two points (the bounds).\n");
      fprintf(stderr,
              "Aborting the program!\n");
      mpq_clear(rat_lower_bound);
      mpq_clear(rat_upper_bound);
      si_vector_free(powers_exponents);
      fclose(inv_matrix_file);
      unlink(argv[1]);
      fclose(points_file);
      unlink(argv[2]);
      mpq_matrix_free(pseudo_vander_matrix);
      mpq_clear(rat_temp_1);
      si_vector_free(drift_directions);
      mpq_matrix_free(inverted_pseudo_vander_matrix);
      mpq_clear(rat_temp_2);
      return(EX_SOFTWARE);
    }
    /* Remember, bounds do not drift, hence i=1 and i < num_exponents - 1. */
    for (i = 1 ; i < num_exponents - 1 ; i++)
    {
      /* Check that the min spacing between points is preserved. */
      mpq_set(rat_temp_2, *(mpq_vector_get_at(points, i - 1)));
      mpq_sub(rat_temp_2, *(mpq_vector_get_at(points, i)), rat_temp_2);
      if (mpq_cmp(rat_temp_2, rat_points_min_spacing) <= 0)
      {
        fprintf(stderr,
                "\n\n%s: Can not find a set of points that gives an invertible ",
                argv[0]);
        fprintf(stderr,
                "matrix.\nAborting the program!\n");
        mpq_clear(rat_lower_bound);
        mpq_clear(rat_upper_bound);
        si_vector_free(powers_exponents);
        fclose(inv_matrix_file);
        unlink(argv[1]);
        fclose(points_file);
        unlink(argv[2]);
        mpq_matrix_free(pseudo_vander_matrix);
        mpq_clear(rat_temp_1);
        si_vector_free(drift_directions);
        mpq_matrix_free(inverted_pseudo_vander_matrix);
        mpq_clear(rat_temp_2);
        return(EX_SOFTWARE);
      }
      mpq_set(rat_temp_2, *(mpq_vector_get_at(points, i + 1)));
      mpq_sub(rat_temp_2, rat_temp_2, *(mpq_vector_get_at(points, i)));
      if (mpq_cmp(rat_temp_2, rat_points_min_spacing) <= 0)
      {
        fprintf(stderr,
                "\n\n%s: Can not find a set of points that gives an invertible ",
                argv[0]);
        fprintf(stderr,
                "matrix.\nAborting the program!\n");
        mpq_clear(rat_lower_bound);
        mpq_clear(rat_upper_bound);
        si_vector_free(powers_exponents);
        fclose(inv_matrix_file);
        unlink(argv[1]);
        fclose(points_file);
        unlink(argv[2]);
        mpq_matrix_free(pseudo_vander_matrix);
        mpq_clear(rat_temp_1);
        si_vector_free(drift_directions);
        mpq_matrix_free(inverted_pseudo_vander_matrix);
        mpq_clear(rat_temp_2);
        return(EX_SOFTWARE);
      }
      mpq_set(rat_temp_1, *(mpq_vector_get_at(points, i)));
      if (*(si_vector_get_at(drift_directions, i)))
      {
        mpq_add(rat_temp_1, rat_temp_1, rat_drift_increment);
      }
      else
      {
        mpq_sub(rat_temp_1, rat_temp_1, rat_drift_increment);
      }
      mpq_canonicalize(rat_temp_1);
      mpq_vector_set_at(points, i, rat_temp_1);
    } /* End for i. */
    /* For each point, compute the power values and add them to the 
       pseudo Vandermonde matrix. */
    for (i = 0 ; i < num_exponents  ; i++)
    {
      powers_values = 
                  mpq_arbitrary_int_powers_of(*(mpq_vector_get_at(points, i)),
                                              powers_exponents);
      if (powers_values == NULL)
      {
        fprintf(stderr,
                "\n\n%s: Can not compute the powers for value ",
                argv[0]);
        mpq_out_str(stderr, 
                    ALL_NUMBERS_BASE,
                    *(mpq_vector_get_at(points, i)));
        fprintf(stderr,
                "Aborting the program!\n");
        mpq_clear(rat_lower_bound);
        mpq_clear(rat_upper_bound);
        si_vector_free(powers_exponents);
        fclose(inv_matrix_file);
        unlink(argv[1]);
        fclose(points_file);
        unlink(argv[2]);
        mpq_matrix_free(pseudo_vander_matrix);
        mpq_clear(rat_temp_1);
        si_vector_free(drift_directions);
        mpq_matrix_free(inverted_pseudo_vander_matrix);
        mpq_clear(rat_temp_2);
        return(EX_SOFTWARE);
      }
      mpq_matrix_set_row_from_vector(pseudo_vander_matrix, 
                                      i,
                                      powers_values);
      mpq_vector_free(powers_values);
    } /* End for i */
    mpq_matrix_free(inverted_pseudo_vander_matrix);
    inverted_pseudo_vander_matrix = mpq_matrix_clone(pseudo_vander_matrix);
  } /* Enb if */
  mpq_matrix_print(stderr, 10, inverted_pseudo_vander_matrix);
  identity_matrix = mpq_matrix_naive_prod(pseudo_vander_matrix,
                                          inverted_pseudo_vander_matrix);
  mpq_matrix_print(stderr, 10, identity_matrix);
  si_vector_print(stderr, drift_directions);
  mpq_matrix_print(inv_matrix_file, ALL_NUMBERS_BASE, inverted_pseudo_vander_matrix);
  mpq_vector_print(points_file, ALL_NUMBERS_BASE, points);
  fclose(inv_matrix_file);
  fclose(points_file);
  return(EX_OK);
} /* End main */

/* Other functions */
