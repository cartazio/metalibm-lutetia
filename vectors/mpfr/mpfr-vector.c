/* mpfr_vector.c */

#include <stdio.h>
#include <stdlib.h>
/*#include <string.>*/
#include <ctype.h>
#include <gmp.h>
#include <mpfr.h>

#include "mpfr-vector.h"
#include "../../misc-utils/std-exit-errors.h"

void
mpfr_vector_errormsg(char *f , char *msgname, char *msg);

/* 
 * Allocate space for vector dimensioned by 'NbElems'.
 */
mpfr_vector* mpfr_vector_alloc(unsigned NbElems, 
                                unsigned precision) 
{
  
  mpfr_vector *Vect;
  mpfr_t *p;
  int i;

  mpfr_set_default_prec(precision);

  /* Deal with pathological parameter */
  if (NbElems < 0)
  {
    mpfr_vector_errormsg("mpfr_vector_alloc", "invalid parameter", "NbElems must be >= 0");
    return 0;
  }
  /* Create the vector struct */
  Vect=(mpfr_vector *)malloc(sizeof(mpfr_vector));
  if(!Vect) {	
    mpfr_vector_errormsg("mpfr_vector_alloc", "outofmem", "out of memory space");
    return 0;
  }
  Vect->NbElems = NbElems;
  /* NbElems == 0: created an "empty" vector. */
  if(NbElems == 0) 
  {
    Vect->p_Init= (mpfr_t *)0;
  }  
  else 
  {
    p = (mpfr_t *)malloc(NbElems * sizeof(mpfr_t));
    if(!p) 
    {
      free(Vect);
      mpfr_vector_errormsg("mpfr_vector_Alloc", "outofmem", "out of memory space");
      return 0;
    }
    Vect->p_Init = p;
    for (i = 0 ; i < NbElems ; i++) 
      {
        mpfr_init_set_ui(*(p+i),0,GMP_RNDN);
      }
  }
  return Vect;
} /* mpfr_vector_alloc */


/**
 * Create a vector from an array.
 */

mpfr_vector* 
mpfr_vector_from_array(mpfr_t* mpfr_array, 
                        unsigned int NbElems,
                        unsigned int precision)
{
  mpfr_vector* new_vector;
  int i;
  
  /* Deal first with pathological values */
  if (NbElems < 0)
  {
    mpfr_vector_errormsg("mpfr_vector_from_array",
                          "invalid number of elements",
                          "negative number of elements for a vector");
    exit(1);
  }
  if (NbElems == 0)
  {
    return((mpfr_vector*)0);
  }
  new_vector = mpfr_vector_alloc(NbElems, precision);
  if (! new_vector)
  {
    mpfr_vector_errormsg("mpfr_vector_from_array", 
                          "outofmem", 
                          "out of memory space");
    return((mpfr_vector*)0);
  }
  for (i = 0 ; i < NbElems ; i++)
  {
    mpfr_vector_set_at(new_vector,
                       i,
                       mpfr_array[i],
                       GMP_RNDN);
  }  
  new_vector->NbElems = NbElems;
  return(new_vector);
} /* End mpfr_vector_from_array */

 
/**
 * Free the memory space occupied by mpfr_vector 'Vect' 
 */
void mpfr_vector_free(mpfr_vector *Vect) 
{
  int i;
  mpfr_t *p;
  
  if (Vect == NULL) return;

  if (Vect->p_Init)
  {
    p = Vect->p_Init;
    for(i = 0 ; i < Vect->NbElems ; i++) 
    {
      mpfr_clear(*p++);
    }
    free(Vect->p_Init);
  }
  free(Vect);
} /* mpfr_vector_Free */

/* 
 * Print the contents of the mpfr_vector 'Vect'.
 */
void mpfr_vector_print(FILE *Dst, 
                        int radix, 
                        size_t digits, 
                        mp_rnd_t rounding_mode,
                        mpfr_vector *Vect) 
{
  mpfr_t *p;
  int i;
  unsigned NbElems;
  
  fprintf(Dst,"%d\n", NbElems=Vect->NbElems);
  p= Vect->p_Init;
  for (i=0 ; i < NbElems ; i++) 
  {
    mpfr_out_str(Dst, radix, digits, *p++, rounding_mode);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* mpfr_vector_print */

/* 
 * Print the contents of the mpfr_vector 'Vect' without printing the
 * the dimensions.
 */
void mpfr_vector_print_no_dims(FILE *Dst, 
                                int radix, 
                                size_t digits, 
                                mp_rnd_t rounding_mode,
                                mpfr_vector *Vect) 
{
  
  mpfr_t *p;
  int i;
  p= Vect->p_Init;
  for (i=0 ; i < Vect->NbElems ; i++) 
  {
    mpfr_out_str(Dst, radix, digits, *p++, rounding_mode);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* mpfr_vector_print_no_dims */

/* 
 * Print the contents of the mpfr_vector 'Vect' in the Matlab
 * input format.
 */
void mpfr_vector_print_matlab(FILE *Dst, 
                                int radix, 
                                size_t digits, 
                                mp_rnd_t rounding_mode,
                                mpfr_vector *Vect) 
{
  mpfr_t *p;
  int i;
  p= Vect->p_Init;
  fprintf(Dst, "[ ");
  for (i=0 ; i < Vect->NbElems ; i++) 
  {
    mpfr_out_str(Dst, radix, digits, *p++, rounding_mode);
    if (i < (Vect->NbElems - 1))
    {
      fprintf(Dst, ", ");
    }
    fprintf(Dst, " ]\n");
  }
} /* mpfr_vector_print_matlab */

/* 
 * Read the elements of the mpfr_vector 'Vect' expressed as
 * strings convertible into mpfr_t or into mpq_t (eventually converted into a
 * mpfr_t).
 */
void mpfr_vector_read_input_from_file(mpfr_vector *Vect, 
                                      int radix, 
                                      mp_rnd_t rounding_mode,
                                      unsigned int precision,
                                      FILE* Src) 
{
  mpfr_t *p;
  mpq_t  rat_temp;
  int i,n;
  char *c, s[1024],str[1024];

  mpfr_set_default_prec(precision);
  
  p = Vect->p_Init;
  /* Move to the line that actually holds the vector elements, skipping
     comments and empty lines */
  do 
  {
    c = fgets(s, 1024, Src);
    while(isspace(*c) && *c != '\n')
	   ++c;
  } while(c && (*c == '#' || *c == '\n'));
  if (!c) 
  {
    mpfr_vector_errormsg("mpfr_vector_read_rat_input_from_file", 
                          "baddim", 
                          "no data" );
  }
  for (i = 0 ; i < Vect->NbElems ; i++) 
  {
    if(!c || *c=='\n' || *c=='#') 
    {
      mpfr_vector_errormsg("mpfr_vector_read_input_from_file", 
                            "baddim", 
                            "not enough elements");
      exit(EX_DATAERR);
    }
    if (sscanf(c,"%s%n", str, &n) == 0) 
    {
	    mpfr_vector_errormsg("mpfr_vector_read_input_from_file", 
                            "baddim", 
                            "not enough elements" );
      exit(EX_DATAERR);
    }
      /* Reading a float:
         - first try to read it as mpfr_t;
         - if this fails, try to read it as a mpq_t and then transform it into
           a mpfr_t */
      /* Try to read the float directly: mpfr_set_str returns 0 if OK. */
      if (mpfr_set_str(*p, str, radix, rounding_mode))
      {
        /* Try to read it as a rational: mpq_set_str returns 0 if OK. */
        mpq_init(rat_temp);
        if (mpq_set_str(rat_temp, str, radix))
        {
          mpfr_vector_errormsg("mpfr_vector_read_input_from_file", 
                                "badformat", 
                                "number in a bad format");
          mpq_clear(rat_temp);
          exit(EX_DATAERR);
        }
        mpfr_set_q(*p, rat_temp, GMP_RNDN);
        mpq_clear(rat_temp);
      }
      p++;
      c += n;
    }
} /* mpfr_vector_read_input_from_file */

/* 
 * Read the contents of the vector 'Vect' from standard input. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of elements, the next 
 * line holds the vector contents.
 */
mpfr_vector *mpfr_vector_read_from_stdin(int radix, 
                                          mp_rnd_t rounding_mode, 
                                          unsigned int precision) 
{
  
  mpfr_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  mpfr_set_default_prec(precision);
  
  while(fgets(s, 1024, stdin)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, stdin);
  Vect = mpfr_vector_alloc(NbElems, precision);
  if(!Vect) {
    mpfr_vector_errormsg("mpfr_vector_read", "outofmem", "out of memory space");
    return(NULL);
  }
  mpfr_vector_read_input_from_file(Vect, radix, rounding_mode, precision, stdin);
  return Vect;
} /* mpfr_vector_read_from_stdin */

/* 
 * Read the contents of the vector 'Vect' from FILE pointed to by Src. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, the next one
 * holds the vector contents.
 * Vector contents is made of strings directly convertible into mpfr_t.
 */
mpfr_vector*
mpfr_vector_read_from_file(int radix, 
                            mp_rnd_t rounding_mode,
                            unsigned int precision, 
                            FILE* Src) 
{
  mpfr_vector *Vect;
  unsigned NbElems;
  char s[1024];
 
  mpfr_set_default_prec(precision);
  
  while(fgets(s, 1024, Src)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, Src);
  Vect = mpfr_vector_alloc(NbElems, precision);
  if(!Vect) {
    mpfr_vector_errormsg("mpfr_vector_read_from_file", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  mpfr_vector_read_input_from_file(Vect, radix, rounding_mode, precision, Src);
  return Vect;
} /* mpfr_vector_read_from_file */

mpfr_t*
mpfr_vector_get_at(mpfr_vector* vect, unsigned int pos)
{

  mpfr_set_default_prec(200);

  if (vect == NULL)
  {
    mpfr_vector_errormsg("mpfr_vector_get_at", "null pointer parameter", "mat parameter is NULL");
    return(NULL);
  }
  if (vect->NbElems == 0)
  {
    mpfr_vector_errormsg("mpfr_vector_get_at", "empty vector", "can't return an element from an empty vector");
    return(NULL);
  }
  if (pos >= vect->NbElems)
  {
    mpfr_vector_errormsg("mpfr_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  if (pos < 0)
  {
    mpfr_vector_errormsg("mpfr_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  return(&(vect->p_Init[pos]));
} /* End mpfr_vector_get_at */

void 
mpfr_vector_set_at(mpfr_vector* vect, 
                    unsigned int pos, 
                    mpfr_t value,
                    mp_rnd_t rounding_mode)
{

  mpfr_set_default_prec(200);
  
  if (vect == NULL)
  {
    mpfr_vector_errormsg("mpfr_vector_set_at", 
                          "null pointer parameter", 
                          "vect parameter is NULL");
  }
  if (vect->NbElems == 0)
  {
    mpfr_vector_errormsg("mpfr_vector_set_at", 
                          "empty vector", 
                          "can't set an element from an empty vector");
  }
  if (pos >= vect->NbElems)
  {
    mpfr_vector_errormsg("mpfr_vector_set_at", 
                          "indexOutOfRange", 
                          "pos index beyond bounds");
  }
  mpfr_set(vect->p_Init[pos], value, rounding_mode);
} /* End mpfr_vector_set_at */

unsigned int
mpfr_vector_get_size(mpfr_vector* vect)
{
  return(vect->NbElems);
} /* End mpfr_vector_get_size */

/**
 * @see mpfr-vector.h#mpfr_vector_has_duplicates
 */
int
mpfr_vector_has_duplicates(mpfr_vector* vect)
{
  unsigned int nb_elems;
  unsigned int i, j;
  mpfr_t   mpfr_temp;
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    return(0);
  }
  if ((nb_elems = mpfr_vector_get_size(vect)) == 0)
  {
    return(0);
  }
  mpfr_init(mpfr_temp);
  for (i = 0 ; i < (nb_elems - 1) ; i++)
  {
    mpfr_set(mpfr_temp, *(mpfr_vector_get_at(vect, i)), GMP_RNDN);
    for (j = i + 1  ; j < nb_elems ; j++)
    {
      if (mpfr_cmp(mpfr_temp, *(mpfr_vector_get_at(vect, j))) == 0)
      {
        mpfr_clear(mpfr_temp);
        return(1);
      }
    } /* End for j. */
  } /* End for i. */
  mpfr_clear(mpfr_temp);
  return(0);
} /* End mpfr_vector_has_duplicates */


/**
 * @see mpfr-vector.h#mpfr_vector_get_abs
 */
mpfr_vector*
mpfr_vector_get_abs(mpfr_vector* vect, unsigned int precision)
{
  mpfr_vector* abs_vect;
  unsigned int num_elems;
  unsigned int i;
  mpfr_t tmp;
  
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    mpfr_vector_errormsg("mpfr_vector_get_abs", 
                          "null pointer parameter", 
                          "vect parameter is NULL");
    return(NULL);
  }
  num_elems = mpfr_vector_get_size(vect);
  if (num_elems == 0)
  {
    mpfr_vector_errormsg("mpfr_vector_get_abs", 
                          "empty vector", 
                          "vect parameter is an empty vector");
    return(NULL);
  }
  abs_vect = mpfr_vector_alloc(num_elems, precision);
  if (abs_vect == NULL)
  {
    mpfr_vector_errormsg("mpfr_vector_get_abs", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  mpfr_init(tmp);
  for (i = 0 ; i < num_elems ; i++)
  {
    mpfr_abs(tmp, *(mpfr_vector_get_at(vect, i)), GMP_RNDN);
    mpfr_vector_set_at(abs_vect, i, tmp, GMP_RNDN);
  }
  mpfr_clear(tmp);
  return(abs_vect);
} /* End mpfr_vector_get_abs */


int
mpfr_vector_scalar_product(mpfr_vector* v1, 
                            mpfr_vector* v2, 
                            mpfr_t* result, 
                            unsigned int precision)
{
  int i = 0;
  mpfr_t part_prod;
  mpfr_t* v1_p_Init;
  mpfr_t* v2_p_Init;

  mpfr_set_default_prec(precision);
  
  /* Check for incompatible vector sizes. */
  if (v1->NbElems != v2->NbElems)
  {
    mpfr_vector_errormsg("mpfr_scalar_product", "vectorSizeConflict", "vectors have different sizes");
    return(1);
  }
  v1_p_Init = v1->p_Init;
  v2_p_Init = v2->p_Init;
  mpfr_init(part_prod);
  mpfr_set_ui(*result, 0, GMP_RNDN);
  for (i = 0 ; i < v1->NbElems ; i++)
  {
    mpfr_mul(part_prod, *v1_p_Init, *v2_p_Init, GMP_RNDN);
    mpfr_add(*result, *result, part_prod, GMP_RNDN);
    v1_p_Init++;
    v2_p_Init++;
  } /* End for */
  mpfr_clear(part_prod);
  return(0);
} /* End mpfr_vector_scalar_product */

int
mpfr_vector_elems_sum(mpfr_vector* vect, 
                            mpfr_t* sum, 
                            unsigned int precision)
{
  unsigned int num_elems;
  unsigned int i;
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    mpfr_vector_errormsg("mpfr_vector_elems_sum", 
                          "nullPointerParameter", 
                          "vect is NULL");
    return(1);
  }
  if (sum == NULL)
  {
    mpfr_vector_errormsg("mpfr_vector_elems_sum", 
                          "nullPointerParameter", 
                          "sum is NULL");
    return(1);
  }
  /* Do something usefull */
  num_elems = mpfr_vector_get_size(vect);
  mpfr_set_ui(*sum, 0, GMP_RNDN);
  for (i = 0 ; i < num_elems ; i++)
  {
    mpfr_add(*sum, 
              *sum,
              *(mpfr_vector_get_at(vect, i)),
              GMP_RNDN); 
  } /* End for i. */
  return(0);
} /* End mpfr_vector_elems_sum */

/**
 * @see mpfr-vector.h#mpfr_vector_set_to_max
 */

mpfr_vector*
mpfr_vector_set_to_max (mpfr_vector* in_out_vect,
                        mpfr_vector* other_vect)
{
  unsigned int num_elems;
  unsigned int i;
  /* Deal with pathological parameters */
  if ((in_out_vect == NULL) || (other_vect == NULL))
  {
    mpfr_vector_errormsg("mpfr_vector_set_to_max", 
                          "nullPointerParameter", 
                          "one of in_out_vect or other_vect is NULL");
    return(NULL);
  }
  num_elems = mpfr_vector_get_size(in_out_vect);
  if (num_elems != mpfr_vector_get_size(other_vect))
  {
    mpfr_vector_errormsg("mpfr_vector_set_to_max", 
                          "vector size conflict", 
                          "vectors have different sizes");
    return(NULL);
  }
  for (i = 0 ; i < num_elems ; i++)
  {
    if (mpfr_cmp(*(mpfr_vector_get_at(in_out_vect, i)),
                  *(mpfr_vector_get_at(other_vect, i))) < 0)
    {
      mpfr_vector_set_at(in_out_vect,
                          i,
                          *(mpfr_vector_get_at(other_vect, i)),
                          GMP_RNDN);
    }
  } /* End for i. */
  return(in_out_vect);
} /* End mpfr_vector_set_to_max. */

/**
 * @see mpfr-vector.h#mpfr_vector_set_to_min
 */

mpfr_vector*
mpfr_vector_set_to_min (mpfr_vector* in_out_vect,
                        mpfr_vector* other_vect)
{
  unsigned int num_elems;
  unsigned int i;
  /* Deal with pathological parameters */
  if ((in_out_vect == NULL) || (other_vect == NULL))
  {
    mpfr_vector_errormsg("mpfr_vector_set_to_max", 
                          "nullPointerParameter", 
                          "one of in_out_vect or other_vect is NULL");
    return(NULL);
  }
  num_elems = mpfr_vector_get_size(in_out_vect);
  if (num_elems != mpfr_vector_get_size(other_vect))
  {
    mpfr_vector_errormsg("mpfr_vector_set_to_max", 
                          "vector size conflict", 
                          "vectors have different sizes");
    return(NULL);
  }
  for (i = 0 ; i < num_elems ; i++)
  {
    if (mpfr_cmp(*(mpfr_vector_get_at(in_out_vect, i)),
                  *(mpfr_vector_get_at(other_vect, i))) > 0)
    {
      mpfr_vector_set_at(in_out_vect,
                          i,
                          *(mpfr_vector_get_at(other_vect, i)),
                          GMP_RNDN);
    }
  } /* End for i. */
  return(in_out_vect);
} /* End mpfr_vector_set_to_min. */

void 
mpfr_vector_errormsg(char *f , char *msgname, char *msg) 
{
  fprintf(stderr, "?%s: %s\n", f, msg);
}


