/* mpq_vector.c */

#include <stdio.h>
#include <stdlib.h>
/*#include <string.>*/
#include <ctype.h>
#include <gmp.h>

#include "mpq-vector.h"
#include "../../misc-utils/std-exit-errors.h"

void
mpq_vector_errormsg(char *f , char *msgname, char *msg);

/* 
 * Allocate space for vector dimensioned by 'NbElems'.
 */
mpq_vector* mpq_vector_alloc(unsigned NbElems) {
  
  mpq_vector *Vect;
  mpq_t *p;
  int i;

  /* Deal with pathological parameter */
  if (NbElems < 0)
  {
    mpq_vector_errormsg("mpq_vector_alloc", "invalid parameter", "NbRows must be >= 0");
    return 0;
  }
  /* Create the vector struct */
  Vect=(mpq_vector *)malloc(sizeof(mpq_vector));
  if(!Vect) {	
    mpq_vector_errormsg("mpq_vector_Alloc", "outofmem", "out of memory space");
    return 0;
  }
  Vect->NbElems = NbElems;
  /* NbRows == 0: created an "empty" vector. */
  if(NbElems == 0) 
  {
    Vect->p_Init= (mpq_t *)0;
  }  
  else 
  {
    p = (mpq_t *)malloc(NbElems * sizeof(mpq_t));
    if(!p) 
    {
      free(Vect);
      mpq_vector_errormsg("mpq_vector_Alloc", "outofmem", "out of memory space");
      return 0;
    }
    Vect->p_Init = p;
    for (i = 0 ; i < NbElems ; i++) 
      {
        /* Init all elements and set them to value 0/1 */
        mpq_init(*(p + i));
        mpq_set_ui(*(p + i),0,1);
      }
  } /* End else */
  Vect->p_Init_size = NbElems;
  return Vect;
} /* mpq_vector_alloc */


/**
 * Create a vector from an array.
 */

mpq_vector* 
mpq_vector_from_array(mpq_t* mpq_array, unsigned int NbElems)
{
  mpq_vector* new_vector;
  if (NbElems == 0)
  {
    return((mpq_vector*)0);
  }
  new_vector = mpq_vector_alloc(0);
  if (! new_vector)
  {
    mpq_vector_errormsg("mpq_vector_from_array", "outofmem", "out of memory space");
    return((mpq_vector*)0);
  }  
  new_vector->p_Init = mpq_array;
  new_vector->p_Init_size = NbElems;
  new_vector->NbElems = NbElems;
  return(new_vector);
} /* End mpq_vector_from_array */

 
/**
 * Free the memory space occupied by mpq_vector 'Vect' 
 */
void mpq_vector_free(mpq_vector *Vect) 
{
  int i;
  mpq_t *p;
  
  if (Vect == NULL) return;

  if (Vect->p_Init)
  {
    p = Vect->p_Init;
    for(i = 0 ; i < Vect->p_Init_size ; i++) 
    {
      mpq_clear(*(p++));
    }
    free(Vect->p_Init);
  }
  free(Vect);
} /* mpq_vector_Free */

/* 
 * Print the contents of the mpq_vector 'Vect'.
 */
void mpq_vector_print(FILE *Dst, 
                        int radix, 
                        mpq_vector *Vect) 
{
  mpq_t *p;
  int i;
  unsigned NbElems;
  
  fprintf(Dst,"%d\n", NbElems=Vect->NbElems);
  p= Vect->p_Init;
  for (i=0 ; i < NbElems ; i++) 
  {
    mpq_out_str(Dst, radix, *p++);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* mpq_vector_print */

/* 
 * Print the contents of the mpq_vector 'Vect' without printing the
 * the dimensions.
 */
void mpq_vector_print_no_dims(FILE *Dst, 
                                int radix, 
                                mpq_vector *Vect) 
{
  
  mpq_t *p;
  int i;
  p= Vect->p_Init;
  for (i=0 ; i < Vect->p_Init_size ; i++) 
  {
    mpq_out_str(Dst, radix,*p++);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* mpq_vector_print_no_dims */

/* 
 * Print the contents of the mpq_vector 'Vect' in the Matlab
 * input format.
 */
void mpq_vector_print_matlab(FILE *Dst, 
                                int radix, 
                                mpq_vector *Vect) 
{
  mpq_t *p;
  int i;
  p= Vect->p_Init;
  fprintf(Dst, "[ ");
  for (i=0 ; i < Vect->p_Init_size ; i++) 
  {
    mpq_out_str(Dst, radix, *p++);
    if (i < (Vect->p_Init_size - 1))
    {
      fprintf(Dst, ", ");
    }
    fprintf(Dst, " ]\n");
  }
} /* mpq_vector_print_matlab */


/* 
 * Read the elements of the mpq_vector 'Vect' expressed as
 * strings convertible into mpq_t or into mpq_t (eventually converted into a
 * mpq_t).
 */
void mpq_vector_read_input_from_file(mpq_vector *Vect, 
                                      int radix, 
                                      FILE* Src) 
{
  mpq_t *p;
  int i,n;
  char *c, s[1024],str[1024];
  p = Vect->p_Init;
  do 
  {
    c = fgets(s, 1024, Src);
    while(isspace(*c) && *c != '\n')
	   ++c;
  } while(c && (*c == '#' || *c == '\n'));
  if (!c) 
  {
    mpq_vector_errormsg("mpq_vector_read_rat_input_from_file", 
                          "baddim", 
                          "no data" );
  }
  for (i = 0 ; i < Vect->NbElems ; i++) 
  {
    if(!c || *c=='\n' || *c=='#') 
    {
      mpq_vector_errormsg("mpq_vector_read_input_from_file", 
                            "baddim", 
                            "not enough elements");
      exit(EX_DATAERR);
    }
    if (sscanf(c,"%s%n", str, &n) == 0) 
    {
	    mpq_vector_errormsg("mpq_vector_read_input_from_file", 
                            "baddim", 
                            "not enough elements" );
      exit(EX_DATAERR);
    }
    /* Try to read the rational: mpq_set_str returns 0 if OK. */
    if (mpq_set_str(*p, str, radix))
    {
      mpq_vector_errormsg("mpq_vector_read_input_from_file", 
                          "badformat", 
                          "number in a bad format");
      exit(EX_DATAERR);
    }
    p++;
    c += n;
    }
} /* mpq_vector_read_input_from_file */

/* 
 * Read the contents of the vector 'Vect' from standard input. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the vector contents, one row per line.
 */
mpq_vector *mpq_vector_read_from_stdin(int radix) 
{
  
  mpq_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  while(fgets(s, 1024, stdin)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, stdin);
  Vect = mpq_vector_alloc(NbElems);
  if(!Vect) {
    mpq_vector_errormsg("mpq_vector_read", "outofmem", "out of memory space");
    return(NULL);
  }
  mpq_vector_read_input_from_file(Vect, radix, stdin);
  return Vect;
} /* mpq_vector_read_from_stdin */

/* 
 * Read the contents of the vector 'Vect' from FILE pointed to by Src. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the vector contents, one row per line.
 * Vectrix contents is made of strings directly convertible into mpq_t.
 */
mpq_vector*
mpq_vector_read_from_file(int radix, 
                            FILE* Src) 
{
  mpq_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  while(fgets(s, 1024, Src)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, Src);
  Vect = mpq_vector_alloc(NbElems);
  if(!Vect) {
    mpq_vector_errormsg("mpq_vector_read_from_file", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  mpq_vector_read_input_from_file(Vect, radix, Src);
  return Vect;
} /* mpq_vector_read_from_file */


mpq_t*
mpq_vector_get_at(mpq_vector* vect, unsigned int pos)
{
  if (vect == NULL)
  {
    mpq_vector_errormsg("mpq_vector_get_at", "null pointer parameter", "vect parameter is NULL");
    return(NULL);
  }
  if (vect->p_Init_size == 0)
  {
    mpq_vector_errormsg("mpq_vector_get_at", "empty vector", "can't return an element from an empty vector");
    return(NULL);
  }
  if (pos >= vect->NbElems)
  {
    mpq_vector_errormsg("mpq_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  if (pos < 0)
  {
    mpq_vector_errormsg("mpq_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  return(&(vect->p_Init[pos]));
} /* End mpq_vector_get_at */

void 
mpq_vector_set_at(mpq_vector* vect, 
                    unsigned int pos, 
                    mpq_t value)
{
  if (vect == NULL)
  {
    mpq_vector_errormsg("mpq_vector_set_at", "null pointer parameter", "mat parameter is NULL");
  }
  if (vect->p_Init_size == 0)
  {
    mpq_vector_errormsg("mpq_vector_set_at", "empty vector", "can't set an element from an empty vector");
  }
  if (pos >= vect->NbElems)
  {
    mpq_vector_errormsg("mpq_vector_set_at", "index out of range", "pos index beyond bounds");
  }
  mpq_set(vect->p_Init[pos], value);
} /* End mpq_vector_set_at */

unsigned int
mpq_vector_get_NbElems(mpq_vector* vect)
{
  return(vect->NbElems);
} /* End mpq_vector_get_NbElems */

unsigned int
mpq_vector_get_size(mpq_vector* vect)
{
  return(mpq_vector_get_NbElems(vect));
} /* End mpq_vector_get_size */

/**
 * @see mpq-vector.c#mpq_vector_has_duplicates
 */
int
mpq_vector_has_duplicates(mpq_vector* vect)
{
  unsigned int nb_elems;
  unsigned int i, j;
  mpq_t rat_temp;
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    return(0);
  }
  if ((nb_elems = mpq_vector_get_size(vect)) == 0)
  {
    return(0);
  }
  mpq_init(rat_temp);
  for (i = 0 ; i < (nb_elems - 1) ; i++)
  {
    mpq_set(rat_temp, *(mpq_vector_get_at(vect, i)));
    for (j = i + 1  ; j < nb_elems ; j++)
    {
      if (mpq_equal(rat_temp, *(mpq_vector_get_at(vect, j))))
      {
        mpq_clear(rat_temp);
        return(1);
      }
    }
  }
  mpq_clear(rat_temp);
  return(0);
  mpq_clear(rat_temp);
}

/**
 * @see mpq-vector.c#mpq_vector_get_abs
 */
mpq_vector*
mpq_vector_get_abs(mpq_vector* vect)
{
  mpq_vector* abs_vect;
  unsigned int num_elems;
  unsigned int i;
  mpq_t tmp;
  
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    mpq_vector_errormsg("mpq_vector_get_abs", 
                          "null pointer parameter", 
                          "vect parameter is NULL");
    return(NULL);
  }
  num_elems = mpq_vector_get_size(vect);
  abs_vect = mpq_vector_alloc(num_elems);
  if (abs_vect == NULL)
  {
    mpq_vector_errormsg("mpq_vector_get_abs", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  mpq_init(tmp);
  for (i = 0 ; i < num_elems ; i++)
  {
    mpq_abs(tmp, *(mpq_vector_get_at(vect, i)));
    mpq_vector_set_at(abs_vect, i, tmp);
  }
  mpq_clear(tmp);
  return(abs_vect);
} /* End mpq_vector_get_abs */

/**
 * see mpq-vector.c#mpq_vector_insertion_sort
 */
void
mpq_vector_insertion_sort(mpq_vector* vi, int direction)
{
  mpq_vector_errormsg("mpq_vector_insertion_sort",
                      "notImplementedYet",
                      "this function is not implemented yet!");
  exit(1);
          
} /* End mpq_vector_insertion_sort */

/**
 * see mpq-vector.h#mpq_vector_are_equal
 */
int 
mpq_vector_are_equal(mpq_vector* vect_1, mpq_vector* vect_2)
{
  unsigned int nb_elems;
  unsigned int i;
  /* Deal with pathological parameters */
  if ((vect_1 == NULL) && (vect_2 == NULL))
  {
    return(1);
  }
  if ((vect_1 == NULL) || (vect_2 == NULL))
  {
    mpq_vector_errormsg("si_vector_are_equal", 
                        "null pointer parameter", 
                        "one of the vectors (and only one) is NULL");
    return(0);
  }
  nb_elems = mpq_vector_get_size(vect_1);
  if ((nb_elems == 0) && (mpq_vector_get_size(vect_2) == 0))
  {
    return(1);
  }
  if (nb_elems != mpq_vector_get_size(vect_2))
  {
    return(0);
  }
  /* Actual check for equality. */
  for (i = 0 ; i < nb_elems ; i++)
  {
    if (mpq_equal(*(mpq_vector_get_at(vect_1, i)),
                  *(mpq_vector_get_at(vect_2, i))))
    {
      return(0);
    }
  } /* End for. */
  return(1);
} /* End mpq_vector_are_equal */

int
mpq_vector_vector_product(mpq_vector* v1, mpq_vector* v2, mpq_t* result)
{
  int i = 0;
  mpq_t part_prod;
  mpq_t* v1_p_Init;
  mpq_t* v2_p_Init;
  /* Check for incompatible vector sizes. */
  if (v1->p_Init_size != v2->p_Init_size)
  {
    mpq_vector_errormsg("mpq_vector_vector_product", "vector size conflict", "vectors have different sizes");
    return(1);
  }
  v1_p_Init = v1->p_Init;
  v2_p_Init = v2->p_Init;
  mpq_init(part_prod);
  mpq_init(*result);
  mpq_set_ui(*result, 0, 0);
  for (i = 0 ; i < v1->p_Init_size ; i++)
  {
    mpq_mul(part_prod, *v1_p_Init, *v2_p_Init);
    mpq_add(*result, *result, part_prod);
    v1_p_Init++;
    v2_p_Init++;
  } /* End for */
  mpq_clear(part_prod);
  return(0);
} /* End mpq_vector_vector_product */

void 
mpq_vector_errormsg(char *f , char *msgname, char *msg) 
{
  fprintf(stderr, "?%s: %s\n", f, msg);
}


