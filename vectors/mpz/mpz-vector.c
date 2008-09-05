/* mpz_vector.c */

#include <stdio.h>
#include <stdlib.h>
/*#include <string.>*/
#include <ctype.h>
#include <gmp.h>

#include "mpz-vector.h"
#include "../../misc-utils/std-exit-errors.h"

void
mpz_vector_errormsg(char *f , char *msgname, char *msg);

/* 
 * Allocate space for vector dimensioned by 'NbElems'.
 */
mpz_vector* mpz_vector_alloc(unsigned NbElems) {
  
  mpz_vector *Vect;
  mpz_t *p;
  int i;

  /* Deal with pathological parameter */
  if (NbElems < 0)
  {
    mpz_vector_errormsg("mpz_vector_alloc", "invalid parameter", "NbRows must be >= 0");
    return 0;
  }
  /* Create the vector struct */
  Vect=(mpz_vector *)malloc(sizeof(mpz_vector));
  if(!Vect) {	
    mpz_vector_errormsg("mpz_vector_Alloc", "outofmem", "out of memory space");
    return 0;
  }
  Vect->NbElems = NbElems;
  /* NbRows == 0: created an "empty" vector. */
  if(NbElems == 0) 
  {
    Vect->p_Init= (mpz_t *)0;
  }  
  else 
  {
    p = (mpz_t *)malloc(NbElems * sizeof(mpz_t));
    if(!p) 
    {
      free(Vect);
      mpz_vector_errormsg("mpz_vector_Alloc", "outofmem", "out of memory space");
      return 0;
    }
    Vect->p_Init = p;
    for (i = 0 ; i < NbElems ; i++) 
      {
        mpz_init(*(p+i));
        mpz_set_ui(*(p+i),0);
      }
  }
  Vect->p_Init_size = NbElems;
  return Vect;
} /* mpz_vector_alloc */


/**
 * Create a vector from an array.
 */

mpz_vector* 
mpz_vector_from_array(mpz_t* mpz_array, unsigned int NbElems)
{
  mpz_vector* new_vector;
  if (NbElems == 0)
  {
    return((mpz_vector*)0);
  }
  new_vector = mpz_vector_alloc(0);
  if (! new_vector)
  {
    mpz_vector_errormsg("mpz_vector_from_array", "outofmem", "out of memory space");
    return((mpz_vector*)0);
  }  
  new_vector->p_Init = mpz_array;
  new_vector->p_Init_size = NbElems;
  new_vector->NbElems = NbElems;
  return(new_vector);
} /* End mpz_vector_from_array */

 
/**
 * Free the memory space occupied by mpz_vector 'Vect' 
 */
void mpz_vector_free(mpz_vector *Vect) 
{
  int i;
  mpz_t *p;
  
  if (Vect == NULL) return;

  if (Vect->p_Init)
  {
    p = Vect->p_Init;
    for(i=0 ; i < Vect->p_Init_size ; i++) 
    {
      mpz_clear(*p++);
    }
    free(Vect->p_Init);
  }
  free(Vect);
} /* mpz_vector_Free */

/* 
 * Print the contents of the mpz_vector 'Vect'.
 */
void mpz_vector_print(FILE *Dst, 
                        int radix, 
                        mpz_vector *Vect) 
{
  mpz_t *p;
  int i;
  unsigned NbElems;
  
  fprintf(Dst,"%d\n", NbElems=Vect->NbElems);
  p= Vect->p_Init;
  for (i=0 ; i < NbElems ; i++) 
  {
    mpz_out_str(Dst, radix, *p++);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* mpz_vector_print */

/* 
 * Print the contents of the mpz_vector 'Vect' without printing the
 * the dimensions.
 */
void mpz_vector_print_no_dims(FILE *Dst, 
                                int radix, 
                                mpz_vector *Vect) 
{
  
  mpz_t *p;
  int i;
  p= Vect->p_Init;
  for (i=0 ; i < Vect->p_Init_size ; i++) 
  {
    mpz_out_str(Dst, radix,*p++);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* mpz_vector_print_no_dims */

/* 
 * Print the contents of the mpz_vector 'Vect' in the Matlab
 * input format.
 */
void mpz_vector_print_matlab(FILE *Dst, 
                                int radix, 
                                mpz_vector *Vect) 
{
  mpz_t *p;
  int i;
  p= Vect->p_Init;
  fprintf(Dst, "[ ");
  for (i=0 ; i < Vect->p_Init_size ; i++) 
  {
    mpz_out_str(Dst, radix, *p++);
    if (i < (Vect->p_Init_size - 1))
    {
      fprintf(Dst, ", ");
    }
    fprintf(Dst, " ]\n");
  }
} /* mpz_vector_print_matlab */


/* 
 * Read the elements of the mpz_vector 'Vect' expressed as
 * strings convertible into mpz_t or into mpz_t (eventually converted into a
 * mpz_t).
 */
void mpz_vector_read_input_from_file(mpz_vector *Vect, 
                                      int radix, 
                                      FILE* Src) 
{
  mpz_t *p;
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
    mpz_vector_errormsg("mpz_vector_read_rat_input_from_file", 
                          "baddim", 
                          "no data" );
  }
  for (i = 0 ; i < Vect->NbElems ; i++) 
  {
    if(!c || *c=='\n' || *c=='#') 
    {
      mpz_vector_errormsg("mpz_vector_read_input_from_file", 
                            "baddim", 
                            "not enough elements");
      exit(EX_DATAERR);
    }
    if (sscanf(c,"%s%n", str, &n) == 0) 
    {
	    mpz_vector_errormsg("mpz_vector_read_input_from_file", 
                            "baddim", 
                            "not enough elements" );
      exit(EX_DATAERR);
    }
    /* Try to read the rational: mpz_set_str returns 0 if OK. */
    if (mpz_set_str(*p, str, radix))
    {
      mpz_vector_errormsg("mpz_vector_read_input_from_file", 
                          "badformat", 
                          "number in a bad format");
      exit(EX_DATAERR);
    }
    p++;
    c += n;
    }
} /* mpz_vector_read_input_from_file */

/* 
 * Read the contents of the vector 'Vect' from standard input. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the vector contents, one row per line.
 */
mpz_vector *mpz_vector_read_from_stdin(int radix) 
{
  
  mpz_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  while(fgets(s, 1024, stdin)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, stdin);
  Vect = mpz_vector_alloc(NbElems);
  if(!Vect) {
    mpz_vector_errormsg("mpz_vector_read", "outofmem", "out of memory space");
    return(NULL);
  }
  mpz_vector_read_input_from_file(Vect, radix, stdin);
  return Vect;
} /* mpz_vector_read_from_stdin */

/* 
 * Read the contents of the vector 'Vect' from FILE pointed to by Src. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the vector contents, one row per line.
 * Vectrix contents is made of strings directly convertible into mpz_t.
 */
mpz_vector*
mpz_vector_read_from_file(int radix, 
                            FILE* Src) 
{
  mpz_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  while(fgets(s, 1024, Src)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, Src);
  Vect = mpz_vector_alloc(NbElems);
  if(!Vect) {
    mpz_vector_errormsg("mpz_vector_read_from_file", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  mpz_vector_read_input_from_file(Vect, radix, Src);
  return Vect;
} /* mpz_vector_read_from_file */


mpz_t*
mpz_vector_get_at(mpz_vector* vect, unsigned int pos)
{
  if (vect == NULL)
  {
    mpz_vector_errormsg("mpz_vector_get_at", "null pointer parameter", "mat parameter is NULL");
    return(NULL);
  }
  if (vect->p_Init_size == 0)
  {
    mpz_vector_errormsg("mpz_vector_get_at", "empty vector", "can't return an element from an empty vector");
    return(NULL);
  }
  if (pos >= vect->NbElems)
  {
    mpz_vector_errormsg("mpz_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  if (pos < 0)
  {
    mpz_vector_errormsg("mpz_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  return(&(vect->p_Init[pos]));
} /* End mpz_vector_get_at */

void 
mpz_vector_set_at(mpz_vector* vect, 
                    unsigned int pos, 
                    mpz_t value)
{
  if (vect == NULL)
  {
    mpz_vector_errormsg("mpz_vector_set_at", "null pointer parameter", "mat parameter is NULL");
  }
  if (vect->p_Init_size == 0)
  {
    mpz_vector_errormsg("mpz_vector_set_at", "empty vector", "can't set an element from an empty vector");
  }
  if (pos >= vect->NbElems)
  {
    mpz_vector_errormsg("mpz_vector_set_at", "index out of range", "pos index beyond bounds");
  }
  mpz_set(vect->p_Init[pos], value);
} /* End mpz_vector_set_at */

int
mpz_vector_vector_product(mpz_vector* v1, mpz_vector* v2, mpz_t* result)
{
  int i = 0;
  mpz_t part_prod;
  mpz_t* v1_p_Init;
  mpz_t* v2_p_Init;
  /* Check for incompatible vector sizes. */
  if (v1->p_Init_size != v2->p_Init_size)
  {
    mpz_vector_errormsg("mpz_vector_vector_product", "vector size conflict", "vectors have different sizes");
    return(1);
  }
  v1_p_Init = v1->p_Init;
  v2_p_Init = v2->p_Init;
  mpz_init(part_prod);
  mpz_init(*result);
  mpz_set_ui(*result, 0);
  for (i = 0 ; i < v1->p_Init_size ; i++)
  {
    mpz_mul(part_prod, *v1_p_Init, *v2_p_Init);
    mpz_add(*result, *result, part_prod);
    v1_p_Init++;
    v2_p_Init++;
  } /* End for */
  mpz_clear(part_prod);
  return(0);
} /* End mpz_vector_vector_product */

/**
 * @see mpz-vector.h#mpz_vector_get_NbElems
 */
unsigned int 
mpz_vector_get_NbElems(mpz_vector* vect)
{
  return(vect->NbElems);
} /* End mpz_vector_get_NbElems. */

/**
 * @see mpz-vector.h#mpz_vector_get_size
 */
unsigned int 
mpz_vector_get_size(mpz_vector* vect)
{
  return(mpz_vector_get_NbElems(vect));
} /* End mpz_vector_get_size. */


void 
mpz_vector_errormsg(char *f , char *msgname, char *msg) 
{
  fprintf(stderr, "?%s: %s\n", f, msg);
}


