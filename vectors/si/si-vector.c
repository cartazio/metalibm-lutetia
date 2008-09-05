/* si_vector.c */

#include <stdio.h>
#include <stdlib.h>
/*#include <string.>*/
#include <ctype.h>

#include "si-vector.h"

void
si_vector_errormsg(char *f , char *msgname, char *msg);

/* 
 * Allocate space for vector dimensioned by 'NbElems'.
 */
si_vector* si_vector_alloc(unsigned NbElems) {
  
  si_vector *Vect;
  signed int *p;
  int i;

  /* Deal with pathological parameter */
  if (NbElems < 0)
  {
    si_vector_errormsg("si_vector_alloc", "invalid parameter", "NbRows must be >= 0");
    return 0;
  }
  /* Create the vector struct */
  Vect=(si_vector *)malloc(sizeof(si_vector));
  if(!Vect) {	
    si_vector_errormsg("si_vector_alloc", "outofmem", "out of memory space");
    return 0;
  }
  Vect->NbElems = NbElems;
  /* NbRows == 0: created an "empty" vector. */
  if(NbElems == 0) 
  {
    Vect->p_Init= (signed int *)0;
  }  
  else 
  {
    p = (signed int *)malloc(NbElems * sizeof(signed int));
    if(!p) 
    {
      free(Vect);
      si_vector_errormsg("si_vector_alloc", "outofmem", "out of memory space");
      return 0;
    }
    Vect->p_Init = p;
    for (i = 0 ; i < NbElems ; i++) 
      {
        *(p+i) = 0;
      }
  }
  Vect->p_Init_size = NbElems;
  return Vect;
} /* si_vector_alloc */


/**
 * Create a vector from an array.
 */

si_vector* 
si_vector_from_array(signed int* si_array, unsigned int NbElems)
{
  si_vector* new_vector;
  if (NbElems == 0)
  {
    return((si_vector*)0);
  }
  new_vector = si_vector_alloc(0);
  if (! new_vector)
  {
    si_vector_errormsg("si_vector_from_array", "outofmem", "out of memory space");
    return((si_vector*)0);
  }  
  new_vector->p_Init = si_array;
  new_vector->p_Init_size = NbElems;
  new_vector->NbElems = NbElems;
  return(new_vector);
} /* End si_vector_from_array */

 
/**
 * Free the memory space occupied by si_vector 'Vect' 
 */
void si_vector_free(si_vector *Vect) 
{
  if (Vect->p_Init)
  {
    free(Vect->p_Init);
  }
  free(Vect);
} /* si_vector_Free */

/* 
 * Print the contents of the si_vector 'Vect'.
 */
void si_vector_print(FILE *Dst, 
                        si_vector *Vect) 
{
  signed int *p;
  int i;
  unsigned NbElems;
  
  fprintf(Dst,"%d\n", NbElems=Vect->NbElems);
  p= Vect->p_Init;
  for (i = 0 ; i < NbElems ; i++) 
  {
    fprintf(Dst, "%d", *p++);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* si_vector_print */

/* 
 * Print the contents of the si_vector 'Vect' without printing the
 * the dimensions.
 */
void si_vector_print_no_dims(FILE *Dst, 
                                si_vector *Vect) 
{
  
  signed int *p;
  int i;
  p= Vect->p_Init;
  for (i = 0 ; i < Vect->p_Init_size ; i++) 
  {
    fprintf(Dst, "%d", *p++);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* si_vector_print_no_dims */

/* 
 * Print the contents of the si_vector 'Vect' in the Matlab
 * input format.
 */
void si_vector_print_matlab(FILE *Dst, 
                                si_vector *Vect) 
{
  signed int *p;
  int i;
  p= Vect->p_Init;
  fprintf(Dst, "[ ");
  for (i=0 ; i < Vect->p_Init_size ; i++) 
  {
    fprintf(Dst, "%d", *p++);
    if (i < (Vect->p_Init_size - 1))
    {
      fprintf(Dst, ", ");
    }
    fprintf(Dst, " ]\n");
  }
} /* si_vector_print_matlab */


/* 
 * Read the contents of the si_vector 'Vect' expressed as
 * strings convertible into signed int.
 */
void si_vector_read_input_from_file(si_vector *Vect, 
                                      int radix, 
                                      FILE* Src) 
{
  signed int *p;
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
    si_vector_errormsg("si_vector_read_rat_input_from_file", 
                          "baddim", 
                          "no data" );
  }
  for (i = 0 ; i < Vect->NbElems ; i++) 
  {
    if(!c || *c=='\n' || *c=='#') 
    {
      si_vector_errormsg("si_vector_read_input_from_file", 
                            "baddim", 
                            "not enough elements");
      break;
    }
    if (sscanf(c,"%s%n", str, &n) == 0) 
    {
	     si_vector_errormsg( "si_vector_Read", "baddim", "not enough elements" );
	     break;
    }
      *(p++) = atoi(str);
      c += n;
    }
} /* si_vector_read_input_from_file */

/* 
 * Read the contents of the vector 'Vect' from standard input. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the vector contents, one row per line.
 */
si_vector *si_vector_read_from_stdin(int radix) 
{
  
  si_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  while(fgets(s, 1024, stdin)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, stdin);
  Vect = si_vector_alloc(NbElems);
  if(!Vect) {
    si_vector_errormsg("si_vector_read", "outofmem", "out of memory space");
    return(NULL);
  }
  si_vector_read_input_from_file(Vect, radix, stdin);
  return Vect;
} /* si_vector_read_from_stdin */

/* 
 * Read the contents of the vector 'Vect' from FILE pointed to by Src. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the vector contents, one row per line.
 * Vectrix contents is made of strings directly convertible into signed int.
 */
si_vector*
si_vector_read_from_file(int radix, 
                            FILE* Src) 
{
  si_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  while(fgets(s, 1024, Src)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, Src);
  Vect = si_vector_alloc(NbElems);
  if(!Vect) {
    si_vector_errormsg("si_vector_read_from_file", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  si_vector_read_input_from_file(Vect, radix, Src);
  return Vect;
} /* si_vector_read_from_file */


signed int*
si_vector_get_at(si_vector* vect, unsigned int pos)
{
  if (vect == NULL)
  {
    si_vector_errormsg("si_vector_get_at", "null pointer parameter", "mat parameter is NULL");
    return(NULL);
  }
  if (vect->p_Init_size == 0)
  {
    si_vector_errormsg("si_vector_get_at", "empty vector", "can't return an element from an empty vector");
    return(NULL);
  }
  if (pos >= vect->NbElems)
  {
    si_vector_errormsg("si_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  if (pos < 0)
  {
    si_vector_errormsg("si_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  return(&(vect->p_Init[pos]));
} /* End si_vector_get_at */

void 
si_vector_set_at(si_vector* vect, 
                    unsigned int pos, 
                    signed int value)
{
  if (vect == NULL)
  {
    si_vector_errormsg("si_vector_set_at", "null pointer parameter", "mat parameter is NULL");
  }
  if (vect->p_Init_size == 0)
  {
    si_vector_errormsg("si_vector_set_at", "empty vector", "can't set an element from an empty vector");
  }
  if (pos >= vect->NbElems)
  {
    si_vector_errormsg("si_vector_set_at", "index out of range", "pos index beyond bounds");
  }
  vect->p_Init[pos] = value;
} /* End si_vector_set_at */

unsigned int
si_vector_get_NbElems(si_vector* vect)
{
  return(vect->NbElems);
}

unsigned int
si_vector_get_size(si_vector* vect)
{
  return(si_vector_get_NbElems(vect));
} /* End si_vector_get_size */

/**
 * @see si-vector.h#si_vector_has_duplicates
 */
int
si_vector_has_duplicates(si_vector* vect)
{
  unsigned int nb_elems;
  unsigned int i, j;
  int int_temp;
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    return(0);
  }
  if ((nb_elems = si_vector_get_size(vect)) == 0)
  {
    return(0);
  }
  for (i = 0 ; i < (nb_elems - 1) ; i++)
  {
    int_temp = *(si_vector_get_at(vect, i));
    for (j = i + 1  ; j < nb_elems ; j++)
    {
      if (int_temp == *(si_vector_get_at(vect, j)))
      {
        return(1);
      }
    } /* End for j. */
  } /* End for i. */
  return(0);
} /* End si_vector_has_duplicates */

/**
 * @see si-vector.h#si_vector_get_abs
 */
si_vector*
si_vector_get_abs(si_vector* vect)
{
  si_vector* abs_vect;
  unsigned int num_elems;
  unsigned int i;
  
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    si_vector_errormsg("si_vector_get_abs", 
                          "null pointer parameter", 
                          "vect parameter is NULL");
    return(NULL);
  }
  num_elems = si_vector_get_size(vect);
  if (num_elems == 0)
  {
    si_vector_errormsg("si_vector_get_abs", 
                          "empty vector", 
                          "vect parameter is an empty vector");
    return(NULL);
  }
  abs_vect = si_vector_alloc(num_elems);
  if (abs_vect == NULL)
  {
    si_vector_errormsg("si_vector_get_abs", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  for (i = 0 ; i < num_elems ; i++)
  {
    si_vector_set_at(abs_vect, i, abs(*(si_vector_get_at(vect, i))));
  }
  return(abs_vect);
} /* End mpq_vector_get_abs */

void
si_vector_shift(si_vector* vect, 
      unsigned int index, 
      int direction)
{
  /* Deal with pathological parameters. */
  /* Do nothing if the pointer is NULL. */
  if (vect == NULL)
  {
    return;
  }
  /* Do nothing if the vect is empty. */
  if (si_vector_get_size(vect) == 0)
  {
    return;
  }
  /* Do nothing if requested to shift up beyond the upper index */
  if ((direction == SI_VECTOR_UP) && (index + 2 > si_vector_get_size(vect)))
  {
    return;
  }
  /* Do nothing if requested to shift down beyond the lower index */
  if ((direction == SI_VECTOR_DOWN) && (index == 0))
  {
    return;
  }
  /* Do some usefull stuff */
  if (direction == SI_VECTOR_UP)
  {
    si_vector_set_at(vect,
                      index + 1,
                      *(si_vector_get_at(vect, index)));
  }
  else
  {
    si_vector_set_at(vect,
                      index - 1,
                      *(si_vector_get_at(vect, index)));
  }
  return;
} /* End si_vector_shift */

/**
 * see si-vector.c#si_vector_insertion_sort
 */
void
si_vector_insertion_sort(si_vector* vect, int direction)
{
  signed int i, j;
  unsigned int size;
  signed int temp;
  /* Deal with pathological parameters. */
  if (vect == NULL)
  {
    si_vector_errormsg("si_vector_insertion_sort", 
                        "null pointer parameter",
                        "vect is a NULL pointer");
    return;
  }
  size = si_vector_get_size(vect);
  if ( size == 0)
  {
    si_vector_errormsg("si_vector_insertion_sort", 
                        "empty vector",
                        "vect is a an empty vector");
    return;
  }
  if (direction == SI_VECTOR_UP)
  {
    for( i = 1 ; i < size ; i++)
    {
      temp = *(si_vector_get_at(vect, i));
      for( j = i - 1 ; j >= 0 && (*(si_vector_get_at(vect, j))> temp) ; j-- )
      {
        si_vector_shift(vect, j, SI_VECTOR_UP);
      } /* End for j */
      si_vector_set_at(vect,
                        j + 1,
                        temp);
    } /* End for i */
  }
  else /* direction == SI_VECTOR_DOWN */
  {
    for( i = 1 ; i < size ; i++)
    {
      temp = *(si_vector_get_at(vect, i));
      for( j = i - 1 ; j >= 0 && (*(si_vector_get_at(vect, j))< temp) ; j-- )
      {
        si_vector_shift(vect, j, SI_VECTOR_UP);
      } /* End for j */
      si_vector_set_at(vect,
                        j + 1,
                        temp);
    } /* End for i */
  }
} /* End si_vector_insertion_sort */

/**
 * see si-vector.h#si_vector_are_equal
 */
int 
si_vector_are_equal(si_vector* vect_1, si_vector* vect_2)
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
    si_vector_errormsg("si_vector_are_equal", 
                        "null pointer parameter", 
                        "one of the vectors (and only one) is NULL");
    return(0);
  }
  nb_elems = si_vector_get_size(vect_1);
  if ((nb_elems == 0) && (si_vector_get_size(vect_2) == 0))
  {
    return(1);
  }
  if (nb_elems != si_vector_get_size(vect_2))
  {
    return(0);
  }
  /* Actual check for equality. */
  for (i = 0 ; i < nb_elems ; i++)
  {
    if (*(si_vector_get_at(vect_1, i)) != *(si_vector_get_at(vect_2, i)))
    {
      return(0);
    }
  } /* End for. */
  return(1);
} /* End si_vector_are_equal */

/**
 * see si-vector.h#si_vector_clone
 */
si_vector*
si_vector_clone(si_vector* orig_vect)
{
  si_vector*    new_vect = NULL;
  unsigned int  nb_elems = 0;
  unsigned int  i;
  /* Deal with pathological parameters */
  if (orig_vect == NULL)
  {
    si_vector_errormsg("si_vector_clone", 
                        "null pointer parameter", 
                        "orig_vect parameter is NULL");
    return(NULL);
  }
  /* Deal with the "classical" case. */
  nb_elems = si_vector_get_size(orig_vect);
  new_vect = si_vector_alloc(nb_elems);
  if (new_vect == NULL)
  {
    si_vector_errormsg("si_vector_alloc", "outofmem", "out of memory space");
    return(NULL);
  }
  for (i = 0; i < nb_elems ; i++)
  {
    si_vector_set_at(new_vect,
                      i,
                      *(si_vector_get_at(orig_vect, i)));
  }
  return(new_vect);
} /* End si_vector_clone */


int
si_vector_vector_product(si_vector* v1, si_vector* v2, signed int* result)
{
  int i = 0;
  signed int part_prod;
  signed int* v1_p_Init;
  signed int* v2_p_Init;
  /* Check for incompatible vector sizes. */
  if (v1->p_Init_size != v2->p_Init_size)
  {
    si_vector_errormsg("si_vector_vector_product", "vector size conflict", "vectors have different sizes");
    return(1);
  }
  v1_p_Init = v1->p_Init;
  v2_p_Init = v2->p_Init;
  *result = 0;
  for (i = 0 ; i < v1->p_Init_size ; i++)
  {
    part_prod = *v1_p_Init * *v2_p_Init;
    *result  += part_prod;
    v1_p_Init++;
    v2_p_Init++;
  } /* End for */
  return(0);
} /* End si_vector_vector_product */

void 
si_vector_errormsg(char *f , char *msgname, char *msg) 
{
  fprintf(stderr, "?%s: %s\n", f, msg);
}


