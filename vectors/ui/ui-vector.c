/* ui_vector.c */

#include <stdio.h>
#include <stdlib.h>
/*#include <string.>*/
#include <ctype.h>

#include "ui-vector.h"

void
ui_vector_errormsg(char *f , char *msgname, char *msg);

/* 
 * Allocate space for vector dimensioned by 'NbElems'.
 */
ui_vector* ui_vector_alloc(unsigned int NbElems) {
  
  ui_vector *Vect;
  unsigned int *p;
  int i;

  /* Deal with pathological parameter */
  if (NbElems < 0)
  {
    ui_vector_errormsg("ui_vector_alloc", "invalid parameter", "NbRows must be >= 0");
    return 0;
  }
  /* Create the vector struct */
  Vect=(ui_vector *)malloc(sizeof(ui_vector));
  if(!Vect) {	
    ui_vector_errormsg("ui_vector_Alloc", "outofmem", "out of memory space");
    return 0;
  }
  Vect->NbElems = NbElems;
  /* NbRows == 0: created an "empty" vector. */
  if(NbElems == 0) 
  {
    Vect->p_Init= (unsigned int *)0;
  }  
  else 
  {
    p = (unsigned int *)malloc(NbElems * sizeof(unsigned int));
    if(!p) 
    {
      free(Vect);
      ui_vector_errormsg("ui_vector_Alloc", "outofmem", "out of memory space");
      return 0;
    }
    Vect->p_Init = p;
    for (i = 0 ; i < NbElems ; i++) 
      {
        *(p+i) = 0;
      }
  }
  return Vect;
} /* ui_vector_alloc */


/**
 * Create a vector from an array.
 */

ui_vector* 
ui_vector_from_array(unsigned int* ui_array, unsigned int NbElems)
{
  ui_vector* new_vector;
  int i;
  if (NbElems == 0)
  {
    return((ui_vector*)0);
  }
  new_vector = ui_vector_alloc(NbElems);
  if (! new_vector)
  {
    ui_vector_errormsg("ui_vector_from_array", "outofmem", "out of memory space");
    return((ui_vector*)0);
  }
  for (i = 0 ; i < NbElems ; i++)
  {
    (new_vector->p_Init)[i] = ui_array[i];
  }
  return(new_vector);
} /* End ui_vector_from_array */

 
/**
 * Free the memory space occupied by ui_vector 'Vect' 
 */
void ui_vector_free(ui_vector *Vect) 
{
  if (Vect->p_Init)
  {
    free(Vect->p_Init);
  }
  free(Vect);
} /* ui_vector_Free */

/* 
 * Print the contents of the ui_vector 'Vect'.
 */
void ui_vector_print(FILE *Dst, 
                        ui_vector *Vect) 
{
  unsigned int *p;
  int i;
  unsigned NbElems;
  
  fprintf(Dst,"%u\n", NbElems=Vect->NbElems);
  p= Vect->p_Init;
  for (i = 0 ; i < NbElems ; i++) 
  {
    fprintf(Dst, "%u", *p++);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* ui_vector_print */

/* 
 * Print the contents of the ui_vector 'Vect' without printing the
 * the dimensions.
 */
void ui_vector_print_no_dims(FILE *Dst, 
                              ui_vector *Vect) 
{
  
  unsigned int *p;
  int i;
  p= Vect->p_Init;
  for (i = 0 ; i < Vect->NbElems ; i++) 
  {
    fprintf(Dst, "%u", *p++);
    fprintf(Dst, " ");
  }
    fprintf(Dst, "\n");
} /* ui_vector_print_no_dims */

/* 
 * Print the contents of the ui_vector 'Vect' in the Matlab
 * input format.
 */
void ui_vector_print_matlab(FILE *Dst, 
                            ui_vector *Vect) 
{
  unsigned int *p;
  int i;
  p= Vect->p_Init;
  fprintf(Dst, "[ ");
  for (i=0 ; i < Vect->NbElems ; i++) 
  {
    fprintf(Dst, "%u", *p++);
    if (i < (Vect->NbElems - 1))
    {
      fprintf(Dst, ", ");
    }
    fprintf(Dst, " ]\n");
  }
} /* ui_vector_print_matlab */


/* 
 * Read the contents of the ui_vector 'Vect' expressed as
 * strings convertible into signed int.
 */
void ui_vector_read_input_from_file(ui_vector *Vect, 
                                      int radix, 
                                      FILE* Src) 
{
  unsigned int *p;
  int i,n;
  char *c, **endptr, s[1024],str[1024];
  p = Vect->p_Init;
  do 
  {
    c = fgets(s, 1024, Src);
    while(isspace(*c) && *c != '\n')
	   ++c;
  } while(c && (*c == '#' || *c == '\n'));
  if (!c) 
  {
    ui_vector_errormsg("ui_vector_read_rat_input_from_file", 
                          "baddim", 
                          "no data" );
  }
  for (i = 0 ; i < Vect->NbElems ; i++) 
  {
    if(!c || *c=='\n' || *c=='#') 
    {
      ui_vector_errormsg("ui_vector_read_input_from_file", 
                          "baddim", 
                          "not enough elements");
      break;
    }
    if (sscanf(c,"%s%n", str, &n) == 0) 
    {
      ui_vector_errormsg( "ui_vector_Read", "baddim", "not enough elements" );
	    break;
    }
      *(p++) = strtoul(str, endptr, radix);
      if (endptr != NULL)
      {
	     ui_vector_errormsg( "ui_vector_Read", "wrongFormat", "string can not be converted into unsigned int" );
	     break;
      }
      c += n;
    }
} /* ui_vector_read_input_from_file */

/* 
 * Read the contents of the vector 'Vect' from standard input. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the vector contents, one row per line.
 */
ui_vector *ui_vector_read_from_stdin(int radix) 
{
  
  ui_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  while(fgets(s, 1024, stdin)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, stdin);
  Vect = ui_vector_alloc(NbElems);
  if(!Vect) {
    ui_vector_errormsg("ui_vector_read", "outofmem", "out of memory space");
    return(NULL);
  }
  ui_vector_read_input_from_file(Vect, radix, stdin);
  return Vect;
} /* ui_vector_read_from_stdin */

/* 
 * Read the contents of the vector 'Vect' from FILE pointed to by Src. 
 * A '#' in the first column denotes a comment line.
 * The first line holds the number of rows and columns, subsequent
 * lines hold the vector contents, one row per line.
 * Vectrix contents is made of strings directly convertible into signed int.
 */
ui_vector*
ui_vector_read_from_file(int radix, 
                            FILE* Src) 
{
  ui_vector *Vect;
  unsigned NbElems;
  char s[1024];
  
  while(fgets(s, 1024, Src)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d", &NbElems) < 1))
    fgets(s, 1024, Src);
  Vect = ui_vector_alloc(NbElems);
  if(!Vect) {
    ui_vector_errormsg("ui_vector_read_from_file", 
                          "outofmem", 
                          "out of memory space");
    return(NULL);
  }
  ui_vector_read_input_from_file(Vect, radix, Src);
  return Vect;
} /* ui_vector_read_from_file */


unsigned int*
ui_vector_get_at(ui_vector* vect, unsigned int pos)
{
  if (vect == NULL)
  {
    ui_vector_errormsg("ui_vector_get_at", "null pointer parameter", "mat parameter is NULL");
    return(NULL);
  }
  if (vect->NbElems == 0)
  {
    ui_vector_errormsg("ui_vector_get_at", "empty vector", "can't return an element from an empty vector");
    return(NULL);
  }
  if (pos >= vect->NbElems)
  {
    ui_vector_errormsg("ui_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  if (pos < 0)
  {
    ui_vector_errormsg("ui_vector_get_at", "index out of range", "pos index beyond bounds");
    return(NULL);
  }
  return(&(vect->p_Init[pos]));
} /* End ui_vector_get_at */

void 
ui_vector_set_at(ui_vector* vect, 
                    unsigned int pos, 
                    unsigned int value)
{
  if (vect == NULL)
  {
    ui_vector_errormsg("ui_vector_set_at", "null pointer parameter", "mat parameter is NULL");
  }
  if (vect->NbElems == 0)
  {
    ui_vector_errormsg("ui_vector_set_at", "empty vector", "can't set an element from an empty vector");
  }
  if (pos >= vect->NbElems)
  {
    ui_vector_errormsg("ui_vector_set_at", "index out of range", "pos index beyond bounds");
  }
  vect->p_Init[pos] = value;
} /* End ui_vector_set_at */

unsigned int
ui_vector_get_size(ui_vector* vect)
{
  return(vect->NbElems);
} /* End ui_vector_get_size */

/**
 * @see ui-vector.h#ui_vector_has_duplicates
 */
int
ui_vector_has_duplicates(ui_vector* vect)
{
  unsigned int nb_elems;
  unsigned int i, j;
  int int_temp;
  /* Deal with pathological parameters */
  if (vect == NULL)
  {
    return(0);
  }
  if ((nb_elems = ui_vector_get_size(vect)) == 0)
  {
    return(0);
  }
  for (i = 0 ; i < (nb_elems - 1) ; i++)
  {
    int_temp = *(ui_vector_get_at(vect, i));
    for (j = i + 1  ; j < nb_elems ; j++)
    {
      if (int_temp == *(ui_vector_get_at(vect, j)))
      {
        return(1);
      }
    } /* End for j. */
  } /* End for i. */
  return(0);
} /* End ui_vector_has_duplicates */

/**
 * @see ui-vector.h#ui_vector_shift
 */
void
ui_vector_shift(ui_vector* vect, 
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
  if (ui_vector_get_size(vect) == 0)
  {
    return;
  }
  /* Do nothing if requested to shift up beyond the upper index */
  if ((direction == UI_VECTOR_UP) && (index + 2 > ui_vector_get_size(vect)))
  {
    return;
  }
  /* Do nothing if requested to shift down beyond the lower index */
  if ((direction == UI_VECTOR_DOWN) && (index == 0))
  {
    return;
  }
  /* Do some usefull stuff */
  if (direction == UI_VECTOR_UP)
  {
    ui_vector_set_at(vect,
                      index + 1,
                      *(ui_vector_get_at(vect, index)));
  }
  else
  {
    ui_vector_set_at(vect,
                      index - 1,
                      *(ui_vector_get_at(vect, index)));
  }
  return;
} /* End ui_vector_shift */

/**
 * see si-vector.c#si_vector_insertion_sort
 */
void
ui_vector_insertion_sort(ui_vector* vect, int direction)
{
  signed int i, j;
  unsigned int size;
  signed int temp;
  /* Deal with pathological parameters. */
  if (vect == NULL)
  {
    ui_vector_errormsg("ui_vector_insertion_sort", 
                        "null pointer parameter",
                        "vect is a NULL pointer");
    return;
  }
  size = ui_vector_get_size(vect);
  if ( size == 0)
  {
    ui_vector_errormsg("ui_vector_insertion_sort", 
                        "empty vector",
                        "vect is a an empty vector");
    return;
  }
  if (direction == UI_VECTOR_UP)
  {
    for( i = 1 ; i < size ; i++)
    {
      temp = *(ui_vector_get_at(vect, i));
      for( j = i - 1 ; j >= 0 && (*(ui_vector_get_at(vect, j))> temp) ; j-- )
      {
        ui_vector_shift(vect, j, UI_VECTOR_UP);
      } /* End for j */
      ui_vector_set_at(vect,
                        j + 1,
                        temp);
    } /* End for i */
  }
  else /* direction == UI_VECTOR_DOWN */
  {
    for( i = 1 ; i < size ; i++)
    {
      temp = *(ui_vector_get_at(vect, i));
      for( j = i - 1 ; j >= 0 && (*(ui_vector_get_at(vect, j))< temp) ; j-- )
      {
        ui_vector_shift(vect, j, UI_VECTOR_UP);
      } /* End for j */
      ui_vector_set_at(vect,
                        j + 1,
                        temp);
    } /* End for i */
  }
} /* End ui_vector_insertion_sort */

int
ui_vector_scalar_product(ui_vector* v1, ui_vector* v2, unsigned int* result)
{
  int i = 0;
  signed int part_prod;
  unsigned int* v1_p_Init;
  unsigned int* v2_p_Init;
  /* Check for incompatible vector sizes. */
  if (v1->NbElems != v2->NbElems)
  {
    ui_vector_errormsg("ui_vector_vector_product", "vector size conflict", "vectors have different sizes");
    return(1);
  }
  v1_p_Init = v1->p_Init;
  v2_p_Init = v2->p_Init;
  *result = 0;
  for (i = 0 ; i < v1->NbElems ; i++)
  {
    part_prod = *v1_p_Init * *v2_p_Init;
    *result  += part_prod;
    v1_p_Init++;
    v2_p_Init++;
  } /* End for */
  return(0);
} /* End ui_vector_vector_product */

void 
ui_vector_errormsg(char *f , char *msgname, char *msg) 
{
  fprintf(stderr, "?%s: %s\n", f, msg);
}


