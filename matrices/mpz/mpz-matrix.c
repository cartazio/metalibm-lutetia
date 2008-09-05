/* mpz_matrix.c */

#include <stdio.h>
#include <stdlib.h>
/*#include <string.>*/
#include <ctype.h>
#include <gmp.h>
#include <mpfr.h>

#include "mpz-matrix.h"

void
mpz_matrix_errormsg(char *f , char *msgname, char *msg);

/* 
 * Allocate space for matrix dimensioned by 'NbRows X NbColumns'.
 */
mpz_matrix* mpz_matrix_alloc(unsigned NbRows,unsigned NbColumns) {
  
  mpz_matrix *Mat;
  mpz_t *p, **q;
  int i,j;

  /* Deal with pathological parameter */
  if (NbRows < 0)
  {
    mpz_matrix_errormsg("mpz_matrix_alloc", "invalid parameter", "NbRows must be >= 0");
    return 0;
  }
  if (NbRows < 0)
  {
    mpz_matrix_errormsg("mpz_matrix_alloc", "invalid parameter", "NbColumns must be >= 0");
    return 0;
  }
  /* Create the matrix struct */
  Mat=(mpz_matrix *)malloc(sizeof(mpz_matrix));
  if(!Mat) {	
    mpz_matrix_errormsg("mpz_matrix_Alloc", "outofmem", "out of memory space");
    return 0;
  }
  Mat->NbRows=NbRows;
  Mat->NbColumns=NbColumns;
  /* NbRows == 0: created an "empty" matrix. */
  if(NbRows==0) 
  {
    Mat->p = (mpz_t **)0;
    Mat->p_Init= (mpz_t *)0;
  }  
  else 
  {
    /* NbColumns == 0: created an "empty" matrix. */
    if(NbColumns==0) 
    {
      Mat->p = (mpz_t **)0;
      Mat->p_Init= (mpz_t *)0;
    }
    else 
    {
      /* Create the elements array that holds the matrix elements */
      q = (mpz_t **)malloc(NbRows * sizeof(*q));
      if(!q) 
      {
        free(Mat);
        mpz_matrix_errormsg("mpz_matrix_Alloc", "outofmem", "out of memory space");
        return 0;
      }
      p = (mpz_t *)malloc(NbRows * NbColumns * sizeof(mpz_t));
      if(!p) 
      {
        free(q);
        free(Mat);
        mpz_matrix_errormsg("mpz_matrix_Alloc", "outofmem", "out of memory space");
        return 0;
      }
      Mat->p = q;
      Mat->p_Init = p;
      for (i=0;i<NbRows;i++) 
      {
        *q++ = p;
        for (j=0;j<NbColumns;j++)   
          mpz_init(*(p+j));
        p += NbColumns;
      }
    }
  }
  p = NULL;
  q = NULL;
  Mat->p_Init_size = NbColumns*NbRows;
  return Mat;
} /* mpz_matrix_Alloc */

/* 
 * Free the memory space occupied by mpz_matrix 'Mat' 
 */
void mpz_matrix_free(mpz_matrix *Mat) {
  
  int i;
  mpz_t *p;

  if (Mat == NULL) return;

  if( Mat->p )
  {
    p = *(Mat->p);
    for(i=0;i<Mat->p_Init_size;i++) {
      mpz_clear(*p++);
    }
  }

  if (Mat->p_Init)
    free(Mat->p_Init);
  if (Mat->p)
    free(Mat->p);
  free(Mat);

} /* mpz_matrix_Free */

/* 
 * Print the contents of the mpz_matrix 'Mat'.
 */
void mpz_matrix_print(FILE *Dst, int radix ,mpz_matrix *Mat) {
  
  mpz_t *p;
  int i, j;
  unsigned NbRows, NbColumns;
  
  fprintf(Dst,"%d %d\n", NbRows=Mat->NbRows, NbColumns=Mat->NbColumns);
  for (i=0;i<NbRows;i++) {
    p=*(Mat->p+i);
    for (j=0;j<NbColumns;j++) 
    {
      mpz_out_str(Dst, radix, *p++);
      fprintf(Dst, " ");
    }
    fprintf(Dst, "\n");
  }
} /* mpz_matrix_print */

/* 
 * Print the contents of the mpz_matrix 'Mat' without printing the
 * the dimensions.
 */
void mpz_matrix_print_no_dims(FILE *Dst, int radix ,mpz_matrix *Mat) {
  
  mpz_t *p;
  int i, j;
  unsigned NbRows, NbColumns;
  
  NbRows    = Mat->NbRows; 
  NbColumns = Mat->NbColumns;
  for (i=0 ; i<NbRows ; i++) 
  {
    p=*(Mat->p+i);
    for (j=0;j<NbColumns;j++) 
    {
      mpz_out_str(Dst, radix, *p++);
      fprintf(Dst, " ");
    }
    fprintf(Dst, "\n");
  }
} /* mpz_matrix_print_no_dims */

/* 
 * Print the contents of the mpz_matrix 'Mat' in the Matlab
 * input format.
 */
void mpz_matrix_print_matlab(FILE *Dst, int radix ,mpz_matrix *Mat) {
  
  mpz_t *p;
  int i, j;
  unsigned NbRows, NbColumns;
  
  NbRows    = Mat->NbRows; 
  NbColumns = Mat->NbColumns;
  fprintf(Dst, "[ ");
  for (i=0 ; i < NbRows ; i++) 
  {
    p=*(Mat->p+i);
    for (j=0 ; j < NbColumns ; j++) 
    {
      mpz_out_str(Dst, radix, *p++);
      if (j < (NbColumns - 1))
      {
        fprintf(Dst, ", ");
      }
      else
      {
        if (i < (NbRows - 1))
        {
          fprintf(Dst, "; ");
        }
        else
        {
          fprintf(Dst, " ]\n");
        }
      }
    }
  }
} /* mpz_matrix_print_matlab */

/* 
 * Read the contents of the mpz_matrix 'Mat' 
 */
void mpz_matrix_read_input(mpz_matrix *Mat, int radix) 
{
  
  mpz_t *p;
  int i,j,n;
  char *c, s[1024],str[1024];
  
  p = Mat->p_Init;
  for (i=0;i<Mat->NbRows;i++) {
    do {
      c = fgets(s, 1024, stdin);
      while(isspace(*c) && *c!='\n')
	++c;
    } while(c && (*c =='#' || *c== '\n'));
    
    if (!c) {
      mpz_matrix_errormsg( "mpz_matrix_Read", "baddim", "not enough rows" );
      break;
    }
    for (j=0;j<Mat->NbColumns;j++) {
      if(!c || *c=='\n' || *c=='#') {
	mpz_matrix_errormsg("mpz_matrix_Read", "baddim", "not enough columns");
	break;
      }
      if (sscanf(c,"%s%n",str,&n) == 0) {
	mpz_matrix_errormsg( "mpz_matrix_Read", "baddim", "not enough columns" );
	break;
      }
      mpz_set_str(*(p++),str, radix);
      c += n;
    }
  }
} /* mpz_matrix_read_input */

/* 
 * Read the contents of the matrix 'Mat' from standard input. 
 * A '#' in the first column is a comment line 
 */
mpz_matrix *mpz_matrix_read(int radix) 
{
  
  mpz_matrix *Mat;
  unsigned NbRows, NbColumns;
  char s[1024];
  
  while(fgets(s, 1024, stdin)==0);
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, "%d %d", &NbRows, &NbColumns)<2))
    fgets(s, 1024, stdin);
  Mat = mpz_matrix_alloc(NbRows,NbColumns);
  if(!Mat) {
    mpz_matrix_errormsg("mpz_matrix_read", "outofmem", "out of memory space");
    return(NULL);
  }
  mpz_matrix_read_input(Mat, radix);
  return Mat;
} /* mpz_matrix_read */

mpz_t*
mpz_matrix_get_at(mpz_matrix* mat, unsigned int row, unsigned int column)
{
  mpz_t** matRow;
  if (mat == NULL)
  {
    mpz_matrix_errormsg("mpz_matrix_get_at", "null pointer parameter", "mat parameter is NULL");
    return(NULL);
  }
  if (mat->p_Init_size == 0)
  {
    mpz_matrix_errormsg("mpz_matrix_get_at", "empty matrix", "can't return an element from an empty matrix");
    return(NULL);
  }
  if (row >= mat->NbRows)
  {
    mpz_matrix_errormsg("mpz_matrix_get_at", "index out of range", "row index beyond bounds");
    return(NULL);
  }
  if (column >= mat->NbColumns)
  {
    mpz_matrix_errormsg("mpz_matrix_get_at", "index out of range", "column index beyond bounds");
    return(NULL);
  }
  matRow = mat->p + row; 
  return(*matRow + column);
} /* End mpz_matrix_get_at */

void 
mpz_matrix_set_at(mpz_matrix* mat, unsigned int row, unsigned int column, mpz_t value)
{
  mpz_t** matRow;
  if (mat == NULL)
  {
    mpz_matrix_errormsg("mpz_matrix_set_at", "null pointer parameter", "mat parameter is NULL");
  }
  if (mat->p_Init_size == 0)
  {
    mpz_matrix_errormsg("mpz_matrix_set_at", "empty matrix", "can't set an element from an empty matrix");
  }
  if (row >= mat->NbRows)
  {
    mpz_matrix_errormsg("mpz_matrix_set_at", "index out of range", "row index beyond bounds");
  }
  if (column >= mat->NbColumns)
  {
    mpz_matrix_errormsg("mpz_matrix_set_at", "index out of range", "column index beyond bounds");
  }
  matRow = mat->p + row; 
  mpz_set(*(*matRow + column), value);
} /* End mpz_matrix_set_at */

void 
mpz_matrix_errormsg(char *f , char *msgname, char *msg) 
{
  fprintf(stderr, "?%s: %s\n", f, msg);
}

#if 0
/* 
 * Basic hermite engine 
 */
static int hermite(mpz_matrix *H,mpz_matrix *U,mpz_matrix *Q) {
  
  int nc, nr, i, j, k, rank, reduced, pivotrow;
  mpz_t pivot,x,aux;
  mpz_t *temp1, *temp2;
  
  /*                     T                     -1   T */
  /* Computes form: A = Q H  and U A = H  and U  = Q  */
  
  if (!H) { 
    mpz_matrix_errormsg("Domlib", "nullH", "hermite: ? Null H");
    return -1;
  }
  nc = H->NbColumns;
  nr = H->NbRows;
  temp1 = (mpz_t *) malloc(nc * sizeof(mpz_t));
  temp2 = (mpz_t *) malloc(nr * sizeof(mpz_t));
  if (!temp1 ||!temp2) {
    mpz_matrix_errormsg("Domlib", "outofmem", "out of memory space");
    return -1;
  }
  
  /* Initialize all the 'mpz_t' variables */
  value_init(pivot); value_init(x); 
  value_init(aux);   
  for(i=0;i<nc;i++)
    value_init(temp1[i]);
  for(i=0;i<nr;i++)
    value_init(temp2[i]);
  
#ifdef DEBUG
  fprintf(stderr,"Start  -----------\n");
  mpz_matrix_Print(stderr,0,H);
#endif
  for (k=0, rank=0; k<nc && rank<nr; k=k+1) {
    reduced = 1;	/* go through loop the first time */
#ifdef DEBUG
    fprintf(stderr, "Working on col %d.  Rank=%d ----------\n", k+1, rank+1);
#endif
    while (reduced) {
      reduced=0;
      
      /* 1. find pivot row */
      value_absolute(pivot,H->p[rank][k]);
      
      /* the kth-diagonal element */
      pivotrow = rank;
      
      /* find the row i>rank with smallest nonzero element in col k */
      for (i=rank+1; i<nr; i++) {
	value_absolute(x,H->p[i][k]);
	if (value_notzero_p(x) &&
	    (value_lt(x,pivot) || value_zero_p(pivot))) {
	  value_assign(pivot,x);
	  pivotrow = i;
	}
      }
      
      /* 2. Bring pivot to diagonal (exchange rows pivotrow and rank) */
      if (pivotrow != rank) {
	Vector_Exchange(H->p[pivotrow],H->p[rank],nc);
	if (U)
	  Vector_Exchange(U->p[pivotrow],U->p[rank],nr);
	if (Q)
	  Vector_Exchange(Q->p[pivotrow],Q->p[rank],nr);

#ifdef DEBUG
	fprintf(stderr,"Exchange rows %d and %d  -----------\n", rank+1, pivotrow+1);
	mpz_matrix_Print(stderr,0,H);
#endif
      }
      value_assign(pivot,H->p[rank][k]);	/* actual ( no abs() ) pivot */
      
      /* 3. Invert the row 'rank' if pivot is negative */
      if (value_neg_p(pivot)) {
	value_oppose(pivot,pivot); /* pivot = -pivot */
	for (j=0; j<nc; j++)
	  value_oppose(H->p[rank][j],H->p[rank][j]);
	
	/* H->p[rank][j] = -(H->p[rank][j]); */
	if (U)
	  for (j=0; j<nr; j++)
	    value_oppose(U->p[rank][j],U->p[rank][j]);
	
	/* U->p[rank][j] = -(U->p[rank][j]); */
	if (Q)
	  for (j=0; j<nr; j++)
	    value_oppose(Q->p[rank][j],Q->p[rank][j]);
	
	/* Q->p[rank][j] = -(Q->p[rank][j]); */
#ifdef DEBUG
	fprintf(stderr,"Negate row %d  -----------\n", rank+1);
	mpz_matrix_Print(stderr,0,H);
#endif

      }      
      if (value_notzero_p(pivot)) {
	
	/* 4. Reduce the column modulo the pivot */
	/*    This eventually zeros out everything below the */
	/*    diagonal and produces an upper triangular matrix */
	
	for (i=rank+1;i<nr;i++) {
	  value_assign(x,H->p[i][k]);
	  if (value_notzero_p(x)) {	    
	    value_modulus(aux,x,pivot);
	    
	    /* floor[integer division] (corrected for neg x) */
	    if (value_neg_p(x) && value_notzero_p(aux)) {
	      
	      /* x=(x/pivot)-1; */
	      value_division(x,x,pivot);
	      value_decrement(x,x);
	    }	
	    else 
	      value_division(x,x,pivot);
	    for (j=0; j<nc; j++) {
	      value_multiply(aux,x,H->p[rank][j]);
	      value_substract(H->p[i][j],H->p[i][j],aux);
	    }
	    
	    /* U->p[i][j] -= (x * U->p[rank][j]); */
	    if (U)
	      for (j=0; j<nr; j++) {
		value_multiply(aux,x,U->p[rank][j]);
		value_substract(U->p[i][j],U->p[i][j],aux);
	      }
	    
	    /* Q->p[rank][j] += (x * Q->p[i][j]); */
	    if (Q)
	      for(j=0;j<nr;j++) {
		value_multiply(aux,x,Q->p[i][j]);
		value_addto(Q->p[rank][j],Q->p[rank][j],aux);
	      }
	    reduced = 1;

#ifdef DEBUG
	    fprintf(stderr,
		    "row %d = row %d - %d row %d -----------\n", i+1, i+1, x, rank+1);
	    mpz_matrix_Print(stderr,0,H);
#endif
	
	  } /* if (x) */
	} /* for (i) */
      } /* if (pivot != 0) */
    } /* while (reduced) */
    
    /* Last finish up this column */
    /* 5. Make pivot column positive (above pivot row) */
    /*    x should be zero for i>k */
    
    if (value_notzero_p(pivot)) {
      for (i=0; i<rank; i++) {
	value_assign(x,H->p[i][k]);
	if (value_notzero_p(x)) { 	  
	  value_modulus(aux,x,pivot);
	  
	  /* floor[integer division] (corrected for neg x) */
	  if (value_neg_p(x) && value_notzero_p(aux)) {
	    value_division(x,x,pivot);
	    value_decrement(x,x);
	    
	    /* x=(x/pivot)-1; */
	  }
	  else
	    value_division(x,x,pivot);
	  
	  /* H->p[i][j] -= x * H->p[rank][j]; */
	  for (j=0; j<nc; j++) {
	    value_multiply(aux,x,H->p[rank][j]);
	    value_substract(H->p[i][j],H->p[i][j],aux);
	  }
	  
	  /* U->p[i][j] -= x * U->p[rank][j]; */
	  if (U)
	    for (j=0; j<nr; j++) {
	      value_multiply(aux,x,U->p[rank][j]);
	      value_substract(U->p[i][j],U->p[i][j],aux);
	    }
	  
	  /* Q->p[rank][j] += x * Q->p[i][j]; */
	  if (Q)
	    for (j=0; j<nr; j++) {
	      value_multiply(aux,x,Q->p[i][j]);
	      value_addto(Q->p[rank][j],Q->p[rank][j],aux);
	    }  
#ifdef DEBUG
	  fprintf(stderr,
		  "row %d = row %d - %d row %d -----------\n", i+1, i+1, x, rank+1);
	  mpz_matrix_Print(stderr,0,H);
#endif
	} /* if (x) */
      } /* for (i) */
      rank++;
    } /* if (pivot!=0) */
  } /* for (k) */
  
  /* Clear all the 'mpz_t' variables */
  value_clear(pivot); value_clear(x); 
  value_clear(aux); 
  for(i=0;i<nc;i++)
    value_clear(temp1[i]);
  for(i=0;i<nr;i++)
    value_clear(temp2[i]);
  free(temp2);
  free(temp1);
  return rank;
} /* Hermite */ 

void right_hermite(mpz_matrix *A,mpz_matrix **Hp,mpz_matrix **Up,mpz_matrix **Qp) {
  
  mpz_matrix *H, *Q, *U;
  int i, j, nr, nc, rank;
  mpz_t tmp;
  
  /* Computes form: A = QH , UA = H */  
  nc = A->NbColumns;
  nr = A->NbRows;
  
  /* H = A */
  *Hp = H = mpz_matrix_Alloc(nr,nc);
  if (!H) { 
    mpz_matrix_errormsg("DomRightHermite", "outofmem", "out of memory space");
    return;
  }
  
  /* Initialize all the 'mpz_t' variables */
  value_init(tmp);
  
  Vector_Copy(A->p_Init,H->p_Init,nr*nc);
  
  /* U = I */
  if (Up) {
    *Up = U = mpz_matrix_Alloc(nr, nr);
    if (!U) {
      mpz_matrix_errormsg("DomRightHermite", "outofmem", "out of memory space");
      value_clear(tmp);
      return;
    }
    Vector_Set(U->p_Init,0,nr*nr);             /* zero's */
    for(i=0;i<nr;i++)                          /* with diagonal of 1's */
      value_set_si(U->p[i][i],1);
  }
  else
    U = (mpz_matrix *)0;
  
  /* Q = I */
  /* Actually I compute Q transpose... its easier */
  if (Qp) {
    *Qp = Q = mpz_matrix_Alloc(nr,nr);
    if (!Q) {
      mpz_matrix_errormsg("DomRightHermite", "outofmem", "out of memory space");
      value_clear(tmp);
      return;
    }
    Vector_Set(Q->p_Init,0,nr*nr);            /* zero's */
    for (i=0;i<nr;i++)                      /* with diagonal of 1's */
      value_set_si(Q->p[i][i],1);
  }
  else
    Q = (mpz_matrix *)0;
  
  rank = hermite(H,U,Q);
  
  /* Q is returned transposed */ 
  /* Transpose Q */
  if (Q) {
    for (i=0; i<nr; i++) {
      for (j=i+1; j<nr; j++) {
	value_assign(tmp,Q->p[i][j]);
	value_assign(Q->p[i][j],Q->p[j][i] );
	value_assign(Q->p[j][i],tmp);
      }
    }
  }
  value_clear(tmp);
  return;
} /* right_hermite */

void left_hermite(mpz_matrix *A,mpz_matrix **Hp,mpz_matrix **Qp,mpz_matrix **Up) {
  
  mpz_matrix *H, *HT, *Q, *U;
  int i, j, nc, nr, rank;
  mpz_t tmp;
  
  /* Computes left form: A = HQ , AU = H , 
     T    T T    T T   T
     using right form  A  = Q H  , U A = H */
  
  nr = A->NbRows;
  nc = A->NbColumns;
  
  /* HT = A transpose */
  HT = mpz_matrix_Alloc(nc, nr);
  if (!HT) {
    mpz_matrix_errormsg("DomLeftHermite", "outofmem", "out of memory space");
    return;
  }
  value_init(tmp);
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      value_assign(HT->p[j][i],A->p[i][j]);
  
  /* U = I */
  if (Up) {
    *Up = U = mpz_matrix_Alloc(nc,nc);
    if (!U) {
      mpz_matrix_errormsg("DomLeftHermite", "outofmem", "out of memory space");
      value_clear(tmp);
      return;
    }
    Vector_Set(U->p_Init,0,nc*nc);            /* zero's */
    for (i=0;i<nc;i++)                        /* with diagonal of 1's */
      value_set_si(U->p[i][i],1);
  }
  else U=(mpz_matrix *)0;
  
  /* Q = I */
  if (Qp) {
    *Qp = Q = mpz_matrix_Alloc(nc, nc);
    if (!Q) {
      mpz_matrix_errormsg("DomLeftHermite", "outofmem", "out of memory space");
      value_clear(tmp);
      return;
    }
    Vector_Set(Q->p_Init,0,nc*nc);            /* zero's */
    for (i=0;i<nc;i++)                        /* with diagonal of 1's */
      value_set_si(Q->p[i][i],1);
  }
  else Q=(mpz_matrix *)0;
  rank = hermite(HT,U,Q);
  
  /* H = HT transpose */
  *Hp = H = mpz_matrix_Alloc(nr,nc);
  if (!H) {
    mpz_matrix_errormsg("DomLeftHermite", "outofmem", "out of memory space");
    value_clear(tmp);
    return;
  }
  for (i=0; i<nr; i++)
    for (j=0;j<nc;j++)
      value_assign(H->p[i][j],HT->p[j][i]);
  mpz_matrix_Free(HT);
  
  /* Transpose U */
  if (U) {
    for (i=0; i<nc; i++) {
      for (j=i+1; j<nc; j++) {
	value_assign(tmp,U->p[i][j]);
	value_assign(U->p[i][j],U->p[j][i] );
	value_assign(U->p[j][i],tmp);
      }
    }
  }
} /* left_hermite */

/*
 * Given a integer matrix 'Mat'(k x k), compute its inverse rational matrix 
 * 'MatInv' k x (k+1). The last column of each row in matrix MatInv is used 
 * to store the common denominator of the entries in a row. The output is 1,
 * if 'Mat' is non-singular (invertible), otherwise the output is 0. Note:: 
 * (1) mpz_matrix 'Mat' is modified during the inverse operation.
 * (2) mpz_matrix 'MatInv' must be preallocated before passing into this function.
 */
int MatInverse(mpz_matrix *Mat,mpz_matrix *MatInv ) {
  
  int i, k, j, c;
  mpz_t x, gcd, piv;
  mpz_t m1,m2;
  
  if(Mat->NbRows != Mat->NbColumns) {
   fprintf(stderr,"Trying to invert a non-square matrix !\n");
    return 0;
  }
  
  /* Initialize all the 'mpz_t' variables */
  value_init(x);  value_init(gcd); value_init(piv);
  value_init(m1); value_init(m2);

  k = Mat->NbRows; 

  /* Initialise MatInv */
  Vector_Set(MatInv->p[0],0,k*(k+1));

  /* Initialize 'MatInv' to Identity matrix form. Each diagonal entry is set*/
  /* to 1. Last column of each row (denominator of each entry in a row) is  */
  /* also set to 1.                                                         */ 
  for(i=0;i<k;++i) {
    value_set_si(MatInv->p[i][i],1);	
    value_set_si(MatInv->p[i][k],1);	/* denum */
  }  
  /* Apply Gauss-Jordan elimination method on the two matrices 'Mat' and  */
  /* 'MatInv' in parallel.                                                */
  for(i=0;i<k;++i) {
    
    /* Check if the diagonal entry (new pivot) is non-zero or not */
    if(value_zero_p(Mat->p[i][i])) {   	
      
      /* Search for a non-zero pivot down the column(i) */
      for(j=i;j<k;++j)      
	if(value_notzero_p(Mat->p[j][i]))
	  break;
      
      /* If no non-zero pivot is found, the matrix 'Mat' is non-invertible */
      /* Return 0.                                                         */
      if(j==k) {
	
	/* Clear all the 'mpz_t' variables */
	value_clear(x);  value_clear(gcd); value_clear(piv);
	value_clear(m1); value_clear(m2);
	return 0;
      }	
      
      /* Exchange the rows, row(i) and row(j) so that the diagonal element */
      /* Mat->p[i][i] (pivot) is non-zero. Repeat the same operations on    */
      /* matrix 'MatInv'.                                                   */
      for(c=0;c<k;++c) {

	/* Interchange rows, row(i) and row(j) of matrix 'Mat'    */
	value_assign(x,Mat->p[j][c]);
	value_assign(Mat->p[j][c],Mat->p[i][c]);
	value_assign(Mat->p[i][c],x);
	
	/* Interchange rows, row(i) and row(j) of matrix 'MatInv' */
	value_assign(x,MatInv->p[j][c]);
	value_assign(MatInv->p[j][c],MatInv->p[i][c]);
	value_assign(MatInv->p[i][c],x);
      }
    }
    
    /* Make all the entries in column(i) of matrix 'Mat' zero except the */
    /* diagonal entry. Repeat the same sequence of operations on matrix  */
    /* 'MatInv'.                                                         */
    for(j=0;j<k;++j) {
      if (j==i) continue;	         /* Skip the pivot */
      value_assign(x,Mat->p[j][i]);
      if(value_notzero_p(x)) {
	value_assign(piv,Mat->p[i][i]);
	Gcd(x,piv,&gcd);
	if (value_notone_p(gcd) ) {
	  value_division(x,x,gcd);
	  value_division(piv,piv,gcd);
	}
	for(c=((j>i)?i:0);c<k;++c) {
	  value_multiply(m1,piv,Mat->p[j][c]);
	  value_multiply(m2,x,Mat->p[i][c]);
	  value_substract(Mat->p[j][c],m1,m2); 
	}
	for(c=0;c<k;++c) {
	  value_multiply(m1,piv,MatInv->p[j][c]);
	  value_multiply(m2,x,MatInv->p[i][c]);
	  value_substract(MatInv->p[j][c],m1,m2);
	}
	      
	/* Simplify row(j) of the two matrices 'Mat' and 'MatInv' by */
	/* dividing the rows with the common GCD.                     */
	Vector_Gcd(&MatInv->p[j][0],k,&m1);
	Vector_Gcd(&Mat->p[j][0],k,&m2);
	Gcd(m1,m2,&gcd);
	if(value_notone_p(gcd)) {
	  for(c=0;c<k;++c) {
	    value_division(Mat->p[j][c],Mat->p[j][c],gcd);
	    value_division(MatInv->p[j][c],MatInv->p[j][c],gcd);
	  }
	}
      }
    }
  }
  
  /* Simplify every row so that 'Mat' reduces to Identity matrix. Perform  */
  /* the same sequence of operations on the matrix 'MatInv'.               */
  for(j=0;j<k;++j) {
    value_assign(MatInv->p[j][k],Mat->p[j][j]);
    
    /* Make the last column (denominator of each entry) of every row greater */
    /* than zero.                                                            */
    Vector_Normalize_Positive(&MatInv->p[j][0],k+1,k); 
  }
  
  /* Clear all the 'mpz_t' variables */
  value_clear(x);  value_clear(gcd); value_clear(piv);
  value_clear(m1); value_clear(m2);

  return 1;
} /* Mat_Inverse */

/*
 * Given (m x n) integer matrix 'X' and n x (k+1) rational matrix 'P', compute
 * the rational m x (k+1) rational matrix  'S'. The last column in each row of
 * the rational matrices is used to store the common denominator of elements
 * in a row.                              
 */
void rat_prodmat(mpz_matrix *S,mpz_matrix *X,mpz_matrix *P) {
  
  int i,j,k;
  int last_column_index = P->NbColumns - 1;
  mpz_t lcm, old_lcm,gcd,last_column_entry,s1,s2;
  mpz_t m1,m2;
  
  /* Initialize all the 'mpz_t' variables */
  value_init(lcm); value_init(old_lcm); value_init(gcd);
  value_init(last_column_entry); value_init(s1); value_init(s2);
  value_init(m1); value_init(m2);

  /* Compute the LCM of last column entries (denominators) of rows */
  value_assign(lcm,P->p[0][last_column_index]);	
  for(k=1;k<P->NbRows;++k) {
    value_assign(old_lcm,lcm);
    value_assign(last_column_entry,P->p[k][last_column_index]);
    Gcd(lcm,last_column_entry,&gcd);
    value_division(m1,last_column_entry,gcd);
    value_multiply(lcm,lcm,m1);
  }
  
  /* S[i][j] = Sum(X[i][k] * P[k][j] where Sum is extended over k = 1..nbrows*/
  for(i=0;i<X->NbRows;++i)
    for(j=0;j<P->NbColumns-1;++j) {
      
      /* Initialize s1 to zero. */
      value_set_si(s1,0);
      for(k=0;k<P->NbRows;++k) {
	
	/* If the LCM of last column entries is one, simply add the products */
	if(value_one_p(lcm)) {
	  value_set_si(s2,0);
	  value_multiply(s2,X->p[i][k],P->p[k][j]);
          value_addto(s1,s1,s2);
	}  
	
	/* Numerator (num) and denominator (denom) of S[i][j] is given by :- */
	/* numerator  = Sum(X[i][k]*P[k][j]*lcm/P[k][last_column_index]) and */
	/* denominator= lcm where Sum is extended over k = 1..nbrows.        */
	else {
	  value_multiply(m1,X->p[i][k],P->p[k][j]);
	  value_division(m2,lcm,P->p[k][last_column_index]);
	  value_multiply(s2,m1,m2);
	  value_addto(s1,s1,s2);
	}
      }	
      value_assign(S->p[i][j],s1);
    }
  
  for(i=0;i<S->NbRows;++i) {
    value_assign(S->p[i][last_column_index],lcm);

    /* Normalize the rows so that last element >=0 */
    Vector_Normalize_Positive(&S->p[i][0],S->NbColumns,S->NbColumns-1);
  }
  
  /* Clear all the 'mpz_t' variables */
  value_clear(lcm); value_clear(old_lcm); value_clear(gcd);
  value_clear(last_column_entry); value_clear(s1); value_clear(s2);
  value_clear(m1); value_clear(m2);
 
  return;
} /* rat_prodmat */

/*
 * Given a matrix 'Mat' and vector 'p1', compute the matrix-vector product 
 * and store the result in vector 'p2'. 
 */
void mpz_matrix_Vector_Product(mpz_matrix *Mat,mpz_t *p1,mpz_t *p2) {

  int NbRows, NbColumns, i, j;
  mpz_t **cm, *q, *cp1, *cp2;
  mpz_t s;
  
  value_init(s);
  NbRows=Mat->NbRows;
  NbColumns=Mat->NbColumns;
  
  cm = Mat->p;
  cp2 = p2;
  for(i=0;i<NbRows;i++) {
    q = *cm++;
    cp1 = p1;
    value_multiply(*cp2,*q,*cp1);
    q++;
    cp1++;
    
    /* *cp2 = *q++ * *cp1++ */
    for(j=1;j<NbColumns;j++) {
      
      value_set_si(s,0);
      value_multiply(s,*q, *cp1);
      value_addto(*cp2,*cp2,s);
      q++;
      cp1++;
    }
    cp2++;
  }
  value_clear(s);
  return;
} /* mpz_matrix_Vector_Product */

/*
 * Given a vector 'p1' and a matrix 'Mat', compute the vector-matrix product 
 * and store the result in vector 'p2'
 */
void Vector_mpz_matrix_Product(mpz_t *p1,mpz_matrix *Mat,mpz_t *p2) {
  
  int NbRows, NbColumns, i, j;
  mpz_t **cm, *cp1, *cp2;
  mpz_t s;
  
  value_init(s);
  NbRows=Mat->NbRows;
  NbColumns=Mat->NbColumns;
  cp2 = p2;
  cm  = Mat->p;
  for (j=0;j<NbColumns;j++) {
    cp1 = p1;
    value_multiply(*cp2,*(*cm+j),*cp1);
    cp1++;
    
    /* *cp2= *(*cm+j) * *cp1++; */
    for (i=1;i<NbRows;i++) {
      
      value_set_si(s,0);
      value_multiply(s,*(*(cm+i)+j),*cp1);
      value_addto(*cp2,*cp2,s);
      cp1++;
    }
    cp2++;
  }
  value_clear(s);
  return;
} /* Vector_mpz_matrix_Product */

/* 
 * Given matrices 'Mat1' and 'Mat2', compute the matrix product and store in 
 * matrix 'Mat3' 
 */
void mpz_matrix_Product(mpz_matrix *Mat1,mpz_matrix *Mat2,mpz_matrix *Mat3) {
  
  int Size, i, j, k;
  unsigned NbRows, NbColumns;
  mpz_t **q1, **q2, *p1, *p3,sum,s;
  
  NbRows    = Mat1->NbRows;
  NbColumns = Mat2->NbColumns;
  
  Size      = Mat1->NbColumns;
  if(Mat2->NbRows!=Size||Mat3->NbRows!=NbRows||Mat3->NbColumns!=NbColumns) {
    fprintf(stderr, "? mpz_matrix_Product : incompatable matrix dimension\n");
    return;
  }     
  value_init(sum); value_init(s);
  p3 = Mat3->p_Init;
  q1 = Mat1->p;
  q2 = Mat2->p;
  
  /* Mat3[i][j] = Sum(Mat1[i][k]*Mat2[k][j] where sum is over k = 1..nbrows */
  for (i=0;i<NbRows;i++) {
    for (j=0;j<NbColumns;j++) {
      p1 = *(q1+i);
      value_set_si(sum,0);
      for (k=0;k<Size;k++) {
	
	/* sum+=*p1++ * *(*(q2+k)+j); */
	value_set_si(s,0);
	value_multiply(s,*p1, *(*(q2+k)+j));
	value_addto(sum,sum,s);
	p1++;
      }
      value_assign(*p3,sum);
      p3++;
    }
  }
  value_clear(sum); value_clear(s);
  return;
} /* mpz_matrix_Product */
  
/*
 * Given a rational matrix 'Mat'(k x k), compute its inverse rational matrix 
 * 'MatInv' k x k. The last column of each row in matrix MatInv is used 
 * to store the common denominator of the entries in a row. The output is 1,
 * if 'Mat' is non-singular (invertible), otherwise the output is 0. Note:: 
 * (1) mpz_matrix 'Mat' is modified during the inverse operation.
 * (2) mpz_matrix 'MatInv' must be preallocated before passing into this function.
 */
int mpz_matrix_Inverse(mpz_matrix *Mat,mpz_matrix *MatInv ) {
  
  int i, k, j, c;
  mpz_t x, gcd, piv;
  mpz_t m1,m2;
  mpz_t *den;
  
  if(Mat->NbRows != Mat->NbColumns) {
   fprintf(stderr,"Trying to invert a non-square matrix !\n");
    return 0;
  }
  
  /* Initialize all the 'mpz_t' variables */
  value_init(x);  value_init(gcd); value_init(piv);
  value_init(m1); value_init(m2);

  k = Mat->NbRows; 

  /* Initialise MatInv */
  Vector_Set(MatInv->p[0],0,k*k);

  /* Initialize 'MatInv' to Identity matrix form. Each diagonal entry is set*/
  /* to 1. Last column of each row (denominator of each entry in a row) is  */
  /* also set to 1.                                                         */ 
  for(i=0;i<k;++i) {
    value_set_si(MatInv->p[i][i],1);	
    // value_set_si(MatInv->p[i][k],1);	/* denum */
  }  
  /* Apply Gauss-Jordan elimination method on the two matrices 'Mat' and  */
  /* 'MatInv' in parallel.                                                */
  for(i=0;i<k;++i) {
    
    /* Check if the diagonal entry (new pivot) is non-zero or not */
    if(value_zero_p(Mat->p[i][i])) {   	
      
      /* Search for a non-zero pivot down the column(i) */
      for(j=i;j<k;++j)      
	if(value_notzero_p(Mat->p[j][i]))
	  break;
      
      /* If no non-zero pivot is found, the matrix 'Mat' is non-invertible */
      /* Return 0.                                                         */
      if(j==k) {
	
	/* Clear all the 'mpz_t' variables */
	value_clear(x);  value_clear(gcd); value_clear(piv);
	value_clear(m1); value_clear(m2);
	return 0;
      }	
      
      /* Exchange the rows, row(i) and row(j) so that the diagonal element */
      /* Mat->p[i][i] (pivot) is non-zero. Repeat the same operations on    */
      /* matrix 'MatInv'.                                                   */
      for(c=0;c<k;++c) {

	/* Interchange rows, row(i) and row(j) of matrix 'Mat'    */
	value_assign(x,Mat->p[j][c]);
	value_assign(Mat->p[j][c],Mat->p[i][c]);
	value_assign(Mat->p[i][c],x);
	
	/* Interchange rows, row(i) and row(j) of matrix 'MatInv' */
	value_assign(x,MatInv->p[j][c]);
	value_assign(MatInv->p[j][c],MatInv->p[i][c]);
	value_assign(MatInv->p[i][c],x);
      }
    }
    
    /* Make all the entries in column(i) of matrix 'Mat' zero except the */
    /* diagonal entry. Repeat the same sequence of operations on matrix  */
    /* 'MatInv'.                                                         */
    for(j=0;j<k;++j) {
      if (j==i) continue;	         /* Skip the pivot */
      value_assign(x,Mat->p[j][i]);
      if(value_notzero_p(x)) {
	value_assign(piv,Mat->p[i][i]);
	Gcd(x,piv,&gcd);
	if (value_notone_p(gcd) ) {
	  value_division(x,x,gcd);
	  value_division(piv,piv,gcd);
	}
	for(c=((j>i)?i:0);c<k;++c) {
	  value_multiply(m1,piv,Mat->p[j][c]);
	  value_multiply(m2,x,Mat->p[i][c]);
	  value_substract(Mat->p[j][c],m1,m2); 
	}
	for(c=0;c<k;++c) {
	  value_multiply(m1,piv,MatInv->p[j][c]);
	  value_multiply(m2,x,MatInv->p[i][c]);
	  value_substract(MatInv->p[j][c],m1,m2);
	}
	      
	/* Simplify row(j) of the two matrices 'Mat' and 'MatInv' by */
	/* dividing the rows with the common GCD.                     */
	Vector_Gcd(&MatInv->p[j][0],k,&m1);
	Vector_Gcd(&Mat->p[j][0],k,&m2);
	Gcd(m1,m2,&gcd);
	if(value_notone_p(gcd)) {
	  for(c=0;c<k;++c) {
	    value_division(Mat->p[j][c],Mat->p[j][c],gcd);
	    value_division(MatInv->p[j][c],MatInv->p[j][c],gcd);
	  }
	}
      }
    }
  }
  
  /* Find common denom for each row */ 
   den = (mpz_t *)malloc(k*sizeof(mpz_t));
   value_set_si(x,1);
   for(j=0 ; j<k ; ++j) {
     value_init(den[j]);
     value_assign(den[j],Mat->p[j][j]);
     
     /* gcd is always positive */
     Vector_Gcd(&MatInv->p[j][0],k,&gcd);
     Gcd(gcd,den[j],&gcd);
     if (value_neg_p(den[j])) 
       value_oppose(gcd,gcd); /* make denominator positive */
     if (value_notone_p(gcd)) {
       for (c=0; c<k; c++) 
	 value_division(MatInv->p[j][c],MatInv->p[j][c],gcd); /* normalize */
       value_division(den[j],den[j],gcd);
     }  
     Gcd(x,den[j],&gcd);
     value_division(m1,den[j],gcd);
     value_multiply(x,x,m1);
   }
   if (value_notone_p(x)) 
     for(j=0 ; j<k ; ++j) {       
       for (c=0; c<k; c++) {
	 value_division(m1,x,den[j]);
	 value_multiply(MatInv->p[j][c],MatInv->p[j][c],m1);  /* normalize */
       }
     }

   /* Clear all the 'mpz_t' variables */
   for(j=0 ; j<k ; ++j) {
     value_clear(den[j]);
   }  
   value_clear(x);  value_clear(gcd); value_clear(piv);
   value_clear(m1); value_clear(m2);
   free(den);
   
   return 1;
} /* mpz_matrix_Inverse */

#endif /* #if 0 */






