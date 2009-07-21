/*
 * Copyright 2008 by 
 * 
 * Laboratoire de l'Informatique du Parallélisme, 
 * UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668
 *
 * Sollya wrapper scripts for vandercoeff techniques by Torres and Brisebarre
 *
 * Contributor: Christoph Quirin Lauter (ENS Lyon) -- christoph.lauter@ens-lyon.fr
 *
 * This file is part of the metalibm library developed by the Arenaire
 * project at Ecole Normale Superieure de Lyon
 * 
 * This file is heavily inspired by the externalproc.c example file in the Sollya 
 * software distribution (http://sollya.gforge.inria.fr) under CeCill-C licence. 
 * It relies on the sollya.h header file of this project.
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

#include <mpfi.h>
#include <mpfr.h>
#include <stdlib.h>
#include "sollya.h"
#include "libVanderCoeffsSparse.h"

/* Example of an external procedure linked to an identifier in sollya

   Procedure vandercoeff will be linked by

   externalproc(vandercoeff, "./vandercoeff", (integer, list of integer, list of constant, list of constant, constant) -> list of range);

   Its signature is 
   
   vandercoeff(absRel, monomials, points, evals, maxerr) -> ranges of [lowerBound, upperBound]

*/


extern mp_prec_t tools_precision;

int noDuplicate(int count, int *array) {
  int i,k;

  for (i=0;i<count;i++) {
    for (k=0;k<count;k++) {
      if ((i != k) && (array[i] == array[k])) return 0;
    }
  }

  return 1;
}

int lengthChain(chain *c) {
  int i;
  chain *curr;

  i = 0;
  curr = c;
  while (curr != NULL) {
    i++;
    curr = curr->next;
  }

  return i;
}

chain *addElement(chain *c, void *elem) {
  chain *newChain;

  newChain = (chain *) malloc(sizeof(chain));
  newChain->next = c;
  newChain->value = elem;
  return newChain;
}

mp_prec_t getToolPrecision() {
  return tools_precision;
}

int myvandercoeff(chain **results, int absRel, chain *monomials, chain* points, chain *evals, mpfr_t maxerr) {
  int monomialsCount, pointsCount, evalsCount;
  mpfr_t *pointsArray;
  mpfr_t *evalsArray;
  int *monomialsArray;
  int res;
  int okay;
  chain *curr;
  int i;
  mpfr_t *upperBounds;
  mpfr_t *lowerBounds; 
  mpfi_t *tempMpfi;
  mp_prec_t prec, p;


  *results = NULL;
  mpfr_set_default_prec(getToolPrecision() * 2); 

  pointsCount = lengthChain(points);
  evalsCount = lengthChain(evals);
  
  if (pointsCount != evalsCount) {
     printf("Error in external procedure vandercoeff: the number of points does not correspond.\n"); 
     return 0; 
  }

  pointsArray = (mpfr_t *) calloc(pointsCount, sizeof(mpfr_t));
  evalsArray = (mpfr_t *) calloc(evalsCount, sizeof(mpfr_t));

  i = 0;
  curr = points;
  while (curr != NULL) {
    mpfr_init2(pointsArray[i],mpfr_get_prec(*((mpfr_t *) (curr->value))));
    mpfr_set(pointsArray[i],*((mpfr_t *) (curr->value)),GMP_RNDN);
    i++;
    curr = curr->next;
  }
  i = 0;
  curr = evals;
  while (curr != NULL) {
    mpfr_init2(evalsArray[i],mpfr_get_prec(*((mpfr_t *) (curr->value))));
    mpfr_set(evalsArray[i],*((mpfr_t *) (curr->value)),GMP_RNDN);
    i++;
    curr = curr->next;
  }
  
  monomialsCount = lengthChain(monomials);

  monomialsArray = (int *) calloc(monomialsCount, sizeof(int));
  i = 0;
  curr = monomials;
  while (curr != NULL) {
    monomialsArray[i] = *((int *) (curr->value));
    i++;
    curr = curr->next;
  }

  if (!noDuplicate(monomialsCount,monomialsArray)) {

    printf("Error in external procedure vandercoeff: duplicate monomial degrees in the list.\n"); 

    free(monomialsArray);

    for (i=0;i<pointsCount;i++) {
      mpfr_clear(pointsArray[i]);
    }
    for (i=0;i<evalsCount;i++) {
      mpfr_clear(evalsArray[i]);
    }

    free(pointsArray);
    free(evalsArray);

    return 0;
  }

  upperBounds = (mpfr_t *) calloc(monomialsCount, sizeof(mpfr_t));
  lowerBounds = (mpfr_t *) calloc(monomialsCount, sizeof(mpfr_t));

  if (absRel) 
    res = vanderCoeffsSparseRel((unsigned int) pointsCount, pointsArray, evalsArray, monomialsCount, (unsigned int *) monomialsArray, maxerr, lowerBounds, upperBounds);
  else 
    res = vanderCoeffsSparseAbs((unsigned int) pointsCount, pointsArray, evalsArray, monomialsCount, (unsigned int *) monomialsArray, maxerr, lowerBounds, upperBounds);

  okay = 0;

  if (res == 0) {
    
    okay = 1;

    for (i=0;i<monomialsCount;i++) {
      if (mpfr_cmp(lowerBounds[i], upperBounds[i]) > 0) {
	printf("Warning: some bounds determined by vanderCoeffsSparse are returned inversed.\n");
	okay = 0;
      }
    }

    if (okay) {
      
      curr = NULL; 
      for (i=0;i<monomialsCount;i++) {
	tempMpfi = (mpfi_t *) malloc(sizeof(mpfi_t));
	prec = mpfr_get_prec(upperBounds[i]);
	p = mpfr_get_prec(lowerBounds[i]);
	if (p > prec) prec = p;
	mpfi_init2(*tempMpfi,prec);
	mpfi_interv_fr(*tempMpfi,lowerBounds[i],upperBounds[i]);
	curr = addElement(curr, tempMpfi);
      }
      *results = curr;

    }

    for (i=0;i<monomialsCount;i++) {
      mpfr_clear(lowerBounds[i]);
      mpfr_clear(upperBounds[i]);
    }

  } else {

    printf("Warning: some error (%d) occured in vanderCoeffsSparse.\n",res);
    
  }

  free(monomialsArray);
  for (i=0;i<pointsCount;i++) {
    mpfr_clear(pointsArray[i]);
  }
  for (i=0;i<evalsCount;i++) {
    mpfr_clear(evalsArray[i]);
  }
  
  free(pointsArray);
  free(evalsArray);

  free(upperBounds);
  free(lowerBounds);
    
  return okay;
}


int vandercoeff(chain **results, void **args) {

  return myvandercoeff(results, *((int *) (args[0])), (chain *) (args[1]), (chain *) (args[2]), (chain *) (args[3]), *((mpfr_t *) (args[4])));
}


