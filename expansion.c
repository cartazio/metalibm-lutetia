/*
 * Copyright 2008-2014 by 
 * 
 * Laboratoire de l'Informatique du Parallélisme, 
 * UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668
 *
 * and by
 *
 * Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
 * UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France.
 *
 * Sollya metalibm test driver helper file.
 *
 * Contributor: Christoph Quirin Lauter (ENS Lyon/ UPMC LIP6 PEQUAN) 
 *              christoph.lauter@ens-lyon.org
 *
 * This file is part of the metalibm library developed by the Arenaire
 * project at Ecole Normale Superieure de Lyon and by the equipe PEQUAN
 * at Universite Pierre et Marie Curie Paris 6.
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

#include <gmp.h>
#include <mpfr.h>
#include <fpu_control.h>
#include "expansion.h"
#include <sys/time.h>
#include <time.h>
#include <stdint.h>

/* A constant always being false */
#define __ALWAYS_FALSE 0


/* Some declarations */

double sollya_mpfr_get_d(mpfr_t, mp_rnd_t);

#if defined(D_TO_D)
extern void POLYNOMIALNAME(double *, double);
#elif defined(D_TO_DD)
extern void POLYNOMIALNAME(double *, double *, double);
#elif defined(D_TO_TD)
extern void POLYNOMIALNAME(double *, double *, double *, double);
#elif defined(DD_TO_DD)
extern void POLYNOMIALNAME(double *, double *, double, double);
#elif defined(DD_TO_TD)
extern void POLYNOMIALNAME(double *, double *, double *, double, double);
#elif defined (TD_TO_TD)
extern void POLYNOMIALNAME(double *, double *, double *, double, double, double);
#endif

void mpfr_to_double(double *dh, mpfr_t op) {
  *dh = sollya_mpfr_get_d(op, GMP_RNDN);
  mpfr_set_d(op, *dh, GMP_RNDN);
}

void mpfr_to_doubledouble(double *dh, double *dm, mpfr_t op) {
  mpfr_t rest, temp;

  mpfr_init2(rest,mpfr_get_prec(op));
  mpfr_init2(temp,mpfr_get_prec(op));

  mpfr_set(rest,op,GMP_RNDN);

  *dh = sollya_mpfr_get_d(rest, GMP_RNDN);
  mpfr_set_d(temp,*dh,GMP_RNDN);
  mpfr_set(op,temp,GMP_RNDN);
  mpfr_sub(rest,rest,temp,GMP_RNDN);

  *dm = sollya_mpfr_get_d(rest, GMP_RNDN);
  mpfr_set_d(temp,*dm,GMP_RNDN);
  mpfr_add(op,op,temp,GMP_RNDN);
  
  mpfr_clear(rest);
  mpfr_clear(temp);
}

void mpfr_to_tripledouble(double *dh, double *dm, double *dl, mpfr_t op) {
  mpfr_t rest, temp;

  mpfr_init2(rest,mpfr_get_prec(op));
  mpfr_init2(temp,mpfr_get_prec(op));

  mpfr_set(rest,op, GMP_RNDN);

  *dh = sollya_mpfr_get_d(rest, GMP_RNDN);
  mpfr_set_d(temp,*dh,GMP_RNDN);
  mpfr_set(op,temp,GMP_RNDN);
  mpfr_sub(rest,rest,temp,GMP_RNDN);

  *dm = sollya_mpfr_get_d(rest, GMP_RNDN);
  mpfr_set_d(temp,*dm,GMP_RNDN);
  mpfr_add(op,op,temp,GMP_RNDN);
  mpfr_sub(rest,rest,temp,GMP_RNDN);

  *dl = sollya_mpfr_get_d(rest, GMP_RNDN);
  mpfr_set_d(temp,*dl,GMP_RNDN);
  mpfr_add(op,op,temp,GMP_RNDN);

  mpfr_clear(rest);
  mpfr_clear(temp);
}

void double_to_mpfr(mpfr_t rop, double dh) {
  mpfr_set_d(rop,dh,GMP_RNDN);
}

void doubledouble_to_mpfr(mpfr_t rop, double dh, double dm) {
  mpfr_t temp;
  mp_prec_t prec;

  prec = mpfr_get_prec(rop);

  if (prec < 53) prec = 53;

  mpfr_init2(temp, prec);

  mpfr_set_d(rop, dh, GMP_RNDN);

  mpfr_set_d(temp, dm, GMP_RNDN);
  mpfr_add(rop, rop, temp, GMP_RNDN);

  mpfr_clear(temp);
}

void tripledouble_to_mpfr(mpfr_t rop, double dh, double dm, double dl) {
  mpfr_t temp;
  mp_prec_t prec;

  prec = mpfr_get_prec(rop);

  if (prec < 53) prec = 53;

  mpfr_init2(temp, prec);

  mpfr_set_d(rop, dh, GMP_RNDN);

  mpfr_set_d(temp, dm, GMP_RNDN);
  mpfr_add(rop, rop, temp, GMP_RNDN);

  mpfr_set_d(temp, dl, GMP_RNDN);
  mpfr_add(rop, rop, temp, GMP_RNDN);


  mpfr_clear(temp);
}


void f(mpfr_t y, mpfr_t xMpfr) {
  unsigned short oldcw, cw;
#if defined(D_TO_D)
  double x;
  double resh;
#elif defined(D_TO_DD)
  double x;
  double resh, resm;
#elif defined(D_TO_TD)
  double x;
  double resh, resm, resl;
#elif defined(DD_TO_DD)
  double xh, xm;
  double resh, resm;
#elif defined(DD_TO_TD)
  double xh, xm;
  double resh, resm, resl;
#elif defined (TD_TO_TD)
  double xh, xm, xl;
  double resh, resm, resl;
#endif

  _FPU_GETCW(oldcw);
  cw = (_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE;
  _FPU_SETCW(cw);

#if defined(D_TO_D)
  mpfr_to_double(&x, xMpfr);
  POLYNOMIALNAME(&resh, x);
  double_to_mpfr(y, resh);
#elif defined(D_TO_DD)
  mpfr_to_double(&x, xMpfr);
  POLYNOMIALNAME(&resh, &resm, x);
  doubledouble_to_mpfr(y, resh, resm);
#elif defined(D_TO_TD)
  mpfr_to_double(&x, xMpfr);
  POLYNOMIALNAME(&resh, &resm, &resl, x);
  tripledouble_to_mpfr(y, resh, resm, resl);
#elif defined(DD_TO_DD)
  mpfr_to_doubledouble(&xh, &xm, xMpfr);
  POLYNOMIALNAME(&resh, &resm, xh, xm);
  doubledouble_to_mpfr(y, resh, resm);
#elif defined(DD_TO_TD)
  mpfr_to_doubledouble(&xh, &xm, xMpfr);
  POLYNOMIALNAME(&resh, &resm, &resl, xh, xm);
  tripledouble_to_mpfr(y, resh, resm, resl);
#elif defined(TD_TO_TD)
  mpfr_to_tripledouble(&xh, &xm, &xl, xMpfr);
  POLYNOMIALNAME(&resh, &resm, &resl, xh, xm, xl);
  tripledouble_to_mpfr(y, resh, resm, resl);
#else
#warning You must define one of the macros for the argument and result formats
  mpfr_set(y,xMpfr,GMP_RNDN);
#endif

  _FPU_SETCW(oldcw);

}

#ifdef ASSEMBLY_TIMING

#warning Assembly timing sequence used

#define READ_TIME_COUNTER(time)                        \
  __asm__ __volatile__(                                \
          "xorl %%eax,%%eax\n\t"                       \
          "cpuid\n\t"                                  \
          "rdtsc\n\t"                                  \
          "movl %%eax,(%0)\n\t"                        \
          "movl %%edx,4(%0)\n\t"                       \
          "xorl %%eax,%%eax\n\t"                       \
          "cpuid\n\t"                                  \
          : /* nothing */                              \
          : "S"((time))                                \
          : "eax", "ebx", "ecx", "edx", "memory")

#else

#define READ_TIME_COUNTER(time) do {                   \
  struct timeval t;                                    \
  gettimeofday(&t,NULL);                               \
  *time = (((uint64_t) t.tv_sec) *                     \
	   ((uint64_t) 1000000)) +                     \
    ((uint64_t) t.tv_usec);                            \
} while ( __ALWAYS_FALSE );

#endif

int timefunc(int *timing, void **args) {
  unsigned short oldcw, cw;
#if defined(D_TO_D)
  double x;
  double resh;
#elif defined(D_TO_DD)
  double x;
  double resh, resm;
#elif defined(D_TO_TD)
  double x;
  double resh, resm, resl;
#elif defined(DD_TO_DD)
  double xh, xm;
  double resh, resm;
#elif defined(DD_TO_TD)
  double xh, xm;
  double resh, resm, resl;
#elif defined (TD_TO_TD)
  double xh, xm, xl;
  double resh, resm, resl;
#endif
  mpfr_t xMpfr;
  mpfr_t a, b;
  mpfr_t h;
  int steps, iterations, actualSteps;
  int i;
  uint64_t start, end;
  uint64_t usecs;
  double overalltime;

  mpfr_init2(xMpfr, 161);
  mpfr_init2(a, mpfr_get_prec(*((mpfr_t *) (args[0]))));
  mpfr_init2(b, mpfr_get_prec(*((mpfr_t *) (args[1]))));
  mpfr_set(a, *((mpfr_t *) (args[0])), GMP_RNDN);
  mpfr_set(b, *((mpfr_t *) (args[1])), GMP_RNDN);
  mpfr_init2(h, 161);

  mpfr_set(xMpfr, a, GMP_RNDU);

  steps = *((int *) (args[2]));
  iterations = *((int *) (args[3]));

  mpfr_sub(h, b, a, GMP_RNDU);
  mpfr_div_si(h, h, steps, GMP_RNDU);

  overalltime = 0;
  actualSteps = 0;

  while (mpfr_cmp(xMpfr,b) <= 0) {

    _FPU_GETCW(oldcw);
    cw = (_FPU_DEFAULT & ~_FPU_EXTENDED)|_FPU_DOUBLE;
    _FPU_SETCW(cw);

#if defined(D_TO_D)
    mpfr_to_double(&x, xMpfr);
    READ_TIME_COUNTER(&start);
    for (i=0;i<iterations;i++) POLYNOMIALNAME(&resh, x);
    READ_TIME_COUNTER(&end);
#elif defined(D_TO_DD)
    mpfr_to_double(&x, xMpfr);
    READ_TIME_COUNTER(&start);
    for (i=0;i<iterations;i++) POLYNOMIALNAME(&resh, &resm, x);
    READ_TIME_COUNTER(&end);
#elif defined(D_TO_TD)
    mpfr_to_double(&x, xMpfr);
    READ_TIME_COUNTER(&start);
    for (i=0;i<iterations;i++) POLYNOMIALNAME(&resh, &resm, &resl, x);
    READ_TIME_COUNTER(&end);
#elif defined(DD_TO_DD)
    mpfr_to_doubledouble(&xh, &xm, xMpfr);
    READ_TIME_COUNTER(&start);
    for (i=0;i<iterations;i++) POLYNOMIALNAME(&resh, &resm, xh, xm);
    READ_TIME_COUNTER(&end);
#elif defined(DD_TO_TD)
    mpfr_to_doubledouble(&xh, &xm, xMpfr);
    READ_TIME_COUNTER(&start);
    for (i=0;i<iterations;i++) POLYNOMIALNAME(&resh, &resm, &resl, xh, xm);
    READ_TIME_COUNTER(&end);
#elif defined(TD_TO_TD)
    mpfr_to_tripledouble(&xh, &xm, &xl, xMpfr);
    READ_TIME_COUNTER(&start);
    for (i=0;i<iterations;i++) POLYNOMIALNAME(&resh, &resm, &resl, xh, xm, xl);
    READ_TIME_COUNTER(&end);
#else
#warning You must define one of the macros for the argument and result formats
#endif

    _FPU_SETCW(oldcw);

    mpfr_add(xMpfr, xMpfr, h, GMP_RNDU);

    if (end >= start) {
      usecs = end - start; 
      actualSteps++;
    } else {
      usecs = 0;
    }

    overalltime += usecs;
    
  }

  if (actualSteps != 0) overalltime /= actualSteps; 

  mpfr_clear(b);
  mpfr_clear(a);
  mpfr_clear(xMpfr);
  mpfr_clear(h);

  *timing = overalltime;

  return 1;
}
