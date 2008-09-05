/*
 * Copyright 2008 by 
 * 
 * Laboratoire de l'Informatique du Parallélisme, 
 * UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668
 *
 * C++ Interface for MPFR (specially patched version)
 *
 * Contributors:
 *  Nathalie Revol (INRIA) -- nathalie.revol@ens-lyon.fr
 *  Serge Torres (ENS Lyon) -- serge.torres@ens-lyon.fr
 *  Christoph Lauter (ENS Lyon) -- christoph.lauter@ens-lyon.fr
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

// C++ Interface for MPFR
// NR 21-11-00 (first steps in C++...)
// NR 23-02-04 (next steps in C++...)
// ST 2006-07-06 (Trying to catch up with MPFR 2.2.0)
// ST 2007-01-03 (Trying to catch up with gcc 4.x.x)

#ifndef _MPFR_CLASS_H_
#define _MPFR_CLASS_H_

#include <stdio.h>

// To be able to use the C++ string
using namespace std;
#include <string>
#define String string // Hmm, what for?

#include <iostream>
#include <math.h>


#include <gmp.h>

extern void * (*__gmp_allocate_func)(size_t);
extern void * (*__gmp_reallocate_func)(void *, size_t, size_t);
extern void   (*__gmp_free_func)(void *, size_t);


#include <mpfr.h>

#ifndef _GIVARO_REF_COUNTER_H_
#define _GIVARO_REF_COUNTER_H_
// ==================================================================== //
// Definition of the Counter class, Counter
// (c) copyright GIVARO 1994                                            
// author: Th. Gautier                                                
// version : 2.7
// date: 1995
// This class definition objects to handle reference
// counter for memory allocation (eg array0). 
// ====================================================================
#include <stddef.h>

class RefCounter {
public:
   // Cstor and Dstor
inline RefCounter( long l = 0) : counter(l) {} 
//inline RefCounter( const RefCounter& ) : counter(C.counter) {} 
inline ~RefCounter() {}

  //  Return the value
inline long  getvalue() const { return counter ; } 
inline long  val() const { return counter ; }
  // Return a ref to the counter
inline long& refvalue() { return counter ; }
  // Increments the counter and returns the new value 
inline long  incr() { return ++counter ; }
  // Decrements the value and returns the new value 
inline long  decr() { return --counter ; }

protected:
  long counter ;
} ;

#endif
// ==================================================================== //
// End of definition of the Counter class, Counter
// ==================================================================== //


#define UNAFFECTED_INEXACT_FLAG 500
#define EXACT_FLAG 0
#define INEXACT_FLAG 25
class InexactFlag {
public:
   // Cstor and Dstor
  inline InexactFlag( int l = UNAFFECTED_INEXACT_FLAG) : inexact_flag(l) {} 
  inline ~InexactFlag() {}

  //  Return the value
  inline int  getvalue() const { return inexact_flag ; } 
  inline int  val() const { return inexact_flag ; }
  // Return a ref to the inexact_flag
  inline int& refvalue() { return inexact_flag ; }

protected:
  int inexact_flag ;
} ;



// Rounding modes
typedef mp_rnd_t RoundingMode;
static const RoundingMode RoundUp=mp_rnd_t(GMP_RNDU); 
static const RoundingMode RoundDown=mp_rnd_t(GMP_RNDD); 
static const RoundingMode RoundNearest=mp_rnd_t(GMP_RNDN); 
static const RoundingMode RoundToZero=mp_rnd_t(GMP_RNDZ);



class MpfrClass {

  protected:
    mpfr_t mpfr_rep;	// representation of the real in a mpfr format
    // Precision
    typedef mp_prec_t PrecisionType;
    InexactFlag *inexact;
    RefCounter *nbref;
    // Current default precision
    static PrecisionType CurrPrecision;
    // Current default rounding mode
    static RoundingMode CurrRndMode;

  public:
    // constructors and destructors
    MpfrClass ();
    MpfrClass (const mpfr_t& mpfr,
                RoundingMode rnd = CurrRndMode,
                PrecisionType prec = CurrPrecision);
    MpfrClass (double d,
              RoundingMode rnd = CurrRndMode,
              PrecisionType prec = CurrPrecision);
    MpfrClass (long double d,
              RoundingMode rnd = CurrRndMode,
              PrecisionType prec = CurrPrecision);
    MpfrClass (int i, 
              RoundingMode rnd = CurrRndMode, 
              PrecisionType prec = CurrPrecision);
    MpfrClass (unsigned int i, 
              RoundingMode rnd = CurrRndMode, 
              PrecisionType prec = CurrPrecision);
    MpfrClass (long int i, 
              RoundingMode rnd = CurrRndMode, 
              PrecisionType prec = CurrPrecision);
   MpfrClass (unsigned long int i, 
              RoundingMode rnd = CurrRndMode, 
              PrecisionType prec = CurrPrecision);
    MpfrClass(string s);
    MpfrClass (mpz_srcptr z,
              RoundingMode rnd = CurrRndMode,
              PrecisionType prec = CurrPrecision);
    MpfrClass (mpq_srcptr q, 
              RoundingMode rnd = CurrRndMode,
              PrecisionType prec = CurrPrecision);
    MpfrClass (mpfr_t r,
              RoundingMode rnd = CurrRndMode,
              PrecisionType prec = CurrPrecision);
    MpfrClass (const MpfrClass& r);
    ~MpfrClass ();
    
    // Accessors
    mpfr_t& getValue();
 
    // Assignment and copy operators
    MpfrClass& operator = (const MpfrClass& r);
    MpfrClass& copy (const MpfrClass& r, 
                    RoundingMode rnd = CurrRndMode, 
                    PrecisionType prec = CurrPrecision);

    // All names of the isxxx functions have been changed
    // because with some c++ implementation (gcc 3.4.x and gcc 4.x.x)
    // the "isinf" and "isnan" functions of cmath are implemented as macros. 
    // In this case we got compilation errors with the original code.
    // To avoid this and to keep the ability to use mpfr++
    // and standard c++ functions in the same code we have decided to 
    // change the names in mpfr++ (and not just undefine the macros). 
    // We changed "iszero" and "isnumber" as well for name consitency among 
    // the isxxx functions family.
    // We has to propagate the changes into the implementation, notably
    // switching many invocation of "iszero" to "isZero".
    
    // Is equal to zero?
    friend bool isZero(const MpfrClass& r) { return ( (mpfr_nan_p(r.mpfr_rep) != 0) && (mpfr_zero_p(r.mpfr_rep) == 0) ); }
    friend bool isInf(const MpfrClass& r) { return (mpfr_inf_p(r.mpfr_rep) != 0); }
    friend bool isNan(const MpfrClass& r) { return (mpfr_nan_p(r.mpfr_rep) != 0); }
    friend bool isNumber(const MpfrClass& r) { return (mpfr_number_p(r.mpfr_rep) != 0); }
    friend int sign(const MpfrClass& r) { return (MPFR_SIGN(r.mpfr_rep)); }

    // Precision and rounding mode
    static void SetDefaultPrecision (PrecisionType newprec);
    void SetPrecision (PrecisionType newprec);
    const static PrecisionType GetDefaultPrecision ();
    PrecisionType GetPrecision () const;
    static void SetDefaultRndMode (RoundingMode newrndmode);
    static const RoundingMode GetDefaultRndMode ();

    // Changing the precision: should be in place or not? In place
    // => to round not in place, copy and round
    void ChangePrec (RoundingMode rnd = CurrRndMode,
                     PrecisionType prec = CurrPrecision);

    // Rounding in the direction of rnd2 when the result is computed 
    // in the direction of rnd1, with a number of wrong bits given as argument.
    // Should be in place or not? In place
    int CanRound (mp_prec_t nb_wrong_bits, RoundingMode rnd1 = CurrRndMode,
                  RoundingMode rnd2 = CurrRndMode,
                  PrecisionType prec = CurrPrecision);

    // "Constants" (but they depend on the precision and the rounding mode)
    static MpfrClass Pi(RoundingMode rnd = CurrRndMode, 
                        PrecisionType prec = CurrPrecision) ;
    static MpfrClass Log2(RoundingMode rnd = CurrRndMode, 
                        PrecisionType prec = CurrPrecision) ;
    static MpfrClass Euler(RoundingMode rnd = CurrRndMode, 
                        PrecisionType prec = CurrPrecision) ;

    // Comparison operators
    friend int compare (const MpfrClass& r1, const MpfrClass& r2);
    friend int compare (const MpfrClass& r1, const double r2);
    friend int compare (const MpfrClass& r1, const int r2);
    friend int compare (const MpfrClass& r1, const unsigned int r2);
    friend int compare (const MpfrClass& r1, const long int r2);
    friend int compare (const MpfrClass& r1, const unsigned long int r2);


    // Arithmetic operators
    // Philosophy: are members only the operations between MpfrClass 
    //static MpfrClass& operator+ (const MpfrClass& r1, const MpfrClass& r2) ;
    MpfrClass operator+ (const MpfrClass& r) const;
    friend MpfrClass operator+ (const MpfrClass& r1, const double r2) ;
    friend MpfrClass operator+ (const MpfrClass& r1, const int r2) ;
    friend MpfrClass operator+ (const MpfrClass& r1, const unsigned int r2) ;
    friend MpfrClass operator+ (const MpfrClass& r1, const long int r2) ;
    friend MpfrClass operator+ (const MpfrClass& r1, const unsigned long int r2) ;
    friend MpfrClass operator+ (const MpfrClass& r1, const mpz_srcptr r2) ;
    friend MpfrClass operator+ (const MpfrClass& r1, const mpq_srcptr r2) ;
    friend MpfrClass operator+ (const double r1, const MpfrClass& r2) ;
    friend MpfrClass operator+ (const int r1, const MpfrClass& r2) ;
    friend MpfrClass operator+ (const unsigned int r1, const MpfrClass& r2) ;
    friend MpfrClass operator+ (const long int r1, const MpfrClass& r2) ;
    friend MpfrClass operator+ (const unsigned long int r1, const MpfrClass& r2) ;
    friend MpfrClass operator+ (const mpz_srcptr r1, const MpfrClass& r2) ;
    friend MpfrClass operator+ (const mpq_srcptr r1, const MpfrClass& r2) ;
    MpfrClass& operator+= (const MpfrClass& r) ;
    MpfrClass& operator+= (const double r) ;
    MpfrClass& operator+= (const int r) ;
    MpfrClass& operator+= (const unsigned int r) ;
    MpfrClass& operator+= (const long int r) ;
    MpfrClass& operator+= (const unsigned long int r) ;
    MpfrClass& operator+= (const mpz_srcptr r) ;
    MpfrClass& operator+= (const mpq_srcptr r) ;
    static void add (MpfrClass& res, 
                          const MpfrClass& r1, const MpfrClass& r2,
                          RoundingMode rnd = CurrRndMode);

    //static MpfrClass& operator- (const MpfrClass& r1, const MpfrClass& r2);
    MpfrClass operator- (const MpfrClass& r) const;
    friend MpfrClass operator- (const MpfrClass& r1, const double r2) ;
    friend MpfrClass operator- (const MpfrClass& r1, const int r2);
    friend MpfrClass operator- (const MpfrClass& r1, const unsigned int r2);
    friend MpfrClass operator- (const MpfrClass& r1, const long int r2);
    friend MpfrClass operator- (const MpfrClass& r1, const unsigned long int r2);
    friend MpfrClass operator- (const MpfrClass& r1, const mpz_srcptr r2);
    friend MpfrClass operator- (const MpfrClass& r1, const mpq_srcptr r2);
    friend MpfrClass operator- (const double r1, const MpfrClass& r2);
    friend MpfrClass operator- (const int r1, const MpfrClass& r2);
    friend MpfrClass operator- (const unsigned int r1, const MpfrClass& r2);
    friend MpfrClass operator- (const long int r1, const MpfrClass& r2);
    friend MpfrClass operator- (const unsigned long int r1, const MpfrClass& r2);
    friend MpfrClass operator- (const mpz_srcptr r1, const MpfrClass& r2);
    friend MpfrClass operator- (const mpq_srcptr r1, const MpfrClass& r2);
    MpfrClass operator- () const;
    MpfrClass& operator-= (const MpfrClass& r) ;
    MpfrClass& operator-= (const double r) ;
    MpfrClass& operator-= (const int r) ;
    MpfrClass& operator-= (const long int r) ;
    MpfrClass& operator-= (const unsigned int r) ;
    MpfrClass& operator-= (const unsigned long int r) ;
    MpfrClass& operator-= (const mpz_srcptr r) ;
    MpfrClass& operator-= (const mpq_srcptr r) ;
    static void neg (MpfrClass& res, 
                          const MpfrClass& r, 
                          RoundingMode rnd = CurrRndMode);
    static void sub (MpfrClass& res, 
                          const MpfrClass& r1, const MpfrClass& r2,
                          RoundingMode rnd = CurrRndMode);

    //static MpfrClass& operator* (const MpfrClass& r1, const MpfrClass& r2) ;
    MpfrClass operator* (const MpfrClass& r) const;
    friend MpfrClass operator* (const MpfrClass& r1, const double r2) ;
    friend MpfrClass operator* (const MpfrClass& r1, const int r2) ;
    friend MpfrClass operator* (const MpfrClass& r1, const unsigned int r2) ;
    friend MpfrClass operator* (const MpfrClass& r1, const long int r2) ;
    friend MpfrClass operator* (const MpfrClass& r1, const unsigned long int r2) ;
    friend MpfrClass operator* (const MpfrClass& r1, const mpz_srcptr r2) ;
    friend MpfrClass operator* (const MpfrClass& r1, const mpq_srcptr r2) ;
    friend MpfrClass operator* (const double r1, const MpfrClass& r2);
    friend MpfrClass operator* (const int r1, const MpfrClass& r2);
    friend MpfrClass operator* (const unsigned int r1, const MpfrClass& r2);
    friend MpfrClass operator* (const long int r1, const MpfrClass& r2);
    friend MpfrClass operator* (const unsigned long int r1, const MpfrClass& r2);
    friend MpfrClass operator* (const mpz_srcptr r1, const MpfrClass& r2);
    friend MpfrClass operator* (const mpq_srcptr r1, const MpfrClass& r2);
    MpfrClass& operator*= (const MpfrClass& r) ;
    MpfrClass& operator*= (const double r) ;
    MpfrClass& operator*= (const int r) ;
    MpfrClass& operator*= (const unsigned int r) ;
    MpfrClass& operator*= (const long int r) ;
    MpfrClass& operator*= (const unsigned long int r) ;
    MpfrClass& operator*= (const mpz_srcptr r) ;
    MpfrClass& operator*= (const mpq_srcptr r) ;
    static void mul (MpfrClass& res, 
                          const MpfrClass& r1, const MpfrClass& r2,
                          RoundingMode rnd = CurrRndMode);

    //static MpfrClass& operator/ (const MpfrClass& r1, const MpfrClass& r2) const;
    MpfrClass operator/ (const MpfrClass& r) const;
    friend MpfrClass operator/ (const MpfrClass& r1, const double r2) ;
    friend MpfrClass operator/ (const MpfrClass& r1, const int r2) ;
    friend MpfrClass operator/ (const MpfrClass& r1, const unsigned int r2) ;
    friend MpfrClass operator/ (const MpfrClass& r1, const long int r2) ;
    friend MpfrClass operator/ (const MpfrClass& r1, const unsigned long int r2) ;
    friend MpfrClass operator/ (const MpfrClass& r1, const mpz_srcptr r2) ;
    friend MpfrClass operator/ (const MpfrClass& r1, const mpq_srcptr r2) ;
    friend MpfrClass operator/ (const double r1, const MpfrClass& r2) ;
    friend MpfrClass operator/ (const int r1, const MpfrClass& r2) ;
    friend MpfrClass operator/ (const unsigned int r1, const MpfrClass& r2) ;
    friend MpfrClass operator/ (const long int r1, const MpfrClass& r2) ;
    friend MpfrClass operator/ (const unsigned long int r1, const MpfrClass& r2) ;
    friend MpfrClass operator/ (const mpz_srcptr r1, const MpfrClass& r2) ;
    friend MpfrClass operator/ (const mpq_srcptr r1, const MpfrClass& r2) ;
    MpfrClass& operator/= (const MpfrClass& r) ;
    MpfrClass& operator/= (const double r) ;
    MpfrClass& operator/= (const int r) ;
    MpfrClass& operator/= (const unsigned int r) ;
    MpfrClass& operator/= (const long int r) ;
    MpfrClass& operator/= (const unsigned long int r) ;
    MpfrClass& operator/= (const mpz_srcptr r) ;
    MpfrClass& operator/= (const mpq_srcptr r) ;
    static void div (MpfrClass& res, 
                          const MpfrClass& r1, const MpfrClass& r2,
                          RoundingMode rnd = CurrRndMode);

    static void fma (MpfrClass& res, 
                          const MpfrClass& r1, const MpfrClass& r2,
                          const MpfrClass& r3, 
                          RoundingMode rnd = CurrRndMode);

    // Input/Output
    ostream& put(ostream& o,
                 RoundingMode rnd = CurrRndMode,
                 PrecisionType prec = CurrPrecision,
                 int base = 10,
                 int nb_digits = 0) const;
                 //PrecisionType prec = PrecisionType(0) ) const;
    friend istream& operator >> (istream &i, MpfrClass& r);
    friend ostream& operator << (ostream &o, const MpfrClass& r);


    // Mathematical functions: exp, sin...
    void random (PrecisionType prec = CurrPrecision);

    friend MpfrClass abs (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass agm (const MpfrClass& r1, const MpfrClass& r2, RoundingMode rnd = CurrRndMode);
    friend MpfrClass sqrt (const MpfrClass& r, RoundingMode rnd = CurrRndMode);

    friend MpfrClass exp (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass expm1 (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass exp2 (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass log (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass log2 (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass log10 (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass log1p (const MpfrClass& r, RoundingMode rnd = CurrRndMode);

    friend MpfrClass sin (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass cos (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend void sin_cos (MpfrClass& res_sin, MpfrClass& res_cos,
                        const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass tan (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
  
    friend MpfrClass acos (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass asin (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass atan (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
  
    friend MpfrClass cosh (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass sinh (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass tanh (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
  
    friend MpfrClass atanh (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass acosh (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass asinh (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
  
    friend MpfrClass pow (const MpfrClass& r, 
                          const int e, 
                          RoundingMode rnd = CurrRndMode);
    friend MpfrClass pow (const MpfrClass& r, 
                          const unsigned long int e, 
                          RoundingMode rnd = CurrRndMode);
    friend MpfrClass pow (const MpfrClass& r, 
                          const long int e, 
                          RoundingMode rnd = CurrRndMode);
    friend MpfrClass pow (const MpfrClass& r, 
                          const double e,
                          RoundingMode rnd = CurrRndMode);
    friend MpfrClass pow (const MpfrClass& r, 
                          const MpfrClass& e, 
                          RoundingMode rnd = CurrRndMode);
    friend MpfrClass cbrt (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass gamma (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass erf (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend MpfrClass factorial (const unsigned long int e, RoundingMode rnd = CurrRndMode);
    friend MpfrClass hypot (const MpfrClass& r1, const MpfrClass& r2, RoundingMode rnd = CurrRndMode);
    friend MpfrClass zeta (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    //end mathematical functions


    // The four rounding modes as in IEEE-754 arithmetic and the fractional part
    friend MpfrClass round (const MpfrClass& r);
    friend MpfrClass floor (const MpfrClass& r);
    friend MpfrClass trunc (const MpfrClass& r);
    friend MpfrClass ceil (const MpfrClass& r);
    friend MpfrClass frac (const MpfrClass& r);

    friend long int to_int (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend unsigned long int to_uint (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend double to_double (const MpfrClass& r, RoundingMode rnd = CurrRndMode);
    friend long double to_ldouble (const MpfrClass& r, RoundingMode rnd = CurrRndMode);

    friend MpfrClass reldiff (const MpfrClass& r1, const MpfrClass& r2, 
                              RoundingMode rnd = CurrRndMode);
    friend MpfrClass nextabove (const MpfrClass& r);
    friend MpfrClass nextbelow (const MpfrClass& r);
    friend MpfrClass nexttoward (const MpfrClass& r, const MpfrClass& dir);
};
//--------------------------
// Comparison operators
//--------------------------
bool operator == (const MpfrClass& r1, const MpfrClass& r2);
bool operator == (const MpfrClass& r1, const double r2);
bool operator == (const MpfrClass& r1, const int r2);
bool operator == (const MpfrClass& r1, const unsigned int r2);
bool operator == (const MpfrClass& r1, const long int r2);
bool operator == (const MpfrClass& r1, const unsigned long int r2);
bool operator == (const double r1, const MpfrClass& r2);
bool operator == (const int r1, const MpfrClass& r2);
bool operator == (const unsigned int r1, const MpfrClass& r2);
bool operator == (const long int r1, const MpfrClass& r2);
bool operator == (const unsigned long int r1, const MpfrClass& r2);

bool operator != (const MpfrClass& r1, const MpfrClass& r2);
bool operator != (const MpfrClass& r1, const double r2);
bool operator != (const MpfrClass& r1, const int r2);
bool operator != (const MpfrClass& r1, const unsigned int r2);
bool operator != (const MpfrClass& r1, const long int r2);
bool operator != (const MpfrClass& r1, const unsigned long int r2);
bool operator != (const double r1, const MpfrClass& r2);
bool operator != (const int r1, const MpfrClass& r2);
bool operator != (const unsigned int r1, const MpfrClass& r2);
bool operator != (const long int r1, const MpfrClass& r2);
bool operator != (const unsigned long int r1, const MpfrClass& r2);

bool operator <  (const MpfrClass& r1, const MpfrClass& r2);
bool operator <  (const MpfrClass& r1, const double r2);
bool operator <  (const MpfrClass& r1, const int r2);
bool operator <  (const MpfrClass& r1, const unsigned int r2);
bool operator <  (const MpfrClass& r1, const long int r2);
bool operator <  (const MpfrClass& r1, const unsigned long int r2);
bool operator <  (const double r1, const MpfrClass& r2);
bool operator <  (const int r1, const MpfrClass& r2);
bool operator <  (const unsigned int r1, const MpfrClass& r2);
bool operator <  (const long int r1, const MpfrClass& r2);
bool operator <  (const unsigned long int r1, const MpfrClass& r2);

bool operator <= (const MpfrClass& r1, const MpfrClass& r2);
bool operator <= (const MpfrClass& r1, const double r2);
bool operator <= (const MpfrClass& r1, const int r2);
bool operator <= (const MpfrClass& r1, const unsigned int r2);
bool operator <= (const MpfrClass& r1, const long int r2);
bool operator <= (const MpfrClass& r1, const unsigned long int r2);
bool operator <= (const double r1, const MpfrClass& r2);
bool operator <= (const int r1, const MpfrClass& r2);
bool operator <= (const unsigned int r1, const MpfrClass& r2);
bool operator <= (const long int r1, const MpfrClass& r2);
bool operator <= (const unsigned long int r1, const MpfrClass& r2);

bool operator >  (const MpfrClass& r1, const MpfrClass& r2);
bool operator >  (const MpfrClass& r1, const double r2);
bool operator >  (const MpfrClass& r1, const int r2);
bool operator >  (const MpfrClass& r1, const unsigned int r2);
bool operator >  (const MpfrClass& r1, const long int r2);
bool operator >  (const MpfrClass& r1, const unsigned long int r2);
bool operator >  (const double r1, const MpfrClass& r2);
bool operator >  (const int r1, const MpfrClass& r2);
bool operator >  (const unsigned int r1, const MpfrClass& r2);
bool operator >  (const long int r1, const MpfrClass& r2);
bool operator >  (const unsigned long int r1, const MpfrClass& r2);

bool operator >= (const MpfrClass& r1, const MpfrClass& r2);
bool operator >= (const MpfrClass& r1, const double r2);
bool operator >= (const MpfrClass& r1, const int r2);
bool operator >= (const MpfrClass& r1, const unsigned int r2);
bool operator >= (const MpfrClass& r1, const long int r2);
bool operator >= (const MpfrClass& r1, const unsigned long int r2);
bool operator >= (const double r1, const MpfrClass& r2);
bool operator >= (const int r1, const MpfrClass& r2);
bool operator >= (const unsigned int r1, const MpfrClass& r2);
bool operator >= (const long int r1, const MpfrClass& r2);
bool operator >= (const unsigned long int r1, const MpfrClass& r2);

//--------------------------------------------------------------------
// Bitwise XOR used as a power operator (to match the "usual" syntax).
//--------------------------------------------------------------------
inline MpfrClass operator^(const MpfrClass& r, 
                            const int e)
{ return pow(r,e); }

inline MpfrClass operator^(const MpfrClass& r, 
                            const unsigned long int e)
{ return pow(r,e); }

inline MpfrClass operator^(const MpfrClass& r, 
                            const long int e)
{ return pow(r,e); }
    
inline MpfrClass operator^(const MpfrClass& r, 
                            const double e)
{ return pow(r,e); }

inline MpfrClass operator^(const MpfrClass& r, 
                            const MpfrClass& e)
{ return pow(r,e); }

//------------
// min and max
//------------
MpfrClass min (const MpfrClass& a, const MpfrClass& b);
MpfrClass max (const MpfrClass& a, const MpfrClass& b);


#endif		// _MPFR_CLASS_H_




