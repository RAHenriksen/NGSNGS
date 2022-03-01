#if defined(__APPLE__) && defined(__MACH__) 
#include "NGSNGS_Random.h"
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <machine/endian.h>

# define weak_alias(name, aliasname) _weak_alias (name, aliasname)
# define _weak_alias(name, aliasname) \
  extern __typeof (name) aliasname __attribute__ ((weak, alias (#name)));

/* Copyright (C) 1995-2019 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@gnu.ai.mit.edu>, August 1995.
   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */

//https://code.woboq.org/userspace/glibc/stdlib/drand48-iter.c.html#__drand48_iterate
int
__drand48_iterate (unsigned short int xsubi[3], struct drand48_data *buffer)
{
  uint64_t X;
  uint64_t result;
  
  /* Do the real work.  We choose a data type which contains at least
     48 bits.  Because we compute the modulus it does not care how
     many bits really are computed.  */
  X = (uint64_t) xsubi[2] << 32 | (uint32_t) xsubi[1] << 16 | xsubi[0];
  result = X * buffer->__a + buffer->__c;
  xsubi[0] = result & 0xffff;
  xsubi[1] = (result >> 16) & 0xffff;
  xsubi[2] = (result >> 32) & 0xffff;
  return 0;
}
//https://code.woboq.org/userspace/glibc/sysdeps/ieee754/ieee754.h.html#ieee754_double
union ieee754_double
  {
    double d;
    /* This is the IEEE 754 double-precision format.  */
    struct
      {
/*#if        __BYTE_ORDER == __BIG_ENDIAN
        unsigned int negative:1;
        unsigned int exponent:11;*/
        /* Together these comprise the mantissa.  */
        /*unsigned int mantissa0:20;
        unsigned int mantissa1:32;
#endif*/                                /* Big endian.  */
#if        __BYTE_ORDER == __LITTLE_ENDIAN
# if        __FLOAT_WORD_ORDER == __BIG_ENDIAN
        unsigned int mantissa0:20;
        unsigned int exponent:11;
        unsigned int negative:1;
        unsigned int mantissa1:32;
# else
        /* Together these comprise the mantissa.  */
        unsigned int mantissa1:32;
        unsigned int mantissa0:20;
        unsigned int exponent:11;
        unsigned int negative:1;
# endif
#endif                                /* Little endian.  */
      } ieee;
    /* This format makes it easier to see if a NaN is a signalling NaN.  */
    struct
      {
#if        __BYTE_ORDER == __BIG_ENDIAN
        unsigned int negative:1;
        unsigned int exponent:11;
        unsigned int quiet_nan:1;
        /* Together these comprise the mantissa.  */
        unsigned int mantissa0:19;
        unsigned int mantissa1:32;
#else
# if        __FLOAT_WORD_ORDER == __BIG_ENDIAN
        unsigned int mantissa0:19;
        unsigned int quiet_nan:1;
        unsigned int exponent:11;
        unsigned int negative:1;
        unsigned int mantissa1:32;
# else
        /* Together these comprise the mantissa.  */
        unsigned int mantissa1:32;
        unsigned int mantissa0:19;
        unsigned int quiet_nan:1;
        unsigned int exponent:11;
        unsigned int negative:1;
# endif
#endif
      } ieee_nan;
  };

#define IEEE754_DOUBLE_BIAS     0x3ff /* Added to exponent.  */

//https://code.woboq.org/userspace/glibc/stdlib/erand48_r.c.html#47
int
erand48_r (unsigned short int xsubi[3], struct drand48_data *buffer,
             double *result)
{
  union ieee754_double temp;
  /* Compute next state.  */
  if (__drand48_iterate (xsubi, buffer) < 0)
    return -1;
  /* Construct a positive double with the 48 random bits distributed over
     its fractional part so the resulting FP number is [0.0,1.0).  */
  temp.ieee.negative = 0;
  temp.ieee.exponent = IEEE754_DOUBLE_BIAS;
  temp.ieee.mantissa0 = (xsubi[2] << 4) | (xsubi[1] >> 12);
  temp.ieee.mantissa1 = ((xsubi[1] & 0xfff) << 20) | (xsubi[0] << 4);
  /* Please note the lower 4 bits of mantissa1 are always 0.  */
  *result = temp.d - 1.0;
  return 0;
}

//https://code.woboq.org/userspace/glibc/stdlib/drand48_r.c.html#drand48_r
int
drand48_r (struct drand48_data *buffer, double *result)
{
  return erand48_r (buffer->__x, buffer, result);
}

//modified the name instead of using the weak alias
//https://code.woboq.org/userspace/glibc/stdlib/srand48_r.c.html#srand48_r
int
srand48_r (long int seedval, struct drand48_data *buffer)
{
  /* The standards say we only have 32 bits.  */
  if (sizeof (long int) > 4)
    seedval &= 0xffffffffl;
  buffer->__x[2] = seedval >> 16;
  buffer->__x[1] = seedval & 0xffffl;
  buffer->__x[0] = 0x330e;
  buffer->__a = 0x5deece66dull;
  buffer->__c = 0xb;
  buffer->__init = 1;
  return 0;
}

#endif