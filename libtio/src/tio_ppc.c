#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>

#include <tell.h>
#include "tio_ppc.h"

/* Adapated from nco_ppc.[ch],
 * NCO utilities for Precision-Preserving Compression (PPC)
 * by Charlie Zender
 * You may redistribute and/or modify this code under the terms of the
 * GNU General Public License (GPL) Version 3 with exceptions described
 * in the NCO LICENSE file
 */

/* Minimum number of explicit significand bits to preserve when
 * zeroing/bit-masking floating point values.  Codes will preserve
 * at least two explicit bits, IEEE significand representation
 * contains one implicit bit. Thus preserve a least three bits which
 * is approximately one sigificant decimal digit. */
#define NCO_PPC_BIT_XPL_NBR_MIN 2

/* Precision-preserving compression: mask-out insignificant bits of significand
 * nsd            number of significant digits
 * ppc_method     bit-masking method
 * sz             number of float values to process
 * val            pointer to array float values
 * pmiss_val      pointer to float value used to indicate missing data;
 *                or NULL if no such value is used.
 */
int _pTIO_ppc_f32_bitmask (int nsd, int ppc_method, int sz, float *val,
                           const float *pmiss_val)
{
   /* 3.32 bits per decimal digit of precision */
   double bit_per_dcm_dgt_prc = M_LN10/M_LN2;

   /* Bits 0-22 of single-precision significands are explicit.
    * Bit 23 is implicitly 1. */
   int bit_xpl_nbr_sgn = 23;

   /* Binary digits of precision, exact */
   double prc_bnr_xct;

   /* Number of explicit bits to zero */
   int bit_xpl_nbr_zro;

   /* Exact binary digits of precision rounded-up */
   unsigned short prc_bnr_ceil;

   /* Explicitly represented binary digits required to retain */
   unsigned short prc_bnr_xpl_rqr;

   long idx;
   unsigned int *u32_ptr;
   unsigned int msk_f32_u32_zro;
   unsigned int msk_f32_u32_one;

   union {
      float *fp;
      unsigned int *uip;
   } pval;

   if (val == NULL)
     {
        Tell_verror (TELL_UNKNOWN_ERROR, "%s: NULL value pointer", __func__);
        return -1;
     }

   /* Disallow unreasonable quantization */
   if (nsd <= 0 || nsd > 7)
     {
        Tell_verror (TELL_USAGE_ERROR,
                     "%s: request for unreasonable number of significant digits, nsd = %d",
                     __func__, nsd);
        return -1;
     }

   /* Use a union so we're twiddling bits on an unsigned int */
   pval.fp = val;
   u32_ptr = pval.uip;

   /* How many bits to preserve? */
   prc_bnr_xct = nsd * bit_per_dcm_dgt_prc;
   /* Be conservative, round upwards */
   prc_bnr_ceil = (unsigned short)ceil(prc_bnr_xct);
   /* The first bit is implicit not explicit but
    * corner cases prevent taking advantage of this */
   prc_bnr_xpl_rqr = prc_bnr_ceil+1;

   if (prc_bnr_xpl_rqr >= bit_xpl_nbr_sgn)
     return 0;

   if (prc_bnr_xpl_rqr < NCO_PPC_BIT_XPL_NBR_MIN)
     return 0;

   bit_xpl_nbr_zro = bit_xpl_nbr_sgn - prc_bnr_xpl_rqr;

   /* Create mask */
   msk_f32_u32_zro = 0u;               /* Zero all bits */
   msk_f32_u32_zro = ~msk_f32_u32_zro; /* Flip all bits to one */

   /* Bit Shave mask for AND: Left shift zeros into bits to be rounded,
    * leave ones in untouched bits */
   msk_f32_u32_zro <<= bit_xpl_nbr_zro;

   /* Bit Set mask for OR:  Put ones into bits to be set,
    * leave zeros in untouched bits */
   msk_f32_u32_one = ~msk_f32_u32_zro;

   switch (ppc_method)
     {
      default:
        Tell_verror (TELL_USAGE_ERROR, "%s: invalid ppc method = %d",
                     __func__, ppc_method);
        return -1;

      case _pTIO_PPC_METHOD_ALT:
        /* Alternately set and unset LSBs */
        if (pmiss_val == NULL)
          {
             for (idx = 0; idx<sz; idx += 2)
               u32_ptr[idx] &= msk_f32_u32_zro;

             for (idx = 1; idx<sz; idx += 2)
               {
                  /* Never quantize upwards floating
                   * point values of zero */
                  if (u32_ptr[idx] != 0U)
                    u32_ptr[idx] |= msk_f32_u32_one;
               }
          }
        else
          {
             float miss_val = *pmiss_val;
             for (idx = 0; idx<sz; idx += 2)
               {
                  if (val[idx] != miss_val)
                    u32_ptr[idx] &= msk_f32_u32_zro;
               }

             for (idx = 1; idx<sz; idx += 2)
               {
                  /* Never quantize upwards floating
                   * point values of zero */
                  if (val[idx] != miss_val && u32_ptr[idx] != 0U)
                    u32_ptr[idx] |= msk_f32_u32_one;
               }
          }
        break;

      case _pTIO_PPC_METHOD_ZERO:
        /* zero LSBs */
        if (pmiss_val == NULL)
          {
             for (idx = 0; idx<sz; idx++)
               u32_ptr[idx] &= msk_f32_u32_zro;
          }
        else
          {
             float miss_val = *pmiss_val;
             for (idx = 0; idx<sz; idx++)
               {
                  if (val[idx] != miss_val)
                    u32_ptr[idx] &= msk_f32_u32_zro;
               }
          }
        break;

      case _pTIO_PPC_METHOD_SET:
        /* set LSBs */
        if (pmiss_val == NULL)
          {
             for (idx = 0; idx<sz; idx++)
               {
                  /* Never quantize upwards floating
                   * point values of zero */
                  if (u32_ptr[idx] != 0U)
                    u32_ptr[idx] |= msk_f32_u32_one;
               }
          }
        else
          {
             float miss_val = *pmiss_val;
             for (idx = 0; idx<sz; idx++)
               {
                  if (val[idx] != miss_val)
                    u32_ptr[idx] |= msk_f32_u32_one;
               }
          }
        break;
     }

   return 0;
}

#ifdef SELF_TEST

/* gcc -o ppc_test -DSELF_TEST -W -Wall tio_ppc.c -lm */

static void uint_char_bits (unsigned int u, char *buf)
{
   int i;
   for (i = 0; i < 32; i++)
     {
        buf[31-i] = (u & (1 << i)) ? '1' : '0';
     }
   buf[32] = 0;
}

static void print_float_bits (float x, const char *s)
{
   char buf[33];
   unsigned int *ux = (unsigned int *)&x;
   uint_char_bits (*ux, buf);
   fprintf (stdout, "%33s  %+0.7e %s\n", buf, x, s ? s : "");
}

static void do_test1 (float test_value)
{
   int nsd;

   for (nsd = 1; nsd < 8; nsd++)
     {
        float x;
        fprintf (stdout, "%d significant digit%s\n",
                 nsd, (nsd > 1) ? "s" : "");

        x = test_value;
        print_float_bits (x, "test value");

        _pTIO_ppc_f32_bitmask (nsd, _pTIO_PPC_METHOD_ZERO, 1, &x, NULL);
        print_float_bits (x, "zero LSB");

        x = test_value;
        _pTIO_ppc_f32_bitmask (nsd, _pTIO_PPC_METHOD_SET, 1, &x, NULL);
        print_float_bits (x, "set LSB");
     }
}

static void do_test2 (int nsd, int method, float test_value, float *miss_value)
{
#define NUM 6
   float ary[NUM];
   int i, n = NUM;
   int sgn = 1;
   char *method_label;

   ary[0] = miss_value ? *miss_value : test_value;
   for (i = 1; i < n; i++)
     {
        ary[i] = sgn * (test_value + i * M_LOG10E) * 1.e11;
        sgn = -1 * sgn;
     }

   fprintf (stdout, "\nInput array values%s\n",
            miss_value ? " (first indicates 'missing data')" : "");
   for (i = 0; i < n; i++)
     {
        print_float_bits (ary[i], NULL);
     }

   switch (method)
     {
      default:
        method = _pTIO_PPC_METHOD_ALT;
        /* drop */
      case _pTIO_PPC_METHOD_ALT:
        method_label = "_pTIO_PPC_METHOD_ALT";
        break;
      case _pTIO_PPC_METHOD_SET:
        method_label = "_pTIO_PPC_METHOD_SET";
        break;
      case _pTIO_PPC_METHOD_ZERO:
        method_label = "_pTIO_PPC_METHOD_ZERO";
        break;
     }

   _pTIO_ppc_f32_bitmask (nsd, method, n, ary, miss_value);
   fprintf (stdout, "%s array values (%d significant digits)\n",
            method_label, nsd);
   for (i = 0; i < n; i++)
     {
        print_float_bits (ary[i], NULL);
     }
}

int main (void)
{
   float test_value = -M_PI * 1.e-8;
   float missing = FLT_MAX;

   do_test1 (test_value);
   do_test2 (2, _pTIO_PPC_METHOD_ALT, test_value, &missing);
   do_test2 (4, _pTIO_PPC_METHOD_ALT, test_value, &missing);
   do_test2 (4, _pTIO_PPC_METHOD_ALT, test_value, NULL);
   do_test2 (2, _pTIO_PPC_METHOD_SET, test_value, &missing);
   do_test2 (4, _pTIO_PPC_METHOD_SET, test_value, NULL);
   do_test2 (2, _pTIO_PPC_METHOD_ZERO, test_value, &missing);
   do_test2 (4, _pTIO_PPC_METHOD_ZERO, test_value, NULL);

   return 0;
}

#endif
