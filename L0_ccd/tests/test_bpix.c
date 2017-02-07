#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <config.h>
#include <bpix.h>

static void set_bitmap (Badpix_Map_Type *m, Badpix_Bitmap_Type mask)
{
   int i, n = m->num_rows * m->num_cols;

   for (i = 0; i < n; i++)
     {
        m->bits[i] = mask;
     }
}

static int assert_bits_set (const Badpix_Map_Type *m, Badpix_Bitmap_Type mask)
{
   int i, n = m->num_rows * m->num_cols;

   for (i = 0; i < n; i++)
     {
        if (0 == (m->bits[i] & mask))
          return 0;
     }

   return 1;
}

static int assert_equal_badpix_maps (const Badpix_Map_Type *m1,
                                     const Badpix_Map_Type *m2)
{
   int i, num_pixels;

   if ((m1->num_rows != m2->num_rows)
       ||(m1->num_cols != m2->num_cols))
     {
        fprintf (stderr, "*** bpix_read failed: read (%d,%d), expected (%d,%d)\n",
                 m2->num_rows, m2->num_cols, m1->num_rows, m1->num_cols);
        return 0;
     }

   num_pixels = m1->num_rows * m1->num_cols;

   for (i = 0; i < num_pixels; i++)
     {
        if (m1->bits[i] != m2->bits[i])
          return 0;
     }

   return 1;
}

static int perform_test (void)
{
   Badpix_Map_Type *m0 = NULL;
   Badpix_Map_Type *m1 = NULL;
   Badpix_Map_Type *m2 = NULL;
   Badpix_Map_Occur_Type *occur = NULL;
   Badpix_Bitmap_Type mask;
   int num_rows = 10, num_cols = 10;
   int num_bits = 8 * sizeof(Badpix_Bitmap_Type);
   const char bpix_outfile[] = "bpix_test.nc";
   int b, num_occur, status = -1;

   if ((NULL == (m0 = bpix_new (num_rows, num_cols)))
       || (NULL == (m1 = bpix_new (num_rows, num_cols))))
     goto return_status;

   /* count all occurences of any bits set */
   mask = -1;
   if (NULL == (occur = bpix_occur_open (num_rows, num_cols, mask)))
     goto return_status;

   for (b = 0; b < num_bits; b++)
     {
        /* Using m1, set 1 bit at a time */
        set_bitmap (m1, 1 << b);
        /* accumulate bits in m0 */
        if (0 != bpix_logical_or (m0, m1->bits))
          goto return_status;
        /* lowest order bit occurs num_bits times,
         * highest order bit occurs only once */
        if (0 != bpix_occur_incr (occur, m0->bits))
          goto return_status;
        set_bitmap (m1, 0);
     }

   /* test assertions about how often each bit was set */
   for (b = 0; b < num_bits; b++)
     {
        set_bitmap (m1, 0);
        mask = 1 << b;
        num_occur = num_bits - b;
        if (0 != bpix_occur_set (occur, num_occur-1, mask, m1->bits))
          goto return_status;
        if (0 == assert_bits_set (m1, mask))
          goto return_status;
     }

   if (0 != bpix_write (m1, bpix_outfile))
     goto return_status;

   if (NULL == (m2 = bpix_read (bpix_outfile)))
     goto return_status;

   if (0 == assert_equal_badpix_maps (m1, m2))
     goto return_status;

   status = 0;
return_status:
   bpix_free (m0);
   bpix_free (m1);
   bpix_free (m2);
   bpix_occur_close (occur);
   return status;
}

int main (void)
{
   int status = EXIT_FAILURE;

   if (0 == perform_test ())
     status = EXIT_SUCCESS;

   if (status)
     fprintf (stderr, "*** ERROR: test_pixelqf failed\n");

   return status;
}

