#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netcdf.h>
#include <tio.h>
#include <tell.h>

#include "config.h"
#include "bpix.h"

#define BITMASK(bit)  (0x01 << (bit))

struct Badpix_Map_Occur_Type
{
   Badpix_Map_Occur_Type *next;
   int *num_occur;
   int num_rows;
   int num_cols;
   int bit;          /* bit=0 is least significant */
};

void bpix_free (Badpix_Map_Type *b)
{
   if (b == NULL)
     return;
   FREE(b->bits);
   FREE(b);
}

Badpix_Map_Type *bpix_new (int num_rows, int num_cols)
{
   Badpix_Map_Type *b = NULL;
   size_t map_size = num_rows * num_cols * sizeof(Badpix_Bitmap_Type);

   if (NULL == (b = (Badpix_Map_Type *)MALLOC (sizeof *b)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   b->bits = NULL;
   if (NULL == (b->bits = (Badpix_Bitmap_Type *)MALLOC (map_size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        bpix_free (b);
        return NULL;
     }

   memset ((char *)b->bits, 0, map_size);
   b->num_rows = num_rows;
   b->num_cols = num_cols;

   return b;
}

int bpix_logical_or (Badpix_Map_Type *a, const Badpix_Bitmap_Type *bits)
{
   Badpix_Bitmap_Type *a_bits = a->bits;
   int i, n = a->num_rows * a->num_cols;

   for (i = 0; i < n; i++)
     {
        a_bits[i] |= bits[i];
     }

   return 0;
}

Badpix_Map_Type *bpix_read (const char *file)
{
   Badpix_Map_Type *b = NULL;
   TIO_Var_Info_Type info;
   int start[2], count[2];
   int ncid;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return NULL;
   tell_vlog (TELL_MSGTYPE_INFO, 1, "reading %s", file);

   if (0 != TIO_inq_var (ncid, "badpix", &info))
     return NULL;

   if (info.ndims != 2)
     {
        tell_verror (TELL_INVALID_DATA_ERROR,
                     "%s: badpix map has %d dimensions, expected 2",
                     __func__, info.ndims);
        return NULL;
     }

   if (NULL == (b = bpix_new (info.dimlens[0], info.dimlens[1])))
     return NULL;

   start[0] = 0;
   start[1] = 0;
   count[0] = b->num_rows;
   count[1] = b->num_cols;
   if (0 != TIO_get_var_section (ncid, "badpix", start, count, BADPIX_BITMAP_TIO_TYPE,
                                 b->bits))
     {
        bpix_free (b);
        return NULL;
     }

   return b;
}

static void free_occur_type1 (Badpix_Map_Occur_Type *ot)
{
   if (ot == NULL)
     return;
   FREE(ot->num_occur);
   FREE(ot);
}

static void free_occur_type (Badpix_Map_Occur_Type *ot)
{
   while (ot != NULL)
     {
        Badpix_Map_Occur_Type *next = ot->next;
        free_occur_type1 (ot);
        ot = next;
     }
}

void bpix_occur_close (Badpix_Map_Occur_Type *ot)
{
   free_occur_type (ot);
}

static Badpix_Map_Occur_Type *new_occur_type (int num_rows, int num_cols,
                                              int bit)
{
   Badpix_Map_Occur_Type *ot = NULL;
   size_t map_size = num_rows * num_cols * sizeof(int);

   if (NULL == (ot = (Badpix_Map_Occur_Type *)MALLOC (sizeof *ot)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   ot->num_occur = NULL;
   if (NULL == (ot->num_occur = (int *)MALLOC (map_size)))
     {
        free_occur_type1 (ot);
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   ot->next = NULL;
   memset ((char *)ot->num_occur, 0, map_size);
   ot->num_rows = num_rows;
   ot->num_cols = num_cols;
   ot->bit = bit;

   return ot;

}

Badpix_Map_Occur_Type *bpix_occur_open (int num_rows, int num_cols,
                                        Badpix_Bitmap_Type mask)
{
   Badpix_Map_Occur_Type *head = NULL;
   int num_bits = sizeof(mask) * CHAR_BIT;
   int bit;

   for (bit = 0; bit < num_bits; bit++)
     {
        Badpix_Map_Occur_Type *ot;

        if (0 == (mask & BITMASK(bit)))
          continue;

        if (NULL == (ot = new_occur_type (num_rows, num_cols, bit)))
          goto return_error;

        if (head != NULL)
          {
             ot->next = head;
             head = ot;
          }
        else head = ot;
     }

   return head;
return_error:
   free_occur_type (head);
   return NULL;
}

static void increment_bit_counters (const Badpix_Map_Occur_Type *ot,
                                    const Badpix_Bitmap_Type *bits)
{
   Badpix_Bitmap_Type mask = BITMASK(ot->bit);
   int *num_occur = ot->num_occur;
   int i, n = ot->num_rows * ot->num_cols;

   for (i = 0; i < n; i++)
     {
        if (bits[i] & mask) num_occur[i] += 1;
     }
}

int bpix_occur_incr (Badpix_Map_Occur_Type *ot, const Badpix_Bitmap_Type *bits)
{
   for ( ; ot != NULL; ot = ot->next)
     {
        increment_bit_counters (ot, bits);
     }

   return 0;
}

static void set_bits_using_count (const Badpix_Map_Occur_Type *ot,
                                  int num_threshold, Badpix_Bitmap_Type *bits)
{
   Badpix_Bitmap_Type mask = BITMASK(ot->bit);
   int *num_occur = ot->num_occur;
   int i, n = ot->num_rows * ot->num_cols;

   for (i = 0; i < n; i++)
     {
        if (num_occur[i] > num_threshold)
          bits[i] |= mask;
     }
}

/* Record selected bpix bits where they occur frequently */
int bpix_occur_set (const Badpix_Map_Occur_Type *ot,
                    int num_threshold, Badpix_Bitmap_Type mask,
                    Badpix_Bitmap_Type *bits)
{
   for ( ; ot != NULL; ot = ot->next)
     {
        if (mask & BITMASK(ot->bit))
          {
             set_bits_using_count (ot, num_threshold, bits);
          }
     }

   return 0;
}
