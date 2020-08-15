#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tell.h>

#include "config.h"
#include "image.h"

void image_free (Image_Type *img)
{
   if (img == NULL)
     return;
   FREE(img->pixels);
   FREE(img->pixel_quality_flags);
   FREE(img);
}

Image_Type *image_new (int num_rows, int num_cols)
{
   Image_Type *img = NULL;
   size_t num_pixels = num_rows * num_cols;

   if ((NULL == (img = (Image_Type *) MALLOC (sizeof *img)))
       || (NULL == (img->pixels = (Image_Pixel_Type *) MALLOC (num_pixels * sizeof(Image_Pixel_Type))))
       || (NULL == (img->pixel_quality_flags = (Image_Pqf_Bitmap_Type *) MALLOC (num_pixels * sizeof(Image_Pqf_Bitmap_Type))))
      )
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        image_free (img);
        return NULL;
     }

   memset ((char *)img->pixels, 0,
           num_pixels * sizeof(Image_Pixel_Type));
   memset ((char *)img->pixel_quality_flags, 0,
           num_pixels * sizeof(Image_Pqf_Bitmap_Type));

   img->num_cols = num_cols;
   img->num_rows = num_rows;
   img->image_type = IMAGE_TYPE_UNKNOWN;

   return img;
}

Image_Type *image_dup (const Image_Type *img)
{
   Image_Type *dup;
   int num_pixels;

   if (img == NULL)
     return NULL;

   if (NULL == (dup = image_new (img->num_rows, img->num_cols)))
     return NULL;

   num_pixels = img->num_rows * img->num_cols;

   dup->image_type = img->image_type;
   memcpy ((char *)dup->pixels,
           (char *)img->pixels, num_pixels * sizeof(Image_Pixel_Type));
   memcpy ((char *)dup->pixel_quality_flags,
           (char *)img->pixel_quality_flags, num_pixels * sizeof(Image_Pqf_Bitmap_Type));

   return dup;
}

int image_copy (const Image_Type *src, Image_Type *dest)
{
   int num_pixels;

   if ((src == NULL) || (dest == NULL))
     return -1;

   /* sizes must match */
   if ((dest->num_rows != src->num_rows)
       || (dest->num_cols != src->num_cols))
     return -1;

   num_pixels = src->num_rows * src->num_cols;

   dest->image_type = src->image_type;
   memcpy ((char *)dest->pixels,
           (char *)src->pixels, num_pixels * sizeof(Image_Pixel_Type));
   memcpy ((char *)dest->pixel_quality_flags,
           (char *)src->pixel_quality_flags, num_pixels * sizeof(Image_Pqf_Bitmap_Type));

   return 0;
}

void image_sqrt (Image_Type *img)
{
   Image_Pixel_Type *pixels;
   int i, num_pixels;

   num_pixels = img->num_rows * img->num_cols;
   pixels = img->pixels;

   for (i = 0; i < num_pixels; i++)
     {
        Image_Pixel_Type pix_i = pixels[i];
        pixels[i] = (pix_i < 0) ? IMAGE_PIXEL_FILL_VALUE : sqrt(pix_i);
     }
}

int image_check_negative_pixels (Image_Type *img, int flag)
{
   Image_Pixel_Type *pixels = img->pixels;
   Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags;
   int i, n, count;

   n = img->num_rows * img->num_cols;

   count = 0;
   for (i = 0; i < n; i++)
     {
        if ((pixels[i] < 0) && (pqf[i] == 0))
          {
             count++;
             if (flag)
               {
                  pqf[i] |= IMAGE_PQF_PROCESSING_ERROR;
               }
          }
     }

   return count;
}

int image_count_mask_pixels (Image_Type *img, unsigned int mask)
{
   Image_Pqf_Bitmap_Type *pqf = img->pixel_quality_flags;
   int i, n, count;

   n = img->num_rows * img->num_cols;

   count = 0;
   for (i = 0; i < n; i++)
     {
        if (pqf[i] & mask) count++;
     }

   return count;
}

int image_transfer_pqf (const Image_Type *a, Image_Type *b)
{
   const Image_Pqf_Bitmap_Type *a_pqf = a->pixel_quality_flags;
   Image_Pqf_Bitmap_Type *b_pqf = b->pixel_quality_flags;
   int i, n;

   if ((a->num_rows != b->num_rows)
       || (a->num_cols != b->num_cols))
     return -1;

   n = a->num_rows * a->num_cols;

   for (i = 0; i < n; i++)
     {
        b_pqf[i] |= a_pqf[i];
     }

   return 0;
}

void image_scale (Image_Type *img, double s)
{
   Image_Pixel_Type *pixels = img->pixels;
   int i, n = img->num_rows * img->num_cols;

   if (isfinite(s))
     {
        for (i = 0; i < n; i++)
          {
             if (pixels[i] != IMAGE_PIXEL_FILL_VALUE)
               pixels[i] *= s;
          }
     }
   else
     {
        for (i = 0; i < n; i++)
          {
             pixels[i] = IMAGE_PIXEL_FILL_VALUE;
          }
     }
}

int image_get_type (const Image_Type *img)
{
   return (img != NULL) ? img->image_type : IMAGE_TYPE_UNKNOWN;
}

void image_set_type (Image_Type *img, int image_type)
{
   if (img == NULL)
     return;

   img->image_type = image_type;
}

Image_Subset_Type *image_new_subsets (int num_subsets)
{
   Image_Subset_Type *img_subset = NULL;
   size_t size = num_subsets * sizeof(Image_Subset_Type);

   if (NULL == (img_subset = (Image_Subset_Type *) MALLOC (size)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return NULL;
     }

   memset ((char *)img_subset, 0, size);

   return img_subset;
}

void image_set_subset (Image_Subset_Type *s,
                       int row_beg, int row_end, int row_step,
                       int col_beg, int col_end, int col_step)
{
   s->row_beg = row_beg;
   s->row_end = row_end;
   s->row_step = row_step;
   s->col_beg = col_beg;
   s->col_end = col_end;
   s->col_step = col_step;
}

/* this routine is intended for debugging only */
int image_write_raw (const Image_Type *img, const char *file_basename)
{
#define BUFSIZE 1024
   char pathbuf[BUFSIZE];
   char mode[] = "w";
   FILE *fp;
   size_t count, num;
   int len;

   if ((img == NULL) || (file_basename == NULL))
     {
        fprintf (stderr, "%s:  received NULL pointer: img=%p file_basename=%p\n",
                 __func__, img, file_basename);
        return -1;
     }

   len = snprintf (pathbuf, BUFSIZE, "%s_pixels.arr", file_basename);
   if (len >= BUFSIZE)
     {
        fprintf (stderr, "%s: file path string too long (%s)\n", __func__, file_basename);
        return -1;
     }

   if (NULL == (fp = fopen (pathbuf, mode)))
     {
        fprintf (stderr, "%s: error opening for writing, %s\n", __func__, pathbuf);
        return -1;
     }

   fprintf (stderr, "writing image (rows,cols) = (%d,%d) to %s\n",
            img->num_rows, img->num_cols, pathbuf);
   num = img->num_rows * img->num_cols;
   count = fwrite (img->pixels, sizeof(Image_Pixel_Type), num, fp);
   if (count != num)
     {
        fprintf (stderr, "%s: error writing file %s\n", __func__, pathbuf);
        return -1;
     }

   fclose(fp);

   len = snprintf (pathbuf, BUFSIZE, "%s_pqf.arr", file_basename);
   if (len >= BUFSIZE)
     {
        fprintf (stderr, "%s: file path string too long (%s)\n", __func__, file_basename);
        return -1;
     }

   if (NULL == (fp = fopen (pathbuf, mode)))
     {
        fprintf (stderr, "%s: error opening for writing, %s\n", __func__, pathbuf);
        return -1;
     }

   fprintf (stderr, "writing image (rows,cols) = (%d,%d) to %s\n",
            img->num_rows, img->num_cols, pathbuf);
   num = img->num_rows * img->num_cols;
   count = fwrite (img->pixel_quality_flags, sizeof(Image_Pqf_Bitmap_Type), num, fp);
   if (count != num)
     {
        fprintf (stderr, "%s: error writing file %s\n", __func__, pathbuf);
        return -1;
     }

   fclose(fp);

   return 0;
}
