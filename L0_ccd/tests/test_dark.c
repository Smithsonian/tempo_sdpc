#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <config.h>
#include <libconfig.h>
#include <image.h>
#include <dark.h>

#include "util.h"

#define NUM_DARKS 100
#define NUM_ROWS 32
#define NUM_COLS 32

#define MAX_TEST_DELTA  0.1

typedef struct
{
   Image_Pixel_Type sdc;
   Image_Pixel_Type fptemp;
   Image_Pixel_Type exptime;
   Image_Pixel_Type pixel_value;
}
Dark_Test_Params;

typedef int Table_Check_Method_Type (Dark_Table_Type *, Dark_Test_Params *,
                                     int, Image_Type *);

static void free_test_params (Dark_Test_Params *tp)
{
   if (tp == NULL)
     return;
   FREE(tp);
}

static Dark_Test_Params *init_test_params (int num_darks)
{
   Dark_Test_Params *tp = NULL;
   double test_delta;
   int i;
   if (NULL == (tp = (Dark_Test_Params *)MALLOC (num_darks * sizeof *tp)))
     return NULL;

   test_delta = MAX_TEST_DELTA / num_darks;

   for (i = 0; i < num_darks; i++)
     {
        Dark_Test_Params *tpi = &tp[i];
        double f = 1.0 + i * test_delta;
        tpi->sdc = f * 1.0;
        tpi->fptemp = 253.0 + i * 0.01;
        tpi->exptime = 2.8;
        tpi->pixel_value = 5.0;
     }

   return tp;
}

static void free_image_array (Image_Type **a, int na)
{
   int i;
   if (a == NULL)
     return;
   for (i = 0; i < na; i++)
     {
        image_free (a[i]);
        a[i] = NULL;
     }
   FREE(a);
}

static Image_Type **alloc_image_array (int num_darks, int num_rows,
                                       int num_cols)
{
   Image_Type **a = NULL;
   size_t array_size = num_darks * sizeof(Image_Type *);
   int i;

   if (NULL == (a = (Image_Type **)MALLOC (array_size)))
     return NULL;
   memset ((char *)a, 0, array_size);

   for (i = 0; i < num_darks; i++)
     {
        if (NULL == (a[i] = image_new (num_rows, num_cols)))
          {
             free_image_array (a, num_darks);
             return NULL;
          }
     }

   return a;
}

static int check_table_fptemp (Dark_Table_Type *dtt, Dark_Test_Params *tp,
                               int num_darks, Image_Type *dark_img)
{
   double fptemp_min, fptemp_max;
   int i;

   if (dtt->dtt_ordering (dtt) != DARK_TABLE_ORDERED_BY_TEMP)
     return -1;

   dtt->dtt_domain (dtt, &fptemp_min, &fptemp_max);
   if ((tp[0].fptemp != (Image_Pixel_Type) fptemp_min)
       || (tp[num_darks-1].fptemp != (Image_Pixel_Type) fptemp_max))
     {
        fprintf (stderr, "*** unexpected dark table fptemp domain limits\n");
        return -1;
     }

   for (i = 0; i < num_darks; i++)
     {
        Dark_Test_Params *tpi = &tp[i];
        if (0 != dtt->dtt_interp (dtt, tpi->fptemp, dark_img))
          return -1;
        if (0 == image_test (dark_img, tpi->pixel_value, 0))
          return -1;
     }

   return 0;
}

static int check_table_sdc (Dark_Table_Type *dtt, Dark_Test_Params *tp,
                            int num_darks, Image_Type *dark_img)
{
   double sdc_min, sdc_max;
   int i;

   if (dtt->dtt_ordering (dtt) != DARK_TABLE_ORDERED_BY_SDC)
     return -1;

   dtt->dtt_domain (dtt, &sdc_min, &sdc_max);
   if ((tp[0].sdc != (Image_Pixel_Type) sdc_min)
       || (tp[num_darks-1].sdc != (Image_Pixel_Type) sdc_max))
     {
        fprintf (stderr, "*** unexpected dark table SDC domain limits\n");
        return -1;
     }

   for (i = 0; i < num_darks; i++)
     {
        Dark_Test_Params *tpi = &tp[i];
        if (0 != dtt->dtt_interp (dtt, tpi->sdc, dark_img))
          return -1;
        if (0 == image_test (dark_img, tpi->pixel_value, 0))
          return -1;
     }

   return 0;
}

static int check_table_exptime (Dark_Table_Type *dtt, Dark_Test_Params *tp,
                                int num_darks, Image_Type *dark_img)
{
   double exptime_min, exptime_max;
   int i;

   if (dtt->dtt_ordering (dtt) != DARK_TABLE_ORDERED_BY_EXPTIME)
     return -1;

   dtt->dtt_domain (dtt, &exptime_min, &exptime_max);
   if ((tp[0].exptime != (Image_Pixel_Type) exptime_min)
       || (tp[num_darks-1].exptime != (Image_Pixel_Type) exptime_max))
     {
        fprintf (stderr, "*** unexpected dark table exptime domain limits\n");
        return -1;
     }

   for (i = 0; i < num_darks; i++)
     {
        Dark_Test_Params *tpi = &tp[i];
        if (0 != dtt->dtt_interp (dtt, tpi->exptime, dark_img))
          return -1;
        if (0 == image_test (dark_img, tpi->pixel_value, 0))
          return -1;
     }

   return 0;
}

static int test_dark (config_t *cfg,
                      Table_Check_Method_Type *table_check_method)
{
   const char table_file[] = "dark_table_test.nc";
   Dark_Array_Type *da = NULL;
   Dark_Test_Params *tp = NULL;
   Dark_Test_Params *tpi = NULL;
   Dark_Config_Type *dcfg = NULL;
   Dark_Table_Type *dtt = NULL;
   Image_Type *dark_img = NULL;
   Image_Type **img_array = NULL;
   int num_darks = NUM_DARKS;
   int num_rows = NUM_ROWS;
   int num_cols = NUM_COLS;
   int is_linearity = 0;
   int i, status = -1;

   if (NULL == (tp = init_test_params (num_darks)))
     goto return_status;

   if (NULL == (img_array = alloc_image_array (num_darks, num_rows, num_cols)))
     goto return_status;

   for (i = 0 ; i < num_darks; i++)
     {
        tpi = &tp[i];
        image_set (img_array[i], tpi->pixel_value, 0);
     }

   if (NULL == (da = dark_array_alloc (num_darks)))
     goto return_status;

   for (i = 0 ; i < num_darks; i++)
     {
        tpi = &tp[i];
        if (0 != dark_array_elem_set (da, i, img_array[i],
                                      tpi->sdc, tpi->fptemp, tpi->exptime))
          goto return_status;
     }

   if (NULL == (dcfg = dark_table_config (cfg, is_linearity)))
     goto return_status;

   if (NULL == (dtt = dark_table_create (dcfg, da)))
     goto return_status;

   if (0 != dtt->dtt_write (dtt, table_file))
     goto return_status;

   dtt->dtt_delete (dtt);
   dtt = NULL;

   if (NULL == (dtt = dark_table_read (table_file)))
     goto return_status;

   if (NULL == (dark_img = image_new (num_rows, num_cols)))
     goto return_status;

   if (0 != (*table_check_method)(dtt, tp, num_darks, dark_img))
     goto return_status;

   status = 0;
return_status:

   image_free (dark_img);
   dark_table_config_free (dcfg);
   if (dtt) dtt->dtt_delete (dtt);
   dark_array_free (da);
   free_image_array (img_array, num_darks);
   free_test_params (tp);

   return status;
}

static int perform_test (config_t *cfg, int sort_order)
{
   Table_Check_Method_Type *table_check_method;

   switch (sort_order)
     {
      case DARK_TABLE_ORDERED_BY_TEMP:
        table_check_method = check_table_fptemp;
        break;
      case DARK_TABLE_ORDERED_BY_SDC:
        table_check_method = check_table_sdc;
        break;
      case DARK_TABLE_ORDERED_BY_EXPTIME:
        table_check_method = check_table_exptime;
        break;
      default:
        fprintf (stderr, "*** unsupported sort order =%d\n", sort_order);
        return -1;
        break;
     }

   return test_dark (cfg, table_check_method);
}

int main (int argc, char **argv)
{
#define BUFSIZE 256
   char config_file[BUFSIZE];
   const char *sort_order_string;
   config_t cfg;
   int sort_order;
   int status = EXIT_FAILURE;

   if (argc != 2)
     {
        fprintf (stderr, "Usage:  %s <sort-order>\n", argv[0]);
        return 0;
     }

   sort_order_string = argv[1];

   if (0 == strcmp (sort_order_string, "fptemp"))
     sort_order = DARK_TABLE_ORDERED_BY_TEMP;
   else if (0 == strcmp (sort_order_string, "sdc"))
     sort_order = DARK_TABLE_ORDERED_BY_SDC;
   else if (0 == strcmp (sort_order_string, "exptime"))
     sort_order = DARK_TABLE_ORDERED_BY_EXPTIME;
   else
     {
        fprintf (stderr, "**** unsupported sort order: %s\n", sort_order_string);
        return 0;
     }

   sprintf (config_file, "test_dark_%s.cfg", argv[1]);

   config_init (&cfg);

   if (0 == config_read_file (&cfg, config_file))
     {
        fprintf (stderr, "%s: Reading %s:%d - %s\n",
                 __func__, config_error_file(&cfg),
                 config_error_line(&cfg), config_error_text(&cfg));
        goto return_status;
     }

   if (0 == perform_test (&cfg, sort_order))
     status = EXIT_SUCCESS;

return_status:
   config_destroy (&cfg);
   if (status)
     fprintf (stderr, "*** ERROR: test_dark failed\n");

   return status;
}

