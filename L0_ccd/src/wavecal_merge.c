#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <limits.h>

#include <ioclib.h>
#include <tell.h>
#include <tio.h>
#include <tio_template.h>

static void usage (void)
{
   fprintf (stderr, "Usage: wavecal_merge -t <target-file> <wavecal-results-dir>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -t | --target FILE   name of netCDF4 file to receive wavecal arrays\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -d | --delete        delete input files and directory after successful merge\n");
   fprintf (stderr, "   -h | --help          print this usage message\n");
   exit (EXIT_SUCCESS);
}

static int compute_wavelength_grid (float *pars, size_t num_pars, int start_pix, int num_pix,
                                    float *waves, size_t num_waves)
{
   float p0, p1, c0;
   int end_pix = start_pix + num_pix - 1;
   int i, num_waves_signed;

   if (num_pars != 2)
     {
        fprintf (stderr, "*** Error:  unsupported polynomial degree = %ld (%s)\n",
                 num_pars, __func__);
        return -1;
     }

   p0 = pars[0];
   p1 = pars[1];
   c0 = 2.0 / (end_pix - start_pix);

   /* we need a signed value of (i-start_pix) */
   num_waves_signed = num_waves;

   for (i = 0; i < num_waves_signed; i++)
     {
        float x =  c0 * (i - start_pix) - 1.0;
        waves[i] = p0 + x * p1;
     }

   return 0;
}

static int replace_wavelength_grid (int grp_target, int mirror_step, int start_pix, int num_pix,
                                    float *param_slab_i, size_t xtrack_dimlen, size_t params_dimlen,
                                    size_t spectral_channel_dimlen,
                                    float *calibrated_waves)
{
   int start[3], count[3];
   size_t xtrack;

   for (xtrack = 0; xtrack < xtrack_dimlen; xtrack++)
     {
        float *pars = param_slab_i + xtrack * params_dimlen;

        if (pars[0] == NC_FILL_FLOAT)
          continue;

        if (0 != compute_wavelength_grid (pars, params_dimlen, start_pix, num_pix,
                                          calibrated_waves, spectral_channel_dimlen))
          return -1;

        start[0] = mirror_step;
        start[1] = xtrack;
        start[2] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = spectral_channel_dimlen;

        if (0 != TIO_put_var_section (grp_target, TEMPO_VAR_WAVELENGTH, start, count,
                                      TIO_FLOAT, calibrated_waves))
          return -1;
     }

   return 0;
}

static int perform_merge (int ncid_target, const char *file)
{
   const char *params_varname = "wavecal_params";
   TIO_Var_Info_Type info = {0};
   char group_name[TIO_MAX_NAME_LEN] = {0};
   int ncid_src, grp_target, varid, start_pix, num_pix;
   int step_dimlen_src, step_dimid, xtrack_dimid, dest_varid, dimid;
   int start[3], count[3];
   size_t step_dimlen, xtrack_dimlen, xtrack_dimlen_src, spectral_channel_dimlen;
   size_t params_dimlen_src, len_params, len_slab, i;
   float *wavecal_params = NULL;
   float *calibrated_waves = NULL;
   int *mirror_step = NULL;
   int status = -1;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid_src))
     return -1;

   /* read group name and verify array dimensions */
   if ((0 != TIO_get_att (ncid_src, NC_GLOBAL, "group_name", NC_CHAR, group_name))
       || (0 != TIO_get_att (ncid_src, NC_GLOBAL, "mirror_step_dimlen", NC_INT, &step_dimlen_src)))
     goto close_and_return;

   if (0 != TIO_inq_var (ncid_src, params_varname, &info))
     goto close_and_return;

   xtrack_dimlen_src = info.dimlens[1];
   params_dimlen_src = info.dimlens[2];

   if (0 != TIO_inq_dim (ncid_target, TEMPO_DIM_STEP, &step_dimid, &step_dimlen))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading step dimension size", __func__);
        goto close_and_return;
     }

   if (0 != TIO_inq_grp (ncid_target, group_name, &grp_target))
     goto close_and_return;

   if (0 != TIO_inq_dim (grp_target, TEMPO_DIM_XTRACK, &xtrack_dimid, &xtrack_dimlen))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading xtrack dimension size", __func__);
        goto close_and_return;
     }

   if (0 != TIO_inq_dim (grp_target, TEMPO_DIM_CHANNEL, &dimid, &spectral_channel_dimlen))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading spectral channel  dimension size", __func__);
        goto close_and_return;
     }
   if (NULL == (calibrated_waves = (float *)MALLOC (spectral_channel_dimlen * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   if (step_dimlen != (size_t) step_dimlen_src)
     {
        tell_vwarn (0, "dimension mismatch (mirror_step) source:%d target:%ld",
                    step_dimlen_src, step_dimlen);
        goto close_and_return;
     }

   if (xtrack_dimlen != xtrack_dimlen_src)
     {
        tell_vwarn (0, "dimension mismatch (xtrack) source:%ld target:%ld",
                    xtrack_dimlen_src, xtrack_dimlen);
        goto close_and_return;
     }

   /* read params */

   len_params = step_dimlen_src * xtrack_dimlen_src * params_dimlen_src;
   if ((NULL == (wavecal_params = MALLOC (len_params * sizeof(float))))
       ||(NULL == (mirror_step = MALLOC (step_dimlen_src * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = info.dimlens[0];
   count[1] = xtrack_dimlen_src;
   count[2] = params_dimlen_src;

   if ((0 != TIO_get_var_section (ncid_src, TEMPO_DIM_STEP, start, count,
                                  TIO_INT, mirror_step))
       ||(0 != TIO_get_var_section (ncid_src, params_varname, start, count,
                                    TIO_FLOAT, wavecal_params)))
     goto close_and_return;

   if ((0 != tio_inq_varid (ncid_src, params_varname, &varid))
       || (0 != TIO_get_att (ncid_src, varid, "start_spectral_channel", NC_INT, &start_pix))
       || (0 != TIO_get_att (ncid_src, varid, "num_spectral_channels", NC_INT, &num_pix)))
     goto close_and_return;

   /* write params, creating target variable if necessary */

   if (0 != tio_inq_varid (grp_target, params_varname, &dest_varid))
     {
        const char *params_dimname = "wavecal_par";
        size_t params_dimlen;
        int params_dimid, params_dimid_list[3];
        if (0 != TIO_inq_dim (grp_target, params_dimname, &params_dimid, &params_dimlen))
          {
             if (0 != TIO_def_dim (grp_target, params_dimname, params_dimlen_src, &params_dimid))
               goto close_and_return;
          }
        params_dimid_list[0] = step_dimid;
        params_dimid_list[1] = xtrack_dimid;
        params_dimid_list[2] = params_dimid;

        if (0 != TIO_def_var (grp_target, params_varname, TIO_FLOAT, 3, params_dimid_list, &dest_varid))
          goto close_and_return;
        if ((0 != TIO_put_att (grp_target, dest_varid, "start_spectral_channel", NC_INT, 1, &start_pix))
            || (0 != TIO_put_att (grp_target, dest_varid, "num_spectral_channels", NC_INT, 1, &num_pix)))
          goto close_and_return;
     }

   len_slab = xtrack_dimlen_src * params_dimlen_src;

   for (i = 0; i < info.dimlens[0]; i++)
     {
        float *param_slab_i = wavecal_params + i * len_slab;

        start[0] = mirror_step[i];
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = xtrack_dimlen_src;
        count[2] = params_dimlen_src;

        if (0 != TIO_put_var_section (grp_target, params_varname, start, count,
                                      TIO_FLOAT, param_slab_i))
          goto close_and_return;

        /* FIXME: presumably drop this when we drop the explicit wavelength grid variable */
        if (0 != replace_wavelength_grid (grp_target, mirror_step[i], start_pix, num_pix,
                                          param_slab_i, xtrack_dimlen_src, params_dimlen_src,
                                          spectral_channel_dimlen, calibrated_waves))
          goto close_and_return;
     }

   status = 0;
close_and_return:
   TIO_close(ncid_src);
   FREE(wavecal_params);
   FREE(mirror_step);
   FREE(calibrated_waves);
   return status;
}

int main (int argc, char **argv)
{
   const char appname[] = "wavecal_merge";
   int status = EXIT_FAILURE;
   const char *target_file = NULL;
   const char *results_dir = NULL;
   char *pattern = NULL;
   IOCLib_Glob_Type *gt = NULL;
   int delete_files = 0;
   int ncid_target;
   size_t i, num_merged;
   static struct option long_options[] =
     {
        {"help",    no_argument, 0, 'h'},
        {"target",  required_argument, 0, 't'},
        {"delete",  no_argument, 0, 'd'},
        {0,0,0,0}
     };

   if (argc < 4)
     usage();

   tell_open (appname, -1, 0);

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "ht:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             fprintf (stderr, "getopt returned character %d??", c);
             goto return_status;
             break;
           case 'h':
             usage();
             break;
           case 'd':
             delete_files++;
             break;
           case 't':
             target_file = optarg;
             break;
          }
     }

   if (optind == argc)
     usage();

   results_dir = argv[optind];

   if (NULL == (pattern = ioclib_pathconcat (results_dir, "*.nc")))
     goto return_status;

   if (NULL == (gt = ioclib_glob (pattern, 0)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "ioclib_glob failed: %s", pattern);
        goto return_status;
     }

   if (0 != TIO_open (target_file, NC_WRITE, &ncid_target))
     goto return_status;

   num_merged = 0;

   for (i = 0; i < gt->num_files; i++)
     {
        if (0 != perform_merge (ncid_target, gt->files[i]))
          continue;

        num_merged++;
        if (delete_files)
          {
             (void) remove (gt->files[i]);
          }
     }

   if (num_merged == gt->num_files)
     {
        status = EXIT_SUCCESS;
        if (delete_files) (void) rmdir (results_dir);
     }
   else
     {
        tell_vwarn (0, "failed merging %ld of %ld wavecal parameter files",
                    gt->num_files - num_merged, gt->num_files);
     }

return_status:
   if (ncid_target)
     {
        (void) TIO_close (ncid_target);
     }
   tell_close();
   ioclib_free (pattern);
   ioclib_glob_free (gt);
   return status;
}
