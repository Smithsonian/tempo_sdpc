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

typedef struct
{
   float *hw1e;
   float *shape;
   float *asym;
   int num_xtrack;
   int num_channel;
}
SlitFun_Type;

static void sf_free (SlitFun_Type *sf)
{
   if (sf == NULL)
     return;
   FREE(sf->hw1e);
   FREE(sf->asym);
   FREE(sf->shape);
   return;
}

static SlitFun_Type *sf_alloc (int num_xtrack, int num_channel)
{
   SlitFun_Type *sf = NULL;
   int len = num_xtrack * num_channel;
   if (NULL == (sf = (SlitFun_Type *)MALLOC (sizeof *sf)))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc error", __func__);
        return NULL;
     }
   memset ((char *)sf, 0, sizeof(*sf));

   if ((NULL == (sf->hw1e = (float *)MALLOC (len * sizeof(float))))
       || (NULL == (sf->shape = (float *)MALLOC (len * sizeof(float))))
       || (NULL == (sf->asym = (float *)MALLOC (len * sizeof(float)))))
     {
        sf_free(sf);
        return NULL;
     }

   return sf;
}

static void usage (void)
{
   fprintf (stderr, "Usage: wavecal_merge -t <target-file> <wavecal-results-dir>\n");
   fprintf (stderr, "  Required:\n");
   fprintf (stderr, "   -t | --target FILE   Name of netCDF4 file to receive wavecal arrays\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -d | --delete        Delete input files and directory after successful merge\n");
   fprintf (stderr, "   -m | --meta          Use this option to indicate that this step\n");
   fprintf (stderr, "                        finalizes the metadata for this data product\n");
   fprintf (stderr, "   -h | --help          Print this usage message\n");
   exit (EXIT_SUCCESS);
}

// added by WHou
/* generate the index of param_slab_i from the index of opt_status_slab_i,
 * here, m = rams_dimlen_src, i is the index of param_slab_i
 */
int generate_vector_index(int *vec_ind, int i, int m)
{
   int j;
   for (j = 0; j < m; j++)
     {
        vec_ind[j] = m * i + j;
     }
   return 0;
}

/* find the neaset index of non-fail pixel for interplation,
 * here, opt_status is vector, index corresponds to the failed pixel,
 * opt_status >= 1 && <= 3 mean the sucessfull wavecal pixel.
 */
int find_nearest_non_fail_index(int *opt_status, int size, int index)
{
    int left_index  = index - 1;
    int right_index = index + 1;

    while (left_index >= 0 || right_index < size)
      {
        if ((left_index >= 0) && (opt_status[left_index] >= 1)
                              && (opt_status[left_index] <= 3))
           return left_index;

        if ((right_index < size) && (opt_status[right_index] >= 1)
                                 && (opt_status[right_index] <= 3))
           return right_index;

        left_index--;
        right_index++;
      }
    // if no failed opt_status, then retun -1
    return -1;
}
// end adding

static int perform_merge (int ncid_target, const char *file)
{
   const char *params_varname = TEMPO_VAR_WAVECAL_PARAM;
   TIO_Var_Info_Type info = {0};
   TIO_Var_Info_Type info1 = {0}; // added by WHou
   char group_name[TIO_MAX_NAME_LEN] = {0};
   int ncid_src, grp_target, varid, start_pix, num_pix, num_coefs;
   int step_dimlen_src, xtrack_dimlen_src;
   int step_dimid, xtrack_dimid, channel_dimid, dest_varid;
   int have_sf, start0, count0, xtrack0;
   int start[3], count[3];
   int have_adjust_att, adjust_att;
   size_t step_dimlen, xtrack_dimlen, channel_dimlen;
   size_t params_dimlen_src, len_params, len_slab;
   size_t channel_dimlen_src;
   size_t i;
   size_t len_niter; // added by WHou
   int *niter = NULL; // added by WHou
   int *opt_status = NULL; // added by WHou
   int *vec_ind_j0 = NULL; // added by WHou
   int *vec_ind_j1 = NULL; // added by WHou
   int *vec_opt_status_flag = NULL; // added by WHou
   float *wavecal_params = NULL;
   SlitFun_Type *sf = NULL;
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

   if (0 == tio_inq_varid (ncid_src, "sf_hw1e", &have_sf))
     {
        TIO_Var_Info_Type sf_info = {0};
        if (0 != TIO_inq_var (ncid_src, "sf_hw1e", &sf_info))
          goto close_and_return;
        channel_dimlen_src = sf_info.dimlens[2];
     }
   else
     {
        have_sf = 0;
        channel_dimlen_src = 0;
     }

   // added by WHou
   if ((0 != TIO_inq_var (ncid_src, "niter", &info1))
       || (0 != TIO_inq_var (ncid_src, "opt_status", &info1)))
     goto close_and_return;
   // end adding

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

   if (0 != TIO_inq_dim (grp_target, TEMPO_DIM_CHANNEL, &channel_dimid, &channel_dimlen))
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: reading xtrack dimension size", __func__);
        goto close_and_return;
     }

   if (step_dimlen != (size_t) step_dimlen_src)
     {
        tell_vwarn (0, "dimension mismatch (mirror_step) source:%d target:%ld",
                    step_dimlen_src, step_dimlen);
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

   // added by WHou
   /* same length for niter & opt_status */
   len_niter = step_dimlen_src * xtrack_dimlen_src;
   if ((NULL == (niter = MALLOC (len_niter * sizeof(int))))
       ||(NULL == (opt_status = MALLOC (len_niter * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }
   // end adding

   if (have_sf)
     {
        if (NULL == (sf = sf_alloc (xtrack_dimlen_src, channel_dimlen_src)))
          goto close_and_return;
     }

   start[0] = 0;
   start[1] = 0;
   start[2] = 0;
   count[0] = info.dimlens[0];
   count[1] = info.dimlens[1];
   count[2] = params_dimlen_src;

   if ((0 != TIO_get_var_section (ncid_src, TEMPO_DIM_STEP, start, count,
                                  TIO_INT, mirror_step))
       ||(0 != TIO_get_var_section (ncid_src, params_varname, start, count,
                                    TIO_FLOAT, wavecal_params)))
     goto close_and_return;

   start0 = 0;
   count0 = 1;
   if (0 != TIO_get_var_section (ncid_src, TEMPO_DIM_XTRACK, &start0, &count0,
                                 TIO_INT, &xtrack0))
     goto close_and_return;

   if ((0 != tio_inq_varid (ncid_src, params_varname, &varid))
       || (0 != TIO_get_att (ncid_src, varid, "num_coefficients", NC_INT, &num_coefs))
       || (0 != TIO_get_att (ncid_src, varid, "start_spectral_channel", NC_INT, &start_pix))
       || (0 != TIO_get_att (ncid_src, varid, "num_spectral_channels", NC_INT, &num_pix)))
     goto close_and_return;

   // added by WHou
   count[2] = 0;

   if ((0 != TIO_get_var_section (ncid_src, "niter", start, count,
                                  TIO_INT, niter))
       || (0 != TIO_get_var_section (ncid_src, "opt_status", start, count,
                                  TIO_INT, opt_status)))
     goto close_and_return;
   // end anding

   /* For back-compatibility, this attribute is optional */
   adjust_att = 0;
   have_adjust_att = 0;
   if (NC_NOERR == nc_get_att_int (ncid_src, varid, "adjust_nominal_wavelength", &adjust_att))
     have_adjust_att = 1;

   if (have_sf)
     {
        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = info.dimlens[0];
        count[1] = info.dimlens[1];
        count[2] = channel_dimlen_src;

        if ((0 != TIO_get_var_section (ncid_src, "sf_hw1e", start, count, TIO_FLOAT, sf->hw1e))
            ||(0 != TIO_get_var_section (ncid_src, "sf_shape", start, count, TIO_FLOAT, sf->shape))
            ||(0 != TIO_get_var_section (ncid_src, "sf_asym", start, count, TIO_FLOAT, sf->asym)))
          goto close_and_return;
     }

   /* write params, creating target variable if necessary */

   if (0 != tio_inq_varid (grp_target, params_varname, &dest_varid))
     {
        const char *params_dimname = "wavecal_par";
        size_t params_dimlen;
        int params_dimid, params_dimid_list[3];
        if (0 != TIO_inq_dim (grp_target, params_dimname, &params_dimid, &params_dimlen))
          {
             int unused_varid;
             if (0 != TIO_def_dim (grp_target, params_dimname, params_dimlen_src, &params_dimid))
               goto close_and_return;
             /* Apparently, python-3.12 netcdf4 needs a variable with this name to exist.
              * If the variable doesn't exist, then attempting to read the resulting file
              * fails with this error:
              *   AttributeError: 'NoneType' object has no attribute 'dimensions'
              * Therefore, don't delete this unused variable unless you're sure that this
              * python issue is no longer relevant.
              */
             if (0 != TIO_def_var (grp_target, params_dimname, TIO_INT, 1, &params_dimid, &unused_varid))
               goto close_and_return;
          }
        params_dimid_list[0] = step_dimid;
        params_dimid_list[1] = xtrack_dimid;
        params_dimid_list[2] = params_dimid;

        if (0 != TIO_def_var (grp_target, params_varname, TIO_FLOAT, 3, params_dimid_list, &dest_varid))
          goto close_and_return;

        if ((0 != TIO_put_att (grp_target, dest_varid, "num_coefficients", NC_INT, 1, &num_coefs))
            || (0 != TIO_put_att (grp_target, dest_varid, "start_spectral_channel", NC_INT, 1, &start_pix))
            || (0 != TIO_put_att (grp_target, dest_varid, "num_spectral_channels", NC_INT, 1, &num_pix)))
          goto close_and_return;
        if (have_adjust_att)
          {
             if ((0 != TIO_put_att (grp_target, dest_varid, "adjust_nominal_wavelength", NC_INT, 1, &adjust_att)))
               goto close_and_return;
          }
     }

   if ((have_sf != 0)
       && (0 != tio_inq_varid (grp_target, "sf_hw1e", &dest_varid)))
     {
        int sf_dimids[3];
        sf_dimids[0] = step_dimid;
        sf_dimids[1] = xtrack_dimid;
        sf_dimids[2] = channel_dimid;
        if (0 != TIO_def_var (grp_target, "sf_hw1e", TIO_FLOAT, 3, sf_dimids, &dest_varid))
          goto close_and_return;
        if (0 != TIO_def_var (grp_target, "sf_shape", TIO_FLOAT, 3, sf_dimids, &dest_varid))
          goto close_and_return;
        if (0 != TIO_def_var (grp_target, "sf_asym", TIO_FLOAT, 3, sf_dimids, &dest_varid))
          goto close_and_return;
     }

   // added by WHou
   if ((0 != tio_inq_varid (grp_target, "wavecal_niter", &dest_varid))
       && (0 != tio_inq_varid (grp_target, "wavecal_opt_status", &dest_varid)))
     {
        int niter_dimids[2];

        niter_dimids[0] = step_dimid;
        niter_dimids[1] = xtrack_dimid;

        if ((0 != TIO_def_var (grp_target, "wavecal_niter", TIO_INT, 2, niter_dimids, &dest_varid))
            || (0 != TIO_def_var (grp_target, "wavecal_opt_status", TIO_INT, 2, niter_dimids, &dest_varid)))
          goto close_and_return;
     }

   /* record the flag of refilled pixel for failure (if opt_status > 3 or < 1) */
   if ((NULL == (vec_opt_status_flag = (int *)MALLOC(info.dimlens[0] * xtrack_dimlen_src * sizeof(int))))
       || (NULL == (vec_ind_j0 = (int *)MALLOC (params_dimlen_src * sizeof(int))))
       || (NULL == (vec_ind_j1 = (int *)MALLOC (params_dimlen_src * sizeof(int)))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        goto close_and_return;
     }

   /* assign the initial value to 0 */
   for (i = 0; i < info.dimlens[0] * xtrack_dimlen_src; i++)
     {
        vec_opt_status_flag[i] = 0;
     }
   // end adding

   len_slab = xtrack_dimlen_src * params_dimlen_src;

   for (i = 0; i < info.dimlens[0]; i++)
     {
        float *param_slab_i = wavecal_params + i * len_slab;

        // added by WHou
        /* niter = -1 means without RAD_wavecal due to the setting of sza_max */
        int *niter_slab_i = niter + i * xtrack_dimlen_src;
        int *opt_status_slab_i = opt_status + i * xtrack_dimlen_src;

        for (int j = 0; j < xtrack_dimlen_src; j++)
          {
             /* if opt_status > 3 or < 1, refilling is required
              * if niter = -1, no refilling is required
              */
             if (((opt_status_slab_i[j] > 3) || (opt_status_slab_i[j] < 1)) && (niter_slab_i[j] > 0))
               {
                  int j1 = find_nearest_non_fail_index(opt_status_slab_i, xtrack_dimlen_src, j);
                  size_t k;

                  if (j1 >= 0)
                    {
                       /* generate vector index of param_slab_i */
                       (void) generate_vector_index(vec_ind_j0, j,  params_dimlen_src);
                       (void) generate_vector_index(vec_ind_j1, j1, params_dimlen_src);

                       for (k = 0; k < params_dimlen_src; k++)
                         {
                            int ind_k0 = vec_ind_j0[k];
                            int ind_k1 = vec_ind_j1[k];

                            /* refill the element's value of param_slab_i
                             * with the generate vector index
                             */
                            param_slab_i[ind_k0] = param_slab_i[ind_k1];
                          }
                        /* record the flag of refill with 1*/
                        vec_opt_status_flag[i * xtrack_dimlen_src + j] = 1;
                    }
               }
          }
        // end adding

        start[0] = mirror_step[i];
        start[1] = xtrack0;
        start[2] = 0;
        count[0] = 1;
        count[1] = xtrack_dimlen_src;
        count[2] = params_dimlen_src;

        if (0 != TIO_put_var_section (grp_target, params_varname, start, count,
                                      TIO_FLOAT, param_slab_i))
          goto close_and_return;

        if (have_sf)
          {
             start[0] = mirror_step[i];
             start[1] = xtrack0;
             start[2] = 0;
             count[0] = 1;
             count[1] = xtrack_dimlen_src;
             count[2] = channel_dimlen_src;
             if ((0 != TIO_put_var_section (grp_target, "sf_hw1e", start, count, TIO_FLOAT, sf->hw1e))
                 || (0 != TIO_put_var_section (grp_target, "sf_shape", start, count, TIO_FLOAT, sf->shape))
                 || (0 != TIO_put_var_section (grp_target, "sf_asym", start, count, TIO_FLOAT, sf->asym)))
               goto close_and_return;
          }
     }

   // added by WHou
   for (i = 0; i < info.dimlens[0]; i++)
     {
        int *niter_slab_i = niter + i * xtrack_dimlen_src;
        int *opt_status_slab_i = opt_status + i * xtrack_dimlen_src;

        /* update the value of opt_status for the refilled pixel */
        for (int j = 0; j < xtrack_dimlen_src; j++)
          {
             if (vec_opt_status_flag[i * xtrack_dimlen_src + j] == 1)
               {
                  /* update the refilled pixel with opt_status = 9 */
                  opt_status_slab_i[j] = 9;
               }
          }

        start[0] = mirror_step[i];
        start[1] = xtrack0;
        start[2] = 0;
        count[0] = 1;
        count[1] = xtrack_dimlen_src;
        count[2] = 0;

        if ((0 != TIO_put_var_section (grp_target, "wavecal_niter", start, count,
                                      TIO_INT, niter_slab_i))
            || (0 != TIO_put_var_section (grp_target, "wavecal_opt_status", start, count,
                                      TIO_INT, opt_status_slab_i)))
          goto close_and_return;
     }
   // end adding

   status = 0;
close_and_return:
   TIO_close(ncid_src);
   FREE(wavecal_params);
   sf_free (sf);
   FREE(mirror_step);
   FREE(niter); // added by WHou
   FREE(opt_status); // added by WHou
   FREE(vec_ind_j0); // added by WHou
   FREE(vec_ind_j1); // added by WHou
   FREE(vec_opt_status_flag);  // added by WHou

   return status;
}

static int read_metadata (TIO_Meta_Type *meta, const char *file)
{
   int ncid, status;

   if (0 != TIO_open (file, NC_NOWRITE, &ncid))
     return -1;

   status = tio_meta_ncinit (meta, ncid, "input_files", TIO_META_TYPE_STRING);

   (void) TIO_close (ncid);

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
   TIO_Meta_Type *meta = NULL;
   int delete_files = 0;
   int ncid_target;
   int finalize_metadata = 0;
   size_t i, num_merged;
   static struct option long_options[] =
     {
        {"help",    no_argument, 0, 'h'},
        {"target",  required_argument, 0, 't'},
        {"delete",  no_argument, 0, 'd'},
        {"meta",    no_argument, 0, 'm'},
        {0,0,0,0}
     };

   if (argc < 4)
     usage();

   tell_open (appname, -1, 0);

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hmt:", long_options, &option_index);
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
           case 'm':
             finalize_metadata++;
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

   tio_set_cmdline (argc, argv);

   if (NULL == (pattern = ioclib_pathconcat (results_dir, "*.nc")))
     goto return_status;

   if (NULL == (gt = ioclib_glob (pattern, 0)))
     {
        tell_verror (TELL_APPLICATION_ERROR, "%s: ioclib_glob failed: %s",
                     __func__, pattern);
        goto return_status;
     }

   if (gt->num_files == 0)
     {
        tell_verror (TELL_RUNTIME_ERROR, "%s: no files match glob pattern: %s",
                     __func__, pattern);
        goto return_status;
     }

   if (NULL == (meta = tio_meta_open ()))
     goto return_status;

   /* Initialize selected metadata entries from target file (to preserve that),
    * and 1st merge block (to append whatever additional metadata it may contain).
    * Assume the additional metadata is the same for all the merge blocks.
    */
   (void) read_metadata (meta, target_file);
   (void) read_metadata (meta, gt->files[0]);

   if (0 != TIO_open (target_file, NC_WRITE, &ncid_target))
     goto return_status;

   if (0 != tio_history_append_cmdline (ncid_target))
     goto return_status;

   /* Write the accumulated metadata entries */
   if (0 != tio_meta_write_ncattr (meta, ncid_target))
     goto return_status;

   if (finalize_metadata == 0)
     {
        tio_meta_set_noexpand (meta, "input_files", 1);
     }

   /* If no template exists, a warning will be printed,
    * but no error will be generated
    */
   if (0 != tio_meta_expand_file (meta, NULL, target_file))
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
   tio_meta_close (meta);
   ioclib_free (pattern);
   ioclib_glob_free (gt);
   return status;
}
