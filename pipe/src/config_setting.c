#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <wordexp.h>
#include <getopt.h>
#include <libconfig.h>

#include <tell.h>

#ifndef DEFAULT_CONFIG_FILENAME
# define DEFAULT_CONFIG_FILENAME "pipe_control.cfg"
#endif

static void usage (void)
{
   fprintf (stderr, "Usage: config_setting [options] <param-path>\n");
   fprintf (stderr, "  Optional:\n");
   fprintf (stderr, "   -c | --config FILE     Configuration file\n");
   fprintf (stderr, "   -h | --help            Print this usage message\n");
   exit (EXIT_SUCCESS);
}

static char *expand_string (const char *s, int quiet)
{
   wordexp_t we = {0};
   char *s_exp = NULL;

   memset ((char *)&we, 0, sizeof (wordexp_t));

   if ((0 != wordexp (s, &we, WRDE_NOCMD | WRDE_UNDEF))
       || (we.we_wordc != 1))
     {
        if (!quiet) tell_verror (TELL_UNKNOWN_ERROR,
                                 "%s: expanding path: %s", __func__, s ? s : "(null)");
        wordfree (&we);
        return NULL;
     }

   s_exp = strdup (we.we_wordv[0]);
   wordfree (&we);

   if (NULL == s_exp)
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: strdup failed", __func__);
     }

   return s_exp;
}

int main(int argc, char **argv)
{
   const char appname[] = "config_setting";
   char *config_file = DEFAULT_CONFIG_FILENAME;
   const char *param_path;
   const char *param_str = NULL;
   int param_type;
   int status = EXIT_FAILURE;
   int malloced_config_file = 0;
   config_t cfg;
   config_setting_t *setting = NULL;
   static struct option long_options[] =
     {
        {"config",  required_argument, 0, 'c'},
	{"help",    no_argument,       0, 'h'},
        {0,0,0,0}
     };

   if (argc < 2)
     usage();

   tell_open (appname, -1, 0);

   config_init (&cfg);

   for (;;)
     {
        int option_index = 0;
        int c = getopt_long (argc, argv, "hc:", long_options, &option_index);
        if (c == -1)
          break;
        switch (c)
          {
           default:
             tell_verror (TELL_INVALID_PARM_ERROR,
                          "%s: getopt returned character %d??",
                          __func__, c);
             goto return_status;
             break;
           case 'c':
             config_file = optarg;
             break;
           case 'h':
	     usage();
             break;
          }
     }

   if (optind == argc)
     usage();

   param_path = argv[optind++];

   if (optind < argc)
     {
        fprintf (stdout, "Remaining arguments ignored:  ");
        while (optind < argc)
          {
             fprintf (stdout, "%s ", argv[optind++]);
          }
        fprintf (stdout, "\n");
     }

   if (0 != access (config_file, F_OK | R_OK))
     {
        char *config_default_paths[] =
          {
             "$SDPC_RUN_DIR_MASTER/etc/" DEFAULT_CONFIG_FILENAME,
             "$SDPC_ROOT/etc/" DEFAULT_CONFIG_FILENAME,
             NULL
          };
        const char *p;
        config_file = NULL;
        for (p = *config_default_paths; p != NULL; p++)
          {
             if (NULL == (config_file = expand_string (p, 1)))
               goto return_status;
             if (0 == access (config_file, F_OK | R_OK))
               {
                  malloced_config_file = 1;
                  break;
               }
             free (config_file);
             config_file = NULL;
          }
     }

   if(! config_read_file(&cfg, config_file))
     {
        tell_verror (TELL_RUNTIME_ERROR, "reading %s:%d - %s\n",
                     config_file, config_error_line(&cfg), config_error_text(&cfg));
        goto return_status;
     }

   if (NULL == (setting = config_lookup (&cfg, param_path)))
     {
        tell_verror (TELL_RUNTIME_ERROR, "reading %s: cannot find %s\n",
                     config_file, param_path);
        goto return_status;
     }

   param_type = config_setting_type (setting);

   switch (param_type)
     {
      case CONFIG_TYPE_STRING:
        param_str = config_setting_get_string (setting);
        (void) fprintf (stdout, "%s", param_str ? param_str : "(null)");
        break;

      case CONFIG_TYPE_BOOL:
        (void) fprintf (stdout, "%d", config_setting_get_bool (setting));
        break;

      case CONFIG_TYPE_INT:
        (void) fprintf (stdout, "%d", config_setting_get_int (setting));
        break;

      case CONFIG_TYPE_FLOAT:
        (void) fprintf (stdout, "%17.15e", config_setting_get_float(setting));
        break;

      default:
        tell_verror (TELL_NOT_IMPLEMENTED_ERROR, "no support for type %d", param_type);
        goto return_status;
        break;
     }

   status = EXIT_SUCCESS;

return_status:
   if (malloced_config_file) free (config_file);
   config_destroy(&cfg);
   tell_close();
   return status;
}
