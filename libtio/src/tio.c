/** @file
 *  @brief Core C functions
 */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <float.h>

#include <cfortran.h>
#include <netcdf.h>
#include <tell.h>

#include "tio.h"
#include "tio_template.h"
#include "_tio.h"
#include "tio_ppc.h"

#define EMPTY()
#define _pTIO_STR(s) #s

#define TIMESTAMP_BUFLEN 26

/*{{{ I/O Tracing facility */

static int _pTIO_TRACING = -1;
static char _pTIO_TRACE_PREFIX[128];

#define TRACE_PATH_BUFSIZE 1024

/* This is linux-specific, but for the current purpose that's ok */
static void program_basename (char *name, size_t namesize)
{
   const char proc_self[] = "/proc/self/exe";
   char buf[1024];
   char *p;
   ssize_t len;

   if ((0 != access (proc_self, F_OK | R_OK))
       || ((len = readlink (proc_self, buf, sizeof(buf))) < 0))
     {
        name[0] = 0;
        return;
     }

   buf[len] = 0;

   if (NULL != (p = strrchr (buf, '/')))
     {
        p++;
     }
   else p = buf;

   (void) snprintf (name, namesize, "%s", p);
}

static void trace_make_prefix (void)
{
   char pgm_basename[64];
   program_basename (pgm_basename, sizeof(pgm_basename));
   if (pgm_basename[0] != 0)
     snprintf (_pTIO_TRACE_PREFIX, sizeof (_pTIO_TRACE_PREFIX), "TIO_TRACE|%d:%s", getpid(), pgm_basename);
   else
     snprintf (_pTIO_TRACE_PREFIX, sizeof (_pTIO_TRACE_PREFIX), "TIO_TRACE|%d", getpid());
}

static void trace_filename (int ncid, char *file, size_t len)
{
   unsigned int n = 2;
   while (n-- > 0)
     {
        int grp;
        if (NC_NOERR == nc_inq_path (ncid, NULL, file))
          return;
        grp = ncid;
        (void) nc_inq_grp_parent (grp, &ncid);
     }
   (void) snprintf (file, len, "unknown");
}

static void trace_init (void)
{
   if (_pTIO_TRACING == -1)
     {
        /* Do this only once for each process id */
        if (NULL == getenv ("TIO_ENABLE_TRACING"))
          {
             _pTIO_TRACING = 1; /* off */
          }
        else
          {
             _pTIO_TRACING = 0; /* on */
             trace_make_prefix();
          }
     }
}

static void trace_print (const char *action, int parent_ncid, const char *file,
                         int grp, const char *name, size_t size, const char *rest)
{
   (void) fprintf (stderr, "%s|%s|%d|%s|%d|%s|%ld|%s",
                   _pTIO_TRACE_PREFIX, action, parent_ncid,
                   file ? file : "null", grp,
                   name ? name : "", size,
                   rest ? rest : "\n");
}

void _pTIO_trace_close (int ncid)
{
   if (_pTIO_TRACING) return;
   trace_print ("CLOSED", ncid, "", 0, NULL, 0, NULL);
}

void _pTIO_trace_open (int ncid, const char *file)
{
   trace_init();
   if (_pTIO_TRACING) return;
   trace_print ("OPENED", ncid, file, 0, NULL, 0, NULL);
}

void _pTIO_trace_create (int ncid, const char *file)
{
   trace_init();
   if (_pTIO_TRACING) return;
   trace_print ("CREATED", ncid, file, 0, NULL, 0, NULL);
}

void _pTIO_trace_group (int parent_ncid, const char *path, int grp)
{
   char file[TRACE_PATH_BUFSIZE];
   if (_pTIO_TRACING) return;
   trace_filename (parent_ncid, file, sizeof(file));
   trace_print ("GROUP", parent_ncid, file, grp, path, 0, NULL);
}

static int trace_var_io (const char *action, int grp, const char *name, int dim, const size_t *start, const size_t *count)
{
   char file[TRACE_PATH_BUFSIZE];
   char grppath[128];
   char varpath[256];
   char rest[128];
   size_t size;
   int i, parent_ncid;
   if (_pTIO_TRACING) return 0;
   if ((NC_NOERR == nc_inq_grpname_full (grp, NULL, grppath))
       && (strlen(grppath) > 1))
     {
        (void) snprintf (varpath, sizeof(varpath), "%s/%s", grppath, name);
     }
   else strncpy (varpath, name, sizeof(varpath)-1);
   if (NC_NOERR != nc_inq_grp_parent (grp, &parent_ncid))
     parent_ncid = grp;
   trace_filename (parent_ncid, file, sizeof(file));
   size = 1;
   for (i = 0 ; i < dim; i++)
     {
        size *= count[i];
     }
   if (dim > 0)
     {
        size_t len = sizeof(rest);
        char *p = rest;
        int n = snprintf (p, len, "dim=%d,start=[%ld", dim, start[0]);
        p   += n;
        len -= n;
        for (i = 1; i < dim; i++)
          {
             n = snprintf (p, len, ",%ld", start[i]);
             p   += n;
             len -= n;
          }
        n = snprintf (p, len, "],count=[%ld", count[0]);
        p   += n;
        len -= n;
        for (i = 1; i < dim; i++)
          {
             n = snprintf (p, len, ",%ld", count[i]);
             p   += n;
             len -= n;
          }
        n = snprintf (p, len, "]\n");
        p   += n;
        len -= n;
     }
   else strcpy (rest, "\n");

   trace_print (action, parent_ncid, file, grp, varpath, size, rest);

   return 0;
}

static void trace_get_var (int grp, const char *name, int dim, const size_t *start, const size_t *count)
{
   (void) trace_var_io ("READ", grp, name, dim, start, count);
}

static void trace_put_var (int grp, const char *name, int dim, const size_t *start, const size_t *count)
{
   (void) trace_var_io ("WROTE", grp, name, dim, start, count);
}

/*}}}*/

int _pTIOMake_Name_UInt_Arrays (_pName_UInt_Pair_Type *array,
                                int *pnum_values, char **pnames,
                                unsigned int **pvalues)
{
   _pName_UInt_Pair_Type *p;
   int total_len, num, i;
   unsigned int *values = NULL;
   char *names = NULL;
   char *s, *end;

   num = total_len = 0;
   p = array;

   while (1)
     {
        if (p->name[0] == 0)
          break;
        total_len += strlen(p->name);
        num += 1;
        p++;
     }

   /* add room for num-1 spaces, plus 1 trailing null character */
   total_len += num;

   if ((NULL == (values = (unsigned int *) TIO_MALLOC (num * sizeof(unsigned int))))
       || (NULL == (names = (char *) TIO_MALLOC (total_len * sizeof(char)))))
     {
        Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        TIO_FREE(values);
        TIO_FREE(names);
        return -1;
     }

   s = names;
   end = names + total_len;

   for (i = 0; i < num; i++)
     {
        int len;
        p = &array[i];

        values[i] = p->value;

        len = strlen(p->name);
        strcpy (s, p->name);
        s += len;

        if (s+1 < end)
          {
             *s = ' ';
             s++;
          }
        else *s = 0;
     }

   *pnum_values = num;
   *pnames = names;
   *pvalues = values;

   return 0;
}

int _pTIO_define_dims_using_offsets (int grp, const _pDim_Offsets_Type *offsets,
                                     _pDim_Table_Type *dim_table)
{
   const _pDim_Offsets_Type *o;
   char *p = (char *)dim_table;
   int status;

   for (o = offsets; o->name != NULL; o++)
     {
        status = nc_def_dim (grp, o->name, *(size_t *)(p + o->len_offset),
                             (int *)(p + o->id_offset));
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dimension %s (%s)",
                          __func__, o->name, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

int TIO_define_enum_table (int grp, const char *name, int base_type,
                           const TIO_Enum_Type *enum_table, int *enum_typeid)
{
   const TIO_Enum_Type *e;
   int status;

   if (NC_NOERR != (status = nc_def_enum (grp, base_type, name, enum_typeid)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining enum %s (%s)",
                     __func__, name, nc_strerror(status));
        return -1;
     }

   for (e = enum_table; e->name != NULL; e++)
     {
        status = nc_insert_enum (grp, *enum_typeid, e->name, &e->value);
        if (status != NC_NOERR)
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: inserting value %s=%d for enum %s (%s)",
                          __func__, e->name ? e->name : "(null)",
                          e->value, name, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_short_attrs (int grp, int varid, const _pShort_Attr_Type *attrs)
{
   const _pShort_Attr_Type *a;
   int status;

   for (a = attrs; a->name != NULL; a++)
     {
        status = nc_put_att_short (grp, varid, a->name, NC_SHORT, 1, &a->value);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: defining short attribute %s (%s)",
                          __func__, a->name, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_int_attrs (int grp, int varid, const _pInt_Attr_Type *attrs)
{
   const _pInt_Attr_Type *a;
   int status;

   for (a = attrs; a->name != NULL; a++)
     {
        status = nc_put_att_int (grp, varid, a->name, NC_INT, a->num_values, a->value);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: defining int attribute %s (%s)",
                          __func__, a->name, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_float_attrs (int grp, int varid, const _pFloat_Attr_Type *attrs)
{
   const _pFloat_Attr_Type *a;
   int status;

   for (a = attrs; a->name != NULL; a++)
     {
        status = nc_put_att_float (grp, varid, a->name, NC_FLOAT, 1, &a->value);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: defining float attribute %s (%s)",
                          __func__, a->name, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_text_attrs (int grp, int varid, const _pText_Attr_Type *attrs)
{
   const _pText_Attr_Type *a;
   int status;

   for (a = attrs; a->name != NULL; a++)
     {
        size_t len = strlen(a->text) + 1;
        status = nc_put_att_text (grp, varid, a->name, len, a->text);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_WRITE_ERROR,
                          "%s: defining text attribute %s (%s)",
                          __func__, a->name, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

int _pTIO_define_var_with_text_attrs (int grp, const char *var_name, nc_type xtype,
                                      int num_dims, const int *dimids,
                                      const _pText_Attr_Type *text_attrs,
                                      int *pvarid)
{
   int status, varid;

   status = nc_def_var (grp, var_name, xtype, num_dims, dimids, &varid);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining variable %s (%s)",
                     __func__, var_name, nc_strerror(status));
        return -1;
     }

   if (text_attrs != NULL)
     {
        if (-1 == _pTIO_define_text_attrs (grp, varid, text_attrs))
          return -1;
     }

   if (pvarid != NULL)
     {
        *pvarid = varid;
     }

   return 0;
}

int _pTIO_put_fillvalue_attr (int grp, int varid, nc_type xtype)
{
   int status;
   void *pfill_value;
   char fill_char = TIO_FILL_CHAR;
   short fill_short = TIO_FILL_SHORT;
   unsigned short fill_ushort = TIO_FILL_USHORT;
   int fill_int = TIO_FILL_INT;
   unsigned int fill_uint = (unsigned int) -1;
   float fill_float = TIO_FILL_FLOAT;
   double fill_double = TIO_FILL_DOUBLE;

   switch (xtype)
     {
      case NC_CHAR: pfill_value = &fill_char;
        break;
      case NC_SHORT: pfill_value = &fill_short;
        break;
      case NC_USHORT: pfill_value = &fill_ushort;
        break;
      case NC_INT: pfill_value = &fill_int;
        break;
      case NC_UINT: pfill_value = &fill_uint;
        break;
      case NC_FLOAT: pfill_value = &fill_float;
        break;
      case NC_DOUBLE: pfill_value = &fill_double;
        break;
      default:
        Tell_verror (TELL_INVALID_PARM,
                     "%s: invalid fill value type xtype=%d", __func__, xtype);
        return -1;
     }

   status = nc_put_att (grp, varid, _FillValue, xtype, 1, pfill_value);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "writing fill value to grp=%d varid=%d (%s)",
                     grp, varid, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_create (const char *path, int cmode, int *ncid)
{
   int status;

   if ((path == NULL) || (ncid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_create (path, cmode, ncid)))
     {
        Tell_verror (TELL_IO_OPEN_ERROR, "%s: creating %s (%s)",
                     __func__, path, nc_strerror(status));
        return -1;
     }
   _pTIO_trace_create (*ncid, path);

   return 0;
}

int TIO_open (const char *path, int omode, int *ncid)
{
   int status;

   if ((path == NULL) || (ncid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_open (path, omode, ncid)))
     {
        Tell_verror (TELL_IO_OPEN_ERROR, "%s: opening %s (%s)",
                     __func__, path, nc_strerror(status));
        return -1;
     }
   _pTIO_trace_open (*ncid, path);

   _pTIO_warn_about_time_reference_mismatch (*ncid);

   return 0;
}

int TIO_close (int ncid)
{
   int status;

   if (NC_NOERR != (status = nc_close (ncid)))
     {
        Tell_verror (TELL_IO_ERROR, "%s: closing file %d (%s)",
                     __func__, ncid, nc_strerror(status));
        return -1;
     }
   _pTIO_trace_close (ncid);

   return 0;
}

int tio_sync (int ncid)
{
   int status;

   if (NC_NOERR != (status = nc_sync (ncid)))
     {
        tell_verror (TELL_IO_ERROR, "%s: syncing file %d (%s)",
                     __func__, ncid, nc_strerror(status));
        return -1;
     }

   return 0;
}

int tio_append_history (int ncid, const char *str)
{
   const char *history_attr = "history";
   char buf[TIMESTAMP_BUFLEN];
   time_t now;
   struct tm tm = {0};
   int status, nctype;
   char *att = NULL;
   size_t len, hlen, entry_len;
   int have_attribute = 0;
   int offset = 0;

   if ((str == NULL) || (*str == 0))
     return 0;

   if (NULL != getenv ("TIO_HISTORY_OFF"))
     return 0;

   now = time(NULL);
   gmtime_r (&now, &tm);
   strftime (buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);

   /* <time>:<str>\n */
   entry_len = strlen(buf) + strlen(str) + 3;
   len = entry_len;

   status = nc_inq_att (ncid, NC_GLOBAL, history_attr, &nctype, &hlen);
   if (status == NC_NOERR)
     {
        if (nctype != NC_CHAR)
          {
             tell_verror (TELL_UNSUPPORTED_ERROR, "%s: attribute %s is not of type NC_CHAR",
                          __func__, history_attr);
             return -1;
          }
        have_attribute++;
        /* plus one in case we need to add \n at the end of the existing attribute */
        len += hlen + 1;
     }

   if (NULL == (att = (char *)TIO_MALLOC (len * sizeof(char))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }
   memset ((char *)att, 0, len * sizeof(char));

   if (have_attribute)
     {
        if (NC_NOERR != (status = nc_get_att_text (ncid, NC_GLOBAL, history_attr, att)))
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: reading history attribute (%s)",
                          __func__, nc_strerror (status));
             goto return_status;
          }
        offset = strlen(att);
        /* ensure that the existing comment ends with \n */
        if (att[offset-1] != '\n')
          {
             att[offset++] = '\n';
          }
     }

   /* history entry ends with \n */
   sprintf (att + offset, "%s:%s\n", buf, str);
   len = strlen (att) + 1;

   if (NC_NOERR != (status = nc_put_att (ncid, NC_GLOBAL, history_attr, NC_CHAR, len, att)))
     {
        tell_verror (TELL_IO_WRITE_ERROR, "%s: updating history attribute (%s)",
                     __func__, nc_strerror (status));
        goto return_status;
     }

   status = 0;
return_status:
   TIO_FREE(att);
   return status ? -1 : 0;
}

char *tio_concat_argv (int argc, char **argv, char *pstr, size_t len_pstr)
{
   char *str = NULL;
   size_t len, offset;
   int i;

   if (argv == NULL)
     {
        tell_verror (TELL_INVALID_PARM_ERROR, "%s: argv != NULL is required", __func__);
        return NULL;
     }

   /* allocate storage to concatenate argc tokens, with a space delimiter,
    * and a terminating null character */

   len = 1;
   for (i = 0; i < argc; i++)
     {
        if (argv[i] == NULL)
          continue;
        len += strlen(argv[i]) + 1;
     }

   if (pstr != NULL)
     {
        if (len_pstr < len)
          {
             tell_verror (TELL_RUNTIME_ERROR, "%s: command line too long (len=%ld, buffer=%ld)",
                          __func__, len, len_pstr);
             return NULL;
          }
        str = pstr;
     }
   else
     {
        if (NULL == (str = (char *)TIO_MALLOC (len * sizeof (char))))
          {
             tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             return NULL;
          }
     }

   /* at this point, str points to at least 'len' characters of storage,
    * which is guaranteed to be enough to hold the command line.
    */

   offset = 0;
   for (i = 0; i < argc; i++)
     {
        if (argv[i])
          {
             sprintf (str + offset, " %s", argv[i]);
             offset += strlen(argv[i]) + 1;
          }
     }

   return str;
}

static int _pTIO_Argc;
static char **_pTIO_Argv;

int tio_set_cmdline (int argc, char **argv)
{
   _pTIO_Argc = argc;
   _pTIO_Argv = argv;
   return 0;
}

int tio_history_append_cmdline (int ncid)
{
   char *str;
   int status;

   if (NULL == (str = tio_concat_argv (_pTIO_Argc, _pTIO_Argv, NULL, 0)))
     return -1;
   status = tio_append_history (ncid, str);
   TIO_FREE(str);
   return status;
}

int tio_inq_varid (int grp, const char *name, int *varid)
{
   return nc_inq_varid (grp, name, varid);
}

int TIO_inq_var (int grp, const char *name, TIO_Var_Info_Type *info)
{
   int status, i;

   if ((name == NULL) || (info == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_varid (grp, name, &info->varid)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing variable %s in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_vartype (grp, info->varid, &info->type)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing variable %s in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_varndims (grp, info->varid, &info->ndims)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing variable %s in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        return -1;
     }

   /* convenient default for the ndims=0 case */
   info->dimlens[0] = 1;

   if (info->ndims == 0)
     return 0;

   if (NC_NOERR != (status = nc_inq_vardimid (grp, info->varid, info->dimids)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing variable %s dimids in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        return -1;
     }
   for (i = 0; i < info->ndims; i++)
     {
        status = nc_inq_dimlen  (grp, info->dimids[i], &info->dimlens[i]);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: accessing dimension %d length in group %d (%s)",
                          __func__, info->dimids[i], grp, nc_strerror(status));
             return -1;
          }
     }

   return 0;
}

#if 0
static int cvt_float_to_type (int num, const float *f, int type, void *v)
{
   double *dbl = (double *)v;
   int i;

   if (type != NC_DOUBLE)
     {
        Tell_verror (TELL_USAGE_ERROR,
                     "%s: Conversion from float to type=%d is not supported",
                     __func__, type);
        return -1;
     }

   for (i = 0; i < num; i++)
     {
        dbl[i] = (double) f[i];
     }

   return 0;
}
#endif

static int cvt_type_to_float (int num, float *f, int type, const void *v)
{
   double *dbl = (double *)v;
   int i;

   if (type != NC_DOUBLE)
     {
        Tell_verror (TELL_USAGE_ERROR,
                     "%s: Conversion from type=%d to float is not supported",
                     __func__, type);
        return -1;
     }

   for (i = 0; i < num; i++)
     {
        f[i] = (float) dbl[i];
     }

   return 0;
}

static int put_float_nsd (int grp, const int *istart, const int *icount,
                          int type, const void *data, const char *name,
                          void *client_data)
{
   TIO_Var_Info_Type info;
   size_t start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS];
   const char nsd_att_name[] = "number_of_significant_digits";
   float *float_error = NULL;
   float missing_value = FLT_MAX;
   int status, i, num, nsd;
   int malloced_temp_float = 0;
   int return_status = -1;

   if (client_data == NULL)
     {
        Tell_verror (TELL_INTERNAL_ERROR, "%s: NULL client_data pointer",
                     __func__);
        return -1;
     }
   nsd = *(int *)client_data;

   if (-1 == TIO_inq_var (grp, name, &info))
     return -1;

   for (i = 0; i < info.ndims; i++)
     {
        start[i] = istart[i];
        count[i] = (icount[i] >= 0) ? (size_t) icount[i] : info.dimlens[i];
     }
   num = count[0];
   for (i = 1; i < info.ndims; i++)
     {
        num *= count[i];
     }

   if (type == NC_FLOAT)
     float_error = (float *) data;
   else
     {
        /* (optionally) convert 'type' to float */
        if (NULL == (float_error = (float *) TIO_MALLOC (num * sizeof(*float_error))))
          {
             Tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
             goto free_and_return;
          }
        malloced_temp_float = 1;
        if (-1 == cvt_type_to_float (num, float_error, type, data))
          goto free_and_return;
     }

   if (-1 == TIO_inq_var_fill (grp, info.varid, NULL, (void *)&missing_value))
     goto free_and_return;

   /* Replace the least-significant bits with either zeros or ones.
    * The modified values will compress better. */
   if (0 != _pTIO_ppc_f32_bitmask (nsd, _pTIO_PPC_METHOD_ALT, num, float_error, &missing_value))
     goto free_and_return;

   if (NC_NOERR != (status = nc_put_vara_float (grp, info.varid, start, count, float_error)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: accessing variable %s in group %d (%s)",
                     __func__, name, grp, nc_strerror(status));
        goto free_and_return;
     }
   if (NC_NOERR != (status = nc_put_att_int (grp, info.varid, nsd_att_name, NC_INT, 1, &nsd)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: writing attribute %s of %s in group %d (%s)",
                     __func__, nsd_att_name, name, grp, nc_strerror(status));
        goto free_and_return;
     }

   trace_put_var (grp, name, info.ndims, start, count);

   return_status = 0;
free_and_return:
   if (malloced_temp_float) TIO_FREE(float_error);
   return return_status;
}

typedef struct
{
   double a;
   double b;
   double *coef;
   int num_coef;
}
Cheb_Type;

static double cheb_expansion_eval (const Cheb_Type *s, double x)
{
   size_t i, order = s->num_coef-1;
   double d1, d2, t, t2, value;

   /* Use Clenshaw recursion to evaluate a truncated Chebyshev series,
    *       f(x) = \sum c(k) T(t;k), k = 0,1,2...order
    * at a particular coordinate x, in the interval [a,b], via the Chebyshev
    * coordinate t, on the interval [-1,1], is: t = (2*x-a-b)/(b-a).
    * In this expansion, T(t;k) is a Chebyshev polynomial of the first kind
    * of order k, and the c(k) are constant coefficients.
    */

   /* Clenshaw recursion: original reference is
    * Clenshaw, C.W., Math. Comp. 9 (1955), 118-120
    */

   /* t is coordinate in [-1,1] interval */
   t = (2.0 * x - s->a - s->b) / (s->b - s->a);
   t2 = 2.0 * t;

   d1 = 0.0;
   d2 = 0.0;

   for (i = order; i >= 1; i--)
     {
        double temp = d1;
        d1 = t2 * d1 - d2 + s->coef[i];
        d2 = temp;
     }

   value = t * d1 - d2 + s->coef[0];

   return value;
}

static int read_wavecal_params (int grp, const int *start, const int *count,
                                Cheb_Type *ct, double **wavecal_params, size_t *params_dimlen)
{
   TIO_Var_Info_Type info = {0};
   size_t len_params, pstart[3], pcount[3];
   double *params = NULL;
   int num_coef, start_pix, num_pix, status;

   *wavecal_params = NULL;
   *params_dimlen = 0;

   if (0 != TIO_inq_var (grp, TEMPO_VAR_WAVECAL_PARAM, &info))
     return -1;

   pstart[0] = start[0];
   pstart[1] = start[1];
   pstart[2] = 0;

   pcount[0] = count[0];
   pcount[1] = count[1];
   pcount[2] = info.dimlens[2];

   len_params = pcount[0] * pcount[1] * pcount[2];
   if (NULL == (params = (double *)TIO_MALLOC (len_params * sizeof(double))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   if ((NC_NOERR != (status = nc_get_vara_double (grp, info.varid, pstart, pcount, params)))
       || (NC_NOERR != (status = nc_get_att_int (grp, info.varid, "num_coefficients", &num_coef)))
       || (NC_NOERR != (status = nc_get_att_int (grp, info.varid, "start_spectral_channel", &start_pix)))
       || (NC_NOERR != (status = nc_get_att_int (grp, info.varid, "num_spectral_channels", &num_pix))))
     {
        tell_verror (TELL_IO_READ_ERROR, "%s: accessing variable %s in group %d (%s)",
                     __func__, TEMPO_VAR_WAVECAL_PARAM, grp, nc_strerror(status));
        TIO_FREE(params);
        return -1;
     }

   trace_get_var (grp, TEMPO_VAR_WAVECAL_PARAM, info.ndims, pstart, pcount);

   ct->a = start_pix;
   ct->b = start_pix + num_pix - 1;
   ct->num_coef = num_coef;

   *params_dimlen = info.dimlens[2];
   *wavecal_params = params;

   return 0;
}

static int get_wavelength_wavecal (int grp, const int *start, const int *count,
                                   int type, void *data, const char *name, void *client_data)
{
   Cheb_Type ct = {0};
   double *wavecal_params = NULL;
   size_t params_dimlen, wavelen_offset, param_offset;
   int step, num_waves, num_step, num_xtrack;

   (void) client_data; (void) name;

   num_step = count[0];
   num_xtrack = count[1];
   num_waves = count[2];

   if (0 != read_wavecal_params (grp, start, count, &ct, &wavecal_params, &params_dimlen))
     return -1;

   wavelen_offset = 0;
   param_offset = 0;

   for (step = 0; step < num_step; step++)
     {
        int xtrack;

        for (xtrack = 0; xtrack < num_xtrack; xtrack++)
          {
             int i;
             ct.coef = wavecal_params + param_offset;

             if (type == NC_FLOAT)
               {
                  float *y_flt = (float *)data + wavelen_offset;
                  for (i = 0; i < num_waves; i++)
                    {
                       y_flt[i] = (float) cheb_expansion_eval (&ct, i + start[2]);
                    }
               }
             else
               {
                  double *y_dbl = (double *)data + wavelen_offset;
                  for (i = 0; i < num_waves; i++)
                    {
                       y_dbl[i] = cheb_expansion_eval (&ct, i + start[2]);
                    }
               }
             wavelen_offset += num_waves;
             param_offset += params_dimlen;
          }
     }

   TIO_FREE(wavecal_params);

   return 0;
}

static int get_wavelength_nominal_1d (int grp, const int *start, const int *count,
                                      int type, void *data, const char *name,
                                      int adjust_nominal, void *client_data)
{
   float *nominal_waves = NULL;
   int step, num_waves, num_step, num_xtrack, varid, status;
   size_t wavelen_offset, wave_start, wave_count;

   (void) client_data; (void) name;

   num_step = count[0];
   num_xtrack = count[1];
   num_waves = count[2];

   if (NULL == (nominal_waves = (float *)TIO_MALLOC (num_waves * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   wave_start = start[2];
   wave_count = count[2];

   if ((NC_NOERR != (status = nc_inq_varid (grp, TEMPO_VAR_WAVELEN_NOMINAL, &varid)))
       || (NC_NOERR != (status = nc_get_vara_float (grp, varid, &wave_start, &wave_count, nominal_waves))))
     {
        TIO_FREE(nominal_waves);
        tell_verror (TELL_IO_READ_ERROR, "%s: failed reading %s (%s)",
                     __func__, TEMPO_VAR_WAVELEN_NOMINAL, nc_strerror (status));
        return -1;
     }

   trace_get_var (grp, TEMPO_VAR_WAVELEN_NOMINAL, 1, &wave_start, &wave_count);

   wavelen_offset = 0;

   for (step = 0; step < num_step; step++)
     {
        int xtrack;

        for (xtrack = 0; xtrack < num_xtrack; xtrack++)
          {
             int i;

             if (type == NC_FLOAT)
               {
                  float *y_flt = (float *)data + wavelen_offset;
                  if (adjust_nominal)
                    {
                       for (i = 0; i < num_waves; i++)
                         {
                            y_flt[i] += nominal_waves[i];
                         }
                    }
                  else
                    {
                       memcpy ((char *)y_flt, (char *)nominal_waves,
                               num_waves * sizeof(float));
                    }
               }
             else
               {
                  double *y_dbl = (double *)data + wavelen_offset;
                  if (adjust_nominal)
                    {
                       for (i = 0; i < num_waves; i++)
                         {
                            y_dbl[i] += (double) nominal_waves[i];
                         }
                    }
                  else
                    {
                       for (i = 0; i < num_waves; i++)
                         {
                            y_dbl[i] = (double) nominal_waves[i];
                         }
                    }
               }
             wavelen_offset += num_waves;
          }
     }

   TIO_FREE(nominal_waves);

   return 0;
}

static int get_wavelength_nominal_2d (int grp, const int *start, const int *count,
                                      int type, void *data, const char *name,
                                      int adjust_nominal, void *client_data)
{
   float *nominal_waves = NULL;
   int step, num_waves, num_step, num_xtrack, varid, status;
   size_t wavelen_offset, nw_offset, nw_start[2], nw_count[2];

   (void) client_data; (void) name;

   num_step = count[0];
   num_xtrack = count[1];
   num_waves = count[2];

   if (NULL == (nominal_waves = (float *)TIO_MALLOC (num_xtrack * num_waves * sizeof(float))))
     {
        tell_verror (TELL_MALLOC_ERROR, "%s: malloc failed", __func__);
        return -1;
     }

   nw_start[0] = start[1];
   nw_start[1] = start[2];
   nw_count[0] = count[1];
   nw_count[1] = count[2];

   if ((NC_NOERR != (status = nc_inq_varid (grp, TEMPO_VAR_WAVELEN_NOMINAL, &varid)))
       || (NC_NOERR != (status = nc_get_vara_float (grp, varid, nw_start, nw_count, nominal_waves))))
     {
        TIO_FREE(nominal_waves);
        tell_verror (TELL_IO_READ_ERROR, "%s: failed reading %s (%s)",
                     __func__, TEMPO_VAR_WAVELEN_NOMINAL, nc_strerror (status));
        return -1;
     }

   trace_get_var (grp, TEMPO_VAR_WAVELEN_NOMINAL, 2, nw_start, nw_count);

   wavelen_offset = 0;

   for (step = 0; step < num_step; step++)
     {
        int xtrack;

        nw_offset = 0;
        for (xtrack = 0; xtrack < num_xtrack; xtrack++)
          {
             int i;

             if (type == NC_FLOAT)
               {
                  float *y_flt = (float *)data + wavelen_offset;
                  float *nw_x = nominal_waves + nw_offset;
                  if (adjust_nominal)
                    {
                       for (i = 0; i < num_waves; i++)
                         {
                            y_flt[i] += nw_x[i];
                         }
                    }
                  else
                    {
                       memcpy ((char *)y_flt, (char *)nw_x, num_waves * sizeof(float));
                    }
               }
             else
               {
                  double *y_dbl = (double *)data + wavelen_offset;
                  if (adjust_nominal)
                    {
                       for (i = 0; i < num_waves; i++)
                         {
                            y_dbl[i] += (double) nominal_waves[nw_offset + i];
                         }
                    }
                  else
                    {
                       for (i = 0; i < num_waves; i++)
                         {
                            y_dbl[i] = (double) nominal_waves[nw_offset + i];
                         }
                    }
               }
             wavelen_offset += num_waves;
             nw_offset += num_waves;
          }
     }

   TIO_FREE(nominal_waves);

   return 0;
}

static int have_wavecal_param_3d (int grp, int *adjust_nominal)
{
   int varid, ndims, dimids[3], dimid_par, adjust_att;

   *adjust_nominal = 0;

   if (NC_NOERR != nc_inq_varid (grp, TEMPO_VAR_WAVECAL_PARAM, &varid))
     return 0;

   if ((NC_NOERR != nc_inq_varndims (grp, varid, &ndims))
       || (ndims != 3))
     return 0;

   if ((NC_NOERR != nc_inq_vardimid (grp, varid, dimids))
       || (NC_NOERR != nc_inq_dimid (grp, TEMPO_DIM_WAVECAL_PARAM, &dimid_par)))
     return 0;

   if (dimids[2] != dimid_par)
     return 0;

   /* For back-compatibility, this attribute is optional */
   if (NC_NOERR == nc_get_att_int (grp, varid, "adjust_nominal_wavelength", &adjust_att))
     {
        *adjust_nominal = adjust_att;
     }

   return 1;
}

static int have_nominal_wavelength (int grp)
{
   int varid, ndims, dimids[2], dimid_spectral_channel, dimid_xtrack;

   if (NC_NOERR != nc_inq_varid (grp, TEMPO_VAR_WAVELEN_NOMINAL, &varid))
     return 0;

   if ((NC_NOERR != nc_inq_varndims (grp, varid, &ndims))
       || (ndims != TIO_NOMINAL_WAVELEN_NUM_DIMS))
     return 0;

   if (NC_NOERR != nc_inq_dimid (grp, TEMPO_DIM_CHANNEL, &dimid_spectral_channel))
     return 0;

   if (ndims == 1)
     {
        if (NC_NOERR != nc_inq_vardimid (grp, varid, dimids))
            return 0;
        return (dimids[0] == dimid_spectral_channel);
     }

   if (ndims == 2)
     {
        if (NC_NOERR != nc_inq_dimid (grp, TEMPO_DIM_XTRACK, &dimid_xtrack))
          return 0;

        if (NC_NOERR != nc_inq_vardimid (grp, varid, dimids))
          return 0;
        return ((dimids[0] == dimid_xtrack) && (dimids[1] == dimid_spectral_channel));
     }

   return 0;
}

static int get_wavelength (int grp, const int *start, const int *count,
                           int type, void *data, const char *name, void *client_data)
{
   int adjust_nominal_wavelength, have_wavecal, have_nominal;
   int status;

   /* When asked to read 'wavelength':
    * First, look for 3D TEMPO_VAR_WAVECAL_PARAM,
    * if that's not present, look for TEMPO_VAR_WAVELEN_NOMINAL.
    * If that's not present, return >0, which means look for 'wavelength'
    * And just to complicate things further, the wavelength calibration parameters
    * can define an adjustment to the nominal wavelength grid, so we have to
    * handle that case too.
    */

   have_wavecal = have_wavecal_param_3d (grp, &adjust_nominal_wavelength);
   have_nominal = have_nominal_wavelength (grp);

   if (have_wavecal)
     {
        if (0 != (status = get_wavelength_wavecal (grp, start, count, type, data, name, client_data)))
          return status;

        /* If the wavelength calibration fully defines the wavelength grid then we're done */
        if (adjust_nominal_wavelength == 0)
          return status;

        /* If the wavelength calibration adjusts missing nominal wavelengths, then that's a problem */
        if (have_nominal == 0)
          {
             tell_verror (TELL_IO_READ_ERROR, "%s: variable not found: %s", __func__, TEMPO_VAR_WAVELEN_NOMINAL);
             return -1;
          }
        /* FALLTHRU when wavelength calibration adjusts the existing nominal grid */
     }

   if (have_nominal)
     {
        if (1 == TIO_NOMINAL_WAVELEN_NUM_DIMS)
          {
             return get_wavelength_nominal_1d (grp, start, count, type, data, name,
                                               adjust_nominal_wavelength, client_data);
          }
        else
          {
             return get_wavelength_nominal_2d (grp, start, count, type, data, name,
                                               adjust_nominal_wavelength, client_data);
          }
     }

   return 1;
}

typedef struct
{
   char *name;
   int (*get)(int, const  int *, const int *,
              int,       void *, const char *, void *);
   void *get_client_data;
   int get_enable;
   int (*put)(int, const  int *, const int *,
              int, const void *, const char *, void *);
   void *put_client_data;
   int put_enable;
}
IO_Methods_Type;
#define IO_METHODS_END {NULL, NULL,NULL,0, NULL,NULL,0}

static
IO_Methods_Type *find_io_methods (const char *name,
                                  IO_Methods_Type *io_methods)
{
   IO_Methods_Type *m;

   if ((name == NULL) || (io_methods == NULL))
     return NULL;

   for (m = io_methods; m->name != NULL; m++)
     {
        if (0 == strcmp (name, m->name))
          return m;
     }

   return NULL;
}

static int Radiance_Error_Num_Significant_Digits = 2;

static IO_Methods_Type IO_Methods[] =
{
   {TEMPO_VAR_WAVELENGTH,
        get_wavelength, NULL, 1,
        NULL, NULL, 0},
   {TEMPO_VAR_RADIANCE_ERROR,
        NULL, NULL, 0,
        put_float_nsd, &Radiance_Error_Num_Significant_Digits, 1},
   {TEMPO_VAR_IRRADIANCE_ERROR,
        NULL, NULL, 0,
        put_float_nsd, &Radiance_Error_Num_Significant_Digits, 1},
   IO_METHODS_END
};

int _TIO_set_io_method_enable (const char *name,
                               int get_enable, int put_enable)
{
   IO_Methods_Type *io_method;

   io_method = find_io_methods (name, IO_Methods);
   if (io_method == NULL)
     return -1;
   io_method->get_enable = get_enable;
   io_method->put_enable = put_enable;
   return 0;
}

#define TIO_IO_VAR_SECTION(action,error_num,const_qual) \
int TIO_##action##_var_section (int grp, const char *name, \
                                const int *istart, const int *icount, int type, \
                                const_qual void *data) \
{ \
   TIO_Var_Info_Type info; \
   int status, varid, ndims, i; \
   size_t start[TIO_MAX_VAR_DIMS], count[TIO_MAX_VAR_DIMS]; \
   IO_Methods_Type *io_method; \
 \
   io_method = find_io_methods (name, IO_Methods); \
   if ((NULL != io_method) && (io_method->action##_enable != 0) \
       && (NULL != io_method->action)) \
     { \
        int io_method_status; \
        io_method_status = io_method->action (grp, istart, icount, type, data, name, \
                                              io_method->action##_client_data); \
        if (io_method_status <= 0) \
          return io_method_status; \
     } \
 \
   if (-1 == TIO_inq_var (grp, name, &info)) \
     return -1; \
 \
   ndims = info.ndims; \
   varid = info.varid; \
 \
   if (ndims < 0) \
     { \
        Tell_verror (TELL_INVALID_PARM, "%s: variable %s has ndims=%d", \
                     __func__, name, ndims); \
        return -1; \
     } \
 \
   for (i = 0; i < TIO_MAX_VAR_DIMS; i++) \
     { \
        start[i] = 0; \
        count[i] = 0; \
     } \
 \
   for (i = 0; i < ndims; i++) \
     { \
        start[i] = (size_t) istart[i]; \
        count[i] = (icount[i] >= 0) ? (size_t) icount[i] : info.dimlens[i]; \
     } \
 \
   switch (type) \
     { \
      case NC_CHAR: \
        status = nc_##action##_vara_text (grp, varid, start, count, (const_qual char *)data); \
        break; \
      case NC_BYTE: \
        status = nc_##action##_vara_schar (grp, varid, start, count, (const_qual signed char *)data); \
        break; \
      case NC_UBYTE: \
        status = nc_##action##_vara_uchar (grp, varid, start, count, (const_qual unsigned char *)data); \
        break; \
      case NC_SHORT: \
        status = nc_##action##_vara_short (grp, varid, start, count, (const_qual short *)data); \
        break; \
      case NC_USHORT: \
        status = nc_##action##_vara_ushort (grp, varid, start, count, (const_qual unsigned short *)data); \
        break; \
      case NC_INT: \
        status = nc_##action##_vara_int (grp, varid, start, count, (const_qual int *)data); \
        break; \
      case NC_UINT: \
        status = nc_##action##_vara_uint (grp, varid, start, count, (const_qual unsigned int *)data); \
        break; \
      case NC_INT64: \
        status = nc_##action##_vara_longlong (grp, varid, start, count, (const_qual long long *)data); \
        break; \
      case NC_UINT64: \
        status = nc_##action##_vara_ulonglong (grp, varid, start, count, (const_qual unsigned long long *)data); \
        break; \
      case NC_FLOAT: \
        status = nc_##action##_vara_float (grp, varid, start, count, (const_qual float *)data); \
        break; \
      case NC_DOUBLE: \
        status = nc_##action##_vara_double (grp, varid, start, count, (const_qual double *)data); \
        break; \
      case NC_STRING: \
        status = nc_##action##_vara_string (grp, varid, start, count, (const_qual char **)data); \
        break; \
      default: \
        status = nc_##action##_vara (grp, varid, start, count, (const_qual void *)data); \
        break; \
     } \
 \
   if (status != NC_NOERR) \
     { \
        Tell_verror (error_num, \
                     "%s: accessing variable %s in group %d (%s)", \
                     __func__, name, grp, nc_strerror(status)); \
        return -1; \
     } \
   \
   trace_##action##_var (grp, name, ndims, start, count); \
 \
   return 0; \
}

TIO_IO_VAR_SECTION(get,TELL_IO_READ_ERROR,EMPTY())
TIO_IO_VAR_SECTION(put,TELL_IO_WRITE_ERROR,const)

#if 0
}
#endif

int TIO_inq_att (int grp, int varid, const char *attname,
                 int *xtype, int *len)
{
   size_t len_size_t;
   int status;

   if (NULL == attname)
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_att (grp, varid, attname, xtype, &len_size_t)))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: accessing attribute %s (varid=%d) (%s)",
                     __func__, attname, varid, nc_strerror(status));
        return -1;
     }

   if (len) *len = len_size_t;

   return 0;
}

int TIO_put_att (int grp, int varid, const char *attname,
                 int xtype, int len, const void *att)
{
   int status;
   int file_atttype;
   size_t len_size_t;

   if ((NULL == attname) || (att == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   /* Enforce attribute type-checking.
    * By default, netCDF will happily overwrite an existing
    * attribute with a different type value, but I don't think
    * we want that behavior.
    */
   status = nc_inq_atttype (grp, varid, attname, &file_atttype);
   if (NC_NOERR == status)
     {
        /* if the attribute exists, make sure we're writing
         * a value of the same type.
         */
        if (xtype != file_atttype)
          {
             Tell_verror (TELL_INVALID_PARM,
                          "%s: type mismatch: writing attribute %s: value type=%d, file value type=%d",
                          __func__, attname, xtype, file_atttype);
             return -1;
          }
     }
   else if (NC_ENOTATT != status)
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: attribute %s (varid=%d) query failed (%s)",
                     __func__, attname, varid, nc_strerror(status));
        return -1;
     }

   len_size_t = len;
   status = nc_put_att (grp, varid, attname, xtype, len_size_t, att);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: writing attribute %s (varid=%d) (%s)",
                     __func__, attname, varid, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_get_att (int grp, int varid, const char *attname,
                 int xtype, void *att)
{
   int status;

   if ((NULL == attname) || (att == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   /* Use the netCDF type-safe interface so that netCDF
    * performs any needed type conversion for built-in types.
    */
   switch (xtype)
     {
      case NC_CHAR:
        status = nc_get_att_text (grp, varid, attname, (char *)att);
        break;
      case NC_BYTE:
        status = nc_get_att_schar (grp, varid, attname, (signed char *)att);
        break;
      case NC_UBYTE:
        status = nc_get_att_uchar (grp, varid, attname, (unsigned char *)att);
        break;
      case NC_SHORT:
        status = nc_get_att_short (grp, varid, attname, (short *)att);
        break;
      case NC_USHORT:
        status = nc_get_att_ushort (grp, varid, attname, (unsigned short *)att);
        break;
      case NC_INT:
        status = nc_get_att_int (grp, varid, attname, (int *)att);
        break;
      case NC_UINT:
        status = nc_get_att_uint (grp, varid, attname, (unsigned int *)att);
        break;
      case NC_INT64:
        status = nc_get_att_longlong (grp, varid, attname, (long long *)att);
        break;
      case NC_UINT64:
        status = nc_get_att_ulonglong (grp, varid, attname, (unsigned long long *)att);
        break;
      case NC_FLOAT:
        status = nc_get_att_float (grp, varid, attname, (float *)att);
        break;
      case NC_DOUBLE:
        status = nc_get_att_double (grp, varid, attname, (double *)att);
        break;
      case NC_STRING:
        status = nc_get_att_string (grp, varid, attname, (char **)att);
        break;
      default:
        if (NC_NOERR != (status = nc_get_att (grp, varid, attname, att)))
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: reading attribute %s in group %d using unrecognized typeid=%d (%s)",
                          __func__, attname, grp, xtype, nc_strerror(status));
             return -1;
          }
        /* drop */
     }

   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_READ_ERROR, "%s: reading attribute %s in group %d (%s)",
                     __func__, attname, grp, nc_strerror(status));
        return -1;
     }

   return 0;
}

int TIO_free_string (size_t len, char **data)
{
   int status;

   if (NC_NOERR != (status = nc_free_string (len, data)))
     {
        Tell_verror (TELL_RUNTIME_ERROR, "%s: freeing string %ld arrays (%s)",
                     __func__, len, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_inq_grp (int parent_ncid, const char *path, int *grp)
{
   int status;

   if ((path == NULL) || (grp == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_grp_full_ncid (parent_ncid, path, grp)))
     {
        Tell_verror (TELL_IO_READ_ERROR, "%s: reading group path %s (%s)",
                     __func__, path, nc_strerror(status));
        return -1;
     }
   _pTIO_trace_group (parent_ncid, path, *grp);

   return 0;
}

int TIO_def_grp (int parent_ncid, const char *path, int *new_ncid)
{
   char delim = '/';
   const char *p;
   unsigned int len;

   if ((path == NULL) || (new_ncid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if ((*path == delim) && (*(path+1) == 0))
     {
        int status = nc_inq_grp_full_ncid (parent_ncid, path, new_ncid);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: accessing id of root group (%s)",
                          __func__, nc_strerror (status));
             return -1;
          }
        _pTIO_trace_group (parent_ncid, path, *new_ncid);
        return 0;
     }

   for (p = path; *p != 0; p += len)
     {
        char buf[TIO_MAX_NAME_LEN];
        int status, ncid;
        const char *end;

        while (*p == delim) p++;

        end = strchr (p, delim);
        if (end == NULL)
          {
             end = p + strlen(p);
          }
        len = end - p;

        if (len == 0) break;

        if (len + 1 > sizeof(buf))
          {
             Tell_verror (TELL_INVALID_PARM, "%s: group name is too long: %s", __func__, p);
             return -1;
          }
        strncpy (buf, p, len);
        buf[len] = 0;

        /* Does this group exist already? */
        status = nc_inq_grp_ncid (parent_ncid, buf, &ncid);
        if (NC_ENOGRP == status)
          {
             /* If the group doesn't exist, create it */
             status = nc_def_grp (parent_ncid, buf, &ncid);
             if (NC_NOERR != status)
               {
                  Tell_verror (TELL_IO_WRITE_ERROR, "%s: creating group %s (%s)",
                               __func__, buf, nc_strerror (status));
                  return -1;
               }
          }
        else if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR, "%s: accessing group %s (%s)",
                          __func__, buf, nc_strerror (status));
             return -1;
          }

        parent_ncid = ncid;
        *new_ncid = ncid;
     }

   _pTIO_trace_group (parent_ncid, path, *new_ncid);

   return 0;
}

int TIO_def_dim (int ncid, const char *name, size_t len, int *id)
{
   int status;

   if ((name == NULL) || (id == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_dim (ncid, name, len, id)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining dimension %s (%s)",
                     __func__, name, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_def_var (int ncid, const char *name, int type,
                 int num_dims, const int *dimids, int *varid)
{
   int status;

   if ((name == NULL) || (varid == NULL)
       || ((dimids == NULL) && (num_dims != 0)))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_def_var (ncid, name, type, num_dims, dimids, varid)))
     {
        Tell_verror (TELL_IO_WRITE_ERROR, "%s: defining variable %s (%s)",
                     __func__, name, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_def_var_fill (int grp, int varid, int no_fill, const void *fill_value)
{
   int status;

   status = nc_def_var_fill (grp, varid, no_fill, fill_value);
   if (status != NC_NOERR)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: setting fill value for varid=%d (%s)",
                     __func__, varid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_inq_var_fill (int grp, int varid, int *no_fill, void *fill_value)
{
   int status;

   /* Note: nc_inq_var_fill doesn't correctly return the fill value
    *       in the file, but nc_get_att does.  Weird.
    */

   if ((NC_NOERR != (status = nc_inq_var_fill (grp, varid, no_fill, NULL)))
       || (NC_NOERR != (status = nc_get_att (grp, varid, _FillValue, fill_value))))
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: getting fill value for varid=%d (%s)",
                     __func__, varid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_def_var_chunking (int ncid, int varid,
                          int storage, size_t *chunksizep)
{
   int status;

   status = nc_def_var_chunking (ncid, varid, storage, chunksizep);
   if (status != NC_NOERR)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: setting chunking for varid=%d (%s)",
                     __func__, varid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_def_var_deflate (int grp, int varid,
                         int shuffle, int deflate, int deflate_level)
{
   int status;

   status = nc_def_var_deflate (grp, varid, shuffle, deflate, deflate_level);
   if (status != NC_NOERR)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: setting deflate values for varid=%d (%s)",
                     __func__, varid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_put_text_attrs (int grp, int varid, const TIO_Attr_Text_Type *attrs)
{
   const TIO_Attr_Text_Type *a;
   int status;

   if (attrs == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   for (a = attrs; a->name != NULL; a++)
     {
        status = nc_put_att_text (grp, varid, a->name, strlen(a->text), a->text);
        if (status != NC_NOERR)
          {
             Tell_verror (TELL_IO_WRITE_ERROR, "%s: writing text attribute %s (%s)",
                          __func__, a->name, nc_strerror (status));
             return -1;
          }
     }

   return 0;
}

int TIO_copy_attrs (int ncid_infile, int id_var_infile,
                    int (*dontcopy_attr)(const char *),
                    int ncid, int id_var)
{
   char attname[TIO_MAX_NAME_LEN];
   int status, attnum, num_atts;

   status = nc_inq_varnatts (ncid_infile, id_var_infile, &num_atts);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_READ_ERROR,
                     "%s: reading number of attributes for variable id=%d (%s)",
                     __func__, id_var_infile, nc_strerror (status));
        return -1;
     }

   for (attnum = 0; attnum < num_atts; attnum++)
     {
        status = nc_inq_attname (ncid_infile, id_var_infile, attnum, attname);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: reading attribute %d=%s (%s)",
                          __func__, attnum, attname,
                          nc_strerror (status));
             return -1;
          }
        if ((NULL != dontcopy_attr) && (1 == dontcopy_attr (attname)))
          continue;
        status = nc_copy_att (ncid_infile, id_var_infile, attname,
                              ncid, id_var);
        if (NC_NOERR != status)
          {
             Tell_verror (TELL_IO_READ_ERROR,
                          "%s: copying attribute %s (%s)",
                          __func__, attname, nc_strerror (status));
             return -1;
          }
     }

   return 0;
}

int TIO_inq_dimname (int grp, int dimid, char *dimname)
{
   int status;

   if (dimname == NULL)
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NC_NOERR != (status = nc_inq_dimname (grp, dimid, dimname)))
     {
        Tell_verror (TELL_IO_READ_ERROR, "%s: reading dimension %d name (%s)",
                     __func__, dimid, nc_strerror (status));
        return -1;
     }

   return 0;
}

int TIO_inq_dimid (int grp, const char *dimname, int *dimid)
{
   if ((dimname == NULL) || (dimid == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   /* For such a low-level operation, it seems better to
    * _not_ generate an error message by default */
   return nc_inq_dimid (grp, dimname, dimid);
}

int TIO_inq_dim (int grp, const char *dimname, int *dimid, size_t *dimlen)
{
   int status;
   if ((dimname == NULL) || (dimid == NULL) || (dimlen == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   /* For such a low-level operation, it seems better to
    * _not_ generate an error message by default */

   status = nc_inq_dimid (grp, dimname, dimid);
   if (NC_NOERR != status)
     return status;

   return nc_inq_dimlen (grp, *dimid, dimlen);
}

int TIO_put_git_commit_hash (int grp, const char *attname)
{
   const char hash[] = GIT_COMMIT_HASH ;
   const char *name = "tio_commit";
   int status;

   /* Use a default attribute name when the supplied attribute
    * name is either NULL, empty, or begins with a space.
    */
   if ((attname != NULL) && (*attname != 0) && (*attname != ' '))
     name = attname;

   status = nc_put_att_text (grp, NC_GLOBAL, name, strlen(hash), hash);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: writing git commit hash (%s)",
                     __func__, nc_strerror(status));
        return -1;
     }

   return 0;
}

#define CONVERT_TO(to_type_name,to_type) \
static int convert_to_##to_type (int from_type, void *from_value, to_type_name *to_value) \
{ \
   to_type_name value; \
 \
   switch (from_type) \
     { \
      case NC_BYTE: value = *(signed char *)from_value; break; \
      case NC_SHORT: value = *(short *)from_value; break; \
      case NC_INT: value = *(int *)from_value; break; \
      case NC_INT64: value = *(long long *)from_value; break; \
      case NC_UBYTE: value = *(unsigned char *)from_value; break; \
      case NC_USHORT: value = *(unsigned short *)from_value; break; \
      case NC_UINT: value = *(unsigned int *)from_value; break; \
      case NC_UINT64: value = *(unsigned long long *)from_value; break; \
      case NC_FLOAT: value = *(float *)from_value; break; \
      case NC_DOUBLE: value = *(double *)from_value; break; \
      default: \
        Tell_verror (TELL_INVALID_PARM, "%s: unsupported type from_type=%d\n", \
                     __func__, from_type); \
        return -1; \
     } \
 \
   *to_value = value; \
 \
   return 0; \
}
CONVERT_TO(signed char,byte)
CONVERT_TO(unsigned char,ubyte)
CONVERT_TO(short,short)
CONVERT_TO(unsigned short,ushort)
CONVERT_TO(int,int)
CONVERT_TO(unsigned int,uint)
CONVERT_TO(long long,int64)
CONVERT_TO(unsigned long long,uint64)
CONVERT_TO(float,float)
CONVERT_TO(double,double)

int TIO_get_fill_value (int grp, const char *name, int type, void *value)
{
   int varid, file_type, status, return_status = -1;
   void *fill = NULL;

   if ((name == NULL) || (value == NULL))
     {
        Tell_verror (TELL_INVALID_PARM, "%s: got a NULL pointer", __func__);
        return -1;
     }

   if (NULL == (fill = TIO_MALLOC (8)))
     {
        Tell_set_error (TELL_MALLOC_ERROR);
        return -1;
     }

   /* Note: nc_inq_var_fill doesn't correctly return the fill value
    *       in the file, but nc_get_att does.  Weird.
    */

   if ((NC_NOERR != (status = nc_inq_varid (grp, name, &varid)))
       || (NC_NOERR != (status = nc_inq_vartype (grp, varid, &file_type)))
       || (NC_NOERR != (status = nc_get_att (grp, varid, _FillValue, fill)))
      )
     {
        Tell_verror (TELL_RUNTIME_ERROR, "%s: getting fill value for %s (%s)\n",
                     __func__, name, nc_strerror(status));
        goto cleanup;
     }

   switch (type)
     {
      case NC_BYTE:
        status = convert_to_byte (file_type, fill, (signed char *)value);
        break;
      case NC_UBYTE:
        status = convert_to_ubyte (file_type, fill, (unsigned char *)value);
        break;
      case NC_SHORT:
        status = convert_to_short (file_type, fill, (short *)value);
        break;
      case NC_USHORT:
        status = convert_to_ushort (file_type, fill, (unsigned short *)value);
        break;
      case NC_INT:
        status = convert_to_int (file_type, fill, (int *)value);
        break;
      case NC_UINT:
        status = convert_to_uint (file_type, fill, (unsigned int *)value);
        break;
      case NC_INT64:
        status = convert_to_int64 (file_type, fill, (long long *)value);
        break;
      case NC_UINT64:
        status = convert_to_uint64 (file_type, fill, (unsigned long long *)value);
        break;
      case NC_FLOAT:
        status = convert_to_float (file_type, fill, (float *)value);
        break;
      case NC_DOUBLE:
        status = convert_to_double (file_type, fill, (double *)value);
        break;
      default:
        Tell_verror (TELL_INVALID_PARM, "%s: unsupported destination type = %d\n",
                     __func__, type);
        goto cleanup;
     }

   if (status)
     {
        Tell_verror (TELL_RUNTIME_ERROR, "%s: converting %s fill value to type = %d\n",
                     __func__, name, type);
        goto cleanup;
     }

   return_status = 0;
cleanup:
   free (fill);
   return return_status;
}

/* Fortran bindings */

FCALLSCFUN6(INT, TIO_get_var_section, TIOF_GET_VAR_SECTION, tiof_get_var_section,
            INT, STRING, PINT, PINT, INT, PVOID)
FCALLSCFUN6(INT, TIO_put_var_section, TIOF_PUT_VAR_SECTION, tiof_put_var_section,
            INT, STRING, PINT, PINT, INT, PVOID)
FCALLSCFUN2(INT, TIO_put_git_commit_hash, TIO_F_PUT_GIT_HASH, tio_f_put_git_hash,
            INT, STRING)
FCALLSCFUN3(INT, TIO_def_grp, TIO_F_DEF_GRP, tio_f_def_grp,
            INT, STRING, PINT)
FCALLSCFUN4(INT, TIO_get_fill_value, TIO_F_GET_FILL_VALUE, tio_f_get_fill_value,
            INT, STRING, INT, PVOID)
