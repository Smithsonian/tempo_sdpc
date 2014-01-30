#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <hdf.h>		       /* define DFNT_* constants */
#include <cfortHdf.h>
#include <HdfEosDef.h>

#include <stdarg.h>
#include "cerr.h"

void cerr_verror (char *fmt, ...)
{
   va_list ap;
   va_start (ap, fmt);

   (void) fprintf (stderr, "**ERROR: ");
   (void) vfprintf (stderr, fmt, ap);
   (void) fprintf (stderr, "\n");
   va_end (ap);
}

/* arbitrary, hopefully overkill.  Unfortunately the underlying library
 * does not specify sizes.  Sigh.
 */
#define MAX_DIMS 32
#define MAX_STRDIMS_LEN	1024

typedef int Swath_Id_Type;
typedef struct _L1b_Info_Type L1b_Info_Type;

typedef struct
{
   int datatype;
   int rank;
   int32 start[MAX_DIMS];
   int32 edge[MAX_DIMS];
   size_t total_num;
}
L1b_Data_Info_Type;

typedef int (*Read_Func_Type)(L1b_Info_Type *, Swath_Id_Type,
			      int32 track, int32 ntracks,
			      int datatype, void **data,
			      L1b_Data_Info_Type *dinfo);

struct _L1b_Info_Type
{
   char *l1bname;		       /* name in L1b file */
   int num_dims;
   Read_Func_Type readfunc;
};

static L1b_Info_Type *find_l1b_info (char *name);

typedef struct
{
   int datatype;
   size_t sizeof_datatype;
}
Type_Info_Type;

static Type_Info_Type Type_Info_Table[] =
{
   {DFNT_INT8, sizeof(char)},
   {DFNT_UINT8, sizeof(unsigned char)},
   {DFNT_INT16, sizeof(int16)},
   {DFNT_UINT16, sizeof(uint16)},
   {DFNT_INT32, sizeof(int32)},
   {DFNT_UINT32, sizeof(uint32)},
   /* {DFNT_INT64, sizeof(int64)}, */
   /* {DFNT_UINT64, sizeof(uint64)}, */
   {DFNT_FLOAT32, sizeof(float)},
   {DFNT_FLOAT64, sizeof(double)},
   {-1, 0}
};

static Type_Info_Type *find_type_info (int datatype)
{
   Type_Info_Type *t = Type_Info_Table;

   while (t->datatype != -1)
     {
	if (t->datatype == datatype) return t;
	t++;
     }
   cerr_verror ("Unable to find type information for type %d", datatype);
   return NULL;
}

static void *alloc_data_array (int type, size_t total_num)
{
   Type_Info_Type *t;
   void *data;

   if (NULL == (t = find_type_info (type)))
     return NULL;

   if (total_num == 0) total_num = 1;
   if (NULL == (data = malloc (t->sizeof_datatype * total_num)))
     cerr_verror ("Out of memory in alloc_data_array");

   return data;
}

typedef int (*Converter_Func_Type)(int, void *, int, void *, size_t);
typedef struct
{
   int from_type, to_type;
   Converter_Func_Type convert;
}
Converter_Type;

static int convert_f32_f64 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   float32 *f32 = (float32 *)ap;
   float64 *f64 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f64[i] = f32[i];
   return 0;
}

static int convert_f64_f32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   float64 *f64 = (float64 *)ap;
   float32 *f32 = (float32 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f32[i] = (float32) f64[i];
   return 0;
}

static int convert_i32_f64 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int32 *i32 = (int32 *)ap;
   float64 *f64 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f64[i] = i32[i];
   return 0;
}

static int convert_i32_f32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int32 *i32 = (int32 *)ap;
   float64 *f32 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f32[i] = i32[i];
   return 0;
}

static int convert_i16_f64 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int16 *i16 = (int16 *)ap;
   float64 *f64 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f64[i] = i16[i];
   return 0;
}

static int convert_i16_f32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int16 *i16 = (int16 *)ap;
   float64 *f32 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f32[i] = i16[i];
   return 0;
}

static int convert_i16_i32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int16 *i16 = (int16 *)ap;
   int32 *i32 = (int32 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) i32[i] = i16[i];
   return 0;
}

static int convert_i8_f64 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int16 *i8 = (int16 *)ap;
   float64 *f64 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f64[i] = i8[i];
   return 0;
}

static int convert_i8_f32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int16 *i8 = (int16 *)ap;
   float64 *f32 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f32[i] = i8[i];
   return 0;
}

static int convert_i8_i32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int16 *i8 = (int16 *)ap;
   int32 *i32 = (int32 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) i32[i] = i8[i];
   return 0;
}


static int convert_i8_i16 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   int8 *i8 = (int8 *)ap;
   int16 *i16 = (int16 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) i16[i] = i8[i];
   return 0;
}

/* unsigned conversions */
static int convert_u32_f64 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint32 *u32 = (uint32 *)ap;
   float64 *f64 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f64[i] = u32[i];
   return 0;
}

static int convert_u32_f32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint32 *u32 = (uint32 *)ap;
   float64 *f32 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f32[i] = u32[i];
   return 0;
}

static int convert_u16_f64 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint16 *u16 = (uint16 *)ap;
   float64 *f64 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f64[i] = u16[i];
   return 0;
}

static int convert_u16_f32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint16 *u16 = (uint16 *)ap;
   float64 *f32 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f32[i] = u16[i];
   return 0;
}

static int convert_u16_i32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint16 *u16 = (uint16 *)ap;
   int32 *i32 = (int32 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) i32[i] = u16[i];
   return 0;
}

static int convert_u8_f64 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint16 *u8 = (uint16 *)ap;
   float64 *f64 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f64[i] = u8[i];
   return 0;
}

static int convert_u8_f32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint16 *u8 = (uint16 *)ap;
   float64 *f32 = (float64 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) f32[i] = u8[i];
   return 0;
}

static int convert_u8_i32 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint16 *u8 = (uint16 *)ap;
   int32 *i32 = (int32 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) i32[i] = u8[i];
   return 0;
}


static int convert_u8_i16 (int from, void *ap, int to, void *bp, size_t n)
{
   size_t i;
   uint8 *u8 = (uint8 *)ap;
   int16 *i16 = (int16 *)bp;
   (void) from; (void) to;
   for (i = 0; i < n; i++) i16[i] = u8[i];
   return 0;
}

static int convert_u8_i8 (int from, void *ap, int to, void *bp, size_t n)
{
   (void) from; (void) to;
   memcpy (bp, ap, n*sizeof(int8));
   return 0;
}

static int convert_u16_i16 (int from, void *ap, int to, void *bp, size_t n)
{
   (void) from; (void) to;
   memcpy (bp, ap, n*sizeof(int16));
   return 0;
}

static int convert_u32_i32 (int from, void *ap, int to, void *bp, size_t n)
{
   (void) from; (void) to;
   memcpy (bp, ap, n*sizeof(int32));
   return 0;
}


static Converter_Type Converter_Table[] =
{
   {DFNT_FLOAT32, DFNT_FLOAT64, convert_f32_f64},
   {DFNT_FLOAT64, DFNT_FLOAT32, convert_f64_f32},

   {DFNT_INT8, DFNT_INT16, convert_i8_i16},
   {DFNT_INT8, DFNT_INT32, convert_i8_i32},
   {DFNT_INT8, DFNT_FLOAT32, convert_i8_f32},
   {DFNT_INT8, DFNT_FLOAT64, convert_i8_f64},

   {DFNT_INT16, DFNT_INT32, convert_i16_i32},
   {DFNT_INT16, DFNT_FLOAT32, convert_i16_f32},
   {DFNT_INT16, DFNT_FLOAT64, convert_i16_f64},

   {DFNT_INT32, DFNT_FLOAT32, convert_i32_f32},
   {DFNT_INT32, DFNT_FLOAT64, convert_i32_f64},

   {DFNT_UINT8, DFNT_INT8, convert_u8_i8},
   {DFNT_UINT8, DFNT_INT16, convert_u8_i16},
   {DFNT_UINT8, DFNT_INT32, convert_u8_i32},
   {DFNT_UINT8, DFNT_FLOAT32, convert_u8_f32},
   {DFNT_UINT8, DFNT_FLOAT64, convert_u8_f64},

   {DFNT_UINT16, DFNT_INT16, convert_u16_i16},
   {DFNT_UINT16, DFNT_INT32, convert_u16_i32},
   {DFNT_UINT16, DFNT_FLOAT32, convert_u16_f32},
   {DFNT_UINT16, DFNT_FLOAT64, convert_u16_f64},

   {DFNT_UINT32, DFNT_INT32, convert_u32_i32},
   {DFNT_UINT32, DFNT_FLOAT32, convert_u32_f32},
   {DFNT_UINT32, DFNT_FLOAT64, convert_u32_f64},
   {-1, -1, NULL}
};

static Converter_Type *find_type_converter (int from_type, int to_type)
{
   Converter_Type *c = Converter_Table;

   while (c->convert != NULL)
     {
	if ((c->to_type == to_type) && (c->from_type == from_type))
	  return c;

	c++;
     }
   cerr_verror ("Unable to find a type converter from %d to %d", from_type, to_type);

   return NULL;
}

static int convert_data_array (int from, void *a, int to, void **bp, size_t n)
{
   void *b = *bp;
   Converter_Type *c;

   if (NULL == (c = find_type_converter (from, to)))
     return -1;
   if (b == NULL)
     {
	if (NULL == (b = alloc_data_array (to, n)))
	  return -1;
     }
   (void) c->convert (from, a, to, b, n);
   if (*bp == NULL) *bp = b;
   return 0;
}

static int get_l1bdata_info (Swath_Id_Type swid, L1b_Info_Type *info,
			     int track, int ntracks,
			     L1b_Data_Info_Type *dinfo)
{
   int status;
   int32 i, rank, numbertype;
   char strdims[MAX_STRDIMS_LEN];
   int32 dims[MAX_DIMS];

   status = SWfieldinfo (swid, info->l1bname, &rank, dims, &numbertype, strdims);
   if (status == -1)
     {
	cerr_verror ("SWfieldinfo failed for %s data", info->l1bname);
	return -1;
     }

   if (info->num_dims != rank)
     {
	cerr_verror ("SWfieldinfo returned rank=%d, expected %d for %s",
		     rank, info->num_dims, info->l1bname);
	return -1;
     }

   if (strncmp (strdims, "nTimes", 6))
     {
	cerr_verror ("Expecting nTimes to be the slowest dimension of L1b %s",
		     info->l1bname);
	return -1;
     }

   /* Note: dims are in C order, with fastest varying last.  The first
    * value must correspond to the nTimes dimension.
    */
   if ((track < 0) || (ntracks < 0)
       || (track + ntracks > dims[0]))
     {
	cerr_verror ("track/num_tracks is invalid for L1b %s\n", info->l1bname);
	return -1;
     }

   dinfo->datatype = numbertype;
   dinfo->rank = rank;
   dinfo->start[0] = track;
   dinfo->edge[0] = ntracks;
   dinfo->total_num = ntracks;

   for (i = 1; i < rank; i++)
     {
	dinfo->start[i] = 0;
	dinfo->edge[i] = dims[i];
	dinfo->total_num *= dinfo->edge[i];	       /* FIXME! overflow possible */
     }

   return 0;
}

static int do_swreadfield (int swid, char *field,
			   int32 *start, int32 *edge,
			   void *data)
{
   if (-1 == SWreadfield (swid, field, start, NULL, edge, data))
     {
	cerr_verror ("Swreadfield for %s failed", field);
	return -1;
     }
   return 0;
}


/* If *datap != NULL, use it, otherwise malloc it */
static int read_l1b_data1 (int swid, char *fieldname,
			   int32 track, int32 ntracks,
			   int datatype, void **datap,
			   L1b_Data_Info_Type *da)
{
   L1b_Info_Type *l;
   void *from_data;
   int status;
   Converter_Type *c;
   void *data;

   if (NULL == (l = find_l1b_info (fieldname)))
     return -1;

   if (l->readfunc != NULL)
     return l->readfunc (l, swid, track, ntracks, datatype, datap, da);

   if (-1 == get_l1bdata_info (swid, l, track, ntracks, da))
     return -1;

   if (ntracks == 0)
     return 0;

   data = *datap;
   if ((data == NULL)
       && (NULL == (data = alloc_data_array (datatype, da->total_num))))
     return -1;

   if (da->datatype == datatype)
     {
	status = do_swreadfield (swid, l->l1bname, da->start, da->edge, data);
	goto free_and_return;
     }

   if ((NULL == (c = find_type_converter (da->datatype, datatype)))
       || (NULL == (from_data = alloc_data_array (datatype, da->total_num))))
     {
	status = -1;
	goto free_and_return;
     }
   status = do_swreadfield (swid, l->l1bname, da->start, da->edge, from_data);

   if (status == 0)
     status = c->convert (da->datatype, from_data, datatype, data, da->total_num);

   free ((char *)from_data);
   /* drop */

free_and_return:
   if (*datap == NULL)
     {
	if (status == -1) free (data);
	else *datap = data;
     }
   return status;
}

static int read_l1b_data (int swid, char *fieldname,
			  int32 track, int32 ntracks,
			  int datatype,
			  void *data)
{
   L1b_Data_Info_Type dinfo;

   fprintf (stdout, "read_l1b_data: swid=%d, track=%d, ntracks=%d, f=%s\n",
	    swid, track, ntracks, fieldname); (void) fflush (stdout);

   return read_l1b_data1 (swid, fieldname, track, ntracks,
			  datatype, &data, &dinfo);
}

static int compute_mant_expon (L1b_Info_Type *l, Swath_Id_Type swid,
			       int32 track, int32 ntracks,
			       int datatype, void **datap,
			       L1b_Data_Info_Type *da,
			       char *mantstr, char *exponstr)
{
   float64 *mant;
   int16 *expon;
   L1b_Data_Info_Type expon_da;
   Converter_Type *c;
   size_t i, total_num;
   void *data;
   float64 mant_fill, val_fill;
   int status;

   /* FIXME!!! */
   val_fill = -1.0 * (pow(2.0,100));
   mant_fill = -32767;

   if (NULL == (c = find_type_converter (DFNT_FLOAT64, datatype)))
     return -1;

   mant = NULL;
   expon = NULL;

   if (-1 == read_l1b_data1 (swid, mantstr, track, ntracks,
			     DFNT_FLOAT64, (void **)&mant, da))
     return -1;

   if (l->num_dims != da->rank)
     {
	free (mant);
	cerr_verror ("%s does not have the same rank as L1B %s and %s",
		     l->l1bname, mantstr, exponstr);
	return -1;
     }

   if (-1 == read_l1b_data1 (swid, exponstr, track, ntracks,
			     DFNT_INT16, (void **)&expon, &expon_da))
     {
	free (mant);
	return -1;
     }

   /* FIXME: I should compare dimension also */
   total_num = expon_da.total_num;
   if (total_num != da->total_num)
     {
	cerr_verror ("%s and %s have different sizes in the l1B file",
		     mantstr, exponstr);
	free (expon);
	free (mant);
	return -1;
     }
   for (i = 0; i < total_num; i++)
     {
	if (mant[i] == mant_fill)
	  {
	     mant[i] = val_fill;
	     continue;
	  }
	mant[i] *= pow(10.0, expon[i]);
     }
   free (expon);

   status = -1;
   data = *datap;
   if (data == NULL)
     {
	if (datatype == DFNT_FLOAT64)
	  {
	     *datap = mant;
	     return 0;
	  }

        if (NULL == (data = alloc_data_array (datatype, total_num)))
	  {
	     free (mant);
	     return -1;
	  }
     }

   status = c->convert (DFNT_FLOAT64, mant, datatype, data, total_num);
   free (mant);
   if (status == -1)
     {
	if (*datap == NULL)
	  free (data);
	return -1;
     }

   *datap = data;		       /* ok to copy onto itself */
   return 0;
}

static int compute_irradiance (L1b_Info_Type *l, Swath_Id_Type swid,
			       int32 track, int32 ntracks,
			       int datatype, void **datap,
			       L1b_Data_Info_Type *da)
{
   return compute_mant_expon (l, swid, track, ntracks, datatype, datap, da,
			      "IrradianceMantissa", "IrradianceExponent");
}

static int compute_radiance (L1b_Info_Type *l, Swath_Id_Type swid,
			     int32 track, int32 ntracks,
			     int datatype, void **datap,
			     L1b_Data_Info_Type *da)
{
   return compute_mant_expon (l, swid, track, ntracks, datatype, datap, da,
			      "RadianceMantissa", "RadianceExponent");
}

static int compute_radiance_precision (L1b_Info_Type *l, Swath_Id_Type swid,
				       int32 track, int32 ntracks,
				       int datatype, void **datap,
				       L1b_Data_Info_Type *da)
{
   return compute_mant_expon (l, swid, track, ntracks, datatype, datap, da,
			      "RadiancePrecisionMantissa", "RadianceExponent");
}

static int do_SWdiminfo (Swath_Id_Type swid, char *name, int32 *val)
{
   int32 v;

   v = SWdiminfo (swid, name);
   if (v == -1)
     {
	cerr_verror ("SWdiminfo failed for name %s", name);
	return -1;
     }
   *val = v;
   return 0;
}

static int compute_wavelength (L1b_Info_Type *l, Swath_Id_Type swid,
			       int32 track, int32 ntracks,
			       int datatype, void **datap,
			       L1b_Data_Info_Type *da)
{
   float32 *wavelengths, *wavelengths_ij, *coeffs, *coeffs_ij;
   int32 i, nwavelengths, ncoeffs, nxtrack;
   int16 *refcol;
   L1b_Data_Info_Type da_refcol;
   size_t num, total_num;

   if ((-1 == do_SWdiminfo (swid, "nWavelCoef", &ncoeffs))
       || (-1 == do_SWdiminfo (swid, "nWavel", &nwavelengths))
       || (-1 == do_SWdiminfo (swid, "nXtrack", &nxtrack)))
     return -1;

   coeffs = NULL;
   if (-1 == read_l1b_data1 (swid, "WavelengthCoefficient", track, ntracks,
			     DFNT_FLOAT32, (void **)&coeffs, da))
     return -1;
   /* The coeffs is an array of dims: [ntracks, nxtrack, ncoeffs] */
   ncoeffs = da->edge[2];
   nxtrack = da->edge[1];

   refcol = NULL;
   if (-1 == read_l1b_data1 (swid, "WavelengthReferenceColumn", track, ntracks,
			     DFNT_INT16, (void **)&refcol, &da_refcol))
     {
	free ((void *) coeffs);
	return -1;
     }
   /* refcol is a 1d array [ntracks] */

   num = (size_t)nxtrack * (size_t)ntracks;
   total_num = num * nwavelengths;

   wavelengths = (float32 *) *datap;
   if ((wavelengths == NULL) || (datatype != DFNT_FLOAT32))
     {
	wavelengths = (float32 *)alloc_data_array (DFNT_FLOAT32, total_num);
	if (wavelengths == NULL)
	  {
	     free ((void *)refcol);
	     free ((void *)coeffs);
	     return -1;
	  }
     }

   coeffs_ij = coeffs;
   wavelengths_ij = wavelengths;
   for (i = 0; i < ntracks; i++)
     {
	int32 j;
	int16 refcol_i = refcol[i];

	for (j = 0; j < nxtrack; j++)
	  {
	     int32 k;
	     for (k = 0; k < nwavelengths; k++)
	       {
		  double w, f;
		  int32 q;

		  f = k - refcol_i;

		  /* Need to evaluate a polynomial:
		   * y_n = c_0 + f*c_1 + f^2*c_2 + ... + f^{n-1}*c_{n-1}
		   *     = c_0 + f*(c_1 + f*(c_2 + ...))
		   * 
		   */
		  w = 0.0;
		  q = ncoeffs;
		  while (q > 0)
		    {
		       q--;
		       w = coeffs_ij[q] + f*w;
		    }
		  wavelengths_ij[k] = w;
	       }
	     wavelengths_ij += nwavelengths;
	     coeffs_ij += ncoeffs;
	  }
     }

   if (0)
     {
	size_t k;
	float minw = wavelengths[0], maxw = minw;
	for (k = 0; k < total_num; k++)
	  {
	     if (wavelengths[k] > maxw) maxw = wavelengths[k];
	     if (wavelengths[k] < minw) minw = wavelengths[k];
	  }
	fprintf (stdout, "Computed wavelengths span %f - %f\n",
		 (double) minw, (double) maxw);
     }

   free ((void*)coeffs);
   free ((void*)refcol);

   if (DFNT_FLOAT32 != datatype)
     {
	int status;
	void *data = *datap;
	/* data could be NULL, causing convert_data_array to allocate it,
	 * which is what is desired for that case.
	 */
	status = convert_data_array (DFNT_FLOAT32, (void*)wavelengths,
				     datatype, (void **)&data, nwavelengths);
	free ((void *)wavelengths);
	if (status == 0) *datap = data;
	return status;
     }

   *datap = (void *) wavelengths;
   return 0;
}

static int my_strcasecmp (char *a, char *b)
{
   while (1)
     {
	char cha = *a++;
	char chb = *b++;

	if ((cha == 0) || (chb == 0))
	  return (int)cha - (int)chb;

	cha |= 0x20; chb |= 0x20;
	if (cha != chb)
	  return cha - chb;
     }
}

static L1b_Info_Type L1b_Info_Table[] =
{
   {"GroundPixelQualityFlags", 2, NULL},
   {"ImageBinningFactor", 1, NULL},
   {"Irradiance", 3, compute_irradiance},
   {"IrradianceExponent", 3, NULL},
   {"IrradianceMantissa", 3, NULL},
   {"Latitude", 2, NULL},
   {"Longitude", 2, NULL},
   {"PixelQualityFlags", 3, NULL},
   {"Radiance", 3, compute_radiance},
   {"RadianceExponent", 3, NULL},
   {"RadianceMantissa", 3, NULL},
   {"RadiancePrecision", 3, compute_radiance_precision},
   {"RadiancePrecisionMantissa", 3, NULL},
   {"SolarAzimuthAngle", 2, NULL},
   {"SolarZenithAngle", 2, NULL},
   {"SpacecraftAltitude", 1, NULL},
   {"TerrainHeight", 2, NULL},
   {"Time", 1, NULL},
   {"ViewingAzimuthAngle", 2, NULL},
   {"ViewingZenithAngle", 2, NULL},
   {"Wavelength", 2, compute_wavelength},
   {"WavelengthCoefficient", 3, NULL},
   {"WavelengthReferenceColumn", 1, NULL},
   {"XTrackQualityFlags", 2, NULL},
   {NULL, 0, NULL}      /* end of table */
};

static L1b_Info_Type *find_l1b_info (char *name)
{
   L1b_Info_Type *l = L1b_Info_Table;

   while (l->l1bname != NULL)
     {
	if (0 == my_strcasecmp (l->l1bname, name))
	  return l;
	l++;
     }
   cerr_verror ("Internal error: find_l1b_info: no entry for %s", name);
   return NULL;
}

#if 0
int main (int argc, char **argv)
{
   char *file, *field, *swname;
   int fid, swid;

   file = argv[1];
   swname = argv[2];
   field = argv[3];

   fid = SWopen (file, DFACC_READ);
   if (fid == -1)
     {
	cerr_verror ("Unable to open %s", file);
	return 1;
     }
   if (-1 == (swid = SWattach (fid, swname)))
     {
	cerr_verror ("Unable to attach to %s", swname);
	(void) SWclose (fid);
	return 1;
     }

   (void) read_l1b_data (swid, field, FLOAT32_TYPE, 0, 20, NULL);
   (void) SWdetach (swid);
   (void) SWclose (fid);
   return 0;
}
#endif
/* FORTRAN bindings */

FCALLSCFUN6(INT, read_l1b_data, NEWGETL1BBLK, newgetl1bblk, INT, STRING,
            INT, INT, INT, PVOID )
