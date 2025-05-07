/*
 * The QF suite of functions unpacks the Quality Flag values in the L1B
 * OMI data and returns the unpacked values in the form of arrays of
 * 4-byte intgers.  The unpacker functions that work on a single pixel
 * (PixelQF, GroundQF, and MeasurementQF) return all the flags in a single
 * array.  The functions that unpack multiple pixel flags (MultiPixelQF,
 * MultiGroundQF, and MultiMeasurementQF) return the different flags in
 * separate integer arrays.  All of the functions assume that space for the
 * unpacked flags has been allocated by the calling program.  Calling
 * programs should also note that the unpackers all skip any flags reserved
 * for later use.  That is, those positions in the packed flags are not
 * unpacked, nor are they returned to the calling program.
 *
 * The GroundQF function unpacks a single GroundPixelQualityFlags element.
 * The PixelQF function unpacks a single PixelQualityFlags element.
 * The MeasurementQF function unpacks a single MeasurementQualityFlags element.
 * The MultGroundQF function unpacks multiple GroundPixelQualityFlags elements.
 * the MultPixelQF function unpacks multiple PixelQualityFlags elements.
 * The MultMeasurementQF function unpacks multiple MeasurementQualityFlags
 * elements.
 *
 * Consequently:
 *
 * The GroundQF function requires an allocated array of 
 * 6 4-byte integers for the unpacked flags.
 *
 * The PixelQF function requires an allocated array of 13 4-byte integers
 * for the unpacked flags.
 *
 * The MeasurementQF function requires an allocated array of 13 4-byte
 * integers for the unpacked flags.
 *
 * All of the QF functions always have a return code of 0.  There is no
 * error detection in these functions, because there are no errors to detect.
 * Also, no messages are output from these functions under any circumstances
 * at any time.
 *
 */

typedef struct {
   unsigned int surface:4;
   unsigned int glint:1;
   unsigned int eclipse:1;
   unsigned int gef:1;
   unsigned int reserve:1;
   unsigned int snowice:7;
   unsigned int nise:1;
} GPQF;

typedef struct {
   unsigned int missing:1;
   unsigned int bad_pix:1;
   unsigned int proc_err:1;
   unsigned int trans_warn:1;
   unsigned int rts_warn:1;
   unsigned int saturation_warn:1;
   unsigned int noise_warn:1;
   unsigned int dark_warn:1;
   unsigned int offset_warn:1;
   unsigned int exposure_warn:1;
   unsigned int stray_light_warn:1;
   unsigned int reserve:3;
   unsigned int dead_pix_id:1;
   unsigned int dead_pix_err:1;
} PQF;

typedef struct {
   unsigned int test_mode:1;
   unsigned int alt_eng_data:1;
   unsigned int alt_seq_read:1;
   unsigned int coadd_err:1;
   unsigned int invalid_coadd_period:1;
   unsigned int coadd_overflow_poss:1;
   unsigned int meas_combo:1;
   unsigned int rebinning:1;
   unsigned int dark_current:1;
   unsigned int detector_smear_calc:1;
   unsigned int saa_poss:1;
   unsigned int sc_manoeuver:1;
   unsigned int geo_err:1;
   unsigned int reserve:3;
} MQF;

int pixelqf_(short *packed_flags, int *unpacked_flags) {
   PQF *pix_test;

   pix_test = (PQF *)packed_flags;
   unpacked_flags[0]  = (int) pix_test->missing;
   unpacked_flags[1]  = (int) pix_test->bad_pix;
   unpacked_flags[2]  = (int) pix_test->proc_err;
   unpacked_flags[3]  = (int) pix_test->trans_warn;
   unpacked_flags[4]  = (int) pix_test->saturation_warn;
   unpacked_flags[5]  = (int) pix_test->noise_warn;
   unpacked_flags[6]  = (int) pix_test->dark_warn;
   unpacked_flags[7]  = (int) pix_test->offset_warn;
   unpacked_flags[8]  = (int) pix_test->exposure_warn;
   unpacked_flags[9]  = (int) pix_test->stray_light_warn;
   unpacked_flags[10] = (int) pix_test->dead_pix_id;
   unpacked_flags[11] = (int) pix_test->dead_pix_err;
   return 0;

}

int groundqf_(short *packed_flags, int *unpacked_flags) {
   GPQF *pix_test;

   pix_test = (GPQF *)packed_flags;
   unpacked_flags[0]  = (int) pix_test->surface;
   unpacked_flags[1]  = (int) pix_test->glint;
   unpacked_flags[2]  = (int) pix_test->eclipse;
   unpacked_flags[3]  = (int) pix_test->gef;
   unpacked_flags[4]  = (int) pix_test->snowice;
   unpacked_flags[5]  = (int) pix_test->nise;
   return 0;

}

int measurementqf_(short *packed_flags, int *unpacked_flags) {
   MQF *pix_test;

   pix_test = (MQF *)packed_flags;
   unpacked_flags[0]   = (int) pix_test->test_mode;
   unpacked_flags[1]   = (int) pix_test->alt_eng_data;
   unpacked_flags[2]   = (int) pix_test->alt_seq_read;
   unpacked_flags[3]   = (int) pix_test->coadd_err;
   unpacked_flags[4]   = (int) pix_test->invalid_coadd_period;
   unpacked_flags[5]   = (int) pix_test->coadd_overflow_poss;
   unpacked_flags[6]   = (int) pix_test->meas_combo;
   unpacked_flags[7]   = (int) pix_test->rebinning;
   unpacked_flags[8]   = (int) pix_test->dark_current;
   unpacked_flags[9]   = (int) pix_test->detector_smear_calc;
   unpacked_flags[10]  = (int) pix_test->saa_poss;
   unpacked_flags[11]  = (int) pix_test->sc_manoeuver;
   unpacked_flags[12]  = (int) pix_test->geo_err;
   return 0;

}

int multpixelqf_(short *packed_flags, int *missing, int *bad_pix, int *proc_err,
      int *trans_warn, int *saturation_warn, int *noise_warn, int *dark_warn,
      int *offset_warn, int *exposure_warn, int *stray_warn, int *dead_pix_id,
      int *dead_pix_err, int *num_flags) {
   PQF *pix_test;
   int  i;

   pix_test = (PQF *)packed_flags;
   for (i=0; i<*num_flags; i++) {
      missing[i]         = (int) pix_test->missing;
      bad_pix[i]         = (int) pix_test->bad_pix;
      proc_err[i]        = (int) pix_test->proc_err;
      trans_warn[i]      = (int) pix_test->trans_warn;
      saturation_warn[i] = (int) pix_test->saturation_warn;
      noise_warn[i]      = (int) pix_test->noise_warn;
      dark_warn[i]       = (int) pix_test->dark_warn;
      offset_warn[i]     = (int) pix_test->offset_warn;
      exposure_warn[i]   = (int) pix_test->exposure_warn;
      stray_warn[i]      = (int) pix_test->stray_light_warn;
      dead_pix_id[i]     = (int) pix_test->dead_pix_id;
      dead_pix_err[i]    = (int) pix_test->dead_pix_err;
      pix_test++;
   }
   return 0;

}

int multgroundqf_(short *packed_flags, int *surface, int *glint, int *eclipse,
      int *geoloc_err, int *snowice, int *nise, int *num_flags) {
   GPQF *pix_test;
   int   i;

   pix_test = (GPQF *)packed_flags;
   for (i=0; i<*num_flags; i++) {
      surface[i]    = (int) pix_test->surface;
      glint[i]      = (int) pix_test->glint;
      eclipse[i]    = (int) pix_test->eclipse;
      geoloc_err[i] = (int) pix_test->gef;
      snowice[i]    = (int) pix_test->snowice;
      nise[i]       = (int) pix_test->nise;
      pix_test++;
   }
   return 0;

}

int multmeasurementqf_(short *packed_flags, int *test_mode, int *alt_eng_data,
      int *alt_seq_read, int *coadd_err, int *invalid_coadd_period,
      int *coadd_overflow_poss, int *meas_combo, int *rebinning,
      int *dark_current, int *detector_smear_calc, int *saa_poss,
      int *sc_manoeuver, int *geo_err, int *num_flags) {
   MQF *pix_test;
   int  i;

   pix_test = (MQF *)packed_flags;
   for (i=0; i<*num_flags; i++) {
      test_mode[i]            = (int) pix_test->test_mode;
      alt_eng_data[i]         = (int) pix_test->alt_eng_data;
      alt_seq_read[i]         = (int) pix_test->alt_seq_read;
      coadd_err[i]            = (int) pix_test->coadd_err;
      invalid_coadd_period[i] = (int) pix_test->invalid_coadd_period;
      coadd_overflow_poss[i]  = (int) pix_test->coadd_overflow_poss;
      meas_combo[i]           = (int) pix_test->meas_combo;
      rebinning[i]            = (int) pix_test->rebinning;
      dark_current[i]         = (int) pix_test->dark_current;
      detector_smear_calc[i]  = (int) pix_test->detector_smear_calc;
      saa_poss[i]             = (int) pix_test->saa_poss;
      sc_manoeuver[i]         = (int) pix_test->sc_manoeuver;
      geo_err[i]              = (int) pix_test->geo_err;
      pix_test++;
   }
   return 0;

}

#if 0

#if defined(DEC_ALPHA) || defined(IRIX) || defined(UNICOS)

#define INT32  INT
#define INT32V INTV
#define PINT32 PINT

#else

#define INT32  LONG
#define INT32V LONGV
#define PINT32 PLONG

#endif

/*
FCALLSCFUN9(INT, OMI_Get_L1B_Data_Block, GETL1BBLK, getl1bblk, STRING, STRING,
            STRING, INT, INT, INT, PINT, INTV, PVOID )
*/

/* Fortran bindings */

FCALLSCFUN9(INT, GroundQF, GROUNDQF, groundqf, PINT, PINT)
FCALLSCFUN9(INT, PixelQF, PINT, PINT)
FCALLSCFUN9(INT, MeasurementQF, PINT, PINT)
FCALLSCFUN9(INT, MultGroundQF, PINT, PINT, PINT, PINT, PINT, PINT, PINT, PINT)
FCALLSCFUN9(INT, MultPixelQF, PINT, PINT, PINT, PINT, PINT, PINT, PINT, PINT,
      PINT, PINT, PINT, PINT, PINT, PINT)
FCALLSCFUN9(INT, MultMeasurementQF, PINT, PINT, PINT, PINT, PINT, PINT, PINT,
      PINT, PINT, PINT, PINT, PINT, PINT, PINT, PINT)

#endif
