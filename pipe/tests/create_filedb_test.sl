#! /usr/bin/env slsh

$1 = path_dirname (__FILE__);
prepend_to_slang_load_path ($1);
prepend_to_slang_load_path ("../src");
require ("pipeutil");

private variable Tempo_Epoch_Time_T = 946728000L;

define make_timet_array (num_days, tt0, num_per_day, delta_sec)
{
   variable day = delta_sec * [[0:num_per_day-1]];

   variable tts = Long_Type[num_days * num_per_day];
   variable i, iday = [[0:num_per_day-1]];

   _for i (0, num_days-1, 1)
     {
        tts[i*num_per_day + iday] = tt0 + day;
        tt0 += 86400L;
     }

   return tts;
}

define make_series_struct (tt0, num_per_day, delta_sec, path_format)
{
   variable s = struct
     {
        tt0 = tt0,
        num_per_day = num_per_day,
        delta_sec = delta_sec,
        path_format = path_format
     };
   return s;
}

private variable Path_Format;

define make_pathname (tt)
{
   variable tm = gmtime (tt);
   return tt, strftime (Path_Format, tm);
}

define make_filename_list (num_days, series)
{
   Path_Format = series.path_format;
   variable tts = make_timet_array (num_days, series.tt0, series.num_per_day, series.delta_sec);
   return array_map (Long_Type, String_Type, &make_pathname, tts);
}

define make_time_tagged_file (tt, path)
{
   variable dir = path_dirname (path);

   if (mkdir_p (dir) != 0)
     throw ApplicationError, "creating $dir"$;
   %print(tt - Tempo_Epoch_Time_T, path);
   print(tt, path);  % use time_t instead of times relative to tempo epoch
}

define slsh_main ()
{
   if (__argc != 3)
     {
        vmessage ("Usage:  %s NUM_DAYS ROOT_DIR_PATH", __argv[0]);
        exit(0);
     }

   variable num_days = eval(__argv[1]);
   variable root_dir = __argv[2];

   message ("Creating files for $num_days days in $root_dir"$);

   variable archive_root = "$root_dir/archive"$;
   variable ancillary_root = "$root_dir/ancillary"$;

   variable hour = 3600L;
   variable year = 86400L * 365L;
   variable tt0 = 1540733409L - 2 * year;

   variable drk = make_series_struct (tt0 - 2*hour, 2, 16*hour,
                                      "$archive_root/L0/drk/%Y/%m/%d/TEMPO_drk_L0_V01_%Y%m%dT%H%M%SZ.nc"$);
   variable irr = make_series_struct (tt0 - 12*hour, 1, 0*hour,
                                      "$archive_root/L1/irr/%Y/%m/%d/TEMPO_irr_L1_V01_%Y%m%dT%H%M%SZ.nc"$);
   variable snow = make_series_struct (tt0, 1, 0*hour,
                                       "$ancillary_root/snow/nsidc/%Y/%m/NISE_SSMISF18_%Y%m%d.HDFEOS"$);
   variable nam227 = make_series_struct (tt0, 12, hour,
                                         "$ancillary_root/met/nam227/%Y/%m/%d/%Y%m%d%H.nam.tffz.conusnest.hiresf%H.tm00.grib2"$);

   variable lst = {};
   list_append (lst, drk);
   list_append (lst, irr);
   list_append (lst, snow);
   list_append (lst, nam227);

   variable s, tts, paths;
   foreach s (lst)
     {
        (tts, paths) = make_filename_list (num_days, s);
        array_map (Void_Type, &make_time_tagged_file, tts, paths);
     }
}
