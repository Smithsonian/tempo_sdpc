#! /usr/bin/env slsh

require ("cmdopt");
require ("curl");
% require ("process");

private variable Root_Url = "https://ftp.ncep.noaa.gov/data/nccf/com";

private define usage ()
{
   variable msg =
`Usage:  dl_met.sl [args]

 Options:
 -d|--dir DIR     Root directory
 -h|--help        print this message
`;
   () = fprintf (stderr, msg);
   exit (1);
}

private define cmdopt_error (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

private define mkdir_p (path)
{
   variable mode = qualifier ("mode", 0777);

   variable st = stat_file (path);
   if (st != NULL)
     {
        if (0 != stat_is ("dir", st.st_mode))
          return 0;
        else return -1;
     }

   ifnot (is_substr (path, "/"))
     {
        if (mkdir (path, mode) != 0)
          return (errno == EEXIST) ? 0 : -1;
     }

   variable dirs = strtok (path, "/");
   if (path[0] == '/') dirs[0] = "/" + dirs[0];

   variable i, n = length(dirs);
   variable s = "";
   _for i (0, n-1, 1)
     {
        s = path_concat (s, dirs[i]);
        if (mkdir(s, mode) != 0)
          {
             if (errno != EEXIST)
               return -1;
          }
     }
   return 0;
}

private define latest_grib2_url_nam227 (when)
{
   % cycle is the model cycle runtime
   % (for NAM 227, docs say 2x per day, but files present on the
   %  web site say 4 times per day:  00,06,12,18)
   % forecast hour 00,01,02,...,60

   % grid 227 is high resolution CONUS nested grid, 3 km resolution(?)

   variable prod_root_url = "${Root_Url}/nam/prod"$;
   variable cycles = [0, 6, 12, 18];   % NAM 227 forecast cycles

   variable forecast_hour = when.hours_past_midnight;
   variable latest_cycle = cycles[wherelast_lt (cycles, forecast_hour)];

   variable dir = sprintf ("nam.%4d%02d%02d", when.year, when.month, when.day);
   variable file = sprintf ("nam.t%02dz.conusnest.hiresf%02d.tm00.grib2",
                            latest_cycle, forecast_hour);
   variable grib2_url = "$prod_root_url/$dir/$file"$;

   return grib2_url;
}

private define latest_grib2_url_rap130 (when)
{
   % RAP 130 updated hourly
   % Filenames: rap.tccz.awp130pgrbfxx.grib2
   % cc is the model cycle runtime
   % xx is the forecast hour

   % grid 130 is Lambert Conformal, 13 km resolution

   variable prod_root_url = "${Root_Url}/rap/prod"$;
   variable cycles = [0:23];   % RAP 130 forecast cycles

   variable forecast_hour = when.hours_past_midnight;
   variable latest_cycle = cycles[wherelast_lt (cycles, forecast_hour)];

   variable dir = sprintf ("rap.%4d%02d%02d", when.year, when.month, when.day);
   variable file = sprintf ("rap.t%02dz.awp130pgrbf%02d.grib2",
                            latest_cycle, forecast_hour);
   variable grib2_url = "$prod_root_url/$dir/$file"$;

   return grib2_url;
}

private define latest_grib2_url_nam221 (when)
{
   % NAM 221 updated every 6 hours
   % Filenames: nam.tccz.awip32xx.tm00.grib2
   % cc is the model cycle runtime
   % xx is the forecast hour

   % grid 221 is Lambert Conformal, 32 km resolution

   % Note:  A rapid-refresh forecast on the 221 grid is also
   % available.  Just replace nam -> rap in the URLs below.
   %

   variable prod_root_url = "${Root_Url}/nam/prod"$;
   variable cycles = [0, 6, 12, 18];   % NAM 221 forecast cycles

   variable forecast_hour = when.hours_past_midnight;
   variable latest_cycle = cycles[wherelast_lt (cycles, forecast_hour)];

   variable dir = sprintf ("nam.%4d%02d%02d", when.year, when.month, when.day);
   variable file = sprintf ("nam.t%02dz.awip32%02d.tm00.grib2",
                            latest_cycle, forecast_hour);
   variable grib2_url = "$prod_root_url/$dir/$file"$;

   return grib2_url;
}

private variable Grib2_Methods = Assoc_Type[];
Grib2_Methods["rap130"] = &latest_grib2_url_rap130;
Grib2_Methods["nam221"] = &latest_grib2_url_nam221;
Grib2_Methods["nam227"] = &latest_grib2_url_nam227;

private define new_stream (src, dest_dir)
{
   variable s = struct
     {
        src_url_method = Grib2_Methods[src],
        dest_dir = dest_dir
     };

   return s;
}

private define parse_byte_offsets (line)
{
   variable byte_offset;
   if (1 != sscanf (line, "%*f:%ld:", &byte_offset))
     throw ApplicationError, "Parsing line: $line"$;
   return byte_offset;
}

private define choose_variables (line)
{
   return (is_substr (line, "PRES:surface")
           || is_substr (line, "PRES:tropopause")
           || string_match (line, "TMP:[0-9]+ mb:"));
}

private define make_byte_range_string (index_file)
{
   variable fp = fopen (index_file, "r");
   if (NULL == fp)
     throw ReadError, "reading $index_file"$;
   variable grib2_inventory = fgetslines (fp);
   if (fclose (fp) != 0)
     throw IOError, "closing $index_file"$;

   grib2_inventory = array_map (String_Type, &strtrim, grib2_inventory, "\n");

   variable byte_offsets = array_map (Int_Type, &parse_byte_offsets,
                                      grib2_inventory);

   variable keep_variables = array_map (Int_Type, &choose_variables,
                                        grib2_inventory);
   variable var_indices = where(keep_variables);

   variable num_lines = length(grib2_inventory);
   variable num_keep = length(var_indices);
   variable k, str, ranges = {};
   _for k (0, num_keep-1, 1)
     {
        variable i = var_indices[k];
        if (i+1 < num_lines)
          str = sprintf ("%ld-%ld", byte_offsets[i], byte_offsets[i+1]);
        else
          str = string(byte_offsets[i]);
        list_append (ranges, str);
     }
   ranges = strjoin (list_to_array (ranges), ",");

   return ranges;
}

private define write_callback (fp, data)
{
   variable len = bstrlen (data);
   if (fwrite (data, fp) != len)
     return -1;
   return 0;
}

private define filter_mime_headers (infile, outfile)
{
   variable fp = fopen (infile, "r");
   if (fp == NULL)
     throw OpenError, "opening $infile for reading"$;
   variable st = stat_file (infile);

   variable num_bytes = st.st_size;
   variable bytes, bytes_read;
   bytes_read = fread_bytes (&bytes, num_bytes, fp);
   if (bytes_read < 0)
     throw ReadError;

   if (fclose (fp) != 0)
     throw IOError, "closing $infile"$;

   % let c = CRLF = 2 bytes 0D0A
   %     b = --   = 2 bytes 2D2D
   %     X = a 16-byte marker
   % file structure is:
   %    block1,block2,block3
   % where
   %  block[0] = <mime-start><mime-header><data-0>
   %  block[i] = <mime-start><mime-header><data-i>
   %  block[N] = <mime-end>
   % and where
   %  <mime-start>  = cb--X
   %  <mime-header> = c<header-info>cc
   %  <data-i>      = GRIB<bytes>7777G
   %  <mime-end>    = cb--X
   %  <header-info> is e.g. something like:
   %       Content-type: text/plain; charset=UTF-8
   %     Content-range: bytes 6195450-7543118/472540423

   fp = fopen (outfile, "w");
   if (NULL == fp)
     throw WriteError, "opening $outfile for writing"$;

   variable num_label_bytes = 20;
   variable label = bytes[[0:num_label_bytes-1]];
   %vmessage ("label %s", label);
   variable grib_beg = typecast ("GRIB", BString_Type);
   variable grib_end = typecast ("7777G\x0D\x0A", BString_Type);
   %vmessage ("grib:  %s - %s", grib_beg, grib_end);

   variable m, b, e, counter = 0;

   variable i = 1L;
   while (i < num_bytes)
     {
        m = is_substrbytes (bytes, label, i);
        if (m == 0)
          {
             %vmessage ("m = 0");
             break;
          }
        %vmessage ("label at m=$m  (%0x)"$, m);
        i = m + num_label_bytes;
        b = is_substrbytes (bytes, grib_beg, i);
        if (b == 0)
          {
             %vmessage ("b = 0");
             break;
          }
        %vmessage ("$grib_beg at b=$b (%0x)"$, b);
        i = b + 4;
        e = is_substrbytes (bytes, grib_end, i);
        if (e == 0)
          {
             %vmessage ("e = 0");
             break;
          }
        %vmessage ("$grib_end at e=$e (%0x)"$, e);
        i = e + 4;
        %vmessage ("data: %0x - %0x", b, e);
        variable len = i - b + 3;
        if (len != fwrite (bytes[b - 1 + [0:len-1]], fp))
          throw WriteError;
        counter += 1;
        %vmessage ("wrote block $counter"$);
     }
   if (fclose (fp) != 0)
     throw IOError, "closing $outfile"$;
}

private define forecast_download (grib2_url, dest_dir)
{
   if (mkdir_p (dest_dir) != 0)
     {
        throw ApplicationError, "***Error: accessing output directory $dest_dir"$;
     }

   variable grib2_url_idx = "${grib2_url}.idx"$;
   variable index_file = path_concat (dest_dir, path_basename (grib2_url_idx));

#iffalse
   variable s, argv;
   argv = ["curl", "-f", "-s", grib2_url_idx, "-o", index_file];
   s = new_process (argv; stderr="idx.log").wait();
   if (s.exit_status != 0)
     throw ApplicationError, "download failed: $grib2_url_idx"$;
#else
   variable ci = curl_new (grib2_url_idx);
   variable fpi = fopen (index_file, "wb");
   if (NULL == fpi)
     throw OpenError, "opening $index_file"$;
   curl_setopt (ci, CURLOPT_WRITEFUNCTION, &write_callback, fpi);
   curl_perform (ci);
   if (fclose (fpi) != 0)
     throw IOError, "closing $index_file"$;
#endif

   variable byte_range_string = make_byte_range_string (index_file);

   variable timestamp = strftime ("%Y%m%d%H", localtime(_time));
   variable grib2_basename =
     sprintf ("%s.%s", timestamp, path_basename (grib2_url));

   variable grib2_file = path_concat (dest_dir, grib2_basename);
   variable tmp_grib2_file = "${grib2_file}.tmp"$;
#iffalse
   argv = ["curl", "-f", "--raw", "-s",
           "-r", byte_range_string, grib2_url,
           "-o", tmp_grib2_file];
   %vmessage (strjoin(argv, " "));
   s = new_process (argv; stderr="grib2.log").wait();
   if (s.exit_status != 0)
     throw ApplicationError, "download failed: $grib2_url"$;
#else
   variable c = curl_new (grib2_url);
   variable fp = fopen (tmp_grib2_file, "wb");
   if (NULL == fp)
     throw OpenError, "opening $tmp_grib2_file"$;
   curl_setopt (c, CURLOPT_WRITEFUNCTION, &write_callback, fp);
   curl_setopt (c, CURLOPT_RANGE, byte_range_string);
   curl_perform (c);
   if (fclose (fp) != 0)
     throw IOError, "closing $tmp_grib2_file"$;
#endif

   variable tmp_grib2_file_filtered = grib2_file + ".filtered";
   filter_mime_headers (tmp_grib2_file, tmp_grib2_file_filtered);

   if (remove (tmp_grib2_file) != 0)
     throw ApplicationError, "removing $tmp_grib2_file"$;
   if (remove (index_file) != 0)
     throw ApplicationError, "removing $index_file"$;

   if (rename (tmp_grib2_file_filtered, grib2_file) != 0)
     throw ApplicationError, "renaming $tmp_grib2_file_filtered"$;

   vmessage ("%s: downloaded %s",
             strftime ("%Y%m%dT%H%M%S %Z", localtime(_time)),
             grib2_file);

   return 0;
}

define slsh_main ()
{
   variable root_dir = ".";
   variable print_usage = 0;

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("h|help", &print_usage; inc);
   opts.add ("d|dir", &root_dir; type="string");

   variable i = opts.process (__argv,1);

   if (print_usage != 0)
     {
        usage();
     }

   variable when = struct
     {
        year, month, day, hours_past_midnight
     };

   variable ts = localtime (_time);
   when.year = ts.tm_year + 1900;
   when.month = ts.tm_mon + 1;
   when.day = ts.tm_mday;
   when.hours_past_midnight = ts.tm_hour;

   variable dest_dir = path_concat (root_dir, "met");

   variable streams = {};
   %list_append (streams, new_stream ("rap130", "$dest_dir/rap130"$));
   list_append (streams, new_stream ("nam221", "$dest_dir/nam221"$));
   list_append (streams, new_stream ("nam227", "$dest_dir/nam227"$));

   variable e, s, fail_count = 0;

   try (e)
     {
        foreach s (streams)
          {
             variable url = s.src_url_method (when);
             variable status = forecast_download (url, s.dest_dir);
             % FIXME - do we care if status != 0?
             % Presumably we would just fall back to the previous forecast.
             if (status != 0) fail_count += 1;
          }
     }
   catch AnyError:
     {
        vmessage ("Caught %s, generated by %s:%d\n",
                  e.descr, e.file, e.line);
        exit (1);
     }

   exit (fail_count);
}
