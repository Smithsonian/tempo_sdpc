#! /usr/bin/env slsh

require ("cmdopt");
require ("glob");
require ("curl");
require ("process");

private variable Root_Dir = ".";

private define usage ()
{
   variable msg =
`Usage: dl_snow_ims.sl [options]
Options:
    -h|--help             Show usage message
    -d|--dir DIR          Root directory
    -w|--when YYYY-MM-DD  Date (year-month-day)
`;
   () = fprintf (stderr, msg);
   exit (0);
}

private define cmdopt_error (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

private define log_entry (log_file, basename)
{
   variable msg = sprintf ("%s: %s\n", time, basename);
   variable fp = fopen (log_file, "a");
   if (fp == NULL)
     throw IOError, "appending to log file $log_file"$;
   () = fprintf (fp, msg);
   () = fclose (fp);
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

private define write_callback (fp, data)
{
   variable len = bstrlen (data);
   if (len != fwrite (data, fp))
     return -1;
   return 0;
}

define download_file (from_url, to_file)
{
   variable dest_dir = "${Root_Dir}/ims/data"$;
   variable log_file = "${Root_Dir}/ims/log"$;

   if (NULL == stat_file (dest_dir))
     {
        if (0 != mkdir_p (dest_dir))
          throw RunTimeError, "Destination directory does not exist: $dest_dir"$;
     }

   % Although curl supports resuming interrupted downloads
   % (-C option), we don't use it because IMS files are small,
   % at least when compressed.
   % The NSIDC ftp server may not support resuming anyway,
   % because using curl -C - to resume download of a file that
   % already exists (locally) causes http error 416
   % "Requested Range not satisfiable".

   variable out_file = path_concat (dest_dir, to_file);
   variable tmp_file = out_file + ".tmp";

#iftrue
   variable argv = ["curl",
                    "-sS",     % show errors, but no progress meter
                    "--fail",  % we want non-zero exit status on http 404, etc
                    "--retry", "10",
                    "--retry-delay", "60",
                    "-L", from_url,
                    "-o", tmp_file];

   vmessage (strjoin (argv, " "));
   variable s = new_process (argv; stdout=">>$log_file"$, dup2=1).wait();
   if (NULL != stat_file (tmp_file))
     {
        if (rename (tmp_file, out_file) != 0)
          throw IOError, "renaming $tmp_file"$;
        log_entry (log_file, "downloaded $to_file"$);
     }
   return s.exit_status;
#else
   variable fp = fopen (out_file, "wb");
   if (fp == NULL)
     throw IOError, "opening $out_file"$;
   vmessage ("downloading $from_url"$);
   variable c = curl_new (from_url);
   curl_setopt (c, CURLOPT_WRITEFUNCTION, &write_callback, fp);
   curl_setopt (c, CURLOPT_FOLLOWLOCATION, 1L);
   curl_setopt (c, CURLOPT_UNRESTRICTED_AUTH, 1L);
   % curl_setopt (c, CURLOPT_FAILONERROR, 1L);  % s-lang module
   % doesn't support this... and that's kind of a show-stopper
   % because '404 not found" is reasonably common.
   % The file dated today often appears late. Maybe we're just
   % looking too early in the day?
   variable e, status = 0;
   try (e)
     {
        curl_perform (c);
     }
   catch CurlError:
     {
        log_entry (log_file,
                   sprintf ("Unable to retrieve %s: %s", from_url, e.message));
        status = 0;
     }
   finally
     {
        () = fclose (fp);
        log_entry (log_file, "downloaded $out_file"$);
     }

   return status;
#endif
}

define record_pending (dir_pending, url, file)
{
   variable pending_file = path_concat (dir_pending, file);
   variable fp = fopen (pending_file, "w");
   if (fprintf (fp, "%s\n", url) < 0)
     throw WriteError, "*** writing $pending_file"$;
   () = fclose (fp);
}

define ims_download_date (dir_pending, year, month, day)
{
   variable ims_v1_url = "ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02156/GIS/1km";
   variable year_dir = sprintf ("%4d", year);

   % Filename contains day of year.  Sigh.
   variable date_string = sprintf ("%4d%02d%02d", year, month, day);
   variable argv = ["date", "-d", date_string, "+%j"];
   variable obj = new_process(argv; write=1);
   variable day_of_year = atoi(fgetslines (obj.fp1)[0]);
   variable s = obj.wait();
   if (s.exit_status != 0)
     {
        throw ApplicationError, "*** Error: computing day of year using: %s", strjoin (argv, " ");
     }
   () = fclose (obj.fp1);
   variable tif_gz_basename = sprintf ("ims%4d%03d_1km_GIS_v1.3.tif.gz", year, day_of_year);

   % Don't re-try files already marked as pending
   if (NULL != stat_file (path_concat (dir_pending, tif_gz_basename)))
     {
        vmessage ("skipping $tif_gz_basename"$);
        return 0;
     }

   variable url = strjoin ([ims_v1_url, year_dir, tif_gz_basename], "/");
   variable status = download_file (url, tif_gz_basename);

   if (status != 0)
     {
	vmessage ("$tif_gz_basename download pending"$);
	record_pending (dir_pending, url, tif_gz_basename);
     }
   else
     {
	vmessage ("$tif_gz_basename downloaded"$);
     }

   % On "connection refused", curl will exit with status=7
   % On "404 not found", curl --fail will exit with status=22
   % Caller should retry later.

   return status;
}

private define ims_download_pending (dir_pending)
{
   variable file_list = glob ("${dir_pending}/ims*tif.gz"$);
   if (length(file_list) == 0)
     return;

   variable file;
   foreach file (file_list)
     {
        variable fp = fopen (file, "r");
        if (fp == NULL)
          throw ReadError, "*** reading $file"$;
        variable url = fgetslines (fp);
        () = fclose (fp);
        url = strtrim (url[0], "\n");
        %vmessage ("try pending url=$url"$);
        variable status = download_file (url, path_basename(file));
        if (status == 0)
          {
             () = remove (file);
          }
     }
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

define slsh_main()
{
   variable show_usage = 0;
   variable dir_pending = NULL;
   variable year_month_day = NULL;

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("h|help", &show_usage; inc);
   opts.add ("d|dir", &Root_Dir; type="string");
   opts.add ("w|when", &year_month_day; type="string");
   variable i = opts.process (__argv,1);

   if (__argc != i || i < 0 || show_usage != 0)
     usage();

   if (0 != mkdir_p (Root_Dir))
     throw UsageError, "Directory $Root_Dir does not exist"$;

   dir_pending = "${Root_Dir}/ims/pending"$;

   if (0 != mkdir_p (dir_pending))
     throw UsageError, "Directory $dir_pending does not exist"$;
   ims_download_pending (dir_pending);

   variable year, month, day;
   if (year_month_day == NULL)
     {
        variable ts = localtime (_time);
        year = ts.tm_year + 1900;
        month = ts.tm_mon + 1;
        day = ts.tm_mday;
     }
   else
     {
        if (3 != sscanf (year_month_day, "%d-%d-%d", &year, &month, &day))
          usage();
        vmessage ("downloading data for %d-%02d-%02d", year, month, day);
     }

   variable status = ims_download_date (dir_pending, year, month, day);
   exit(status);
}
