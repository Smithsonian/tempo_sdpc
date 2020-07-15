#! /usr/bin/env slsh

%  To download NSIDC data, it's necessary to have an Earthdata
%  account, and to provide the userid and password for that account
%  in a netrc file containing a line that looks like this:
%
%  machine urs.earthdata.nasa.gov login @USER@ password @PASSWORD@
%

require ("cmdopt");
require ("glob");
require ("curl");
require ("process");

private variable Root_Dir = ".";

private define usage ()
{
   variable msg =
`Usage: dl_nise.sl [options]
Options:
    -h|--help             Show usage message
    -d|--dir DIR          Root directory
    -w|--when YYYY/MM/DD  Date (year/month/day)
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
   variable cookie_file = "${Root_Dir}/ancillary/nise/etc/nsidc.cookies"$;
   variable netrc_file = "${HOME}/.netrc"$;
   variable dest_dir = "${Root_Dir}/ancillary/nise/data"$;
   variable log_file = "${Root_Dir}/ancillary/nise/log"$;

   if (NULL == stat_file (dest_dir))
     {
        if (0 != mkdir_p (dest_dir))
          throw RunTimeError, "Destination directory does not exist: $dest_dir"$;
     }

   if (NULL == stat_file (cookie_file))
     {
	if (0 != mkdir_p (path_dirname (cookie_file)))
	  throw RunTimeError, "Cookie file directory does not exist: $cookie_file"$;
     }

   % Although curl supports resuming interrupted downloads
   % (-C option), we don't use it because NISE files are small.
   % The NISE ftp server may not support resuming anyway,
   % because using curl -C - to resume download of a file that
   % already exists (locally) causes http error 416
   % "Requested Range not satisfiable".

#iftrue
   variable argv = ["curl",
                    "-sS",     % show errors, but no progress meter
                    "--fail",  % we want non-zero exit status on http 404, etc
                    "--retry", "10",
                    "--retry-delay", "60",
                    "-b", cookie_file,
                    "-c", cookie_file,
                    "--netrc-file", netrc_file,
                    "-L", from_url,
                    "-o", path_concat (dest_dir, to_file)];

   vmessage (strjoin (argv, " "));
   variable s = new_process (argv; stdout=">>$log_file"$, dup2=1).wait();
   log_entry (log_file, "downloaded $to_file"$);
   return s.exit_status;
#else
   variable out_file = path_concat (dest_dir, to_file);
   variable fp = fopen (out_file, "wb");
   if (fp == NULL)
     throw IOError, "opening $out_file"$;
   vmessage ("downloading $from_url"$);
   variable c = curl_new (from_url);
   curl_setopt (c, CURLOPT_WRITEFUNCTION, &write_callback, fp);
   curl_setopt (c, CURLOPT_COOKIEFILE, cookie_file);
   curl_setopt (c, CURLOPT_COOKIEJAR, cookie_file);
   curl_setopt (c, CURLOPT_NETRC_FILE, netrc_file);
   % support authorization redirect:
   curl_setopt (c, CURLOPT_NETRC, CURL_NETRC_REQUIRED);
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

define nise_download_date (dir_pending, year, month, day)
{
   % Available files are:
   %  NISE_SSMISF18_yyyymmdd.1.jpg
   %  NISE_SSMISF18_yyyymmdd.2.jpg
   %  NISE_SSMISF18_yyyymmdd.HDFEOS
   %  NISE_SSMISF18_yyyymmdd.HDFEOS.xml

   variable nise_v5_url = "https://n5eil01u.ecs.nsidc.org/OTHR/NISE.005";
   variable ymd_dir = sprintf ("%4d.%02d.%02d", year, month, day);
   variable hdfeos_basename
     = sprintf ("NISE_SSMISF18_%4d%02d%02d.HDFEOS", year, month, day);

   % Don't re-try files already marked as pending
   if (NULL != stat_file (path_concat (dir_pending, hdfeos_basename)))
     {
        vmessage ("skipping $hdfeos_basename"$);
        return 0;
     }

   variable url = strjoin ([nise_v5_url, ymd_dir, hdfeos_basename], "/");
   variable status = download_file (url, hdfeos_basename);

   if (status != 0)
     {
	vmessage ("$hdfeos_basename download pending"$);
	record_pending (dir_pending, url, hdfeos_basename);
     }
   else
     {
	vmessage ("$hdfeos_basename downloaded"$);
     }

   % On "connection refused", curl will exit with status=7
   % On "404 not found", curl --fail will exit with status=22
   % Caller should retry later.

   return status;
}

private define nise_download_pending (dir_pending)
{
   variable file_list = glob ("${dir_pending}/NISE_*.HDFEOS"$);
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

   dir_pending = "${Root_Dir}/ancillary/nise/pending"$;

   if (0 != mkdir_p (dir_pending))
     throw UsageError, "Directory $dir_pending does not exist"$;
   nise_download_pending (dir_pending);

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
        if (3 != sscanf (year_month_day, "%d/%d/%d", &year, &month, &day))
          usage();
        vmessage ("downloading data for $year/$month/$day"$);
     }

   variable status = nise_download_date (dir_pending, year, month, day);
   exit(status);
}
