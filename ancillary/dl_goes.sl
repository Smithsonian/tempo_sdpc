#! /usr/bin/env slsh

%  To download GOES data, it's necessary to have an PDA service
%  account and to provide the userid and password for that account
%  in a netrc file containing a line that looks like this:
%
%   machine 140.90.190.143 login @SERVICE_USERID@ password @SERVICE_PASSWORD@
%
%  For TEMPO, the service account userid is currently AO_TEMPO_SERV1.
%  We're required to change the password every 90 days.
%
%  Note that the PDA IP address can only be accessed by local machines
%  that have been granted access through the PDA firewall.

require ("cmdopt");
require ("curl");
% require ("process");

private variable Root_Dir = ".";

private define usage ()
{
   variable msg =
`Usage: dl_goes.sl [options]
Options:
    -h|--help             Show usage message
    -d|--dir DIR          Root directory
`;
   () = fprintf (stderr, msg);
   exit (0);
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

private define log_entry (log_file, s)
{
   if (NULL == stat_file (path_dirname (log_file)))
     {
        if (0 != mkdir_p (path_dirname (log_file)))
          throw IOError, "Cannot create $log_file"$;
     }

   variable msg = sprintf ("%s: %s\n", time, s);
   variable fp = fopen (log_file, "a");
   if (fp == NULL)
     throw IOError, "appending to log file $log_file"$;
   () = fprintf (fp, msg);
   () = fclose (fp);
}

private define write_callback (fp, data)
{
   variable len = bstrlen (data);
   if (fwrite (data, fp) != len)
     return -1;
   return 0;
}

private define parse_list_entry (lst)
{
   lst = strtrim (lst, "\n");
   variable t = strtok (lst, "\\s");
   return t[-1];
}

private define read_file_list (file_list)
{
   variable fp = fopen (file_list, "r");
   if (fp == NULL)
     {
        vmessage("*** Failed to open $file_list for reading"$);
        return NULL;
     }
   variable lst = fgetslines (fp);
   () = fclose (fp);
   if (length(lst) == 0)
     {
        %vmessage("*** Empty file list: $file_list"$);
        return NULL;
     }

   return array_map (String_Type, &parse_list_entry, lst);
}

private define curl_ssl_get (url, file)
{
   variable netrc_file = "${Root_Dir}/ancillary/etc/netrc_pda"$;
   variable log_file = "${Root_Dir}/ancillary/goes/log"$;

   if (NULL == stat_file (netrc_file))
     {
        throw UsageError, "Cannot access netrc file: $netrc_file"$;
     }

   variable dest_dir = path_dirname (file);

   variable fp = fopen (file, "wb");
   if (NULL == fp)
     throw OpenError, "opening $file"$;

   variable c = curl_new (url);
   curl_setopt (c, CURLOPT_WRITEFUNCTION, &write_callback, fp);
   curl_setopt (c, CURLOPT_USE_SSL, CURLUSESSL_CONTROL);
   curl_setopt (c, CURLOPT_SSL_VERIFYPEER, 0);
   curl_setopt (c, CURLOPT_NETRC, CURL_NETRC_REQUIRED);
   curl_setopt (c, CURLOPT_NETRC_FILE, netrc_file);

   variable e;
   try (e)
     {
        curl_perform (c);
     }
   catch CurlError:
     {
        log_entry (log_file, sprintf ("Error retrieving %s: %s", url, e.message));
     }
   finally
     {
        log_entry (log_file, sprintf ("downloaded %s", file));
        return fclose (fp);
     }
}

private define pda_download (src_url, dest_dir)
{
   if (mkdir_p (dest_dir) != 0)
     throw ApplicationError, "***Error: accessing output directory $dest_dir"$;

   variable dotlis_file =
     sprintf ("%d.%s.pda.lis", getpid(), strftime ("%Y%m%d%H%M%S"));
   dotlis_file = path_concat (dest_dir, dotlis_file);

   if (0 != curl_ssl_get (src_url + "/", dotlis_file))
     return -1;

   variable lst = read_file_list (dotlis_file);
   () = remove (dotlis_file);
   if (lst == NULL)
     return -1;

   variable has_sha1_substr = array_map (Integer_Type, &is_substr, lst, ".sha1");
   variable is_nc, is_sha1 = where (has_sha1_substr, &is_nc);

   % When a GOES data product file is downloaded, both the product
   % file and the corresponding sha1 checksum file are automatically
   % deleted.  For this reason, it's important to download the sha1
   % checksum before downloading the corresponding data file.

   variable file, file_list = lst[[is_sha1, is_nc]];

   foreach file (file_list)
     {
        variable file_url = path_concat (src_url, file);
        variable dest_file = path_concat (dest_dir, file);

        if (0 != curl_ssl_get (file_url, dest_file))
          {
             vmessage ("*** Error downloading $file"$);
             continue;
          }
     }

   % FIXME: Because of PDA's auto-deletion behavior, human
   % intervention is probably required if a checksum fails.
   % We only get one chance to download the file, so there's
   % no option to automatically try again.  We'll worry about
   % this later.

   return 0;
}

define slsh_main()
{
   variable show_usage = 0;

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("h|help", &show_usage; inc);
   opts.add ("d|dir", &Root_Dir; type="string");
   variable i = opts.process (__argv,1);

   if (__argc != i || i < 0 || show_usage != 0)
     usage();

   variable pda_root_url = "ftp://140.90.190.143/PDAFileLinks/";
   variable dest_root_dir = path_concat (Root_Dir, "ancillary/goes");

   variable status;
   variable d, subdir_list = ["cmi_ch01", "cmi_ch02"];

   foreach d (subdir_list)
     {
        variable src_url = path_concat (pda_root_url, d);
        variable dest_dir = path_concat (dest_root_dir, d);
        status = pda_download (src_url, dest_dir);
        if (status != 0)
          break;
     }

   exit (status);
}
