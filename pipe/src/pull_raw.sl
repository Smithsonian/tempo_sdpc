#!/usr/bin/env slsh
require ("cmdopt");
require ("timestamp");

$1 = path_dirname (__FILE__);
prepend_to_slang_load_path (path_concat ($1, "../share/slsh/local-packages"));
prepend_to_slang_load_path ($1);
require ("pipeutil");

private variable Script_Version_String = "0.1.0";
private variable Register_With_Symlink = 0;
private variable Archive_Root_Dir = getenv ("SDPC_ARCHIVE_DIR");
private variable Registry_Subdir = "registry/incoming";

private define exit_version ()
{
   () = fprintf (stdout, "Version: %S\n", Script_Version_String);
   exit (0);
}

private define exit_usage ()
{
   variable fp = stderr;
   () = fprintf (fp, "Usage: %s [options] indir outdir\n", __argv[0]);
   variable opts =
     [
      "Options:\n",
      " -v|--version               Print version\n",
      " -h|--help                  This message\n",
      " -r|--register              Put symlink in archive $Registry_Subdir \n"$,
     ];
   foreach (opts)
     {
        variable opt = ();
        () = fputs (opt, fp);
     }
   exit (1);
}

private define make_target_path (infile, outdir)
{
   variable basename = path_basename (infile);

   % tempo_YYYYMMDDThhmmssZ_ccccc_type.raw
   % where type=grddp_XXX
   %    or type=ccsds_XXX
   variable tok = strtok (basename, "_");
   variable tt = timestamp_parse (tok[1]);
   variable tm = gmtime(tt);

   % tm_year is years since 1900, tm_yday is numbered from 0
   % Target directory is YYYY/ddd, where ddd is [001-366].
   variable target_dir = sprintf ("%s/%4d/%03d", outdir, 1900 + tm.tm_year, 1 + tm.tm_yday);

   variable timestamp = tok[1];
   variable counter_string = tok[2];
   variable type_string = strup(tok[3]);

   % For example:
   %   TEMPO_GRDDP_YYYYMMDDThhmmssZ_Cccccc.raw
   %   TEMPO_CCSDS_YYYYMMDDThhmmssZ_Cccccc.raw

   variable asdc_basename = "TEMPO_${type_string}_${timestamp}_C${counter_string}.raw"$;

   return path_concat (target_dir, asdc_basename);
}

% Return -1 on error, 0 on success, +1 when target file exists
private define pull_file (infile, outdir_root)
{
   variable target_path = make_target_path (infile, outdir_root);
   if (target_path == NULL)
     return -1;

   % Do we already have this file?  If so, we're done:
   variable st = stat_file (target_path);
   if (NULL != st)
     return 1;

   variable outdir = path_dirname (target_path);

   % If necessary, create the target directory:
   if (0 != mkdir_p (outdir))
     return -1;

   % Copy it because we may not own it.
   variable tmp = path_concat (outdir, sprintf (".%X.%X", _time(), getpid()));
   variable fp = fopen (tmp, "wb");

   if (fp == NULL)
     {
	() = fprintf (stderr, "Unable to open %S for writing: %S\n", tmp, errno_string());
	return -1;
     }
   variable fpin = fopen (infile, "rb");
   if (fpin == NULL)
     {
	() = fprintf (stderr, "Unable to open %S for reading: %S\n", infile, errno_string());
	() = fclose (fp);
	() = remove (tmp);
	return -1;
     }

   vmessage ("copying %s to %s", infile, outdir);

   variable buf;
   while (-1 != fread_bytes (&buf, 8192, fpin))
     {
	if (-1 == fwrite (buf, fp))
	  {
	     () = fprintf (stderr, "Write error: %S\n");
	     () = fclose (fpin);
	     () = fclose (fp);
	     () = remove (tmp);
	     return -1;
	  }
     }
   () = fclose (fpin);
   () = fclose (fp);

   () = rename (tmp, target_path);

   if (Register_With_Symlink != 0)
     {
        variable symlink_dir = path_concat (Archive_Root_Dir, Registry_Subdir);
        variable symlink_target = path_concat (symlink_dir, path_basename (target_path));
        if (0 != symlink (target_path, symlink_target))
          {
             () = fprintf (stderr, "*** Warning: could not create symlink: %s\n", symlink_target);
          }
     }

   return 0;
}

define slsh_main ()
{
   variable c = cmdopt_new ();

   c.add("h|help", &exit_usage);
   c.add("v|version", &exit_version);
   c.add("r|register", &Register_With_Symlink; inc);

   variable i = c.process (__argv, 1);

   if (i +2 != __argc)
     exit_usage ();

   variable indir = __argv[i];
   i++;
   variable outdir = __argv[i];
   variable status;

   if (Register_With_Symlink != 0)
     {
        if (NULL == Archive_Root_Dir)
          throw ApplicationError,
          "*** Error: Archive root directory not specified (SDPC_ARCHIVE_DIR not set)";
     }

   % The monitored directory should normally contain
   % no more than a few hundred files, at most.
   % (A few days of data at 144 files/day, or one file every 10 min)

   forever
     {
	variable files = glob ("$indir/tempo_*.raw"$);
	if (length (files) == 0)
	  {
	     sleep (1);
	     continue;
	  }

	files = files [array_sort (files)];
	foreach (files)
	  {
	     variable file = ();
	     () = pull_file (file, outdir);
	  }
     }
}
