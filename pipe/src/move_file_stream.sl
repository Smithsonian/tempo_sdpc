#!/usr/bin/env slsh
require ("cmdopt");

private variable Script_Version_String = "0.1.0";

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
     ];
   foreach (opts)
     {
        variable opt = ();
        () = fputs (opt, fp);
     }
   exit (1);
}

private define move_file (infile, outdir)
{
   % Copy it because we may not own it.
   variable tmp = path_concat (outdir, sprintf (".%X.%X", _time(), getpid()));
   variable fp = fopen (tmp, "wb");
   variable outfile = path_concat (outdir, path_basename (infile));

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
	     () = fprintf (stderr, "Write error: %S\n", tmp);
	     () = fclose (fpin);
	     () = fclose (fp);
	     () = remove (tmp);
	     return -1;
	  }
     }
   () = fclose (fpin);

   if (-1 == fclose (fp))
     {
        () = fprintf (stderr, "Close error: %S\n", tmp);
        () = remove (tmp);
        return -1;
     }

   if (-1 == rename (tmp, outfile))
     {
        () = fprintf (stderr, "Rename error: %S -> %S\n", tmp, outfile);
        return -1;
     }

   () = remove (infile);
   return 0;
}

define slsh_main ()
{
   variable c = cmdopt_new ();

   c.add("h|help", &exit_usage);
   c.add("v|version", &exit_version);

   variable i = c.process (__argv, 1);

   if (i +2 != __argc)
     exit_usage ();

   variable indir = __argv[1];
   variable outdir = __argv[2];

   forever
     {
	variable files = glob ("$indir/tempo_*"$);
	if (length (files) == 0)
	  {
	     sleep (1);
	     continue;
	  }
	files = files [array_sort (files)];
	foreach (files)
	  {
	     variable file = ();
	     () = move_file (file, outdir);
	  }
     }
}
