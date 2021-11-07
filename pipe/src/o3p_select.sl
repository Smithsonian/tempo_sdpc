#! /usr/bin/env slsh

require ("cmdopt");

private define usage ()
{
   variable msg =
`Usage: o3p_select.sl [opts] <l1-radiance-path>
Options:
    -s|--step N      step=N means "process every Nth"
    -o|--offset M    start sequence with scan M<N"
`;
   () = fprintf (stderr, msg);
   exit (0);
}

private define cmdopt_error (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

define slsh_main ()
{
   variable step = 1;
   variable offset = 0;

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("s|step", &step; type="int");
   opts.add ("o|offset", &offset; type="int");
   variable i = opts.process (__argv,1);
   if (__argc - i < 1)
     usage();

   if (step <= 0)
     {
        () = fprintf (stdout, "yes");
        exit (0);
     }

   variable path = __argv[i];
   variable tok = strtok (path_basename_sans_extname(path), "_");

   % TEMPO_RAD_L1_V01_YYYYMMDDThhmmssZ_SsssGgg.nc
   variable scan_num, granule_num;
   variable n = sscanf (tok[5], "S%dG%d", &scan_num, &granule_num);
   if (n != 2)
     {
        throw ApplicationError, "Parsing file name: $path"$;
     }

   variable remainder = (scan_num - 1 - offset) mod step;

   if (remainder == 0)
     () = fprintf (stdout, "yes");
   else
     () = fprintf (stdout, "");
}
