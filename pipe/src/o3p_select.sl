#! /usr/bin/env slsh

require ("cmdopt");

private define usage ()
{
   variable msg =
`Usage: o3p_select.sl [opts] <l1-radiance-path>
Options:
    -l|--list STRING   Scans for which O3PROF should be generated.
                       STRING = "all" | "none" | comma-separated scan_num list
`;
   () = fprintf (stderr, msg);
   exit (0);
}

private define cmdopt_error (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

define answer (s)
{
   () = fprintf (stdout, s);
   exit(0);
}

define slsh_main ()
{
   variable scan_list;
   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("l|list", &scan_list; type="string");
   variable i = opts.process (__argv,1);
   if (__argc - i < 1)
     usage();

   scan_list = strlow (scan_list);

   if (scan_list == "all" || scan_list == "*")
     answer ("yes");
   else if (scan_list == "none" || scan_list == "")
     answer ("");

   variable path = __argv[i];
   variable tok = strtok (path_basename_sans_extname(path), "_");

   % TEMPO_RAD_L1_V01_YYYYMMDDThhmmssZ_SsssGgg.nc
   variable scan_num, granule_num;
   variable n = sscanf (tok[5], "S%dG%d", &scan_num, &granule_num);
   if (n != 2)
     {
        throw ApplicationError, "Parsing file name: $path"$;
     }

   variable list = eval(scan_list);

   if (any(list == scan_num))
     answer ("yes");
   else
     answer ("");
}
