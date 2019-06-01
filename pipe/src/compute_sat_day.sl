#! /usr/bin/env slsh

$1 = path_dirname (__FILE__);
prepend_to_slang_load_path (path_concat ($1, "../share/slsh/local-packages"));
prepend_to_slang_load_path ($1);
require ("pipeutil");

define slsh_main ()
{
   variable sdpc_root_dir = getenv ("SDPC_ROOT");
   if (sdpc_root_dir == NULL)
     throw ApplicationError, "*** Error: SDPC_ROOT is not set";

   if (__argc != 2)
     {
        vmessage ("Usage: %s SEC", __argv[0]);
        exit(1);
     }

   % taix = TAI seconds since the TEMPO epoch
   variable taix = eval(__argv[1]);
   variable sc_timezone = read_sc_timezone (sdpc_root_dir);
   variable sat_day = satellite_day (taix, sc_timezone);

   () = fprintf (stdout, "%d\n", int(sat_day));
}
