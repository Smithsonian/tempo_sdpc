#! /usr/bin/env slsh

require ("csv");
require ("process");

define compute_sza (sza_str)
{
   variable argv = ["../src/gnuobjs/plan", "--Zenith", sza_str];
   variable obj = new_process (argv; write=1);
   variable result = fgetslines (obj.fp1);
   variable s = obj.wait();
   if (s.exit_status != 0) throw ApplicationError;
   () = fclose (obj.fp1);
   return eval(result[0]);
}

define slsh_main ()
{
   () = fprintf (stdout, "Checking SZA calculation...");
   variable tol = 0.02;  % degrees
   variable s = csv_readcol ("sza_tests.dat"; fields=["sza_str", "sza"], type="sd");
   variable sza = array_map (Double_Type, &compute_sza, s.sza_str);
   variable i = where(abs(sza - s.sza) > tol);
   if (length(i) > 0)
     {
        vmessage ("*** FAIL: unexpected solar zenith angle values:");
        csv_writecol (stdout, s.sza_str[i], sza[i], s.sza[i]; names=["input", "sza", "sza_expected"]);
        throw ApplicationError;
     }
   () = fprintf (stdout, "OK\n");
}
