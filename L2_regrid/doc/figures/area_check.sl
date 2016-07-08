#! /usr/bin/env isis

require ("readascii");

define read_area_check (file)
{
   variable s = struct
     {
        box_size, albers_area, exact_area, area_err
     };
   variable n = readascii (file, &s.box_size, &s.albers_area,
                           &s.exact_area, &s.area_err;
                           format = "%*g %le %le %le %*g %le",
                           comment="#");
   return s;
}

define slsh_main ()
{
   variable dir = "";
   if (__argc > 1)
     {
        dir = __argv[1];
     }

   variable bs = read_area_check (path_concat (dir, "area_check_bs.dat"));
   variable ne = read_area_check (path_concat (dir, "area_check_ne.dat"));

   variable id = plot_open ("area_check.ps/cps");
   set_frame_line_width (3);
   set_line_width(3);
   charsize(1.25);
   apj_viewport();

   yrange (1.e-8, 0.5);
   xlog;
   ylog;
   %point_style (-4);

   xlabel (latex2pg("Pixel angular size [\mu rad]"R));
   ylabel (latex2pg("\Delta A/A"R));

   plot (ne.box_size, abs(ne.area_err));
   oplot (bs.box_size, abs(bs.area_err), red);

   charsize(1.2);

   color(1);
   xylabel (100, 1.e-1,
            latex2pg("(\lambda,\phi)_{NE} = (-60.2847, +50.1414)"R), 0, -0.25);
   color(red);
   xylabel (100, 3.e-2,
            latex2pg("(\lambda,\phi)_{aim} = (-96.9657, +33.9239)"R), 0, -0.25);

   line_style ("-");
   point_style (1);
   oplot ([120, 120], [1.e-8, 1], blue);
   charsize (1);
   color(blue);
   xylabel (125, 3.e-8, latex2pg("slit width = 120 \mu rad"R));

   _pgsci(1);
   _pgsch(1.25);
   _pgmtxt ("T", 2, 0.5, 0.5, "Fractional error in pixel area");
   _pgsch(1.0);
   _pgmtxt ("T", 1, 0.5, 0.5, "Albers-projected quadrilateral vs. \"exact\"");

   plot_close(id);
}
