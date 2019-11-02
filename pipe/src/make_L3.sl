#! /usr/bin/env slsh

require ("csv");
require ("cmdopt");
require ("glob");
require ("process");
require ("setfuns");

_debug_info=1;
_traceback=1;

$1 = path_dirname (__FILE__);
prepend_to_slang_load_path (path_concat ($1, "../share/slsh/local-packages"));
prepend_to_slang_load_path ($1);
require ("pipeutil");

private variable Clobber_Output_Files = 0;

private define config_file_case (cases)
{
   return struct
     {
        config_file = cases[0],
        product_list = list_to_array (cases[1])
     };
}

private define applicable_config_file (cfg_case, product_list)
{
   return length (intersection (cfg_case.product_list, product_list));
}

private variable Config_File_Cases =
  array_map (Struct_Type, &config_file_case,
             [{"${SDPC_ROOT}/etc/l3.cfg"$,     {"HCHO", "NO2", "O3T"}},
              {"${SDPC_ROOT}/etc/l3_cldrr.cfg"$, {"CLDRR"}},
              {"${SDPC_ROOT}/etc/l3_o3p.cfg"$, {"O3P"}}
             ]);

private define usage ()
{
   variable argv0 = __argv[0];
   variable s =
`Usage:  $argv0 [options] -p prod1,prod2,prod3 DIR1 [DIR2 DIR3 ...]
   Options:
         -a|--archive_root DIR
         -c|--clobber        If present, will overwrite pre-existing files
         -h|--help
         -p|--products p1,p2,p3

    At run-time, environment variable SDPC_ROOT must be set.
`$;
   vmessage (s);
   exit(0);
}

private define error_routine (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

define parse_string_option (str, pval)
{
   @pval=strtok (str, ",");
}

private define assert_nonexistent (path)
{
   if (Clobber_Output_Files)
     return;
   if (NULL != stat_file (path))
     throw ApplicationError, "*** Error: exists $path"$;
}

private define is_directory (file)
{
   variable st = stat_file (file);
   if (st == NULL) return 0;
   return stat_is ("dir", st.st_mode);
}

private define set_eval_struct_field (g, name, val)
{
   switch (name)
     {
      case "product_type":
        set_struct_field (g, name, val);
     }
     {
        % default:
        set_struct_field (g, name, eval(val));
     }
}

private define read_ident_file (csv_ident_file)
{
   variable _g = csv_readcol (csv_ident_file);
   variable g = @Struct_Type (_g.col1);
   array_map (Void_Type, &set_eval_struct_field, g, _g.col1, _g.col2);
   return g;
}

private define get_ident_time_structs (ident)
{
   variable tm_start = struct
     {
        tm_year = ident.tstart_year - 1900,
        tm_mon = ident.tstart_month - 1,
        tm_mday = ident.tstart_mday,
        tm_hour = ident.tstart_hour,
        tm_min = ident.tstart_min,
        tm_sec = ident.tstart_sec,
        tm_wday = ident.tstart_wday,
        tm_yday = ident.tstart_yday,
        tm_isdst = 0
     };
   variable tm_end = struct
     {
        tm_year = ident.tend_year - 1900,
        tm_mon = ident.tend_month - 1,
        tm_mday = ident.tend_mday,
        tm_hour = ident.tend_hour,
        tm_min = ident.tend_min,
        tm_sec = ident.tend_sec,
        tm_wday = ident.tend_wday,
        tm_yday = ident.tend_yday,
        tm_isdst = 0
     };

   return tm_start, tm_end;
}

private define tstart_since_epoch (tm)
{
   return tm.time_coverage_start_since_epoch;
}

private define make_l3_filename_format (idents)
{
   variable tstart = array_map (Double_Type, &tstart_since_epoch, idents);
   variable i = array_sort (tstart);

   variable
     beg = i[0],
     end = i[-1];

   variable tbeg;
   (tbeg, ) = get_ident_time_structs (idents[beg]);
   variable tstart_str = strftime ("%Y%m%dT%H%M%SZ", tbeg);

   variable
     scan_num = idents[beg].scan_num,
     processing_version = idents[beg].processing_version;

   variable filename_format =
     sprintf ("TEMPO_%%s_L3_V%02d_%s_S%03d.nc",
              processing_version, tstart_str, scan_num);

   return filename_format;
}

private define scan_subdir (g)
{
   % Derive the destination archive directory from the contents of
   % the granule_ident CSV file.
   variable subdir_seq = [g.processing_version,
                          g.sat_local_day_start,
                          g.scan_num];
   subdir_seq = array_map (String_Type, &string, subdir_seq);
   return strjoin (subdir_seq, "/");
}

private define write_filename_list (output_dir, l3_output_file,
                                    prod, file_list)
{
   variable path = path_concat (output_dir, "l2_${prod}.lst"$);
   assert_nonexistent (path);

   variable fp = fopen (path, "w");
   if (fp == NULL)
     throw IOError, "*** Error: opening $path for writing"$;

   () = fprintf (fp, "%s\n", l3_output_file);
   variable i, n = length(file_list);
   _for i (0, n-1, 1)
     {
        () = fprintf (fp, "%s\n", file_list[i]);
     }

   if (0 != fclose (fp))
     {
        throw IOError, sprintf ("*** Error: closing ${path} (%s)",
                                errno_string(errno));
     }
}

private define perform_regridding (output_dir, cfg_path)
{
   variable regrid_exec = "${SDPC_ROOT}/bin/L2_regrid"$;
   variable s = new_process ([regrid_exec, cfg_path]; dir=output_dir).wait();
   if (s.exit_status != 0)
     throw ApplicationError, "*** Error: L2_regrid failed";
}

define process_scan_granules (scan_dir, archive_root_dir, products)
{
   variable file, dir, prod, prod_files = Assoc_Type[];
   variable granule_dir_list = scan_dir + "/" + listdir (scan_dir);
   foreach prod (products)
     {
        variable lst = {};
        foreach dir (granule_dir_list)
          {
             variable st = stat_file (dir);
             if (0 == stat_is ("dir", st.st_mode))
               continue;
             file = glob ("${dir}/${prod}/TEMPO_${prod}_L2_*.nc"$);
             if (length(file) == 0) continue;
             list_append (lst, file[0]);
          }
        if (length(lst) == 0) continue;
        lst = list_to_array (lst);
        % Sorting the filenames is equivalent to time-ordering
        % because the basenames have an embedded time stamp.
        % L2_regrid doesn't need them time-ordered, but sorting
        % the list seems like a good idea.
        variable basenames = array_map (String_Type, &path_basename, lst);
        variable i = array_sort (basenames);
        prod_files[prod] = lst[i];
     }

   variable ident_files = array_map (String_Type, &path_concat,
                                     granule_dir_list, "granule_ident.csv");
   variable idents = array_map (Struct_Type, &read_ident_file, ident_files);
   variable l3_outfile_fmt = make_l3_filename_format (idents);

   variable level3_root_dir = path_concat (archive_root_dir, "L3");
   variable output_dir = path_concat (level3_root_dir,
                                      scan_subdir(idents[0]));

   if (0 != mkdir_p (output_dir))
     {
        throw ApplicationError, "*** Error creating directory $output_dir"$;
     }

   foreach prod (products)
     {
        if (0 == assoc_key_exists (prod_files, prod))
          continue;
        variable l3_output_file = sprintf (l3_outfile_fmt, prod);
        assert_nonexistent (l3_output_file);
        write_filename_list (output_dir, l3_output_file,
                             prod, prod_files[prod]);
     }

   variable cfg_case;
   foreach cfg_case (Config_File_Cases)
     {
        if (applicable_config_file (cfg_case, products) != 0)
          {
             perform_regridding (output_dir, cfg_case.config_file);
          }
     }
}

define slsh_main()
{
   variable archive_root_dir = getenv ("SDPC_ARCHIVE_DIR");
   variable sdpc_root_dir = getenv ("SDPC_ROOT");
   variable products = NULL;

   variable c = cmdopt_new (&error_routine);
   c.add ("h|help", &usage);
   c.add ("a|archive_root_dir", &archive_root_dir; type="string");
   c.add ("c|clobber", &Clobber_Output_Files; inc);
   c.add ("p|products", &parse_string_option, &products; type="string");
   variable __i = c.process (__argv, 1);

   if (__argc - __i < 1)
     usage();

   if (products == NULL)
     throw ApplicationError, "*** Error: Level 2 product list not specified";

   if (archive_root_dir == NULL)
     throw ApplicationError,
     "*** Error: Archive root directory not specified (SDPC_ARCHIVE_DIR not set)";

   if (sdpc_root_dir == NULL)
     throw ApplicationError, "*** Error: SDPC_ROOT is not set";

   variable scan_dir, scan_dir_list = __argv[[__i:]];

   foreach scan_dir (scan_dir_list)
     {
        process_scan_granules (scan_dir, archive_root_dir, products);
     }
}
