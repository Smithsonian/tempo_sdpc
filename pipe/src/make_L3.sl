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
             [{"${SDPC_ROOT}/etc/l3.cfg"$,     {"HCHO", "NO2", "O3TOT"}},
              {"${SDPC_ROOT}/etc/l3_cldrr.cfg"$, {"CLDRR"}},
              {"${SDPC_ROOT}/etc/l3_o3p.cfg"$, {"O3PROF"}}
             ]);

private define usage ()
{
   variable argv0 = __argv[0];
   variable s =
`Usage:  $argv0 [options] SCAN_DIR1 [SCAN_DIR2 SCAN_DIR3 ...]
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

private define write_filename_list (output_dir, l3_output_file,
                                    prod, file_list)
{
   variable path = path_concat (output_dir, "TEMPO_${prod}_L2.lis"$);

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

private define make_l3_filename (product_files)
{
   variable name = path_basename (product_files[0]);
   variable tok = strchopr (name, 'G', 0);
   return strreplace (tok[1], "_L2_", "_L3_") + ".nc";
}

define insert_fixed_metadata (path)
{
   variable argv = ["insert_fixed_metadata.py", path];
   variable s = new_process (argv; dup2=1).wait();
   if (s.exit_status != 0)
     throw ApplicationError, "*** Error: inserting fixed metadata: $path"$;
}

define process_scan_granules (scan_dir, archive_root_dir, products)
{
   if (products == NULL)
     return;
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
             file = glob ("${dir}/${prod}/TEMPO_${prod}_L2_*_S???G??.nc"$);
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

   variable output_dir = strtrans (scan_dir, "/L2/", "/L3/");

   if (0 != mkdir_p (output_dir))
     {
        throw ApplicationError, "*** Error creating directory $output_dir"$;
     }

   variable l3_output_file;

   foreach prod (products)
     {
        if (0 == assoc_key_exists (prod_files, prod))
          continue;
        l3_output_file = make_l3_filename (prod_files[prod]);
        assert_nonexistent (path_concat (output_dir, l3_output_file));
        write_filename_list (output_dir, l3_output_file, prod, prod_files[prod]);
     }

   variable cfg_case;
   foreach cfg_case (Config_File_Cases)
     {
        if (applicable_config_file (cfg_case, products) != 0)
          {
             perform_regridding (output_dir, cfg_case.config_file);
          }
     }

   variable register_dir = "$SDPC_ARCHIVE_DIR/registry/incoming"$;
   if (is_directory (register_dir))
     {
        % Register L3 products produced
        foreach prod (products)
          {
             if (0 == assoc_key_exists (prod_files, prod))
               continue;

             l3_output_file = make_l3_filename (prod_files[prod]);

             variable from_path = path_concat (output_dir, l3_output_file);
             if (NULL == stat_file (from_path))
               continue;

             insert_fixed_metadata (from_path);

             variable to_path = path_concat (register_dir, l3_output_file);
             if (NULL == stat_file (to_path))
               {
                  if (0 != symlink (from_path, to_path))
                    {
                       throw IOError, "*** Error creating symbolic link: $from_path to $to_path"$;
                    }
               }
          }
     }
}

define find_product_dirs (scan_dir)
{
   variable granule_dir_list = scan_dir + "/" + listdir (scan_dir);
   variable granule_dir;

   variable lst = {};

   foreach granule_dir (granule_dir_list)
     {
        variable st = stat_file (granule_dir);
        variable f;
        if (0 == stat_is ("dir", st.st_mode))
          continue;
        variable files = glob ("${granule_dir}/*/TEMPO_*_L2_*_S???G??.nc"$);
        if (length(files) == 0)
          continue;
        foreach f (files)
          {
             list_append (lst, path_basename (path_dirname (f)));
          }
     }

   if (length(lst) == 0)
     return NULL;

   variable lst_arr = list_to_array (lst);

   return lst_arr[unique(lst_arr)];
}

define slsh_main()
{
   variable archive_root_dir = getenv ("SDPC_ARCHIVE_DIR");
   variable sdpc_root_dir = getenv ("SDPC_ROOT");
   variable cmdline_products = NULL;

   variable c = cmdopt_new (&error_routine);
   c.add ("h|help", &usage);
   c.add ("a|archive_root_dir", &archive_root_dir; type="string");
   c.add ("c|clobber", &Clobber_Output_Files; inc);
   c.add ("p|products", &parse_string_option, &cmdline_products; type="string");
   variable __i = c.process (__argv, 1);

   if (__argc - __i < 1)
     usage();

   if (archive_root_dir == NULL)
     throw ApplicationError,
     "*** Error: Archive root directory not specified (SDPC_ARCHIVE_DIR not set)";

   if (sdpc_root_dir == NULL)
     throw ApplicationError, "*** Error: SDPC_ROOT is not set";

   variable scan_dir, scan_dir_list = __argv[[__i:]];
   variable products;

   foreach scan_dir (scan_dir_list)
     {
        if (cmdline_products != NULL)
          {
             products = cmdline_products;
          }
        else
          {
             products = find_product_dirs (scan_dir);
          }

        process_scan_granules (scan_dir, archive_root_dir, products);
     }
}
