#! /usr/bin/env slsh

require ("process");
require ("chksum");
require ("cmdopt");

private variable Node_Name = "$HOST"$;

define make_file_entry (path, data_type, st, extname_is_nc)
{
   variable s = struct
     {
        path = path,
        data_type,
        file_id,
        file_type, file_size,
        file_chksum, file_chksum_type
     };

   s.data_type = data_type;
   s.file_type = extname_is_nc ? "SCIENCE" : "METADATA";
   s.file_id = path_basename (path);
   s.file_size = st.st_size;
   s.file_chksum_type = "MD5";
   s.file_chksum = md5sum (path);

   return s;
}

define process_file (types, path_nc)
{
   variable extname = path_extname (path_nc);
   if (extname != ".nc")
     {
        () = fprintf (stderr, "*** skipping file: %s\n", path_nc);
        return;
     }

   variable st_nc = stat_file (path_nc);
   if (st_nc == NULL)
     return;

   variable basename_nc = path_basename (path_nc);
   variable tok = strtok (basename_nc, "_");
   variable data_type = strjoin (tok[[0:2]], "_");
   variable product_type = strjoin (tok[[1:3]], "_");
   variable data_version = atoi(strtrim_beg (tok[3], "V"));
   variable nc_entry = make_file_entry (path_nc, data_type, st_nc, 1);

   variable path_met = path_nc + ".met";
   variable st_met = stat_file (path_met);
   variable met_entry = NULL;
   if (st_met != NULL)
     met_entry = make_file_entry (path_met, data_type, st_met, 0);

   variable group = struct
     {
        data_version = data_version,
        nc_entry = nc_entry,
        met_entry = met_entry
     };

   if (assoc_key_exists (types, product_type))
     {
        list_append (types[product_type], group);
        return;
     }

   variable lst = {};
   list_append (lst, group);
   types[product_type] = lst;
}

define read_file_list (list_file)
{
   variable fp = fopen (list_file, "r");
   if (fp == NULL)
     throw IOError, "opening ${list_file} for reading"$;

   variable lst = fgetslines (fp);
   () = fclose (fp);

   return array_map (String_Type, &strtrim, lst, "\n\t");
}

define entry_string (entry, target_dir)
{
   variable
     id = entry.file_id,
     size = entry.file_size,
     type = entry.file_type,
     chksum = entry.file_chksum,
     chksum_type = entry.file_chksum_type;

   variable str =
` OBJECT = FILE_SPEC;
    DIRECTORY_ID = $target_dir;
    FILE_ID = $id;
    FILE_TYPE = $type;
    FILE_SIZE = $size;
    FILE_CKSUM_TYPE = $chksum_type;
    FILE_CKSUM_VALUE = "$chksum";
  END_OBJECT = FILE_SPEC;`$;

   return str;
}

define write_file_group (fp, g, target_dir)
{
   variable data_version = g.data_version;
   variable data_type = g.nc_entry.data_type;

   variable entries = entry_string (g.nc_entry, target_dir);

   if (g.met_entry != NULL)
     {
        variable met_str = entry_string (g.met_entry, target_dir);
        entries = sprintf ("%s\n%s", entries, met_str);
     }

   variable str =
`OBJECT = FILE_GROUP;
  DATA_TYPE = $data_type;
  DATA_VERSION = $data_version;
  NODE_NAME = $Node_Name;
  $entries
END_OBJECT = FILE_GROUP;
`$;

   () = fputs (str, fp);
}

define write_manifest (dest, lst, filename)
{
   variable g, num_files = 0;

   foreach g (lst)
     {
        if (g.met_entry == NULL)
          num_files += 1;
        else
          num_files += 2;
     }

   variable hdr =
`ORIGINATING_SYSTEM = TEMPO;
TOTAL_FILE_COUNT = $num_files
`$;

   variable fp = fopen (filename, "w");
   if (fp == NULL)
     throw IOError, "opening $filename for writing"$;

   () = fputs (hdr, fp);

   foreach g (lst)
     {
        write_file_group (fp, g, dest.target_dir);
     }

   if (0 != fclose (fp))
     throw IOError, "closing $filename"$;
}

define make_target_dir ()
{
   return strftime ("tempo_test_submission_%Y-%m-%d", gmtime(_time()));
}

define write_lftp_script (dest, types, pdr_files, script_file)
{
   variable fp = fopen (script_file, "w");
   if (NULL == fp)
     throw IOError, "opening $script_file for writing"$;

   () = fprintf (fp, "open --user %s --password %s sftp://%s\n",
                 dest.user, dest.password, dest.host);
   () = fprintf (fp, "mkdir %s\n", dest.target_dir);
   () = fprintf (fp, "cd %s\n", dest.target_dir);

   variable f, g, t;

   foreach f (pdr_files)
     {
        () = fprintf (fp, "put %s\n", f);
     }

   foreach t (types) using ("keys")
     {
        foreach g (types[t])
          {
             () = fprintf (fp, "put %s\n", g.nc_entry.path);
             if (g.met_entry != NULL)
               {
                  () = fprintf (fp, "put %s\n", g.met_entry.path);
               }
          }
     }

   () = fputs ("exit\n", fp);

   if (0 != fclose (fp))
     throw IOError, "closing $script_file"$;
}

define make_manifest_filename (type)
{
   return strftime ("TEMPO_${type}_%Y%m%dT%H%M%SZ.PDR"$, gmtime(_time));
}

define process_file_list (dest, nc_file_list, script_file)
{
   variable nc_paths = read_file_list (nc_file_list);

   variable types = Assoc_Type[];
   variable path;

   foreach path (nc_paths)
     {
        process_file (types, path);
     }

   variable type_list = assoc_get_keys (types);
   variable mf_files = array_map (String_Type, &make_manifest_filename, type_list);

   variable i, num_types = length(type_list);

   _for i (0, num_types-1, 1)
     {
        write_manifest (dest, types[type_list[i]], mf_files[i]);
     }

   write_lftp_script (dest, types, mf_files, script_file);
}

private define read_netrc_file (file)
{
   variable fp = fopen (file, "r");
   if (fp == NULL)
     throw IOError, "Error opening $file for reading"$;
   variable s = fgetslines (fp);
   () = fclose (fp);

   variable dest = struct
     {
        host, user, password, target_dir
     };

   variable i, n = length(s);
   _for i (0, n-1, 1)
     {
        if (s[i][0] == '#')
          continue;
        variable num = sscanf (s[i], "machine %s login %s password %s dir %s",
                               &dest.host,
                               &dest.user,
                               &dest.password,
                               &dest.target_dir);
        if (num != 4)
          throw IOError, "Error parsing netrc file: $file"$;
     }

   return dest;
}

private define usage ()
{
   variable msg =
`Usage: asdc_mkscript.sl [options] <files.lis>
Options:
    -n|--netrc FILE       Read host,user,passwd,target_dir from FILE
    -o|--output FILE      Write lftp script to FILE
    -h|--help             Show usage message
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
   variable show_usage = 0;
   variable netrc_file = NULL;
   variable script_file = "script.lftp";

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("h|help", &show_usage; inc);
   opts.add ("n|netrc", &netrc_file; type="string");
   opts.add ("o|output", &script_file; type="string");
   variable i = opts.process (__argv,1);

   if (__argc != i+1 || show_usage != 0)
     usage();

   variable dest = read_netrc_file (netrc_file);
   dest.target_dir = path_concat (dest.target_dir, make_target_dir ());

   variable nc_file_list_file = __argv[i];

   if (NULL == stat_file (nc_file_list_file))
     throw IOError, "Product list file not found: $nc_file_list_file"$;

   if (NULL != stat_file (script_file))
     throw IOError, "Target file exists: $script_file"$;

   process_file_list (dest, nc_file_list_file, script_file);
}
