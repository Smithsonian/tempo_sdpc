#! /usr/bin/env slsh

require ("process");
require ("chksum");
require ("cmdopt");
require ("pcre");

private variable Ancillary_Type_List = ["GEOSCF", "CMIEAST", "CMIWEST", "IMS"];

private variable GOES_Path_Pattern = "/20\d{2}/\d{3}/(east|west)_cmi/OR_ABI-"R;
private variable GOES_Path_Regex = pcre_compile (GOES_Path_Pattern);

private variable Node_Name_Entry;
private variable Dest_Target_Dir;

define make_file_entry (path, data_type, st, file_type)
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
   s.file_type = file_type;
   s.file_id = path_basename (path);
   s.file_size = st.st_size;
   s.file_chksum_type = "MD5";
   s.file_chksum = md5sum_file (path);

   return s;
}

define process_file_corrfile (types, path)
{
   variable st = stat_file (path);
   if (st == NULL)
     return;

   variable basename = path_basename (path);
   variable tok = strtok (basename, "_");
   % example:  TEMPO_RADREF_L1_V01_YYYYMMDD_S123456789_E123456789_S003.nc
   % example:  TEMPO_DSTRHCHO_L2_V01_YYYYMMDD_S123456789_E123456789_S003.nc
   variable product_type = tok[1];
   variable version_string = strtrim_beg (tok[3], "V0");  % e.g. 1
   variable data_type = "TEMPO_NONORDERABLE";
   variable entry = make_file_entry (path, data_type, st, "SCIENCE");

   variable group = struct
     {
        data_version = version_string,
        entry = entry,
        met_entry = NULL,
        json_entry = NULL
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

define process_file_nc (types, path, is_aws_upload)
{
   variable st = stat_file (path);
   if (st == NULL)
     return;

   variable basename = path_basename (path);
   variable tok = strtok (basename, "_");
   variable product_type = strjoin (tok[[1:3]], "_");
   variable data_type, data_version;
   if (0 != is_substr (basename, "_L0_"))
     {
        data_type = "TEMPO_NONORDERABLE";
        data_version = strtrim_beg (tok[3], "V0");   % e.g. 1
     }
   else
     {
        data_type = strjoin (tok[[0:2]], "_");       % e.g. TEMPO_RAD_L1
        data_version = tok[3];                       % e.g. V01
     }
   variable entry = make_file_entry (path, data_type, st, "SCIENCE");

   variable path_met = path + ".met";
   variable st_met = stat_file (path_met);
   variable met_entry = NULL;
   if (st_met != NULL and is_aws_upload == 0)
     met_entry = make_file_entry (path_met, data_type, st_met, "METADATA");

   variable path_json = path + ".cmr.json";
   variable st_json = stat_file (path_json);
   variable json_entry = NULL;
   if (st_json != NULL and is_aws_upload != 0)
     json_entry = make_file_entry (path_json, data_type, st_json, "METADATA");

   variable group = struct
     {
        data_version = data_version,
        entry = entry,
        met_entry = met_entry,
        json_entry = json_entry
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

define process_file_raw (types, path)
{
   variable st = stat_file (path);
   if (st == NULL)
     return;

   variable basename = path_basename (path);
   variable tok = strtok (basename, "_");
   variable product_type = "RAW";
   variable data_type    = "TEMPO_NONORDERABLE";
   variable data_version = "1";
   variable entry = make_file_entry (path, data_type, st, "SCIENCE");

   variable group = struct
     {
        data_version = data_version,
        entry = entry,
        met_entry = NULL,
        json_entry = NULL
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

define process_file_ancillary (types, path)
{
   variable st = stat_file (path);
   if (st == NULL)
     return;

   variable basename = path_basename (path);

   variable product_type = NULL;
   if (0 == strncmp ("GEOS-CF", basename, 7))
     {
        product_type = "GEOSCF";
     }
   else if (0 == strncmp ("OR_ABI", basename, 6))
     {
        % path matches /east_cmi/ tree => product_type=CMIEAST
        % path matches /west_cmi/ tree => product_type=CMIWEST
        variable tok = pcre_matches (GOES_Path_Regex, path);
        product_type = "CMI" + strup(tok[1]);
     }
   else if (0 == strncmp ("ims", basename, 3))
     {
        product_type = "IMS";
     }
   if (product_type == NULL)
     {
        () = fprintf (stderr, "*** unrecognized file type: skipping file: %s\n", path);
        return;
     }

   variable data_type    = "TEMPO_NONORDERABLE";
   variable data_version = "1";
   variable entry = make_file_entry (path, data_type, st, "SCIENCE");

   variable group = struct
     {
        data_version = data_version,
        entry = entry,
        met_entry = NULL,
        json_entry = NULL
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

define entry_string (entry)
{
   variable
     id = entry.file_id,
     size = entry.file_size,
     type = entry.file_type,
     chksum = entry.file_chksum,
     chksum_type = entry.file_chksum_type;

   variable str =
` OBJECT = FILE_SPEC;
    DIRECTORY_ID = $Dest_Target_Dir;
    FILE_ID = $id;
    FILE_TYPE = $type;
    FILE_SIZE = $size;
    FILE_CKSUM_TYPE = $chksum_type;
    FILE_CKSUM_VALUE = "$chksum";
  END_OBJECT = FILE_SPEC;`$;

   return str;
}

define write_file_group (fp, g)
{
   variable data_version = g.data_version;
   variable data_type = g.entry.data_type;

   variable entries = entry_string (g.entry);

   if (g.met_entry != NULL)
     {
        variable met_str = entry_string (g.met_entry);
        entries = sprintf ("%s\n%s", entries, met_str);
     }

   if (g.json_entry != NULL)
     {
        variable json_str = entry_string (g.json_entry);
        entries = sprintf ("%s\n%s", entries, json_str);
     }

   variable str =
`OBJECT = FILE_GROUP;
  DATA_TYPE = $data_type;
  DATA_VERSION = $data_version;
  $Node_Name_Entry$entries
END_OBJECT = FILE_GROUP;
`$;

   () = fputs (str, fp);
}

define write_manifest (lst, filename)
{
   variable g, num_files = 0;

   foreach g (lst)
     {
        num_files += 1;
        if (g.met_entry != NULL) num_files += 1;
        if (g.json_entry != NULL) num_files += 1;
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
        write_file_group (fp, g);
     }

   if (0 != fclose (fp))
     throw IOError, "closing $filename"$;
}

define print_mput (fp, ary)
{
   () = fprintf (fp, "mput -P 4 \\\n");
   () = fprintf (fp, "%s\n", strjoin (ary, " \\\n"));
}

define print_script_entry (fp, ary)
{
   variable has_met = array_map (Int_Type, &is_substr, ary, ".nc.met");
   variable i, j = where (has_met, &i);
   if (length(i) > 0) print_mput (fp, ary[i]);
   if (length(j) > 0) print_mput (fp, ary[j]);
}

define print_list_entry (fp, ary)
{
   variable has_met = array_map (Int_Type, &is_substr, ary, ".nc.met");
   variable has_json = array_map (Int_Type, &is_substr, ary, ".nc.cmr.json");
   variable i = where (has_met == 0 and has_json == 0);
   variable j = where (has_json);
   variable k = where (has_met);
   if (length(i) > 0) () = fprintf (fp, "%s\n", strjoin (ary[i], "\n"));
   if (length(j) > 0) () = fprintf (fp, "%s\n", strjoin (ary[j], "\n"));
   if (length(k) > 0) () = fprintf (fp, "%s\n", strjoin (ary[k], "\n"));
}

define write_lftp_script (dest, types, type_list, pdr_files, script_file)
{
   variable fp = fopen (script_file, "w");
   if (NULL == fp)
     throw IOError, "opening $script_file for writing"$;

   () = fprintf (fp, "open --user %s --password %s sftp://%s\n",
                 dest.user, dest.password, dest.host);
   %() = fprintf (fp, "set xfer:log-file lftp_log.%s\n", strftime ("%Y%m%dT%H%M%SZ", gmtime(_time)));
   () = fprintf (fp, "set xfer:use-temp-file yes\n");
   () = fprintf (fp, "set xfer:log-file lftp.log\n");
   % Use a relative path to work around some confusion at ASDC
   () = fprintf (fp, "cd %s || exit\n", dest.subdir);

   % Transfer data files first, then the manifest (PDR) file.
   % The ASDC ingest system assumes the PDR file is uploaded last.
   % When a PDR file is detected, the system immediately starts
   % looking for the files it refers to, and if those files aren't
   % available, the ingest system fails with "FILE NOT FOUND".

   % "mput -P N" means transfer N files in parallel.
   % In testing data transfer of 2GB files from SAO to ASDC,
   % N=4 reduced the transfer time from 296 sec to 84 sec,
   % a factor of 3.5 speedup.

   variable i, num_types = length(type_list);

   _for i (0, num_types-1, 1)
     {
        variable this_type = type_list[i];
        variable g, lst = {};
        foreach g (types[this_type])
          {
             list_append (lst, g.entry.path);
             if (g.met_entry != NULL)
               list_append (lst, g.met_entry.path);
             % cmr.json metadata not used with direct upload to ASDC
          }
        if (length(lst) > 0)
          {
             variable ary = list_to_array (lst);
             % To avoid long lines, continue lines with a trailing '\\n'
             print_script_entry (fp, ary);
             () = fprintf (fp, "put %s\n", pdr_files[i]);
          }
     }

   () = fputs ("exit\n", fp);

   if (0 != fclose (fp))
     throw IOError, "closing $script_file"$;
}

define write_aws_upload_sequence (dest, types, type_list, pdr_files, script_file)
{
   variable fp = fopen (script_file, "w");
   if (NULL == fp)
     throw IOError, "opening $script_file for writing"$;

   % Transfer data files first, then the manifest (PDR) file.
   % The ASDC ingest system assumes the PDR file is uploaded last.
   % When a PDR file is detected, the system immediately starts
   % looking for the files it refers to, and if those files aren't
   % available, the ingest system fails with "FILE NOT FOUND".

   variable i, num_types = length(type_list);

   _for i (0, num_types-1, 1)
     {
        variable this_type = type_list[i];
        variable g, lst = {};
        foreach g (types[this_type])
          {
             list_append (lst, g.entry.path);
             % if (g.met_entry != NULL)
             %   list_append (lst, g.met_entry.path);
             if (g.json_entry != NULL)
               list_append (lst, g.json_entry.path);
          }
        if (length(lst) > 0)
          {
             variable ary = list_to_array (lst);
             % To avoid long lines, continue lines with a trailing '\\n'
             print_list_entry (fp, ary);
             () = fprintf (fp, "%s\n", pdr_files[i]);
          }
     }

   if (0 != fclose (fp))
     throw IOError, "closing $script_file"$;
}

define make_manifest_filename (type)
{
   % Ensure that two separate pipelines cannot generate the same PDR filename.
   % SDPC_PIPE_ID uniquely identifies each pipeline.
   variable pipe_id_env = getenv ("SDPC_PIPE_ID");
   variable pipe_id = (pipe_id_env == NULL) ? "" : "_" + pipe_id_env;

   if (any (type == Ancillary_Type_List))
     return strftime ("${type}_%Y%m%dT%H%M%SZ${pipe_id}.PDR"$, gmtime(_time));
   else
     return strftime ("TEMPO_${type}_%Y%m%dT%H%M%SZ${pipe_id}.PDR"$, gmtime(_time));
}

define upload_priority (product_type)
{
   variable tok = string_matches (product_type, "_L[0-3]_");
   if (NULL == tok)
     {
        % assign Level 0 priority (raw and ancillary data)
        return 0;
     }
   variable level;
   () = sscanf (tok[0], "_L%d_", &level);
   return level;
}

define priority_ordered_type_list (types)
{
   variable type_list = assoc_get_keys (types);
   if (length(type_list) == 0)
     return type_list;
   variable priority = array_map (Integer_Type, &upload_priority, type_list);
   % Sort the types by descending order of upload priority
   return type_list [reverse(array_sort(priority))];
}

define write_pdr_list (pdr_names, list_file)
{
   if (list_file == NULL)
     return;

   variable fp = fopen (list_file, "w");
   if (NULL == fp)
     throw IOError, "opening $list_file for writing"$;

   if (1 != fputslines (strjoin (pdr_names, "\n"), fp))
     throw IOError, "writing $list_file"$;

   if (0 != fclose (fp))
     throw IOError, "closing $list_file"$;
}

define process_file_list (dest, file_list, script_file, is_aws_upload, pdr_file_list)
{
   variable path_list = read_file_list (file_list);

   variable path, types = Assoc_Type[];

   foreach path (path_list)
     {
        if (0 == strlen(strtrim(path)))
          continue;
        variable st = stat_file (path);
        if (st == NULL)
          {
             vmessage ("file not found:  $path"$);
             continue;
          }
        variable basename = path_basename (path);
        % Filter correction files by their prefixes.
        % Ancillary data files do not have a "TEMPO" prefix.
        % TEMPO data products are prefixed with "TEMPO_"
        % TEMPO raw tar files are prefixed with "tempo_"
        if ((0 == strncmp ("TEMPO_RADREF", basename, 12))
            || (0 == strncmp ("TEMPO_DSTR", basename, 10)))
          {
             process_file_corrfile (types, path);
             continue;
          }
        else if (0 != strncmp ("TEMPO", strup(basename), 5))
          {
             process_file_ancillary (types, path);
             continue;
          }
        variable extname = path_extname (path);
        if (extname == ".met" or 0 != string_match (path, ".cmr.json$"))
          {
             % Silently skip .met files and .cmr.json files  They get handled
             % along with the corresponding data file.
             continue;
          }
        else if (extname == ".nc")
          {
             process_file_nc (types, path, is_aws_upload);
          }
        else if (extname == ".tar")
          {
             process_file_raw (types, path);
          }
        else
          {
             % Complain about unrecognized file extensions
             () = fprintf (stderr, "*** unrecognized file type: skipping file: %s\n", path);
          }
     }

   variable type_list = priority_ordered_type_list (types);
   if (length(type_list) == 0)
     return;

   variable mf_files = array_map (String_Type, &make_manifest_filename, type_list);

   write_pdr_list (mf_files, pdr_file_list);

   variable i, num_types = length(type_list);

   _for i (0, num_types-1, 1)
     {
        write_manifest (types[type_list[i]], mf_files[i]);
     }

   if (is_aws_upload == 0)
     {
        write_lftp_script (dest, types, type_list, mf_files, script_file);
     }
   else
     {
        write_aws_upload_sequence (dest, types, type_list, mf_files, script_file);
     }
}

private define usage ()
{
   variable msg =
`Usage: asdc_mkscript.sl [options] <files.lis>
Options:
    -d|--dest USER@HOST:dirpath     Destination account/host/path
    -o|--output FILE                Write lftp script to FILE
    -b|--bucket Bucket:Bucket_dir   Generate output script for AWS S3 upload
    -p|--pdr FILE                   Write PDR filenames to FILE
    -h|--help                       Show usage message
`;
   () = fprintf (stderr, msg);
   exit (0);
}

private define cmdopt_error (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

private define get_hostname ()
{
   variable name = "$HOST"$;
   if (0 != strlen(name))
     return name;

   % In a cron job environment, HOST may not be set.
   % In that case, we run /bin/hostname, which should exist.

   variable obj, s;

   obj = new_process ("/bin/hostname"; write=1);
   name = fgetslines (obj.fp1);
   s = obj.wait();
   () = fclose (obj.fp1);

   if (s.exit_status != 0)
     {
        name = "unknown";
     }

   return strtrim (name[0], "\n\t");
}

define slsh_main ()
{
   variable show_usage = 0;
   variable script_file = "script.lftp";
   variable pdr_file_list = NULL;
   variable user_at_host = "tempo@xfr140.larc.nasa.gov:ingest/tempo";
   variable s3_bucket = NULL;

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("h|help", &show_usage; inc);
   opts.add ("b|bucket", &s3_bucket; type="string");
   opts.add ("d|dest", &user_at_host; type="string");
   opts.add ("o|output", &script_file; type="string");
   opts.add ("p|pdr", &pdr_file_list; type="string");
   variable i = opts.process (__argv,1);

   if (__argc != i+1 || show_usage != 0)
     usage();

   variable dest = struct
     {
        host, user, password = "DUMMY", subdir
        % ASDC (password DUMMY causes lftp to use the ssh-agent for authentication)
     };
   variable tok = strtok (user_at_host, "@");
   dest.user = tok[0];
   variable htok = strtok (tok[1], ":");
   dest.host = htok[0];
   dest.subdir = htok[1];

   variable file_list_file = __argv[i];

   if (NULL == stat_file (file_list_file))
     throw IOError, "Product list file not found: $file_list_file"$;

   if (NULL != stat_file (script_file))
     throw IOError, "Target file exists: $script_file"$;

   % PDR files differ somewhat between ASDC direct ingest and cloud ingest:
   variable is_aws_upload = 0;
   if (s3_bucket == NULL)
     {
        Node_Name_Entry = sprintf ("NODE_NAME = %s;\n", get_hostname());
        Dest_Target_Dir = ".";
     }
   else
     {
        variable bckt = strtok (s3_bucket, ":");
        if (length(bckt) != 2) usage();
        Node_Name_Entry = "";
        Dest_Target_Dir = bckt[1];
        is_aws_upload = 1;
     }

   process_file_list (dest, file_list_file, script_file, is_aws_upload, pdr_file_list);
}
