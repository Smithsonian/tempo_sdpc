#! /usr/bin/env slsh

require ("process");
require ("chksum");
require ("cmdopt");

private variable Node_Name = "$HOST"$;

private variable Dest_Subdir   = "ingest/tempo";
private variable Dest_Target_Dir = ".";
% Weirdly, absolute directory paths viewed from inside ASDC are
% different from those viewed from outside ASDC. To avoid
% confusion, we use "."

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

define process_file_nc (types, path)
{
   variable st = stat_file (path);
   if (st == NULL)
     return;

   variable basename = path_basename (path);
   variable tok = strtok (basename, "_");
   variable data_type = strjoin (tok[[0:2]], "_");
   variable product_type = strjoin (tok[[1:3]], "_");
   variable data_version = tok[3];
   variable entry = make_file_entry (path, data_type, st, "SCIENCE");

   variable path_met = path + ".met";
   variable st_met = stat_file (path_met);
   variable met_entry = NULL;
   if (st_met != NULL)
     met_entry = make_file_entry (path_met, data_type, st_met, "METADATA");

   variable group = struct
     {
        data_version = data_version,
        entry = entry,
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

define process_file_raw (types, path)
{
   variable st = stat_file (path);
   if (st == NULL)
     return;

   variable basename = path_basename (path);
   variable tok = strtok (basename, "_");
   variable data_type    = strjoin (tok[[0:1]], "_");  % e.g. TEMPO_GRDDP
   variable product_type = tok[1];                     % e.g.       GRDDP
   variable data_version = "V01";    % no versions for unprocessed telemetry
   variable entry = make_file_entry (path, data_type, st, "SCIENCE");

   variable group = struct
     {
        data_version = data_version,
        entry = entry,
        met_entry = NULL
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
        write_file_group (fp, g);
     }

   if (0 != fclose (fp))
     throw IOError, "closing $filename"$;
}

define write_lftp_script (dest, types, pdr_files, script_file)
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
   () = fprintf (fp, "cd %s\n", Dest_Subdir);

   variable f, g, t;

   foreach t (types) using ("keys")
     {
        foreach g (types[t])
          {
             () = fprintf (fp, "put %s\n", g.entry.path);
             if (g.met_entry != NULL)
               {
                  () = fprintf (fp, "put %s\n", g.met_entry.path);
               }
          }
     }

   % The ASDC ingest system assumes the PDR files are uploaded last.
   % When the PDR is detected, the system immediately starts
   % looking for files, and if the files aren't available, the
   % ingest system fails with "FILE NOT FOUND".
   % Yet another rule not mentioned in the ICD...
   foreach f (pdr_files)
     {
        () = fprintf (fp, "put %s\n", f);
     }

   () = fputs ("exit\n", fp);

   if (0 != fclose (fp))
     throw IOError, "closing $script_file"$;
}

define make_manifest_filename (type)
{
   return strftime ("TEMPO_${type}_%Y%m%dT%H%M%SZ.PDR"$, gmtime(_time));
}

define process_file_list (dest, file_list, script_file)
{
   variable path_list = read_file_list (file_list);

   variable types = Assoc_Type[];
   variable path, extname;

   foreach path (path_list)
     {
        extname = path_extname (path);
        if (extname == ".met")
          {
             % Silently skip .met files.  They get handled
             % along with the corresponding data file.
             continue;
          }
        else if (extname == ".nc")
          {
             process_file_nc (types, path);
          }
        else if (extname == ".raw")
          {
             process_file_raw (types, path);
          }
        else
          {
             % Complain about unrecognized file extensions
             () = fprintf (stderr, "*** unrecognized file type: skipping file: %s\n", path);
          }
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

private define usage ()
{
   variable msg =
`Usage: asdc_mkscript.sl [options] <files.lis>
Options:
    -d|--dest USER@HOST   Destination account
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
   variable script_file = "script.lftp";
   variable user_at_host = "tempo@xfr140.larc.nasa.gov";

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("h|help", &show_usage; inc);
   opts.add ("d|dest", &user_at_host; type="string");
   opts.add ("o|output", &script_file; type="string");
   variable i = opts.process (__argv,1);

   if (__argc != i+1 || show_usage != 0)
     usage();

   variable dest = struct
     {
        host, user, password = "DUMMY"
        % ASDC (password DUMMY causes lftp to use the ssh-agent for authentication)
     };
   variable tok = strtok (user_at_host, "@");
   dest.user = tok[0];
   dest.host = tok[1];

   variable file_list_file = __argv[i];

   if (NULL == stat_file (file_list_file))
     throw IOError, "Product list file not found: $file_list_file"$;

   if (NULL != stat_file (script_file))
     throw IOError, "Target file exists: $script_file"$;

   process_file_list (dest, file_list_file, script_file);
}
