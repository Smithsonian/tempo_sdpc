#! /usr/bin/env slsh

require ("cmdopt");
require ("chksum");
require ("csv");
require ("process");

$1 = path_dirname (__FILE__);
prepend_to_slang_load_path (path_concat ($1, "../share/slsh/local-packages"));
prepend_to_slang_load_path ($1);
require ("pipeutil");

private variable Delete_Tarfiles = 0;
private variable Clobber_Output_Files = 0;
private variable Archive_Subdir_File = "archive_subdir";

private variable Archive_Root_Dir;

private define usage ()
{
   variable argv0 = __argv[0];
   variable s =
`Usage:  $argv0 [options] DIR
     or $argv0 [options] tarfile1 tarfile2 ...
      Options:
         -a|--archive_root DIR
         -c|--clobber        If present, will overwrite pre-existing files
         -d|--delete         If present, will delete input tar files
                             after archiving the contents.
         -l|--level LEVEL    Data product level (L1 | L2 | L3)
         -h|--help
`$;
   vmessage (s);
   exit(0);
}

private define error_routine (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

private define assert_nonexistent (path)
{
   if (Clobber_Output_Files)
     return;
   if (NULL != stat_file (path))
     throw ApplicationError, "*** Error: exists $path"$;
}

private define compare_checksums (src_path, cpy_path)
{
   variable sha1_src = sha1sum_file (src_path);
   variable sha1_cpy = sha1sum_file (cpy_path);
   return strcmp (sha1_cpy, sha1_src);
}

define copy_file (file_path, dest_dir)
{
   variable basename = strtrim_beg (path_basename(file_path), ".");
   variable dest_path = path_concat (dest_dir, basename);
   assert_nonexistent (dest_path);

   % Create the destination directory
   if (0 != mkdir_p (dest_dir))
     {
        throw ApplicationError, "*** Error: creating $dest_dir"$;
     }

   variable s = new_process (["/bin/cp", file_path, dest_path]).wait();
   if (s.exit_status != 0)
     throw ApplicationError, "*** Error: copy failed: $file_path"$;

   % Verifying checksums at this point might be excessive.
   if (0 != compare_checksums (file_path, dest_path))
     {
        throw ApplicationError,
          "*** Error: file copy checksum mismatch: $file_path"$;
     }

   return dest_path;
}

private define get_tarfile_archive_subdir (tar_file)
{
   % Extract the archive_subdir file.

   variable dir = path_dirname (tar_file);
   variable tar_file_basename = path_basename (tar_file);

   variable name_fields = strtok (tar_file_basename, ".");
   variable granule_name = name_fields[0];

   variable argv = ["tar", "-x", "-f", tar_file, "--to-stdout",
                    path_concat (granule_name, Archive_Subdir_File)];

   variable obj = new_process(argv; write=1);
   variable archive_subdir = fgetslines (obj.fp1);
   variable s = obj.wait();
   if (s.exit_status != 0)
     {
        throw ApplicationError,
          "*** Error: extracting ${Archive_Subdir_File} from tar archive: $tar_file"$;
     }
   () = fclose (obj.fp1);

   archive_subdir = array_map (String_Type, &strtrim, archive_subdir, "\n");
   return archive_subdir[0];
}

define insert_fixed_metadata (path)
{
   variable argv = ["insert_fixed_metadata.py", path];
   variable s = new_process (argv; dup2=1).wait();
   if (s.exit_status != 0)
     throw ApplicationError, "*** Error: inserting fixed metadata: $path"$;
}

define fix_met_file_format (path)
{
   variable argv = ["fix_met_format.py", path];
   variable s = new_process (argv; dup2=1).wait();
   if (s.exit_status != 0)
     throw ApplicationError, "*** Error: fixing met file format $path"$;
}

define make_public_mirror_symlink (oldpath)
{
   % If the mirror directory does not exist, do nothing.
   variable mirror_dir = "$SDPC_RUN_DIR_MASTER/public_mirror"$;
   if (NULL == stat_file (mirror_dir))
     return;

   variable basename = path_basename (oldpath);
   variable argv, obj, s, result;

   % Filter out selected data products, either because we don't
   % release them at all, or because we release them elsewhere.
    if (is_substr (basename, "TEMPO_NO2_L2"))
     {
        % TEMPO_NO2_L2 files are released after L2_split runs
        return;
     }
   else if (is_substr (basename, "TEMPO_RAD_L1"))
     {
        % TEMPO_RAD_L1 files are released after INR processing
        argv = ["radiance_inr_status.py", oldpath];
        obj = new_process (argv; write=1);
        result = fgetslines (obj.fp1);
        s = obj.wait();
        if (s.exit_status != 0)
          {
             throw ApplicationError, "*** Error: examining L1 radiance file: $oldpath"$;
          }
        variable inr_status = eval(result[0]);
        if (inr_status < 1)
          return;
     }

   % Construct the target directory path:
   argv = ["level1_info", "--localday", oldpath];
   obj = new_process (argv; write=1);
   result = fgetslines (obj.fp1);
   s = obj.wait();
   if (s.exit_status != 0)
     {
        throw ApplicationError, "*** Error: constructing target directory path: $oldpath"$;
     }
   variable sat_local_day = result[0];
   variable name_fields = strtok (basename, "_");
   % name_fields[2] = L1 | L2 | L3
   variable target_dir = sprintf ("%s/%s/%s", mirror_dir, sat_local_day, name_fields[2]);
   if (0 != mkdir_p (target_dir))
     {
        throw ApplicationError, "*** Error: creating $target_dir"$;
     }

   % Create the symbolic link
   variable newpath = path_concat (target_dir, basename);
   if (0 != symlink (oldpath, newpath))
     throw ApplicationError, "*** Error: creating symlink $newpath"$;
}

define register_using_symlink (tar_file, archive_dest_subdir)
{
   % Extract partial paths to archived data products.
   % Some tar files contain block_??? subdirectories with .nc files
   % that will eventually be merged to generate the final product.
   % These block .nc files should not be entered in the product
   % database, so we filter them out of this query.
   % Note that these tar options will yield a list of partial paths
   % that contain basenames like *TEMPO_*.nc so we will need to filter
   % out any files with undesired prefixes.
   variable argv = ["tar", "tf", tar_file,
		    "--exclude=block_*", "--strip-components=1",
                    "--show-transformed-names", "--wildcards",
                    "--no-anchored", "TEMPO_*.nc"];
   %vmessage (strjoin (argv, " "));
   variable obj = new_process(argv; write=1);
   variable partial_paths = fgetslines (obj.fp1);
   variable s = obj.wait();
   if (s.exit_status != 0)
     {
        throw ApplicationError, "*** Error: reading product file names from: $tar_file"$;
     }
   () = fclose (obj.fp1);

   partial_paths = array_map (String_Type, &strtrim, partial_paths, "\n");

   % "partial_path" means something like
   %      HCHO/TEMPO_HCHO_L2_V01_20130715T165956Z_S002G01.nc

   % For each product file, trigger registration in the product
   % database by making a symbolic link in $incoming_dir
   % in the archive directory (usually on the master node)
   variable incoming_dir = path_concat (Archive_Root_Dir, "registry/incoming");
   if (NULL == stat_file (incoming_dir))
     {
        if (0 != mkdir_p (incoming_dir))
          {
             throw ApplicationError, "*** Error: creating $incoming_dir"$;
          }
     }

   variable pp, product, oldpath, newpath;
   foreach pp (partial_paths)
     {
        % Some files are to be excluded.
        %  => skip basenames with unwanted prefixes:
        variable pp_basename = path_basename (pp);
        if (0 != strncmp (pp_basename, "TEMPO_", 6))
          continue;
        variable exclude_substrs = ["TEMPO_INR", "_diag.nc"];
        if (any(array_map (Integer_Type, &is_substr, pp_basename, exclude_substrs)))
          continue;

        oldpath = path_concat (archive_dest_subdir, pp);
	if (NULL == stat_file (oldpath))
	  continue;

        % insert fixed metadata
        insert_fixed_metadata (oldpath);

        % fix formatting of .met file
        variable metfile_path = oldpath + ".met";
        if (NULL != stat_file (metfile_path))
          fix_met_file_format (metfile_path);

        % create symbolic link to trigger product registration
        newpath = path_concat (incoming_dir, pp_basename);
        if (0 != symlink (oldpath, newpath))
          throw ApplicationError, "*** Error: creating symlink $newpath"$;

        % create symbolic link for public mirror
        make_public_mirror_symlink (oldpath);
     }

   return 0;
}

define unpack_and_archive (tar_file, archive_level_dir)
{
   % The unpack and archive process involves the following steps:
   %  1. Extract the ${Archive_Subdir_File} file
   %  2. Use ${Archive_Subdir_File} to define the destination sub-directory
   %  3. Unpack the contents in the appropriate destination
   %     directory, avoiding file collisions.

   variable subdir = get_tarfile_archive_subdir (tar_file);
   variable archive_dest_subdir = path_concat (archive_level_dir, subdir);

   % Create the destination directory
   if (0 != mkdir_p (archive_dest_subdir))
     {
        throw ApplicationError, "*** Error: creating $archive_dest_subdir"$;
     }

   variable argv = ["tar", "-x"];

   if (Clobber_Output_Files == 0)
     {
        argv = [argv, "--keep-old-files"];
     }

   %if (NULL != stat_file (path_concat (archive_dest_subdir, Archive_Subdir_File)))
   argv = [argv, "--exclude=${Archive_Subdir_File}"$];

   argv = [argv, "-f", tar_file, "-C", archive_dest_subdir,
           "--strip-components=1"];

   %vmessage (strjoin (argv, " "));
   variable unpack_log = "${tar_file}.unpack"$;
   variable s = new_process(argv; stdout=unpack_log, dup2=1).wait();
   if (s.exit_status != 0)
     {
        throw ApplicationError, "*** Error: unpacking $tar_file"$;
     }
   else
     {
        () = remove (unpack_log);
        () = register_using_symlink (tar_file, archive_dest_subdir);
     }

   % It's now safe to delete this copy
   if (0 != remove (tar_file))
     throw ApplicationError, "*** Error: removing $tar_file"$;

   return archive_dest_subdir;
}

define process_tar_file (tar_file, archive_incoming_dir,
                         archive_level_dir)
{
   if (-1 == access (tar_file, F_OK | R_OK))
     {
        vmessage ("*** Error: cannot access tar file $tar_file"$);
        return -1;
     }

   % This process is usually running on a compute node and the
   % archive incoming directory is probably on a different host.
   % Therefore, we copy the tar file to the remote archive incoming
   % directory before unpacking it. Once unpacked, we delete the
   % remote archive's copy of the tar file, and once all that has
   % succeeded, we delete the local copy of the tar file.

   variable e, tar_file_cpy;
   try (e)
     {
        tar_file_cpy = copy_file (tar_file, archive_incoming_dir);
        () = unpack_and_archive (tar_file_cpy, archive_level_dir);
     }
   catch AnyError:
     {
        vmessage (e.message);
        return -1;
     }

   if (Delete_Tarfiles == 0)
     {
        vmessage ("Skipped remove $tar_file"$);
     }
   else if (0 != remove (tar_file))
     {
        vmessage ("*** Error: removing %s (%s)", tar_file,
                  errno_string(errno));
        return -1;
     }

   return 0;
}

private define is_directory (file)
{
   variable st = stat_file (file);
   if (st == NULL) return 0;
   return stat_is ("dir", st.st_mode);
}

define slsh_main ()
{
   Delete_Tarfiles = 0;
   variable archive_root_dir = getenv ("SDPC_ARCHIVE_DIR");
   variable archive_level = NULL;

   variable c = cmdopt_new (&error_routine);
   c.add ("h|help", &usage);
   c.add ("a|archive_root_dir", &archive_root_dir; type="string");
   c.add ("c|clobber", &Clobber_Output_Files; inc);
   c.add ("d|delete", &Delete_Tarfiles; inc);
   c.add ("l|level", &archive_level; type="string");
   variable __i = c.process (__argv, 1);

   if (__argc - __i < 1)
     usage();

   if (archive_root_dir == NULL)
     throw ApplicationError,
     "*** Error: Archive root directory not specified (SDPC_ARCHIVE_DIR not set)";

   Archive_Root_Dir = archive_root_dir;

   % The tar file basename must be prefixed by the granule name.
   % The tar file must unpack into a directory with the granule
   % name, and must contain the archive_subdir file at the top
   % level:
   %      $granule_name/${Archive_Subdir_File}

   variable archive_level_dir = path_concat (archive_root_dir, archive_level);
   variable archive_incoming_dir = path_concat (archive_level_dir, "incoming");

   if (NULL == stat_file (archive_level_dir))
     {
        throw ApplicationError, "*** Error: cannot stat $archive_level_dir"$;
     }

   variable path_list = __argv[[__i:]];

   variable tar_file_list;
   if ((__argc - __i == 1)
       && is_directory (path_list[0]))
     {
        tar_file_list = glob (sprintf ("%s/*.tar", path_list[0]));
     }
   else tar_file_list = path_list;

   variable tar_file, status;
   foreach tar_file (tar_file_list)
     {
        status = process_tar_file (tar_file, archive_incoming_dir,
                                   archive_level_dir);
        if (status != 0) exit(1);
     }

   return 0;
}
