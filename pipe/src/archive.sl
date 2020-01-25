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
private variable Granule_Ident_File = "granule_ident.csv";

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
   % Extract the ${Granule_Ident_File} file.
   %
   % To avoid granule-specific directory name collisions when
   % multiple processes simultaneously archive data products,
   % we dump the data via stdout to a unique temporary file
   % Tar options are chosen to ensure that no existing files are
   % overwritten

   variable dir = path_dirname (tar_file);
   variable tar_file_basename = path_basename (tar_file);

   variable name_fields = strtok (tar_file_basename, ".");
   variable granule_name = name_fields[0];

   variable temp_granule_ident_file =
     path_concat (dir, tar_file_basename + "_" + Granule_Ident_File);

   variable argv = ["tar", "-x", "-f", tar_file, "--to-stdout",
                    path_concat (granule_name, Granule_Ident_File)];
   variable s = new_process(argv; stdout=temp_granule_ident_file).wait();
   if (s.exit_status != 0)
     {
        throw ApplicationError,
          "*** Error: extracting ${Granule_Ident_File} from tar archive: $tar_file"$;
     }

   % Determine the destination archive directory from the
   % contents of the granule_ident CSV file.
   variable _g = csv_readcol (temp_granule_ident_file);
   variable g = @Struct_Type (_g.col1);
   array_map (Void_Type, &set_struct_field, g, _g.col1, _g.col2);

   % Now delete the temporary file
   if (0 != remove (temp_granule_ident_file))
     throw ApplicationError, "*** Error removing $temp_granule_ident_file"$;

   variable sat_day = atof(g.sat_local_day_start);
   variable sat_day_dir = "D" + string(int(sat_day));

   % Construct the archive subdirectory path:
   variable subdir_seq = ["V" + g.processing_version,
                          g.product_type,
                          sat_day_dir];

   % Radiances will have assigned scan/granule number values,
   % dark, and irradiance will not.
   % If anything else comes through, the granule name should be
   % unique for the given day.
   if (is_substr(g.product_type, "RAD"))
     {
        subdir_seq = [subdir_seq,
                      "S" + g.scan_num,
                      "G" + g.granule_num];
     }
   else
     {
        subdir_seq = [subdir_seq, granule_name];
     }

   variable dest_subdir = strjoin (subdir_seq, "/");

   return dest_subdir;
}

define register_using_symlink (tar_file, archive_dest_subdir)
{
   % Dump partial paths to archived data products into a temporary
   % file on the local machine (usually a compute node), ideally on
   % a RAM disk.
   % Some tar files contain block_??? subdirectories with .nc files
   % that will eventually be merged to generate the final product.
   % These block .nc files should not be entered in the product
   % database, so we filter them out of this query.
   variable tmpfile_dir = "/var/tmp/$USER/$SDPC_PIPE_NAME"$;
   () = mkdir_p (tmpfile_dir);
   variable tmpfile = sprintf ("%s/register_symlink.%d", tmpfile_dir, getpid());
   variable argv = ["tar", "tf", tar_file,
		    "--exclude=block_*", "--strip-components=1",
                    "--show-transformed-names", "--wildcards",
                    "--no-anchored", "TEMPO_*.nc"];
   %vmessage (strjoin (argv, " "));
   variable s = new_process(argv; stdout=tmpfile, dup2=1).wait();
   if (s.exit_status != 0)
     {
        throw ApplicationError, "*** Error: creating file:$tmpfile"$;
     }

   % Read the partial paths from the temporary file.
   % "partial_path" means something like
   %      HCHO/TEMPO_HCHO_L2_V01_20130715T165956Z_S002G01.nc
   % (and remember that this temporary file may not contain anything
   % relevant).
   variable partial_paths;
   variable fp = fopen (tmpfile, "r");
   if (fp == NULL)
     throw ApplicationError, "*** Error: reading file:$tmpfile"$;
   partial_paths = fgetslines (fp);
   () = fclose(fp);
   () = remove (tmpfile);
   partial_paths = array_map (String_Type, &strtrim, partial_paths, "\n");

   % For each product file, make a symbolic link in $incoming_dir
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
        oldpath = path_concat (archive_dest_subdir, pp);
	if (NULL == stat_file (oldpath))
	  continue;
        newpath = path_concat (incoming_dir, path_basename(pp));
        if (0 != symlink (oldpath, newpath))
          throw ApplicationError, "*** Error: creating symlink $newpath"$;
     }

   return 0;
}

define unpack_and_archive (tar_file, archive_dest_dir)
{
   % The unpack and archive process involves the following steps:
   %  1. Extract the ${Granule_Ident_File} file
   %  2. Use ${Granule_Ident_File} to determine the destination directory
   %  3. Unpack the contents in the appropriate destination
   %     directory, avoiding file collisions and taking care
   %     to ensure consistency between ${Granule_Ident_File} files.

   variable subdir = get_tarfile_archive_subdir (tar_file);
   variable archive_dest_subdir = path_concat (archive_dest_dir, subdir);

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

   if (NULL != stat_file (path_concat (archive_dest_subdir, Granule_Ident_File)))
     {
        argv = [argv, "--exclude=${Granule_Ident_File}"$];
        % At this point, we could also check that the two
        % ${Granule_Ident_File} files are identical...
     }
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
                         archive_dest_dir)
{
   if (-1 == access (tar_file, F_OK | R_OK))
     {
        vmessage ("*** Error: cannot access tar file $tar_file"$);
        return -1;
     }

   variable e, tar_file_cpy;
   try (e)
     {
        tar_file_cpy = copy_file (tar_file, archive_incoming_dir);
        () = unpack_and_archive (tar_file_cpy, archive_dest_dir);
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
   variable sdpc_root_dir = getenv ("SDPC_ROOT");
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

   if (sdpc_root_dir == NULL)
     throw ApplicationError, "*** Error: SDPC_ROOT is not set";

   if (archive_root_dir == NULL)
     throw ApplicationError,
     "*** Error: Archive root directory not specified (SDPC_ARCHIVE_DIR not set)";

   Archive_Root_Dir = archive_root_dir;

   % The tar file basename must be prefixed by the granule name.
   % The tar file must unpack into a directory with the granule
   % name, and must contain the granule ident file at the top
   % level:
   %      $granule_name/${Granule_Ident_File}

   variable archive_dest_dir = path_concat (archive_root_dir, archive_level);
   variable archive_incoming_dir = path_concat (archive_dest_dir, "incoming");

   if (NULL == stat_file (archive_dest_dir))
     {
        throw ApplicationError, "*** Error: cannot stat $archive_dest_dir"$;
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
                                   archive_dest_dir);
        if (status != 0) exit(1);
     }

   return 0;
}
