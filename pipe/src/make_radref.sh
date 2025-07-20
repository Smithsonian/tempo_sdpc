# This file is intended to contain only function definitions,
# so running it alone doesn't do anything.
# The functions are intended to be called from a parent script
# which imports these definitions.

: "${SDPC_ARCHIVE_DBFILE:?SDPC_ARCHIVE_DBFILE not set}"
: "${SDPC_ARCHIVE_DBFILE_L1:?SDPC_ARCHIVE_DBFILE_L1 not set}"

: "${SDPC_RADREF_VERSION:=1}"

# Format the YAML file lists
# The result will have a trailing newline, but
# AFAIK, blank lines in yaml files are ok.
print_yaml_list()
{
   files="$1"
   lines=""
   for f in $files ; do
      line="  - $f\n"
      lines="${lines}${line}"
   done
   echo "$lines"
}

PROGNAME="$(basename $0)"
mr_error_exit()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error: line $2: $1"
  fi
}
trap 'mr_error_exit $? on $LINENO' EXIT

trigger_level2()
{
  radref_basename="$1"
  rad_files="$2"

  level2_input_dir="$SDPC_PIPE_DIR/stage/granules/level2_input"
  radref_wait_dir="$level2_input_dir/radref_pending"
  if ! test -d $radref_wait_dir ; then
     return
  fi

  for rf in $rad_files ; do
      basename_sans_extname=$(basename $rf .nc)
      tar_hostfile_path="$radref_wait_dir/${basename_sans_extname}_radref.tar"
      if test -f $tar_hostfile_path ; then
         echo "radref_file=$radref_basename" >> $tar_hostfile_path
         /bin/mv $tar_hostfile_path $level2_input_dir
      fi
  done
}

make_radref()
{
   l2_cloud_paths="$1"

   # Get the L1 radiance path for each cloud product.
   # If we were certain RAD_L1 and CLDO4_L2 were in the same sqlite dbfile, we could
   # use an SQL inner join to get the matching list of radiance files in one SQL query.
   # Unfortunately, during reprocessing, RAD_L1 and CLDO4_L2 may be in different dbfiles,
   # so two queries are needed, one to get scan_id from the CLDO4_L2 entry, and one to
   # get all the RAD_L1 files for that scan_id.  The input CLDO4_L2 path list should be
   # time-ordered, so as long as the RAD_L1 paths are time-ordered, the two lists should
   # match up.
   first_cloud_path=$(echo $l2_cloud_paths | cut -d' ' -f1)
   first_cloud_basename="$(basename $first_cloud_path)"
   scan_id=$(sqlite3 -readonly -cmd ".timeout 10000" $SDPC_ARCHIVE_DBFILE "select scan_id from CLDO4_L2 where filename = \"$first_cloud_basename\"")
   rad_files=$(sqlite3 -readonly -cmd ".timeout 10000" $SDPC_ARCHIVE_DBFILE_L1 "select path from RAD_L1 where scan_id = $scan_id order by istart")

   # Get the start/end time for each granule
   tbeg_lis=""
   tend_lis=""
   for rpath in $rad_files ; do
      tstart=$(global_attribute.py --attr time_coverage_start_since_epoch $rpath)
      tend=$(global_attribute.py --attr time_coverage_end_since_epoch $rpath)
      tbeg_lis="${tbeg_lis}${tstart}\n"
      tend_lis="${tend_lis}${tend}\n"
   done

   l2_yaml_list=$(print_yaml_list "$l2_cloud_paths")
   rad_yaml_list=$(print_yaml_list "$rad_files")

   first_radiance_path=$(echo $rad_files | cut -d' ' -f1)
   first_radiance_filename_sans_extname=$(basename $first_radiance_path .nc)

   # When reprocessing, the radiance file may have come from $SDPC_ARCHIVE_DBFILE_L1,
   # which does not point to $SDPC_ARCHIVE_DIR.  To be sure that we put the radref into
   # $SDPC_ARCHIVE_DIR, we construct such a directory path explicitly:
   archive_subdir_for_granule=$(level1_info --dir $first_radiance_path)
   radref_dir="$SDPC_ARCHIVE_DIR/L1/$(dirname $archive_subdir_for_granule)/radref"

   if ! test -d $radref_dir ; then
      mkdir -p $radref_dir
   fi

   # make radref filename
   ymd=$(global_attribute.py --attr time_coverage_start $first_radiance_path | cut -dT -f1 | tr -d '-')
   tbeg="$(echo $tbeg_lis | sort -n | head -1 | cut -d'.' -f1)"
   tend="$(echo $tend_lis | sort -n | tail -1 | cut -d'.' -f1)"
   scan_label="$(echo $first_radiance_filename_sans_extname | cut -d_ -f6 | cut -dG -f1)"
   #version="$(echo $first_radiance_filename_sans_extname | cut -d_ -f4)"
   version="$(printf V%02d $SDPC_RADREF_VERSION)"
   radref_filename="TEMPO_RADREF_L1_${version}_${ymd}_S${tbeg}_E${tend}_${scan_label}.nc"

   radref_path="$radref_dir/$radref_filename"
   config_file="$radref_dir/make_radref.yml"
   log_file="$radref_dir/make_radref.log"

   # edit the control file template
   sed -e s,'@RADIANCE_FILE_PATHS@',"$rad_yaml_list", \
       -e s,'@CLOUD_FILE_PATHS@',"$l2_yaml_list", \
       -e s,'@RADREF_PATH@',"$radref_path", \
       -e s,'@SDPC_RADREF_VERSION@',"$SDPC_RADREF_VERSION", \
       $SDPC_ROOT/etc/make_radref.yml.in > $config_file

   # Generate the radiance reference file
   make_radref.py $config_file > $log_file 2>&1 || mr_error_exit "make_radref.py failed (see $log_file)" $LINENO

   # Register the file in the sqlite database.
   ln -s $radref_path $SDPC_ARCHIVE_DIR/registry/incoming

   # Trigger Level 2 processing of any granules that were waiting for this.
   if test $SDPC_RADREF_CRON_TRIGGER_LEVEL2 -eq 0 ; then
      trigger_level2 "$radref_filename" "$rad_files"
   fi
}
