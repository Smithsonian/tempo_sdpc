# This file is intended to contain only function definitions,
# so running it alone doesn't do anything.
# The functions are intended to be called from a parent script
# which imports these definitions.

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
md_error_exit()
{
  if test "$1" != "0" ; then
    echo "*** ${PROGNAME}: Error: line $2: $1"
    exit 1
  fi
}
trap 'md_error_exit $? $LINENO' EXIT

make_path_for_scan_bkgcorr_file()
{
   first_path="$1"

   product_type="$(global_attribute.py --attr product_type $first_path)"
   bkgcorr_filename="$(basename $first_path | sed -E -e s,${product_type},BACK${product_type}, -e 's,G[0-9]+,,')"

   product_dir=$(dirname $first_path)
   granule_dir=$(dirname $product_dir)
   bkgcorr_dir="$(dirname $granule_dir)/bkgcorr/$product_type"

   echo "$bkgcorr_dir/$bkgcorr_filename"
}

make_path_for_day_bkgcorr_file() # currently unused
{
   first_path="$1"

   product_type="$(global_attribute.py --attr product_type $first_path)"
   bkgcorr_filename="$(basename $first_path | sed -E -e s,${product_type},BACK${product_type}, -e 's,T[0-9]+Z_S[0-9]+G[0-9]+,,')"

   day_dir=$(level1_info -l $first_path)
   bkgcorr_dir="$SDPC_ARCHIVE_DIR/L2/RAD/$day_dir/bkgcorr/$product_type"

   echo "$bkgcorr_dir/$bkgcorr_filename"
}

bkgcorr_scan()
{
   l2_paths="$1"

   # If any file has been background corrected previously, silently do nothing.
   for p in $l2_paths ; do
       res="$(variable_exists.py --var /support_data/background_correction $p)"
       if test x"$res" == x"yes" ; then
          return
       fi
   done

   first_path=$(echo $l2_paths | cut -d' ' -f1)
   processing_version="$(global_attribute.py --attr processing_version $first_path)"

   bkgcorr_path="$(make_path_for_scan_bkgcorr_file $first_path)"
   bkgcorr_dir=$(dirname $bkgcorr_path)
   if ! test -d $bkgcorr_dir ; then
      mkdir -p $bkgcorr_dir
   fi

   config_file="$bkgcorr_dir/make_background.yml"
   log_file="$bkgcorr_dir/make_background.log"

   l2_yaml_list=$(print_yaml_list "$l2_paths")

   # edit the control file template
   sed -e s,'@LEVEL2_PRODUCT_PATHS@',"$l2_yaml_list", \
       -e s,'@BKGCORR_FILE_PATH@',"$bkgcorr_path", \
       -e s,'@SDPC_BKGCORR_VERSION@',"$processing_version", \
       $SDPC_ROOT/etc/trace_gas/make_background.yml.in > $config_file

   # Generate the background correction file
   # (We don't register this file in the archive
   #  because it won't be used for anything else)
   make_background.py $config_file > $log_file 2>&1 || md_error_exit "make_background.py failed (see $log_file)" $LINENO

   # Apply background correction
   apply_log="$bkgcorr_dir/bkgcorr.log"
   background.py --corrfile "$bkgcorr_path" $l2_paths > $apply_log 2>&1 || md_error_exit "background.py failed (see $apply_log)" $LINENO
}

make_day_bkgcorr_file() # currently unused
{
   l2_path_list_file="$1"

   first_path=$(head -1 $l2_path_list_file)
   processing_version="$(global_attribute.py --attr processing_version $first_path)"

   bkgcorr_path="$(make_path_for_day_bkgcorr_file $first_path)"
   bkgcorr_dir=$(dirname $bkgcorr_path)
   if ! test -d $bkgcorr_dir ; then
      mkdir -p $bkgcorr_dir
   fi

   config_file="$bkgcorr_dir/make_background.yml"
   log_file="$bkgcorr_dir/make_background.log"

   l2_yaml_list=$(print_yaml_list "$(cat $l2_path_list_file)")

   # edit the control file template
   sed -e s,'@LEVEL2_PRODUCT_PATHS@',"$l2_yaml_list", \
       -e s,'@BKGCORR_FILE_PATH@',"$bkgcorr_path", \
       -e s,'@SDPC_BKGCORR_VERSION@',"$processing_version", \
       $SDPC_ROOT/etc/trace_gas/make_background.yml.in > $config_file

   # Generate the background correction file
   make_background.py $config_file > $log_file 2>&1 || md_error_exit "make_background.py failed (see $log_file)" $LINENO

   # Register the file in the sqlite database.
   ln -s $bkgcorr_path $SDPC_ARCHIVE_DIR/registry/incoming
}
