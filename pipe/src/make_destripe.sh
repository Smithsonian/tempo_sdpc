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

make_path_for_scan_destripe_file()
{
   first_path="$1"

   product_type="$(global_attribute.py --attr product_type $first_path)"
   destripe_filename="$(basename $first_path | sed -E -e s,${product_type},DSTR${product_type}, -e 's,G[0-9]+,,')"

   product_dir=$(dirname $first_path)
   granule_dir=$(dirname $product_dir)
   destripe_dir="$(dirname $granule_dir)/destripe/$product_type"

   echo "$destripe_dir/$destripe_filename"
}

make_path_for_day_destripe_file()
{
   first_path="$1"

   product_type="$(global_attribute.py --attr product_type $first_path)"
   destripe_filename="$(basename $first_path | sed -E -e s,${product_type},DSTR${product_type}, -e 's,T[0-9]+Z_S[0-9]+G[0-9]+,,')"

   day_dir=$(level1_info -l $first_path)
   destripe_dir="$SDPC_ARCHIVE_DIR/L2/RAD/$day_dir/destripe/$product_type"

   echo "$destripe_dir/$destripe_filename"
}

destripe_scan()
{
   l2_paths="$1"

   # If any file has been destriped previously, silently do nothing.
   for p in $l2_paths ; do
       res="$(variable_exists.py --var /support_data/destriping_correction $p)"
       if test x"$res" == x"yes" ; then
          return
       fi
   done

   first_path=$(echo $l2_paths | cut -d' ' -f1)
   product_type="$(global_attribute.py --attr product_type $first_path)"
   processing_version="$(global_attribute.py --attr processing_version $first_path)"

   destripe_path="$(make_path_for_scan_destripe_file $first_path)"
   destripe_dir=$(dirname $destripe_path)
   if ! test -d $destripe_dir ; then
      mkdir -p $destripe_dir
   fi

   config_file="$destripe_dir/make_destripe.${product_type}.yml"
   log_file="$destripe_dir/make_destripe_${product_type}.log"

   l2_yaml_list=$(print_yaml_list "$l2_paths")

   # edit the control file template
   sed -e s,'@LEVEL2_PRODUCT_PATHS@',"$l2_yaml_list", \
       -e s,'@DESTRIPE_FILE_PATH@',"$destripe_path", \
       -e s,'@SDPC_DSTR_VERSION@',"$processing_version", \
       $SDPC_ROOT/etc/trace_gas/make_destripe.${product_type}.yml.in > $config_file

   # Generate the destriping correction file
   # (We don't register this file in the archive
   #  because it won't be used for anything else)
   make_destripe.py $config_file > $log_file 2>&1 || md_error_exit "make_destripe.py failed (see $log_file)" $LINENO

   # Apply destriping correction
   apply_log="$destripe_dir/destripe_${product_type}.log"
   destripe.py --corrfile "$destripe_path" $l2_paths > $apply_log 2>&1 || md_error_exit "destripe.py failed (see $apply_log)" $LINENO
}

make_day_destripe_file()
{
   l2_path_list_file="$1"

   first_path=$(head -1 $l2_path_list_file)
   product_type="$(global_attribute.py --attr product_type $first_path)"
   processing_version="$(global_attribute.py --attr processing_version $first_path)"

   destripe_path="$(make_path_for_day_destripe_file $first_path)"
   destripe_dir=$(dirname $destripe_path)
   if ! test -d $destripe_dir ; then
      mkdir -p $destripe_dir
   fi

   config_file="$destripe_dir/make_destripe.${product_type}.yml"
   log_file="$destripe_dir/make_destripe.${product_type}.log"

   l2_yaml_list=$(print_yaml_list "$(cat $l2_path_list_file)")

   # edit the control file template
   sed -e s,'@LEVEL2_PRODUCT_PATHS@',"$l2_yaml_list", \
       -e s,'@DESTRIPE_FILE_PATH@',"$destripe_path", \
       -e s,'@SDPC_DSTR_VERSION@',"$processing_version", \
       $SDPC_ROOT/etc/trace_gas/make_destripe.${product_type}.yml.in > $config_file

   # Generate the destriping correction file
   make_destripe.py $config_file > $log_file 2>&1 || md_error_exit "make_destripe.py failed (see $log_file)" $LINENO

   # Register the file in the sqlite database.
   if test -f $destripe_path ; then
      ln -s $destripe_path $SDPC_ARCHIVE_DIR/registry/incoming
   fi
}
