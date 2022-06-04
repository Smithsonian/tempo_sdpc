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

make_destripe()
{
  l2_paths="$1"

   # Get timestamps from each data product
   tbeg_lis=""
   tend_lis=""
   for f in $l2_paths ; do
      b=$(basename $f)
      tstart=$(global_attribute.py --attr time_coverage_start_since_epoch $f)
      tend=$(global_attribute.py --attr time_coverage_end_since_epoch $f)
      tbeg_lis="${tbeg_lis}${tstart}\n"
      tend_lis="${tend_lis}${tend}\n"
   done

   l2_yaml_list=$(print_yaml_list "$l2_paths")

   # put destripe in destripe/$product_type subdir of granule scan directory
   first_path=$(echo $l2_paths | cut -d' ' -f1)
   first_filename_sans_extname=$(basename $first_path .nc)
   product_dir=$(dirname $first_path)
   granule_dir=$(dirname $product_dir)
   product_type="$(global_attribute.py --attr product_type $first_path)"
   destripe_dir="$(dirname $granule_dir)/destripe/$product_type"

   if ! test -d $destripe_dir ; then
      mkdir -p $destripe_dir
   fi

   # make destripe filename
   tbeg="$(echo $tbeg_lis | sort -n | head -1 | cut -d'.' -f1)"
   tend="$(echo $tend_lis | sort -n | tail -1 | cut -d'.' -f1)"
   scan_label="$(echo $first_filename_sans_extname | cut -d_ -f6 | cut -dG -f1)"
   version="$(echo $first_filename_sans_extname | cut -d_ -f4)"
   destripe_filename="TEMPO_DSTR${product_type}_L2_${version}_S${tbeg}_E${tend}_${scan_label}.nc"

   destripe_path="$destripe_dir/$destripe_filename"
   config_file="$destripe_dir/make_destripe.yml"
   log_file="$destripe_dir/make_destripe.log"

   # edit the control file template
   sed -e s,'@LEVEL2_PRODUCT_PATHS@',"$l2_yaml_list", \
       -e s,'@DESTRIPE_FILE_PATH@',"$destripe_path", \
       $SDPC_ROOT/etc/make_destripe.yml.in > $config_file

   # Generate the destriping correction file
   make_destripe.py $config_file > $log_file 2>&1 || md_error_exit "make_destripe.py failed (see $log_file)" $LINENO

   # Register the file in the sqlite database.
   ln -s $destripe_path $SDPC_ARCHIVE_DIR/registry/incoming
}
