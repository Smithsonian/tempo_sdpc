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
   ymd=$(global_attribute.py --attr time_coverage_start $first_path | cut -dT -f1 | tr -d '-')
   tbeg="$(echo $tbeg_lis | sort -n | head -1 | cut -d'.' -f1)"
   tend="$(echo $tend_lis | sort -n | tail -1 | cut -d'.' -f1)"
   scan_label="$(echo $first_filename_sans_extname | cut -d_ -f6 | cut -dG -f1)"
   version="$(echo $first_filename_sans_extname | cut -d_ -f4)"
   destripe_filename="TEMPO_DSTR${product_type}_L2_${version}_${ymd}_S${tbeg}_E${tend}_${scan_label}.nc"

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

   # Apply the destriping correction
   apply_destripe=$(config_setting destripe.HCHO.apply)
   if test $apply_destripe -ne 0 ; then
      # If destripe_search==True, and the search succeeded, then
      # the files have already been destriped, and could be marked with
      # asdc_status=new|pending|uploaded|accepted|problem, depending on
      # what happened afterward. Otherwise (if the search failed, or if
      # destripe_search==False) then destriping was delayed until now,
      # and the products will be marked with asdc_status='defer'=100.
      # Now that we've generated the necessary destriping correction,
      # destriping can proceed, but we apply the correction only to
      # products marked 'defer'.

      scan_num="$(global_attribute.py --attr scan_num $first_path)"
      sql="select HCHO_L2.path from RAD_L1 inner join HCHO_L2 on RAD_L1.istart = HCHO_L2.istart and RAD_L1.scan_num = $scan_num and HCHO_L2.asdc_status = 100 order by HCHO_L2.istart"
      needs_destripe=$(sqlite3 -cmd ".timeout 5000" $SDPC_ARCHIVE_DBFILE "$sql")

      if test -n "$needs_destripe" ; then
         apply_log="$destripe_dir/destripe.log"
         destripe.py --corrfile "$destripe_path" $needs_destripe > $apply_log 2>&1 || md_error_exit "destripe.py failed (see $apply_log)" $LINENO
         # Change asdc_status of HCHO_L2 products from 'defer' to 'new'
         tmpfile=$(mktemp)
         printf "%s\n" $needs_destripe > $tmpfile
         asdc_track_uploads.py --stat --set new $tmpfile || error_exit "asdc_track_uploads failed: changing HCHO_L2 asdc_status defer to new"
         printf "%s.met\n" $needs_destripe > $tmpfile
         asdc_track_uploads.py --set new $tmpfile || error_exit "asdc_track_uploads failed: changing HCHO_L2 met asdc_status defer to new"
         /bin/rm -f $tmpfile
      fi
   fi

}
