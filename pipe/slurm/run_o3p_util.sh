#! /bin/sh -x
#SBATCH --output=/dev/null

# exit on error
set -e
# exit upon any usage of an undefined variable
set -u
ulimit -s unlimited

# 1. Processing will run in the subdirectory provided on the command line,
#    which already contains all necessary inputs.
# 2. Processing ultimately stores all results in a tar file in an
#    appropriate destination directory.
# 3. When processing ends, remove the processing directory.

num_blocks=0
range_spec="0 0 0 0"

mode="$1"

case $mode in
  init)
    num_blocks="$2"
    range_spec="$3"
    ;;
  cleanup)
    ;;
  *)
    echo "*** $0: unsupported mode = $mode"
    exit 1
    ;;
esac

work_dir="o3p"
run_dir=$(pwd)
parent_dir=$(dirname $run_dir)
cd $work_dir

l2_out_dir="$SDPC_RUN_DIR/L2/out"
l2_repro_dir="$SDPC_RUN_DIR/L2/repro"
etc_dir="$SDPC_ROOT/etc"

# get input file names
. ./files.lis
rad_file=$RAD
irr_file=$IRR
cld_file=$CLD

rad_basename=$(basename $rad_file .nc)
irr_basename=$(basename $irr_file .nc)

# Define product file name template
#
lev1_file_fmt=$(mkgranule_name -p %s ${rad_basename}.nc)
lev1_base_fmt=$(basename $lev1_file_fmt .nc)

out_basename=$(printf "$lev1_base_fmt" l2_o3p)

tar_product_to_dest_dir()
{
   dest_dir=$1

   cd $run_dir

   /bin/rm $work_dir/${rad_basename}.nc
   /bin/rm $work_dir/${irr_basename}.nc
   /bin/rm $work_dir/$cld_file

   work_dir_tarfile="${rad_basename}.${work_dir}.tar"
   granule_dir=$(basename $run_dir)

   /bin/mkdir -p $dest_dir
   /bin/mv $work_dir/granule_ident.csv .
   cd $parent_dir
   tar c --remove-files -f $dest_dir/.$work_dir_tarfile \
         $granule_dir/granule_ident.csv $granule_dir/$work_dir
   /bin/mv $dest_dir/.$work_dir_tarfile $dest_dir/$work_dir_tarfile
}

decide_cleanup_dest_dir()
{
   dirs="$(ls -d block_[0-9]*)"
   cleanup_dest_dir="$l2_out_dir"
   for d in $dirs ; do
     status_file="$d/exit_status"
     if ! test -r $status_file ; then
        cleanup_dest_dir="$l2_repro_dir"
        break
     fi
     s=$(cat $status_file)
     if test x"$s" != x0 ; then
        cleanup_dest_dir="$l2_repro_dir"
        break
     fi
   done
}

perform_merge()
{
   product_file="${out_basename}.nc"
   input_files=$(ls block_*/$product_file)
   num_files=$(echo $input_files | wc -w)
   quoted_input_files=$(echo $input_files | sed -e "s/[^ ]*/\'&\'/g")

cat << EOF > merge_o3p_iolist.nml
&merge_o3p_iolist
ninput=$num_files
input_files=$quoted_input_files
outfile='$product_file'
EOF

  merge_o3p_files > merge.log 2>&1

  /bin/rm $input_files
}

case $mode in
  cleanup)
    decide_cleanup_dest_dir
    if test "$cleanup_dest_dir" != "$l2_out_dir" ; then
        tar_product_to_dest_dir "$cleanup_dest_dir"
    else
        perform_merge
        tar_product_to_dest_dir "$cleanup_dest_dir"
        tarfile_path="$l2_out_dir/${rad_basename}.o3p.tar"
        archive.sl --delete -a $SDPC_ARCHIVE_DIR -l L2 $tarfile_path
    fi
    exit
    ;;
  *)
    ;;
esac

#--------------------------------------------------
#  The code below is used only when mode != cleanup
#--------------------------------------------------

finish()
{
   tar_product_to_dest_dir "$l2_repro_dir"
}
trap finish EXIT ERR

error_exit(){
  echo $1 >&2
  exit 1
}

define_blocking()
{
   # Get sizes as separate values:
   # mirror step = [min_ms:max_ms]
   #      xtrack = [min_xt:max_xt]
read min_ms max_ms min_xt max_xt << EOF
   $range_spec
EOF
   size_ms=$(($max_ms - $min_ms + 1))
   size_xt=$(($max_xt - $min_xt + 1))

   # For now, default bin factors are mirror_step=2, xtrack=4
   # File to be broken up along north-south axis,
   # check xtrack size is divisible by block size
   if test $(($size_xt % $num_blocks)) -ne 0 ; then
     error_exit "size_xt must be divisible by block"
   fi
   size_xt_block=$(($size_xt / $num_blocks))
}

create_subdir()
{
   blk=$1
   subdir_name=$2

   k=$(($blk - 1))
   beg=$(($min_xt + $k \* $size_xt_block))
   end=$(($beg + $size_xt_block - 1))

   if test $end -gt $max_xt ; then
      end=$max_xt
   fi

   mkdir $subdir_name

cat << EOF > "$subdir_name/block.txt"
   xt_beg=$beg
   xt_end=$end
   ms_beg=$min_ms
   ms_end=$max_ms
EOF
}

config_subdir()
{
   subdir_name=$1

   # Copy control file to product directory
   control_file="${etc_dir}/o3_profile/default_main_control.inp"
   /bin/cp $control_file $subdir_name

   # Load default config parameters
   config_file="$SDPC_ROOT/etc/o3_profile/o3_profile.rc"
   . $config_file

   # file names
   radiance_file="../${rad_basename}.nc"
   irradiance_file="../${irr_basename}.nc"
   cloud_file="../${cld_file}"
   product_file="${out_basename}.nc"
   control_file_basename=$(basename $control_file)

   # Read the block parameters:
   . ./$subdir_name/block.txt

   template_pcf="$etc_dir/o3_profile/default.pcf.in"
# Edit the PCF file template:
   sed \
    -e s,@refdata_dir@,$refdata_dir,g \
    -e s,@spectra_dir@,$spectra_dir,g \
    -e s,@cloud_dir@,$cloud_dir,g \
    -e s,@product_dir@,$product_dir,g \
    -e s,@radiance_file@,$radiance_file,g \
    -e s,@irradiance_file@,$irradiance_file,g \
    -e s,@cloud_file@,$cloud_file,g \
    -e s,@product_file@,$product_file,g \
    -e s,@control_file@,$control_file_basename,g \
    -e s,@line_sample_extent@,$ms_beg\ $ms_end\ $xt_beg\ $xt_end,g \
    $template_pcf > $subdir_name/$pcf_file
}

define_blocking

for blk in $(seq $num_blocks) ; do
   subdir_name="block_$blk"
   create_subdir $blk $subdir_name
   config_subdir $subdir_name
done
trap - EXIT
