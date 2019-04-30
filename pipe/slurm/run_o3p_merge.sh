#! /bin/sh
#SBATCH --output=/dev/null

set -e
set -u

if test $# -ne 2 ; then
   echo "Usage: $0 path <L2-incoming-tar-file> OR"
   echo "Usage: $0 merge <granule-arch-dir-path>"
   exit 1
fi

#. ../ctrl/sdpc_setup.sh
export PATH="$SDPC_ROOT/bin:$PATH"

get_granule_arch_dir_path()
{
  tarfile_path="$1"
  tarfile_basename=$(basename $tarfile_path .tar)

  granule_ident=$(tar -x -f $tarfile_path --to-stdout $tarfile_basename/granule_ident.csv | tr ',' '=')

  #echo "$granule_ident"
  eval "$granule_ident"

  arch_type="$SDPC_ARCHIVE_DIR/L2/$processing_version/$product_type"
  ymd_scan="$tstart_year/$tstart_month/$tstart_mday/$scan_num"

  granule_arch_dir_path="${arch_type}/${ymd_scan}/${granule_num}"

  printf "$granule_arch_dir_path"
}

perform_merge()
{
  granule_arch_dir_path="$1"

  o3p_dir="$granule_arch_dir_path/o3p"

  if ! test -d "$o3p_dir" ; then
      printf "*** Error: directory not found: $o3p_dir"
      exit 1
  fi

  cd $o3p_dir
  #echo "chdir $o3p_dir"

  input_files=$(ls block_*/TEMPO_o3p*.nc)

  # The blocks all have the same filename.
  # Use the same filename for the merged result.
  first_file=$(echo $input_files | cut -d' ' -f 1)
  product_file=$(basename $first_file)

  num_files=$(echo $input_files | wc -w)
  quoted_input_files=$(echo $input_files | sed -e "s/[^ ]*/\'&\'/g")

cat << EOF > merge_o3p_iolist.nml
&merge_o3p_iolist
ninput=$num_files
input_files=$quoted_input_files
outfile='$product_file'
EOF

  srun --ntasks=1 --output=merge.log merge_o3p_files

  # make sure the merged product gets added to the product registry
  if test -f $product_file ; then
     ln -s $o3p_dir/$product_file $SDPC_ARCHIVE_DIR/registry/incoming
  fi
}

mode="$1"

case "$mode" in
  path )
   get_granule_arch_dir_path "$2"
  ;;

  merge )
   perform_merge "$2"
  ;;

  * )
  echo "*** Error: unrecognized mode = $mode"
  exit 1
  ;;
esac
