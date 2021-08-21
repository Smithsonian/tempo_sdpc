#! /bin/sh
#SBATCH --output=/dev/null

set -e
set -u

if test $# -ne 1 ; then
   echo "Usage: $0 <granule-arch-dir-path>"
   exit 1
fi

error_exit()
{
   printf "*** Error: o3prof_merge.sh: $1"
   exit 1
}

public_mirror_symlink()
{
   src=$1

   mirror_dir="$SDPC_RUN_DIR_MASTER/public_mirror"

   # Do nothing when the mirror directory is absent
   if ! test -d $mirror_dir ; then
      return 0
   fi

   bn=$(basename $src)
   day=$(level1_info --localday $src)
   target_dir="$mirror_dir/$day/L2"
   if ! test -d $target_dir ; then
      mkdir -p $target_dir
   fi
   ln -s $src $target_dir/$bn || error_exit "public_mirror_symlink failed: $bn"
}

perform_merge()
{
  granule_arch_dir_path="$1"

  o3p_dir="$granule_arch_dir_path/O3PROF"

  if ! test -d "$o3p_dir" ; then
      printf "*** Error: directory not found: $o3p_dir"
      exit 1
  fi

  cd $o3p_dir
  #echo "chdir $o3p_dir"

  input_files=$(ls block_*/TEMPO_O3PROF*.nc)

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

  # The .met files should be equivalent, so any one will do
  first_input_file=$(echo $input_files | cut -d' ' -f1)
  /bin/cp ${first_input_file}.met .

  met_file=$(basename ${first_input_file}.met)
  fix_met_format.py $met_file

  # insert fixed metadata and, make sure the merged product
  # gets added to the product registry, and the public mirror
  if test -f $product_file ; then
     insert_fixed_metadata.py $o3p_dir/$product_file
     ln -s $o3p_dir/$product_file $SDPC_ARCHIVE_DIR/registry/incoming
     public_mirror_symlink $o3p_dir/$product_file
  fi

  block_tarfile_name=$(basename $product_file .nc)
  tar cf ${block_tarfile_name}.tar --remove-files --transform s,\^,${block_tarfile_name}/, block_*
}

perform_merge "$1"
