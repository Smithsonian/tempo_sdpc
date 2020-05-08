#! /bin/sh
#SBATCH --output=/dev/null

set -e
set -u

if test $# -ne 1 ; then
   echo "Usage: $0 <granule-arch-dir-path>"
   exit 1
fi

#. ../ctrl/sdpc_setup.sh
export PATH="$SDPC_ROOT/bin:$PATH"

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

  # insert fixed metadata and, make sure the merged product
  # gets added to the product registry
  if test -f $product_file ; then
     insert_fixed_metadata.py $o3p_dir/$product_file
     ln -s $o3p_dir/$product_file $SDPC_ARCHIVE_DIR/registry/incoming
  fi

  block_tarfile_name=$(basename $product_file .nc)
  tar cf ${block_tarfile_name}.tar --remove-files --transform s,\^,${block_tarfile_name}/, block_*
}

perform_merge "$1"
