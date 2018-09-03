
# The script assumes that the environment is properly set (e.g. PATH, SDPC_ROOT, etc.)

run_wavecal()
{
   input_file="$1"
   chunks="$2"

   bname=$(basename $input_file .nc)
   result_dir="wavecal_${bname}"
   mkdir $result_dir

   jid_list=$(sbatch -w $SLURMD_NODENAME --parsable --array=$chunks wavecal_block.sh band_290_490_nm $input_file $result_dir)

   # FIXME - vis band narrow cal window not yet specified
   #jid_vis=$(sbatch -w $SLURMD_NODENAME --parsable --array=$chunks wavecal_block.sh band_540_740_nm $input_file $result_dir)
   #jid_list=${jid_list}:${jid_vis}

   sbatch -w $SLURMD_NODENAME --wait --dependency=afterany:$jid_list \
          run_wavecal_merge.sh $input_file $result_dir
}
