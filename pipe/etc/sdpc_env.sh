
# SDPC_ROOT = root directory where pipeline processing software is installed
export SDPC_ROOT="/soft/tempo/sdpc/install/v1_gnu/sdpc"

# SDPC_OTS_ROOT = root directory where OTS software is installed
export SDPC_OTS_ROOT="/soft/tempo/sdpc/install/v1_gnu/ots"

# SDPC_RUN_DIR = root directory where pipeline processing takes place on slave nodes
export SDPC_RUN_DIR="/scratch/sdpc_test/sdpc_run_dir"

# SDPC_RUN_DIR_MASTER = root directory where pipeline processing takes place on the master node
export SDPC_RUN_DIR_MASTER="/home/houck/sdpc_test"

# SDPC_ARCHIVE_DIR = root directory for the mission archive
export SDPC_ARCHIVE_DIR="/home/houck/sdpc_test/archive"

export PATH="$SDPC_ROOT/bin:$SDPC_OTS_ROOT/bin:$PATH"

