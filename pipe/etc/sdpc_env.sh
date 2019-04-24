#---------------------------
#  Installed software paths:
#---------------------------

# SDPC_ROOT = root directory where pipeline processing software is installed
export SDPC_ROOT="/soft/tempo/sdpc/install/v1_gnu/sdpc"

# SDPC_INRSW_ROOT = root directory where INR SW is installed
export SDPC_INRSW_ROOT="/soft/tempo/sdpc/install/v1_gnu/inr"

# SDPC_OTS_ROOT = root directory where OTS software is installed
export SDPC_OTS_ROOT="/soft/tempo/sdpc/install/v1_gnu/ots"

# SDPC_S6_ROOT = root directory where s6 supervisor is installed
export SDPC_S6_ROOT="/soft/tempo/sdpc/skarnet"

#---------------------------------
#  Installed reference data paths:
#---------------------------------

# ANCILLARY_DIR contains time-sensitive reference data,
# shared by all processing nodes
# (e.g. GOES data, meteorology data, snow & ice cover)
export SDPC_ANCILLARY_ROOT="/data/tempo/sdpc/ancillary"

# SDPC_REFDATA_DIR contains static reference data,
# shared by all processing nodes
# (e.g. reference spectra, cross-sections, climatologies, etc.)
export SDPC_REFDATA_DIR="/data/tempo/sdpc/refdata"

#-----------------------------------------
#  Environment-dependent pipeline defaults
#-----------------------------------------

: "${SDPC_PIPE_HOME:=$HOME}"
export SDPC_PIPE_HOME

: "${SDPC_PIPE_NAME:=sdpc_test}"
export SDPC_PIPE_NAME

: "${SDPC_PROCESSING_VERSION:=1}"
export SDPC_PROCESSING_VERSION

#---------------------------------
#  Live processing directory paths
#---------------------------------

# SDPC_RUN_DIR_MASTER = root directory where pipeline processing takes place on the master node
export SDPC_RUN_DIR_MASTER="${SDPC_PIPE_HOME}/${SDPC_PIPE_NAME}"

# SDPC_INR_RUN_DIR = root directory where INR processing takes place on the master node
export SDPC_INR_RUN_DIR="${SDPC_RUN_DIR_MASTER}/inr"

# SDPC_ARCHIVE_DIR = root directory for the mission archive
export SDPC_ARCHIVE_DIR="${SDPC_RUN_DIR_MASTER}/archive"

# SDPC_RUN_DIR = root directory where pipeline processing takes place on slave nodes
export SDPC_RUN_DIR="/scratch/${SDPC_PIPE_NAME}/sdpc_run_dir"

#---------------------------
#  Command search path, etc.
#---------------------------

export PATH="$SDPC_ROOT/bin:$SDPC_OTS_ROOT/bin:$PATH"
