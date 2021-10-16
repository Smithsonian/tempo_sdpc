#*****************************************************************************
# BEGIN_FILE_PROLOG:
#
# FILENAME:
#       OMI_52090.t
#
# DESCRIPTION:
#       This file contains PGS_SMF standard return code definitions for the
#       OMI_MET group of tools.
#
#       The file is intended to be used as input by the smfcompile utility,
#       which generates the message file, and the C and FORTRAN header files.
#
#       usage   : $PGSBIN/smfcompile -f OMIE_52090.t -r -i
# AUTHOR:
#       Ellyne Kinney / SSAI
#
# HISTORY:
#       08-APR-2003  EKK  Initial version
#
# END_FILE_PROLOG:
#*****************************************************************************

%INSTR  = OMI
%LABEL  = OMI
%SEED   = 52090

# General
OMI_S_SUCCESS            PGE finishes normally with non-fatal errors, exit code = 0
OMI_E_FATAL              Execution aborted due to a fatal error, exit code = 1

OMI_W_GENERAL            General warning message
OMI_E_GENERAL            General error message

# Memory allocation
OMI_E_MEM_ALLOC          Memory allocation error
OMI_E_MEM_REF            Memory reference error   

# HDF and HDF-EOS
OMI_E_HDF                HDF failure
OMI_E_HDFEOS             HDF-EOS failure

# File I/O
OMI_E_FILE_OPEN          Failed to open file
OMI_E_FILE_READ          Failed to read file
OMI_E_FILE_WRITE         Failed to write file
OMI_E_FILE_CLOSE         Failed to close file

# PCF interface
OMI_E_PCF_REFID          Failed to get Reference ID
OMI_E_PCF_SPLIT          Failed to split Reference ID

# Swath interface
OMI_E_SWATH_OPEN         Open swath failed
OMI_E_SWATH_CREATE       Create swath failed
OMI_E_SWATH_ATTACH       Attach swath failed
OMI_E_SWATH_DEF_DIM      Define swath geo fields failed
OMI_E_SWATH_DEF_GEO      Define swath data fields failed
OMI_E_SWATH_DEF_DATA     Define swath dimensions failed
OMI_E_SWATH_WRITE_DATA   Write swath data failed
OMI_E_SWATH_WRITE_ATTR   Write swath attribute failed
OMI_E_SWATH_DETTACH      Detach swath failed
OMI_E_SWATH_CLOSE        Close swath failed

# Metadata Interface
OMI_E_MET_INIT           Metadata initalization failure
OMI_E_MET_SD_START       Metadata science data interface start failure
OMI_E_MET_READ           Metadata read failure
OMI_E_MET_WRITE          Metadata write failure
OMI_E_MET_SD_END         Metadata science data interface end failure
OMI_E_MET_REMOVE         Metadata remove failure

OMI_E_INVALID_CHANNEL    Invalid channel identifier
OMI_E_VALID_RANGE        Variable out of valid ranges 
OMI_E_TYPE_INVALID       Invalid input type
OMI_E_TYPE_NOT_SET       Input type not set
OMI_E_INPUT              Error in input data
OMI_E_TOO_MANY_PIXELS    Too many pixels in image
OMI_E_GRAN_PROC_LEVEL    Processing Level for Granule ID Invalid
OMI_E_GRAN_ID            Local Granule ID could not be set

# Informational messages
OMI_S_MET_READ           Read Metadata successfully
OMI_S_GRAN_ID            Local Granule ID set successfully
