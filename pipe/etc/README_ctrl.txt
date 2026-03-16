Function of the control directory: $SDPC_PIPE_DIR/ctrl
======================================================
Existence of files in this directory is used to enable/disable
selected pipeline functions as follows:

FILENAME: disable
FUNCTION: When $SDPC_PIPE_DIR/ctrl/disable exists,
          all cachemon-based services are disabled according to
          the control strings defined in $SDPC_PIPE_DIR/etc/disable.conf.

FILENAME: disable-asdc
FUNCTION: When $SDPC_PIPE_DIR/ctrl/disable-asdc exists,
          the ASDC service(s) are disabled.

FILENAME: disable-asdc-transfer
FUNCTION: When $SDPC_PIPE_DIR/ctrl/disable-asdc-transfer exists,
          the ASDC service(s) prepare files for upload by
          generating file manifests, but no data files are
          uploaded or downloaded.

FILENAME: disable-destripe-CLDO4
FUNCTION: When $SDPC_PIPE_DIR/ctrl/disable-destripe-CLDO4 exists,
          destriping of CLDO4 data products is disabled

FILENAME: disable-nrt
FUNCTION: When $SDPC_PIPE_DIR/ctrl/disable-nrt exists, data products
          will not be delivered to the NRT processing pipeline, thus
          disabling generation of NRT products.

FILENAME: disable-radt-level3
FUNCTION: When $SDPC_PIPE_DIR/ctrl/disable-radt-level3 exists,
          post-processing of RADT_L1 granules is disabled.

FILENAME: save-inr-input
FUNCTION: When $SDPC_PIPE_DIR/ctrl/save-inr-input exists, the level1a
          service will keep a copy of INR input files by creating a
          hard link to each file in the 'save' subdirectory of the INR
          input cache.
