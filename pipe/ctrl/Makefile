# -*- sh -*-

default:	inr-targets pipe-targets

inr-targets:
	@echo "INR targets:  inr-install inr-feed clean"

#-----------------------------------------------------
# Feed data into pipeline

INR_FEED = $(SDPC_INRSW_ROOT)/bin/TempoMockPipeline

inr-feed:
	/bin/cp $(INR_MOCK_FEED_CFG) $(SDPC_INR_RUN_DIR)
	cd $(SDPC_INR_RUN_DIR) ; \
	exec $(INR_FEED) $(INR_MOCK_FEED_CFG) > feed.log 2>&1 &

#-----------------------------------------------------
# Initialize INR directory tree

inr-install:	inr-mkdirs inr-install-refdata

INR_DIRS = config logs scantailoring Output CloudProducts \
           Staging/Granules Staging/Left Staging/Right Staging/rsr

inr-mkdirs:
	mkdir -p $(SDPC_INR_RUN_DIR)
	for d in $(INR_DIRS) ; do \
	  mkdir -p $(SDPC_INR_RUN_DIR)/$$d ; \
	done

# INR_GOES_REFDATA_DIR contains:
#   GOES imagery
#   GOES rsr files
INR_GOES_REFDATA_DIR=/data/tempo/sdpc/inr/synth_granules_v2/FullGranules

inr-install-refdata:
	-ln -s -t $(SDPC_INR_RUN_DIR)/Staging/rsr $(INR_GOES_REFDATA_DIR)/rsr/*.txt
	-ln -s -t $(SDPC_INR_RUN_DIR)/Staging/Left $(INR_GOES_REFDATA_DIR)/Reference/Left/*.nc
	-ln -s -t $(SDPC_INR_RUN_DIR)/Staging/Right $(INR_GOES_REFDATA_DIR)/Reference/Right/*.nc
	# ln -s -t $(SDPC_INR_RUN_DIR)/CloudProducts $(CLOUD_DATA_DIR)/*.nc

#-----------------------------------------------------

include ./pipe_paths.sh

pipe-targets:
	@echo "Pipe targets:  pipe-setup-root pipe-setup-nodes delete-node-dirs"

pipe-setup-root:	pipe-setup-archive pipe-setup-processing-root

SDPC_ARCHIVE_DIR=$(PROCESSING_ROOT)/archive

pipe-setup-archive:
	-mkdir -p $(SDPC_ARCHIVE_DIR)/L0/incoming \
	          $(SDPC_ARCHIVE_DIR)/L1/incoming \
                  $(SDPC_ARCHIVE_DIR)/L2/incoming \
                  $(SDPC_ARCHIVE_DIR)/L3/incoming

pipe-setup-processing-root:
	-mkdir -p $(PROCESSING_ROOT)/L0/incoming/ccd \
                  $(PROCESSING_ROOT)/L0/incoming/telem \
                  $(PROCESSING_ROOT)/L0/repro \
                  $(PROCESSING_ROOT)/L1/irradiance \
                  $(PROCESSING_ROOT)/L1/repro \
	          $(PROCESSING_ROOT)/L2/incoming \
                  $(PROCESSING_ROOT)/L2/repro \
                  $(PROCESSING_ROOT/status \
                  $(NODE_ROOT_DIR)
	-ln -s -n $(ANCILLARY_DIR) $(PROCESSING_ROOT)/ancillary
	-ln -s -n $(IRRADIANCE_DIR) $(PROCESSING_ROOT)/L1/irradiance/out

pipe-setup-nodes:
	-srun -N$(NUM_NODES) mkdir -p $(NODE_ROOT_DIR)/L0/incoming \
	                              $(NODE_ROOT_DIR)/L0 \
	                              $(NODE_ROOT_DIR)/L1 \
	                              $(NODE_ROOT_DIR)/L2
	-srun -N$(NUM_NODES) ln -s -n $(REFDATA_DIR) $(NODE_ROOT_DIR)/refdata
	-srun -N$(NUM_NODES) ln -s -n $(PROCESSING_ROOT)/ancillary \
	                              $(NODE_ROOT_DIR)/ancillary
	-srun -N$(NUM_NODES) ln -s -n $(PROCESSING_ROOT)/L0/incoming/telem \
	                              $(NODE_ROOT_DIR)/L0/incoming/telem
	-srun -N$(NUM_NODES) ln -s -n $(PROCESSING_ROOT)/L0/repro \
	                              $(NODE_ROOT_DIR)/L0/repro
	-srun -N$(NUM_NODES) ln -s -n $(PROCESSING_ROOT)/L1/irradiance \
	                              $(NODE_ROOT_DIR)/L1/irradiance
	-srun -N$(NUM_NODES) ln -s -n $(PROCESSING_ROOT)/L1/repro \
	                              $(NODE_ROOT_DIR)/L1/repro
	-srun -N$(NUM_NODES) ln -s -n $(SDPC_INR_RUN_DIR)/Staging/Granules \
	                              $(NODE_ROOT_DIR)/L1/radiance_inr_staging
	-srun -N$(NUM_NODES) ln -s -n $(PROCESSING_ROOT)/L2/incoming \
	                              $(NODE_ROOT_DIR)/L2/incoming
	-srun -N$(NUM_NODES) ln -s -n $(PROCESSING_ROOT)/L2/repro \
	                              $(NODE_ROOT_DIR)/L2/repro

NODE_ROOT_DIR_SYMLINKS = \
        $(NODE_ROOT_DIR)/refdata \
        $(NODE_ROOT_DIR)/ancillary \
        $(NODE_ROOT_DIR)/L1/irradiance \
        $(NODE_ROOT_DIR)/L2/incoming
NODE_ROOT_DIR_SUBDIRS = \
        $(NODE_ROOT_DIR)/L0/out \
        $(NODE_ROOT_DIR)/L0 \
        $(NODE_ROOT_DIR)/L1/out \
        $(NODE_ROOT_DIR)/L1 \
        $(NODE_ROOT_DIR)/L2/out \
        $(NODE_ROOT_DIR)/L2/repro \
        $(NODE_ROOT_DIR)/L2

delete-node-dirs:
	-srun -N$(NUM_NODES) /bin/rm -f $(NODE_ROOT_DIR_SYMLINKS)
	-srun -N$(NUM_NODES) /bin/rmdir $(NODE_ROOT_DIR_SUBDIRS)
