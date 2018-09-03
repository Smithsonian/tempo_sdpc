# -*- sh -*-

include ./inr_paths.sh

default:	inr-targets pipe-targets

inr-targets:
	@echo "INR targets:  inr-config inr-install inr-start inr-feed clean"

LOGGING_CONF_FILE  = inr_config/TempoLogging.conf
PIPELINE_CONF_FILE = inr_config/TempoPipelineInterfaceSAO.conf

INR_AUTOGEN_FILES = $(LOGGING_CONF_FILE) $(PIPELINE_CONF_FILE)

INR_DIRS = config logs scantailoring Output CloudProducts \
           Staging/Granules Staging/Left Staging/Right Staging/rsr

INR_FEED = $(INRSW_INSTALL_DIR)/bin/TempoMockPipeline
INR_PIPE = $(INRSW_INSTALL_DIR)/bin/PipelineInterface

inr-start:
	cd $(INR_PROCESSING_ROOT) ; \
	exec $(INR_PIPE) $(PIPELINE_CONF_FILE) > run.log 2>&1 &

inr-feed:
	/bin/cp $(INR_MOCK_FEED_CFG) $(INR_PROCESSING_ROOT)
	cd $(INR_PROCESSING_ROOT) ; \
	exec $(INR_FEED) $(INR_MOCK_FEED_CFG) > feed.log 2>&1 &

inr-config:
	# INR logging
	sed -e s,@INR_PROCESSING_ROOT@,$(INR_PROCESSING_ROOT),g \
	       $(LOGGING_CONF_FILE).in > $(LOGGING_CONF_FILE)
	# INR pipeline
	sed -e s,@INR_PROCESSING_ROOT@,$(INR_PROCESSING_ROOT),g \
	    -e s,@INR_IERS_DIR@,$(INR_IERS_DIR),g \
	    -e s,@INR_REFDATA_DIR@,$(INR_REFDATA_DIR),g \
	   $(PIPELINE_CONF_FILE).in > $(PIPELINE_CONF_FILE)

inr-install:	inr-mkdirs inr-install-config inr-install-refdata

inr-mkdirs:
	mkdir -p $(INR_PROCESSING_ROOT)
	for d in $(INR_DIRS) ; do \
	  mkdir -p $(INR_PROCESSING_ROOT)/$$d ; \
	done

inr-install-config:
	/bin/cp -a inr_config/*.conf $(INR_PROCESSING_ROOT)/config

inr-install-refdata:
	-ln -s -t $(INR_PROCESSING_ROOT)/Staging/rsr $(INR_GOES_REFDATA_DIR)/rsr/*.txt
	-ln -s -t $(INR_PROCESSING_ROOT)/Staging/Left $(INR_GOES_REFDATA_DIR)/Reference/Left/*.nc
	-ln -s -t $(INR_PROCESSING_ROOT)/Staging/Right $(INR_GOES_REFDATA_DIR)/Reference/Right/*.nc
	# ln -s -t $(INR_PROCESSING_ROOT)/CloudProducts $(CLOUD_DATA_DIR)/*.nc

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
	-srun -N$(NUM_NODES) ln -s -n $(INR_PROCESSING_ROOT)/Staging/Granules \
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

clean:
	/bin/rm -f $(INR_AUTOGEN_FILES) .##* inr_config/.##*
