# -*- sh -*-

IGNORE := $(shell sh -c ". ../sdpc_env.sh; env | grep SDPC_ | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
include makeenv

# INR_REFDATA_DIR contains:
#   ellipsoid altitude map
#   GNA polygon
#   JPL ephemeris
INR_REFDATA_DIR=$(SDPC_INRSW_ROOT)/config

# INR_IERS_DIR contains:
#   IERS bulletin A
INR_IERS_DIR=$(SDPC_INRSW_ROOT)/config

LOGGING_CONF_FILE = TempoLogging.conf
PIPELINE_CONF_FILE = TempoPipelineInterfaceSAO.conf

INR_AUTOGEN_FILES = $(LOGGING_CONF_FILE) $(PIPELINE_CONF_FILE)

inr-config:
	# INR logging
	sed -e s,@INR_PROCESSING_ROOT@,$(SDPC_INR_RUN_DIR),g \
	       $(LOGGING_CONF_FILE).in > $(LOGGING_CONF_FILE)
	# INR pipeline
	sed -e s,@INR_PROCESSING_ROOT@,$(SDPC_INR_RUN_DIR),g \
	    -e s,@INR_IERS_DIR@,$(INR_IERS_DIR),g \
	    -e s,@INR_REFDATA_DIR@,$(INR_REFDATA_DIR),g \
	   $(PIPELINE_CONF_FILE).in > $(PIPELINE_CONF_FILE)

clean:
	/bin/rm -f $(INR_AUTOGEN_FILES) makeenv
