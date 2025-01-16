#! /usr/bin/env python3

import sys
import os
import re
import argparse
import datetime
import shutil
import glob
import signal

from threading import Event
from subprocess import run

def log_entry (s):
    print("{}: handle_inr_crash.py: {}".format (datetime.datetime.now().isoformat(), s), flush=True)

class Signal_Catcher:

  exit = None
  signum = None

  def __init__(self):
    self.exit = Event()
    signal.signal(signal.SIGINT, self.handler)
    signal.signal(signal.SIGTERM, self.handler)

  def wait(self, delay):
      self.exit.wait(delay)

  def caught(self):
      return self.exit.is_set()

  def handler(self,signum, frame):
    self.exit.set()
    self.signum = signum

def basename_timestamp_field (path):
    regex = r'TEMPO_RAD[T]?_L1_V\d{2}_(\d{8}T\d{6}Z)_S\d{3}G\d{2}'
    return re.search (regex, os.path.basename(path)).group(1)

def sorted_granule_list (dir, globexpr):
    path_list = glob.glob(os.path.join (dir, globexpr))
    path_list.sort (key = basename_timestamp_field)
    return path_list

def redirect_failed_scan_granules (inr_input_dir, problem_dir):
    # Get a list of granules in the INR input cache
    inr_input_files = sorted_granule_list (inr_input_dir, "TEMPO_RAD*_L1*.nc")
    if len(inr_input_files) == 0:
        return -1

    # The crash must have occurred in the first scan.
    # Parse the scan number, keeping it as a string with leading zeros
    scan_num_regex = r'S(\d{3})G\d{2}'
    scan_num = re.search (scan_num_regex, inr_input_files[0]).group(1)
    log_entry ("moving scan {} to {}".format(scan_num, problem_dir))

    sig = Signal_Catcher()

    # Move all granules from the problem scan into the problem directory
    # until a granule from a different scan is detected.
    while not sig.caught():
        granules = sorted_granule_list (inr_input_dir, "TEMPO_RAD*_L1*S???G??.nc")
        if len(granules) > 0:
            for path in granules:
                fname = os.path.basename (path)
                this_scan_num = re.search (scan_num_regex, fname).group(1)
                if this_scan_num == scan_num:
                    to_path = os.path.join (problem_dir, fname)
                    log_entry ("rename {} to {}".format(path, to_path))
                    os.rename (path, to_path)
                else:
                    log_entry ("new scan detected: {}".format(path))
                    return 0
        log_entry ("waiting for new scan {} granules".format(scan_num))
        sig.wait(60)

    print ("Exiting: caught signal = {}".format(sig.signum))
    return -1

def replace_inr_config_file (sdpc_pipe_dir):
    src_inr_config_path = os.path.join (sdpc_pipe_dir, "inr/config/TempoPipelineInterfaceSAO.conf.cold_start")
    dst_inr_config_path = os.path.join (sdpc_pipe_dir, "inr/config/TempoPipelineInterfaceSAO.conf")
    log_entry ("replacing INR config file with: {}".format(src_inr_config_path))
    shutil.copy (src_inr_config_path, dst_inr_config_path)

def main ():

    sdpc_pipe_dir = os.getenv ("SDPC_PIPE_DIR")
    if sdpc_pipe_dir is None:
        print ("*** SDPC_PIPE_DIR is not defined")
        sys.exit(1)
    if not os.path.isdir (sdpc_pipe_dir):
        print ("*** nonexistent directory: SDPC_PIPE_DIR={}".format(sdpc_pipe_dir))
        sys.exit(1)

    inr_input_dir = os.path.join (sdpc_pipe_dir, "stage/granules/inr_input")
    problem_dir = os.path.join (inr_input_dir, "problem")

    # Create the problem subdirectory if it doesn't already exist
    os.makedirs(problem_dir, exist_ok=True)

    # Shut down the INR service
    log_entry("stopping INR service")
    run (["sdpcsrv.sh", "down", "inr"], check=True)

    # Replace INR config file to force a cold restart
    # [obsolete as of INR SW r3.0.0]
    # replace_inr_config_file (sdpc_pipe_dir)

    # Until granules from the next scan arrive,
    # redirect granules from failed scan into problem subdirectory
    status = redirect_failed_scan_granules (inr_input_dir, problem_dir)
    if status != 0:
        sys.exit(1)

    # Restart the INR service
    log_entry("starting INR service")
    #inr_service_dir = os.path.join (sdpc_pipe_dir, "services/inr")
    #run (["s6-svdt-clear", inr_service_dir], check=True)
    run (["sdpcsrv.sh", "up", "inr"], check=True)

if __name__ == '__main__':
    main()
