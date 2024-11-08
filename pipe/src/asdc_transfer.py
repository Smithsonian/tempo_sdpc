#! /usr/bin/env python3

import os, sys
import signal
import time
import subprocess
import argparse
from threading import Event

Prefix = "asdc:"

# python3 will provide file= redirection to stderr
# noninteractive stderr is line-buffered so use explicit flush.
def eprint(*args, **kwargs):
    print(Prefix, *args, file=sys.stderr, flush=True, **kwargs)

def logprint(*args, **kwargs):
    print(Prefix, *args, file=sys.stdout, flush=True, **kwargs)

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

class Enable_Switch:

    dir = os.getenv ("SDPC_PIPE_DIR")
    disable_file = "disable-asdc"
    disable_path = os.path.join (dir, disable_file)

    def __init__(self, name):
        self.name = name
        self.name_file = "{}-{}".format(self.disable_file, name)
        self.name_path = os.path.join (self.dir, self.name_file)

    def enable (self):
        if os.path.isfile (self.disable_path):
            logprint ("ASDC {} service disabled ({} exists)".format(self.name, self.disable_file))
            return False
        elif os.path.isfile (self.name_path):
            logprint ("ASDC {} service disabled ({} exists)".format(self.name, self.name_file))
            return False
        else:
            return True

def main():
    parser = argparse.ArgumentParser(description='Manage ASDC interface')
    parser.add_argument('--wait', type=int, default=300,
                        help="time interval [sec] between updates")
    parser.add_argument('user_at_host', default=None,
                         metavar='USER@HOST:dirpath', help="ASDC dropbox account")
    args = parser.parse_args()

    wait = abs(args.wait)
    user_at_host = args.user_at_host

    sig = Signal_Catcher()

    logprint ("Started")

    enable_s3 = Enable_Switch ("s3")
    enable_direct = Enable_Switch ("direct")

    while not sig.caught():
        if enable_s3.enable():
            obj = subprocess.run (["asdc_pull_s3.sh"])
            obj = subprocess.run (["asdc_push_s3.sh"])
        if enable_direct.enable():
            obj = subprocess.run (["asdc_pull.sh", user_at_host])
            obj = subprocess.run (["asdc_push.sh", user_at_host])
        sig.wait(wait)

    logprint ("Exiting: caught signal = {}".format(sig.signum))

if __name__ == "__main__":
    main()
