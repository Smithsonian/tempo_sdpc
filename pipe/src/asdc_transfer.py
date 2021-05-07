#! /usr/bin/env python3

import os, sys
import signal
import time
import subprocess
import argparse
from threading import Event

class Signal_Catcher:

  exit = None
  signum = None

  def __init__(self):
    self.exit = Event()
    signal.signal(signal.SIGINT, self.handler)
    signal.signal(signal.SIGHUP, self.handler)
    signal.signal(signal.SIGTERM, self.handler)

  def wait(self, delay):
      self.exit.wait(delay)

  def caught(self):
      return self.exit.is_set()

  def handler(self,signum, frame):
    self.exit.set()
    self.signum = signum

def main():
    parser = argparse.ArgumentParser(description='Manage ASDC interface')
    parser.add_argument('--wait', type=int, default=None,
                        help="time interval [sec] between updates")
    parser.add_argument('user_at_host', nargs='?', default=None,
                         metavar='USER@HOST', help="ASDC dropbox account")
    args = parser.parse_args()

    if args.wait:
        wait = args.wait
    else:
        wait = os.getenv ("SDPC_ASDC_INTERVAL")
        if wait == None:
            wait = 600
        else:
            wait = int(wait)

    wait = abs(wait)

    if args.user_at_host:
        user_at_host = args.user_at_host
    else:
        user_at_host = os.getenv ("SDPC_ASDC_DROPBOX")
        if user_at_host == None:
            print ("Exiting: SDPC_ASDC_DROPBOX not set")
            sys.exit(1)

    sig = Signal_Catcher()

    while not sig.caught():
        obj = subprocess.run (["asdc_pull.sh", user_at_host])
        obj = subprocess.run (["asdc_push.sh", user_at_host])
        sig.wait(wait)

    print ("Exiting: caught signal = {}".format(sig.signum))

if __name__ == "__main__":
    main()
