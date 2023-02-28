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
    print(Prefix, *args, file=sys.stdout, **kwargs)

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

def main():
    parser = argparse.ArgumentParser(description='Manage ASDC interface')
    parser.add_argument('--wait', type=int, default=300,
                        help="time interval [sec] between updates")
    parser.add_argument('user_at_host', default=None,
                         metavar='USER@HOST', help="ASDC dropbox account")
    args = parser.parse_args()

    wait = abs(args.wait)
    user_at_host = args.user_at_host

    sig = Signal_Catcher()

    logprint ("Started", flush=True)

    while not sig.caught():
        obj = subprocess.run (["asdc_pull.sh", user_at_host])
        obj = subprocess.run (["asdc_push.sh", user_at_host])
        sig.wait(wait)

    logprint ("Exiting: caught signal = {}".format(sig.signum))

if __name__ == "__main__":
    main()
