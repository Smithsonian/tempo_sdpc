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

def main():
    parser = argparse.ArgumentParser(description='Manage ASDC interface')
    parser.add_argument('--wait', type=int, default=300,
                        help="time interval [sec] between updates")
    parser.add_argument('--script-fmt', default=None,
                        help="Push/pull script basename format, e.g. asdc_%%s.sh where %%s is push|pull")
    parser.add_argument('dest', default=None,
                        help="Destination string, e.g. bucket:dir or user@host:dir")
    if len(sys.argv) == 1:
        parser.print_usage(sys.stderr)
        sys.exit(0)
    args = parser.parse_args()

    wait = abs(args.wait)
    dest = args.dest

    pull_script = args.script_fmt % ("pull")
    push_script = args.script_fmt % ("push")

    sig = Signal_Catcher()

    logprint ("Started")

    while not sig.caught():
        obj = subprocess.run ([pull_script, dest])
        obj = subprocess.run ([push_script, dest])
        sig.wait(wait)

    logprint ("Exiting: caught signal = {}".format(sig.signum))

if __name__ == "__main__":
    main()
