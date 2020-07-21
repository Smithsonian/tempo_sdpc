#! /usr/bin/env python3

import os, sys
import signal
import time
import subprocess

class Signal_Catcher:
  kill_now = False
  signum = None
  def __init__(self):
    signal.signal(signal.SIGINT, self.exit_gracefully)
    signal.signal(signal.SIGHUP, self.exit_gracefully)
    signal.signal(signal.SIGTERM, self.exit_gracefully)

  def exit_gracefully(self,signum, frame):
    self.kill_now = True
    self.signum = signum

def main():

    sig = Signal_Catcher()

    argv = ["update_public_mirror.sh"]

    while not sig.kill_now:
        obj = subprocess.run (argv)
        time.sleep (240)

    print ("Exiting: caught signal = {}".format(sig.signum))

if __name__ == "__main__":
    main()
