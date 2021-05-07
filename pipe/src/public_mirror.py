#! /usr/bin/env python3

import os, sys
import signal
import time
import subprocess
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

    sig = Signal_Catcher()

    argv = ["update_public_mirror.sh"]

    while not sig.caught():
        obj = subprocess.run (argv)
        sig.wait(240)

    print ("Exiting: caught signal = {}".format(sig.signum))

if __name__ == "__main__":
    main()
