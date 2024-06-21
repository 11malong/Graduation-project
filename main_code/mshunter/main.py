#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import argparse
import os
import sys

curpath = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(os.path.dirname(curpath))

from genotype import *
from qc import *
from paras import args_process
from global_dict import *

logger.info(" ".join(sys.argv))

def main():

    global_init()
    arg = args_process()

    if arg:
        parase = arg.parse_args()

        if parase.command == "genotype":
            genotype(parase)

        if parase.command == "qc":
            qc(parase)

if __name__ == "__main__":
    main()
