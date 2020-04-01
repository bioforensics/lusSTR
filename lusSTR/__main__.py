#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import lusSTR
import sys


def main(args=None):  # pragma: no cover
    if args is None:  # pragma: no cover
        if len(sys.argv) == 1:
            lusSTR.cli.get_parser().parse_args(['-h'])
        args = lusSTR.cli.get_parser().parse_args()
    mainmethod = lusSTR.cli.mains[args.subcmd]
    mainmethod(args)
