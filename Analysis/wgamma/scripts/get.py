#!/usr/bin/env python
import sys

vals = {
    "lumi_2016": "35.9 fb^{-1} (13 TeV, 2016)",
    "lumi_2017": "41.5 fb^{-1} (13 TeV, 2017)",
    "lumi_2018": "59.5 fb^{-1} (13 TeV, 2018)",
}

print vals[sys.argv[1]],
