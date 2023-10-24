#!/usr/bin/env python

from __future__ import print_function

import argparse
import json

def process_timing_file(file):

    with open(file) as json_file:
        timing_dict = json.load(json_file)

    print("Min time (μs) Function")
    print("              (mean,       s.d.,  max, scaled by minimum)")
    for sub in timing_dict.keys():
        tmin = timing_dict[sub]['min']
        tmean = timing_dict[sub]['mean']
        tstd = timing_dict[sub]['std']
        tmax = timing_dict[sub]['max']

        fmn = (tmean-tmin)/tmin
        fmx = (tmax-tmin)/tmin
        fsd = tstd/tmin
        print("%13.4E %s" % (tmin*1e6, sub))
        print("              (mean=+%4.2f%% ±%4.2f%% max=+%3.1f%%)" % (fmn*100, fsd*100, fmx*100))

def compare_timing_files(file, ref):

    with open(file) as json_file:
        timing_dict = json.load(json_file)

    with open(ref) as json_file:
        ref_dict = json.load(json_file)

    print("Delta (μs)  %-change Function")
    print("                     new: tmin (mean,       s.d.,  max, scaled by minimum)")
    print("                     old: tmin (mean,       s.d.,  max, scaled by minimum)")
    for sub in {**timing_dict, **ref_dict}.keys():
        tmin = timing_dict[sub]['min']
        tmean = timing_dict[sub]['mean']
        tstd = timing_dict[sub]['std']
        tmax = timing_dict[sub]['max']

        fmn = (tmean-tmin)/tmin
        fmx = (tmax-tmin)/tmin
        fsd = tstd/tmin
        print("%13.4E %s" % (tmin*1e6, sub))
        print("              (mean=+%4.2f%% ±%4.2f%% max=+%3.1f%%)" % (fmn*100, fsd*100, fmx*100))

# Parse arguments
parser = argparse.ArgumentParser(
    description="Beautify timing output from MOM6 timing tests."
)
parser.add_argument(
    'file',
    help="File to process."
)
parser.add_argument(
    '-r', '--reference',
    help="Reference file to compare against."
)
args = parser.parse_args()

# Do the thing
if args.reference is None:
    process_timing_file(args.file)
else:
    compare_timing_files(args.file, args.reference)
