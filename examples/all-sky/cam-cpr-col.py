#! /usr/bin/env python3

# This script compares SW fluxes for 2 cam_test runs, column by column.

import sys

import numpy as np
import xarray as xr

if (len(sys.argv) != 3):
    sys.stderr.write('Usage: cam-cpr.py file1 file2\n')
    sys.exit(1)

case1 = sys.argv[1]
case2 = sys.argv[2]

f1 = xr.open_dataset(case1)
f2 = xr.open_dataset(case2)

# the flux fields have dimensions (plev_rad,lon)
ncols = len(f1['lon'])

for v in ['sw_flux_up', 'sw_flux_dn', 'sw_flux_dir']:

    print('Comparing columns for %s' % (v))

    for i in range(ncols):

        v1 = f1[v].data[:,i]
        v2 = f2[v].data[:,i]
        abs_diff = abs(v1 - v2)
        avg = 0.5 * (v1 + v2)
        abs_rel_diff = abs_diff / avg
        print(
            'Column %d: max abs difference: %e; '
            'max abs rel difference: %e' % (
                i, abs_diff.max(), abs_rel_diff.max()))
