#! /usr/bin/env python3

# This script compares SW fluxes for 2 cam_test runs.

import numpy as np
import xarray as xr

case1='cam-profile.nc'
case2='cam-test.nc'
#case2='/work/test-src/rp-cam-test/examples/all-sky/cam-profile.nc'

f1 = xr.open_dataset(case1)
f2 = xr.open_dataset(case2)

for v in ['sw_flux_up', 'sw_flux_dn', 'sw_flux_dir']:

    v1 = f1[v].data[:,0]
    v2 = f2[v].data[:,0]
    abs_diff = abs(v1 - v2)
    avg = 0.5 * (v1 + v2)
    # Division raises a runtime warning when we divide by zero even if the
    # values in those locations will be ignored.
    with np.errstate(divide='ignore', invalid='ignore'):
        abs_rel_diff = np.abs(
            np.where((avg > 2. * np.finfo(float).eps), abs_diff / avg, 0))

        print(
            'Variable %s: max abs difference: %e; '
            'max abs rel difference: %e' % (
                v, abs_diff.max(), abs_rel_diff.max()))
