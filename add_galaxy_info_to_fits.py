'''

Reads in a MISTY fits file and adds information about the galaxy to the header

'''

from __future__ import print_function

import glob
import os
import argparse

os.sys.path.insert(0, '/Users/molly/Dropbox/foggie/foggie')

import numpy as np

from astropy.io import fits
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii

def add_galaxy_info_to_fits(filename, haloinfo):
    hdu = fits.open(filename)
    zsnap = hdu[0].header['REDSHIFT']

    # need to find this redshift in the halo_info file
    t = ascii.read(haloinfo, format='fixed_width')
    thisid = t['redshift'] ==  zsnap
    print('adding physical information from \n', t[thisid])
    assert len(t[thisid]) == 1

    hdu[0].header['MVIR'] = (t['Mvir'][thisid][0], 'Msun')
    hdu[0].header['RVIR'] = (t['Rvir'][thisid][0], 'kpc')
    hdu[0].header['MSTAR'] = (t['Mstar'][thisid][0], 'Msun')
    hdu[0].header['MISM'] = (t['Mism'][thisid][0], 'Msun')
    hdu[0].header['SFR'] = (t['SFR'][thisid][0], 'Msun/yr')

    hdu.writeto(filename, overwrite=True, output_verify='fix')


if __name__ == "__main__":

    long_dataset_list = glob.glob(os.path.join(".", 'hlsp_misty_foggie_halo008508_nref11n_nref10f*v6_los.fits.gz'))
    haloinfo = '/Users/molly/Dropbox/foggie/foggie/halo_infos/008508/orig_nref11n/nref11n_nref10f/halo_info'
    dataset_list = long_dataset_list

    for filename in dataset_list:
        new_filename = filename
        print('adding galaxy info to ', filename)
        add_galaxy_info_to_fits(filename, haloinfo)
