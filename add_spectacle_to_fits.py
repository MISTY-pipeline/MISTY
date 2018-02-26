'''

Reads in a MISTY fits file that may or may not have spectacle run on it yet;
runs spectacle and creates new fits file with new header information

'''

from __future__ import print_function

import MISTY

import numpy as np
from astropy.io import fits

def add_spectacle_to_fits(old_fits_name, new_fits_name):
    orig_hdu = fits.open(old_fits_name)

    new_hdu = fits.HDUList([orig_hdu[0]])
    new_hdu.append(orig_hdu[1])

    keys_to_copy = ('LINENAME',
                    'RESTWAVE',
                    'F_VALUE',
                    'GAMMA',
                    'SIM_TAU_HDENS',
                    'SIM_TAU_TEMP',
                    'SIM_TAU_METAL',
                    'TOT_COLUMN',
                    'EXTNAME',
                    'XTENSION',     ## hack-y since don't want to not delete these
                    'BITPIX',
                    'NAXIS',
                    'NAXIS1',
                    'NAXIS2',
                    'PCOUNT',
                    'GCOUNT',
                    'TFIELDS',
                    'TTYPE1',
                    'TFORM1',
                    'TTYPE2',
                    'TFORM2',
                    'TUNIT2',
                    'TTYPE3',
                    'TFORM3',
                    'TTYPE4',
                    'TFORM4',
                    'TTYPE5',
                    'TFORM5',
                    'TTYPE6',
                    'TFORM6',
                    'TTYPE7',
                    'TFORM7',
                    'TTYPE8',
                    'TFORM8',
                    'TTYPE9',
                    'TFORM9',
                    'TTYPE10',
                    'TFORM10')

    ## now for the individual lines
    nlines = np.int(orig_hdu[0].header['NLINES'])
    for line_num in np.arange(nlines):
        key = 'LINE_'+str(line_num+1)
        line_name = orig_hdu[0].header[key]
        new_ext = orig_hdu[line_name]
        for k in orig_hdu[line_name].header:
            if k not in keys_to_copy:
                print("deleting ", k)
                del new_ext.header[k]

        print("~~~~> now trying to run spectacle ~~~~~~>")
        zsnap = np.median(orig_hdu[line_name].data['redshift'])  ## hack
        lines_properties = MISTY.get_line_info(orig_hdu[line_name].data['wavelength'], \
                                        orig_hdu[line_name].data['flux'], \
                                        tau=orig_hdu[line_name].data['tau'], \
                                        redshift=zsnap, \
                                        lambda_0=orig_hdu[line_name].header['RESTWAVE'], \
                                        f_value=orig_hdu[line_name].header['F_VALUE'], \
                                        gamma=orig_hdu[line_name].header['GAMMA'], \
                                        ion_name=line_name)

        for line_key in lines_properties:
            new_ext.header[line_key] = lines_properties[line_key]

        new_hdu.append(new_ext)

    print("writing out to .... " + new_fits_name)
    new_hdu.writeto(new_fits_name, overwrite=True, output_verify='fix')
