'''

Reads in a MISTY fits file that may or may not have spectacle run on it yet;
runs spectacle and creates new fits file with new header information

'''

from __future__ import print_function

import glob
import os

import MISTY
from spectacle.analysis import Resample

os.sys.path.insert(0, '/Users/molly/Dropbox/foggie/foggie')
from plot_misty_spectra import plot_misty_spectra

import numpy as np

from astropy.io import fits
import astropy.units as u
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.table import Table

from scipy.signal import argrelextrema

def add_spectacle_to_fits(old_fits_name, new_fits_name, **kwargs):
    threshold = kwargs.get('threshold', 0.01)
    plot = kwargs.get('plot', False)

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
        print('~~~~> trying',line_name,'~~~~~>>>>')

        if any([x.name.upper() == line_name.upper() for x in orig_hdu]):
            new_ext = orig_hdu[line_name]
            for k in orig_hdu[line_name].header:
                if k not in keys_to_copy:
                    print("deleting ", k)
                    del new_ext.header[k]

            lambda_0 = orig_hdu[line_name].header['RESTWAVE']
            try:
                disp = orig_hdu[line_name].data['wavelength']
                flux = orig_hdu[line_name].data['flux']
                tau = orig_hdu[line_name].data['tau']
                redshift = orig_hdu[line_name].data['redshift']
            except:
                disp = orig_hdu[line_name].data['disp_obs']
                flux = orig_hdu[line_name].data['flux_obs']
                tau = orig_hdu[line_name].data['tau_obs']
                redshift = orig_hdu[line_name].data['redshift_obs']

            zsnap = np.median(redshift)

            ## we want Nmin
            Nmin = np.size(np.where(flux[argrelextrema(flux, np.less)[0]] < (1-threshold)))
            new_ext.header['Nmin'] = Nmin

            print("~~~~> now trying to run spectacle on line ",line_name, "~~~~~~>")
            lines_properties = MISTY.get_line_info(disp, flux, \
                                            tau=tau, \
                                            redshift=zsnap, \
                                            lambda_0=lambda_0, \
                                            f_value=orig_hdu[line_name].header['F_VALUE'], \
                                            gamma=orig_hdu[line_name].header['GAMMA'], \
                                            ion_name=line_name, \
                                            threshold = threshold)
            print(lines_properties)


            for line_key in lines_properties:
                if isinstance(lines_properties[line_key], tuple):
                    if np.isnan(lines_properties[line_key][0]):
                        lines_properties[line_key] = -99.
                new_ext.header[line_key] = lines_properties[line_key]


            new_hdu.append(new_ext)
            print('~~~~> all done with',line_name,'~~~~~<<<')
        else:
            print('<<<<<~~~~ ',line_name,' not found :-(  ~~~~~<<<<<<')

    print("writing out to .... " + new_fits_name)
    new_hdu.writeto(new_fits_name, overwrite=True, output_verify='fix')

    plotname = '.' + new_fits_name.strip('.fits.gz') + '.png'
    print('plotting to... ' + plotname)
    plot_misty_spectra(new_hdu, overplot=True, outname=plotname)



if __name__ == "__main__":

    long_dataset_list = glob.glob(os.path.join(".", 'hlsp*v4_lsf.fits.gz'))
    dataset_list = long_dataset_list

    for filename in dataset_list:
        new_filename = '.' + filename.strip('lsf.fits.gz') + 'lsf.fits.gz'
        print('adding spectacle to ', filename, ' and saving as ', new_filename)
        add_spectacle_to_fits(filename, new_filename, threshold=0.005)
