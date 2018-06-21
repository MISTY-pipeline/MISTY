'''

Reads in a MISTY fits file that may or may not have spectacle run on it yet;
runs spectacle and creates new fits file with new header information

'''

from __future__ import print_function

import MISTY
from spectacle.analysis import Resample

import glob
import os
import argparse

os.sys.path.insert(0, '/Users/molly/Dropbox/foggie/foggie')
from plot_misty_spectra import plot_misty_spectra

import numpy as np

from astropy.io import fits
import astropy.units as u
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.table import Table

from scipy.signal import argrelextrema

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="resamples and convovles LSF to MISTY spectra; optionally runs spectacle")

    parser.add_argument('--lsf', metavar='lsf', type=float, action='store',
                        help='FWHM of LSF? if zero; skips. default: 7km/s')
    parser.set_defaults(lsf=7.)

    parser.add_argument('--resample', metavar='resample', type=float, action='store',
                        help='resampling of spectrum? in km/s. if zero; skips. default: 2km/s')
    parser.set_defaults(resample=2.)

    parser.add_argument('--spectacle', dest='spectacle', action='store_true',
                            help='run spectacle? default is yes')
    parser.set_defaults(spectacle=True)

    parser.add_argument('--appendix', dest='appendix',  type=str, action='store',
                            help='what appendix to add? default is lsf')
    parser.set_defaults(appendix='lsf')


    args = parser.parse_args()
    return args


def add_spectacle_to_fits(old_fits_name, new_fits_name, **kwargs):
    resample = kwargs.get('resample', 0.0)  ## km/s of pixels
    fwhm = kwargs.get('fwhm', 0.0)  ## km/s of lsf
    threshold = kwargs.get('threshold', 0.01)
    plot = kwargs.get('plot', False)
    use_spectacle = kwargs.get('use_spectacle', True)
    plotname = kwargs.get('plotname', 'temp.png')

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
            disp = orig_hdu[line_name].data['wavelength']
            flux = orig_hdu[line_name].data['flux']
            old_flux = flux
            tau = orig_hdu[line_name].data['tau']
            redshift = orig_hdu[line_name].data['redshift']
            data = Table([disp, redshift, flux, tau], names=('disp', 'redshift', 'flux', 'tau'))
            data.sort('disp')
            zsnap = np.median(redshift)
            velocity = ((data['disp'] / orig_hdu[line_name].header['RESTWAVE']) - 1 - zsnap) * 299792.458 * u.km / u.s
            vcen = np.median(velocity.value)
            vmin, vmax = vcen - 1200, vcen + 1200
            print('flux for original line ',line_name,': ',np.min(flux),np.max(flux),'with Nmin=',np.size(np.where(old_flux[argrelextrema(old_flux, np.less)[0]] < (1-threshold))))
            if resample > 0:
                print('Resampling at',resample,' km/s:')
                pixsize = resample
                print(vmax, vmin, pixsize)
                nstep = np.round((vmax - vmin)/pixsize)
                new_vel = np.linspace(vmin, vmax, nstep) * u.km / u.s
                new_flux = Resample(new_vel)(velocity, flux)
                new_tau = Resample(new_vel)(velocity, data['tau'])
                new_redshift = Resample(new_vel)(velocity, data['disp'] / lambda_0 - 1)
                new_disp = lambda_0 * (((new_vel.value / 299792.458)) + (1 + zsnap))
                disp = new_disp
                flux = new_flux
                redshift = new_redshift
                tau = new_tau
                print('resampled, Nmin now = ',np.size(np.where(flux[argrelextrema(flux, np.less)[0]] < (1-threshold))))
                print(len(old_flux), len(new_flux))
            if fwhm > 0:
                print('Adding LSF with FWHM = ',fwhm,'km/s:')
                if resample > 0:
                    dv = nstep
                else:
                    dv = (np.max(velocity) - np.min(velocity))/len(velocity) ## number of velocity pixels
                sigma= (fwhm/dv) / (2*np.sqrt(2.*np.log(2.)))
                print('dv = ', dv, 'fwhm = ', fwhm, 'sigma = ', sigma)
                gaussker = Gaussian1DKernel(stddev=sigma) ## sigma is in units of pixels
                new_flux = convolve(flux, gaussker)
                new_tau = convolve(tau, gaussker)
                tau = new_tau
                print(len(flux), len(new_flux))
                flux = new_flux


            if resample > 0 or fwhm > 0:
                print('lengths: ', len(redshift), len(disp), len(tau), len(flux))
                z_col = fits.Column(name='redshift_obs', format='E', array=redshift)
                wavelength = fits.Column(name='disp_obs', format='E',
                                         array=disp, unit='Angstrom')
                tau_col = fits.Column(name='tau_obs', format='E', array=tau)
                flux_col = fits.Column(name='flux_obs', format='E', array=flux)
                col_list = [z_col, wavelength, tau_col, flux_col]

            ## we want Nmin
            Nmin = np.size(np.where(new_flux[argrelextrema(new_flux, np.less)[0]] < (1-threshold)))
            new_ext.header['Nmin'] = Nmin

            print('flux for new line ',line_name,': ',np.min(flux),np.max(flux), 'with Nmin=',Nmin)

            if use_spectacle:
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
                    new_ext.header[line_key] = lines_properties[line_key]

            cols = fits.ColDefs(col_list)
            final_new_ext = fits.BinTableHDU.from_columns(cols, header=new_ext.header, name=line_name)

            new_hdu.append(final_new_ext)
            print('~~~~> all done with',line_name,'~~~~~<<<')
        else:
            print('<<<<<~~~~ ',line_name,' not found :-(  ~~~~~<<<<<<')

    print('plotting...')
    plot_misty_spectra(new_hdu, overplot=True, outname=plotname)

    print("writing out to .... " + new_fits_name)
    new_hdu.writeto(new_fits_name, overwrite=True, output_verify='fix')




if __name__ == "__main__":

    args = parse_args()

    long_dataset_list = glob.glob(os.path.join(".", 'hlsp*rd0020*v5_los.fits.gz'))
    dataset_list = long_dataset_list

    for filename in dataset_list:
        fileroot = '.' + filename.strip('los.fits.gz')
        new_filename = fileroot + args.appendix + '.fits.gz'
        plotname = fileroot + args.appendix + '.png'
        print('adding spectacle to ', filename, ' and saving as ', new_filename, ' and plotting to ', plotname)
        add_spectacle_to_fits(filename, new_filename, resample=args.resample, \
                        fwhm=args.lsf, use_spectacle=args.spectacle, plotname=plotname)
