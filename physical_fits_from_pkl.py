import glob
import os

os.sys.path.insert(0, '/Users/molly/Dropbox/foggie/foggie')

import numpy as np

import pickle

from astropy.io import fits
import astropy.units as u
from astropy.table import Table


### these are for the FOGGIE I ion list
ion_keys = ['H_p0_number_density', 'Si_p1_number_density', 'Si_p3_number_density', \
        'C_p3_number_density', 'O_p5_number_density']
ion_dict = { 'H_p0_number_density':'HI',
                "Si_p1_number_density":'SiII' ,
                "Si_p3_number_density":'SiIV' ,
                 'C_p3_number_density':'CIV',
                'O_p5_number_density':'OVI'}

def physical_fits_from_pkl(pkl_name, orig_fits_name, new_fits_name, **kwargs):
    orig_hdu = fits.open(orig_fits_name)
    new_hdu = fits.HDUList(orig_hdu)
    pkl_dict = pickle.load( open( pkl_file, "rb" ) )
    tbl = Table.from_pandas(pkl_dict['ray_df'])

    hdr = fits.Header()
    hdr['losfile'] = orig_fits_name
    phys_keys = ['density', 'temperature', 'metallicity']
    for ion in ion_dict:
        column = tbl[ion] * tbl['dx']
        sumcol = np.sum(column)
        density = np.sum(column * tbl['density']) / sumcol
        temperature = np.sum(column * tbl['temperature']) / sumcol
        metallicity = np.sum(column * tbl['metallicity']) / sumcol
        key = ion_dict[ion] + '_density'
        hdr[key] = (density, 'g cm^-3')
        key = ion_dict[ion] + '_temperature'
        hdr[key] = (temperature, 'K')
        key = ion_dict[ion] + '_metallicity'
        hdr[key] = (metallicity, 'Solar')


    col_list = []
    ## position
    col = fits.Column(name='los_position', format='E', array=tbl['x_ray'], unit='ckpc/h')
    col_list = np.append(col_list, col)

    ## velocity
    if 'axx' in pkl_name:
        velkey = 'x-velocity'
    elif 'axy' in pkl_name:
        velkey = 'y-velocity'
    elif 'axz' in pkl_name:
        velkey = 'z-velocity'
    col = fits.Column(name='los_velocity', format='E', array=tbl[velkey], unit='km s^-1')
    col_list = np.append(col_list, col)

    col = fits.Column(name='density', format='E', array=tbl['density'], unit='g cm^-3')
    col_list = np.append(col_list, col)
    col = fits.Column(name='temperature', format='E', array=tbl['temperature'], unit='K')
    col_list = np.append(col_list, col)
    col = fits.Column(name='metallicity', format='E', array=tbl['metallicity'], unit='Solar')
    col_list = np.append(col_list, col)
    for key in ion_keys:
        col = fits.Column(name=key, format='E', array=tbl[key], unit='cm^-3')
        col_list = np.append(col_list, col)
    cols = fits.ColDefs(col_list)

    hdu = fits.BinTableHDU.from_columns(cols, header=hdr)
    hdu.writeto(new_fits_name, overwrite=True, output_verify='fix')


if __name__ == "__main__":

    long_dataset_list = glob.glob(os.path.join(".", 'hlsp*v6_los.fits.gz'))
    ##long_dataset_list = ['./hlsp_misty_foggie_halo008508_nref11n_nref10f_rd0020_axy_dx042.9_dz082.4_v6_lsf.fits.gz']
    dataset_list = long_dataset_list

    for filename in dataset_list:
        pkl_file = '../' + filename.strip('_los.fits.gz').replace('.','') + '_sizes.pkl'
        out_fits_name =  '.' + filename.strip('_los.fits.gz') + '_physical.fits.gz'
        print('reading in pkl file ', pkl_file, ' and creating fits file ', out_fits_name)
        physical_fits_from_pkl(pkl_file, filename, out_fits_name)
