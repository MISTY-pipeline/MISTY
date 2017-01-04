import yt
import trident
import numpy as np

from astropy.io import fits

fn = 'WindTest/DD0010/DD0010'
ds = yt.load(fn)

ray_start = [1,0,1]
ray_end = [1,2,1]

## line_list = ['H','C','N','O','Mg','S','Si','Ne']
ldb = trident.LineDatabase('lines.txt')
line_list = ['H I 1216', 'C IV 1548', 'O VI 1032']
ll = ldb.parse_subset(line_list)

ray = trident.make_simple_ray(ds,start_position=ray_start,end_position=ray_end,
                                data_filename="ray.h5",lines=line_list,ftype='gas')


# sg = trident.SpectrumGenerator('COS-G130M')
sg = trident.SpectrumGenerator(lambda_min=1150, lambda_max=1250, dlambda=0.01)

sg.make_spectrum(ray,lines=line_list)

filespecout = 'spectrum_'+ds.basename+'.png'
sg.plot_spectrum(filespecout,flux_limits=(0.9,1.0))

sg.save_spectrum('spec_raw.txt')

prihdr = fits.Header()
prihdr['RAY_START'] = str(ray_start)
prihdr['RAY_END'] = str(ray_end)
prihdr['SIMULATION_NAME'] = fn
prihdu = fits.PrimaryHDU(header=prihdr)
sghdulist = fits.HDUList([prihdu])


## want separate sg's with different lines --> need different rays
ray = trident.make_simple_ray(ds,start_position=ray_start,end_position=ray_end,
                                data_filename="ray.h5",lines=line_list,ftype='gas')
ar = ray.all_data()
for line in ll:
    lambda_rest = line.wavelength
    lambda_min = lambda_rest * (1+min(ar['redshift_eff'])) - 1
    lambda_max = lambda_rest * (1+max(ar['redshift_eff'])) + 1
    print lambda_min, lambda_max
    sg = trident.SpectrumGenerator(lambda_min=lambda_min.value, lambda_max=lambda_max.value, dlambda=0.01)
    sg.make_spectrum(ray,lines=line.name)

    # make this sg an hdu
    col1 = fits.Column(name='wavelength', format='E', array=sg.lambda_field)
    col2 = fits.Column(name='tau', format='E', array=sg.tau_field)
    col3 = fits.Column(name='flux', format='E', array=sg.flux_field)
    cols = fits.ColDefs([col1, col2, col3])
    sghdr = fits.Header()
    sghdr['LINE'] = line.name
    sghdu = fits.BinTableHDU.from_columns(cols,header=sghdr)

    sghdulist.append(sghdu)



sghdulist.writeto('spec_lots.fits', overwrite=True)
