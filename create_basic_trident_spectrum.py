import yt
import trident
import numpy as np

from astropy.io import fits

## read in example data
fn = 'WindTest/DD0010/DD0010'
ds = yt.load(fn)

## create ray start and end
ray_start = [1,0,1]
ray_end = [1,2,1]


## lines to calculate spectrum fo
## line_list = ['H','C','N','O','Mg','S','Si','Ne']
ldb = trident.LineDatabase('lines.txt')
line_list = ['H I 1216', 'C IV 1548', 'O VI 1032']
ll = ldb.parse_subset(line_list)


## begin making fits header
prihdr = fits.Header()
prihdr['AUTHOR'] = "Molly Peeples"
prihdr['DATE'] = time.strftime("%c") ## doesn't have time zone
prihdr['RAY_START'] = str(ray_start)
prihdr['RAY_END'] = str(ray_end)
prihdr['SIMULATION_NAME'] = fn
i = 1
for line in ll:
    keyword = 'LINE_'+str(i)
    prihdr[keyword] = line.name
    i += 1
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

    # let's plot it!
    filespecout = 'spectrum_'+ds.basename+'_'+line.identifier.replace(" ", "_")+'.png'
    sg.plot_spectrum(filespecout,flux_limits=(0.0,1.0))

    # make this sg an hdu
    col1 = fits.Column(name='wavelength', format='E', array=sg.lambda_field)
    col2 = fits.Column(name='tau', format='E', array=sg.tau_field)
    col3 = fits.Column(name='flux', format='E', array=sg.flux_field)
    col_list = [col1,col2,col3]

    for key in sg.line_observables[line.identifier].keys():
        col = fits.Column(name='sim_'+key,format='E',array=sg.line_observables[line.identifier][key])
        col_list = np.append(col_list,col)    
 
    cols = fits.ColDefs(col_list)
    sghdr = fits.Header()
    sghdr['LINE_NAME'] = line.identifier
    sghdr['LINE_REST_WAVELENGTH'] = line.wavelength
    sghdr['LINE_F_VALUE'] = line.f_value
    sghdr['LINE_GAMMA'] = line.gamma
    sghdu = fits.BinTableHDU.from_columns(cols,header=sghdr)

    sghdulist.append(sghdu)



sghdulist.writeto('spec_lots.fits', overwrite=True)
