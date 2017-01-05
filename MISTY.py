import numpy as np
from astropy.io import fits
import time
import trident

ldb = trident.LineDatabase('lines.txt')


def write_header(ray,start_pos=None,end_pos=None,lines=None,author='NAME'):
	## begin making fits header
	prihdr = fits.Header()
	prihdr['AUTHOR'] = author
	prihdr['DATE'] = time.strftime("%c") ## doesn't have time zone
	prihdr['RAY_START'] = str(start_pos)
	prihdr['RAY_END'] = str(end_pos)
	prihdr['SIMULATION_NAME'] = ray.basename
	
	lines = ldb.parse_subset(lines)
	i = 1
	for line in lines:
	    keyword = 'LINE_'+str(i)
    	    prihdr[keyword] = line.name
    	    i += 1
	prihdu = fits.PrimaryHDU(header=prihdr)
	sghdulist = fits.HDUList([prihdu])
	return sghdulist

def generate_line(ray,line,write=False,hdulist=None):
    if (write == True) & ((type(hdulist) != fits.hdu.hdulist.HDUList)):
       raise ValueError('Must pass HDUList in order to write. Call write_header first.')

    if not isinstance(line,trident.Line):
    	line = ldb.parse_subset(line)
    	line = line[0]

    ar = ray.all_data()
    lambda_rest = line.wavelength
    lambda_min = lambda_rest * (1+min(ar['redshift_eff'])) - 1
    lambda_max = lambda_rest * (1+max(ar['redshift_eff'])) + 1

    sg = trident.SpectrumGenerator(lambda_min=lambda_min.value, lambda_max=lambda_max.value, dlambda=0.01)
    sg.make_spectrum(ray,lines=line.name)

    if write:
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

    	hdulist.append(sghdu)

    return sg

def write_out(hdulist,filename='spectrum.fits'):
	hdulist.writeto(filename, overwrite=True) 
	return




