import numpy as np
from astropy.io import fits
import time
import trident
import spectacle

from spectacle.core.spectra import Spectrum1D
from spectacle.modeling.models import Absorption1D
from spectacle.core.lines import Line
from astropy.modeling.fitting import LevMarLSQFitter

import getpass
import datetime
import os.path

ldb = trident.LineDatabase('lines.txt')


def write_header(ray,start_pos=None,end_pos=None,lines=None,**kwargs):
	## begin making fits header
    prihdr = fits.Header()
    prihdr['AUTHOR'] = kwargs.get("author",getpass.getuser())
    prihdr['DATE'] = datetime.datetime.now().isoformat()
    prihdr['RAYSTART'] = str(start_pos)
    prihdr['RAYEND'] = str(end_pos)
    prihdr['SIM_NAME'] = ray.basename
    prihdr['NLINES'] = str(len(np.array(lines)))
    prihdr['DOI'] = "doi.corlies2017.paper.thisistotesnotmadeup"
    prihdr['PAPER'] = "Corlies et al. (2017) ApJ, ###, ###"
    prihdr['EUVB'] = "HM12" ## probably shouldn't be hardcoded
    prihdr['IMPACT'] = kwargs.get("impact","undef")
    prihdr['ANGLE'] = kwargs.get("angle","undef")

    lines = ldb.parse_subset(lines)

    i = 1
    for line in lines:
        keyword = 'LINE_'+str(i)
        prihdr[keyword] = line.name
        i += 1
    prihdu = fits.PrimaryHDU(header=prihdr)
    sghdulist = fits.HDUList([prihdu])
    return sghdulist

def write_parameter_file(ds,filename=None,hdulist=None):
    if type(hdulist) != fits.hdu.hdulist.HDUList:
        raise ValueError('Must pass HDUList in order to write. Call write_header first.')

    # is a filename given? then use that
    if filename != None and os.path.isfile(filename):
        param_file = np.genfromtxt(filename,delimiter='=',dtype=str,autostrip=True)
        col1 = fits.Column(name='PARAMETERS',format='A50',array=param_file[:,0])
        col2 = fits.Column(name='VALUES',format='A50',array=param_file[:,1])
    else:
        #  use ds.parameters
        col1 = fits.Column(name='PARAMETERS',format='A50',array=ds.parameters.keys())
        col2 = fits.Column(name='VALUES',format='A50',array=[str(v) for v in ds.parameters.values()])

    col_list = [col1,col2]
    cols = fits.ColDefs(col_list)
    sghdr = fits.Header()
    ds.index
    sghdr['SIM_CODE'] = ds.dataset_type
    print "---> SIM_CODE set to ", ds.dataset_type, "if you don't like this, change it!"
    sghdr['COMPUTER'] = 'pleiades'
    print "---> ASSUMING PLEIADES FOR NOW BUT SHOULD BE PASSSSSSED IN"

    sghdu = fits.BinTableHDU.from_columns(cols,header=sghdr,name='PARAMS')
    hdulist.append(sghdu)

    return sghdu

def generate_line(ray,line,write=False,hdulist=None):
    '''
    input: a lightray and a line; writes info to extension of hdulist
    '''
    if (write == True) & ((type(hdulist) != fits.hdu.hdulist.HDUList)):
       raise ValueError('Must pass HDUList in order to write. Call write_header first.')

    if not isinstance(line,trident.Line):
        ldb = trident.LineDatabase('lines.txt')
        line_out = ldb.parse_subset(line)
        print line, line_out
        line_out = line_out[0]

    ar = ray.all_data()
    lambda_rest = line_out.wavelength
    lambda_min = lambda_rest * (1+min(ar['redshift_eff'])) - 1.5
    lambda_max = lambda_rest * (1+max(ar['redshift_eff'])) + 1.5

    sg = trident.SpectrumGenerator(lambda_min=lambda_min.value, \
        lambda_max=lambda_max.value, dlambda=0.0001)
    sg.make_spectrum(ray,lines=line_out.name)

    if write:
    	col1 = fits.Column(name='wavelength', format='E', array=sg.lambda_field,unit='Angstrom')
    	col2 = fits.Column(name='tau', format='E', array=sg.tau_field)
    	col3 = fits.Column(name='flux', format='E', array=sg.flux_field)
    	col_list = [col1,col2,col3]

    	for key in sg.line_observables[line_out.identifier].keys():
    	    col = fits.Column(name='sim_'+key,format='E',array=sg.line_observables[line_out.identifier][key])
    	    col_list = np.append(col_list,col)

    	cols = fits.ColDefs(col_list)
    	sghdr = fits.Header()
    	sghdr['LINENAME'] = line_out.identifier
        print "----->>>>using ", line_out.identifier, "as LINENAME, whereas ", line, " was passed. Change?"
    	sghdr['RESTWAVE'] = line_out.wavelength
    	sghdr['F_VALUE'] = line_out.f_value
    	sghdr['GAMMA'] = line_out.gamma


        ## want to leave blank spaces now for values that we're expecting to generate for MAST
        ## first let's add some spaces for the simulated, tau-weighted values!
        sghdr['SIM_TAU_HDENS'] = -9999.
        sghdr['SIM_TAU_TEMP'] = -9999.
        sghdr['SIM_TAU_METAL'] = -9999.
        sghdr['TOT_COLUMN'] = np.log10(np.sum(sg.line_observables[line_out.identifier]['column_density'].value))


        ## we're also going to want data from Nick's fitting code
        ## it's going to give values for all of it's components
        ## for now, let's give it five and assume that many are going to be empty
        sghdr['NCOMPONENTS'] = 5.
        # names = ['fitEW','fitcol','fitvcen','fitb','fitv90']
        line_properties = get_line_info(sg)
        for key in line_properties:
            stringin = key+str(0)
            if np.isnan(line_properties[key]):
                sghdr[stringin] = "NaN"
            else:
                sghdr[stringin] = line_properties[key]
        names = ['fitEW','fitcol','fitvcen','fitb']
        ncomponent_standard = 5
        j = 1
        while j < ncomponent_standard:
            for name in names:
                stringin = name+str(j)
                sghdr[stringin] = -9999.
            j = j + 1

    	sghdu = fits.BinTableHDU.from_columns(cols,header=sghdr,name=line_out.name)

    	hdulist.append(sghdu)

    return sg

def get_line_info(sg):
    '''
    runs spectacle on a trident spectrum object and returns absorber properties
    '''
    disp = sg.lambda_field
    flux = sg.flux_field
    spectrum = Spectrum1D(flux, dispersion=disp)
    line = Line(name=sg.line_list[0]['label'], \
            lambda_0=sg.line_list[0]['wavelength'].value, \
            f_value=sg.line_list[0]['f_value'], \
            gamma=sg.line_list[0]['gamma'], fixed={'lambda_0': False,
                                                 'f_value': True,
                                                 'gamma': True,
                                                 'v_doppler': False,
                                                 'column_density': False,
                                                })
    spec_mod = Absorption1D(lines=[line])

    # Create a fitter. The default fitting routine is a LevMarLSQ.
    fitter = LevMarLSQFitter()
    fit_spec_mod = fitter(spec_mod, spectrum.dispersion, spectrum.data, maxiter=500)
    fit_y = fit_spec_mod(spectrum.dispersion.value)

    # OK now we want line properties
    line_properties = {'fitcol' : fit_spec_mod.column_density_1.value,
                       'fitb': fit_spec_mod.v_doppler_1.value,
                       'fitvcen' : fit_spec_mod.lambda_0_1.value,
                       'fitEW' : fit_y.equivalent_width(x_0=fit_spec_mod.lambda_0_1)[0]}

    return line_properties

def write_out(hdulist,filename='spectrum.fits'):
	hdulist.writeto(filename, overwrite=True)
	return ""
