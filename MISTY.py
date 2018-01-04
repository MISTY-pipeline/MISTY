import numpy as np
from astropy.io import fits
import time
import trident
import spectacle

from spectacle.core.spectra import Spectrum1D
from spectacle.modeling.models import Absorption1D
from spectacle.core.lines import Line
from spectacle.modeling.fitting import DynamicLevMarFitter

import getpass
import datetime
import os.path

## ldb = trident.LineDatabase('lines.txt')
ldb = trident.LineDatabase('atom_wave_gamma_f.dat')

def write_header(ray,start_pos=None,end_pos=None,lines=None,**kwargs):
	## begin making fits header
    prihdr = fits.Header()
    prihdr['AUTHOR'] = kwargs.get("author",getpass.getuser())
    prihdr['DATE'] = datetime.datetime.now().isoformat()
    prihdr['RAYSTART'] = str(start_pos[0]) + "," + str(start_pos[1]) + "," + str(start_pos[2])
    prihdr['RAYEND'] = str(end_pos[0]) + "," + str(end_pos[1]) + "," + str(end_pos[2])
    prihdr['SIM_NAME'] = ray.basename
    prihdr['NLINES'] = str(len(np.array(lines)))
    prihdr['DOI'] = "doi.corlies2017.paper.thisistotesnotmadeup"
    prihdr['PAPER'] = "Corlies et al. (2017) ApJ, ###, ###"
    prihdr['EUVB'] = "HM12" ## probably shouldn't be hardcoded
    prihdr['IMPACT'] = (kwargs.get("impact","undef"), "impact parameter, kpc")
    prihdr['ANGLE'] = (kwargs.get("angle","undef"), "radians")

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

def generate_line(ray,line,write=False,use_spectacle=True,hdulist=None):
    '''
    input: a lightray and a line; writes info to extension of hdulist
    '''
    if (write == True) & ((type(hdulist) != fits.hdu.hdulist.HDUList)):
       raise ValueError('Must pass HDUList in order to write. Call write_header first.')

    if not isinstance(line,trident.Line):
        ## ldb = trident.LineDatabase('lines.txt')
        ldb = trident.LineDatabase('atom_wave_gamma_f.dat')
        line_out = ldb.parse_subset(line)
        print line, line_out
        line_out = line_out[0]

    ar = ray.all_data()
    lambda_rest = line_out.wavelength
    if line_out.name == "H I 1216":
        padding = 5.
    else:
        padding = 5.
    lambda_min = lambda_rest * (1+min(ar['redshift_eff'])) - padding
    lambda_max = lambda_rest * (1+max(ar['redshift_eff'])) + padding

    sg = trident.SpectrumGenerator(lambda_min=lambda_min.value, \
        lambda_max=lambda_max.value, dlambda=0.0001, \
        line_database = 'atom_wave_gamma_f.dat')
    sg.make_spectrum(ray,lines=line_out.name, min_tau=1.e-5, store_observables=True)

    if write and str(line_out) in sg.line_observables_dict:
    	redshift = fits.Column(name='redshift', format='E', array=(sg.lambda_field.value/lambda_rest - 1))
    	wavelength = fits.Column(name='wavelength', format='E', array=sg.lambda_field,unit='Angstrom')
    	tau = fits.Column(name='tau', format='E', array=sg.tau_field)
    	flux = fits.Column(name='flux', format='E', array=sg.flux_field)
    	col_list = [redshift, wavelength, tau, flux]

    	for key in sg.line_observables_dict[str(line_out)].keys():
    	    col = fits.Column(name='sim_'+key,format='E',array=sg.line_observables_dict[str(line_out)][key])
    	    col_list = np.append(col_list,col)

    	cols = fits.ColDefs(col_list)
    	sghdr = fits.Header()
    	sghdr['LINENAME'] = line_out.name
        print "----->>>>using ", line_out.name, "as LINENAME, whereas ", line, " was passed. Change?"
    	sghdr['RESTWAVE'] = (line_out.wavelength, "Angstroms")
    	sghdr['F_VALUE'] = line_out.f_value
    	sghdr['GAMMA'] = line_out.gamma


        ## want to leave blank spaces now for values that we're expecting to generate for MAST
        ## first let's add some spaces for the simulated, tau-weighted values!
        sghdr['SIM_TAU_HDENS'] = -9999.
        sghdr['SIM_TAU_TEMP'] = -9999.
        sghdr['SIM_TAU_METAL'] = -9999.
        sghdr['TOT_COLUMN'] = (np.log10(np.sum(sg.line_observables_dict[line_out.identifier]['column_density'].value)), "log cm^-2")

        ## we're also going to want data from spectacle
        if use_spectacle:
            lines_properties = get_line_info(sg)
            for key in lines_properties:
                sghdr[key] = lines_properties[key]

    	sghdu = fits.BinTableHDU.from_columns(cols,header=sghdr,name=line_out.name)

    	hdulist.append(sghdu)

    return sg

def get_line_info(sg):
    '''
    runs spectacle on a trident spectrum object and returns absorber properties
    '''
    disp = sg.lambda_field / (1+ 2.)
    flux = sg.flux_field
    sg_line = sg.line_list[0]
    spectrum = Spectrum1D(np.array(list(flux)), dispersion=np.array(list(disp)))
    tot_ew = spectrum.equivalent_width()[0]
    print "FYI, tot_ew = ", tot_ew
    # This process will find lines in the trident spectrum
    # and assign the values set in the `defaults` dict to
    # the new lines found.
    try:
        lines = spectrum.find_lines(threshold=0.02/max(1-spectrum.data),
                            min_dist=10,
                            smooth=True,
                            interpolate=True,
                            defaults=dict(
                                lambda_0=sg_line['wavelength'].value,
                                f_value=sg_line['f_value'],
                                gamma=sg_line['gamma'],
                                fixed={
                                    'lambda_0': True,
                                    'delta_v': True,
                                    'delta_lambda': False}
                           ))

        # Create absorption Spectrum1D from line information
        spec_mod = Absorption1D(lines=lines)

        # Create a Spectrum1D object from the Absorption1D model
        gen_spec = spec_mod(spectrum.dispersion)

        # Create a fitter. The default fitting routine is a LevMar.
        fitter = DynamicLevMarFitter()
        fit_spec_mod = fitter(spec_mod, spectrum.dispersion, spectrum.data,
                      maxiter=5, initialize=False)

        # Get the results of the fit
        fit_spec = fit_spec_mod(spectrum.dispersion)

        # OK now we want line properties
        NCOMP = len(fit_spec.lines)
        lines_properties = {'NCOMP' : (NCOMP, "number of fitted components")}
        lines_properties['totEW'] = (tot_ew, "Angstrom")
        for i in np.arange(1,NCOMP+1):
            lines_properties['fitcol'+str(i)] = (fit_spec_mod[i].column_density[0], "log cm^-2")
            lines_properties['fitb'+str(i)] = (fit_spec_mod[i].v_doppler[0]/100000., "km s^-1")
            lines_properties['fitlcen'+str(i)] = (fit_spec_mod[i].lambda_0[0] + fit_spec_mod[i].delta_v[0], "Angstrom, center of component, observed wavelength")
            lines_properties['fitEW'+str(i)] = (fit_spec.equivalent_width(x_0=fit_spec_mod.lambda_0_1)[0], "Angstrom, total equivalent width")
    except Exception:
        print "***** --->> line finding SO did not work ****"
        lines_properties = {'NCOMP': 0}


    return lines_properties

def write_out(hdulist,filename='spectrum.fits'):
    print "printing to .... " + filename
    hdulist.writeto(filename, overwrite=True)
    return ""
