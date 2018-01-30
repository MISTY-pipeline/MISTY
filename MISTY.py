from __future__ import print_function

import datetime
import getpass
import os.path
import time

import numpy as np
from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter

import trident
from spectacle.analysis.line_finder import LineFinder
from spectacle.core.spectrum import Spectrum1D

# ldb = trident.LineDatabase('lines.txt')
ldb = trident.LineDatabase('atom_wave_gamma_f.dat')


def write_header(ray, start_pos=None, end_pos=None, lines=None, **kwargs):
    # begin making fits header
    prihdr = fits.Header()
    prihdr['AUTHOR'] = kwargs.get("author", getpass.getuser())
    prihdr['DATE'] = datetime.datetime.now().isoformat()
    prihdr['RAYSTART'] = str(start_pos[0]) + "," + \
        str(start_pos[1]) + "," + str(start_pos[2])
    prihdr['RAYEND'] = str(end_pos[0]) + "," + \
        str(end_pos[1]) + "," + str(end_pos[2])
    prihdr['SIM_NAME'] = ray.basename
    prihdr['NLINES'] = str(len(np.array(lines)))
    prihdr['DOI'] = "doi.corlies2017.paper.thisistotesnotmadeup"
    prihdr['PAPER'] = "Corlies et al. (2017) ApJ, ###, ###"
    prihdr['EUVB'] = "HM12"  # probably shouldn't be hardcoded
    prihdr['IMPACT'] = (kwargs.get("impact", "undef"), "impact parameter, kpc")
    prihdr['ANGLE'] = (kwargs.get("angle", "undef"), "radians")

    lines = ldb.parse_subset(lines)

    i = 1
    for line in lines:
        keyword = 'LINE_' + str(i)
        prihdr[keyword] = line.name
        i += 1
    prihdu = fits.PrimaryHDU(header=prihdr)
    sghdulist = fits.HDUList([prihdu])
    return sghdulist


def write_parameter_file(ds, filename=None, hdulist=None):
    if type(hdulist) != fits.hdu.hdulist.HDUList:
        raise ValueError(
            'Must pass HDUList in order to write. Call write_header first.')

    # is a filename given? then use that
    if filename != None and os.path.isfile(filename):
        param_file = np.genfromtxt(filename, delimiter='=', dtype=str,
                                   autostrip=True)
        col1 = fits.Column(name='PARAMETERS', format='A50',
                           array=param_file[:, 0])
        col2 = fits.Column(name='VALUES', format='A50', array=param_file[:, 1])
    else:
        #  use ds.parameters
        col1 = fits.Column(name='PARAMETERS', format='A50',
                           array=list(ds.parameters.keys()))
        col2 = fits.Column(name='VALUES', format='A50',
                           array=[str(x) for x in ds.parameters.values()])

    col_list = [col1, col2]
    cols = fits.ColDefs(col_list)
    sghdr = fits.Header()

    sghdr['SIM_CODE'] = ds.dataset_type
    print("---> SIM_CODE set to ", ds.dataset_type,
          "if you don't like this, change it!")
    sghdr['COMPUTER'] = 'pleiades'
    print("---> ASSUMING PLEIADES FOR NOW BUT SHOULD BE PASSSSSSED IN")

    # primary_hdu = fits.PrimaryHDU(header=sghdr)

    sghdu = fits.BinTableHDU.from_columns(cols, header=sghdr)
    hdulist.append(sghdu)

    return sghdu


def generate_line(ray, line, write=False, use_spectacle=True, hdulist=None):
    '''
    input: a lightray and a line; writes info to extension of hdulist
    '''
    if write and type(hdulist) != fits.hdu.hdulist.HDUList:
        raise ValueError(
            'Must pass HDUList in order to write. Call write_header first.')

    if not isinstance(line, trident.Line):
        # ldb = trident.LineDatabase('lines.txt')
        ldb = trident.LineDatabase('atom_wave_gamma_f.dat')
        line_out = ldb.parse_subset(line)
        print(line, line_out)
        line_out = line_out[0]

    ar = ray.all_data()
    lambda_rest = line_out.wavelength
    if line_out.name == "H I 1216":
        padding = 5.
    else:
        padding = 5.
    lambda_min = lambda_rest * (1 + min(ar['redshift_eff'])) - padding
    lambda_max = lambda_rest * (1 + max(ar['redshift_eff'])) + padding

    sg = trident.SpectrumGenerator(lambda_min=lambda_min.value,
                                   lambda_max=lambda_max.value,
                                   dlambda=0.0001,
                                #    line_database='lines.txt'
                                   line_database='atom_wave_gamma_f.dat'
                                   )
    sg.make_spectrum(ray, lines=line_out.name, min_tau=1.e-5,
                     store_observables=True)

    if write and str(line_out) in sg.line_observables_dict:
        redshift = fits.Column(name='redshift', format='E',
                               array=(sg.lambda_field.value / lambda_rest - 1))
        wavelength = fits.Column(name='wavelength', format='E',
                                 array=sg.lambda_field, unit='Angstrom')
        tau = fits.Column(name='tau', format='E', array=sg.tau_field)
        flux = fits.Column(name='flux', format='E', array=sg.flux_field)
        col_list = [redshift, wavelength, tau, flux]

        for key in sg.line_observables_dict[str(line_out)].keys():
            col = fits.Column(name='sim_' + key, format='E',
                              array=sg.line_observables_dict[str(line_out)][key])
            col_list = np.append(col_list, col)

        cols = fits.ColDefs(col_list)
        sghdr = fits.Header()
        sghdr['LINENAME'] = line_out.name
        print("----->>>>using ", line_out.name,
              "as LINENAME, whereas ", line, " was passed. Change?")
        sghdr['RESTWAVE'] = (line_out.wavelength, "Angstroms")
        sghdr['F_VALUE'] = line_out.f_value
        sghdr['GAMMA'] = line_out.gamma

        # want to leave blank spaces now for values that we're expecting to generate for MAST
        # first let's add some spaces for the simulated, tau-weighted values!
        sghdr['SIM_TAU_HDENS'] = -9999.
        sghdr['SIM_TAU_TEMP'] = -9999.
        sghdr['SIM_TAU_METAL'] = -9999.
        sghdr['TOT_COLUMN'] = (np.log10(np.sum(
            sg.line_observables_dict[line_out.identifier][
                'column_density'].value)), "log cm^-2")

        # we're also going to want data from spectacle
        if use_spectacle:
            lines_properties = get_line_info(sg)
            for key in lines_properties:
                sghdr[key] = lines_properties[key]

        sghdu = fits.BinTableHDU.from_columns(
            cols, header=sghdr, name=line_out.name)

        hdulist.append(sghdu)

    return sg


def get_line_info(sg):
    '''
    runs spectacle on a trident spectrum object and returns absorber properties
    '''
    import astropy.units as u
    from spectacle.analysis.line_finder import LineFinder

    disp = sg.lambda_field / (1 + 2.) * u.Unit('Angstrom')
    flux = sg.flux_field
    tau = sg.tau_field
    sg_line = sg.line_list[0]

    # This process will find lines in the trident spectrum
    # and assign the values set in the `defaults` dict to
    # the new lines found.
    # try:

    # Create a dictionary to hold the default values we want 
    # the lines to have
    default_values = dict(
        lambda_0=sg_line['wavelength'].value * u.Unit('Angstrom'),
        f_value=sg_line['f_value'],
        gamma=sg_line['gamma'],
        fixed={'lambda_0': True,
                'delta_v': True,
                'delta_lambda': False})

    # Have the line finder attempt to find absorption features. Fit the
    # result to the data.
    spec_mod = LineFinder(disp, flux,
                         ion_name=sg_line,
                         redshift=0,  # This could be tied to sg, e.g.
                         data_type='flux',
                         defaults=default_values,
                         threshold=0.1).fit()
    fitter = LevMarLSQFitter()
    fit_spec_mod = fitter(spec_mod, disp, flux, maxiter=2000)

    # This will be more user-friendly in the future: get all the Voigt
    # profiles that make up this spectrum model.
    line_mods = [x for x in fit_spec_mod if hasattr(x, 'lambda_0')]

    line_properties = {
        'NCOMP': len(line_mods)
    }

    for i, line in enumerate(line_mods):
        dv90 = line.dv90()
        fwhm = line.fwhm()

        line_properties.update({
            'fitcol' + str(i): (line.column_density.value, line.column_density.unit.to_string()),
            'fitb' + str(i): (line.v_doppler.value, line.v_doppler.unit.to_string()),
            'fitlcen' + str(i): (line.lambda_0.value + line.delta_lambda.value, line.lambda_0.unit.to_string()),
            # 'fitEW' + str(i): (line.equivalent_width.value, line.equivalent_width.unit.to_string())
            'fitdv90' + str(i): (dv90.value, dv90.unit.to_string()),
            'fitfwhm' + str(i): (fwhm.value, fwhm.unit.to_string())
        })

    # except Exception:
    #     print("***** --->> line finding SO did not work ****")
    #     lines_properties = {'NCOMP': 0}

    return line_properties


def write_out(hdulist, filename='spectrum.fits'):
    print("printing to .... " + filename)
    hdulist.writeto(filename, overwrite=True, output_verify='fix')
    return ""
