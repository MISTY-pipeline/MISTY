import yt
import trident
import numpy as np
import argparse

from astropy.io import fits

import sys
import time
import MISTY

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="makes a MISTY spectrum")

    parser.add_argument('--ray', dest='ray', type=str, action='store',
                        default="",
                        help='a yt ray object')

    parser.add_argument('--sim', dest='sim', type=str, action='store',
                        default='WindTest/DD0010/DD0010',
                        help='name of a simulation to be loaded by yt')

    parser.add_argument('--author', dest='author', type=str, action='store',
                        default='A MISTY User',
                        help='name of person running this script')

    parser.add_argument('--param_file',dest='paramfile',type=str,action='store',                        default='a parameter text file')
    # parser.error("blah")

    args = parser.parse_args()
    return args

def create_misty_spectrum():

    my_line_list = ['H I 1216', 'C II 1335', 'C IV 1548', 'O VI 1032']

    if len(args.ray) == 0:
        ## read in example data
        ds = yt.load(args.sim)

        print "ray start and end are hardcoded and should be passed in !!!!!!"
        ## create ray start and end
        ray_start = [1,0,1]
        ray_end = [1,2,1]
        ray = trident.make_simple_ray(ds,start_position=ray_start,end_position=ray_end,
                        lines=my_line_list,ftype='gas')
        filespecout_base = 'spectrum_'+ds.basename
    else:
        ray = yt.load(args.ray)
        ray_start = ray.light_ray_solution[0]['start']
        ray_end = ray.light_ray_solution[0]['end']
        filespecout_base = 'spectrum_'+ray.light_ray_solution[0]['filename']
        print ray_start, ray_end, filespecout_base

    hdulist = MISTY.write_header(ray,start_pos=ray_start,end_pos=ray_end,
                      lines=my_line_list, author=args.author)

    ## put parameter file into the fits file
    MISTY.write_parameter_file(args.paramfile,hdulist=hdulist)

    for line in my_line_list:
        sg = MISTY.generate_line(ray,line,write=True,hdulist=hdulist)
        filespecout = filespecout_base+'_'+line.replace(" ", "_")+'.png'
        sg.plot_spectrum(filespecout,flux_limits=(0.0,1.0))

    MISTY.write_out(hdulist,filename='lauren_advance_spectrum.fits')


if __name__ == "__main__":

    args = parse_args()
    create_misty_spectrum()
    sys.exit("~~~*~*~*~*~*~all done!!!! spectra are fun!")
