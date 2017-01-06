import yt
import trident
import numpy as np

from astropy.io import fits

import time
import MISTY

## read in example data
fn = 'WindTest/DD0010/DD0010'
ds = yt.load(fn)

## create ray start and end
ray_start = [1,0,1]
ray_end = [1,2,1]

my_line_list = ['H I 1216', 'C II 1335', 'C IV 1548', 'O VI 1032']

ray = trident.make_simple_ray(ds,start_position=ray_start,end_position=ray_end,
                               data_filename="ray.h5",lines=my_line_list,ftype='gas')


hdulist = MISTY.write_header(ray,start_pos=ray_start,end_pos=ray_end,
                      lines=my_line_list, author='Lauren')


for line in my_line_list:
    sg = MISTY.generate_line(ray,line,write=True,hdulist=hdulist)
    filespecout = 'spectrum_'+ds.basename+'_'+line.replace(" ", "_")+'.png'
    sg.plot_spectrum(filespecout,flux_limits=(0.0,1.0))

MISTY.write_out(hdulist,filename='advance_spectrum.fits')
