import yt
import trident


fn = 'WindTest/DD0010/DD0010'
ds = yt.load(fn)

ray_start = [1,0,1]
ray_end = [1,2,1]

line_list = ['H','C','N','O','Mg','S','Si','Ne']

ray = trident.make_simple_ray(ds,start_position=ray_start,end_position=ray_end,
                                data_filename="ray.h5",lines=line_list,ftype='gas')


sg = trident.SpectrumGenerator('COS-G130M')
sg.make_spectrum(ray,lines=line_list)

filespecout = 'spectrum_'+ds.basename+'.png'
sg.plot_spectrum(filespecout,flux_limits=(0.998,1.0))

sg.save_spectrum('spec_raw.txt') 



