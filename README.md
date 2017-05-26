# The MAST Interface to Synthetic Telescopes with yt (MISTY)
## for observing simulations of the circumgalactic medium in absorption
### Molly Peeples, Lauren Corlies, Nicholas Earl
### May 26, 2017

### Acknowledgments:
* Supported by [HST AR #13919](http://adsabs.harvard.edu/abs/2014hst..prop13919P)
* Makes heavy use of [yt](yt-project.org), [Trident](http://trident-project.org/), and [astropy](http://www.astropy.org/).

#### File name convention:
* hlsp_misty_sim_halo_output_geometry_version_type.fits
* hlsp_misty_sim_halo_output_impactparameter-angle_version_type.fits
* hlsp_misty_foggie_halo008508_rd0042_i177.9-a1.30_v1_los.fits
These follow, as closely as possible, the [HSLP file naming conventions from MAST](https://archive.stsci.edu/hlsp/hlsp_guidelines.html#Filenames)


#### FITS files contain the following HDUs:
0. Primary Header, documented below. Includes any simulation or sightline based parameters to be searched over.
1. Binary table of the parameter file used to run the simulation
2. Extension name is ionic transitions as listed in the primary header, following element+ion+wavelength convention (e.g., "H I 1216"). Header contains absorber information, arrays contain spectra, documented below.
3. (and more) any other absorption lines included

So, e.g., a given file might look like:
Filename: spectra/hlsp_misty_foggie_halo008508_rd0042_i101.9-a1.94_v1_los.fits

No.    Name         Type      Cards   Dimensions   Format
0  PRIMARY     PrimaryHDU      22   ()      
1  PARAMS      BinTableHDU     15   608R x 2C   [50A, 50A]   
2  H I 1216    BinTableHDU     57   108837R x 9C   [E, E, E, E, E, E, E, E, E]   
3  SI II 1260  BinTableHDU     57   109162R x 9C   [E, E, E, E, E, E, E, E, E]   
4  C II 1335   BinTableHDU     57   109701R x 9C   [E, E, E, E, E, E, E, E, E]   
5  C III 977   BinTableHDU     57   107102R x 9C   [E, E, E, E, E, E, E, E, E]   
6  SI III 1207  BinTableHDU     57   108771R x 9C   [E, E, E, E, E, E, E, E, E]   
7  C IV 1548   BinTableHDU     57   111254R x 9C   [E, E, E, E, E, E, E, E, E]   
8  O VI 1032   BinTableHDU     57   107501R x 9C   [E, E, E, E, E, E, E, E, E]   

A thing to be done is to combine extensions 0 and 1 in a clean way, so the following extensions correspond to the LINE_no in the primary header.

In HDU 0, the header contains the following useful cards, for example (and these need to be cleaned up...):

AUTHOR  = 'molly   '                                                            
DATE    = '2017-05-25T18:36:03.110700'                                          
RAYSTART= '[ 0.48954376  0.47206     0.5078625 ] unitary'                       
RAYEND  = '[ 0.48954376  0.47206     0.5113375 ] unitary'                       
SIM_NAME= 'RD0042_ray_z_imp101.9_ang1.94_tri.h5'                                
NLINES  = '7       '                                                            
DOI     = 'doi.corlies2017.paper.thisistotesnotmadeup'                          
PAPER   = 'Corlies et al. (2017) ApJ, ###, ###'                                 
EUVB    = 'HM12    '                                                            
IMPACT  =    101.8704386189004 / kpc                                            
ANGLE   =    1.941119505343093 / radians                                        
LINE_1  = 'H I 1216'                                                            
LINE_2  = 'Si II 1260'                                                          
LINE_3  = 'C II 1335'                                                           
LINE_4  = 'C III 977'                                                           
LINE_5  = 'Si III 1207'                                                         
LINE_6  = 'C IV 1548'                                                           
LINE_7  = 'O VI 1032'                                                                   

HDU 2-N are then the individual lines, with header formats like:

LINENAME= 'H I 1216'                                                            
RESTWAVE=              1215.67                                                  
F_VALUE =                0.416                                                  
GAMMA   =          469000000.0                                                  
HIERARCH SIM_TAU_HDENS = -9999.0                                                
HIERARCH SIM_TAU_TEMP = -9999.0                                                 
HIERARCH SIM_TAU_METAL = -9999.0                                                
HIERARCH TOT_COLUMN = 14.29506725344586                                         
HIERARCH NCOMPONENTS = 5                                                        
FITB0   =    3355253.270032837                                                  
FITCOL0 =    14.13695228140041                                                  
FITEW0  = 'NaN     '                                                            
FITVCEN0=    1215.673395410955                                                  
FITEW1  =              -9999.0                                                  
FITCOL1 =              -9999.0                                                  
FITVCEN1=              -9999.0                                                  
FITB1   =              -9999.0                                                  
FITEW2  =              -9999.0                                                  
FITCOL2 =              -9999.0                                                  
FITVCEN2=              -9999.0                                                  
FITB2   =              -9999.0                                                  
FITEW3  =              -9999.0                                                  
FITCOL3 =              -9999.0                                                  
FITVCEN3=              -9999.0                                                  
FITB3   =              -9999.0                                                  
FITEW4  =              -9999.0                                                  
FITCOL4 =              -9999.0                                                  
FITVCEN4=              -9999.0                                                  
FITB4   =              -9999.0                                                  
EXTNAME = 'H I 1216'           / extension name                                 
TTYPE1  = 'wavelength'                                                          
TFORM1  = 'E       '                                                            
TUNIT1  = 'Angstrom'                                                            
TTYPE2  = 'tau     '                                                            
TFORM2  = 'E       '                                                            
TTYPE3  = 'flux    '                                                            
TFORM3  = 'E       '                                                            
TTYPE4  = 'sim_delta_lambda'                                                    
TFORM4  = 'E       '                                                            
TTYPE5  = 'sim_lambda_obs'                                                      
TFORM5  = 'E       '                                                            
TTYPE6  = 'sim_thermal_b'                                                       
TFORM6  = 'E       '                                                            
TTYPE7  = 'sim_column_density'                                                  
TFORM7  = 'E       '                                                            
TTYPE8  = 'sim_thermal_width'                                                   
TFORM8  = 'E       '                      

The "fit" parameters correspond to the column densities, equivalent widths, Doppler parameters, etc., as found by spectacle; the different numbers are supposed to be for different detected components. This is still being worked on. Note that the TTYPEs at the end are the arrays included in the extension; the wavelength, tau, and flux arrays are the synthetic data of most interest. The flux is normalized to have continuum = 1.
