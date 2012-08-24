''' This generates separate ascii files for all templates, by band (J, H, and K). The files contain five columns: wavelength, average flux, average flux variance, min flux, max flux.'''

import nir_opt_comp_strip as nocs
import astrotools as at

FOLDER_OUT = '/Users/alejo/Dropbox/KCData/Output/templates/'
TYPES = ['L0','L1','L2','L3','L4','L5','L6','L7','L8']
GRAVS = ['f','g','b']
BANDS = ['J','H','K']

for sptp in TYPES:
    print sptp
    for grav in GRAVS:
        templ = nocs.main(sptp, grav, templ=True, plot=False)
        if templ is None:
            continue
        
        print ' ' + grav
        for bdidx, band in enumerate(templ):
            # Create template spectrum file
            # columns are: wavelength, mean flux, mean flux variance, min flux, max flux
            at.create_ascii(band, FOLDER_OUT + sptp + BANDS[bdidx] + '_' + grav)

