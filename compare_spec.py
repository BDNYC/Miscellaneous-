''' 
This procedure compares two or more individual spectra. It plots their NIR spectra normalized by band (J, H, and K).
'''

def get_spec(unum, separate=False, bandnames=None, bandlimits=None, bandnorms=None, retCoord=False, retST=False):
    '''
    This function finds and pulls spectrum from the BDNYC database. Specifically, it looks for low-res, NIR SpeX Prism spectra. It can split and normalize it by the NIR bands (J, H, and K).
    
    *unum*
      String with the U-number of the target (e.g. U20268).
    *separate*
      Boolean, whether to split and normalize the spectrum by NIR bands (J, H, and K)
    *bandnames*
      List with strings containing the band names used to separate the spectrum (when separate=True)
    *bandlimits*
      Dictionary with keys *bandnames*, each key containing a list with float numbers specifying the bands limits (e.g. {'J': [0.8,1.4], 'H': [1.4,1.9], 'K': [1.9,2.4]}
    *bandnorms*
      Same structure as *bandlimits*, this time the float numbers specify the wavelength ranges used to normalize the bands (when separate=True)
    *retCoord*
      Boolean, whether to pull from the database the target coordinates
    *retST*
      Boolean, whether to pull from the database the target spectral type
    '''
    
    import pickle
    import BDNYC
    import astrotools as at
    
    # 1. Initialize variables ---------------------------------------
    FOLDER_DB = '/Users/alejo/Dropbox/Python/Python_Database/'
    FILE_DB = 'BDNYCData.txt'
    
    # 2. Load database ----------------------------------------------
    f = open(FOLDER_DB + FILE_DB,'rb')
    bdnyc = pickle.load(f)
    f.close()
    
    # 3. Check data available for target ----------------------------
    availData = bdnyc.show_data(unum, dump=True)
    if availData is None:
        return
    
    # 4. Find Spex Prism data ---------------------------------------
    spexFound = False
    for row in availData:
        # Check that row is data row
        try:
            row[0] + 1
        except TypeError:
            continue
        # Check that row is nir row
        try:
            loc = row.index('nir')
        except ValueError:
            continue
        # Check that row is low res row
        try:
            loc = row.index('low')
        except ValueError:
            continue
        # Check that row is Spex Prism row
        instr = row[3].lower()
        loc = instr.find('spex')
        if loc != -1:
            spexFound = True
            specIdx = row[0]
            break
    if not spexFound:
        print 'Spectrum for target not found.'
        return
    
    # 5. Fetch target parameters if requested -----------------------
    params = []
    if retCoord:
        ra = availData[4][1][0:5]
        dec = availData[5][1][0:6]
        coord = ra + dec
        coord = coord.replace(' ','')
        params.append(coord)
    if retST:
        st = availData[3][1]
        params.append(st)
    
    # 6. Get spectrum -----------------------------------------------
    specRaw = bdnyc.get_data(unum, specIdx)
    
    # 7. Separate spectrum by bands ---------------------------------
    if separate:
        spec = [None] * 3
        for bIdx, band in enumerate(bandnames):
            bLim = bandlimits[band][0]
            bMax = bandlimits[band][1]
            
            idx1 = np.where(specRaw[0,:] >= bLim)
            idx2 = np.where(specRaw[0,:] <= bMax)
            idx = np.intersect1d(idx1[0], idx2[0])
            if len(idx) == 0:
                print 'Error in spectrum range.'
                return
            spec[bIdx] = specRaw[:,idx]
        
        # 8. Normalize spectrum -------------------------------------
        specNorm = [None] * 3
        for bIdx, band in enumerate(bandnames):
            specNorm[bIdx] = at.norm_spec(spec[bIdx], bandnorms[band])[0]
        
        if params != []:
            return specNorm, params
        else:
            return specNorm
    else:
        if params != []:
                return specRaw, params
        else:
            return specRaw


def plotspec(specData, bandNames, limits, objID, plotInput=None, templ=True):
    
    import numpy
    import matplotlib.pyplot as plt
    import types
    import pdb
    
    # 1) Check data consistency ===============================================
    try:
        specData.keys()
    except AttributeError:
        print 'Spectra not received as dictionaries.'
        return
    try:
        limits.keys()
    except AttributeError:
        print 'Limits not received as dictionaries.'
        return
    
    # 2) Initialize variables and color sets to use in plots ==================
    COLOR_SET = numpy.array(['#CC3333','#FF0000','#CC0000','#990000','#CC3300', \
                             '#FF3333','#FF6666','#FF3399','#CC0099','#FF0066', \
                             '#663300','#CC9900','#FFCC33','#666600','#669966', \
                             '#666666','#99CC99','#66CC99','#CCFF00','#66FF33', \
                             '#009933','#006600','#003300','#000066','#3333FF', \
                             '#33CCFF','#00FFFF','#9999FF','#3399CC','#0000CC'])
                # 0-plum, 1-red, 2-indian red, 3-maroon, 4-brick,
                # 5-tomato, 6-salmon, 7-fuchsia, 8-deep pink, 9-pink,
                # 10-brown, 11-chocolate, 12-wheat, 13-dk olive, 14-olive,
                # 15-silver, 16-lt green, 17-aquamarine, 18-yellow green, 19-lime,
                # 20-green, 21-forest, 22-dk green, 23-navy, 24-blue
                # 25-sky blue, 26-lt blue, 27-orchid, 28-steel blue, 29-royal blue
    colors = [None] * 16
    colors[15] = COLOR_SET[[1,3,6,7,10,11,12,18,19,20,21,24,25,27,29]].tolist()
    colors[14] = COLOR_SET[[1,3,6,7,11,12,18,19,20,21,24,25,27,29]].tolist()
    colors[13] = COLOR_SET[[1,3,6,7,11,12,19,20,21,24,25,27,29]].tolist()
    colors[12] = COLOR_SET[[1,3,6,7,11,12,19,20,21,25,27,29]].tolist() 
    colors[11] = COLOR_SET[[1,3,6,11,12,19,20,21,25,27,29]].tolist()
    colors[10] = COLOR_SET[[1,6,11,12,19,20,21,25,27,29]].tolist()
    colors[9]  = COLOR_SET[[1,6,11,12,19,20,25,27,29]].tolist()
    colors[8]  = COLOR_SET[[1,6,11,12,19,20,25,29]].tolist()
    colors[7]  = COLOR_SET[[1,6,11,19,20,25,29]].tolist()
    colors[6]  = COLOR_SET[[1,6,11,20,25,29]].tolist()
    colors[5]  = COLOR_SET[[1,6,11,20,29]].tolist()
    colors[4]  = COLOR_SET[[1,11,20,29]].tolist()
    colors[3]  = COLOR_SET[[1,20,29]].tolist()
    colors[2]  = COLOR_SET[[1,29]].tolist()
    
    numColors = len(specData['J'])
    plotColors = colors[numColors][:]
    plotColors.reverse()    
    BLACK = '#000000'
    GRAY  = '#CCCCCC'
    WHITE = '#FFFFFF'
    RED   = '#FF0000'
    X_LABEL = 'Wavelength ($\mu$m)'
    Y_LABEL = 'Normalized Flux (F$_{\lambda}$)'
    
    # 3) Initialize Figure ====================================================
    plt.close()
    plt.rc('font', size=9)
    fig = plt.figure(1, figsize=(9, 4.5))
    plt.clf()
    
    # 4) Generate Subplots ====================================================
    bandNames.reverse()
    for bandIdx, band in enumerate(bandNames):
        
        # 4.1) If band data is only one set, convert it into array of sets ----
        if specData[band][0] is not None:
            if len(specData[band][0]) > 3:
                specData[band] = [specData[band],]
        
        # 4.2) Initialize variables -------------------------------------------
        spLines = []
        minPlot = 1
        maxPlot = 1
        copyColors = list(plotColors)
        if band == 'J':
            textColors = []
        
        # 4.3) Initialize Subplot ---------------------------------------------
        tmpLeft = 0.06 + (2 - bandIdx) * 0.32
        subPlot = plt.figure(1).add_subplot(1,3,3 - bandIdx, \
                            position=[tmpLeft,0.085,0.265,0.85])
                                   # [left,bottom,width,height]
        subPlot.set_autoscale_on(False)
        
        # Create dummy axes instance to be able to later manipulate upper axis
        ax2 = subPlot.axes.twiny()
        
        # Set figure and axes labels (on left-most subplot)
        if bandIdx == 2:
            subPlot.set_xlabel(X_LABEL, position=(1.65,0.08), fontsize=10)
            subPlot.set_ylabel(Y_LABEL, fontsize=10)
        
        # 4.4) Plot spectra ---------------------------------------------------
        for specIdx, spec in enumerate(specData[band]):
            if spec is None:
                continue
            
            # Define plot parameters
            lnStyle = '-'
            objLabel = objID[specIdx]
            lnWidth = 0.8
            plotColor = copyColors.pop()
            
            if band == 'J':
                textColors.append(plotColor)
            wls = np.array(spec[0])
            fluxes = np.array(spec[1])
            
            subPlot.plot(wls, fluxes, color=plotColor, linestyle=lnStyle, \
                         dash_joinstyle='round', linewidth=lnWidth, label=objLabel)
            
            # Plot a dummy line on secondary axis to later modify upper x-axis
            if specIdx == 0:
                ax2.plot(wls,[-0.5] * len(wls), color=WHITE)
                
            # Track the highest & lowest y-axis values to fix y-axis limits later            
            tmpMin = numpy.nanmin(fluxes)
            if tmpMin < minPlot:
                minPlot = tmpMin
            tmpMax = numpy.nanmax(fluxes)
            if tmpMax > maxPlot:
                maxPlot = tmpMax
        
        # 4.5) Fix axes limits ------------------------------------------------
        minPlot = minPlot - minPlot * 0.1
        maxOff = 0.02
        maxPlot = maxPlot + maxPlot * maxOff
        
        plt.ylim(ymin=minPlot, ymax=maxPlot)
        subPlot.set_xlim(xmin=limits[band][0], \
                         xmax=limits[band][1] * 1.001)
        ax2.set_xlim(xmin=limits[band][0], \
                         xmax=limits[band][1] * 1.001)
        
        # 4.6) Customize y axis -----------------------------------------------
        subPlot.spines['left'].set_color('none')
        subPlot.spines['right'].set_color('none')
        subPlot.yaxis.set_ticks([])
        
        # 5) Add legend =======================================================
        if bandIdx == 2:
            objLegends = subPlot.legend(handlelength=0, handletextpad=0.1, \
                                      loc='lower right', numpoints=1, \
                                      labelspacing=0.2) #, \
                                      #bbox_to_anchor=(-0.05,0.98)
            objLegends.draw_frame(False)
            
            for legendIdx, legendText in enumerate(objLegends.get_texts()):                
                plt.setp(legendText, color=textColors[legendIdx], \
                           fontsize=7, fontname='Andale Mono')
        
    return fig


# ============================= PROCEDURE =====================================

# 1. LOAD RELEVANT MODULES ----------------------------------------------------
import astrotools as at
import numpy as np
import matplotlib.pyplot as plt
import os

# 2. SET UP VARIABLES ---------------------------------------------------------
# Enter U-numbers with the "U"
UNUMS = ['U50125','U50188']

FOLDER_NIR = '/Users/alejo/KCData/NIR/'
FOLDER_OUT = '/Users/alejo/KCData/Output/compare/'

# Variables used to plot result
BANDS = ['J','H','K']
BAND_LIMS = {}.fromkeys(BANDS)
BAND_NORMS = {}.fromkeys(BANDS)
for band in BANDS:
    BAND_LIMS[band] = [None] * 2
    BAND_NORMS[band] = [None] * 2
BAND_LIMS['J'][0] = 0.8
BAND_LIMS['J'][1] = 1.4 
BAND_LIMS['H'][0] = 1.4
BAND_LIMS['H'][1] = 1.9
BAND_LIMS['K'][0] = 1.9
BAND_LIMS['K'][1] = 2.4
BAND_NORMS['J'][0] = 0.87
BAND_NORMS['J'][1] = 1.39
BAND_NORMS['H'][0] = 1.41
BAND_NORMS['H'][1] = 1.89
BAND_NORMS['K'][0] = 1.91
BAND_NORMS['K'][1] = 2.39

# 3. GET SPECTRA --------------------------------------------------------------
# Initialize holders of consolidated spectra and parameters
all_params = []
spectra = {}.fromkeys(BANDS)
for band in BANDS:
    spectra[band] = []

# Loop through provided U-numbers and consolidate them
for uIdx, unum in enumerate(UNUMS):
    specN, parameters = get_spec(unum, separate=True, bandnames=BANDS, \
                             bandlimits=BAND_LIMS, bandnorms=BAND_NORMS, \
                             retCoord=True, retST=True)
    if specN is None:
        continue
    else:
        all_params.append([unum, parameters[0], parameters[1]])
        for bdIdx, band in enumerate(BANDS):
            spectra[band].append(specN[bdIdx])
    
    # if len(unum) > 3:
    #     specRaw, lblRaw = noc.main(unum, plot=False, lbl=True)
    #     lblRaw[0] = lblRaw[0] + ' peculiar'
    # else:
    #     templRaw = noc.main(unum, 'f', plot=False, templ=True)
    #     lblRaw = [unum + ' field template']
    #     specRaw = {}.fromkeys(BANDS)
    #     specRaw['J'] = [templRaw[0]]
    #     specRaw['H'] = [templRaw[1]]
    #     specRaw['K'] = [templRaw[2]]
    # 
    # J_spectra.append(specRaw['J'][0])
    # H_spectra.append(specRaw['H'][0])
    # K_spectra.append(specRaw['K'][0])
    # labels.append(lblRaw[0])

# 4. CREATE LABELS USING PARAMETERS -------------------------------------------
labels = []
for pIdx, param in enumerate(all_params):
    tmplbl = param[0] + ' ' + param[1] + ' ' + param[2]
    labels.append(tmplbl)

# 5. PLOT ALL SPECTRA ---------------------------------------------------------
figure = plotspec(specData=spectra, bandNames=BANDS, limits=BAND_LIMS, \
                  objID=labels)

plt.savefig(FOLDER_OUT + '/comp_spec.pdf', dpi=600)