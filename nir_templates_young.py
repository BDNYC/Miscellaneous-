''' 
This procedure plots the NIR L template for the gamma objects of type COMPTYPE
against the NIR L templates calculated using nir_opt_comp_strip.

'''
def plotspec(specData, bandNames, limits, objID, plotInstructions, compName, plotInput=None, figNum=1):
# Plots set of spectral data and saves plots in a PDF file.
# specData and limits must be dictionaries.
    
    import numpy
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
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
    GRAYS = ['#585858', '#686868', '#707070', '#808080', '#909090', \
             '#A0A0A0', '#B0B0B0', '#C0C0C0', '#D0D0D0', '#E0E0E0']
    # Color order goes from reds to blues
    colors = ['#FF0000','#FF6699','#CC9900','#FFCC33', \
              '#66FF33','#009933','#99FFFF','#33CCFF','#0066FF']
    colors.reverse()
    BLACK = '#000000'
    GRAY  = '#CCCCCC'
    WHITE = '#FFFFFF'
    RED   = '#FF0000'
    X_LABEL = 'Wavelength ($\mu$m)'
    Y_LABEL = 'Normalized Flux (F$_{\lambda}$) + constant'
    templName = compName + '$\gamma$'
    
    # 3) Initialize Figure ====================================================
    plt.close()
    plt.rc('font', size=9)
    fig = plt.figure(figNum, figsize=(7.33,8.6))
    plt.clf()
    
    # 4) Generate Subplots ====================================================
    bandNames.reverse()
    for bandIdx, band in enumerate(bandNames):
        
        # 4.1) If band data is only one set, convert it into array of sets ----
        if specData[band][0] is not None:
            if len(specData[band][0]) > 6:
                specData[band] = [specData[band],]
        
        # 4.2) Initialize variables -------------------------------------------
        spLines = []
        minPlot = 1
        maxPlot = 1
        copyColors = list(colors)
        
        # 4.3) Initialize Subplot ---------------------------------------------
        tmpLeft = 0.06 + (2 - bandIdx) * 0.32
        subPlot = plt.figure(figNum).add_subplot(1,3,3 - bandIdx, \
                            position=[tmpLeft,0.05,0.265,0.92])
                                   # [left,bottom,width,height]
        subPlot.set_autoscale_on(False)
        
        # Create dummy axes instance to be able to later manipulate upper axis
        ax2 = subPlot.axes.twiny()
        
        # Set figure and axes labels
        if bandIdx == 2:
            subPlot.set_xlabel(X_LABEL, position=(1.65,0.08), fontsize=10)
            subPlot.set_ylabel(Y_LABEL, fontsize=10)
        
        # 4.4) Plot spectra ---------------------------------------------------
        offset = 0
        lastColor = BLACK
        for specIdx, spec in enumerate(specData[band]):
            if spec is None:
                continue
                
            plotType = plotInstructions[specIdx]
            objLabel = objID[specIdx]
            
            # Determine which field templates to highlight as similar to young ones
            similar = False # to indicate young template that resembles field one
            if plotType == 'template':
                if compName == 'L0' and specIdx:
                    if bandIdx == 0:
                        if objLabel == 'L0':
                            similar = True
                elif compName == 'L2':
                    if bandIdx == 0:
                        if objLabel == 'L0' or objLabel == 'L1' or objLabel == 'L2':
                            similar = True
                    elif bandIdx == 1:
                        if objLabel == 'L7':
                            similar = True
                elif compName == 'L3':
                    if bandIdx == 0:
                        if objLabel == 'L1' or objLabel == 'L2' or objLabel == 'L3':
                            similar = True
                    elif bandIdx == 2:
                        if objLabel == 'L3' or objLabel == 'L4':
                            similar = True
                elif compName == 'L4':
                    if bandIdx == 0:
                        if objLabel == 'L1' or objLabel == 'L2' or objLabel == 'L3':
                            similar = True
                    elif bandIdx == 1:
                        if objLabel == 'L7':
                            similar = True
                    elif bandIdx == 2:
                        if objLabel == 'L4' or objLabel == 'L5':
                            similar = True
            
            # Define plot parameters
            lnStyle = '-'
            if plotType == 'young':
                plotColor = BLACK
                lnWidth = 0.7
            elif plotType == 'template':
                plotColor = copyColors.pop()
                lastColor = plotColor
                if similar:
                    lnWidth = 1.6
                else:
                    lnWidth = 0.7
            if specIdx > 0 and specIdx % 2 == 0:
                if bandIdx == 0:
                    offset = offset + 0.55
                elif bandIdx == 1:
                    offset = offset + 0.75
                else:
                    offset = offset + 1
            
            wls = np.array(spec[0])
            fluxes = np.array(spec[1])
            
            # Plot spectral strip when available
            if plotType == 'template' and len(spec) > 3:
                errs = np.array(spec[2])
                mins = np.array(spec[3])
                maxs = np.array(spec[4])
                
                for wlIdx, wl in enumerate(wls):
                    # Skip first and last points
                    if wlIdx == 0:
                        continue
                    elif wlIdx == len(wls) - 1:
                        continue
                    elif not numpy.isfinite(mins[wlIdx]):
                        continue
                    
                    # Set location of lower left corner of rectangle
                    rect_x = wl - ((wl - wls[wlIdx - 1]) / 2)
                    rect_y = mins[wlIdx] + offset
                    
                    # Set dimensions of rectangle
                    rect_width = ((wl - wls[wlIdx - 1]) / 2) + \
                                 ((wls[wlIdx + 1] - wl) / 2)
                    rect_height = maxs[wlIdx] - mins[wlIdx]
                    
                    # Set color fill of rectangle
                    err = errs[wlIdx]
                    if err > 0.19:
                        grayIdx = 9
                    elif err > 0.17:
                        grayIdx = 8
                    elif err > 0.16:
                        grayIdx = 7
                    elif err > 0.15:
                        grayIdx = 6
                    elif err > 0.13:
                        grayIdx = 5
                    elif err > 0.11:
                        grayIdx = 4
                    elif err > 0.09:
                        grayIdx = 3
                    elif err > 0.07:
                        grayIdx = 2
                    elif err > 0.06:
                        grayIdx = 1
                    else:
                        grayIdx = 0
                    rect_color = GRAYS[grayIdx]
                    
                    # Add rectangle to plot
                    rect_patch = mpatches.Rectangle(xy=(rect_x, rect_y), \
                                                    width=rect_width, edgecolor='none', \
                                                    height=rect_height, color=rect_color)
                    subPlot.add_patch(rect_patch)
            
            # Plot spectral lines
            subPlot.plot(wls, fluxes + offset, color=plotColor, linestyle=lnStyle, \
                    dash_joinstyle='round', linewidth=lnWidth, label=objLabel, \
                    drawstyle='steps-mid')
            
            # Plot a dummy line on secondary axis to later modify upper x-axis
            if specIdx == 0:
                ax2.plot(wls, [-0.5] * len(wls), color=WHITE)
                
            # Track the highest & lowest y-axis values to fix y-axis limits later            
            tmpMin = numpy.nanmin(fluxes)
            if tmpMin < minPlot:
                minPlot = tmpMin
            tmpMax = numpy.nanmax(fluxes + offset)
            if tmpMax > maxPlot:
                maxPlot = tmpMax
            
            # Track the highest y-axis value for the tail of each band
            if specIdx % 2 == 0:
                tailMax = numpy.nanmax(fluxes[-8:-1])
            else:
                currMax = numpy.nanmax(fluxes[-8:-1])
                if currMax > tailMax:
                    tailMax = currMax
            
            # Add annotation to template plot
            if specIdx % 2 == 1:
                if bandIdx == 0:
                    yTextLoc = 13
                else:
                    yTextLoc = 8
                textLoc = (0, yTextLoc)
                annotLoc = (wls[-10], tailMax + offset)
                subPlot.annotate(objLabel, xy=annotLoc, xycoords='data', \
                                color=lastColor, xytext=textLoc, \
                                textcoords='offset points')
                textLoc = (0, yTextLoc + 7)
                subPlot.annotate(templName, xy=annotLoc, xycoords='data', \
                                color=BLACK, xytext=textLoc, \
                                textcoords='offset points')
                                
        # 4.5) Fix axes limits ------------------------------------------------
        minPlot = minPlot - minPlot * 0.1
        maxOff = 0.02
        maxPlot = maxPlot + maxPlot * maxOff
        
        plt.ylim(ymin=minPlot, ymax=maxPlot)
        subPlot.set_xlim(xmin=limits[band]['lim'][0], \
                         xmax=limits[band]['lim'][1] * 1.001)
        ax2.set_xlim(xmin=limits[band]['lim'][0], \
                         xmax=limits[band]['lim'][1] * 1.001)
        
        # 4.6) Customize y axis -----------------------------------------------
        subPlot.spines['left'].set_color('none')
        subPlot.spines['right'].set_color('none')
        subPlot.yaxis.set_ticks([])
    
    return fig

# ============================= PROCEDURE =====================================

# 1. LOAD RELEVANT MODULES ----------------------------------------------------
import asciidata as ad
import astrotools as at
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pdb


# 2. SET UP VARIABLES ---------------------------------------------------------
COMPTYPE = raw_input('Enter young spectral type (e.g. L0): ').upper()
grav = raw_input('Enter gravity (b or g): ').upper()
DELL_CHAR = '\t' # Delimiter character
FOLDER_OUT = '/Users/alejo/KCData/Output/special'
FOLDER_TEMPL = '/Users/alejo/KCData/Output/templates/'
SP_TYPES = ['L0','L1','L2','L3','L4','L5','L6','L7','L8']

BANDS = ['J','H','K']
BAND_LIMS = {}.fromkeys(BANDS)
for bandKey in BANDS:
    BAND_LIMS[bandKey] = dict(lim = [None] * 2)
BAND_LIMS['J']['lim'][0] = 0.8
BAND_LIMS['J']['lim'][1] = 1.4 
BAND_LIMS['H']['lim'][0] = 1.4
BAND_LIMS['H']['lim'][1] = 1.9
BAND_LIMS['K']['lim'][0] = 1.9
BAND_LIMS['K']['lim'][1] = 2.4

# 3. GET YOUNG SPECTRAL NIR TEMPLATE ------------------------------------------
# Fetch template from ascii files generated by make_templ.py
youngTempl = {}.fromkeys(BANDS)
for band in BANDS:
    youngTempl[band] = []
for bdIdx, band in enumerate(BANDS):
        fileNm = COMPTYPE + band + '_' + grav + '.txt'
        try:
            templRaw = ad.open(FOLDER_TEMPL + fileNm, delimiter=DELL_CHAR)
        except IOError:
            templRaw = None
        
        if templRaw is None:
            print 'No template for type requested.'
            sys.exit(0)
        else:
            templLs = np.array(templRaw).tolist()
            templRaw = ''
            youngTempl[band].append(templLs)

# 4. GET SPECTRAL NIR FIELD TEMPLATES -----------------------------------------
spTypes = []
spectra = {}.fromkeys(BANDS)
for band in BANDS:
    spectra[band] = []

for idxTp, spTp in enumerate(SP_TYPES):
    # Attach copy of young template
    for bdIdx, band in enumerate(BANDS):
        spectra[band].append(youngTempl[band][0])
    
    # Fetch template from ascii files generated by make_templ.py
    for bdIdx, band in enumerate(BANDS):
        fileNm = spTp + band + '_f.txt'
        try:
            templRaw = ad.open(FOLDER_TEMPL + fileNm, delimiter=DELL_CHAR)
        except IOError:
            print fileNm + ' file could not be opened. Check problem and try again.'
            sys.exit(0)
        
        templLs = np.array(templRaw).tolist()
        templRaw = ''
        spectra[band].append(templLs)
    
    # Append plot label list
    spTypes.append(SP_TYPES[idxTp])
    spTypes.append(SP_TYPES[idxTp])

spTypes.reverse()
for band in BANDS:
    spectra[band].reverse()

# 5. PLOT ALL SPECTRA ---------------------------------------------------------
plotInstructions = ['template','young'] * len(spTypes)

figure = plotspec(spectra, BANDS, BAND_LIMS, spTypes, \
                  plotInstructions, COMPTYPE)

plt.savefig(FOLDER_OUT + '/' + COMPTYPE + grav + '-field.pdf', dpi=600)