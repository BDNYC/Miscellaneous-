''' 
This procedure plots NIR spectra fanned out. Reads fits files directly. This procedures does not split the NIR spectra by bands, it only normalizes them and draws them in black. specData must be a list, even if it is just one object.
'''

def plotspec(specData, limits, objID, wide=False):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import types
    
    # 1) Initialize variables =================================================
    numSpecs = len(specData)
    BLACK = '#000000'
    GRAY  = '#CCCCCC'
    WHITE = '#FFFFFF'
    RED   = '#FF0000'
    X_LABEL = 'Wavelength ($\mu$m)'
    if numSpecs == 1:
        Y_LABEL = 'Normalized Flux (F$_{\lambda}$)'
    else:
        Y_LABEL = 'Normalized Flux (F$_{\lambda}$) + constant'
    
    # 2) Initialize Figure ====================================================
    plt.close()
    if numSpecs < 2:
        height = 5
    else:
        height = 8.6
    plt.rc('font', size=9)
    fig = plt.figure(1, figsize=(7.33, height))
    plt.clf()
    
    # 3) Generate Subplots ====================================================
    # Initialize variables -----------------------------------------------
    spLines = []
    minPlot = 1
    maxPlot = 1
    
    # Initialize Subplot -------------------------------------------------
    subPlot = plt.figure(1).add_subplot(1,1,1, \
                                        position=[0.08,0.08,0.85,0.87])
                                        # [left,bottom,width,height]
    subPlot.set_autoscale_on(False)
    
    # Create dummy axes instance to be able to later manipulate upper axis
    ax2 = subPlot.axes.twiny()
    
    # Set figure and axes labels
    subPlot.set_xlabel(X_LABEL, position=(0.5,0.08), fontsize=9)
    subPlot.set_ylabel(Y_LABEL, fontsize=9)
    
    # 4.4) Plot spectra --------------------------------------------------
    offset = 0
    tailMax = 0
    for specIdx, spec in enumerate(specData):
        if spec is None:
            continue
        
        # Define plot parameters
        lnStyle = '-'
        objLabel = objID[specIdx]
        plotColor = BLACK
        lnWidth = 0.8
        
        wls = np.array(spec[0])
        fluxes = np.array(spec[1])
        
        # Plot spectral lines
        subPlot.plot(wls, fluxes + offset, color=plotColor, linestyle=lnStyle, \
                dash_joinstyle='round', linewidth=lnWidth, label=objLabel, \
                drawstyle='steps-mid')
        
        # Plot a dummy line on secondary axis to later modify upper x-axis
        if specIdx == 0:
            ax2.plot(wls,[-0.5] * len(wls), color=WHITE)
            
        # Track the highest & lowest y-axis values to fix y-axis limits later            
        tmpMin = np.nanmin(fluxes)
        if tmpMin < minPlot:
            minPlot = tmpMin
        tmpMax = np.nanmax(fluxes + offset)
        if tmpMax > maxPlot:
            maxPlot = tmpMax
        
        # Track the highest y-axis value for the tail of each band
        currMax = np.nanmax(fluxes[-40:-1])
        if currMax > tailMax:
            tailMax = currMax
        
        # Add annotation to template plot
        textLoc = (0, 6)
        annotLoc = (wls[-40], tailMax + offset)
        annotTxt = objLabel
        subPlot.annotate(annotTxt, xy=annotLoc, xycoords='data', \
                         color=BLACK, xytext=textLoc, textcoords='offset points')
        
        # Fix axes limits ------------------------------------------------
        minPlot = minPlot - minPlot * 0.1
        maxOff = 0.02
        maxPlot = maxPlot + maxPlot * maxOff
        
        plt.ylim(ymin=minPlot, ymax=maxPlot)
        subPlot.set_xlim(xmin=limits[0], xmax=limits[1] * 1.001)
        ax2.set_xlim(xmin=limits[0], xmax=limits[1] * 1.001)
        
        # Customize y axis -----------------------------------------------
        subPlot.spines['left'].set_color('none')
        subPlot.spines['right'].set_color('none')
        subPlot.yaxis.set_ticks([])
        
        if wide:
            offset = offset + 3.
        else:
            offset = offset + 1.
    
    return fig

# ============================= PROCEDURE =====================================

# 1. LOAD RELEVANT MODULES ----------------------------------------------------
import matplotlib.pyplot as plt
import astrotools as at
import numpy as np

# 2. SET UP VARIABLES ---------------------------------------------------------
# emission line galaxies
# DATA = [['2243-1525','U12153'],
#         ['2015-1252','U11982'],
#         ['2349+1833','U21012']]
# KIND = 'emln_galaxies'
# galaxies
DATA = [['0256+2013','U20138.fits'],
        ['0256+1935','U20137.fits'],
        ['1721-0619','U20732.fits'],
        ['055741-1333','U10474.fits'],
        ['055742-1333','U10475.fits'],
        ['0558-1339','U10480.fits'],
        ['2240+3848','U13203.fits'],
        ['0515-0656','U20216.fits'],
        ['1357-3946','U20534_050323.fits'],
        ['0421+1528','U20186.fits'],
        ['0011+5149','U13005.fits'],
        ['1506+2759','U20599.fits']]
FILE_LBL = 'galaxies'
FOLDER_OUT = '/Users/alejo/KCData/Output/special/'
FOLDER_NIR = '/Users/alejo/KCData/NIR/'
BAND = 'NIR'
BAND_LIMS = [0.8, 2.4]

# 3. ORGANIZE TARGETS ---------------------------------------------------------
dataOrg = sorted(DATA)

tgtNames = []
fitsNames = []
for tgt in dataOrg:
    tgtNames.append(tgt[0])
    fitsNames.append(FOLDER_NIR + tgt[1])

# Break into several plots if too many objects
numData = len(dataOrg)
if numData > 6:
    numPlots = int(np.ceil(numData / 6.))
else:
    numPlots = 1

# 3. GET SPECTRA --------------------------------------------------------------
spectra = at.read_spec(fitsNames, aToMicron=True, negToZero=True, errors=False, plot=False)

# 4. CLEAN AND NORMALIZE SPECTRA ----------------------------------------------
spectraC = at.sel_band(spectra, BAND_LIMS)
spectraN = at.norm_spec(spectraC, BAND_LIMS)

# 6. PLOT ALL SPECTRA ---------------------------------------------------------
# Parameter to space spectra more
if FILE_LBL == 'emln_galaxies':
    broad = True
else:
    broad = False
start = 0
stop = 6
for plIdx in range(numPlots):
    if stop > numData:
        plSpecs = spectraN[start:]
        plNames = tgtNames[start:]
    else:
        plSpecs = spectraN[start:stop]
        plNames = tgtNames[start:stop]
    figure = plotspec(specData=plSpecs, limits=BAND_LIMS, objID=plNames, wide=broad)
    plt.savefig(FOLDER_OUT + 'notM' + '_' + FILE_LBL + str(plIdx + 1) + '.pdf', dpi=600)
    
    start = start + 6
    stop = stop + 6
