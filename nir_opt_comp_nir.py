''' 
The main() procedure plots normalized spectral data in the NIR band (sorted by J-K magnitudes) for a given spectral type.

NEEDED: 1) FILE_IN: ASCII tab-delimited txt file with data for each object
           (Access query is "nir_spex_prism_with_optical")
           (columns are specified under HDR_FILE_IN).
        2) FILE_IN_STD: ASCII tab-delimited txt file with data for standard NIR objects
           (columns are specified under HDR_FILE_IN_STD).
        3) EXCL_FILE: ASCII tab-delimite txt file with list of unums of objects to exclude
        4) FOLDER_ROOT: Folder containing (1)-(3) above and all .fits files (which are stored in two folders: OPT and NIR. It also contains (5) below.
        5) FOLDER_OUT: Folder to store output, within (4) above.

INPUT:  1) spInput: Spectral type to select (e.g. L0); it can also be a single object, identified by unum (e.g. U20268).
        2) grav: All young: y, Gamma: g, Beta: b, Field: f, All: leave blank.

OUTPUT: PDF file with four plots for selected spectral type.
'''

def addannot(specData, subPlot, bandName, classType):
# Adds annotations to indicate spectral absorption lines for optical and 
# near infrared spectra.
    
    from scipy.stats import nanmean
    
    # Initialize strings
    TXT_SIZE = 9
    H2O      = 'H' + '$\sf_2$' + 'O'
    CH4      = 'CH' + '$\sf_4$'
    COH2O    = 'CO + ' + H2O
    
    # Define the spectral lines to annotate
    if bandName == 'OPT':
        if classType <= 'L5':
            offVO = 0.4
        else:
            offVO = 3.5
        if classType == 'L4':
            offCrH = 60
        else:
            offCrH = 47
        
        ANNOT = [None] * 12
        # [Name, wl of absorption(s), offset of annotation from plot, absorption type]
        # Offset: For Line/Doublet, if < 0 then annot below line; 
        #         For Band, if < 1 then annot below line
        ANNOT[0]  = ['Li I',  0.6707,            20, 'Line']
        ANNOT[1]  = ['VO',   (0.7300,0.7550), offVO, 'Band']
        ANNOT[2]  = ['K I',  (0.7665,0.7699),    65, 'Doublet']
        ANNOT[3]  = ['Rb I',  0.7800,            60, 'Line']
        ANNOT[4]  = ['VO',   (0.7850,0.8000), offVO, 'Band']
        ANNOT[5]  = ['Rb I',  0.7948,            62, 'Line']
        ANNOT[6]  = ['Na I', (0.8178,0.8200),    45, 'Doublet']
        ANNOT[7]  = ['TiO',   0.8432,            35, 'Line']
        ANNOT[8]  = ['Cs I',  0.8521,           -25, 'Line']
        ANNOT[9]  = ['CrH',   0.8611,        offCrH, 'Line']
        ANNOT[10] = ['FeH',   0.8692,           -50, 'Line']
        ANNOT[11] = ['Cs I',  0.8943,           -30, 'Line']
    
    elif bandName == 'J':
        ANNOT = [None] * 12
        ANNOT[0]  = [H2O,   (0.890,0.990),  1.6, 'Band']
        ANNOT[1]  = ['FeH', (0.990,1.007),  0.7, 'Band']
        ANNOT[2]  = ['VO',  (1.050,1.080),  1.2, 'Band']
        ANNOT[3]  = [H2O,   (1.090,1.200),  0.4, 'Band']
        ANNOT[4]  = ['?',   (1.085,1.123),  1.1, 'Band']
        ANNOT[5]  = [CH4,   (1.100,1.240),  0.3, 'Band']
        ANNOT[6]  = ['Na I', 1.141,          40, 'Line']
        ANNOT[7]  = ['K I',  1.170,         -40, 'Line']
        ANNOT[8]  = ['VO',  (1.160,1.200), 1.15, 'Band']
        ANNOT[9]  = ['FeH', (1.194,1.239),  0.7, 'Band']
        ANNOT[10] = ['K I',  1.250,         -25, 'Line']
        ANNOT[11] = [H2O,   (1.300,1.390), 0.45, 'Band']
    
    elif bandName == 'H':
        if classType == 'L5':
            offH2O = 1.2
        else:
            offH2O = 1.46
        
        ANNOT = [None] * 4
        ANNOT[0] = [H2O,   (1.410,1.510),   1.35, 'Band']
        ANNOT[1] = ['K I',  1.517,           -35, 'Line']
        ANNOT[2] = ['FeH', (1.583,1.750),   0.65, 'Band']
        ANNOT[3] = [H2O,   (1.750,1.890), offH2O, 'Band']
    
    elif bandName == 'K':
        if classType <= 'L7':
            offCH4 = 1.05
        else:
            offCH4 = 0.87
        if classType == 'L4' or classType == 'L5':
            offCOH2O = 1.2
        else:
            offCOH2O = 1.3
            
        ANNOT = [None] * 4
        ANNOT[0] = [H2O,    (1.910,2.050),      0.8, 'Band']
        ANNOT[1] = [CH4,    (2.150,2.390),   offCH4, 'Band']
        ANNOT[2] = ['Na I',  2.210,             -25, 'Line']
        ANNOT[3] = [COH2O,  (2.293,2.390), offCOH2O, 'Band']
    
    else:
        return
    
    # Add annotation for each absorption band/line
    for annotation in ANNOT:
        
        # Determine distances between annotated point and annotation's objects
        offLine = annotation[2]     # Distance betw. annotation line & plot
        if offLine > 0:
            offText = offLine + 10  # Distance betw. text & plot
        else:
            offText = offLine - 15
        
        # Determine annotation line style
        if annotation[1] == 0.8943:  # Exception OPT-Band, Cs I
            annLineType = dict(arrowstyle='-', \
                               connectionstyle='angle,angleA=0,angleB=90,rad=0', \
                               shrinkB=offLine, shrinkA=0.5)
        else:
            annLineType = dict(arrowstyle='-', shrinkB=offLine, shrinkA=0.5)
        
        annotType = annotation[3]
        
        if annotType == 'Line':
        # For Line absorption: Add annotation with vertical connector
            
            # Initialize variables
            objsFluxIdxs = [numpy.nan] * len(specData)
            objsFluxes   = [numpy.nan] * len(specData)
            
            # Find spectrum with the highest/lowest flux @ absorption wl
            for objIdx, objSpec in enumerate(specData):
                wlRange = numpy.where(objSpec[0] <= annotation[1])
                if len(wlRange[0]) == 0:
                    objsFluxIdxs[objIdx] = 0
                else:
                    objsFluxIdxs[objIdx] = wlRange[0][-1]
                objsFluxes[objIdx] = objSpec[1][objsFluxIdxs[objIdx]]
            
            if offLine > 0:
                xtremeObj = numpy.array(objsFluxes).argmax()
                
            else:
                xtremeObj = numpy.array(objsFluxes).argmin()
            
            # Set the coordinate location for the annotated point
            annotWL  = specData[xtremeObj][0][objsFluxIdxs[xtremeObj]]
            annotLoc = (annotWL, objsFluxes[xtremeObj])
            
            # Set the coordinate location for the annotation's text
            if annotation[1] == 0.8943:  # Exception OPT-Band, Cs I
                textLoc = (-5, offText)
            else:
                textLoc = (0, offText)
            
            # Add annotation
            subPlot.annotate(annotation[0], xy=annotLoc, xycoords='data', \
                             xytext=textLoc, textcoords='offset points', \
                             fontsize=TXT_SIZE, ha='center', arrowprops=annLineType)
            
        elif annotType == 'Band':  # Draw a horizontal line
        # For band absorption: Add horizontal line AND annotation with no connector
            
            # Initialize variables
            objsFluxIdxs = [numpy.nan] * len(specData)
            objsFluxAvgs = [numpy.nan] * len(specData)
            xPos = numpy.zeros([len(specData),2])
            xPos.fill(numpy.nan)
            
            # Find spectrum with the highest/lowest flux average @ absorption wls
            for objIdx, objSpec in enumerate(specData):
                xLoRange = numpy.where(objSpec[0] <= annotation[1][0])
                if len(xLoRange[0]) == 0:
                    xPos[objIdx,0] = objSpec[0][0]
                else:
                    xPos[objIdx,0] = xLoRange[0][-1]
                
                xHiRange = numpy.where(objSpec[0] >= annotation[1][1])
                if len(xHiRange[0]) == 0:
                    xPos[objIdx,1] = objSpec[0][-1]
                else:
                    xPos[objIdx,1] = xHiRange[0][0]
                
                # Set up limits of section with which to calculate average flux
                if annotation[1] == (2.150,2.390): # Exception K-Band, CH4
                    if offCH4 > 1:
                        firstxPos = xPos[objIdx,0]
                        lastxPos  = xPos[objIdx,0] + \
                                    (xPos[objIdx,1] - xPos[objIdx,0]) * 1 / 4
                    else:
                        firstxPos = xPos[objIdx,0] + \
                                    (xPos[objIdx,1] - xPos[objIdx,0]) * 3 / 4
                        lastxPos  = xPos[objIdx,1]
                else:
                    firstxPos = xPos[objIdx,0]
                    lastxPos  = xPos[objIdx,1]
                
                objsFluxAvgs[objIdx] = nanmean(objSpec[1][firstxPos:lastxPos])
                
            if offLine > 1:
                textLoc   = (0,1)
                xtremeObj = numpy.array(objsFluxAvgs).argmax()
            else:
                textLoc   = (0,-11)
                xtremeObj = numpy.array(objsFluxAvgs).argmin()
            
            # Set the coordinate locations for the horizontal line
            # and the annotated point
            xMin   = specData[xtremeObj][0][xPos[xtremeObj][0]]
            xMax   = specData[xtremeObj][0][xPos[xtremeObj][1]]
            xMid   = xMin + (xMax - xMin) / 2
            annotY = objsFluxAvgs[xtremeObj] * offLine
            
            annotLoc = (xMid, annotY)
                
            # Add horizontal line AND annotation
            subPlot.plot([xMin,xMax],[annotY,annotY], color='k', \
                         linewidth=1, label='_ann')
            subPlot.annotate(annotation[0], xy=annotLoc, xycoords='data', \
                             xytext=textLoc, textcoords='offset points', \
                             fontsize=TXT_SIZE, ha='center')
        
        elif annotType == 'Doublet':  # Draw two vertical lines
        # For Doublet absorption: Add two annotations with vertial connectors
            
            # Find spectrum with the highest/lowest flux @ doublet's first absorption wl
            for objIdx, objSpec in enumerate(specData):
                wlRange = numpy.where(objSpec[0] <= annotation[1][0])
                if len(wlRange[0]) == 0:
                    objsFluxIdxs[objIdx] = 0
                else:
                    objsFluxIdxs[objIdx] = wlRange[0][-1]
                objsFluxes[objIdx] = objSpec[1][objsFluxIdxs[objIdx]]
            
            if offLine > 0:
                xtremeObj = numpy.array(objsFluxes).argmax()
            else:
                xtremeObj = numpy.array(objsFluxes).argmin()
            
            # Set the coordinate location of the first annotated point
            loc1      = numpy.where(specData[xtremeObj][0] <= annotation[1][0])
            annotLoc1 = (specData[xtremeObj][0][loc1[0][-1]], \
                         specData[xtremeObj][1][loc1[0][-1]])
            
            # Set the coordinate location of the first annotations' text
            txtLoc = (0, offText)
            
            # Add first annotation
            subPlot.annotate(annotation[0], xy=annotLoc1, xycoords='data', \
                             xytext=txtLoc, textcoords='offset points', \
                             fontsize=TXT_SIZE, ha='center', arrowprops=annLineType)
            
            # Set the coordinate location of the second annotated point
            loc2      = numpy.where(specData[xtremeObj][0] <= annotation[1][1])
            annotLoc2 = (specData[xtremeObj][0][loc2[0][-1]], annotLoc1[1])
            txtLoc    = (0,offText)
            
            # Add second annotation (with no text)
            subPlot.annotate(' ', xy=annotLoc2, xycoords='data', xytext=txtLoc, \
                             textcoords='offset points', ha='center', \
                             arrowprops=annLineType)
    return


def plotspec(specData, bandNames, limits, objID, classType, grav=None, plotInstructions=None, figNum=1):
# Plots set of spectral data and saves plots in a PDF file.
# specData and limits must be dictionaries.
    
    import matplotlib.pyplot as plt
    import types
    import numpy
    
    # 1) Check data consistency ===============================================
    # Stop if specData or limits are not dictionaries
    try:
        specData.keys()
        limits.keys()
    except AttributeError:
        print 'PLOTSPEC: Data not received as dictionaries.'
        return
    
    # 2) Initialize variables & color sets (hex codes) ========================
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
    
    colors = [None] * 31
    colors[30] = COLOR_SET.copy().tolist()
    colors[29] = COLOR_SET[[0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, \
                            19,20,21,22,23,24,25,26,27,28,29]].tolist()
    colors[28] = COLOR_SET[[0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, \
                            19,20,21,23,24,25,26,27,28,29]].tolist()
    colors[27] = COLOR_SET[[0,1,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18, \
                            19,20,21,23,24,25,26,27,28,29]].tolist()
    colors[26] = COLOR_SET[[0,1,3,4,5,6,7,8,10,11,12,14,15,16,17,18, \
                            19,20,21,23,24,25,26,27,28,29]].tolist()
    colors[25] = COLOR_SET[[1,3,4,5,6,7,8,10,11,12,14,15,16,17,18, \
                            19,20,21,23,24,25,26,27,28,29]].tolist()
    colors[24] = COLOR_SET[[1,3,4,5,6,7,8,10,11,12,14,15,16,17,18, \
                            19,20,21,24,25,26,27,28,29]].tolist()
    colors[23] = COLOR_SET[[1,3,4,5,6,7,8,10,11,12,14,15,16,18, \
                            19,20,21,24,25,26,27,28,29]].tolist()
    colors[22] = COLOR_SET[[1,3,4,6,7,8,10,11,12,14,15,16,18, \
                            19,20,21,24,25,26,27,28,29]].tolist()
    colors[21] = COLOR_SET[[1,3,4,6,7,8,10,11,12,14,15,16,18, \
                            19,20,21,24,25,27,28,29]].tolist()
    colors[20] = COLOR_SET[[1,3,4,6,7,8,10,11,12,14,16,18, \
                            19,20,21,24,25,27,28,29]].tolist()
    colors[19] = COLOR_SET[[1,3,4,6,7,8,10,11,12,16,18, \
                            19,20,21,24,25,27,28,29]].tolist()
    colors[18] = COLOR_SET[[1,3,4,6,7,10,11,12,16,18, \
                            19,20,21,24,25,27,28,29]].tolist()
    colors[17] = COLOR_SET[[1,3,4,6,7,10,11,12,18,19,20,21,24,25,27,28,29]].tolist()
    colors[16] = COLOR_SET[[1,3,6,7,10,11,12,18,19,20,21,24,25,27,28,29]].tolist()
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
    colors[1]  = COLOR_SET[[29]].tolist()
    
    BLACK = '#000000'
    GRAY  = '#999999'
    WHITE = '#FFFFFF'
    X_LABEL = 'Wavelength ($\mu$m)'
    Y_LABEL = 'Normalized Flux (F$_{\lambda}$)'
    
    # 3) Initialize Figure ====================================================
    plt.close()
    plt.rc('font', size=10)
    fig = plt.figure(figNum, figsize=(9,6))
    plt.clf()
    
    # 4) Generate Subplots ====================================================
    for bandIdx, band in enumerate(bandNames):
        
        # 4a) If band data is only one set, convert into array of sets --------
        if specData[band][0] is not None:
            if len(specData[band][0]) > 3:
                specData[band] = [specData[band],]
        
        # 4b) Initialize variables --------------------------------------------
        spLines  = []
        minPlot  = 1
        maxPlot  = 1
        
        # Count the number of plots in order to select color set
        if plotInstructions is not None:
            tmpFld  = numpy.where(numpy.array(plotInstructions) == 'field')
            tmpYng  = numpy.where(numpy.array(plotInstructions) == 'young')
            numFld  = len(tmpFld[0])
            numYng  = len(tmpYng[0])
            specNum = numFld + numYng
        else:
            specNum  = len(filter(None, specData[band]))
        
        # Select color set based on count above
        if specNum > len(COLOR_SET):
            plotColors = colors[len(COLOR_SET)][:]
        elif specNum == 0:
            plotColors = None
        else:
            plotColors = colors[specNum][:]
        
        
        # 4c) Initialize Subplot ----------------------------------------------
        # Determine position of Subplot
        multH = 0
        multV = 1
        
        gapHoriz   = 0.03
        gapVertic  = 0.04
        plotHeight = 0.90  # proportion of total height (11 inches)
        plotWidth  = 0.92  # proportion of total width (8.5 inches)
        edgeLeft   = 0.06 #+ (plotWidth + gapVertic) * multH
        edgeBottom = 0.08 #+ (plotHeight + gapHoriz) * multV
        plotLoc    = [edgeLeft, edgeBottom, plotWidth, plotHeight]
        
        subPlot = plt.figure(figNum).add_axes(plotLoc)
        
        subPlot.set_autoscale_on(False)
        
        # Set figure and axes labels
        if grav == 'Y':
            plotType = ' young'
        elif grav == 'B':
            plotType = r'$\beta$'
        elif grav == 'G':
            plotType = r'$\gamma$'
        elif grav == 'F':
            plotType = ' field'
        else:
            plotType = ''
        title = classType
        if plotType != '':
            title = title + plotType
        
        subPlot.set_xlabel(X_LABEL)
        subPlot.set_ylabel(Y_LABEL)
        subPlot.set_title(title, fontsize=20, fontweight='bold', \
                          position=(0.01,0.9), ha='left')
        
        # 4d) Determine order of spectra plotting -----------------------------
        zOrders = [None] * len(plotInstructions)
        countColor = specNum
        for plotIdx,plot in enumerate(plotInstructions):
            zOrders[plotIdx] = specNum - countColor
            countColor = countColor - 1
        
        # 4e) Plot spectral data in Subplot -----------------------------------
        countColors = specNum - 1
        textColors = [] # For legend purposes only
        for specIdx, spec in enumerate(specData[band]):
            if spec is None:
                continue
            if plotInstructions[specIdx] == 'exclude':
                continue
            
            # Determine line parameters
            lnStyle = '-'
            lnWidth = 0.5
            objLabel  = objID[specIdx]
            
            # Consolidate color plot and legend designation
            plotColor   = plotColors[countColors] # Color for plot line
            legColor    = plotColor               # Color for legend text
            countColors = countColors - 1
            
            textColors.append(legColor) # Colors for legend labels
            
            subPlot.plot(spec[0], spec[1], color=plotColor, linestyle=lnStyle, \
                        dash_joinstyle='round', linewidth=lnWidth, label=objLabel, \
                        drawstyle='steps-mid', zorder=zOrders[specIdx])
            
            # Track the highest & lowest y-axis values to fix y-axis limits later            
            if plotInstructions[specIdx] != 'exclude':
                tmpMin = numpy.nanmin(spec[1])
                if tmpMin < minPlot:
                    minPlot = tmpMin
                tmpMax = numpy.nanmax(spec[1])
                if tmpMax > maxPlot:
                    maxPlot = tmpMax
        
        # 4f) Fix axes limits -------------------------------------------------
        minPlot = minPlot - minPlot * 0.1
        maxOff = 0.15
        maxPlot = maxPlot + maxPlot * maxOff
        plt.ylim(ymin=minPlot, ymax=maxPlot)
        plt.xlim(xmin=limits[band]['lim'][0], \
                 xmax=limits[band]['lim'][1] * 1.001)
        
        # 4g) Customize y axis ------------------------------------------------
        subPlot.spines['left'].set_color('none')
        subPlot.spines['right'].set_color('none')
        subPlot.yaxis.set_ticks([])
        
        # 4h) Create and format legend
        objLegends = plt.legend(handlelength=0, handletextpad=0.1, loc='upper right', \
                                bbox_to_anchor=(1.012,0.97), labelspacing=0.04, \
                                borderpad=0.2, numpoints=1)
        objLegends.draw_frame(True)
        
        for legendIdx, legendText in enumerate(objLegends.get_texts()):
            plt.setp(legendText, color=textColors[legendIdx], \
                     fontsize=7, fontname='Andale Mono')
        
        # Add Titles for the legends
        legendTitles1 = '                  Optical'
        legendTitles2 = 'Coords.      SpType   J-K'
        xCoord = 0.87
        yCoord1 = 0.98
        yCoord2 = 0.96
        
        subPlot.text(xCoord, yCoord1, legendTitles1, fontsize=6, \
                     transform=subPlot.transAxes, zorder=20)
        subPlot.text(xCoord, yCoord2, legendTitles2, fontsize=6, \
                     transform=subPlot.transAxes, zorder=20)
        
        # 4i) Add absorption annotations to Subplots
        #addannot(filter(None, specData[band]), subPlot, band, classType)
        
    return fig


def main(spInput, grav=''):
    # 1. LOAD RELEVANT MODULES ---------------------------------------------------------
    import astrotools as at
    import asciidata
    import pyfits
    import matplotlib.pyplot as plt
    import numpy
    import sys
    import pdb
    
    
    # 2. SET UP VARIABLES --------------------------------------------------------------
    FOLDER_ROOT = '/Users/alejo/KCData/'  # Location of NIR and OPT folders
    FOLDER_OUT  = 'Output/NOCN/'
    OPTNIR_KEYS  = ['OPT', 'NIR']
    BAND_NAME  = ['NIR']
    data       = ''
    dataRaw    = ''
    specFiles  = ''
    spectraRaw = ''
    spectra    = ''
    
    # For TXT objects file (updatable here directly)
    FILE_IN     = 'nir_spex_prism_with_optical_12aug15.txt' # ASCII file w/ data
    HDR_FILE_IN = ('Ref','Designation','J','H','K','SpType','SpType_T','NIRFobs',\
                   'NIRFtel','NIRfile','OPTobs','OPTtel','OPTinst','OPTfile',\
                   'Young?','Dusty?','Blue?','Multiple?','Pec?')
    
    colNameRef   = HDR_FILE_IN[0]
    colNameDesig = HDR_FILE_IN[1]
    colNameJ     = HDR_FILE_IN[2]
    colNameK     = HDR_FILE_IN[4]
    colNameJK    = 'J-K'
    colNameType  = HDR_FILE_IN[6]
    colNameYng   = HDR_FILE_IN[14]
    colNameDust  = HDR_FILE_IN[15]
    colNameBlue  = HDR_FILE_IN[16]
    colNamePec   = HDR_FILE_IN[18]
    
    # For TXT exclude-objects file
    EXCL_FILE = 'Exclude_Objs.txt'   # ASCII file w/ U#s of objects to exclude
    
    
    # 3. READ DATA FROM INPUT FILES ----------------------------------------------------
    NULL_CHAR = ''   # Null character
    DELL_CHAR = '\t' # Delimiter character
    COMM_CHAR = '#'  # Comment character
    
    # File with objects (query in Access)
    dataRaw = asciidata.open(FOLDER_ROOT + FILE_IN, NULL_CHAR, DELL_CHAR, COMM_CHAR)
    
    # Store data in a dictionary-type object
    data = {}.fromkeys(HDR_FILE_IN)
    for colIdx,colData in enumerate(dataRaw):
        data[HDR_FILE_IN[colIdx]] = colData.tonumpy()
    
    
    # 4. FORMAT SOME ASCII COLUMNS -----------------------------------------------------
    # 4.1 Convert into unicode the Spectral Type-Text column
    uniSpType = [None] * len(data[colNameType])
    for sIdx,sType in enumerate(data[colNameType]):
        uniSpType[sIdx] = sType.decode('utf-8')
    
    data[colNameType] = numpy.array(uniSpType)
    
    # 4.2 Calculate J-K Color And Add J-K Column
    data[colNameJK] = data[colNameJ] - data[colNameK]
    
    # 4.3 Format Designation Number from Designation Column
    for desigIdx,desig in enumerate(data[colNameDesig]):
        desig = ''.join(desig.split())
        signType = '+'
        signPos = desig.find(signType)
        if signPos == -1:
            signType = '-'
            signPos  = desig.find(signType)
        
        desigProper = desig[:4] + signType + desig[signPos+1:signPos+5]
        data[colNameDesig][desigIdx] = desigProper
    
    
    # 5. FILTER DATA BY USER INPUT IN spInput ------------------------------------------
    # Find all spectra of same spectral type
    specIdx = []
    for spIdx,spType in enumerate(data[colNameType]):
        if spType.upper().startswith(spInput.upper()):
            specIdx.append(spIdx)
    
    if not specIdx:
        print 'No target found for given input.'
        return
    spTypeInput = spInput.upper()
    
    # Sort relevant objects by JKmag value
    specIdx     = numpy.array(specIdx)
    specSortIdx = data[colNameJK][specIdx].argsort()
    
    
    # 6. READ SPECTRAL DATA FROM SPECTRAL FILES ----------------------------------------
    spectraRaw    = {}.fromkeys(OPTNIR_KEYS) # Used to store the raw data from fits files
    specFilesDict = {}.fromkeys(OPTNIR_KEYS) # Used for reference purposes
    
    for key in OPTNIR_KEYS:
        specFiles = [None] * len(specSortIdx)
        
        for sortIdx,specSort in enumerate(specSortIdx):
            tmpFullName = FOLDER_ROOT + key + '/' + data[key + 'file'][specIdx[specSort]]
            specFiles[sortIdx] = tmpFullName
            specFilesDict[key] = specFiles
        
        spectraRaw[key] = at.read_spec(specFiles, atomicron=True, negtonan=True, \
                                       errors=True, verbose=False)
    
    # Clear out spectral data for objects missing either OPT or NIR data
    allNone = True
    for spIdx in range(0,len(spectraRaw['OPT'])):
        if spectraRaw['OPT'][spIdx] is None:
            spectraRaw['NIR'][spIdx] = None
        elif spectraRaw['NIR'][spIdx] is None:
            spectraRaw['OPT'][spIdx] = None
        else:
            allNone = False
    
    if allNone:
        print 'No spectral data found for objects of the given spectral type.'
        return
    
    # Convert spectraRaw contents into lists if only one spectral data
    for key in spectraRaw.keys():
        if spectraRaw[key][0] is not None:
            if len(spectraRaw[key][0]) > 3:
                spectraRaw[key] = [spectraRaw[key],]
    
    
    # 7. GATHER OBJECTS' NAMES----------------------------------------------------------
    # Filtered objects
    refs = [None] * len(specSortIdx)
    for idx,spIdx in enumerate(specSortIdx):
        tmpRef    = data[colNameRef][specIdx[spIdx]]
        refs[idx] = str(int(tmpRef))
    
    
    #8. SMOOTH SPECTRA -----------------------------------------------------------------
    # Smooth the flux data to a reasonable resolution
    spectraS = at.smooth_spec(spectraRaw['NIR'], specFile=specFilesDict['NIR'], \
                              winWidth=0)
    
    
    # 9. SET LIMITS FOR BAND AND NORMALIZING SECTION------------------------------------
    # Initialize dictionary to store limits
    BAND_LIMS = {}.fromkeys(BAND_NAME)
    for bandKey in BAND_NAME:
        BAND_LIMS[bandKey] = dict(lim = [None] * 2, limN = [None] * 2)
    
    # Set wl limits for band
    # Limits are in microns
    BAND_LIMS['NIR']['lim'][0] = 0.8
    BAND_LIMS['NIR']['lim'][1] = 2.4
    
    # Set wl limits for normalizing sections; this is the peak of the J band
    # Limits are in microns
    BAND_LIMS['NIR']['limN'][0] = 1.28
    BAND_LIMS['NIR']['limN'][1] = 1.32
    
    
    # 10. SELECT SPECTRAL DATA FOR NIR BAND---------------------------------------------
    # Initialize variables
    spectraN = {}.fromkeys(BAND_NAME)
    
    # Gather reference numbers of objects
    objRef = data[colNameRef][specIdx[specSortIdx]]
    
    # Select band
    spectra = at.sel_band(spectraS, BAND_LIMS['NIR']['lim'], objRef)
    
    # Normalize band
    spectraN['NIR'] = at.norm_spec(spectra, BAND_LIMS['NIR']['limN'])
    
    
    # 11. CHARACTERIZE TARGETS (i.e. identify young, blue, to exclude...)---------------
    # Determine which targets to exclude using the "Exclude_Objs" file
    toExclude = [False] * len(refs)
    dataExcl = asciidata.open(FOLDER_ROOT + EXCL_FILE, NULL_CHAR, DELL_CHAR, COMM_CHAR)
    if len(dataExcl[0]) > 0:
        # Extract data from "Exclude_Objs" file
        excludeObjs = [None] * len(dataExcl[0])
        for rowIdx, rowData in enumerate(dataExcl[0]):
            excludeObjs[rowIdx] = str(rowData)
        
        # Find intersection of exclude-obj list and filtered targets list
        setExclude = set(excludeObjs).intersection(set(refs))
        
        # Create list with intersection targets
        if len(setExclude) != 0:
            for exclIdx in setExclude:
                tmpExclIdx = numpy.where(numpy.array(refs) == exclIdx)
                toExclude[tmpExclIdx[0]] = True
    
    # Determine which targets are blue
    blueObjs = [False] * len(refs)
    for idx,spIdx in enumerate(specIdx[specSortIdx]):
        if data[colNameBlue][spIdx].upper() == 'YES':
            blueObjs[idx] = True
    
    # Determine which targets are dusty
    dustyObjs = [False] * len(refs)
    for idx,spIdx in enumerate(specIdx[specSortIdx]):
        if data[colNameDust][spIdx].upper() == 'YES':
            dustyObjs[idx] = True
    
    # Determine which targets are peculiar
    pecObjs = [False] * len(refs)
    for idx,spIdx in enumerate(specIdx[specSortIdx]):
        if data[colNamePec][spIdx].upper() == 'YES':
            pecObjs[idx] = True
    
    # Determine which plots are young objects
    youngObjs = [False] * len(refs)
    for idx,spIdx in enumerate(specSortIdx):
        if data[colNameYng][specIdx[spIdx]].upper() == 'YES':
            youngObjs[idx] = True
    
    # Determine which targets are GAMMA
    gammaObjs = [False] * len(refs)
    for idx,spIdx in enumerate(specIdx[specSortIdx]):
        tmpType = data[colNameType][spIdx].encode('utf-8')
        tmpLen  = len(tmpType)
        utcA = tmpType[tmpLen - 2]
        utcB = tmpType[tmpLen - 1]
        # GAMMA in utf-8 code is "\xce\xb3"
        if utcA == '\xce' and utcB == '\xb3':
            gammaObjs[idx] = True
    
    # Determine which targets are BETA
    betaObjs = [False] * len(refs)
    for idx,spIdx in enumerate(specIdx[specSortIdx]):
        tmpType = data[colNameType][spIdx].encode('utf-8')
        tmpLen  = len(tmpType)
        utcA = tmpType[tmpLen - 2]
        utcB = tmpType[tmpLen - 1]
        # GAMMA in utf-8 code is "\xce\xb2"
        if utcA == '\xce' and utcB == '\xb2':
            betaObjs[idx] = True
    
    # Determine which targets to include in plots (based on user input)
    # Consolidate plotting instructions
    grav = grav.upper()
    plotInstructions = ['exclude'] * len(refs)
    if grav == 'Y': # If plot request is Young, include gamma, beta & young targets
        for plotIdx in range(len(refs)):
            if toExclude[plotIdx]:
                continue
            if gammaObjs[plotIdx] or betaObjs[plotIdx] or youngObjs[plotIdx]:
                if blueObjs[plotIdx] or dustyObjs[plotIdx] or pecObjs[plotIdx]:
                    continue
                plotInstructions[plotIdx] = 'young'
    
    elif grav == 'G': # If plot request is Gamma, include only gamma targets
        for plotIdx in range(len(plotInstructions)):
            if toExclude[plotIdx]:
                continue
            if gammaObjs[plotIdx]:
                if blueObjs[plotIdx] or dustyObjs[plotIdx] or pecObjs[plotIdx]:
                    continue
                plotInstructions[plotIdx] = 'young'
    
    elif grav == 'B': # If plot request is Beta, include only beta targets
        for plotIdx in range(len(plotInstructions)):
            if toExclude[plotIdx]:
                continue
            if betaObjs[plotIdx]:
                if blueObjs[plotIdx] or dustyObjs[plotIdx] or pecObjs[plotIdx]:
                    continue
                plotInstructions[plotIdx] = 'young'
    
    elif grav == 'F': # If plot request is Field, include Field & Standard targets
        for plotIdx in range(len(plotInstructions)):
            if toExclude[plotIdx]:
                continue
            if betaObjs[plotIdx] or gammaObjs[plotIdx] or youngObjs[plotIdx]:
                continue
            if blueObjs[plotIdx] or dustyObjs[plotIdx] or pecObjs[plotIdx]:
                continue
            plotInstructions[plotIdx] = 'field'
    
    else:   # Otherwise, print Field, gamma, beta, young & Standard targets
        for plotIdx in range(len(plotInstructions)):
            if toExclude[plotIdx]:
                continue
            if blueObjs[plotIdx] or dustyObjs[plotIdx] or pecObjs[plotIdx]:
                continue
            if youngObjs[plotIdx]:
                plotInstructions[plotIdx] = 'young'
            else:
                plotInstructions[plotIdx] = 'field'
    
    # If all plot instructions are "exclude", then stop procedure
    allExcl = True
    for instr in plotInstructions:
        if instr != 'exclude':
            allExcl = False
    if allExcl:
        if not uniqueSpec:
            print 'No spectral data to plot based on your request.'
            return
    
    
    # 12. PLOT DATA --------------------------------------------------------------------
    # Gather info on each object (for legend purposes)
    objInfo = [None] * len(refs)
    for posIdx,spIdx in enumerate(specIdx[specSortIdx]):
        tmpDesig  = data[colNameDesig][spIdx]
        tmpJK     = data[colNameJK][spIdx]
        tmpSPtype = data[colNameType][spIdx]
        tmpSPtype = tmpSPtype + ' ' * (5 - len(tmpSPtype))  # For alignment purposes
    
        objInfo[posIdx] = (tmpDesig + ' ' + tmpSPtype + ' ' + '%.2f' %tmpJK)
    
    # Create Figure with Subplots
    figObj = plotspec(spectraN, BAND_NAME, BAND_LIMS, objInfo, spTypeInput, grav, \
                        plotInstructions)
    
    figObj.savefig(FOLDER_ROOT + FOLDER_OUT + spTypeInput + grav + '_fan.pdf', \
                   dpi=800)

