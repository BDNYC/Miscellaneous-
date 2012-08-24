''' 
This procedure reads and calculates J-K colors for all objects in Paper XII, and then plot their average and extreme points as a function of spectral type. It plots them together with Faherty-12 average and extreme points.
Until we get a NIR spectrum of U50171 (0835+1953, L5 standard) with uncertainties, it will be manually excluded from this plot in step (10) of the procedure, since it is not used to calculate the L5 template.
'''

def makeplot(dataEB1, dataEB2):
    
    import numpy
    import matplotlib.pyplot as plt
    import pdb
    
    try:
        dataEB1.shape
        dataEB2.shape
    except AttributeError:
        return
    
    # 1) Initialize variables
    BLACK = '#000000'
    GRAY  = '#999999'
    RED   = '#FF0000'
    WHITE = '#FFFFFF'
    X_LABEL = 'Spectral Type'
    Y_LABEL = 'J-K$_s$ (2MASS)'
    X_OFFSET = 0.08
    
    # 2) Initialize Figure
    plt.close()
    plt.rc('font', size=12)
    fig = plt.figure(1, figsize=(6,5))
    plt.clf()
    
    # 3) Create bins for Kelle's data
    bins = []
    labels = []
    for tp in range(11, 20):
        idxHi = np.where(dataEB1[0,:] < tp)
        idxLo = np.where(dataEB1[0,:] >= tp - 1)
        if len(idxHi[0]) > 0 and len(idxLo[0]) > 0:
            idxtp = np.intersect1d_nu(idxHi[0], idxLo[0])
        else:
            idxtp = []
        if len(idxtp) != 0:
            bins.append(dataEB1[1,idxtp])
            labels.append('L' + str(tp - 11))
        # Rememeber the smallest J-K value of L6 to annotate label there later
        if tp == 16:
            lblymin = np.min(dataEB1[1,idxtp])        
    
    # 4) Initialize subplot
    ax = fig.add_subplot(1,1,1, position=[0.12,0.11,0.85,0.84])
                                        # [left,bottom,width,height]
    ax.set_autoscale_on(False)
    
    # 5) Plot the dots for Jackie's data
    xs = np.array(range(1, 10))
    avgs2 = dataEB2[0,:]
    extremes2 = np.row_stack((avgs2 - dataEB2[1,:], dataEB2[2,:] - avgs2))
    counts2 = dataEB2[3,:]
    ebars2 = plt.errorbar(xs - X_OFFSET, avgs2, yerr=extremes2, ecolor='r', \
                         linestyle='', marker='o', markersize=10, markerfacecolor='r', \
                         markeredgecolor='r', markeredgewidth=1.1, label='Faherty-12')
    
    # Annotate top of bars
    for idxc, count in enumerate(counts2):
        loc = (xs[idxc] - 0.08, dataEB2[2,idxc])
        loctext = (-7, 5)
        linetype = dict(arrowstyle='-', shrinkB=4, shrinkA=2, color=GRAY) 
        plt.annotate(str(int(count)), xy=loc, xytext=loctext, xycoords='data', \
                     textcoords='offset points', fontsize=8, fontstyle='italic', \
                     color=GRAY)#, arrowprops=linetype)
    
    # 6) Plot the dots for Kelle's data
    avgs1 = []
    maxs1 = []
    mins1 = []
    counts1 = []
    for bn in bins:
        avg = np.average(bn)
        max1 = np.max(bn)
        min1 = np.min(bn)
        avgs1.append(avg)
        maxs1.append(max1)
        mins1.append(min1)
        counts1.append(len(bn))
    avgs1 = np.array(avgs1)
    mins1 = np.array(mins1)
    maxs1 = np.array(maxs1)
    extremes1 = np.row_stack((avgs1 - mins1, maxs1 - avgs1))
    ebars1 = plt.errorbar(xs + X_OFFSET, avgs1, yerr=extremes1, ecolor='k', \
                         linestyle='', marker='s', markersize=8, markerfacecolor='k', \
                         markeredgecolor='k', markeredgewidth=1.1, label='Paper XII')
    
    # Annotate top of bars
    for idxc, count in enumerate(counts1):
        loc = (xs[idxc] + 0.08, maxs1[idxc])
        loctext = (-1, 5)
        linetype = dict(arrowstyle='-', shrinkB=4, shrinkA=2, color=GRAY) 
        plt.annotate(str(int(count)), xy=loc, xytext=loctext, xycoords='data', \
                     textcoords='offset points', fontsize=8, fontstyle='italic', \
                     color=GRAY)#, arrowprops=linetype)
    
    # 7) Format plot (make lots of things disappear)
    # Hide axes lines
    ax.spines['top'].set_color(WHITE)
    ax.spines['bottom'].set_color(WHITE)
    ax.spines['left'].set_color(WHITE)
    ax.spines['right'].set_color(WHITE)
    
    # Hide ticker lines
    tcks = ax.xaxis.get_ticklines()
    for tl in tcks:
        tl.set_color(WHITE)
    tcks = ax.yaxis.get_ticklines()
    for tl in tcks:
        tl.set_color(WHITE)
    
    # Add axes labels
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_ylim(0.9, 2.3)
    ax.set_xlim(0.5, 9.5)
    
    # Add tick labels
    plt.xticks(range(1,10), labels, size=14)
    
    # Add horizontal grid
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.8)
    ax.set_axisbelow(True)
    
    # 8) Annotate Data Labels
    # Annotate Jackie data set
    loc = (6. - X_OFFSET, dataEB2[1,6]) # The min J-K value at L6
    loctext = (-10, -65)
    linetype = dict(arrowstyle='-', shrinkB=4, shrinkA=2, color=RED, relpos=(1,0))
    plt.annotate('Faherty et al. \'12', xy=loc, xytext=loctext, xycoords='data', \
                 textcoords='offset points', fontsize=13, color=RED, \
                 ha='right', va='bottom', arrowprops=linetype) 
    
    # Annotate Paper XII data set
    loc = (6. + X_OFFSET, lblymin) # The min J-K value at L6
    loctext = (10, -60)
    linetype = dict(arrowstyle='-', shrinkB=4, shrinkA=2, color=BLACK, relpos=(0,0))
    plt.annotate('Paper XII', xy=loc, xytext=loctext, xycoords='data', \
                 textcoords='offset points', fontsize=13, color=BLACK, \
                 ha='left', va='bottom', arrowprops=linetype)
    
    return fig


# 1. LOAD RELEVANT MODULES ---------------------------------------------------------
import asciidata
import astrotools as at
import matplotlib.pyplot as plt
import numpy as np
import sys
import pdb

# 2. SET UP VARIABLES --------------------------------------------------------------
# Faherty-12 J-K averages, extremes, and counts, L0 to L8 in sequential order
AVGJK    = np.array([1.3,1.35,1.48,1.64,1.69,1.72,1.84,1.75,1.85])
AVGJKMIN = np.array([1.,1.,1.2,1.2,1.3,1.4,1.4,1.4,1.4])
AVGJKMAX = np.array([1.8,1.8,1.8,2.,2.,2.,2.,2.,2.])
COUNTS   = np.array([10,10,10,10,10,10,10,10,10])

# General variables
GRAV = 'f'
FOLDER_ROOT = '/Users/alejo/KCData/'  # Location of NIR and OPT folders
FOLDER_OUT  = 'Output/special/'
dataRaw = ''

# For TXT objects file (updatable here directly)
FILE_EXCL = 'Exclude_Objs.txt'   # ASCII file w/ objects to exclude
FILE_IN = 'NIR_Spex_Prism_with_optical_12Aug15.txt' # ASCII file w/ all objects
HDR_FILE_IN = ('Ref','Designation','J','H','K','SpType','SpType_T','NIRFobs',\
               'NIRFtel','NIRfile','OPTobs','OPTtel','OPTinst','OPTfile',\
               'Young?','Dusty?','Blue?', 'Multiple?','Pec?')
colRef = HDR_FILE_IN[0]
colJ   = HDR_FILE_IN[2]
colK   = HDR_FILE_IN[4]
colJK  = 'J-K'
colTypeN = HDR_FILE_IN[5]
colType  = HDR_FILE_IN[6]
colYng  = HDR_FILE_IN[14]
colDust = HDR_FILE_IN[15]
colBlue = HDR_FILE_IN[16]
colBin  = HDR_FILE_IN[17]
colPec  = HDR_FILE_IN[18]

# 3. READ DATA FROM INPUT FILES-----------------------------------------------------
NULL_CHAR = ''   # Null character
DELL_CHAR = '\t' # Delimiter character
COMM_CHAR = '#'  # Comment character

# File with objects (query in Access)
dataRaw = asciidata.open(FOLDER_ROOT + FILE_IN, NULL_CHAR, DELL_CHAR, COMM_CHAR)

# Store data in a dictionary-type object
dataDict = {}.fromkeys(HDR_FILE_IN)
for colIdx,colData in enumerate(dataRaw):
    dataDict[HDR_FILE_IN[colIdx]] = colData.tonumpy()

numRows = len(dataDict[colRef])

# 4. FORMAT SOME ASCII COLUMNS -----------------------------------------------------
# Convert into unicode the Spectral Type-Text column
uniSpType = [None] * numRows
for sIdx,sType in enumerate(dataDict[colType]):
    uniSpType[sIdx] = sType.decode('utf-8')
dataDict[colType] = np.array(uniSpType)

# Calculate J-K Color And Add J-K Column
dataDict[colJK] = dataDict[colJ] - dataDict[colK]

# 5. CREATE PYTHON LIST WITH RELEVANT INFO ON OBJECTS ------------------------------
dataLs = [dataDict[colRef], dataDict[colJK], dataDict[colType], \
          dataDict[colTypeN], dataDict[colDust], dataDict[colBlue], \
          dataDict[colBin], dataDict[colPec]]

# 6. SELECT TARGETS BASED ON GRAVITY PARAMETER -------------------------------------
# Determine which targets are young
gravObjs = ['f'] * numRows
for objIdx, obj in enumerate(dataLs[2]):
    tmpType = dataLs[2][objIdx].encode('utf-8')
    tmpLen  = len(tmpType)
    utcA = tmpType[tmpLen - 2]
    utcB = tmpType[tmpLen - 1]
    # GAMMA in utf-8 code is "\xce\xb3"
    if utcA == '\xce' and utcB == '\xb3':
        gravObjs[objIdx] = 'g'
    # BETA in utf-8 code is "\xce\xb2"
    elif utcA == '\xce' and utcB == '\xb2':
        gravObjs[objIdx] = 'b'

# Determine which targets to include in Kelle's data set
inclEB1 = np.array([False] * numRows)
if GRAV.upper() == 'G':
    incIdx = np.where(np.array(gravObjs) == 'g')
    if len(incIdx[0]) != 0:
        inclEB1[incIdx] = True
elif GRAV.upper() == 'B':
    incIdx = np.where(np.array(gravObjs) == 'b')
    if len(incIdx[0]) != 0:
        inclEB1[incIdx] = True
elif GRAV.upper() == 'F':
    incIdx = np.where(np.array(gravObjs) == 'f')
    if len(incIdx[0]) != 0:
        inclEB1[incIdx] = True

# 7. EXCLUDE OBJECTS BASED ON EXCLUDE.TXT FILE -------------------------------------
allRefs = dataLs[0]
dataExcl = asciidata.open(FOLDER_ROOT + FILE_EXCL, NULL_CHAR, DELL_CHAR, COMM_CHAR)
if len(dataExcl[0]) > 0:
    # Extract data from "Exclude_Objs" file
    excludeObjs = [None] * len(dataExcl[0])
    for rowIdx, rowData in enumerate(dataExcl[0]):
        excludeObjs[rowIdx] = rowData
    
    # Find intersection of exclude-obj list and all objects
    setExclude = set(excludeObjs).intersection(set(allRefs))
    
    # Exclude objects if present in exclude txt file
    if len(setExclude) != 0:
        for exclIdx in setExclude:
            tmpExclIdx = np.where(np.array(allRefs) == exclIdx)
            inclEB1[tmpExclIdx] = False

# 8. EXCLUDE BINARIES --------------------------------------------------------------
for binIdx, bin in enumerate(dataLs[6]):
    if bin.upper() == 'YES':
        inclEB1[binIdx] = False

# 9. EXCLUDE SPECIAL OBJECTS FROM THE KELLE'S DATA SET -----------------------------
for refIdx, ref in enumerate(dataLs[0]):
    # Exclude dusty, blue, and pec objects
    if dataLs[4][refIdx] == 'Yes' or dataLs[5][refIdx] == 'Yes' or \
                                                dataLs[7][refIdx] == 'Yes':
        inclEB1[refIdx] = False

# 10. MANUALLY EXCLUDE U50171 (0835+1953) (L5 STANDARD, MISSING UNCERTAINTIES) -----
for refIdx, ref in enumerate(dataLs[0]):
    if ref == 50171:
        inclEB1[refIdx] = False

# 11. PREPARE DATA TO PLOT ----------------------------------------------------------
toplotEB1 = np.array([np.array(dataLs[3])[inclEB1], np.array(dataLs[1])[inclEB1]])
toplotEB2 = np.array([AVGJK, AVGJKMIN, AVGJKMAX, COUNTS])

# 12. PLOT DATA --------------------------------------------------------------------
figObj = makeplot(toplotEB1, toplotEB2)
figObj.savefig(FOLDER_ROOT + FOLDER_OUT + 'JK_JF12_' + GRAV.lower() +  \
              '.pdf', dpi=600)
