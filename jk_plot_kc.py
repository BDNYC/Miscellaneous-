''' 
This procedure reads and calculates J-K colors for all objects Paper XII, and then plot them as a function of spectral type. It also calculates the average for each spectral type and plots those too. It can also plot excluded peculiar objects (blue, dusty, pec) (the code for this is commented out inside makeplot()).

'''

def makeplot(dataSP):
    
    import numpy
    import matplotlib.pyplot as plt
    import pdb
    
    try:
        dataSP.shape
    except AttributeError:
        return
    
    # 1) Initialize variables
    BLACK = '#000000'
    GRAY  = '#999999'
    RED   = '#FF0000'
    WHITE = '#FFFFFF'
    X_LABEL = 'Spectral Type'
    Y_LABEL = 'J-K$_s$ (2MASS)'
    
    # 2) Initialize Figure
    plt.close()
    plt.rc('font', size=12)
    fig = plt.figure(1, figsize=(6,5))
    plt.clf()
    
    # 3) Create spectral bins
    bins = []
    labels = []
    for tp in range(11, 20):
        idxHi = np.where(dataSP[0,:] < tp)
        idxLo = np.where(dataSP[0,:] >= tp - 1)
        if len(idxHi[0]) > 0 and len(idxLo[0]) > 0:
            idxtp = np.intersect1d_nu(idxHi[0], idxLo[0])
        else:
            idxtp = []
        if len(idxtp) != 0:
            bins.append(dataSP[1,idxtp])
            labels.append('L' + str(tp - 11))
        # Rememeber the smallest J-K value of L5 to annotate label there later
        if tp == 16:
            lblymin = np.min(dataSP[1,idxtp])
        
    # 4) Initialize subplot
    ax = fig.add_subplot(1,1,1, position=[0.12,0.11,0.85,0.84], zorder=10)
                                        # [left,bottom,width,height]
    ax.set_autoscale_on(False)
    
    # 5) Plot the boxplot bins
    xs = np.floor(dataSP[0,:] - 10)
    ys = dataSP[1,:]
    sc = plt.scatter(xs, ys, marker='o', s=20, edgecolor=GRAY, facecolor=GRAY, \
                     label='individual L dwarfs')
    
    # 6) Plot the averages
    # Generate four point with same J-K and spectral type right around the center
    # to re-create a bar type of marker
    XS2 = np.array([-0.09,-0.04,0.,0.04,0.09])
    xs2 = XS2.copy()
    for shift in range(1,9):
        tmpx = XS2 + shift
        xs2 = np.concatenate((xs2,tmpx))
    avgs = []
    for bn in bins:
        avg = np.average(bn)
        avgs.append(avg)
        avgs.append(avg)
        avgs.append(avg)
        avgs.append(avg)
        avgs.append(avg)
    ax.scatter(xs2, avgs, marker='s', s=5, edgecolor='k', facecolor='k', \
               label='average')
    
    # # 7) Plot the scatter peculiar dots
    # xs = np.floor(dataSP[0,:] - 10)
    # ys = dataSP[1,:]
    # plt.scatter(xs, ys, c=GRAY, s=20, marker='x', alpha=0.5, zorder=1)
    
    # 8) Format plot (make lots of things disappear)
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
    ax.set_ylim(0.9,2.3)
    ax.set_xlim(-0.5, 8.5)
    
    # Add tick labels
    plt.xticks(range(0,9), labels, size=14)
    
    # Add horizontal grid
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.8)
    ax.set_axisbelow(True)
    
    # Annotate average label
    loc = (6. + XS2[4], avgs[6 * len(XS2)])
    loctext = (55, -50)
    linetype = dict(arrowstyle='-', shrinkB=4, shrinkA=2, color=BLACK, relpos=(0,0))
    plt.annotate('averages', xy=loc, xytext=loctext, xycoords='data', \
                 textcoords='offset points', fontsize=14, color=BLACK, \
                 ha='center', arrowprops=linetype)
    
    # Annotate individual object label
    loc = (5., lblymin)
    loctext = (-10, -60)
    linetype = dict(arrowstyle='-', shrinkB=6, shrinkA=4, color=GRAY)
    plt.annotate('individual L dwarfs', xy=loc, xytext=loctext, xycoords='data', \
                 textcoords='offset points', fontsize=14, color=GRAY, \
                 ha='center', arrowprops=linetype)
    
    return fig


# 1. LOAD RELEVANT MODULES ---------------------------------------------------------
import asciidata
import astrotools as at
import matplotlib.pyplot as plt
import numpy as np
import sys
import pdb

# 2. SET UP VARIABLES --------------------------------------------------------------
# General variables
GRAV = 'f'
FOLDER_ROOT = '/Users/alejo/KCData/'  # Location of NIR and OPT folders
FOLDER_OUT  = 'Output/special/'
dataRaw = ''

# For TXT objects file (updatable here directly)
FILE_IN = 'NIR_Spex_Prism_with_optical_12Aug15.txt' # ASCII file w/ all objects
FILE_EXCL = 'Exclude_Objs.txt'   # ASCII file w/ objects to exclude
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

# Determine which targets to include in the data set
inclSP = np.array([False] * numRows)
if GRAV.upper() == 'G':
    incIdx = np.where(np.array(gravObjs) == 'g')
    if len(incIdx[0]) != 0:
        inclSP[incIdx] = True
elif GRAV.upper() == 'B':
    incIdx = np.where(np.array(gravObjs) == 'b')
    if len(incIdx[0]) != 0:
        inclSP[incIdx] = True
elif GRAV.upper() == 'F':
    incIdx = np.where(np.array(gravObjs) == 'f')
    if len(incIdx[0]) != 0:
        inclSP[incIdx] = True

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
            inclSP[tmpExclIdx] = False

# 8. EXCLUDE BINARIES --------------------------------------------------------------
for binIdx, bin in enumerate(dataLs[6]):
    if bin == 'YES':
        inclSP[binIdx] = False

# 9. EXCLUDE SPECIAL OBJECTS FROM THE DATA SET ------------------------------------
for spIdx, sp in enumerate(dataLs[0]):
    if dataLs[4][spIdx] == 'Yes' or dataLs[5][spIdx] == 'Yes' or \
                                                dataLs[7][spIdx] == 'Yes':
        inclSP[spIdx] = False

# 10. PREPARE DATA TO PLOT ----------------------------------------------------------
toplotSP = np.array([np.array(dataLs[3])[inclSP], np.array(dataLs[1])[inclSP]])

# 11. PLOT DATA --------------------------------------------------------------------
figObj = makeplot(toplotSP)
figObj.savefig(FOLDER_ROOT + FOLDER_OUT + 'JK_XII_' + GRAV.lower() +  \
              '.pdf', dpi=600)
