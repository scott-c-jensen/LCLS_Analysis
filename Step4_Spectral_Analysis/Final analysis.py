"""This function will further bin the data, combining spectra and plotting the data"""

#!/usr/bin/python
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter as sf
import scottLib as sl
import scipy.signal as sci
from scipy.interpolate import interp1d
import glob
import re
reload(sl)

mpl.rcParams['axes.formatter.useoffset'] = False

"""
I could make a 2d plot with intesity for IO vs Pulse intensity
"""
#-------------------------------------------------------------------------------
#Functions for Solids

def shiftData(xIn, yIn, shift):
    """
    Fits the xIn and yIn data with a cubic spline with the data shifted by "shift" in the same units as xIn. 
    It then resamples the shifted fit using the original unshifted x values
    xIn, yIn: List or arrays for the x and y data to be shifted
    shift: Amount to shift the x-axis of the data 
    """
    spline=interp1d(xIn+shift,yIn, kind='cubic')
    outSpec = np.zeros(np.shape(xIn))
    outSpec[10:-10] = spline(xIn[10:-10])
    return(outSpec)

def sortInputs(folder, sortNum=-2):
    """
    Lists all the xes files in order of run number
    folder: string path to the folder of interest
    sortnum: The 
    """
    fileList = glob.glob(folder+'*.h5')
    fileListOut = sl.sortFiles(fileList,sortNum)
    return(fileListOut)
   
def summedSpec(inputFile):
    """Sums data in a single file
    inputFile: file path str"""
    with h5py.File(inputFile,'r') as fIn:
        data = np.array(fIn['BinnedMatrix']).sum(axis=0).sum(axis=0)
        xData = np.array(fIn['SpectraBins'])
        dataNew = shiftData( xData, data, 0.6)
    return(xData, dataNew)

def loadSolids(inputFolder):
    """Loads all the data from solid samples tested
    inputFolder: path to folder"""
    useFiles = sortInputs(inputFolder)
    xData = summedSpec(useFiles[0])[0]
    yOut = [summedSpec(fileIn)[1] for fileIn in useFiles]
    summedData = [yOut[0],np.add(yOut[1],yOut[2]),yOut[3],np.add(yOut[4],yOut[5])]
    return(xData,summedData)
    
def makeFM(xData, yList):
    """Finds the first moment of the data given the range
    xData = ev data 
    yList = List of y arrays"""
    fms = [sl.firstMoment(xData, yData, *fmRange) for yData in yList]
    return(fms)

def plotSolids(xData, yList):
    """Plots solid data
    xData: array x values
    yList: list of arrays of y values"""
    plt.figure('Oxides')
    plt.clf()
    for i in range(len(yList)):
        plt.plot(xData,sl.normPeak(yList[i]), label = solidNames[i])
    plt.legend()
    plt.show()
    return()

#-------------------------------------------------------------------------
#Load MnCl2
def sortFiles(fileList):
    """Sorts files by the second to last number in file path
    fileList: list of file paths"""
    p = re.compile('\d+')
    runNum = [p.findall(tempFile)[-2] for tempFile in fileList]
    sortList = [fileList[i] for i in np.argsort(runNum)]
    return(sortList, np.sort(runNum))
   
def loadMatrix(inputFile, shift =0):
    """Loads in a matrix of data from inputFile and shifts it
    inputFile: file path str
    shift: amount to shift data (based on standards)"""
    with h5py.File(inputFile, 'r') as fIn:
        data = np.asarray(fIn['BinnedMatrix'])
        ioRange = np.asarray(fIn['IOBins'] )
        eRange = np.asarray(fIn['SpectraBins'])+shift
        if 'PulseDurationBins' in fIn:
            pulseRange = np.asarray(fIn['PulseDurationBins'])
        else:
            print('No xtCav Data')
            pulseRange = None
    return(data, eRange, ioRange, pulseRange)

def loadAllMatrix(inputFolder, shiftList):
    """loads all data in folder
    inputFolder: folder path str containing files to load
    shiftList: list of shift amounts for each file"""
    dataList, eList, ioList, pulseList=[],[],[],[]
    inputFileList = glob.glob(inputFolder+'*.h5')
    sortedFiles, runs = sortFiles(inputFileList)
    for i, fileIn in enumerate(sortedFiles):
        dArray, e, io, pul = loadMatrix(fileIn, shiftList[i])
        dataList.append(dArray)
        eList.append(e)
        ioList.append(io)
        pulseList.append(pul)
    return(dataList, eList, ioList, pulseList, runs)
    
def get2016Data():
    """Loads in parameters for 2016 and loads data"""
    inFolder16 = 'C:/Scott/LCLS/Data - 10-18-2017/MnCl2 2016/'
    #'C:/Scott/LCLS/Data - 10-18-2017/MnCl2 2015/Run_0006_NOHITFINDING.h5']
    shiftAll = 0.6
    shift305 = 1.1
    shift310 = 0.6
    shiftList16 = [shiftAll]*7 + [shift305] + [shiftAll]*2 + [shift310] + [shiftAll]
    dataList, eList, ioList, pulseList, runs = loadAllMatrix(inFolder16, shiftList16)
    return(dataList, eList, ioList, pulseList, runs)

def get2015Data():
    """Loads in parameters for 2016 and loads data"""
    inFolder15 = 'C:/Scott/LCLS/Data - 10-18-2017/MnCl2 2015/'
    #'C:/Scott/LCLS/Data - 10-18-2017/MnCl2 2015/Run_0006_NOHITFINDING.h5']
    shift1 = 1.3
    shift2 = 2.525
    shift3 = 2.1
    shiftList15 = [shift1]*5 + [shift2] + [shift3]*3
    dataList, eList, ioList, pulseList, runs = loadAllMatrix(inFolder15, shiftList15)
    return(dataList, eList, ioList, pulseList,  runs)

#-----------------------------------------------------------------------------
#Data splitting, bad point removal, FM calculation..
def runNumToIdx(runNum, runNumList):
    '''
    Converts the run num to the index
    '''
    runNumList = [int(xx) for xx in runNumList]
    idx=int(np.where(np.asarray(runNumList) == runNum)[0])
    return(idx)

def manualRm(dataX, dataY, year, runNum):
    '''
    Removes bad points from noisey/dead pixels that couldn't be interpolated previous
    dataX, dataY: lists of x and y data
    year: int of experiment
    runNum: int run identification
    '''
    if year==2016:
        if runNum in range(113,130)+[308,311]:
            badPoints = [x / 10. for x in range(64830, 64838)]
        elif runNum == 305:
            badPoints = [x / 10. for x in range(64834, 64841)]
        elif runNum == 310:
            badPoints = [x / 10. for x in range(64830, 64841)]
        elif runNum == 306:
            badPoints = [x / 10. for x in range(64829, 64843)]
        else:
            raise Exception(str(runNum)+' is not a valid runNumber')
        
    elif year == 2015:
        if runNum in [6,8,9,13]:
            badPoints = [x / 10. for x in range(64901, 64906)]
        elif runNum == 10:
            bP1 = [x / 10. for x in range(64901, 64906)]
            badPoints = [6487.5, 6487.6]+bP1
        elif runNum == 60:
            badPoints = [6492.2,6492.3]+[x / 10. for x in range(64985, 64990)]
        elif runNum in [155,156,157,158]:
            badPoints = [6493.6,6493.7,6493.8]+[x / 10. for x in range(64998, 65003)]
        else:
            raise Exception(str(runNum)+' is not a valid runNumber')

    else:
        raise Exception(str(year)+' is not a valid year. Only 2015 and 16')

    xNew = np.round(dataX, decimals=1)
    yNew = dataY
    for xVal in badPoints:
        idx = np.where(xNew==xVal)[0]
        if len(idx)==1:
            xNew=np.delete(xNew, idx)
            yNew=np.delete(yNew, idx)
        else:
            print"bad point"
    xOut, yOut = sl.rebin(xNew, yNew,xStart, xEnd, xStep)
    return(xOut, yOut)
    
def getMean(values, intensityVector):
    '''
    This function will give the average IO or pulse range for the matrix
    values: array of values
    intensityVector: vector for each point in the spectrum
    '''
    if len(values)!=len(intensityVector):
        raise Exception ('Lenghts for the two vectors are not equal, no mean calculated')
    binVals = np.array(values)
    inten = np.array(intensityVector)
    normVector = (binVals*inten).sum()/(inten.sum()+0.0000001)   
    return(normVector)

def findMeans(matrix, io, ioRange, pulse=None, pulseRange=None):
    """Finds the means for io and pulse if it is given for each of the final binned data
    matrix: combination of all data
    io: x-ray intensity
    ioRange: range of values to bin
    pulse: x-ray pulse duration
    pulseRange: range to integrate into one spectrum"""
    ioDigi = np.digitize(ioRange, io)
    if type(pulse) == type(None):
        ioMean = getMean(io[ioDigi[0]:ioDigi[1]],matrix[ioDigi[0]:ioDigi[1],:].sum(axis=1))
        return (ioMean, None)
    else:
        pulseDigi = np.digitize(pulseRange, pulse)
        tempM = matrix[ioDigi[0]:ioDigi[1],pulseDigi[0]:pulseDigi[1],:]
        ioMean = getMean(io[ioDigi[0]:ioDigi[1]],tempM.sum(axis=2).sum(axis=1))
        pulseMean = getMean(pulse[pulseDigi[0]:pulseDigi[1]],tempM.sum(axis=2).sum(axis=0))
        return(ioMean, pulseMean)
    
def collapseMatrix(xData, data, io, pulse, ioRange, pulseRange, year, runNum):
    """
    This will collapse data for one set given the ioRange and pulseRange desired
    Also runs the manual Rm which is important for rebinning
    It only gives back the collapsed Data
    xData: x-ray energy as list
    data: array of all data
    io: pulse intensity list
    pulse: pulse duration list
    ioRange, pulseRange: ranges to integrate the respective variables
    year: int year for date first captured data
    runNum: unique identifier for the run
    """
    ioDigi = np.digitize(ioRange, io)
    if pulse is None:
        tempData = data[ioDigi[0]:ioDigi[1],:].sum(axis=0)
    else:
        pulseDigi = np.digitize(pulseRange, pulse)
        tempData = data[ioDigi[0]:ioDigi[1], pulseDigi[0]:pulseDigi[1],:].sum(axis=0).sum(axis=0)
    xData, outData = manualRm(xData, tempData, year, runNum)
    return(xData, sl.normPeak(outData))

def getSingleFM(emissionE, data, io, pulse, pulseRange, ioRange, runNum, year):
    """
    Finds the FM and the average IO or pulse duration for each dataSet give the IO and pulseRange desired
    emissionE: energy of emission
    data: array of data
    io: ion chamber intensities
    pulse: duration of x-ray pulse
    pulseRange, ioRange: Range of values to integrate
    runNum: unique identifier for a run
    year: int year of experiment
    """
    #Get mean values for pulse and IO
    meanIO, meanPulse = findMeans(data, io, ioRange, pulse, pulseRange)
    #Getting FM (remove bad points)
    xData, yData = collapseMatrix(emissionE, data, io, pulse, ioRange, pulseRange, year, runNum)
    fm =sl.firstMoment(xData, yData,*fmRange)
    return(xData, yData, fm, meanIO, meanPulse)


def getFM(useRuns,years,ioRange,pulseRange):
    '''
    Finds all first moments for the data considered
    Combines data based on the ranges chosen
    useRuns: list of run numbers to use
    years: years of same length as run numbers (one for each run)
    ioRange: range of values to bin by in x-ray intensity
    pulseRange: range of values to bin by in x-ray pulse duration
    '''
    xDataList,yDataList,meanIOList,meanPulseList,fmList = [],[],[],[],[]
    for i, useRun in enumerate(useRuns):
        if years[i] == 2016:
            idx = runNumToIdx(useRun, runs16)
            xData, yData, fm, meanIO, meanPulse = getSingleFM(e16[idx], data16[idx], io16[idx], pulse16[idx], pulseRange[i], ioRange[i], useRun, 2016)
        elif years[i] == 2015:
            idx = runNumToIdx(useRun, runs15)
            xData, yData, fm, meanIO, meanPulse = getSingleFM(e15[idx], data15[idx], io15[idx], pulse15[idx], pulseRange[i], ioRange[i], useRun, 2015)
        xDataList.append(xData)
        yDataList.append(yData)
        meanIOList.append(meanIO)
        meanPulseList.append(meanPulse)
        fmList.append(fm)

    return (xDataList,yDataList,fmList,meanIOList,meanPulseList)
    
#-------------------------------------------------------------------------------
#Final Plots and Normalization
def getPulseParams(runs, years):
    """Get pulse parameters based on the run number and year
    runs: list of runs to get parameters for
    years: list of years corresponding to runs (must be same length)"""
    photonEArray, attenArray = np.empty(len(runs)),np.empty(len(runs))
    for i, run in enumerate(runs):
        if years[i] == 2016:
            photonEnergy = 8690
            runs = range(113,133)+range(305,314)
            attens = [1]*5+[0.0056]*3+[0.507]*2+[0.017]*3+[0.2]*3+[0.053]*2+[0.0056]*2+[1]*3+[0.507]*2+[1]+[0.017]*3
        elif years[i] ==2015:
            photonEnergy = 6900
            runs = range(6,14)+[60]+range(155,161)
            attens = [1]*2+[0.045]+[0.019]+[0.012]*3+[0.267]+[0.019]+[0.267]*2+[1]+[0.012]*3
        else:
            raise Exception("Must use a valid year")
        idx = runNumToIdx(run, runs)
        photonEArray[i]=photonEnergy
        attenArray[i]=attens[idx]
    return(photonEArray, attenArray)

def photonsPerPulse(ioAverage, runs, years):
    """Gets the total number of photons in a single x-ray pulse
    ioAverage: the average value read out by the gas chamber
    runs: is needed to estimate the attenuation factor used for that data set
    Years: is the years each run was taken (this is used for photon Energy)"""
    photonE, atten = getPulseParams(runs,years)
    mJ_per_eV = 1.6021773e-16 #This is for mJ not J
    photonArray = np.array(ioAverage)*atten/(photonE*mJ_per_eV)
    return(photonArray)

def getPhotonParams2(runs, years):
    """
    Returns the parameters for each run to make normalizations
    In each list of atoms in the sample the Z is high to low (so Mn is the first in each array)
    runs: is needed to estimate the attenuation factor used for that data set
    Years: is the years each run was taken (this is used for photon Energy:
    """
    photonEnergyA, crossSectionA, molarA, mwA = np.zeros(len(runs)), np.zeros((len(runs))), np.zeros((len(runs),5)), np.zeros((len(runs),5))
    for i, year in enumerate(years):
        if year == 2015:
            photonEnergy =  6900
            if runs[i] < 15: #Use H2O as the solvent
                crossSection = np.array([3.17880*10**-12] )
            else: #Use ethanol as solvent
                crossSection = np.array([3.17880*10**-12] )
        elif year == 2016:
            photonEnergy = 8690
            crossSection = np.array([1.75774*10**-12] )

        photonEnergyA[i], crossSectionA[i]= photonEnergy, crossSection
    return(photonEnergyA, crossSectionA)

def getMnNorm2(ioAverage,useRuns, years):
    """
    This normalized by number of photons absorbed per Mn atom
        incident photons * probability of absorption / # MnAtoms         (use this)
    or for absorption probability
        cross section * density * pathlength
        cross section * ((molarity*volume*molarMass)/Volume)*pathlength
        cross section * molarity*molarMass*pathlength                    (use this)
    The number of Mn atoms are
        molarity * avegadros number * volume probed                      (use this)  
    So the final equation is
        incident photons * cross section * molarMass / (avegadros number * beam area)
    ioAverage: average intensity for binned specra
    useRuns, years: lists to uniquely identify which run and year the data was captured
    """
    
    radius = 1 #in um
    area = np.pi*radius**2
    photonE, crossS = getPhotonParams2(useRuns, years)
    print(np.shape(crossS))
    norm = photonsPerPulse(ioAverage, useRuns, years)*crossS/area #Use only the cross sections for Mn
    print norm
    normStr = 'Photons absorbed per Mn atom'
    return(norm, normStr)
################################################################################
################################################################################
################################################################################
    
def getDoseNorm(ioAverage,useRuns, years):
    """
    Here I find the dose delivered for each incident photon so that this normalization factor will 
    To get Dose (energy deposited per area) I use the following steps:
        incident photons * probability of absorption * photon energy / volume 
    So I can simplify the probability of absorption as:
        incident photons *cross section * density * pathlength
        incident photons *cross section * ((molar*volume*molarMass)/Volume)*pathlength
    So I can rewrite the whole thing as:
        incident photons *cross section * ((molar*molarMass)*pathlength * photon energy / volume 
        incident photons *cross section * molar * molarMass * (photon energy in J)/ area
    ioAverage: average intensity for binned specra
    useRuns, years: lists to uniquely identify which run and year the data was captured
    """
    radius = 0.0001 #in cm
    area = np.pi*radius**2
    eConversion = 1.6021*10**-19
    photonE, crossS, molarity, mws = getPhotonParams(useRuns, years)
    norm = crossS*molarity*mws*(photonE*eConversion)/area*photonsPerPulse(ioAverage, useRuns, years)
    normStr = 'Dose delivered (J/cm^3)'
    return(norm, normStr)

def getPhotosAbsorbed(ioAverage, useRuns, years):
    """
    This normalization is for the number of photoelectrons
    To get Dose (energy deposited per area) I use the following steps:
        incident photons * probability of absorption 
        incident photons *cross section * ((molar*volume*molarMass)/Volume)*pathlength
        incident photons *cross section * molar * molarMass * pathlength
    ioAverage: average intensity for binned specra
    useRuns, years: lists to uniquely identify which run and year the data was captured
    """
    pathlength = 0.05  # 50 um diameter jet
    photonE, crossS, molarity, mws = getPhotonParams(useRuns, years)
    norm = crossS*molarity*mws*pathlength*photonsPerPulse(ioAverage, useRuns, years)
    normStr = 'Number of photons absorbed per pulse'
    return(norm, normStr)
    
def makeDataSet(useRuns, years):
    """
    This function separates the list into the appropriate groups
    This also lists the reference run at the beginning
    useRuns, years: lists to uniquely identify which run and year the data was captured
    """
    refRuns = [10,158,118,311]
    setNum = [[],[],[],[]]
    for i, year in enumerate(years):
        if year == 2015:
            if useRuns[i] <= 15:
                if useRuns[i]== refRuns[0]:
                    setNum[0].insert(0,i)
                else:
                    setNum[0].append(i)
            else:
                if useRuns[i]== refRuns[1]:
                    setNum[1].insert(0,i)
                else:
                    setNum[1].append(i)
        elif year == 2016:
            if useRuns[i] <= 150:
                if useRuns[i]== refRuns[2]:
                    setNum[2].insert(0,i)
                else:
                    setNum[2].append(i)
            else:
                if useRuns[i]== refRuns[3]:
                    setNum[3].insert(0,i)
                else:
                    setNum[3].append(i)
    labels = ['40 fs pulse (6900 eV)','20 fs pulse (6900 eV)','16 fs pulse (8690 eV)','34 fs pulse (8690 eV)']
    return(setNum, labels)
    
def getDataSets(useRuns, years, xIn, yIn):
    """gets all data sets for runs/years and the change in spectra from that and the standard
    useRuns: which runs to use
    years: which year data taken
    xIn, yIn: reference spectrum data"""
    sets, labels = makeDataSet(useRuns, years)
    xList = [np.array(xIn)[np.array(sets[idx])] for idx in range(len(sets))]
    yList = [np.array(yIn)[np.array(sets[idx])] for idx in range(len(sets))]
    yChangeList = [yData-yData[0] for yData in yList]
    return(xList, yList, yChangeList, labels)
    
def plotFMTrend(normFn, runs, years, fmList,ioAverage, solidFMChange):
    """
    Plots the FM as a function of the normalization spectrum
    normFn: function to normalize spectra
    runs, years: identifies for which dataset
    fmList: list of first moments of each input run
    ioAverage: average x-ray intensity for each dataset
    solidFMChange: change in the FM for solid standards
    """
    xData, xLabel = normFn(ioAverage, runs, years)
    xDataSets, _, yChangeList, labels = getDataSets(runs, years, xData, fmList)
    
    plt.figure('First Moment Change')
    plt.clf()
        
    for i in [2,3,1,0]:
        #plt.plot(photonsListOfLists[i], fmListOfLists[i], linestyle = '', marker = 'o', label=labels[i], markersize=11)
        plt.plot(xDataSets[i], yChangeList[i], linestyle = '-', marker = 'o', label=labels[i], markersize=8)
        
        #plt.errorbar(fmList[i], fitFm(fmList[i]), xerr=seFmList[i], color='black')
    #oxColor = ['r','b','g','black']
    oxLines = ['-.','--','-']
    for i in range(len(solidFMChange)):
        #plt.axhline(y=solidFMChange[i], color=oxColor[i], linestyle='-', label=oxLabel[i])
        plt.axhline(y=solidFMChange[i], color='black', linestyle=oxLines[i], label=solidNames[i])
    plt.xlabel(xLabel, fontsize='20')
    plt.ylabel('Change in First Moment (eV)', fontsize= '20')
    plt.legend(bbox_to_anchor=(.5, 1), loc=2, borderaxespad=0., frameon=False,  fontsize = '16')
    plt.tight_layout()
    #plt.xscale('log')
    ax = plt.gca()
    #ax.get_xaxis().get_major_formatter().set_useOffset(False)
    #x_formatter = sf(useOffset=False)
    #ax.xaxis.set_major_formatter(x_formatter)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.tick_params( labelsize=16)
    plt.tight_layout()
    #plt.figtext(x=.1, y=.05, fontdict=None, s=description,  fontsize = '16') #add s="figure 1 etc
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15)
    plt.rcParams.update({'font.size': 13.5})
    plt.show() 
    return()
    
def plotSpectraSets(xIn, yIn, photonsPPulse):
    """Plot Spectral sets together
    xIn, yIn: data for all sets
    photonsPPulse: number of photons in a single x-ray pulse"""
    setNames = ['2015 - 60fs', '2015 - 20fs', '2016 - 16fs', '2016 - 34fs']
    sets = [[0,6],[6,10],[10,16],[16,20]]
    for j in range(len(setNames)):
        plt.figure(setNames[j])
        plt.clf()
        for i in range(sets[j][1]-1,sets[j][0]-1,-1):
            sl.makePlot(xIn[i], yIn[i], figure= setNames[j], label='{:.1e}'.format(photonsPPulse[i])+' Photons/pulse')

        plt.xlim(6473,6497)
    return() 
 
def runFM(normFn, runs, years, ioRange, pulseRange):
    """Runs first moment analysis on the files in runs/year and normalizes them
    normFn: function to normalize the spectra
    runs, years: identifiers for the run/year the data was taken
    ioRange: range of values to bin by in x-ray intensity
    pulseRange: range of values to bin by in x-ray pulse duration"""
    xDataList, yDataList, fmList, meanIOList, meanPulseList = getFM(useRuns, years,ioRange,pulseRange)
    plotFMTrend(normFn, runs, years, fmList, meanIOList, solidFMChange)
    
    phNum = photonsPerPulse(meanIOList, useRuns, years)
    plotSpectraSets(xDataList, yDataList, phNum)
    return(xDataList, yDataList, fmList, meanIOList, meanPulseList)    


################################################################################
    
#-------------------------------------------------------------------------------
#Load Solid Oxides
xStart, xEnd, xStep = [6460, 6512, 0.1]
fmRange =[6485., 6495.]
plotRange = [6472, 6498]
solidFolder = "C:\\Users\\Pushkar\\Desktop\\Beamtime_Data\\LCLS\\12-11-2017\\processed\\Solids\\"
solidNames = ['MnO - MnO'] +[ 'Mn2O3 - MnO']+ ['MnO2 - MnO']+['KMnO4']
solidSpin = [2.5,2,1.5,0]
solidFM = makeFM(*loadSolids(solidFolder))
plotSolids(*loadSolids(solidFolder))
solidFMChange = [ x-solidFM[0] for x in solidFM[:-1]]
#-------------------------------------------------------------------------------
#Note that the range for IO is 0.1-5 while pulse is 0-50

#Load MnCl2 Data
data16, e16, io16, pulse16, runs16 = get2016Data()
data15, e15, io15, pulse15, runs15 = get2015Data()

#Data to use and IO/pulse duration to use
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#I dont sort the values for plotting so go from least to highest intensity for plotting for runs

useRuns2015 = [10,9,8,13,6,6,158,155,157,157] #60,
useRuns2016 = [118,123,129,126,121,117,311,308,306,306] #305,310 #113  ###Reasonably exclude 129 or is something off?
useRuns = useRuns2015+useRuns2016
years = [2015]*len(useRuns2015)+[2016]*len(useRuns2016)

#Manually set the IO and pulse ranges used
rangeIO15 = [[0.1,15]]*4+[[0.1,2.5]]+ [[2.5,5]]+[[0.1,5]]*2+[[0.1,1.5]]+ [[2,4.5]]
rangeIO16 = [[0.1,5]]*(len(useRuns2016)-2)+[[0.1,2]]+ [[2,5]] 
rangeIO =  rangeIO15+rangeIO16
rangePulse15 = [[None]]*len(useRuns2015)
rangePulse16 = [[10,20]]*6+[[30,40]]*2+[[30,36]]+[[31,40]]
#rangePulse16 = [[1,50]]*10
rangePulse = rangePulse15+rangePulse16

xData, yData, fm, meanIOList, meanPulseList = runFM(getMnNorm2, useRuns, years, rangeIO, rangePulse)
