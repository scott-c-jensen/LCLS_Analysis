"""
Author: Scott C. Jensen 2016
This program just sums spectra within the same bin region (ie bins by pulse duration, incident pulse energy, and emission energy)

This will also rebin the spectra so they all have comparable x-axis values


"""
import numpy as np
import h5py
import glob
import os
import argparse

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def sortFiles(files, whichNumToIndex=-1):
    """
    This function sorts files by a number in the file name. 
    Typically used to sort files that are increasing but do not have leading zeros
    files: list of files which have numbers in their strings to sort by
    whichNumToIndex: the index of the number sorted by, defaults to the last number (-1)
    """
    import re
    
    p = re.compile('\d+')
    unsort = [float(p.findall(tempFile)[whichNumToIndex]) for tempFile in files]
    sortList = [files[i] for i in np.argsort(unsort)]
    return(sortList)

def getVals(h5File, xtCav=True):
    """Loads up the H5 file with all the data
    h5file: path to file as string
    xtCav: Bool if it was recorded"""
    #Read 
    with h5py.File(h5File, 'r') as fileIn:
        
        #This checks if I0 exists or if 
        if 'I0' not in fileIn.keys() :
            return( None, None,None,None,None)
        if len(fileIn['I0'].value) ==0:
            return( None, None,None,None,None)

        #Get Spectrum
        spec = fileIn['spectra'].value
        specX = fileIn['energy'].value

        #Bin every 0.1 for IO
        io = fileIn['I0'].value
        if xtCav ==True:
            #Bin for every 1 from 0 to 45 for pulse duration
            pulseDur = fileIn['XTCAV_Pulse_Duration'].value
            pulseAccuracy = fileIn['XTCAV_Reconstruction_Agreement']
        else:
            pulseDur = None
            pulseAccuracy = None

    return(spec, specX, io, pulseDur, pulseAccuracy)
    
def getFilesFromRun(inDir, runNum, curClass):
    """Get the files for the run/class
    inDir: directory to get files
    runNum: int umber assigned to the run
    curClass: 4 assignations given by crystallographers (laser on/off and crystal diffraction present/not present)
    """
    #Get the files for the run/classs
    files = glob.glob('{Dir}*r{run:04d}*{curClass}*'.format(run=runNum, Dir=inDir, curClass=curClass))
    return files

def rebin (xIn, yIn, Start, Stop, Step, interp = 11):
    """Rebins xIn and yIn and interpolates to a new grid
    xIn, yIn: x and y data as array
    start, stop, step: values to create a new grid to interpolate to
    interp: number of surrounding points used to interpolate to any given point"""
    if interp%2==0:
        raise Exception('Interpolation number not odd')
        
    halfStep = Step/float(interp)*(interp-1)/2.
    xInterp =(Start-halfStep, Stop+halfStep, Step/interp)
    yInterp = np.interp(xInterp, xIn, yIn)
    yOut = np.reshape(yInterp,(-1,interp)).sum(axis=-1)
    return(np.arange(Start,Stop,Step),yOut)
    
def stepNum(start, stop, step):
    """Finds an approximate even number of steps between start and stop
    start: begining value for range
    stop: ending value for range
    step: step size to take
    
    returns number of steps and the range"""
    stepNum = np.round((stop-start)/step)+1
    newStop = start+step*(stepNum-1)
    stepRange = np.linspace(start, newStop, stepNum)
    return (stepNum, stepRange)
    
def makePath(mkDir):
    if not os.path.exists(mkDir):
        os.makedirs(mkDir)
    return()
    
def xtCavStat(inDir, runNum, curClass):
    """Gets stats for XTCav for each run and class
    inDir: input directory
    runNum: int assigned to the run
    curClass: current class"""
    inFiles = getFilesFromRun(inDir, runNum, curClass)
    with h5py.File(inFiles[0]) as fTest:
        if 'XTCAV_Pulse_Duration' in fTest:
            xtCavBool = True
        else:
            xtCavBool = False
    return xtCavBool
    
def makeList(arg):
    if isinstance(arg, list):
        return arg
    else:
        return([arg])
        
def makeClass(useClasses):
    """Which class the data were binned into (not the python class but group of data)
    There are four classes, laser on/off and crystal diffraction present/not present
    NOHITFINDING is for non-crystal data which didn't require a laser (ie one class)
    """
    if useClasses == True:
        return ['class' + str(i) for i in range(0,4)]
    else:
        return ['NOHITFINDING']

def newDigi(listIn, rangeIn, binMax):
    """
    Rebins the list
    listIn: list of data
    rangeIn: range of values to use
    binMax: max value to consider a bin
    """
    digiTemp = np.digitize(listIn, rangeIn)
    digiOut = [binCur for binCur in digiTemp if binCur<=binMax]
    return(digiOut)

def get_pulse_info(xtCavBool, inFile, pulseRange, pulseBinMax):
    """Returns the pulse characteristics
    xtCavBool: Bool indicating if the XTCav was recorded
    inFile: H5 file where the data was stored previously"""
    specList, specXList, ioList, pulseDurList, pulseAccuracyList = getVals(inFile, xtCav=xtCavBool)
    if specXList is None:
        print('This file has no data')
        return None #skip this iteration since there is no data
    elif xtCavBool is True:
        pulseDigi = newDigi(pulseDurList, pulseRange, pulseBinMax)
    else:
        pulseDigi = [None]*6
    return (pulseDigi, specList, specXList, ioList, pulseDurList, pulseAccuracyList)
    
parser = argparse.ArgumentParser(description='Bin spectra for different pulse durations and I0')
parser.add_argument('-s', '--sortA', help='Sort by accuracy (must be equal or higher to this value)', type=float, default = 0.5)
parser.add_argument('-i', '--inputDir', help='Input directory', type=str)
parser.add_argument('-o', '--outputDir', help='Output directory', type=str)
parser.add_argument('-c', '--useClass', help='If classes (ie hits (1 and 3) misses (0 and 2) and laser on (0 and 1) and laser off (2 and 3)', type=str)
parser.add_argument('-r', '--runList', nargs='+', help='Run numbers to make sum folder from (Use as last input)', type=int)
parser.add_argument('-x', '--xtCav', help='use xtCav', type=str)

args = parser.parse_args()
    
def main(args):
    """Takes all the individual spectra and combines them based on bins
    Bin the input x-ray intensity, duration and photon energy
    This will also rebin the spectra so they all have comparable x-axis values"""

    inputDir = args.inputDir
    outputDir = args.outputDir
    useClasses = str2bool(args.useClass)
    xtCavBool = str2bool(args.xtCav)

    #Add spectra into binned format of Binned Matrix
    runList = makeList(args.runList)
    notes = 'Souce: ' + inputDir+ '\nRuns: '+str(args.runList) + ' \nAccepted Accuracy for xtCav agreement: '+str(args.sortA)

    #Spectrum Rebin Factor
    eStep = 0.1
    eStart = 6460
    eStop = 6512

    #Bin IO
    ioStep = 0.1
    ioStart = 0.1
    ioStop = 5
    ioBinMax = ioStop/ioStep

    #Bin Duration
    pulseStep = 1
    pulseStart = 0
    pulseStop = 50
    pulseBinMax =pulseStop/pulseStep

    #Num of Bins for each variable
    ioBinNum, ioRange = stepNum(ioStart, ioStop, ioStep)
    pulseBinNum, pulseRange = stepNum(pulseStart, pulseStop, pulseStep)
    specBinNum, eRange = stepNum(eStart, eStop, eStep)

    #Create the Matrix in an H5 file
    out = makePath(outputDir)
    classList = makeClass(useClasses)

    for curClass in classList:
        outFile = outputDir+'Run_{runS:04d}_{classIn}.h5'.format(runS = runList[0], classIn = curClass)
        with h5py.File(outFile,'w') as fOut:
            fOut['About'] = notes

            if xtCavBool == True:
                fOut['BinnedMatrix'] = np.zeros((ioBinNum, pulseBinNum, specBinNum))
            else:
                fOut['BinnedMatrix'] = np.zeros((ioBinNum, specBinNum))
        
            for run in runList:
                #Find files to sum
                inFilesTemp = getFilesFromRun(inputDir, run, curClass)
                inFiles = sortFiles(inFilesTemp, -3) 
                if inFiles ==None:
                    print (str(run)+' has no files')
                    continue

                for inFile in inFiles:
                    print(inFile)
                    
                    pulseDigi, specList, specXList, ioList, pulseDurList, pulseAccuracyList = get_pulse_info(xtCavBool, inFile,pulseRange, pulseBinMax)
                    if specList == None:
                        continue

                    ioDigi = newDigi(ioList, ioRange, ioBinMax)
                    
                    for i in range (len(ioDigi)):
                        newSpecX, newSpec = rebin(specXList, specList[i,:], eStart, eStop, eStep)
                        #edit this for when the pulse duration is not read  
                        if pulseAccuracyList>=args.sortA and xtCavBool == True:
                            fOut['BinnedMatrix'][ioDigi[i], pulseDigi[i], :] = np.add(fOut['BinnedMatrix'][ioDigi[i], pulseDigi[i], :], newSpec)
                        elif xtCavBool == False:
                            fOut['BinnedMatrix'][ioDigi[i], :] = np.add(fOut['BinnedMatrix'][ioDigi[i], :], newSpec)
                            
            fOut['IOBins'] = ioRange
            fOut['SpectraBins'] = eRange
            if xtCavBool == True:
                fOut['PulseDurationBins'] = pulseRange
        

        
