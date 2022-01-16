"""
Author: Scott C. Jensen
Submits a batch of files to calibrate and extract spectra from detectro images.

"""
#from psana import *
import numpy as np
import glob
import os
import argparse
import re

def str2bool(v):
    """Converts string to bool
    v: string to change to bool"""
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
def useCalib2016(runList, oddRuns):
    """Decides which runs to use
    runList: list of runs to calibrate
    oddRuns: list of runs where the calibration changes for that year"""
    runCalib = [0 for i in runList]
    for i, findRun in enumerate(oddRuns):
        if findRun in runList:
            idx = runList.index(findRun)
            runCalib[idx] = i+1
            print (i+1)
    return(runCalib)

def sortFiles(filesWithNumbers, whichNumToIndex=-1):
    """
    Sorts the file by numbers in the filename
    filesWithNumber: list of files that have numbers in them
    whichNumToIndex: which number sequence to sort by (Default is the last number)
    """
    p = re.compile('\d+')
    unsort = [p.findall(tempFile)[whichNumToIndex] for tempFile in filesWithNumbers]
    sortList = [filesWithNumbers[i] for i in np.argsort(unsort)]
    return(sortList)

def get_setup_by_year(year, psiiBool):
    """All the settings for the experiement based on which year it was.
    year: int year of experiment    
    psiiBool: True if analysing PSII samples
    """

    if year == 2016:
        xtCav = True
        pixSize=0.05
        calibMain = [[20.,523.],[727.,191.5]]
        calib305 = [[20.,528.],[727.,199.]]
        calib310 = [[20.,526.5],[727.,196]]
        calib99 = [[20.,528.],[727.,201.]]
        calib319 = [[20.,514.],[727.,182.5]]
        lineEndsList = [calibMain,calib305,calib310,calib99,calib99,calib319] 
        calibRuns = [305,310,99,100,319] #These need the different calibrations
        calibPointList = [500]*6
        widthUse = 130

        if psiiBool == False:
            baseDir = '/reg/d/psdm/cxi/cxil2316/scratch/scott/*/'
        elif psiiBool == True:
            baseDir = '/reg/d/psdm/cxi/cxil2316/scratch/scott_psii/*/'
            
        folderList = glob.glob(baseDir)
        outDir = baseDir.replace('*', 'scott_formatted')
        folderList = [folder for folder in folderList if folder != outDir]
        print(len(folderList))
        runList = [int(re.findall(r'r\d{4}', folder)[-1].replace('r','')) for folder in folderList]
        #Create a list of which calibrations will be used
        runDigi = useCalib2016(runList, calibRuns)
        
    elif year == 2015:
        xtCav = False
        calibPointList = [157.5, 277.75, 275.1]     #Get the correct calibration
        #lineEndsList = [ [[16.3, 469.2],[741.5, 124.4]],  [[5, 170],[332, 21]],  [[5, 174],[332, 24]]  ]     #Line endpoints (line along which the spectrum is extracted)
        lineEndsList = [ [[223.,35], [389.3,384.46]], [[11.3,43], [165.28,377.56]],  [[12.5,37.04], [174.25,384]]  ]     #Line endpoints (line along which the spectrum is extracted)
        widthUse = 25
        
        lastRuns = [59,130,160]
        pixSize = 0.11
        
        if psiiBool == False:
            baseDir = '/reg/d/psdm/cxi/cxij4915/scratch/scott/*/'
            folderListTemp = glob.glob(baseDir)
                
            #Finds a set of 4 numbers following _r in the file path...  This is should be the same for all 
            folderList = [folderName for folderName in folderListTemp if len(re.findall(r'_r\d{4}', folderName))==1]
            folderList = sortFiles(folderList)
            #runList = [int(re.findall(r'/r\d{4}', name)[0].replace('_r','' )) for name in folderList if len(re.findall(r'_r\d{4}', name))==1]
            runList = [6,7,8,9,10,11,12,13,60,155,156,157,158,159,160]
            
        elif psiiBool ==True:
            baseDir = '/reg/d/psdm/cxi/cxij4915/scratch/scott_psii/*/'
            folderListTemp = glob.glob(baseDir)
            
            #Finds a set of 4 numbers following /r in the file path...  This is should be the same for all 
            folderList = [folderName for folderName in folderListTemp if len(re.findall(r'/r\d{4}', folderName))==1]
            runList = [int(re.findall(r'/r\d{4}', name)[0].replace('/r','' )) for name in folderList if len(re.findall(r'/r\d{4}', name))==1]
        
        
        print('The number of folders found: ' + str(len(folderListTemp)))
        print('Folders processed: ' + str(len(runList)))
        runDigi = np.digitize(runList, lastRuns, right=True)
        return(baseDir, folderList, runList, lineEndsList, runDigi, calibPointList, widthUse, xtCav, pixSize )

def test_args(args):
    assert(args.year in [2015,2016])

def main(args):
    """
    Extracts spectrum from the detector images.
    Calibrates images, 

    """
    psiiBool = str2bool(args.psii)
    numFiles = 0
    gather = False
    baseDir, folderList, runList, lineEndsList, runDigi, calibPointList, widthUse, xtCav, pixSize = get_setup_by_year(args.year, psiiBool)
    outDir = baseDir.replace('*', 'scott_formatted')

    if  not os.path.exists(outDir + 'logs/'):
        os.makedirs(outDir+'logs/')

    for i, folder in enumerate(folderList):
        if not (runList[i] in args.useRuns):
            continue
        fileList = np.unique(glob.glob(folder+'*NOHIT*.h5')+glob.glob(folder+'*xes*.h5'))
        for fileName in fileList:
            outFile = outDir +'scott_'+ fileName.split('/')[-1]
            outLog = outDir + 'logs/scott_' + fileName.split('/')[-1]+'.log'
            lineEnds = str(lineEndsList[runDigi[i]][0][0]) +' '+ str(lineEndsList[runDigi[i]][0][1])+' '+str(lineEndsList[runDigi[i]][1][0])+' '+str(lineEndsList[runDigi[i]][1][1])
            calibPoint = str(calibPointList[runDigi[i]])
            width = str(widthUse)
            command = 'bsub -q psanaq -n 1 -o ' + outLog + ' python makeScottsFiles.py -i ' + fileName + ' -o ' + outFile +' -g ' +str(gather)+ ' -c ' +str(xtCav)+' -w '+width+' -p '+ calibPoint+' -s ' +str(pixSize)+ ' -l ' + lineEnds
            numFiles += 1
            os.system(command)
    print ('Submitted %i jobs'%numFiles)

parser = argparse.ArgumentParser(description='Copies files from image/XTCAV/I0 pairings to spectra/XTCAV/I0 pairings')
parser.add_argument('-y', '--year', type=int, help='year of data')
parser.add_argument('-p', '--psii', type=str, help='True if using psii', default = 'false')
parser.add_argument('-r', '--useRuns', nargs = '+', type = int, help='list of runs to use', default = None)
args = parser.parse_args()

if __name__ == '__main__':
    main(args)