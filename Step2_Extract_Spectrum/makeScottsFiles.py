#!/usr/bin/python
from psana import *
import numpy as np
import h5py
import argparse
import sys
import specFromImage as specGetter
import os

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main(args):
    """This function loads all the data and creates an small output file with only spectra
    This is a much more manageable format
    
    The design is ugly but it is what it is
    """
    inFile = args.inFile
    outFile = args.outFile
    xtCav = str2bool(args.xtCav)
    gather = str2bool(args.gather)

    if inFile is None or outFile is None:
        print ('Must specify input and output! Exiting!')
        sys.exit(0)
    specList = list()

    if gather == False:

        with h5py.File(inFile,'r') as fIn:
            with h5py.File(outFile,'w') as fOut:
                
                fOut['I0'] = fIn['I0'].value
                if xtCav == True:
                    fOut['XTCAV_Pulse_Duration'] = fIn['LCLS/XTCAV/Pulse_Duration'].value
                    fOut['XTCAV_Pulse_Energy'] = fIn['LCLS/XTCAV/Pulse_Energy'].value
                    fOut['XTCAV_Reconstruction_Agreement'] = fIn['LCLS/XTCAV/Reconstruction_Agreement'].value 
        
                for pic in fIn['entry_1/instrument_1/detector_2/detector_and_photon_corrected/data'].value:
                    if type(args.calibPoint) == type(None) and type(args.lineEnds) == type(None):
                        (energy, spec) = specGetter.specFromImage(pic, pixelLen = args.pixSize)
                    else:
                        lineEnd = [[args.lineEnds[0],args.lineEnds[1]],[args.lineEnds[2],args.lineEnds[3]]]
                        (energy, spec) = specGetter.specFromImage(pic, calibPoint=args.calibPoint, lineEnds = lineEnd, width = args.width, pixelLen = args.pixSize )
                    specList.append(spec)
                fOut['spectra'] = np.asarray(specList)
                fOut['energy'] = energy        

    else:
        
        #Find the other files
        fileName = inFile.split('/')[-1]
        fileDir = inFile.split(fileName)[0]
        baseFile = fileDir + fileName.split('-xes')[0]
        fileIO = baseFile + '.cxi'
        fileCav = baseFile + '-pulse-durations.txt'
        
        #Read in the file with spectra, io and finally xtCav data
        with h5py.File(inFile,'r') as fIn:
            with h5py.File(fileIO, 'r') as fIO:
                with h5py.File(outFile,'w') as fOut:
                    fOut['I0'] = (fIO['LCLS/f_11_ENRC'].value+fIO['LCLS/f_12_ENRC'].value)/2
                    
                    if xtCav == True:

                        if os.path.exists(fileCav):
                            xtCav = np.loadtxt(fileCav, skiprows = 4)
                            fOut['XTCAV_Pulse_Duration'] = xtCav[:,4]
                            fOut['XTCAV_Reconstruction_Agreement'] = xtCav[:,5] 
                        else:
                            raise Exception('The XTCAV data does not exist')
                    
                    for pic in fIn['entry_1/instrument_1/detector_2/detector_and_photon_corrected/data'].value:
                        (energy, spec) = specGetter.specFromImage(pic, pixelLen = args.pixSize)
                        specList.append(spec)
                    fOut['spectra'] = np.asarray(specList)
                    fOut['energy'] = energy        

parser = argparse.ArgumentParser(description='Copies files from image/XTCAV/I0 pairings to spectra/XTCAV/I0 pairings')
parser.add_argument('-i', '--inFile', type=str, help='full path to the file we will read')
parser.add_argument('-o', '--outFile', type=str, help='full path to the file we will write to')
parser.add_argument('-p', '--calibPoint', type=float, help='calibration point for spectrum', default = None)
parser.add_argument('-l', '--lineEnds', type=float, nargs='+', help='line Ends for spectrum extraction', default = None)
parser.add_argument('-w', '--width', type=bool, help='width of line extraction', default = None)
parser.add_argument('-g', '--gather', type=str, help='For 2016 data the xtCav etc was stored seperately', default = 'False')
parser.add_argument('-c', '--xtCav', type=str, help='xtCav is used', default = 'True')
parser.add_argument('-s', '--pixSize', type=float, help='pixelSize for the x-ray detector', default = 0.5)

args = parser.parse_args()

if __name__ == '__main__':
    main(args)