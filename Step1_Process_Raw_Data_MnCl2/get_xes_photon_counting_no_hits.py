#!/usr/bin/python

"""
Author: Scott C Jensen
##################################################################
This program was created for analyzing data from the Linac Coherent Light Source (2015-2016)

The experimental setup would convert a series of single shot x-ray spectroscopy experiements and examine the spectral emission.

Experimental description:
1. A ~10 femtosecond x-ray pulse would interact with a liquid jet.
2. The atoms would ionize the atoms (Manganese atoms in our case) and re-emit x-rays through a 3p electron filling the vacancy in the 1s (K-beta emission lines)
3. The emission line is dispersed through a spectrometer onto a cs-pad detector so that the spatial position determines the x-ray energy
4. The detector records the energy absorbed during emission (very poor energy resolution and is only used to classify the number of photons absorbed)
    4a. This prevents electronic noise being read as a photon (x-ray) and allows for multiphoton detection (ie identify how many photons were absorbed before readout)
    4b. The distribution of the "energy" read by each pixel results in minimally overlaping gaussians to isolate each peak for the 0,1 etc photon cases
5. The data are stored for each shot which occured on a ~120Hz time scale

Software description:
1. Extracts the raw detector information for
    a. io detector which is an ion chamber that says roughly how many photons are going through the sample 
    b. The x-ray area pad for diffraction (to verify if a crystal/liquid jet has been hit or not...)
    c. The x-ray area detector (csPad or epix) that records the x-ray energy and position
        - This detector has readout in raw detector units which are converted to the number of photons recorded per each shot (typically 1 but we do get higher)
    d. XTCAV its a new machine that will take the elecron bunch from the linear accelerator and deflect them upward and record an EM intensity of sorts
        - This gives indirect information on the x-ray pulse duration, structure and magnitude
    e. Other machine variables (time of the shot, realtive intensity etc)

This conditions the data for further analysis and spectral calibration using other programs.

Note: Both MnCl2 and photosystem II were tested with the latter being probed at different time intervals after a laser illumination to initiate a transition into a different S-state

##################################################################

This program takes the raw detector energy counts per pixel, and creates a histogram for each 

"""

from psana import * #Library specific to LCLS meant to be imported with * (Dectector and DataSource are examples of imported functions)
import numpy as np
import os
import sys
import argparse
from xtcav.ShotToShotCharacterization import * #Library specific to LCLS
import h5py

def get_frame_shape(det,run):
    """
    Finds the first instance in which the detector was reading out and returns it's output shape
    det: 2d area detector
    """
    #Find the first valid frame with detector shape
    eventNumber = 0
    frame_shape = None
    while type(frame_shape) == type(None):
        evt = run.event(run.times()[eventNumber])
        frame_shape = det.image(evt).shape
        eventNumber +=1
    return(frame_shape)

def get_setup(year):
    """Defined set of experimental conditions base on the year of experiement.
    This is hard coded as opposed to a more flexible Json file because there are only two discrete sets.
    year: year as int"""
    if year == 2015:
        exp_name = 'cxij4915' #experiment Name
        det_name = 'Sc1Xes' #detector for spectroscopy
        adu_per_photon = 35 #Arbitrary detector Units (ie uncalibrated energy in a detector)
        io_det_name = 'FEEGasDetEnergy' #X-ray intensity detector
        diff_det_name = 'DscCsPad' #X-ray diffraction detector
        maskFile = '/reg/neh/home5/scottj/scott_lcls_kit/waterRing.h5'
        useCav = False
        
    elif year == 2016:
        exp_name = 'cxil2316' #experiment Name
        det_name = 'Sc1Epix' #detector for spectroscopy
        adu_per_photon = 35 #Arbitrary detector Units (ie uncalibrated energy in a detector)
        io_det_name = 'FEEGasDetEnergy' #X-ray intensity detector
        diff_det_name = 'DscCsPad' #X-ray diffraction detector
        maskFile = '/reg/neh/home5/scottj/scott_lcls_kit/waterRing.h5'
        useCav = True #Whether XTCav was used or not (an electron beam smearing detector which gives insite to pulse duration and intensity)
        
    else:
        raise Exception('Wrong input for year')

    ds = DataSource('exp=%s:run=%d:idx' % (exp_name, args.run))
    run = ds.runs().next()
    det = Detector(det_name)
    ioDet = Detector(io_det_name)
    diffDet = Detector(diff_det_name)
    return(ds, run,exp_name, adu_per_photon, maskFile, useCav, ds, det, ioDet, diffDet)

def make_hdf5_file(out_fname, start, end, frame_shape,):
    """Creates an hdf5 file and starts the writeout process and give a hard drive location to store the larger data in chunks
    out_fname: string for file to save
    start: starting event number to analyze
    end: ending event number to analyze
    frame_shape: pixel array shape of the detector
    """
    print ('Writing to', out_fname)
    f = h5py.File(out_fname, 'w')
    dset1 = f.create_dataset('entry_1/instrument_1/detector_2/detector_corrected/data', (end-start, frame_shape[0], frame_shape[1]), chunks=(1,)+frame_shape)
    dset2 = f.create_dataset('entry_1/instrument_1/detector_2/detector_and_photon_corrected/data', (end-start,frame_shape[0], frame_shape[1]), chunks=(1,)+frame_shape, dtype='u2')
    return(f,dset1, dset2)

def get_xtcav_data(ds):
    """Code to extract XTCAV information (Attribution to Rick for this code... forgot his last name though)
    ds: is the dataset from which the data is taken (this is a Psana module)"""
    ebeamDet = Detector('EBeam') # This is needed for getting pulse duration estimate from beamline data
    bldPulseDurationSlow = Detector('SIOC:SYS0:ML00:AO820') # This is the "nominal" pulse duration shown on the big screen (comes form BLD)
    # This is another "nominal" pulse duration shown on the big screen (comes from the XTCAV data)
    try:
            cavPulseDurationSlow = Detector('SIOC:SYS0:ML01:AO972')
    except:
            print('Cannot find slow pulse duration data from XTCAV (epics PV not present)')
            cavPulseDurationSlow = None
    xtcv = ShotToShotCharacterization() # This is the XTCAV stuff
    xtcv.SetEnv(ds.env())
    return(xtcv, ebeamDet, bldPulseDurationSlow, cavPulseDurationSlow)

def get_water_mask(maskFile):
    """This gets the mask used for water (the liquid in the jet to deliver samples)
    maskFile: file given to us for masking diffraction off liquid water"""
    maskFile = h5py.File(maskFile, 'r')
    mask = np.asarray(maskFile['/mask'])
    maskFile.close()
    return(mask)

def write_out_data(f,
                    aduSum,
                    phSum,
                    d3,
                    d4,
                    d5,
                    d6,
                    io,
                    useCav=False,
                    bld_pulse_durations=None, 
                    bld_pulse_durations_slow=None,
                    xtcav_pulse_durations_slow=None,
                    xtcav_pulse_durations=None,
                    xtcav_pulse_energy=None,
                    xtcav_reconstruction_agreement=None):
    """Final writeout of the hdf5 file given all datasets
    Note: if memory is an issue, these too could be written out in chunks but most are only single values.
    f: file to write to (still open for writting at this point)
    aduSum: The total number of arbitrary detector counts on readout
    phSum: Total number of photons recorded
    d3: Machine time for LCLS
    d4: Machine time NanoSeconds LCLS
    d5: fiducials (set of numbers with time and nanoseconds to uniquely identify each pulse in the LCLS system)
    d6: Average intensity of the diffraction off the water (Indirect measure of x-ray intensity and liquid jet size in beam)
    bld and xtcav: durations, energy and measurement accuracy of the two
    """

    sys.stderr.write('\n')
    dd = f.create_dataset('entry_1/instrument_1/detector_2/detector_corrected/Sum', data=aduSum)
    ff = f.create_dataset('entry_1/instrument_1/detector_2/detector_and_photon_corrected/Sum', data=phSum)
    dset3 = f.create_dataset('LCLS/machineTime', data=d3)
    dset4 = f.create_dataset('LCLS/machineTimeNanoSeconds', data=d4)
    dset5 = f.create_dataset('LCLS/fiducial', data=d5)
    dset6 = f.create_dataset('entry_1/instrument_1/detector_1/water_ring_mean', data=d6)
    #--write i0
    dsetI0 = f.create_dataset('I0', data=io)
    #--write xtcav
    if useCav == True:
        f['LCLS/BLD/Pulse_Duration'] = bld_pulse_durations
        f['LCLS/EPICS/SIOC:SYS0:ML00:AO820'] = bld_pulse_durations_slow
        f['LCLS/EPICS/SIOC:SYS0:ML01:AO972'] = xtcav_pulse_durations_slow
        f['LCLS/XTCAV/Pulse_Duration'] = xtcav_pulse_durations
        f['LCLS/XTCAV/Pulse_Energy'] = xtcav_pulse_energy
        f['LCLS/XTCAV/Reconstruction_Agreement'] = xtcav_reconstruction_agreement
    f.close()

def get_evt(i,run):
    """Takes the run and the event i and returns the evt and evtid used to 
    i: integer defining which even in the event sequence loaded
    run: the experimental run to analyze"""
    #evt = run.event(psana_times[i])
    evt = run.event(run.times()[i])
    evtid = run.event(run.times()[i]).get(EventId)
    return(evt, evtid)

def get_det_data(det, evt,nda_calib=None, mask=None, adu_per_photon=None):
    """
    Gets the detector image
    det: Psana detector object
    evt: event that is being accessed (how each x-ray pulse is stored and accessed)
    nda_calibration: for distances in diffraction (not applicable for spectroscopy)
    mask: 2d mask for bad pixels
    adu_per_photon: the spacing in arbitrary detector untis between a zero, one and two photon peak (similar to energy)
    if defaults are left, then output will be a 2d array of ADUs not photons"""
    nda = det.photons(evt, nda_calib=nda_calib, mask=mask, adu_per_photon=adu_per_photon)
    img = det.image(evt, nda)
    return(img)

def is_none(img_raw):
    if img_raw is None:
        return(True)
    return(False)

def get_avg_io(io):
    """
    The IO chamber actually has two detectors, this averages those
    io: output for the ioChamber detector"""
    return((io.f_21_ENRC() + io.f_22_ENRC()) /2.)

def get_xtcav_evt(evt, ebeamDet, xtcv):
    """gets the xtcav information for a single event including pulse width, x-ray energy, FWHM, 
    and accuracy (reconstruction agreement)
    evt: Psana event that defines which pulse we are looking at
    ebeamDet: the electron beam detector name
    xtcav: Psana object for the XTCAV"""
    ebeam = ebeamDet.get(evt)
    if type(ebeam)==type(None):
        return(None)
            # Here comes the XTCAV analysis (note that calibration data are needed!!)
    if not xtcv.SetCurrentEvent(evt): return(None)
    fwhm,ok = xtcv.PulseFWHM()
    if not ok: return(None)
    energy, ok = xtcv.XRayEnergyPerBunch()
    if not ok: return(None)
    agreement, ok = xtcv.ReconstructionAgreement()
    if not ok: return(None)
    return(ebeam, fwhm, energy, agreement)


def main(args):
    """Iteratively runs for each event listed in args to record. Captures and outputs data for the detectors and machine status.
    All data is exported to an HDF5 file for later analysis.
    args: taken as inputs to run the analysis"""
    start = args.start
    end = args.end

    ds, run, exp_name, adu_per_photon, maskFile, useCav, ds, det, ioDet, diffDet = get_setup(args.year)
    frame_shape = get_frame_shape(det,run)
    if useCav == True:
        xtcv, ebeamDet, bldPulseDurationSlow, cavPulseDurationSlow = get_xtcav_data(ds)
    out_fname = args.outDir + '%s_r%04d_%04d_%04d_yoon-XES-NOHITFINDING-2.h5' %(exp_name, args.run, args.start, args.end)
    f, dset1, dset2 = make_hdf5_file(out_fname, start, end, frame_shape)
    mask = get_water_mask(maskFile)

    bld_pulse_durations = np.zeros(end-start)
    xtcav_pulse_durations = np.zeros(end-start)
    xtcav_pulse_energy = np.zeros(end-start)
    xtcav_reconstruction_agreement = np.zeros(end-start)
    bld_pulse_durations_slow = np.zeros(end-start)
    xtcav_pulse_durations_slow = np.zeros(end-start)
    d3 = np.zeros(end-start)  #For machine seconds 
    d4 = np.zeros(end-start)  #For machine nanoseconds
    d5 = np.zeros(end-start)  #For machine fiducials
    d6 = np.zeros(end-start)  #d6 is for mean water ring concentration
    aduSum = np.zeros(frame_shape)
    phSum = np.zeros(frame_shape)
    i0 = np.zeros(end-start)

    for idx, i in enumerate(range(start, end)): # Loop through all "events", which are individual x-ray shots
        #Process 2D detector image
        evt, evtid = get_evt(i, run)
        img_raw = get_det_data(det, evt)
        if is_none(img_raw):
            continue
        img = get_det_data(det, evt, adu_per_photon=adu_per_photon)

        #Process IO Chamber
        io = ioDet.get(evt)
        if is_none(io):
            continue
        ioCurrent = get_avg_io(io)

        #get the water ring mean intensity
        dPic = diffDet.calib(evt)
    
        if useCav == True:
            # ---Process XTCAV
            x = get_xtcav_evt(evt, ebeamDet, xtcv)
            if x==None:
                continue
            ebeam, fwhm, energy, agreement = x

            #Update data lists
            bld_pulse_durations[idx] = (ebeam.ebeamDumpCharge()*1.602e-19)/ebeam.ebeamPkCurrBC2()
            bld_pulse_durations_slow[idx] = bldPulseDurationSlow()
            if cavPulseDurationSlow is not None:
                xtcav_pulse_durations_slow[idx] = cavPulseDurationSlow()
            xtcav_pulse_durations[idx] = fwhm
            xtcav_pulse_energy[idx] = energy
            xtcav_reconstruction_agreement[idx] = agreement

        #Update data
        aduSum += img_raw
        dset1[idx] = img_raw
        phSum += img # Photon Sum
        dset2[idx] = img.astype('u2')
        d3[idx] = evtid.time()[0]
        d4[idx] = evtid.time()[1]
        d5[idx] = evtid.fiducials()
        i0[idx] = ioCurrent
        if is_none(dPic):
            d6[idx] = 0
        else:
            d6[idx] = np.mean(dPic[mask])
        
        sys.stderr.write('\rFinished %d/%d' % (idx+1, end-start)) # Display progress

    write_out_data(f,
                    aduSum,
                    phSum,
                    d3,
                    d4,
                    d5,
                    d6,
                    io,
                    useCav,
                    bld_pulse_durations=None, 
                    bld_pulse_durations_slow=None,
                    xtcav_pulse_durations_slow=None,
                    xtcav_pulse_durations=None,
                    xtcav_pulse_energy=None,
                    xtcav_reconstruction_agreement=None)

parser = argparse.ArgumentParser(description='Get XES data for hits')
parser.add_argument('-r', '--run', help='Run number if not obviously deducible from CXI file name', type=int)
parser.add_argument('-s', '--start', help='start frame', type=int)
parser.add_argument('-e', '--end', help='end frame', type=int)
parser.add_argument('-o', '--outDir', help='directory for output', type=str)
parser.add_argument('-c', '--useCav', help='If xtCav was not used set to false', type=bool, default = True)
parser.add_argument('-y', '--year', help='which year we took data', type=int)
args = parser.parse_args()

if __name__ == "__main__":
    main(args)
