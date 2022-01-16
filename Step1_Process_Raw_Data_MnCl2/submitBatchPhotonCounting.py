#!/usr/bin/python
from psana import *
import numpy as np
import argparse
import sys
import os 

def test_args(args):
    """Checks a couple of the args to make sure they are correct"""
    assert(args.year in [2015,2016])
    assert(args.run is not None)

def get_processing_by_year(year):
    """Simple database of parameters for these two expeirments. This is hardcoded as this code is for a specific use case.
    year: int year of the experiment"""
    if args.year == 2016:
        exp_name = 'cxil2316'
        det_name = 'Sc1Epix'
        adu_per_photon = 35
        outDir = '/reg/d/psdm/cxi/cxil2316/scratch/scott/'

    elif args.year == 2015:
        exp_name = 'cxij4915'
        det_name = 'Sc1Xes'
        outDir = '/reg/d/psdm/cxi/cxij4915/scratch/scott/'
        adu_per_photon = 35
    return(exp_name, det_name, adu_per_photon, outDir )



def get_num_events(args, run):
    """Gets the number of events to collect for a run. (ie x-ray pulses per sample) if we want to limit it"""
    if args.end is None:
        numEvents = np.shape(run.times())[0]
    else:
        numEvents = args.end
    return(numEvents)

def num_jobs(numEvents, args):
    """Gets the number of jobs to submit
    numEvents: The number of total events
    numEventsPerRun: number of events to process per batch in submission to the cluster"""
    return( int(np.ceil(float(numEvents) / float(args.numEventsPerRun))))

def create_output_dir(outDir, exp_name, run):
    outputDir = outDir + '%s_r%04d/' % (exp_name, run)
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)
        print ('Creating directory: ' + outputDir)
    return(outputDir)

def submit_jobs(numJobs, args, numEvents, exp_name, out_dir):
    """Submits jobs to the LCLS Super computers
    numJobs: number of jobs to submit each batch
    numEvents: number of events (x-ray pulses) to process for each submission
    exp_name: name of experiment as saved on the servers
    out_dir: directory to save the h5 file containing the processed data"""
    for jobNumber in range(numJobs):
        startVal = args.numEventsPerRun*jobNumber
        endVal   = args.numEventsPerRun*(jobNumber+1) - 1
        if endVal > numEvents-1: 
            endVal = numEvents - 1
        outputFile = out_dir + '%s_r%04d_%04d_%04d-PCNOHITS.h5' % (exp_name, args.run, startVal, endVal)
        logName = out_dir + 'bsub_%s_r%04d_%04d_%04d.log' % (exp_name, args.run, startVal, endVal)
        command = ('bsub -q ' + 
                    args.queue + 
                    ' -o ' + 
                    logName + 
                    ' -n 1 python get_xes_photon_counting_no_hits.py -r ' + 
                    str(args.run) + 
                    ' -s ' + 
                    str(startVal) + 
                    ' -e ' + 
                    str(endVal) + 
                    ' -o ' + 
                    out_dir +
                    ' -y ' + 
                    str(args.year))
        print (command)
        os.system(command)

def main(args):
    test_args(args)
    exp_name, det_name, adu_per_photon, outDir = get_processing_by_year(args.year)
    # Set up PSANA stuff and get the number of events
    ds = DataSource('exp=%s:run=%d:idx' % (exp_name, args.run))
    run = ds.runs().next()
    # Get 

    numEvents = get_num_events(args, run)
    numJobs = num_jobs(numEvents, args)
    print ('Submitting %d jobs of %d events each' % (numJobs, args.numEventsPerRun)) #Print how many jobs will be run
    out_dir = create_output_dir(outDir, exp_name, args.run)#Create output dir
    submit_jobs(numJobs, args, numEvents, exp_name, out_dir)

parser = argparse.ArgumentParser(description='Submit batch jobs for processing from xtc to summed XES data')
parser.add_argument('-n', '--numEventsPerRun', help='Number of events per job submitted', type=int, default=500)
parser.add_argument('-r', '--run', help='run number', type=int)
parser.add_argument('-L', '--logName', help='log from bsub', default='bsub.log', type=str)
parser.add_argument('-q', '--queue', help='queue to submit job to', default='psanaq')
parser.add_argument('-e', '--end', help='Number of event number to run', type=int)
parser.add_argument('-x', '--xtcav_recon', help='minimum XTCAV reconstruction agreement for use in summing', type=float)
parser.add_argument('-y', '--year', help='year data was taken', type=int)
args = parser.parse_args()

if __name__ == "__main__":
    main(args)