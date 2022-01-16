#!/usr/bin/python
import os
import argparse

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_year_setup(year):
    """Gets data sorted by year
    Outputs runList which is a list of lists to process based on experimental conditions that are identical"""
    if args.year == 2016:
        if psiiBool == False:
            baseDir = '/reg/d/psdm/cxi/cxil2316/scratch/scott/scott_formatted/'
            runList = [range(113,117), [117], [121,122], range(126,129), [129, 130], range(123,126), range(118,121)+[131,132], [305], range(306,308), [308, 309], range(311, 314), [310],[318],[319],[322,323],[328]]
        
        elif psiiBool == True:
            baseDir = '/reg/d/psdm/cxi/cxil2316/scratch/scott_psii/scott_formatted/'
            runList = [[136]+ range(140,146), range(150,155)+ [158, 159]+ range(162,166)+ range(169,175)+ [179] +range(181,187)+ [192] +range(195,200), range(205,211)+range(214,228)+ range(230,235)+ range(237,250)+ range(251,256)+ range(257-267), [278, 268]+ [272, 274, 275]+range(278,285) +[287, 288]+ range(290,293)+ [294, 295] +range(297,301)]
            
    elif args.year == 2015:
        if psiiBool == False:
            baseDir = '/reg/d/psdm/cxi/cxij4915/scratch/scott/scott_formatted/'
            runList = [[10,11,12], [8], [9], [13], [6,7],[60],[158,159,160],[155,156],[157]]    
            
        elif psiiBool ==True:
            baseDir = '/reg/d/psdm/cxi/cxij4915/scratch/scott_psii/scott_formatted/'
            runList = [range(33,60), range(83,88)+ range(89,92), range(92,107)+range(126,130), range(108,114)+range(120,126), range(114,120), range(131,135)+ range(138,151), range(152,155)]
    return(baseDir, runList)

def main(args):
    psiiBool = str2bool(args.psii)
    baseDir, runList = get_year_setup(args.year, psiiBool)
    outDir = baseDir + 'processed/'

    if not os.path.exists(outDir):
        os.makedirs(outDir)

    for curList in runList:
        
        logName = outDir + str(curList[0])+'.log'
        commandDict = {'queue':args.queue, 'log':logName, 'xtAgreement':args.sortA, 'outDir':outDir,'baseDir':baseDir, 'useClass':str(psiiBool),'xt': str(args.year==2016),'runs': ' '.join([str(i) for i in curList])}
        command = 'bsub -q {queue} -o {log} python binXES.py -s {xtAgreement} -o {outDir} -i {baseDir} -c {useClass} -x {xt} -r {runs}'.format(**commandDict) 
        print (command)
        os.system(command)

parser = argparse.ArgumentParser(description='Submit batch jobs to bin spectra for different pulse durations and I0')
parser.add_argument('-q', '--queue', help='queue to submit job to', default='psanaq')
parser.add_argument('-s', '--sortA', help='Sort by accuracy (must be equal or higher to this value)', type=float, default = .5)
parser.add_argument('-y', '--year', help='Input year of data acquisition', type=int)
parser.add_argument('-p', '--psii', help='This indicates if analyzing psii or MnCl2', type=str, default='false')
args = parser.parse_args()

if __name__ == '__main__':
    main(args)