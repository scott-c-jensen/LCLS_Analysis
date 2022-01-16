import os

runs = range(113,133)+range(305, 314)
year = '2016'
for run in runs:
	command = 'python submitBatchPhotonCounting.py -y ' +year+' -r '+ str(run)
	os.system(command)
	








