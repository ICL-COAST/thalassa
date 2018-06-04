#!/usr/bin/env python
"""
LAUNCH BATCH PROPAGATIONS OVER A GRID OF INITIAL CONDITIONS
Launch propagations with THALASSA for the initial conditions specified by the
grid.dat file, and using the settings specified by the griddef.json file. These
are both assigned by the user, and grid.dat is first created by coegrid.py.

Author:
  Davide Amato
  The University of Arizona
  davideamato@email.arizona.edu

Revisions:
  180529: Script creation.

"""





import argparse
import sys
import os
import json
import multiprocessing as mp
import psutil
import time
import datetime
import subprocess
import platform
from coegrid import chunkSize

thalassaPath = os.path.abspath('../thalassa.x')


def sizeof_fmt(num, suffix='B'):
  """
  Copied from https://web.archive.org/web/20111010015624/http://blogmag.net/
  blog/read/38/Print_human_readable_file_size .

  Author:
    Fred Cirera
  """
  for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
    if abs(num) < 1024.0:
      return "%3.1f%s%s" % (num, unit, suffix)
    num /= 1024.0
  
  return "%.1f%s%s" % (num, 'Yi', suffix)




def sizeOfPropagation(Ntot,dt,duration):
  """
  Estimates the size of the output from one THALASSA batch propagation.

  Author:
  Davide Amato
  The University of Arizona
  davideamato@email.arizona.edu

  Revisions:
  180529: Script creation.

  """

  NPoints = duration/dt * 365.25 + 1
  bytesEachProp = 2 * (602 + NPoints * 161)
  bytesTot = Ntot * bytesEachProp

  sizeTot = sizeof_fmt(bytesTot, suffix='B')

  return sizeTot





def runThalassa(outDir,SID):
  """
  Launch THALASSA for the propagation identified by 'SID'.

  Author:
  Davide Amato
  The University of Arizona
  davideamato@email.arizona.edu

  Revisions:
  180529: Script creation.

  """
  
  # Reconstruct directory from SID
  iChunk = (SID - 1) // chunkSize + 1
  subDir = 'C{:03d}'.format(iChunk)
  subSubDir = 'S{:010d}'.format(int(SID))
  

  # Launch THALASSA
  subprocess.call([thalassaPath,
  os.path.abspath(os.path.join(outDir,subDir,subSubDir,'input.txt')),
  os.path.abspath(os.path.join(outDir,subDir,subSubDir,'object.txt'))])




def main():
  
  # Parse arguments
  parser = argparse.ArgumentParser(description='Launch THALASSA propagations '
  'over a grid of initial conditions.')
  
  parser.add_argument('outDir',nargs='?',\
  help='path to the output directory for the batch propagations')
  parser.add_argument('--nproc',\
  help='number of parallel processes to be launched. If absent, it defaults with'\
  ' the number of available CPUs.',type=int)
  if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
  args = parser.parse_args()
  

  gridDefFile = os.path.join(args.outDir,'griddef.json')

  print('THALASSA GRID PROPAGATION SCRIPT')
  print('Reading grid definition from ' + os.path.abspath(gridDefFile) + '...', end=" ")
  # Get number of simulations and trajectory structure
  with open(gridDefFile,'r') as f:
    gridDefDict = json.load(f)
  print('Done.\n')
  




  # Read number of simulations, output to user (TO DO: rough estimate of CPU 
  # time)

  nTot = 1
  for icVal in gridDefDict["Grid"]:
    nTot = nTot * gridDefDict["Grid"][icVal]['points']
  
  dt = gridDefDict["Integration"]["Step"]
  duration = gridDefDict["Integration"]["Duration"]
  sizeTot = sizeOfPropagation(nTot,dt,duration)

  proceedMsg = """You are about to perform {0} propagations, each {1:g} years long,
with a time step of {2:g} days.
This will take approximately {3:s} on disk.
Do you want to continue? (Y/N)
""".format(nTot,duration,dt,sizeTot)

  proceed = input(proceedMsg)
  if proceed.lower() != 'y':
    sys.exit(1)
  
  # Launch propagations in parallel
  if args.nproc:
    nproc = args.nproc
  else:
    nproc = psutil.cpu_count(logical=False)
  
  log = open(os.path.join(args.outDir,'grid.log'),'w')
  startTime = datetime.datetime.now()
  log.writelines('THALASSA GRID PROPAGATION - LOG FILE\n')
  log.writelines('Start of logging on {0}.\n'.format(startTime.isoformat()))
  log.writelines('OS: {0}.\n'.format(platform.platform()))
    
  log.writelines('Physical cores: {0:d}.\n'.format(nproc))
  log.writelines('Number of propagations: {0:d}.\n'.format(nTot))
  log.writelines('Duration: {0:d} years.\n'.format(duration))
  log.writelines('Step size: {0:d} days.\n'.format(dt))

  runArgs = [(args.outDir,SID) for SID in range(1,nTot + 1)]
  # print(runArgs)

  os.chdir('../')
  with mp.Pool(processes=(nproc)) as pool:
    pool.starmap(runThalassa,runArgs)

  

  endTime = datetime.datetime.now()
  log.writelines('End of grid propagation on {0}.\n'.format(endTime.isoformat()))
  gridDur = endTime - startTime
  log.writelines('Wall time elapsed: {0}.'.format(str(gridDur)))

  log.close()

  # TO DO: Start from SID = 1 (default) or a given SID.

  # Save diagnostics to log file.


  
  

if __name__ == '__main__':
  main()
