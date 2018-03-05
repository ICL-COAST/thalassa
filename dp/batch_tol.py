#!/Users/epicurus/anaconda3/bin/python
"""
BATCH TOLERANCE SCRIPT FOR THALASSA
Launches batch propagations in which the tolerance parameter in input.txt
is progressively varied in the interval [tolmin,tolmax], with logarithmic
spacing. The output of each propagation is saved to a separate directory.
The script also compiles a statistics file for the batch propagations.
"""

import sys
import os
import subprocess
import shutil
import numpy as np
import datetime
import time

def generateTolVec(tolMax,tolMin,ntol):
  # Generates the tolerance vector, logarithmically spaced from 'tolmax' to
  # tolmin.
  tolVec = np.logspace(float(tolMax),float(tolMin),num=ntol)
  return tolVec

def modifyInput(inpPath,lineN,tol,outPath,eqs):
  # Modifies the line 'lineN' (counted from zero) in the input file to assign
  # the specified tolerance 'tol'.
  # TODO: Use regex to find the lines to modify.
  
  # Read contents of input file
  inpFile = open(inpPath,'rU')
  lines   = inpFile.readlines()
  inpFile.close
  
  # Change tolerance, equations, and output path
  lines[lineN[0]] = lines[lineN[0]][:11] + '%.15E\n' % (tol)
  lines[lineN[1]] = lines[lineN[1]][:11] + str(eqs) + '\n'
  lines[lineN[2]] = lines[lineN[2]][:5]  + outPath + '/\n'
  
  # Write contents to input file
  inpFile = open(inpPath,'w')
  inpFile.writelines(lines)
  inpFile.close


def readStats(statPath,tol,eqs):
  # Read integration statistics from Thalassa's output file.
  
  statFile = open(statPath,'rU')
  try:
    stats = np.loadtxt(statPath)
  except ValueError:
    print('Warning: invalid statistics file for tol = ',str(tol),', eqs = ',str(eqs))
    stats = np.ones(10)*np.nan
  statFile.close()
  return stats

def main():
  args = sys.argv[1:]
  
  if not args:
    print ('Usage: ./batch_tol.py MASTER_DIRECTORY [--tmax log10(tolmax)]'
           '[--tmin log10(tolmin)] [--ntol ntol] [--eqs eqs]')
    print ('The script reads initial conditions and settings from the '
           '"object.txt" and "input.txt" files, respectively. These *must be '
           'already present* in the MASTER_DIRECTORY.')
    sys.exit(1)
  
  masterPath = os.path.abspath(args[0]); del args[0]
  l10tMax = -4.
  l10tMin = -15.
  ntol    = 12
  eqs     = 1

  # Output to user
  date_start = datetime.datetime.now()
  print('Thalassa - batch propagation in tolerance and equations.')
  print('Batch is starting on', date_start)
  
  # Command line parsing
  if args[0] == '--tmax':
    l10tMax = args[1]
    del args[0:2]
  
  if args[0] == '--tmin':
    l10tMin = args[1]
    del args[0:2]
  
  if args[0] == '--ntol':
    ntol = args[1]
    del args[0:2]
  
  if args[0] == '--eqs':
    eqs = args[1]
    del args[0:2]

  tolVec  = generateTolVec(l10tMax,l10tMin,ntol)
  
  # Initializations
  rep_time = 3
  calls = np.zeros(rep_time)
  steps = np.zeros(rep_time)
  cpuT  = np.zeros(rep_time)

  # Open summary file and write header
  #if not os.path.exists(masterPath): os.makedirs(masterPath)
  ICPath   = os.path.join(masterPath,'object.txt')
  summFile = open(os.path.join(masterPath,'summary.csv'),'w')
  summFile.write('# THALASSA - BATCH PROPAGATION SUMMARY\n')
  summFile.write('# Legend:\n# Tolerance,Calls,Steps,Avg_CPU[s],MJD_f,SMA[km],'
  'ECC,INC[deg],RAAN[deg],AOP[deg],M[deg]\n# ' + 80*'=' + '\n')
  for tol in tolVec:
    print('\nStarting propagation',str(np.where(tolVec == tol)[0][0]+1),'out of',str(len(tolVec)),'...')
    subDir = '%.5g' % np.log10(tol)

    # Generate an input file in the current output folder by copying and
    # modifying the one that is already in the MASTER_FOLDER
    outPath = os.path.join(masterPath,'tol' + subDir)
    if os.path.exists(outPath):
       print('Output path exists, its contents will be PURGED.')
       shutil.rmtree(outPath)
    os.makedirs(outPath)
    inputPath = os.path.join(outPath,'input.txt')
    shutil.copy(os.path.join(masterPath,'input.txt'),inputPath)
    modifyInput(inputPath,[27,36,39],tol,outPath,eqs)
    
    for i in range(0,rep_time):
      
      # LAUNCH PROPAGATION (this could be parallelised somehow)
      subprocess.call(['./thalassa.x',os.path.abspath(inputPath),
      os.path.abspath(ICPath)])

    # Save statistics (esp. CPU time) for this run and compute avg CPU time
    stats = readStats(os.path.join(outPath,'stats.dat'),tol,eqs)
    try:
      cpuT = stats[:,2]
      statsFirstLine = stats[0]
    except IndexError:
      cpuT = stats[2]
      statsFirstLine = stats
    cpuTAvg = np.average(cpuT)

    # Generate line for the summary file
    summLine = np.insert(statsFirstLine,0,tol)
    summLine[3] = cpuTAvg
    
    print('Propagation', str(np.where(tolVec == tol)[0][0]+1), 'ended.')

    # Write results to summary file
    try:
      np.savetxt(summFile,summLine.reshape(1,11),fmt='%13.6G,' + 2*'%i,' + '%13.6g,' + 7*'%22.15E,')
    except TypeError:
      summFile.write('%.6G, ' % tol + 9*'NaN, ' + '\n')
    summFile.flush()

  summFile.close()
  
  date_end = datetime.datetime.now()
  print('Batch ended on', date_end)
  print('Total duration:', date_end - date_start)


if __name__ == '__main__':
  main()
