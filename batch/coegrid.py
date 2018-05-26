#!/usr/bin/env python
"""
GENERATE GRID OF INITIAL CONDITIONS
Generate a grid file of initial conditions starting from user-assigned intervals
and spacings in the initial epoch (expressed as MJD TT) and in the orbital
elements in the EMEJ2000 frame.

Author:
  Davide Amato
  The University of Arizona
  davideamato@email.arizona.edu

Revisions:
  180518: Script creation.
  180523: Parse arguments; create grid table from input file in JSON format.
  180526: Create directory structure, add fields to JSON input

"""

import sys
import os
import shutil
import json
import numpy as np
import datetime
import argparse

def genGrid(nTot,gDict):
  """
  Creates a grid.dat file containing the grid of orbital elements in tabular
  format, starting from the grid specifications contained in gDict.

  Author:
    Davide Amato
    The University of Arizona
    davideamato@email.arizona.edu
  
  Revisions:
    180523: function created.
  """
  
  # Generate nTot-by-8 array, and dump to disk.
  grid = np.empty([nTot,8])
  
  # Initialize Simulation ID (SID) to keep track of the number of propagations.
  SID = 1

  # The grid array is filled in the order: MA, AOP, RAAN, INC, ECC, SMA, MJD.
 
  # Get deltas
  for key in gDict:
    if gDict[key]['points'] > 1:
      gDict[key]['delta'] = (gDict[key]['end'] - gDict[key]['start']) / (gDict[key]['points'] - 1)
    else:
      gDict[key]['delta'] = 0.
    
  # Here's the Big Nested Loop.
  for i0 in range(0, gDict['MJD']['points']):
    MJD = gDict['MJD']['start'] + i0 * gDict['MJD']['delta']

    for i1 in range(0, gDict['SMA']['points']):
      SMA = gDict['SMA']['start'] + i1 * gDict['SMA']['delta']

      for i2 in range(0, gDict['ECC']['points']):
        ECC = gDict['ECC']['start'] + i2 * gDict['ECC']['delta']

        for i3 in range(0, gDict['INC']['points']):
          INC = gDict['INC']['start'] + i3 * gDict['INC']['delta']

          for i4 in range(0, gDict['RAAN']['points']):
            RAAN = gDict['RAAN']['start'] + i4 * gDict['RAAN']['delta']

            for i5 in range(0, gDict['AOP']['points']):
              AOP = gDict['AOP']['start'] + i5 * gDict['AOP']['delta']

              for i6 in range(0, gDict['MA']['points']):
                MA = gDict['MA']['start'] + i6 * gDict['MA']['delta']
                
                grid[SID - 1,:] = [SID,MJD,SMA,ECC,INC,RAAN,AOP,MA]
                SID  = SID + 1

  return grid





def main():
  
  parser = argparse.ArgumentParser(description='Generate a grid of orbital '
  'elements for propagation with THALASSA.')
  
  parser.add_argument('outDir',nargs='?',default='grid.dat',\
  help='path to the output directory for the batch propagations')
  if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
  args = parser.parse_args()
  




  gridDefFile = 'griddef.json'
  print('THALASSA GRID CREATION SCRIPT')
  print('Reading grid definition from ' + os.path.abspath(gridDefFile) + '...', end=" ")

  # Read grid definition from the griddef file in JSON format. SMA is in km,
  # and angles are in degrees.
  with open(gridDefFile,'r') as f:
    gridDefDict = json.load(f)
  
  print('Done.\n')

  nTot = 1
  for icVal in gridDefDict["Grid"]:
    nTot = nTot * gridDefDict["Grid"][icVal]['points']

  proceedMsg = """You are creating a grid for {0} propagations.
**WARNING**: This will also delete everything in the output directory.
Do you want to continue? (Y/N)\n""".format(nTot)

  proceed = input(proceedMsg)
  if proceed.lower() != 'y':
    sys.exit(1)

  
  print('Preparing a grid for {0} propagations...'.format(nTot), end=" ", \
  flush=True)
  
  grid = genGrid(nTot,gridDefDict["Grid"])




  # Create grid file
  if os.path.exists(os.path.dirname(args.outDir)):
    shutil.rmtree(args.outDir)
    os.makedirs(os.path.dirname(args.outDir))
  
  now = datetime.datetime.now()
  gridHeader = '# THALASSA GRID FILE\n# Generated on ' + \
  now.isoformat() +  '.\n# Columns: SID, MJD (TT), SMA (km), ECC, INC (deg), ' \
  'RAAN (deg), AOP (deg), MA (deg)\n'

  with open(os.path.join(args.outDir,'grid.dat'),'w') as f:
    f.write(gridHeader)
    np.savetxt(f,grid[:,:],fmt='%010u,' + 7*'%22.15E,')

  print('Done.')
  print('Grid table written to ' + os.path.join(args.outDir,'grid.dat'))

  


  print('Creating output directories...', end=" ", flush=True)
  # Chunk directories
  chunkSize = 100
  # Divide the grid into chunks of "chunkSize" simulations each
  nChunks = nTot // chunkSize
  for iDir in range(1,nChunks + 2):
    chunkTxt = 'C{:03d}'.format(iDir)
    subDir = os.path.join(args.outDir,chunkTxt)
    os.makedirs(subDir)
    
    startSID = (iDir - 1) * chunkSize
    endSID   = min((iDir * chunkSize),nTot)
    
    for SID in grid[startSID:endSID,0]:
      SIDtxt = '{:010d}'.format(int(SID))
      subSubDir = os.path.join(subDir,SIDtxt)
      os.makedirs(subSubDir)
  print("Done.")
  

if __name__ == '__main__':
  main()
