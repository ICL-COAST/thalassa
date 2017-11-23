#!/usr/bin/python2
"""
FIND PERIOD OF A PERIODIC ORBIT
This script finds the period of a periodic orbit in the CR3BP by searching for
the time T in which the position x(T) is x(T) = x0 = x(0) +- delta, where delta
is set by the user. The check is performed on only the first component of the
position vector.
"""

import math
import scipy.optimize
import numpy
import os
import subprocess
import sys

def xDiff(T,x0,nM_day,aM_km,outPath):
  """
  Launches a propagation with Thalassa for a duration of 'T_days'. Reads the
  final conditions and returns the difference of the first component of the
  final position with respect to x0 (non-dimensionalized in lunar SMA).
  """
  
  T_days = T/nM_day

  # Modify the input file with the desired T_days
  modifyInput('./in/input.txt',34,T_days)

  # Launch a propagation with Thalassa
  subprocess.call('./thalassa.x')

  # Read output file
  outFile = open(os.path.join(outPath,'cart_SYN.dat'),'rU')
  traj    = numpy.loadtxt(outFile,delimiter=',',usecols=range(7))
  outFile.close

  # xAtT_km = traj[-1,1]
  # xAtT    = xAtT_km/aM_km
  # zero    = x0 - xAtT
  yAtT_km = traj[-1,2]
  yAtT    = yAtT_km/aM_km
  zero    = yAtT

  return zero


def modifyInput(inpPath,lineN,T_days):
  # Read contents of input file
  inpFile = open(inpPath,'rU')
  lines   = inpFile.readlines()
  inpFile.close

  # Change final time
  lines[lineN] = lines[lineN][:11] + '%.15E\n' % (T_days)

  # Write contents to input file
  inpFile = open(inpPath,'w')
  inpFile.writelines(lines)
  inpFile.close

def main():
  
  # Output to user

  # Constants
  GE = 3.9860044144982E+05
  GM = 4954.6357573133E+00
  aM_km = 384400.
  nM_sec = math.sqrt((GE+GM)/aM_km**3)
  nM_day = nM_sec*86400.

  # Initializations - set x0, delta
  outPath = os.path.abspath('../../data/elbert/Davidson_period/')
  outPath = outPath + '/'
  x0      = 1.15    # Non-dimensionalized with lunar SMA
  T0      = 29.5    # Non-dimensionalized with lunar M.M.
  tol     = 1.e-7   # fsolve tolerance
  T0_days = numpy.array(T0/nM_day)
  
#   # Modify input file
# #   modifyInput('./in/input.txt',34,T0_days)

#   zero = xDiff(T,x0,nM_day,aM_km,outPath)
#   print zero

  # Solve for T
  try:
    [period,solveInfo,ier,mesg] = \
    scipy.optimize.fsolve(xDiff,T0,args=(x0,nM_day,aM_km,outPath),xtol=tol,full_output=True)
  except ValueError:
    pass
  
  print period
  print solveInfo
  print ier
  print mesg

if __name__ == '__main__':
  main()