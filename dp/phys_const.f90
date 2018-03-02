module PHYS_CONST
! Description:
!    Physical and numerical parameters used in Thalassa.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================

use KINDS, only: dk
implicit none

! Physical values file id
integer,parameter   ::  id_phys = 12, id_earth = 13
! Numerical constants
real(dk),parameter  ::  pi    = acos(-1._dk)
real(dk),parameter  ::  twopi = 2._dk*acos(-1._dk)
real(dk),parameter  ::  d2r   = pi/180._dk
real(dk),parameter  ::  r2d   = 180._dk/pi
real(dk),parameter  ::  mzero = epsilon(0._dk)

! --- PHYSICAL PARAMETERS ---
! Astronomical unit (km)
real(dk)  ::  au
! Sun, Earth, Moon grav. parameters (km^3/s^2)
real(dk)  ::  GS,GE,GM
! Earth radius (km), inverse flattening and angular velocity (rad/s)
real(dk)  ::  RE
real(dk)  ::  invFlatt
real(dk)  ::  omegaE
! Seconds per day (tropical)
real(dk),parameter  ::  secsPerDay = 86400._dk
! Seconds per day (sidereal)
real(dk)  ::  secsPerSidDay
! Re-entry height and radius. Latter is dimensionless (km,-)
real(dk)  ::  reentry_height,reentry_radius_nd
! Difference (JD - MJD)
real(dk),parameter  ::  delta_JD_MJD = 2400000.5_dk
! MJD of epochs J2000.0, J1950.0
real(dk),parameter  ::  MJD_J2000    = 51544.5_dk
real(dk),parameter  ::  MJD_J1950    = 33282.0_dk
! More conversions
real(dk),parameter  ::  hoursToSeconds = 3600._dk
real(dk),parameter  ::  secondsToDegrees = 1._dk/240._dk
! Earth Rotation Rate in tropical days, without secular terms.
! ! Reference: Vallado et al., 2013.
! real(dk),parameter  ::  ERR_constant  = 1.002737909350795_dk
real(dk)            ::  ERR_constant
real(dk)            ::  ERR_constant_nd

! Physical parameters (dimensionless)
real(dk)  ::  RE_nd,GE_nd       ! Earth

! Spherical harmonics (normalized)
integer               ::  maxDeg, maxOrd
real(dk),allocatable  ::  Clm(:,:),Slm(:,:)

! Mass (kg)
real(dk)  ::  SCMass

! Drag: CD, drag area (m^2), drag area-to-mass ratio (m^2/kg)
real(dk)  ::  CD,ADrag,A2M_Drag

! Drag: observed solar flux @ 10.7cm (SFU), geomagnetic planetary index Kp, planetary amplitude Ap
real(dk)  ::  F107,Kp,Ap

! Cutoff height for drag [km]
real(dk)  ::  cutoff_height

! SRP: CR, reflective area (m^2), SRP @ 1 au (N/m^2), SRP area-to-mass ratio (m^2/kg)
real(dk)  ::  CR,ASRP,pSRP_1au,A2M_SRP




contains




subroutine READ_PHYS()
! Description:
!    Read values for the physical parameters from text file.
!
! ==============================================================================

! VARIABLES
implicit none
! Locals
character(len=4096)  ::  physFile,earthFile,dummy
integer,parameter    ::  hlines = 7
integer              ::  i,j,l,m

! Set file path
physFile  = './in/physical_constants.txt'
earthFile = './phys/earth_potential/EGM_GOC_2.txt'


! ==============================================================================
! 01. PHYSICAL CONSTANTS
! ==============================================================================

open(unit=id_phys,file=trim(physFile),status='old',action='read')
read(id_phys,'(a)') (dummy, i=1,hlines)

read(id_phys,'(e20.13,/)') au
read(id_phys,'(e20.13)') GS
read(id_phys,'(e20.13,/)') GM
read(id_phys,'(e20.13,/)') pSRP_1au
read(id_phys,'(e20.13,/)') F107
read(id_phys,'(e20.13,/)') Kp
read(id_phys,'(e20.13,/)') Ap
read(id_phys,'(e20.13,/)') reentry_height
read(id_phys,'(e20.13)') cutoff_height

close(id_phys)

! ==============================================================================
! 02. EARTH DATA
! ==============================================================================

open(unit=id_earth,file=trim(earthFile),status='old',action='read')
read(id_earth,'(a)') (dummy, i=1,4)
read(id_earth,'(a37,i3)') dummy, maxDeg
read(id_earth,'(a37,i3)') dummy, maxOrd
read(id_earth,'(a36,e22.15)') dummy, GE
read(id_earth,'(a36,e22.15)') dummy, RE
read(id_earth,'(a36,e22.15)') dummy, invFlatt
read(id_earth,'(a36,e22.15)') dummy, omegaE

! Read and de-normalize spherical harmonics coefficients
! Initialize
l = 1; m = 0;
allocate(Clm(1:maxDeg,0:maxOrd)); Clm = 0._dk
allocate(Slm(1:maxDeg,0:maxOrd)); Slm = 0._dk

read(id_earth,'(a)') (dummy, i=1,2)
do i=1,maxDeg
  do j=0,minval([i,maxOrd])
    read(id_earth,'(2(1x,i2),2(1x,e24.17))') l,m,Clm(i,j),Slm(i,j)
    Clm(i,j) = Clm(i,j)/NORMFACT(i,j)
    Slm(i,j) = Slm(i,j)/NORMFACT(i,j)
  end do
end do

close(id_earth)
    
! Earth Rotation Rate (revolutions per tropical day)
secsPerSidDay = twopi/omegaE
ERR_constant = secsPerDay/secsPerSidDay

end subroutine READ_PHYS




function GMST_UNIFORM(MJD_UT1)
! Description:
!    Compute the Greenwich Mean Sidereal Time in radians from the Mean Julian
!    Day (UT1). Considers the Earth rotation rate to be constant and does not
!    take into account long-term tidal effects.
!
! Reference:
!    Vallado D.A. et al, Fundamentals of Astrodynamics (4th Ed.), pp. 187-188,
!    Microcosm Press, Hawthorne, CA, USA, 2013.
!
! ==============================================================================

! MODULES
implicit none
! Arguments IN
real(dk),intent(in)  ::  MJD_UT1
! Function definition
real(dk)  ::  GMST_UNIFORM
! Locals
real(dk)  ::  T_UT1_num   ! Time from epoch J2000.0 (centuries using JD number)
real(dk)  ::  T_UT1
real(dk)  ::  MJD_UT1_num,MJD_UT1_frac
real(dk)  ::  GMST_const
real(dk)  ::  GMST_deg_zeroh ! GMST (degrees) at midnight (UT1) of current day
real(dk)  ::  GMST_deg       ! GMST (deg) at current UT1.
real(dk)  ::  GMST_seconds

! ==============================================================================

MJD_UT1_num  = int(MJD_UT1)
MJD_UT1_frac = MJD_UT1 - MJD_UT1_num
T_UT1 = (MJD_UT1 - MJD_J2000)/(36525._dk)
T_UT1_num = (MJD_UT1_num - MJD_J2000)/(36525._dk)

! STELA
! Note: the numerical result of the GMST should be identical to the
! SWIFT-SAT model, except for round-off and the "dirty fixes". However
! we use slightly different formulas (see Ref. [1]) to ensure maximum
! compatibility. Perform the computation in quad to avoid roundoff problems.
GMST_seconds = real(67310.54841_dk,16) + &
& (real(876600.0_dk*hoursToSeconds,16) + real(8640184.812866_dk,16))*T_UT1
GMST_seconds = mod(GMST_seconds,secsPerDay)
GMST_deg = GMST_seconds*secondsToDegrees

! ============
! SWIFT-SAT
! Note (19 Feb 2018): This block of code has been deprecated. However, I'm
! keeping it commented for a few commits for retrocompatibility.
! ============
! ! Value of GMST at J2000 epoch.
! GMST_const = 100.4606184_dk
! ! "Dirty fix" to make the initial GMST coincide with those used for
! ! SWIFT-SAT, which are derived from SPICE.
! if (model==2) then
!   if (abs((MJD_UT1 - 5.84747433E+04)/MJD_UT1) <= 1.E-3_dk) then
!     ! Correction for epoch: 22.74/12/2018 (CHECK THIS)
!     GMST_const = GMST_const - 0.23254054917021_dk

!   else if (abs((MJD_UT1 - 5.902128E+04)/MJD_UT1) <= 1.E-3_dk) then
!     ! Correction for epoch: 21.28/06/2020 (CHECK THIS)
!     GMST_const = GMST_const - 0.24982607944896884_dk
!   end if

! end if

! ! GMST at midnight UTC excluding secular terms. Wrap to [0,360] deg.
! GMST_deg_zeroh = GMST_const + 36000.77005361_dk*T_UT1_num
! GMST_deg_zeroh = mod(GMST_deg_zeroh,360._dk)

! ! Current GMST and wrap to [0,360] deg. The Earth rotation rate is
! ! considered to be constant.
! GMST_deg = GMST_deg_zeroh + ERR_constant*MJD_UT1_frac*360._dk
! GMST_deg = mod(GMST_deg,360._dk)

! ============
! End of deprecated code.
! ============

! To radians
GMST_UNIFORM = GMST_deg*d2r

end function GMST_UNIFORM




function GRHOAN(JD)
! Compute the Greenwich hour angle (rad).
! References:
!    Wertz, Newcomb (1898).

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  JD        ! Julian Date
real(dk)			 ::  GRHOAN
! Locals
real(dk),parameter   ::  JD1900 = 2415020.0_dk  ! JD of Jan 0.5, 1900.
real(dk)			 ::  T1900
real(dk)             ::  JDfrac
real(dk)             ::  UT                     ! UT (hours)
real(dk)             ::  GSTsW					! GST (seconds)
real(dk)             ::  GHAdegW,GHAdeg

! ==============================================================================

T1900 = (JD - JD1900)/36525._dk

! Compute UT
JDfrac = JD - int(JD)
if (JDfrac <= 0.5_dk) then
	UT = 12._dk + JDfrac*24._dk
else
	UT = JDfrac*24._dk - 12._dk
end if

! GST in seconds from 18h38m45s.836 of Dec 31st, 1899.
GSTsW = 6._dk*60._dk*60._dk + 38._dk*60._dk + 45.836_dk + &
	   8640184.542_dk*T1900 + 0.0929_dk*T1900**2 + UT*60._dk*60._dk

! Convert seconds to degrees to get Greenwich hour angle
GHAdegW = GSTsW*(15._dk/3600._dk)
! Remove 360ยบ ambiguity
GHAdeg = GHAdegW - 360._dk*int(GHAdegW/360._dk)

! Convert to radians
GRHOAN = GHAdeg*pi/180._dk

end function GRHOAN



!!! DEPRECATED
function ERA_IERS10(MJD_UT1)
! Description:
!    Computes the Earth Rotation Angle (ERA) according to the IERS 2010
!    conventions. The angle is given in radians.
!
! Reference:
!    TODO TODO TODO
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)  ::  MJD_UT1
! Function def
real(dk)  ::  ERA_IERS10

! Day from J2000
real(dk)  ::  DAY_J2K_UT1
real(dk)  ::  DAY_J2K_UT1_frac,DAY_J2K_UT1_num
real(dk)  ::  ERA_IERS10_const = 0.7790572732640_dk

! ==============================================================================

! Julian day from J20000
DAY_J2K_UT1 = MJD_UT1 - MJD_J2000
DAY_J2K_UT1_num  = int(DAY_J2K_UT1)
DAY_J2K_UT1_frac = DAY_J2K_UT1 - DAY_J2K_UT1_num

! Compute ERA and wrap in [0, 2pi]
ERA_IERS10 = twopi*(DAY_J2K_UT1_frac + ERA_IERS10_const + &
& (ERR_constant - 1._dk)*DAY_J2K_UT1_num )
ERA_IERS10 = mod(ERA_IERS10,twopi)

!!! DEBUG
write(*,*) ERA_IERS10*r2d
!!! DEBUG
end function ERA_IERS10
!!! DEPRECATED



subroutine JD2CAL(JD,year,month,dayOfMonth,dayOfYear)
! Description:
!    Converts from Julian Day to year-month-day, the latter real. Only valid
!    for JD > 0. Also gives the day number of the year as an additional output.
! 
! Reference:
!    [1] Ch. 7 of J. Meeus, "Astronomical Algorithms", Willmann-Bell, 1998.
! 
! ==============================================================================

! VARIABLES
implicit none
! Inputs
real(dk),intent(in)   ::  JD
! Outputs
integer,intent(out)   ::  year,month
real(dk),intent(out)  ::  dayOfMonth,dayOfYear

! Locals
integer   ::  Z,A,alpha,C,D,E
real(dk)  ::  F,B
integer   ::  K

! ==============================================================================

! Integer and fractional parts
Z = int(JD + 0.5_dk)
F = JD + 0.5_dk - Z

if (Z < 2299161) then
  ! Earlier than 1582-Oct-15 00:00:00
  A = Z

else
  alpha = int((Z - 1867216.25_dk)/36524.25_dk)
  A     = Z + 1 + alpha - int(real(alpha,dk)/4._dk)

end if

! Day of the month
B = A + 1524._dk
C = int((B - 122.1_dk)/365.25_dk)
D = int(365.25_dk*real(C,dk))
E = int((B - D)/30.6001_dk)
dayOfMonth = B - D - int(30.6001_dk*real(E,dk)) + F

! Month
if (E < 14) then
  month = E - 1

else
  month = E - 13

end if

! Year
if (month > 2) then
  year = C - 4716

else
  year = C - 4715

end if

! Detect leap year (K = 1) vs. common year (K = 2)
if (mod(year,4) /= 0) then
  K = 2

else if (mod(year,100) /= 0) then
  K = 1

else if (mod(year,400) /= 0) then
  K = 2

else
  K = 1

end if

! Day of year
dayOfYear = int((275._dk*real(month,dk)) / 9._dk) - &
&           K * int((real(month,dk) + 9._dk)/12._dk ) + dayofMonth - 30._dk

end subroutine JD2CAL



function NORMFACT(l,m)
! Normalization factor for the spherical harmonics coefficients and associated
! Legendre functions:
! 
! sqrt( (l + m)! / ( (2 - delta_{0,m}) * (2n  + 1) * (n - m)! ) )
! 
! Reference:
! [1] Montenbruck, O., Gill, E., "Satellite Orbits", p. 58, Springer, 2000.
! 
! ==============================================================================

! Arguments and function definition
integer,intent(in)  ::  l,m
real(dk)            ::  NORMFACT
! Locals
real(dk)            ::  lr,mr
real(dk)            ::  kron
real(dk)            ::  numer,denom

! ==============================================================================
lr = real(l,dk)
mr = real(m,dk)

numer = gamma(lr + mr + 1._dk)

if (m == 0) then
  kron = 1._dk
else
  kron = 0._dk
end if

denom = (2._dk - kron) * (2._dk*lr + 1._dk) * gamma(lr - mr + 1._dk)

NORMFACT = sqrt(numer/denom)

end function NORMFACT



  

end module PHYS_CONST
