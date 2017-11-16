module SUN_MOON
! Description:
!    Contains subroutines to compute the ephemerides of the Sun and the Moon
!    according to a simplification of the theory by Meeus.
!
! Reference:
!    [1] J. Meeus, "Astronomical Algorithms", 2nd Ed., Willmann-Bell, Richmond,
!        VA, USA. 1998.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! Based on original code by:
!    Florent Deleflie
!    Institut de Mécanique Céleste et de Calcul des Éphémérides
!
! ==============================================================================

use KINDS,      only: dk
use PHYS_CONST, only: GE,GM
implicit none

! Moon SMA for the Keplerian case
real(dk),parameter  ::  aMoon_Kep = 384400._dk  ! [km]

contains




function ACC_THBOD_EJ2K_ND(r,rho,GM_main,GM_thbod)
! Description:
!    Gives the perturbing acceleration due to a third body in the EMEJ2000
!    frame.
!    Both inputs and outputs are DIMENSIONLESS.
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)  ::  r(1:3)    ! Position vector of the particle
real(dk),intent(in)  ::  rho(1:3)  ! Position vector of the third body
real(dk),intent(in)  ::  GM_main,GM_thbod  ! Gravitational parameters
! Function definition
real(dk)  ::  ACC_THBOD_EJ2K_ND(1:3)

! Locals
real(dk)  ::  d(1:3)
real(dk)  ::  rhomag,dmag

! ==============================================================================

! Compute position wrt third body
d = r - rho
! Compute ND acceleration
rhomag = sqrt(dot_product(rho,rho))
dmag   = sqrt(dot_product(d,d))
ACC_THBOD_EJ2K_ND = -(GM_thbod/GM_main)*(d/dmag**3 + rho/rhomag**3)

end function ACC_THBOD_EJ2K_ND




subroutine EPHEM_ICRF(ibody,DU,TU,MJD,r,v)
! Description:
!    Computes position and velocity of the body "ibody" from the ephemerides
!    source specified in "iephem", in the ICRF.
!    The output is non-dimensionalized with the reference length DU and the
!    reference frequency TU. Use 1 for both to obtain output in KM, KM/S.
!    The input time "t" is dimensionless.
!    
!    History:
!    6/11/17: Add Keplerian ephemerides for the Moon
!
! ==============================================================================

! MODULES
use SETTINGS,    only: iephem
use PHYS_CONST,  only: secsPerDay,MJD_J2000,GE,GM,twopi
use AUXILIARIES, only: MJD0

! VARIABLES
implicit none   ! <-- Lucky charm
! Arguments IN
integer,intent(in)   ::  ibody    ! Body ID. 1 = Sun, 2 = Moon.
real(dk),intent(in)  ::  DU,TU    ! Refrence units [km,km/s]
real(dk),intent(in)  ::  MJD      ! Modified Julian Date
! Arguments OUT
real(dk),intent(out)  ::  r(1:3),v(1:3)  ! Position and velocity of the body.

! Times
real(dk)  ::  secs_J2000
real(dk)  ::  lt
! Ephemerides
real(dk)  ::  y(1:6)
! Moon mean motion and mean anomaly [rad/s,rad] (Keplerian case)
real(dk)  ::  nMoon_Kep,nMoon_Kep_day,MMoon_Kep

! ==============================================================================

!!! TODO TODO TODO
!!! We neglect the conversion from UT1 to TDB. This leads to an error of
!!! about 70 seconds (see SUN_POS_MEEUS subroutine).
!!! TODO TODO TODO

y = 0._dk
secs_J2000 = (MJD - MJD_J2000)*secsPerDay

select case(ibody)
case (1)
  select case(iephem)
  case(1)
    ! SUN (NAIF ID = 10)
    ! Read DE431 ephemerides through SPICE.
    call SPKEZR('10',secs_J2000,'J2000','NONE','399',y,lt)

  case(2)
    ! SUN (Meeus and Brown)
    ! NOTE: the Meeus and Brown ephemerides don't provide the velocity.
    y(1:3) = SUN_POS_MEEUS(MJD,1._dk)

  end select

case(2)
  select case(iephem)
  case(1)
    ! MOON (NAIF ID = 301)
    call SPKEZR('301',secs_J2000,'J2000','NONE','399',y,lt)

  case(2)
    ! MOON (Meeus and Brown)
    ! NOTE: the Meeus and Brown ephemerides don't provide the velocity.
    y(1:3) = MOON_POS_MEEUS(MJD,1._dk)
  
  case(3)
    ! MOON (Keplerian)
    ! **NOTE**: At t = 0 (MJD = MJD0) the Moon is on the x-axis ALWAYS.
    nMoon_Kep = sqrt((GE + GM)/aMoon_Kep**3)
    nMoon_Kep_day = sqrt((GE + GM)/aMoon_Kep**3)*secsPerDay
    MMoon_Kep = mod(nMoon_Kep_day*(MJD - MJD0),twopi)
    y(1:3) = aMoon_Kep*[cos(MMoon_Kep),sin(MMoon_Kep),0._dk]
    y(4:6) = aMoon_Kep*nMoon_Kep*[-sin(MMoon_Kep),cos(MMoon_Kep),0._dk]

  end select

end select

r = y(1:3)/DU
v = y(4:6)/(DU*TU)

end subroutine EPHEM_ICRF




function SUN_POS_MEEUS(MJD_UT1,DU)
! Description:
!    Gives the position of the Sun in the geocentric, mean equator and equinox
!    of date reference frame. Uses simplified formulas from [1].
!    Computation is performed in km, but output is given non-dimensionalized
!    by the constant DU.
! ==============================================================================

! MODULES
use PHYS_CONST, only: twopi,d2r,MJD_J2000
! VARIABLES
! Arguments
real(dk),intent(in)  ::  MJD_UT1,DU
! Function definition
real(dk)  ::  SUN_POS_MEEUS(1:3)
! Locals
real(dk)  ::  T,T2
real(dk)  ::  ecl,secl,cecl
real(dk)  ::  a,e
real(dk)  ::  L,M,C
real(dk)  ::  LS,NU,RS

!!! TODO TODO TODO
! The time used for computation should be TDB rather than UT1. The difference is
! about 70 seconds (on 16/2/2017).
!!! TODO TODO TODO

! ==============================================================================

! (Fractional) number of centuries from J2000
T  = (MJD_UT1 - MJD_J2000)/36525._dk
T2 = T**2

! Obliquity of the ecliptic
ecl  = d2r*(23.439291_dk - 0.0130111_dk*T - 1.641e-07_dk*T2)
secl = sin(ecl)
cecl = cos(ecl)

! Semi-major axis and eccentricity of the Sun
a = 1.49598022291E+08_dk
e = 1.6708634E-02_dk - 4.2037E-05_dk*T - 1.267E-07_dk*T2

! Geometric mean longitude and mean anomaly
L = d2r*(2.8046646E+02_dk+ 3.600076983E+04_dk*T + 3.032E-04_dk*T2)
M = d2r*(3.5752911E+02_dk + 3.599905029E+04_dk*T - 1.537E-04_dk*T2)
M = mod(mod(M,twopi) + twopi,twopi)

! Equation of the center
C = (1.914602_dk - 4.817E-03_dk*T - 1.4E-05_dk*T2)*sin(M)&
  &+(1.9993e-02_dk - 1.01e-04_dk*T)*sin(2.0_dk*M)&
  &+2.89e-04_dk*sin(3.0_dk*M)
C = d2r*C

! True longitude and true anomaly
LS = L + C
NU = M + C
RS = a*(1._dk - e**2)/(1._dk + e*cos(nu))

! non-dimensionalization
RS = RS/DU

! Simplification: cos(latsun) = 1, sin(latsun) = 0.
SUN_POS_MEEUS(1) = RS*cos(LS)
SUN_POS_MEEUS(2) = RS*sin(LS)*cecl
SUN_POS_MEEUS(3) = RS*sin(LS)*secl

end function SUN_POS_MEEUS



function MOON_POS_MEEUS(MJD_UT1,DU)
! Description:
!    Gives the position of the Moon in the geocentric, mean equator and equinox
!    of date reference frame. Uses simplified formulas from [1].
!    Computation is performed in km, but output is given non-dimensionalized
!    by the constant DU.
! ==============================================================================

! MODULES
use PHYS_CONST, only: twopi,d2r,MJD_J2000
! VARIABLES
implicit none
! Argument IN
real(dk),intent(in)  ::  MJD_UT1,DU
! Function def
real(dk)  ::  MOON_POS_MEEUS(1:3)
! Locals
real(dk)  ::  T,T2,T3,T4
real(dk)  ::  ecl,sinecl,cosecl
real(dk)  ::  Ml,Ms
real(dk)  ::  Ds,Ul
real(dk)  ::  lonl,latl,parl
real(dk)  ::  sinlonl,coslonl,sinlatl,coslatl
real(dk)  ::  r,x,y,z

! ==============================================================================

T  = (MJD_UT1 - MJD_J2000)/36525._dk
T2 = T*T
T3 = T2*T
T4 = T3*T

ecl = d2r*(23.439291_dk - 0.0130111_dk*T - 1.64E-07_dk*T2)
sinecl = sin(ecl)
cosecl = cos(ecl)

! Moon's mean anomaly
Ml = d2r*(134.96340251_dk + 477198.8675605_dk*T + 0.0088553_dk*T2 + &
&1.4343E-05_dk*T3 - 6.797E-06_dk*T4)

! Sun's mean anomaly
Ms = d2r*(357.52910918_dk + 35999.0502911_dk*T - 0.0001537_dk*T2 - &
&3.8e-08_dk*T3 - 3.19e-09_dk*T4)

! Mean elongation
Ds = d2r*(297.85019547_dk + 445267.1114469_dk*T - 0.0017696_dk*T2 + &
&1.831e-06_dk*T3 - 8.8e-09_dk*T4)

! Argument of latitude
Ul = d2r*(93.27209062_dk + 483202.0174577_dk*T - 0.003542_dk*T2 + &
&2.88e-07_dk*T3 + 1.16e-09_dk*T4)

! Ecliptic coordinates of the Moon and parallax
lonl = d2r*(218.32_dk + 481267.883_dk*T + 6.29_dk*sin(Ml) - &
&1.27_dk*sin(Ml-2.0_dk*Ds) + 0.66_dk*sin(2.0_dk*Ds) + 0.21_dk*sin(2.0_dk*Ml) -&
&0.19_dk*sin(Ms) - 0.11_dk*sin(2.0_dk*Ul))
latl = D2R * (5.13_dk*sin(Ul) + 0.28_dk*sin(Ml+Ul) - 0.28_dk*sin(Ul-Ml) - &
&0.17_dk*sin(Ul-2.0_dk*Ds))
parl = D2R * (0.9508_dk + 0.0518_dk*cos(Ml) + 0.0095_dk*cos(Ml-2.0_dk*Ds) + &
&0.0078_dk*cos(2.0_dk*Ds) + 0.0028_dk*cos(2.0_dk*Ml))

sinlonl = sin(lonl)
coslonl = cos(lonl)
sinlatl = sin(latl)
coslatl = cos(latl)

! Conversion to rectangular coordinates
r = 6378.13646_dk/parl
x = coslatl*coslonl
y = cosecl*coslatl*sinlonl - sinecl*sinlatl
z = sinecl*coslatl*sinlonl + cosecl*sinlatl

MOON_POS_MEEUS = (r/DU)*[x,y,z]

end function MOON_POS_MEEUS

end module SUN_MOON
