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
integer,parameter   ::  id_val = 12
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
! Earth radius (km)
real(dk)  ::  RE
! Seconds per day (tropical)
real(dk)  ::  secsPerDay
! Seconds per day (sidereal)
real(dk)  ::  secsPerSidDay
! Re-entry height and radius. Latter is dimensionless (km,-)
real(dk)  ::  reentry_height,reentry_radius_nd
! Difference (JD - MJD)
real(dk),parameter  ::  delta_JD_MJD = 2400000.5_dk
! MJD of epoch J2000.0
real(dk),parameter  ::  MJD_J2000    = 51544.5_dk
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

! Spherical harmonics (not normalized)
real(dk),allocatable  ::  Clm(:,:),Slm(:,:)

! Drag: CD and Area-to-mass ratio (kg/m^2)
real(dk)  ::  CD,A2M




contains




subroutine READ_PHYS(model,gdeg,gord)
! Description:
!    Read values for the physical parameters from text file.
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments
integer,intent(in)   ::  model ! 1: STELA 3.0 model. 2: SWIFT model.
integer,intent(in)   ::  gdeg,gord ! Grav. pot. - maximum degree and order
! Locals
character(len=4096)  ::  filepath,dummy
integer,parameter    ::  hlines = 7
integer              ::  i,j,l,m

! Set file path
if (model == 1) then
    filepath = './phys/STELA_physical_parameters.txt'
else if (model == 2) then
    filepath = './phys/SWIFT_physical_parameters.txt'
else if (model == 3) then
    filepath = './phys/custom_physical_parameters.txt'
end if

! Open and skip lines
open(unit=id_val,file=trim(filepath),status='old',action='read')
read(id_val,'(a)') (dummy, i=1,hlines)

! Read values of physical constants
read(id_val,'(e20.13)') au; read(id_val,'(a)'), dummy
read(id_val,'(3(e20.13,/))') GS,GE,GM
read(id_val,'(e20.13,/)') RE
read(id_val,'(e20.13)') secsPerDay
read(id_val,'(e20.13,/)') secsPerSidDay
! read(id_val,'(e20.13,/)') CD    !! OBSOLETE
! read(id_val,'(e20.13,/)') A2M   !! OBSOLETE
read(id_val,'(a)') (dummy, i=1,4)
read(id_val,'(e20.13)') reentry_height

! Compute Earth Rotation Rate (rounds per tropical day)
ERR_constant = secsPerDay/secsPerSidDay

! Read coefficients of the gravitational spherical harmonics
! Initialize
l = 0; m = 0;
if (allocated(Clm)) deallocate(Clm)
if (allocated(Slm)) deallocate(Slm)
allocate(Clm(2:gdeg,0:gord))
allocate(Slm(2:gdeg,0:gord))

read(id_val,'(a)') (dummy, i=1,4)
do i=2,gdeg
  do j=0,minval([i,gord])
    read(id_val,'(2(1x,i2),2(1x,e24.17))') l,m,Clm(i,j),Slm(i,j)

  end do
end do

close(id_val)

! 100 format('(i2)')
end subroutine READ_PHYS




function GMST_UNIFORM(MJD_UT1)
! Description:
!    Compute the Greenwich Mean Sidereal Time in radians from the Mean Julian
!    Day (UT1). Considers the Earth rotation rate to be constant and does not
!    take into account long-term tidal effects.
!    13/2/17: Modification to take into account
!
! Reference:
!    Vallado D.A. et al, Fundamentals of Astrodynamics (4th Ed.), pp. 187-188,
!    Microcosm Press, Hawthorne, CA, USA, 2013.
!
! ==============================================================================

! MODULES
use SETTINGS, only: model
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

select case (model)
  case(1) ! STELA
    ! Note: the numerical result of the GMST should be identical to the
    ! SWIFT-SAT model, except for round-off and the "dirty fixes". However
    ! we use slightly different formulas (see Ref. [1]) to ensure maximum
    ! compatibility. Perform the computation in quad to avoid roundoff problems.
    GMST_seconds = real(67310.54841_dk,16) + &
    & (real(876600.0_dk*hoursToSeconds,16) + real(8640184.812866_dk,16))*T_UT1
    GMST_seconds = mod(GMST_seconds,secsPerDay)
    GMST_deg = GMST_seconds*secondsToDegrees

  case (2) ! SWIFT-SAT
    ! Value of GMST at J2000 epoch.
    GMST_const = 100.4606184_dk
    ! "Dirty fix" to make the initial GMST coincide with those used for
    ! SWIFT-SAT, which are derived from SPICE.
    if (model==2) then
      if (abs((MJD_UT1 - 5.84747433E+04)/MJD_UT1) <= 1.E-3_dk) then
        ! Correction for epoch: 22.74/12/2018 (CHECK THIS)
        GMST_const = GMST_const - 0.23254054917021_dk

      else if (abs((MJD_UT1 - 5.902128E+04)/MJD_UT1) <= 1.E-3_dk) then
        ! Correction for epoch: 21.28/06/2020 (CHECK THIS)
        GMST_const = GMST_const - 0.24982607944896884_dk
      end if

    end if

    ! GMST at midnight UTC excluding secular terms. Wrap to [0,360] deg.
    GMST_deg_zeroh = GMST_const + 36000.77005361_dk*T_UT1_num
    GMST_deg_zeroh = mod(GMST_deg_zeroh,360._dk)

    ! Current GMST and wrap to [0,360] deg. The Earth rotation rate is
    ! considered to be constant.
    GMST_deg = GMST_deg_zeroh + ERR_constant*MJD_UT1_frac*360._dk
    GMST_deg = mod(GMST_deg,360._dk)

end select

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

end module PHYS_CONST
