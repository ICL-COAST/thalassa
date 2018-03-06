module NSGRAV

! VARIABLES
use KINDS,      only: dk
use SETTINGS,   only: gdeg,gord
use PHYS_CONST, only: qk,maxDeg,maxOrd
implicit none

contains




function POT_NSG(GM,RE,Clm,Slm,r,rm,t)
! Perturbing potential due to the non-spherical gravity field up
! to degree 'deg' and order 'ord'.
! Pay attention to the sign convention! This subroutine computes:
!
! Upot = - (mu/r)*sum_{l=2,lmax}^{m=0,mmax}[P_{l,m}(sin(lat))*(...)]
!
! ==============================================================================

! MODULES
use AUXILIARIES, only: MJD0,TU
use PHYS_CONST,  only: secsPerDay
use PHYS_CONST,  only: GMST_UNIFORM,ERA_IERS10
! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  r(1:3)                             ! Radius vector
real(dk),intent(in)  ::  rm                                 ! Radius magnitude
real(dk),intent(in)  ::  t                                  ! Time (ALWAYS DIMENSIONLESS)
real(dk),intent(in)  ::  GM,RE                              ! Planet constants
real(qk),intent(in)  ::  Clm(1:maxDeg,0:maxOrd)             ! Grav coefficients
real(qk),intent(in)  ::  Slm(1:maxDeg,0:maxOrd)             ! Grav coefficients
! Function definition
real(dk)  ::  POT_NSG
! Locals
real(dk)  ::  sinph,alpha,lam,ERA
real(dk)  ::  MJD_UT1
real(dk)  ::  Plm(0:gdeg,0:gord)
real(dk)  ::  CLAM(0:gord),SLAM(0:gord)
integer   ::  l,m,mmax
real(dk)  ::  insum,outsum

! ==============================================================================

!!! TODO TODO TODO
!!! Efficiency improvement: compute sinph, tanph, cos(lam), sin(lam) only ONCE,
!!! possibly avoiding unnecessary evaluations of trigonometric functions.
!!! Also it is possible to avoid evaluating the Legendre functions twice.
!!! TODO TODO TODO

! ==============================================================================
! 01. Latitude and longitude in the planet-fixed frame
! ==============================================================================

! Latitude
sinph = r(3)/rm

! Longitude
alpha = atan2(r(2),r(1))
! Current GMST
MJD_UT1 = MJD0 + t/TU/secsPerDay
ERA = GMST_UNIFORM(MJD_UT1)      !!! TODO TODO TODO: STELA modifications.
! ERA = ERA_IERS10(MJD_UT1)
lam = alpha - ERA

! ==============================================================================
! 02. Evaluation of perturbing potential
! ==============================================================================

! Associated Legendre functions
Plm = 0._dk
Plm = ALF(gdeg,gord,sinph)
! Compute trig functions of lambda through recurrence relations
SLAM(0) = 0._dk; CLAM(0) = 1._dk
if (gord > 0 ) then
	SLAM(1) = sin(lam)
	CLAM(1) = cos(lam)
end if
do m=2,gord
	SLAM(m) = 2._dk*CLAM(1)*SLAM(m-1) - SLAM(m-2)
	CLAM(m) = 2._dk*CLAM(1)*CLAM(m-1) - CLAM(m-2)
end do

! Non-spherical gravity potential function.
outsum = 0._dk
do l = gdeg,2,-1
	insum = 0._dk
	mmax = minval([l,gord])
	do m=mmax,0,-1
		insum = insum + Plm(l,m)*(Clm(l,m)*CLAM(m) + Slm(l,m)*SLAM(m))
	end do
	outsum = outsum + (RE/rm)**l*insum
end do
POT_NSG = -(GM/rm)*outsum

end function POT_NSG




function ACC_NSG(r,rmag,dU)
! Compute the perturbing acceleration due to the non-spherical gravity field
! in the inertial frame given the position vector "r" and the partials of the
! potential in spherical coordinates "dU".

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  r(1:3)	  ! Position vector
real(dk),intent(in)  ::  rmag	  	! Position vector magnitude
real(dk),intent(in)  ::  dU(1:3)	! Partials of the potential in spherical coordinates
! Function def
real(dk)  ::  ACC_NSG(1:3)

! Locals
real(dk)  ::  xr,yr,zr,rxy,rxyrt
real(dk)  ::  dUdr,dUdp,dUdl

! ==============================================================================

! PARTIALS OF THE POTENTIAL:
! dU(1) = dU/dr
! dU(2) = dU/d(phi)
! du(3) = dU/d(lam)

! ACCELERATION COMPONENTS:
! ACC_NSG(1) = a_x in the inertial or Earth-fixed reference frame
! ACC_NSG(2) = a_y in the inertial or Earth-fixed reference frame
! ACC_NSG(3) = a_z in the inertial or Earth-fixed reference frame

! Unpack:
dUdr = dU(1); dUdp = dU(2); dUdl = dU(3)
! Direction cosines
xr = r(1)/rmag; yr = r(2)/rmag; zr = r(3)/rmag
rxy = r(1)**2 + r(2)**2
rxyrt = sqrt(rxy)

! Acceleration
ACC_NSG(1) = xr*(dUdr - zr/rxyrt*dUdp) - r(2)/rxy*dUdl
ACC_NSG(2) = yr*(dUdr - zr/rxyrt*dUdp) + r(1)/rxy*dUdl
ACC_NSG(3) = zr*dUdr + (rxyrt/rmag**2)*dUdp

end function ACC_NSG



function POTPAR(GM,RE,Clm,Slm,r,rm,t)
! Description:
!    Compute the partials of the non-spherical perturbing gravitational
!    potential in spherical coordinates.
!
! References:
!    Vallado D. A., Fundamentals of Astrodynamics 4th Ed., 2013.
! ==============================================================================

! MODULES
use AUXILIARIES, only: MJD0,TU
use PHYS_CONST,  only: secsPerDay
use PHYS_CONST,  only: GMST_UNIFORM,GRHOAN
! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)  ::  GM			               ! Gravitational parameter
real(dk),intent(in)  ::  RE			               ! Equatorial radius
real(dk),intent(in)  ::  r(1:3)		               ! Position vector in inertial frame
real(dk),intent(in)  ::  rm		                   ! Position vector magnitude
real(dk),intent(in)  ::  t                         ! Current time
real(qk),intent(in)  ::  Clm(1:maxDeg,0:maxOrd)    ! Grav coefficients
real(qk),intent(in)  ::  Slm(1:maxDeg,0:maxOrd)    ! Grav coefficients
! Function definition
real(dk)	::	POTPAR(1:3)

! LOCALS
! Current time
real(dk)  ::  MJD_UT1
! Coordinates in the Earth-fixed frame
real(dk)  ::  sinph,phi,alpha,lam,ERA
! Trig and ALFs
real(dk)  ::  PLexp(0:gdeg,0:gord+1)    ! Associated Legendre functions matrix (expanded)
! Indices
integer   ::  l,m,mmax
! Trigonometric functions
real(dk)  ::  slam(0:gord),CLAM(0:gord)
real(dk)  ::  tanph
! Auxiliaries
real(dk)  ::  mur,outsum,insum,REr

! ==============================================================================

! PARTIALS OF THE POTENTIAL:
! POTPAR(1) = dU/dr, partial wrt radius
! POTPAR(2) = dU/dphi, partial wrt geographic latitude
! POTPAR(3) = dU/dlam, partial wrt geographic longitude

!!! TODO TODO TODO
!!! The spherical harmonics have to be normalized to be used properly at high
!!! orders without round-off problems.
!!! TODO TODO TODO

! ==============================================================================
! 01. Latitude and longitude in the planet-fixed frame
! ==============================================================================

! Latitude
sinph = r(3)/rm
phi = asin(sinph)
PLexp = ALF(gdeg,gord+1,sinph)

! Right ascension
alpha = atan2(r(2),r(1))

! Current GMST
MJD_UT1 = MJD0 + t/TU/secsPerDay
ERA = GMST_UNIFORM(MJD_UT1)
! ERA = GRHOAN(MJD_UT1+2400000.5_dk)
lam = alpha - ERA

! ==============================================================================
! 02. Trigonometric functions
! ==============================================================================

! Compute trigonometric functions recursively
slam(0) = 0._dk;    clam(0) = 1._dk
if (gord > 0 ) then
	! This is the only evaluation of trigonometric functions of lambda ;)
	slam(1) = sin(lam)
	clam(1) = cos(lam)

end if

do m=2,gord
	! Compute sin(m*lambda), cos(m*lambda)
	slam(m) = 2._dk*clam(1)*slam(m-1) - slam(m-2)
	clam(m) = 2._dk*clam(1)*clam(m-1) - clam(m-2)

end do
tanph = tan(phi)
mur  = GM/rm
REr  = RE/rm

POTPAR = 0._dk

! ==============================================================================
! 03. dU/dr
! ==============================================================================

outsum = 0._dk
do l=gdeg,2,-1
	insum = 0._dk
	mmax = minval([l,gord])
	do m=mmax,0,-1
		insum = insum + PLexp(l,m)*(Clm(l,m)*clam(m) + Slm(l,m)*slam(m))
	end do
	outsum = outsum + (REr**l)*(l+1)*insum
end do
POTPAR(1) = -(mur/rm)*outsum

! ==============================================================================
! 04. dU/dphi
! ==============================================================================

outsum = 0._dk
do l=gdeg,2,-1
	insum = 0._dk
	mmax = minval([l,gord])
	do m=mmax,0,-1
		insum = insum + &
		(PLexp(l,m+1) - m*tanph*PLexp(l,m))*(Clm(l,m)*clam(m) + Slm(l,m)*slam(m))
	end do
	outsum = outsum + (REr**l)*insum
end do
POTPAR(2) = mur*outsum

! ==============================================================================
! 05. dU/dlam
! ==============================================================================

outsum = 0._dk
do l=gdeg,2,-1
	insum = 0._dk
	mmax = minval([l,gord])
	do m=mmax,0,-1
		insum = insum + m*PLexp(l,m)*(Slm(l,m)*clam(m) - Clm(l,m)*slam(m))
	end do
	outsum = outsum + (REr**l)*insum
end do
POTPAR(3) = mur*outsum

end function POTPAR



function ALF(deg,ord,x)
! Compute the values of the associated Legendre functions for argument x, for
! degree = 0, 1, ..., deg and order = 0, 1, ...,ord <= deg.

! VARIABLES
implicit none
! Arguments
integer    ::  deg,ord		! Degree and order
real(dk)   ::  x			! Argument
! Function def
real(dk)   ::  ALF(0:deg,0:ord)
real(dk)   ::  xrt
integer    ::  l,m,k

! ==============================================================================
ALF = 0._dk
xrt = sqrt(1._dk - x**2)

if (ord > 0) then
	! Initialize top left corner
	ALF(0:1,0:1) = reshape([1._dk,x,0._dk,xrt],[2,2])
else
	! ALF is a column array
	ALF(0:1,0) = [1._dk,x]
end if

! First column (Legendre polynomials)
do l=2,deg
	ALF(l,0) = ((2._dk*l - 1._dk)*x*ALF(l-1,0) - (l - 1._dk)*ALF(l-2,0))/l
end do

if (ord > 0) then
	! Diagonal
	do k=2,minval([ord,deg])
		ALF(k,k) = (2._dk*k - 1._dk)*xrt*ALF(k-1,k-1)
	end do

	! Sub-diagonal elements
	do m=1,ord
		do l=m+1,deg
			ALF(l,m) = ALF(l-2,m) + (2._dk*l - 1._dk)*xrt*ALF(l-1,m-1)
		end do
	end do
end if


end function ALF

!
! function DDGRHOAN(JD)
! ! Compute the angular rate of the Greenwich Hour Angle, i.e. the Earth's rotational
! ! velocity, using the derivative of the expression by Newcomb (1898).
! ! Output units are rad/s.
!
! ! VARIABLES
! use GLOBALS, only: pi
! implicit none
! ! Arguments
! real(dk),intent(in)  ::  JD        ! Julian Date
! real(dk)			 ::  DDGRHOAN
! ! Locals
! real(dk),parameter   ::  JD1900 = 2415020.0_dk  ! JD of Jan 0.5, 1900.
! real(dk)			 ::  T1900
!
! ! ==============================================================================
!
! T1900 = (JD - JD1900)/36525._dk
!
! DDGRHOAN = (pi/(12._dk*3600._dk))*(8640184.542_dk/(36525._dk*86400._dk) + &
! & 0.1858_dk/(36525._dk*86400._dk)*T1900 + 1._dk)
!
! end function DDGRHOAN

end module NSGRAV
