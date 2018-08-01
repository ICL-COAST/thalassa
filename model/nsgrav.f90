module NSGRAV

! MODULES
use KINDS,      only: dk
use SETTINGS,   only: gdeg,gord
use IO,         only: id_earth
use PHYS_CONST, only: qk,GE,RE,flatt,omegaE,secsPerDay,secsPerSidDay,twopi,&
&ERR_constant
implicit none

! VARIABLES
! Spherical harmonics (unnormalized) 
integer               ::  maxDeg, maxOrd
real(qk),allocatable  ::  Cnm(:,:),Snm(:,:)

! Pines algorithm matrices
real(dk),allocatable  ::  Anm(:,:),Dnm(:,:),Enm(:,:),Fnm(:,:),Rm(:),Im(:),Pn(:)
real(dk),allocatable  ::  Aux1(:),Aux2(:),Aux3(:),Aux4(:,:)

contains



subroutine INITIALIZE_NSGRAV(earthFile)
! Description:
!    Reads Earth gravity data from a text file. Initializes the gravitational
!    parameter, equatorial radius, flattening, rotational velocity, spherical
!    harmonics coefficients, and auxiliary matrices for Pines' algorithm. The
!    latter computes the Earth potential in Cartesian coordinates, avoiding
!    singularities due to the spherical harmonics.
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Revisions:
!    180730: Subroutine created from part of READ_PHYS().
! 
! ==============================================================================

! Arguments
character(len=*),intent(in)  ::  earthFile

! Locals
character(len=72)  ::  dummy
real(dk)           ::  invFlatt
integer            ::  i,j,l,m,n


! ==============================================================================

open(unit=id_earth,file=trim(earthFile),status='old',action='read')
read(id_earth,'(a)') (dummy, i=1,4)
read(id_earth,'(a37,i3)') dummy, maxDeg
read(id_earth,'(a37,i3)') dummy, maxOrd
read(id_earth,'(a36,e22.15)') dummy, GE
read(id_earth,'(a36,e22.15)') dummy, RE
read(id_earth,'(a36,e22.15)') dummy, invFlatt
read(id_earth,'(a36,e22.15)') dummy, omegaE

flatt    = 1._dk/invFlatt

! Read and de-normalize spherical harmonics coefficients and auxiliary matrices
! for Pines' algorithm.
! Initialize
l = 1; m = 0;
allocate(Cnm(1:maxDeg,0:maxDeg)); Cnm = 0._dk
allocate(Snm(1:maxDeg,0:maxDeg)); Snm = 0._dk
allocate(Anm(0:maxDeg+2, 0:maxDeg+2))
allocate(Dnm(1:maxDeg,   0:maxDeg))
allocate(Enm(1:maxDeg,   1:maxDeg))
allocate(Fnm(1:maxDeg,   1:maxDeg))
allocate(Rm(0:maxDeg)    )
allocate(Im(0:maxDeg)    )
allocate(Pn(0:maxDeg + 1))
allocate(Aux1(1:maxDeg+1))
allocate(Aux2(1:maxDeg+1))
allocate(Aux3(1:maxDeg+1))
allocate(Aux4(1:maxDeg+1, 0:maxDeg+1))

read(id_earth,'(a)') (dummy, i=1,2)
do i=1,maxDeg
  do j=0,minval([i,maxOrd])
    read(id_earth,'(2(1x,i2),2(1x,e24.17))') l,m,Cnm(i,j),Snm(i,j)
    Cnm(i,j) = Cnm(i,j)/NORMFACT(i,j)
    Snm(i,j) = Snm(i,j)/NORMFACT(i,j)
  end do
end do

close(id_earth)

! Fill coefficient arrays for Pines algorithm
Anm(:,:) = 0._dk
Dnm(:,:) = 0._dk
Enm(:,:) = 0._dk
Fnm(:,:) = 0._dk
Anm(0,0) = 1._dk
Anm(1,1) = 1._dk
do n = 1,maxDeg + 1
  Aux1(n) = 2._dk*n + 1._dk
  Aux2(n) = Aux1(n) / (n+1._dk)
  Aux3(n) = n / (n+1._dk)
  do m = 0, n-1
    Aux4(n,m) = (n+m+1._dk)
  end do
end do

! Earth Rotation Rate (revolutions per tropical day)
secsPerSidDay = twopi/omegaE
ERR_constant = secsPerDay/secsPerSidDay

end subroutine INITIALIZE_NSGRAV




function NORMFACT(l,m)
! Description:
!    Normalization factor for the spherical harmonics coefficients and
!    associated Legendre functions:
! 
!    sqrt( (l + m)! / ( (2 - delta_{0,m}) * (2n  + 1) * (n - m)! ) )
! 
! Reference:
!    [1] Montenbruck, O., Gill, E., "Satellite Orbits", p. 58, Springer, 2000.
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! ==============================================================================

! Arguments and function definition
integer,intent(in)  ::  l,m
real(qk)            ::  NORMFACT
! Locals
real(qk)            ::  lr,mr
real(qk)            ::  kron
real(qk)            ::  numer,denom

! ==============================================================================
lr = real(l,qk)
mr = real(m,qk)

numer = gamma(lr + mr + 1._qk)

if (m == 0) then
	kron = 1._qk
else
	kron = 0._qk
end if

denom = (2._qk - kron) * (2._qk*lr + 1._qk) * gamma(lr - mr + 1._qk)

NORMFACT = sqrt(numer/denom)

end function NORMFACT



subroutine PINES_NSG(GM,RE,rIn,F,pot,dPot)
! Description:
!    Compute the perturbing acceleration, perturbing potential (optional), and
!    time derivative of the perturbing potential in the body-fixed frame
!    (optional), given the position vector wrt the non-spherical body, its
!    gravitational parameter and its potential coefficients.
!    Uses the method described in the Reference to perform the calculation in
!    Cartesian coordinates, thereby avoiding associated Legendre functions and
!    polynomials. The formulation is non-singular everywhere for r > RE, where
!    RE is the radius of the non-spherical body.
! 
!    Adapted from code developed by Ho dei Urrutxua (Universidad Rey Juan Carlos,
!    Madrid, Spain) and Claudio Bombardelli (Universidad Politécnica de Madrid,
!    Madrid, Spain).
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Reference:
!    S. Pines, "Uniform Representation of the Gravitational Potential and its
!    derivatives," AIAA J. 11 (11), pp. 1508-1511, 1973.
! 
! Revisions:
!
! ==============================================================================

! Use associations
use AUXILIARIES, only: MJD0,TU
use PHYS_CONST,  only: GMST_UNIFORM

! Arguments
real(dk),intent(in)   ::  GM             ! Gravitational parameter
real(dk),intent(in)   ::  RE             ! Equatorial radius
real(dk),intent(in)   ::  rIn(1:3)       ! Position in the inertial frame
real(dk),intent(out),optional  ::  F(1:3)         ! Perturbing acceleration
real(dk),intent(out),optional  ::  pot   ! Perturbing potential
real(dk),intent(out),optional  ::  dPot  ! Time derivative of the potential in body-fixed frame
! Locals
real(dk)  ::  rNorm,rNormSq  ! Norm of position vector and square
real(dk)  ::  rEqNorm ! Norm of projection of position vector on equatorial frame
real(dk)  ::  cosRA,sinRA ! Direction cosines in the inertial frame
real(dk)  ::  s,t,u  ! Direction cosines in the body-fixed frame
real(dk)  ::  MJD_TT ! MJD in Terrestrial Time
real(dk)  ::  ERA,cosERA,sinERA    ! Earth Rotation Angle and its trig functions
real(dk)  ::  rho    ! = equatorial radius / r 
real(dk)  ::  a1,a2,a3,a4  ! Acceleration components 
integer   ::  n,m    ! Harmonic indices

! ==============================================================================

rNormSq = dot_product(rIn,rIn)
rNorm = sqrt(rNormSq)
rEqNorm = sqrt(rNormSq - rIn(3)**2)

! ==============================================================================
! 01. Transform from inertial to body-fixed frame
! ==============================================================================

! z_Inertial = z_Body (main body rotates along this axis)
u = rIn(3)/rNorm

! Direction cosines in the inertial frame
cosRA = rIn(1)/rEqNorm
sinRA = rIn(2)/rEqNorm

! Get Earth Rotation Angle 
MJD_TT = MJD0 + t/TU/secsPerDay ! NOTE: This will have to be modified to compute UT1 from TT in a next version.
ERA = GMST_UNIFORM(MJD_TT)
cosERA = cos(ERA); sinERA = sin(ERA)

! Direction cosines in the body-fixed frame
s = cosRA * cosERA + sinRA * sinERA
t = sinRA * cosERA - cosRA * sinERA

rho = RE/rNorm

! ==============================================================================
! 02. Fill in coefficient matrices and auxiliary vectors
! ==============================================================================

! Fill A Matrix
! TODO: Reduce the number of calculations if gord < gdeg.
Anm(0,0) = 1._dk
Anm(1,1) = 1._dk
Anm(1,0) = u
do n = 1, gdeg + 1
  Anm(n+1,n+1) = Aux1(n) * Anm(n,n) ! Fill the diagonal
  Anm(n+1,0) = Aux2(n) * u * Anm(n,0) - Aux3(n) * Anm(n-1,0) ! Fill the 1st column
  Anm(n+1,n) = u * Anm(n+1,n+1)    ! Fill the subdiagonal

end do
do n = 2, gdeg + 1 ! Fill remaining elements
  do m = 0, n - 2
    Anm(n+1,m+1) = Aux4(n,m) * Anm(n,m) + u * Anm(n,m+1)
  end do

end do

! Fill R, I, and P vectors
Rm(0) = 1._dk
Im(0) = 0._dk
Pn(0) = GM / rNorm
Pn(1) = rho * Pn(0)
do n = 1, gdeg
  Rm(n)   = s  * Rm(n-1) - t * Im(n-1)
  Im(n)   = s  * Im(n-1) + t * Rm(n-1)
  Pn(n+1) = rho * Pn(n)

end do

! Fill D, E, and F matrices
! TODO: add auxiliary matrix for the time derivative of pot.
! TODO: avoid computing E, F when the acceleration is not requested.
do m = 1, gord
  do n = m, gord
    Dnm(n,m) = Cnm(n,m)*Rm(m)   + Snm(n,m)*Im(m)
    Enm(n,m) = Cnm(n,m)*Rm(m-1) + Snm(n,m)*Im(m-1)
    Fnm(n,m) = Snm(n,m)*Rm(m-1) - Cnm(n,m)*Im(m-1)

  end do

end do
do n = 1, gdeg
  Dnm(n,0) = Cnm(n,0)*Rm(0)  !+ S(n,0)*I(0) = 0

end do

! ==============================================================================
! 03. Calculate perturbing function
! ==============================================================================

if (present(pot)) then
  pot = 0._dk
  do m = 1, gord
    do n = m, gord
      pot = pot + Pn(n) * Anm(n,m) * Dnm(n,m)
		
    end do
	
  end do
  do n = 1, gdeg
	pot = pot + Pn(n) * Anm(n,0) * Dnm(n,0)
	
  end do
  ! pot = pot + Pn(0) -> avoid this since we are only computing the perturbing potential
  ! pot = - 1.E-6_dk * pot -> non-dimensionalization is not needed necessarily
end if

! ==============================================================================
! 04. Calculate perturbing acceleration
! ==============================================================================

if (present(F)) then
  a1 = 0._dk; a2 = 0._dk; a3 = 0._dk; a4 = 0._dk
  do m = 1, gord
    do n = m, gord
      a1 = a1 + Pn(n+1) * Anm(n,m) * m * Enm(n,m)
      a2 = a2 + Pn(n+1) * Anm(n,m) * m * Fnm(n,m)
      a3 = a3 + Pn(n+1) * Anm(n,m+1)   * Dnm(n,m)
      a4 = a4 - Pn(n+1) * Anm(n+1,m+1) * Dnm(n,m)
    end do
  end do
  do n = 1, gord
    a3 = a3 + Pn(n+1) * Anm(n,1)     * Dnm(n,0)
    a4 = a4 - Pn(n+1) * Anm(n+1,1)   * Dnm(n,0)

  end do
  F = [a1, a2, a3] + [s, t, u] * a4
  F = F / RE
  ! F = F - Pn(0) * [s, t, u] / rNorm -> should not add the Keplerian term
  ! F = F * 1.E3_dk -> should not dimensionalize
end if


end subroutine PINES_NSG




function POT_NSG(GM,RE,Cnm,Snm,r,rm,t)
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
real(qk),intent(in)  ::  Cnm(1:maxDeg,0:maxOrd)             ! Grav coefficients
real(qk),intent(in)  ::  Snm(1:maxDeg,0:maxOrd)             ! Grav coefficients
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
		insum = insum + Plm(l,m)*(Cnm(l,m)*CLAM(m) + Snm(l,m)*SLAM(m))
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



function POTPAR(GM,RE,Cnm,Snm,r,rm,t)
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
real(qk),intent(in)  ::  Cnm(1:maxDeg,0:maxOrd)    ! Grav coefficients
real(qk),intent(in)  ::  Snm(1:maxDeg,0:maxOrd)    ! Grav coefficients
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
		insum = insum + PLexp(l,m)*(Cnm(l,m)*clam(m) + Snm(l,m)*slam(m))
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
		(PLexp(l,m+1) - m*tanph*PLexp(l,m))*(Cnm(l,m)*clam(m) + Snm(l,m)*slam(m))
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
		insum = insum + m*PLexp(l,m)*(Snm(l,m)*clam(m) - Cnm(l,m)*slam(m))
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
