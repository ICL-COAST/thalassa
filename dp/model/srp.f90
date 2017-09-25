module SRP
! Description:
!    Contains subroutines to compute the perturbing acceleration due to SRP.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
! 
! ==============================================================================

! MODULES
use KINDS, only: dk
implicit none




contains




function SRP_ACC(pSRP_1au,au,CR,A2M_SRP,r,r_sun)
! Description:
!    Computes the perturbing acceleration due to solar radiation pressure.
!    Input units have to be consistent.
! 
! ==============================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)  ::  pSRP_1au,au,CR,A2M_SRP,r(1:3),r_sun(1:3)
! Function definition
real(dk)  ::  SRP_ACC(1:3)

! Locals
real(dk)  :: pSRP
real(dk)  :: dSun(1:3),dSunNorm

! ==============================================================================

! Sun->Spacecraft vector
dSun = r - r_sun
dSunNorm = sqrt(dot_product(dSun,dSun))

! SRP constant at the current distance from the Sun
pSRP = pSRP_1au*(au/dSunNorm)**2

! Perturbing acceleration
SRP_ACC = pSRP*CR*A2M_SRP*(dSun/dSunNorm)

end function SRP_ACC

end module SRP