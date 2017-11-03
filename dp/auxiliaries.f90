module AUXILIARIES
! Description:
!    Contains auxiliary variables to be communicated between different
!    procedures.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================

! VARIABLES
use KINDS, only: dk
implicit none
! Reference quantities for non-dimensionalization
real(dk)  ::  DU,TU
! Initial epoch, next epoch during integration, and final epoch (MJD-UT1)
real(dk)  ::  MJD0,MJDnext,MJDf
! Current coordinate frame and center
character(len=12)  :: coordSyst



contains




subroutine SET_UNITS(R)
! Description:
!    Contains reference units for non-dimensionalization. DU is a unit of
!    distance, while TU is a unit of frequency. With this choice of reference
!    units, GM = 1 (DU^3)*(TU^2).

! MODULES
use PHYS_CONST, only: CURRENT_MU
implicit none
real(dk),intent(in)  :: R(1:3)
real(dk)             :: mu

! ==============================================================================

mu = CURRENT_MU(coordSyst)

DU = sqrt(dot_product(R,R))
TU = sqrt(mu/DU**3)

end subroutine SET_UNITS



function T2MJD(t)
! Description:
!    Converts dimensionless time to MJD.

! ==============================================================================

! MODULES
use PHYS_CONST, only: secsPerDay

! VARIABLES
implicit none
real(dk),intent(in)  ::  t
real(dk)  ::  T2MJD

! ==============================================================================

T2MJD = t/TU/secsPerDay + MJD0

end function T2MJD




end module AUXILIARIES
