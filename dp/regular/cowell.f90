module COWELL
! Description:
!    Contains procedures necessary to integrate the Cowell formulation.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================

use KINDS, only: dk
implicit none


contains


subroutine COWELL_RHS(neq,t,y,ydot)
! Description:
!    Computes the value of the right-hand side of the 1st-order equations of
!    motion of the Cowell formulation.
!
! ==============================================================================

! MODULES
use AUXILIARIES,   only: DU,TU,coordSyst
use SETTINGS,      only: insgrav,isun,imoon,idrag,iSRP
use PERTURBATIONS, only: PACC_ICRF,PACC_MMEIAUE

! VARIABLES
implicit none

! Arguments
integer,intent(in)     ::  neq             ! Number of equations.
real(dk),intent(in)    ::  t               ! Time, ND.
real(dk),intent(in)    ::  y(1:neq)        ! Cartesian state vector, ND.
real(dk),intent(out)   ::  ydot(1:neq)     ! RHS of EoM's, ND.

! Local variables
real(dk)          ::  rMag                 ! Magnitude of position vector. [-]
real(dk)          ::  p(1:3)          ! Perturbation acceleration. [-]

! ==============================================================================

rMag = sqrt(dot_product(y(1:3),y(1:3)))

! ==============================================================================
! 01. COMPUTE PERTURBATIONS
! ==============================================================================

p = 0._dk
select case(trim(coordSyst))
  case('ICRF')
    p = PACC_ICRF(insgrav,isun,imoon,idrag,iSRP,y(1:3),y(4:6),rMag,t)
  
  case('MMEIAUE')
    p = PACC_MMEIAUE(insgrav,isun,imoon,iSRP,y(1:3),y(4:6),rMag,t)

end select

! ==============================================================================
! 02. EVALUATE RIGHT-HAND SIDE
! ==============================================================================

ydot(1:3) = y(4:6)
ydot(4:6) = -y(1:3)/rMag**3 + p

end subroutine COWELL_RHS


subroutine COWELL_EVT(neq,t,y,ng,roots)
! Description:
!    Finds roots to stop the integration for the Cowell formulation.
!    Quad precision.
!
! ==============================================================================

! MODULES
use AUXILIARIES, only: MJD0,MJDnext,MJDf,TU,DU,coordSyst,T2MJD
use PHYS_CONST,  only: secsPerDay,reentry_radius_nd,GE,GM
use SETTINGS,    only: iswitch,RSw_Hill
use SUN_MOON,    only: aMoon_Kep
use SUN_MOON,    only: EPHEM_ICRF
use COORD_SYST,  only: R_SOI,POS_VEL_ICRF

! VARIABLES
implicit none
! Arguments IN
integer,intent(in)       ::  neq
integer,intent(in)       ::  ng
real(dk),intent(in)      ::  t
real(dk),intent(in)      ::  y(1:neq)
! Arguments OUT
real(dk),intent(out)     ::  roots(1:ng)
! Locals
real(dk)  ::  MJD
real(dk)  ::  rmag
real(dk)  ::  RSw
real(dk)  ::  r_ICRF(1:3),v_ICRF(1:3)
real(dk)  ::  rMoon_ICRF(1:3),vMoon_ICRF(1:3)
real(dk)  ::  RHillMoon,d(1:3),dmag

! ==============================================================================

roots = 1._dk
! ==============================================================================
! 01. Next timestep (DISABLED)
! ==============================================================================

!roots(1) = t - MJDnext*secsPerDay*TU
roots(1) = 1._dk

! ==============================================================================
! 02. Stop integration
! ==============================================================================

roots(2) = t - (MJDf-MJD0)*secsPerDay*TU
!roots(2) = 1._dk

! ==============================================================================
! 03. Re-entry (temporarily disabled)
! ==============================================================================

! rmag = sqrt(dot_product(y(1:3),y(1:3)))
! roots(3) = rmag - reentry_radius_nd
roots(3) = 1._dk

! ==============================================================================
! 04. Distance-based switch of coordinate system
! ==============================================================================

if (iswitch == 1) then
  RHillMoon = aMoon_Kep*(GM/(3._dk*GE))**(1._dk/3._dk)
  RSw = RSw_Hill*RHillMoon

  ! Position of the particle in ICRF
  MJD = T2MJD(t)
  call POS_VEL_ICRF(coordSyst,MJD,DU,TU,y(1:3),y(4:6),r_ICRF,v_ICRF)
  call EPHEM_ICRF(2,1._dk,1._dk,MJD,rMoon_ICRF,vMoon_ICRF)
  
  ! Position of the particle wrt Moon
  d = r_ICRF - rMoon_ICRF
  dmag = sqrt(dot_product(d,d))
  roots(4) = dmag - RSw

else
  roots(4) = 1._dk
  
end if

end subroutine COWELL_EVT


end module COWELL
