module COORD_SYST
! Description:
!    Contains procedures to switch between different coordinate systems.
! 
! Author:
!    Davide Amato
!    d.amato@upm.es
! 
! ==============================================================================

use KINDS, only: dk
implicit none

! Quadruple kind
integer,parameter  ::  qk = selected_real_kind(33)



contains




subroutine SWITCH_CS(coordSyst,t,R_from,V_from,R_to,V_to)
! Description:
!    Changes the coordinate system of the position and velocity 'R_from',
!    'V_from', giving 'R_to', 'V_to' as output. If 'coordSyst == ICRF' it
!    switches to MEIAUE, and vice versa.
!    The operation is carried out in quadruple precision to avoid issues due to
!    the change of coordinate origins.
! 
! ==============================================================================

use SUN_MOON, only: EPHEM_ICRF
implicit none
! VARIABLES
character(len=12),intent(inout)  ::  coordSyst
real(dk),intent(in)   ::  t                        ! Time [-]
real(dk),intent(in)   ::  R_from(1:3),V_from(1:3)  ! Original pos and vel [km,km/s]
real(dk),intent(out)  ::  R_to(1:3),V_to(1:3)      ! New pos and vel [km,km/s]

! Locals
real(dk)  ::  rMoon(1:3),vMoon(1:3)  ! Moon position and velocity [km,km/s]

! ==============================================================================

! Call Moon ephemerides (dimensional)
call EPHEM_ICRF(2,1._dk,1._dk,t,rMoon,vMoon)

select case (trim(coordSyst))
  case ('ICRF')
    ! Switch to MMEIAUE.
    R_to = real(R_from,qk) - real(rMoon,qk)
    V_to = real(V_from,qk) - real(vMoon,qk)
    coordSyst = 'MMEIAUE'
  
  case ('MMEIAUE')
    ! Switch to ICRF.
    R_to = real(R_from,qk) + real(rMoon,qk)
    V_to = real(V_from,qk) + real(vMoon,qk)
    coordSyst = 'ICRF'
  
end select

end subroutine SWITCH_CS

end module COORD_SYST