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

use SUN_MOON,    only: EPHEM_ICRF
use AUXILIARIES, only: DU,TU,SET_UNITS
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
    call SET_UNITS(R_to)
  
  case ('MMEIAUE')
    ! Switch to ICRF.
    R_to = real(R_from,qk) + real(rMoon,qk)
    V_to = real(V_from,qk) + real(vMoon,qk)
    coordSyst = 'ICRF'
    call SET_UNITS(R_to)
    
end select

end subroutine SWITCH_CS




function R_SOI(a_sat,mu_main,mu_sat)
! Description:
!    Computes the radius of the sphere of influence of the body of semi-major
!    axis 'a_sat' and gravitational parameter 'mu_sat', orbiting around the
!    body of mass 'mu_main'.
!    Units have to be consistent.
! 
! ==============================================================================

! VARIABLES
implicit none
real(dk),intent(in)  ::  a_sat    ! Semi-major axis of the satellite
real(dk),intent(in)  ::  mu_sat   ! Gravitational parameter of the satellite
real(dk),intent(in)  ::  mu_main  ! Gravitational parameter of the main body
real(dk)             ::  R_SOI

! ==============================================================================

R_SOI = a_sat * (mu_sat / mu_main)**(2._dk/5._dk)

end function R_SOI


subroutine POS_VEL_ICRF(coordSyst,t,DU,TU,rCurr,vCurr,rICRF,vICRF)
! Description:
!    Computes the position and velocity in the ICRF coordinate system, starting
!    from position and velocity in the coordinate system specified by
!    'coordSyst'.
! 
! ==============================================================================

use SUN_MOON, only: EPHEM_ICRF

! VARIABLES
implicit none
character(len=12),intent(in)  ::  coordSyst
real(dk),intent(in)           ::  t,DU,TU,rCurr(1:3),vCurr(1:3)
real(dk),intent(out)          ::  rICRF(1:3),vICRF(1:3)

! LOCALS
real(dk)  ::  rMoon_ICRF(1:3),vMoon_ICRF(1:3)

! ==============================================================================

select case(trim(coordSyst))
  case ('ICRF')
    rICRF = rCurr
    vICRF = vCurr
  
  case ('MMEIAUE')
    call EPHEM_ICRF(2,DU,TU,t,rMoon_ICRF,vMoon_ICRF)
    rICRF = rMoon_ICRF*DU    + rCurr
    vICRF = vMoon_ICRF*DU*TU + vCurr

end select

end subroutine POS_VEL_ICRF




function ICRF2SYN(rICRF,vICRF,r2,v2,m2,m1,n)
! Description:
!     Transforms from ICRF coordinates rICRF,vICRF [km,km/s] to synodic
!     coordinates.
!
! ==============================================================================

! MODULES
use SUN_MOON, only: EPHEM_ICRF

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  rICRF(1:3),vICRF(1:3)    ! ICRF state [km,km/s].
real(dk),intent(in)  ::  r2(1:3),v2(1:3)          ! State of secondary body (Moon) [km,km/s].
real(dk),intent(in)  ::  m2,m1,n
! Function definition
real(dk)  ::  ICRF2SYN(1:6)             ! Synodic state [km,km/s]
! Locals
real(dk)  ::  CM(1:6)                   ! CM position and velocity [km,km/s]
real(dk)  ::  m2_red                    ! Reduced mass of secondary
real(dk)  ::  y_CM(1:6)                 ! Position and velocity wrt CM [km,km/s]
real(dk)  ::  theta                     ! Secondary's true longitude
real(dk)  ::  Rot(3,3),dRot(3,3)        ! Rotation matrices

! ==============================================================================

m2_red = m2/(m2+m1)
! Compute position and velocity of the CM
CM = m2_red*[r2,v2]

! Compute position and velocity WRT CM in helio coord.
y_CM(1:3) = rICRF - CM(1:3)
y_CM(4:6) = vICRF - CM(4:6)

! Secondary's true longitude theta
theta = atan2(r2(2),r2(1))

! Rotation matrices
Rot  = RESHAPE([cos(theta),-sin(theta),0._dk&   ! 1st col
             &,sin(theta),cos(theta),0._dk&     ! 2nd col
             &,0._dk,0._dk,1._dk],[3,3])        ! 3rd col

dRot = RESHAPE([-sin(theta),-cos(theta),0._dk&  ! 1st col
             &,cos(theta),-sin(theta),0._dk&    ! 2nd col
             &,0._dk,0._dk,0._dk],[3,3])        ! 3rd col

! Synodic position
ICRF2SYN(1:3) = MATMUL(Rot,y_CM(1:3))
! Synodic velocity
ICRF2SYN(4:6) = MATMUL(Rot,y_CM(4:6)) + n*MATMUL(dRot,y_CM(1:3))

end function ICRF2SYN



function SYN2ICRF(rSYN,vSYN,r2,v2,m2,m1,n)
! Description:
!     Transforms from synodic coordinates rSYN,vSYN [km,km/s] to synodic
!     coordinates.
!
! ==============================================================================

! MODULES
use SUN_MOON, only: EPHEM_ICRF

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  rSYN(1:3),vSYN(1:3)    ! ICRF state [km,km/s].
real(dk),intent(in)  ::  r2(1:3),v2(1:3)          ! State of secondary body (Moon) [km,km/s].
real(dk),intent(in)  ::  m2,m1,n
! Function definition
real(dk)  ::  SYN2ICRF(1:6)             ! Synodic state [km,km/s]
! Locals
real(dk)  ::  CM(1:6)                   ! CM position and velocity [km,km/s]
real(dk)  ::  m2_red                    ! Reduced mass of secondary
real(dk)  ::  y_CM(1:6)                 ! Position and velocity wrt CM [km,km/s]
real(dk)  ::  theta                     ! Secondary's true longitude
real(dk)  ::  Rot(3,3),dRot(3,3)        ! Rotation matrices

! ==============================================================================

! Secondary's true longitude theta
theta = atan2(r2(2),r2(1))

! Rotation matrices
Rot  = RESHAPE([cos(theta),sin(theta),0._dk&   ! 1st col
             &,-sin(theta),cos(theta),0._dk&     ! 2nd col
             &,0._dk,0._dk,1._dk],[3,3])        ! 3rd col

dRot = RESHAPE([-sin(theta),cos(theta),0._dk&  ! 1st col
             &,-cos(theta),-sin(theta),0._dk&    ! 2nd col
             &,0._dk,0._dk,0._dk],[3,3])        ! 3rd col

! Barycentric coordinates
y_CM(1:3)  =  MATMUL(Rot,rSYN)
y_CM(4:6)  =  MATMUL(Rot,vSyn) + n*MATMUL(dRot,rSyn)

! CM position and velocity in ICRF
m2_red = m2/(m2+m1)
CM = m2_red*[r2,v2]

! ICRF position and velocity
SYN2ICRF = y_CM + CM

end function SYN2ICRF




end module COORD_SYST