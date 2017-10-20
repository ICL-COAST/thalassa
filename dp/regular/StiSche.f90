module STISCHE
! Description:
!    Contains subroutines necessary for the Stiefel-Scheifele formulation.
!
! References:
! [1] Stiefel & Scheifele, "Linear and Regular Celestial Mechanics", 1971
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


! ==============================================================================
! 01. INTEGRATION PROCEDURES
! ==============================================================================


subroutine STISCHE_RHS(neq,phi,z,zdot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of the
!    EDromo formulation.
!
! ==============================================================================

! MODULES
use SETTINGS,      only: eqs,insgrav,isun,imoon,idrag
use PHYS_CONST,    only: GE_nd,RE_nd,ERR_constant_nd
use PERTURBATIONS, only: PPOTENTIAL,PACC_EJ2K

! VARIABLES
implicit none

! Arguments
integer,intent(in)    ::  neq              ! Number of equations
real(dk),intent(in)   ::  phi              ! StiSche independent variable
real(dk),intent(in)   ::  z(1:neq)         ! StiSche state vector, ND
real(dk),intent(out)  ::  zdot(1:neq)      ! RHS of EoM's, ND

! Auxiliary quantities
real(dk)  ::  sph2,cph2
real(dk)  ::  aux0,aux1,aux2,aux3,aux4,aux5,aux6
integer   ::  flag_time

! K-S parameters and  their derivatives
real(dk)  ::  u_vec(1:4)
real(dk)  ::  du_vec(1:4)

! State
real(dk)  ::  rV(1:3),vV(1:3),t
real(dk)  ::  x_vec(1:3),y_vec(1:3)
real(dk)  ::  rmag,vmag,vsq

! Perturbations
real(dk)  ::  Upot,dUdr(1:3),dUdu(1:4),dUdt  ! Perturbing potential and its derivatives
real(dk)  ::  gradU_sph(1:3)
real(dk)  ::  p(1:3),Lp_vec(1:4)             ! Non-conservative perturbing accelerations

! ==============================================================================

! INDEX OF STATE VECTORS

! z(1) = zeta1
! z(2) = zeta2
! z(3) = zeta3
! z(4) = zeta4
! z(5) = zeta5
! z(6) = zeta6
! z(7) = zeta7
! z(8) = zeta8
! z(9) = zeta9
! z(10) = zeta10 (physical time,
!                 constant time element,
!                 linear time element,
!                 depending on flag_time)


! Safety check
!if (any(z/=z)) then
!	write(*,*) 'NaNs detected in the RHS, stopping execution.'
!	stop
!end if

! ==============================================================================
! 01. AUXILIARY QUANTITIES (1)
! ==============================================================================

! Store trig functions
sph2 = sin(phi/2._dk)
cph2 = cos(phi/2._dk)

! ==============================================================================
! 02. POSITION IN INERTIAL FRAME, TIME
! ==============================================================================

! K-S parameters
u_vec = [z(2)*cph2 + z(6)*sph2,  &
         z(3)*cph2 + z(7)*sph2,  &
         z(4)*cph2 + z(8)*sph2,  &
         z(5)*cph2 + z(9)*sph2]
! Derivatives of u_vec wrt the independent variable
du_vec = .5_dk*[-z(2)*sph2 + z(6)*cph2,  &
                -z(3)*sph2 + z(7)*cph2,  &
                -z(4)*sph2 + z(8)*cph2,  &
                -z(5)*sph2 + z(9)*cph2]

! Position in inertial frame
rV = [u_vec(1)**2 - u_vec(2)**2 - u_vec(3)**2 + u_vec(4)**2,  &
      2._dk*(u_vec(1)*u_vec(2) - u_vec(3)*u_vec(4)),  &
      2._dk*(u_vec(1)*u_vec(3) + u_vec(2)*u_vec(4))]
rmag = u_vec(1)**2 + u_vec(2)**2 + u_vec(3)**2 + u_vec(4)**2

flag_time = eqs - 2

! Get time
if ( flag_time == 0 ) then
    ! Physical time
    t = z(10)
elseif  ( flag_time == 1 ) then
    ! Constant time element
    ! I risultati non sono migliori rispetto ad un elemento tempo lineare
elseif  ( flag_time == 2 ) then
    ! Linear time element
    t = z(10) - dot_product(u_vec,du_vec)/z(1)
end if

! ==============================================================================
! 03. PERTURBING POTENTIAL
! ==============================================================================

! Initialize
Upot = 0._dk; dUdt = 0._dk; dUdr = 0._dk
! Evaluate potential
Upot = PPOTENTIAL(insgrav,GE_nd,RE_nd,rV,rmag,t)
! Evaluate time and spatial derivatives (note that velocity is not needed here)
dUdr = PACC_EJ2K(insgrav,0,0,0,0,rV,vV,rmag,t,gradU_sph)
dUdt = gradU_sph(3)*ERR_constant_nd

! ==============================================================================
! 04. VELOCITY IN THE INERTIAL FRAME
! ==============================================================================

vV = 4._dk*z(1)/rmag*[u_vec(1)*du_vec(1) - u_vec(2)*du_vec(2) - u_vec(3)*du_vec(3) + u_vec(4)*du_vec(4),  &
                      u_vec(2)*du_vec(1) + u_vec(1)*du_vec(2) - u_vec(4)*du_vec(3) - u_vec(3)*du_vec(4),  &
                      u_vec(3)*du_vec(1) + u_vec(4)*du_vec(2) + u_vec(1)*du_vec(3) + u_vec(2)*du_vec(4)]
vsq = vV(1)**2 + vV(2)**2 + vV(3)**2
vmag = sqrt(vsq)

! ==============================================================================
! 05. PERTURBING ACCELERATIONS
! ==============================================================================

! Initializations
p = 0._dk; f = 0._dk
p = PACC_EJ2K(0,isun,imoon,idrag,0,rV,vV,rmag,t)

! ==============================================================================
! 06. COMPUTE AUXILIARY QUANTITIES (2)
! ==============================================================================

Lp_vec = [ u_vec(1)*p(1) + u_vec(2)*p(2) + u_vec(3)*p(3),  &
          -u_vec(2)*p(1) + u_vec(1)*p(2) + u_vec(4)*p(3),  &
          -u_vec(3)*p(1) - u_vec(4)*p(2) + u_vec(1)*p(3),  &
           u_vec(4)*p(1) - u_vec(3)*p(2) + u_vec(2)*p(3)]
dUdu = 2._dk*[ u_vec(1)*dUdr(1) + u_vec(2)*dUdr(2) + u_vec(3)*dUdr(3),  &
              -u_vec(2)*dUdr(1) + u_vec(1)*dUdr(2) + u_vec(4)*dUdr(3),  &
              -u_vec(3)*dUdr(1) - u_vec(4)*dUdr(2) + u_vec(1)*dUdr(3),  &
               u_vec(4)*dUdr(1) - u_vec(3)*dUdr(2) + u_vec(2)*dUdr(3)]
aux0 = 2._dk/z(1)*zdot(1)
aux1 = .5_dk/z(1)**2*(.5_dk*Upot*u_vec(1) + rmag/4._dk*(dUdu(1) - 2._dk*Lp_vec(1))) + aux0*du_vec(1)
aux2 = .5_dk/z(1)**2*(.5_dk*Upot*u_vec(2) + rmag/4._dk*(dUdu(2) - 2._dk*Lp_vec(2))) + aux0*du_vec(2)
aux3 = .5_dk/z(1)**2*(.5_dk*Upot*u_vec(3) + rmag/4._dk*(dUdu(3) - 2._dk*Lp_vec(3))) + aux0*du_vec(3)
aux4 = .5_dk/z(1)**2*(.5_dk*Upot*u_vec(4) + rmag/4._dk*(dUdu(4) - 2._dk*Lp_vec(4))) + aux0*du_vec(4)
aux5 = dot_product(u_vec,dUdu)
aux6 = dot_product(u_vec,Lp_vec)

! ==============================================================================
! 07. COMPUTE RIGHT-HAND SIDE
! ==============================================================================

zdot(1) = -rmag/(8._dk*z(1)**2)*dUdt - .5_dk/z(1)*dot_product(du_vec,Lp_vec)
zdot(2) = aux1*sph2
zdot(3) = aux2*sph2
zdot(4) = aux3*sph2 
zdot(5) = aux4*sph2
zdot(6) = -aux1*cph2
zdot(7) = -aux2*cph2
zdot(8) = -aux3*cph2 
zdot(9) = -aux4*cph2 

! Time / Time Element
if (flag_time == 0) then
    zdot(10) = rmag/(2._dk*z(1))
else if (flag_time == 1) then   ! Constant Time Element
    ! NO
else if (flag_time == 2) then   ! Linear Time Element
   GE = 1._dk
   zdot(10) = 1._dk/(8._dk*z(1)**3)*(GE - 2._dk*rmag*Upot - 0.5_dk*rmag*(aux5 - 2._dk*aux6)) -  &
              2._dk/z(1)**2*zdot(1)*aux
end if

end subroutine STISCHE_RHS


subroutine STISCHE_EVT(neq,phi,z,ng,roots)

! MODULES
use AUXILIARIES, only: MJD0,MJDnext,MJDf,DU,TU
use PHYS_CONST,  only: secsPerDay,RE,reentry_radius_nd
use SETTINGS,    only: eqs

! VARIABLES
implicit none
! Arguments IN
integer,intent(in)    ::  neq
integer,intent(in)    ::  ng
real(dk),intent(in)   ::  phi
real(dk),intent(in)   ::  z(1:neq)
! Arguments OUT
real(dk),intent(out)  ::  roots(1:ng)

! Locals
integer   ::  flag_time
real(dk)  ::  t  ! Current time [-]
real(dk)  ::  sph2,cph2,rmag
real(dk)  ::  u_vec(1:4)

! ==============================================================================

roots = 1._dk

! Get time
flag_time = eqs - 2
t = STISCHE_TE2TIME(z,phi,flag_time)

! ==============================================================================
! 01. Next timestep
! ==============================================================================

roots(1) = t - (MJDnext - MJD0)*secsPerDay*TU

! ==============================================================================
! 02. Stop integration
! ==============================================================================

roots(2) = t - (MJDf - MJD0)*secsPerDay*TU

! ==============================================================================
! 03. Re-entry
! ==============================================================================
sph2 = sin(phi/2._dk)
cph2 = cos(phi/2._dk)
u_vec = [z(2)*cph2 + z(6)*sph2,  &
         z(3)*cph2 + z(7)*sph2,  &
         z(4)*cph2 + z(8)*sph2,  &
         z(5)*cph2 + z(9)*sph2]
rmag = u_vec(1)**2 + u_vec(2)**2 + u_vec(3)**2 + u_vec(4)**2

roots(3) = rmag - reentry_radius_nd

end subroutine STISCHE_EVT


function STISCHE_PHI0(R,V,pot,GM,DU,TU)
! Description:
!    Computes initial value of the Stiefel-Scheifele fictitious time
!
! ==============================================================================

implicit none
! Arguments
real(dk),intent(in)  ::  R(1:3),V(1:3),pot,GM,DU,TU
real(dk)  ::  STISCHE_PHI0
! Locals
real(dk)  ::  R_nd(1:3),V_nd(1:3),GM_nd,Rmag_nd
real(dk)  ::  rvdot
real(dk)  ::  pot_nd,totEn

R_nd = R/DU; V_nd = V/(DU*TU); GM_nd = GM/(DU**3*TU**2); pot_nd = pot/(DU*TU)**2
Rmag_nd = sqrt(dot_product(R_nd,R_nd))
rvdot = dot_product(R_nd,V_nd)
totEn = .5_dk*dot_product(V_nd,V_nd) - GM_nd/Rmag_nd - pot_nd
STISCHE_PHI0 = atan2(rvdot*sqrt(-2._dk*totEn), 1._dk + 2._dk*totEn*Rmag_nd)

end function STISCHE_PHI0


! ==============================================================================
! 02. TRANSFORMATIONS AND PROCESSING PROCEDURES
! ==============================================================================


function STISCHE_TE2TIME(z,phi,flag_time)
! Description:
!    Gets the value of physical time from the Stiefel-Scheifele state vector "z" and
!    fictitious time "phi".
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)  ::  z(1:10),phi
integer,intent(in)   ::  flag_time
! Function definition
real(dk)  ::  STISCHE_TE2TIME

! ==============================================================================

if ( flag_time == 0 ) then
    ! Physical time
    STISCHE_TE2TIME = z(10)

else if  ( flag_time == 1 ) then
    ! Constant time element
    ! NO

else if  ( flag_time == 2 ) then
    ! Linear time element
    STISCHE_TE2TIME = z(10) - dot_product(u_vec,du_vec)/z(1)

end if

end function STISCHE_TE2TIME


subroutine CART2STISCHE(R,V,t0,DU,TU,z,phi,W,flag_time)
! Description:
!    Converts from Cartesian coordinates to Stiefel-Scheifele elements. It requires
!    the current values of the fictitious time, perturbing potential and initial
!    time.
!
!===============================================================================

! VARIABLES
implicit none

! Arguments
integer,parameter     ::  neq = 10	  ! Number of elements of Stiefel-Scheifele state vector
real(dk),intent(in)   ::  R(1:3),V(1:3)	  ! Dimensional position and velocity [km,km/s]
real(dk),intent(in)   ::  DU,TU           ! Ref. quantities for non-dimensionalization
real(dk),intent(in)   ::  phi             ! Initial phi value
real(dk),intent(in)   ::  t0              ! Initial time value [s]
real(dk),intent(in)   ::  W               ! Potential energy [km^2/s^2]
integer,intent(in)    ::  flag_time	  ! = 0 physical time
         				  ! = 1 constant time element
					  ! = 2 linear time element
real(dk),intent(out)  ::  z(1:neq)        ! Stiefel-Scheifele state vector

! Local variables
real(dk)  ::  y(1:6)           		        ! Cartesian state vector, ND
real(dk)  ::  Rmag      	                ! Radius and velocity magnitudes, ND
real(dk)  ::  pot0           		        ! Potential, ND
real(dk)  ::  totEn0           		        ! Initial total energy, ND
real(dk)  ::  cph2,sph2
real(dk)  ::  u_vec(1:4),du_vec(1:4)            ! K-S parameters and their derivatives
real(dk)  ::  GE                	        ! Auxiliary variable
real(dk)  ::  zero             		        ! Reference machine zero

! ==============================================================================

! ==============================================================================
! 01. NON-DIMENSIONALIZATION
! ==============================================================================
y(1:3) = R/DU
y(4:6) = V/(DU*TU)
pot0   = W/(DU*TU)**2

! ==============================================================================
! 02. AUXILIARY QUANTITIES
! ==============================================================================
! Compute machine zero for comparisons
zero = epsilon(0._dk)

Rmag = sqrt(dot_product(y(1:3),y(1:3)))
cph2 = cos(phi/2._dk)
sph2 = sin(phi/2._dk)

! ==============================================================================
! 03. COMPUTE z(1) - z(9)
! ==============================================================================
! Total energy
GE = 1._dk
totEn0  = .5_dk*Vmag**2 - GE/Rmag + pot0
z(1) = sqrt(-totEn0/2._dk)

! K-S parameters
if ( y(1).ge.0._dk ) then
   u_vec(1) = 0._dk
   u_vec(4) = sqrt(.5_dk*(Rmag + y(1)) - u_vec(1)**2)
   u_vec(2) = (y(2)*u_vec(1) + y(3)*u_vec(4))/(Rmag + y(1))
   u_vec(3) = (y(3)*u_vec(1) - y(2)*u_vec(4))/(Rmag + y(1))
else
   u_vec(2) = 0._dk
   u_vec(3) = sqrt(.5_dk*(Rmag - y(1)) - u_vec(2)**2)
   u_vec(1) = (y(2)*u_vec(2) + y(3)*u_vec(3))/(Rmag - y(1))
   u_vec(4) = (y(3)*u_vec(2) - y(2)*u_vec(3))/(Rmag - y(1))
end if
! Derivatives of the K-S parameters wrt the independent variable
du_vec(1) = ( u_vec(1)*y(4) + u_vec(2)*y(5) + u_vec(3)*y(6))/(4._dk*z(1))
du_vec(2) = (-u_vec(2)*y(4) + u_vec(1)*y(5) + u_vec(4)*y(6))/(4._dk*z(1))
du_vec(3) = (-u_vec(3)*y(4) - u_vec(4)*y(5) + u_vec(1)*y(6))/(4._dk*z(1))
du_vec(4) = ( u_vec(4)*y(4) - u_vec(3)*y(5) + u_vec(2)*y(6))/(4._dk*z(1))

z(2:5) = cph2*u_vec - 2._dk*sph2*du_vec
z(6:9) = sph2*u_vec + 2._dk*cph2*du_vec

! ==============================================================================
! 04. TIME / TIME ELEMENT
! ==============================================================================
if ( flag_time == 0 ) then
    ! Physical time
    z(10) = t0*TU
elseif  ( flag_time == 1 ) then
    ! Constant time element
    ! NO
elseif  ( flag_time == 2 ) then
    ! Linear time element
    z(10) = t0*TU + dot_product(u_vec,du_vec)/z(1)
end if

contains

end subroutine CART2STISCHE


subroutine STISCHE2CART(phi,z,r_vec,v_vec)
! Description:
!    Transforms from the Stiefel-Scheifele state vector "z", fictitious time "phi" and
!    potential "Upot" to Cartesian position and velocity "r_vec", "v_vec".
!    **All quantities are dimensionless, unlike in CART2STISCHE**.
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)   ::  z(1:10),phi
! Arguments OUT
real(dk),intent(out)  ::  r_vec(1:3),v_vec(1:3)

! Auxiliaries
real(dk)  ::  sph2,cph2
real(dk)  ::  rmag
real(dk)  ::  u_vec(1:4),du_vec(1:4)

! ==============================================================================

! ==============================================================================
! 01. AUXILIARY QUANTITIES
! ==============================================================================

! Store trig functions
sph2 = sin(phi/2._dk)
cph2 = cos(phi/2._dk)

! ==============================================================================
! 02. POSITION IN INERTIAL FRAME
! ==============================================================================

! K-S parameters
u_vec = [z(2)*cph2 + z(6)*sph2,  &
         z(3)*cph2 + z(7)*sph2,  &
         z(4)*cph2 + z(8)*sph2,  &
         z(5)*cph2 + z(9)*sph2]
! Derivatives of the K-S parameters wrt the independent variable
du_vec = .5_dk*[-z(2)*sph2 + z(6)*cph2,  &
                -z(3)*sph2 + z(7)*cph2,  &
                -z(4)*sph2 + z(8)*cph2,  &
                -z(5)*sph2 + z(9)*cph2]

! Position in inertial frame
r_vec = [u_vec(1)**2 - u_vec(2)**2 - u_vec(3)**2 + u_vec(4)**2,  &
         2._dk*(u_vec(1)*u_vec(2) - u_vec(3)*u_vec(4)),  &
         2._dk*(u_vec(1)*u_vec(3) + u_vec(2)*u_vec(4))]
rmag = u_vec(1)**2 + u_vec(2)**2 + u_vec(3)**2 + u_vec(4)**2

! ==============================================================================
! 03. VELOCITY IN INERTIAL FRAME
! ==============================================================================

v_vec = 4._dk*z(1)/rmag*[u_vec(1)*du_vec(1) - u_vec(2)*du_vec(2) - u_vec(3)*du_vec(3) + u_vec(4)*du_vec(4),  &
                         u_vec(2)*du_vec(1) + u_vec(1)*du_vec(2) - u_vec(4)*du_vec(3) - u_vec(3)*du_vec(4),  &
                         u_vec(3)*du_vec(1) + u_vec(4)*du_vec(2) + u_vec(1)*du_vec(3) + u_vec(2)*du_vec(4)]

end subroutine STISCHE2CART


end module STISCHE
