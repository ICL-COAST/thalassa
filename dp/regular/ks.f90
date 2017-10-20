module KUST_STI
! Description:
!    Contains the subroutines necessary for the Kustaanheimo-Stiefel
!    formulation.
! 
! References:
! [1] Stiefel, E. L. and Scheifele, G. "Linear and Regular Celestial Mechanics",
!     Springer-Verlag, 1971.
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

subroutine KS_RHS(neq,s,u,udot)
! Description:
!    Computes the value of the right-hand side of the equations of motion of the
!    Kustaanheimo-Stiefel formulation.
!
! ==============================================================================

! MODULES

! VARIABLES
implicit none
! Arguments
integer,intent(in)       ::  neq                  ! Number of equations.
real(dk),intent(in)      ::  s                    ! Value of fictitious time.
real(dk),intent(in)      ::  u(1:neq)             ! KS state vector.
real(dk),intent(out)     ::  udot(1:neq)          ! RHS of EoM's, ND.
! Local variables
! -- State variables
real(dk)    ::  x(1:4),xdot(1:4)                  ! Radius and velocity in R^4, ND
real(dk)    ::  r                                 ! Radius magnitude, ND
real(dk)    ::  t                                 ! Physical time, ND
! -- Perturbations
real(dk)    ::  rVpot               ! Perturbing potential term
real(dk)    ::  drVdu(1:4)          ! Derivatives of the perturbing potential term


! ==============================================================================

! STATE VECTOR DICTIONARY
! u(1:4)        u1,...,u4; KS-position, in R^4
! u(5:8)        u1',...,u4'; KS-velocity, in R^4
! u(9)          h; (-total energy) = (-Keplerian energy) + (-potential)
! u(10)         t; non-dimensional physical time

! ==============================================================================
! 01. COMPUTE CARTESIAN COORDINATES, ND
! ==============================================================================

call KS2CART(u,x,xdot)
r = sqrt(dot_product(x,x))
t = u(10)

! ==============================================================================
! 02. COMPUTE PERTURBING POTENTIAL, ND
! ==============================================================================

! Initialize
rVpot = 0._dk; drVdu = 0._dk

end subroutine KS_RHS


function KSMAT(u)
! Description:
!    Computes the KS-matrix from the K-S state vector.
! 
! ==============================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)  ::  u(:)                ! K-S state vector
real(dk)             ::  KSMAT(1:4,1:4)      ! K-S-matrix

! ==============================================================================

KSMAT = RESHAPE([ u(1),  u(2),  u(3),  u(4)&
                &,-u(2),  u(1),  u(4), -u(3)&
                &,-u(3), -u(4),  u(1),  u(2)&
                &, u(4), -u(3),  u(2), -u(1)],[4,4])

end function KSMAT 


subroutine KS2CART(u,x,xdot)
! Description:
!    Computes the Cartesian state (dimensionless) from the KS state vector.
! 
! ==============================================================================

! VARIABLES
implicit none
! Arguments
real(dk),intent(in)   ::  u(:)                  ! KS (extended) state vector
real(dk),intent(out)  ::  x(:),xdot(:)          ! Position and velocity, ND
! Locals
real(dk)              ::  r                     ! Radius magnitude, ND

! ==============================================================================

! Position
x(1) = u(1)**2 - u(2)**2 - u(3)**2 + u(4)**2
x(2) = 2._dk*(u(1)*u(2) - u(3)*u(4))
x(3) = 2._dk*(u(1)*u(3) + u(2)*u(4))
r    = u(1)**2 + u(2)**2 + u(3)**2 + u(4)**2

! Velocity
xdot(1) = 2._dk*(u(1)*u(5) - u(2)*u(6) - u(3)*u(7) + u(4)*u(8))/r
xdot(2) = 2._dk*(u(2)*u(5) + u(1)*u(6) - u(4)*u(7) - u(3)*u(8))/r
xdot(3) = 2._dk*(u(3)*u(5) + u(4)*u(6) + u(1)*u(7) + u(2)*u(8))/r

! If x, xdot are in R^4 their fourth component is null.
if (size(x,1) > 3)    x(4:) = 0._dk
if (size(xdot,1) > 3) xdot(4:) = 0._dk

end subroutine KS2CART


end module KUST_STI