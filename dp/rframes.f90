module RFRAMES
! Description:
!    Procedures to convert to the output reference frames.
!
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu

use KINDS, only: dk
implicit none

! SPICE routines
external SXFORM


contains



subroutine EQ2EC(npts,cart_EQ,cart_EC)
! Description:
!    Transforms a Cartesian trajectory array in the equatorial frame (EJ2K) to 
!    the ecliptic frame (ECJ2K).
! 
! ==============================================================================

implicit none
! Arguments IN/OUT
integer,intent(in)    ::  npts
real(dk),intent(in)   ::  cart_EQ(1:npts,1:7)
real(dk),intent(out)  ::  cart_EC(1:npts,1:7)

! Locals
integer   ::  i
real(dk)  ::  R_EQ2EC(6,6)

! ==============================================================================

call SXFORM('J2000','ECLIPJ2000',0._dk,R_EQ2EC)

do i=1,npts
  cart_EC(i,1) = cart_EQ(i,1)
  cart_EC(i,2:7) = matmul(R_EQ2EC,cart_EQ(i,2:7))

end do

end subroutine EQ2EC

end module RFRAMES