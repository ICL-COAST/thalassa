module INITIALIZE
! Description:
!    Contains wrapper initialization procedures for Thalassa.
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


subroutine INIT_STATE(eqs,R0,V0,t0,neq,y0,x0,mu)
! Description:
!    Initializes the state vector "y0" and the independent variable "x0" from
!    Cartesian coordinates and initial epoch.
!
! ==============================================================================

! MODULES
use AUXILIARIES,   only: DU,TU
use PHYS_CONST,    only: GE,RE,Clm,Slm
use SETTINGS,      only: insgrav
use PERTURBATIONS, only: PPOTENTIAL
use EDROMO,        only: EDROMO_PHI0,CART2EDROMO
use KUST_STI,      only: CART2KS

! VARIABLES
implicit none
! Arguments
integer,intent(in)               ::  eqs
real(dk),intent(in)              ::  R0(1:3),V0(1:3)
real(dk),intent(in)              ::  t0     ! Initial time [s]
real(dk),intent(in)              ::  mu     ! Grav parameter [km^3/s^2]
integer,intent(out)              ::  neq
real(dk),intent(out),allocatable ::  y0(:)
real(dk),intent(out)             ::  x0
! Locals
real(dk)                         ::  Rm
real(dk)                         ::  U0     ! Perturbing potential [km^2/s^2]
integer                          ::  ftime  ! Flag for time-like variable

! ==============================================================================

! Deallocate x0 if already allocated for some reason
if (allocated(y0)) deallocate(y0)

select case (eqs)

    case(1)   ! Cowell, 1st order
        neq = 6
        allocate(y0(1:neq))
        y0(1:3) = R0/DU
        y0(4:6) = V0/(DU*TU)
        x0 = t0*TU

    case(2:4) ! EDromo
        neq = 8
        ftime = eqs - 2
        allocate(y0(1:neq))
        Rm = sqrt(dot_product(R0,R0))
        U0 = PPOTENTIAL(insgrav,GE,RE,R0,Rm,t0*TU)  ! TBD
        x0 = EDROMO_PHI0(R0,V0,U0,mu,DU,TU)

        call CART2EDROMO(R0,V0,t0,DU,TU,y0,x0,U0,ftime)
    
    case(5:6) ! KS
        neq = 10
        ftime = eqs - 5
        allocate(y0(1:neq))
        Rm = sqrt(dot_product(R0,R0))
        U0 = PPOTENTIAL(insgrav,mu,RE,R0,Rm,t0*TU)  ! TBD
        call CART2KS(R0,V0,t0*TU,mu,DU,TU,y0,U0)    ! Input time to CART2KS is dimensionless
        x0 = 0._dk

end select

end subroutine INIT_STATE


end module INITIALIZE
