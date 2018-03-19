program THALASSA
! Description:
! Thalassa propagates Earth-bound orbits from an initial to a final epoch. It
! uses either Newtonian or regularized equations of motion written in the
! EMEJ2000 frame. It takes into account the following perturbations:
! - Non-spherical Earth gravitational potential
! - Third-body perturbations from Sun and Moon
! - Drag
! - SRP
!
! The choice of the mathematical model and constants reflect those used in
! STELA, a widely-known semi-analytical propagator currently in usage for studying
! the LEO, MEO and GEO regions.
! Also, the user can choose between the following numerical integrators:
! - LSODAR
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================


! MODULES
use KINDS,       only: dk
use AUXILIARIES, only: MJD0
use IO,          only: READ_IC,CREATE_OUT,DUMP_TRAJ
use CART_COE,    only: COE2CART,CART2COE
use PHYS_CONST,  only: READ_PHYS,GMST_UNIFORM
use PROPAGATE,   only: DPROP_REGULAR
use SETTINGS,    only: READ_SETTINGS,input_path
use IO,          only: id_cart,id_orb,id_stat,object_path
use SETTINGS,    only: gdeg,gord,outpath,tol,eqs
use PHYS_CONST,  only: GE,d2r,r2d,secsPerDay,secsPerSidDay,twopi
implicit none

! VARIABLES
! Initial conditions (EMEJ2000)
real(dk)  ::  COE0(1:6),COE0_rad(1:6)
real(dk)  ::  R0(1:3),V0(1:3)
real(dk)  ::  GMST0,aGEO
real(dk)  ::  period
! Integration span and dt [solar days]
real(dk)  ::  tspan,tstep
! Trajectory
integer               ::  npts,ipt
real(dk),allocatable  ::  cart(:,:),orb(:,:)
real(dk)              ::  R(1:3),V(1:3)
! Measurement of CPU time
integer  ::  rate,tic,toc
real(dk) ::  cputime
! Function calls and integration steps
integer  ::  int_steps,tot_calls
! Command arguments
integer  ::  command_arguments


! ==============================================================================

! ==============================================================================
! 01. COMMAND LINE PARSING
! ==============================================================================

command_arguments = COMMAND_ARGUMENT_COUNT()
input_path  = './in/input.txt'
object_path = './in/object.txt'
if (command_arguments > 0) call GET_COMMAND_ARGUMENT(1,input_path)
if (command_arguments > 1) call GET_COMMAND_ARGUMENT(2,object_path)

! ==============================================================================
! 02. INITIALIZATIONS
! ==============================================================================

! Start clock
call SYSTEM_CLOCK(tic,rate)

! Read initial conditions, settings and physical model data.
call READ_IC(MJD0,COE0)
call READ_SETTINGS(tspan,tstep)
call READ_PHYS()

! Initialize output (if not done by python script already)
if (command_arguments == 0) then
  call SYSTEM('mkdir -p '//trim(outpath))
  call SYSTEM('cp in/*.txt '//trim(outpath))
end if

! Load SPICE kernels
call FURNSH('./data/kernels_to_load.furnsh')

! ==============================================================================
! 03. TEST PROPAGATION
! ==============================================================================

! Convert to Cartesian coordinates
COE0_rad = [COE0(1:2),COE0(3:6)*real(d2r,dk)]
call COE2CART(COE0_rad,R0,V0,GE)

! Output to user
GMST0 = GMST_UNIFORM(MJD0)
aGEO  = (GE*(secsPerSidDay/twopi)**2)**(1._dk/3._dk)
period = twopi*sqrt(COE0(1)**3/GE)/secsPerSidDay

write(*,'(a,g15.8)') 'Tolerance: ',tol
write(*,'(a,i2)') 'Equations: ',eqs

call DPROP_REGULAR(R0,V0,tspan,tstep,cart,int_steps,tot_calls)

! ==============================================================================
! 04. PROCESSING AND OUTPUT
! ==============================================================================

! Initialize orbital elements array
npts = size(cart,1)
if (allocated(orb)) deallocate(orb)
allocate(orb(1:npts,1:7))
orb = 0._dk

! End timing BEFORE converting back to orbital elements
call SYSTEM_CLOCK(toc)
cputime = real((toc-tic),dk)/real(rate,dk)
write(*,'(a,g9.2,a)') 'CPU time: ',cputime,' s'

! Convert to orbital elements.
! orb(1): MJD,  orb(2): a,  orb(3): e, orb(4): i
! orb(5): Om,   orb(6): om, orb(7): M
do ipt=1,npts
    orb(ipt,1) = cart(ipt,1)    ! Copy MJD
    R = cart(ipt,2:4)
    V = cart(ipt,5:7)
    call CART2COE(R,V,orb(ipt,2:7),GE)
    orb(ipt,4:7) = orb(ipt,4:7)/d2r

end do

! Dump output and copy input files to the output directory
call CREATE_OUT(id_cart)
call CREATE_OUT(id_orb)
call CREATE_OUT(id_stat)
call DUMP_TRAJ(id_cart,npts,cart)
call DUMP_TRAJ(id_orb,npts,orb)

! Write statistics line: calls, steps, CPU time, final time and orbital elements
write(id_stat,100) tot_calls, int_steps, tol, cputime, orb(npts,:)

100 format((2(i10,1x),9(es22.15,1x)))
end program THALASSA
