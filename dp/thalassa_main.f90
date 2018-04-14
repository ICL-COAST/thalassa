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
use IO,          only: id_cEQUAT,id_oEQUAT,id_cECLIP,id_oECLIP,id_stat,object_path
use SETTINGS,    only: gdeg,gord,outpath,tol,eqs,iEQUAT,iECLIP
use PHYS_CONST,  only: GE,d2r,r2d,secsPerDay,secsPerSidDay,twopi
use RFRAMES,     only: EQ2EC
implicit none

! VARIABLES
! Initial conditions (EMEJ2000)
real(dk)  ::  COE0(1:6),COE0_rad(1:6)
real(dk)  ::  R0_EQUAT(1:3),V0_EQUAT(1:3)
real(dk)  ::  GMST0,aGEO
real(dk)  ::  period
! Integration span and dt [solar days]
real(dk)  ::  tspan,tstep
! Trajectory
integer               ::  npts,ipt
real(dk),allocatable  ::  cart_EQUAT(:,:),orb_EQUAT(:,:)
real(dk),allocatable  ::  cart_ECLIP(:,:),orb_ECLIP(:,:)
real(dk)              ::  R_EQUAT(1:3),V_EQUAT(1:3)
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
call COE2CART(COE0_rad,R0_EQUAT,V0_EQUAT,GE)

! Output to user
GMST0 = GMST_UNIFORM(MJD0)
aGEO  = (GE*(secsPerSidDay/twopi)**2)**(1._dk/3._dk)
period = twopi*sqrt(COE0(1)**3/GE)/secsPerSidDay

write(*,'(a,g15.8)') 'Tolerance: ',tol
write(*,'(a,i2)') 'Equations: ',eqs

call DPROP_REGULAR(R0_EQUAT,V0_EQUAT,tspan,tstep,cart_EQUAT,int_steps,tot_calls)

! ==============================================================================
! 04. PROCESSING AND OUTPUT
! ==============================================================================

! Initialize orbital elements array
npts = size(cart_EQUAT,1)
if (allocated(orb_EQUAT)) deallocate(orb_EQUAT)
allocate(orb_EQUAT(1:npts,1:7))
orb_EQUAT = 0._dk

! End timing BEFORE converting back to orbital elements
call SYSTEM_CLOCK(toc)
cputime = real((toc-tic),dk)/real(rate,dk)
write(*,'(a,g9.2,a)') 'CPU time: ',cputime,' s'

! Convert to orbital elements.
! orb_EQUAT(1): MJD,  orb_EQUAT(2): a,  orb_EQUAT(3): e, orb_EQUAT(4): i
! orb_EQUAT(5): Om,   orb_EQUAT(6): om, orb_EQUAT(7): M
do ipt=1,npts
    orb_EQUAT(ipt,1) = cart_EQUAT(ipt,1)    ! Copy MJD
    R_EQUAT = cart_EQUAT(ipt,2:4)
    V_EQUAT = cart_EQUAT(ipt,5:7)
    call CART2COE(R_EQUAT,V_EQUAT,orb_EQUAT(ipt,2:7),GE)
    orb_EQUAT(ipt,4:7) = orb_EQUAT(ipt,4:7)/d2r

end do

! If additional output reference frames are selected, transform R,V
if (iEQUAT == 1) then
  call CREATE_OUT(id_cEQUAT)
  call CREATE_OUT(id_oEQUAT)
  call DUMP_TRAJ(id_cEQUAT,npts,cart_EQUAT)
  call DUMP_TRAJ(id_oEQUAT,npts,orb_EQUAT)

end if

if (iECLIP == 1) then
  if (allocated(cart_ECLIP)) deallocate(cart_ECLIP)
  if (allocated(orb_ECLIP)) deallocate(orb_ECLIP)
  allocate(cart_ECLIP,source=cart_EQUAT)
  allocate(orb_ECLIP,source=orb_EQUAT)
  
  ! Transform R_EQUAT, V_EQUAT into R_ECLIP, V_ECLIP
  call EQ2EC(npts,cart_EQUAT,cart_ECLIP)
  
  orb_ECLIP(:,1) = cart_EQUAT(:,1)
  do ipt=1,npts
    call CART2COE(cart_ECLIP(ipt,2:4),cart_ECLIP(ipt,5:7),orb_ECLIP(ipt,2:7),GE)

  end do
  orb_ECLIP(:,4:7) = orb_ECLIP(:,4:7)/d2r
  
  call CREATE_OUT(id_cECLIP)
  call CREATE_OUT(id_oECLIP)
  call DUMP_TRAJ(id_cECLIP,npts,cart_ECLIP)
  call DUMP_TRAJ(id_oECLIP,npts,orb_ECLIP)

end if

! Write statistics line: calls, steps, CPU time, final time and orbital elements
call CREATE_OUT(id_stat)
write(id_stat,100) tot_calls, int_steps, tol, cputime, orb_EQUAT(npts,:)

100 format((2(i10,1x),9(es22.15,1x)))
end program THALASSA
