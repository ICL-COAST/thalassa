program THALASSA
! Description:
! Thalassa propagates Earth-bound orbits from an initial to a final epoch. It
! uses either Newtonian or regularized equations of motion written in the
! EMEJ2000 frame. It takes into account the following perturbations:
! - Non-spherical Earth gravitational potential
! - Third-body perturbations from Sun and Moon
! - Drag (TBD)
! - SRP  (TBD)
!
! The choice of the mathematical model and constants reflect those used in
! SWIFT and STELA, two widely-known propagators currently in usage for studying
! the MEO and GEO regions.
! Also, the user can choose between the following numerical integrators:
! - LSODAR
! - CVODE     (TBD)
! - DOPRI853  (TBD)
! - RKG       (TBD)
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================


! MODULES
use KINDS,       only: dk
use IO,          only: READ_IC,CREATE_OUT,DUMP_TRAJ
use CART_COE,    only: COE2CART,CART2COE
use COORD_SYST,  only: CHOOSE_CS
use PHYS_CONST,  only: READ_PHYS,GMST_UNIFORM,CURRENT_MU
use PROPAGATE,   only: DPROP_REGULAR
use SETTINGS,    only: READ_SETTINGS
use AUXILIARIES, only: SET_UNITS
use AUXILIARIES, only: MJD0,coordSyst
use IO,          only: id_cart,id_orb,id_stat
use SETTINGS,    only: model,gdeg,gord,outpath,mxstep
use PHYS_CONST,  only: GE,d2r,r2d,secsPerDay,secsPerSidDay,twopi
use PROPAGATE,   only: trajectory
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
! Current gravitational parameter
real(dk)              ::  mu
! Measurement of CPU time
integer  ::  rate,tic,toc
real(dk) ::  cputime
! Function calls and integration steps
integer  ::  int_steps,tot_calls

! Results
type(trajectory),allocatable  ::  traj_ICRF(:)


! ==============================================================================


! ==============================================================================
! 01. INITIALIZATIONS
! ==============================================================================

! Start clock
call SYSTEM_CLOCK(tic,rate)

! Read initial conditions, settings and physical model data.
call READ_IC(MJD0,R0,V0,coordSyst)
call READ_SETTINGS(tspan,tstep)
call READ_PHYS(model,gdeg,gord)

! Load SPICE kernels
call FURNSH('in/kernels_to_load.furnsh')

! ==============================================================================
! 02. TEST PROPAGATION
! ==============================================================================

! Select gravitational parameter
mu = CURRENT_MU(coordSyst)

! Check on initial coordinate system
! NOTE: CHOOSE_CS needs DU, TU to be set. SET_UNITS will be called again later.
call SET_UNITS(R0)
call CHOOSE_CS(0._dk,coordSyst,R0,V0,coordsyst,R0,V0)

! ! Output to user (TBD)
! GMST0 = GMST_UNIFORM(MJD0)
! aGEO  = (GE*(secsPerSidDay/twopi)**2)**(1._dk/3._dk)
! period = twopi*sqrt(COE0(1)**3/GE)/secsPerSidDay

! write(*,'(a)') 'INITIAL COORDINATES (aGEO, aGEO/sidDay):'
! write(*,'(a,g22.15)') 'X = ',R0(1)/aGEO
! write(*,'(a,g22.15)') 'Y = ',R0(2)/aGEO
! write(*,'(a,g22.15)') 'Z = ',R0(3)/aGEO
! write(*,'(a,g22.15)') 'VX = ',V0(1)/aGEO*(secsPerDay/1.0027379093508_dk)
! write(*,'(a,g22.15)') 'VY = ',V0(2)/aGEO*(secsPerDay/1.0027379093508_dk)
! write(*,'(a,g22.15)') 'VZ = ',V0(3)/aGEO*(secsPerDay/1.0027379093508_dk)
! write(*,'(a,g22.15)') 'Initial GMST (deg): ',GMST0*r2d
! write(*,'(a,g22.15)') 'Initial orbital period (sid. days): ',period

call DPROP_REGULAR(coordSyst,R0,V0,tspan,tstep,cart,int_steps,tot_calls,traj_ICRF)

! ==============================================================================
! 03. PROCESSING AND OUTPUT
! ==============================================================================

! Initialize orbital elements array
npts = size(traj_ICRF,1)
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
    orb(ipt,1) = traj_ICRF(ipt)%MJD    ! Copy MJD
    R = traj_ICRF(ipt)%RV(1:3)
    V = traj_ICRF(ipt)%RV(4:6)
!    mu = CURRENT_MU(traj_ICRF(ipt)%CS)
    call CART2COE(R,V,orb(ipt,2:7),CURRENT_MU('MMEIAUE'))
    orb(ipt,4:7) = orb(ipt,4:7)/d2r

end do

! Dump output and copy input files to the output directory
call SYSTEM('mkdir -p '//trim(outpath))
call SYSTEM('cp in/*.txt '//trim(outpath))
call CREATE_OUT(id_cart)
call CREATE_OUT(id_orb)
call CREATE_OUT(id_stat)
!call DUMP_TRAJ(id_cart,npts,cart)
call DUMP_TRAJ(id_orb,npts,orb)

! Write statistics line: calls, steps, CPU time, final time and orbital elements
write(id_stat,100) tot_calls, int_steps, cputime, orb(npts,:)

100 format((2(i10,'',''),8(es22.15,'','')))
end program THALASSA
