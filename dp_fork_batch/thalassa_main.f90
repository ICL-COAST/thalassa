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
use AUXILIARIES, only: MJD0
use IO,          only: READ_IC,CREATE_OUT,DUMP_TRAJ
use CART_COE,    only: COE2CART,CART2COE
use PHYS_CONST,  only: READ_PHYS,GMST_UNIFORM
use PROPAGATE,   only: DPROP_REGULAR
use SETTINGS,    only: READ_SETTINGS
use IO,          only: id_cartF,id_orbF,id_stat
use SETTINGS,    only: model,gdeg,gord,outpath,statpath,tol,tol_lim,ntol,eqs
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
real(dk)              ::  MJD_fin,R_fin(1:3),V_fin(1:3)
real(dk)              ::  orb_fin(1:6)
! Measurement of CPU time
integer,parameter  ::  rep_time = 3
integer     ::  itime
integer(8)  ::  rate,tic,toc
real(dk) ::  cputime_sum,cputime_avg
! Function calls and integration steps
integer  ::  int_steps,tot_calls
! Batch tolerance loop
real(dk)  ::  delta_tol
real(dk),allocatable  ::  tolarr(:)
integer   ::  itol

! ==============================================================================


! ==============================================================================
! 01. INITIALIZATIONS
! ==============================================================================

! Read initial conditions, settings and physical model data.
call READ_IC(MJD0,COE0)
call READ_SETTINGS(tspan,tstep)
call READ_PHYS(model,gdeg,gord)

! Load SPICE kernels
call FURNSH('in/kernels_to_load.furnsh')

! ==============================================================================
! 02. TEST PROPAGATION
! ==============================================================================

! Convert to Cartesian coordinates
COE0_rad = [COE0(1:2),COE0(3:6)*real(d2r,dk)]
call COE2CART(COE0_rad,R0,V0,GE)

! Output to user
GMST0 = GMST_UNIFORM(MJD0)
aGEO  = (GE*(secsPerSidDay/twopi)**2)**(1._dk/3._dk)
period = twopi*sqrt(COE0(1)**3/GE)/secsPerSidDay

write(*,'(a)') 'INITIAL COORDINATES (aGEO, aGEO/sidDay):'
write(*,'(a,g22.15)') 'X = ',R0(1)/aGEO
write(*,'(a,g22.15)') 'Y = ',R0(2)/aGEO
write(*,'(a,g22.15)') 'Z = ',R0(3)/aGEO
write(*,'(a,g22.15)') 'VX = ',V0(1)/aGEO*(secsPerDay/1.0027379093508_dk)
write(*,'(a,g22.15)') 'VY = ',V0(2)/aGEO*(secsPerDay/1.0027379093508_dk)
write(*,'(a,g22.15)') 'VZ = ',V0(3)/aGEO*(secsPerDay/1.0027379093508_dk)
! write(*,'(a,g22.15)') 'Initial GMST (deg): ',GMST0*r2d
write(*,'(a,g22.15)') 'Initial orbital period (sid. days): ',period
write(*,*) 'Formulation: ',eqs
write(*,*) 'Tolmin, Tolmax: ',tol_lim
write(*,*) 'ntol: ',ntol

! Allocations
allocate(tolarr(1:ntol))

delta_tol = (log10(tol_lim(2)) - log10(tol_lim(1)))/(ntol - 1)
tolarr(1) = tol_lim(1)
do itol=2,ntol
    tolarr(itol) = 10._dk**(log10(tolarr(1)) + (itol - 1)*delta_tol)
end do

! Create output files and copy input files to the output directory
call SYSTEM('mkdir -p '//trim(outpath));    call SYSTEM('mkdir -p '//trim(statpath))
call SYSTEM('cp in/*.txt '//trim(outpath)); call SYSTEM('cp in/*.txt '//trim(statpath))
call CREATE_OUT(id_cartF)
call CREATE_OUT(id_orbF)
call CREATE_OUT(id_stat)

tol_loop: do itol=1,ntol
    tol = tolarr(itol)
    cputime_sum = 0._dk
    time_loop: do itime=1,rep_time
        ! Start clock
        call SYSTEM_CLOCK(tic,rate)
    
        call DPROP_REGULAR(R0,V0,tspan,tstep,cart,int_steps,tot_calls)

        ! End timing BEFORE converting back to orbital elements
        call SYSTEM_CLOCK(toc)
        cputime_sum = cputime_sum + real((toc-tic),dk)/real(rate,dk)
    
    end do time_loop
    cputime_avg = cputime_sum/real(rep_time,dk)
    write(*,*) 'Run ',itol,', CPU time (avg.): ',cputime_avg,'s'

    ! Save results and convert to orbital elements
    npts = size(cart,1)
    MJD_fin = cart(npts,1)
    R_fin  = cart(npts,2:4)
    V_fin  = cart(npts,5:7)
    call CART2COE(R_fin,V_fin,orb_fin,GE)

    ! Dump output
    write(id_cartF,'(2(g15.8,1x),7(es23.15,1x))') tol, cputime_avg, MJD_fin, R_fin, V_fin
    write(id_orbF,'(2(g15.8,1x),7(es23.15,1x))') tol, cputime_avg, MJD_fin,&
    & orb_fin(1:2), orb_fin(3:6)*r2d
    write(id_stat,'(2(g15.8,1x),es23.15,1x,2(i11,1x))')&
    & tol, cputime_avg, MJD_fin, int_steps, tot_calls

    ! Close and reopen files to write buffer to disk
    close(id_cartF); close(id_orbF); close(id_stat)
    open(id_cartF,file=trim(adjustl(trim(statpath)//'cart_fin.dat')),&
    &action='write',status='old',position='append')
    open(id_orbF,file=trim(adjustl(trim(statpath)//'orb_fin.dat')),&
    &action='write',status='old',position='append')
    open(id_stat,file=trim(adjustl(trim(statpath)//'stats.dat')),&
    &action='write',status='old',position='append')

end do tol_loop

! ==============================================================================
! 03. PROCESSING AND OUTPUT
! ==============================================================================    

! call CREATE_OUT(id_cart)
! call CREATE_OUT(id_orb)
! call CREATE_OUT(id_stat)
! call DUMP_TRAJ(id_cart,npts,cart)
! call DUMP_TRAJ(id_orb,npts,orb)

! ! Write statistics line: calls, steps, CPU time, final time and orbital elements
! write(id_stat,100) tot_calls, int_steps, cputime_avg, orb(npts,:)

! 100 format((2(i10,'',''),8(es22.15,'','')))

! Cleanup
deallocate(tolarr)
close(id_cartF)
close(id_orbF)
close(id_stat)

end program THALASSA
