module SETTINGS
! Description:
!    Settings for Thalassa.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================

! MODULES
use KINDS, only: dk

implicit none
! Settings file id
integer,parameter   ::  id_set = 13
! Physical model
integer  ::  model              ! 1 = STELA model, 2 = SWIFT model, 3 = custom.
integer  ::  insgrav            ! Non-spherical gravity field flag.
integer  ::  isun               ! 0 = no Sun perturbation, 1 = otherwise.
integer  ::  imoon              ! 0 = no Moon perturbation, 1 = otherwise.
integer  ::  idrag              ! 0 = no drag, 1 = exponential.
integer  ::  iephem             ! Ephemerides source. 1 = DE431 ephemerides.
integer  ::  gdeg,gord          ! Gravitational potential - maximum degree and order
! Integrator settings
integer  ::  integ              ! Integrator type. 1 = LSODAR, 2 = CVODE.
integer  ::  mxstep             ! Max. number of integration/output steps.
integer  ::  ntol               ! Number of runs in tolerance
real(dk) ::  tol                ! Integrator tolerance. (TODO)
real(dk) ::  tolarr(1:2)        ! Max and min tolerances for batch propagation.
! Equations of motion settings
integer  ::  eqs                ! Equations of motion type. 1 = Cowell,
!                                 2 = EDromo(t), 3 = EDromo(c), 4 = EDromo(l)
! Output directories
character(len=512)  ::  outpath
character(len=512)  ::  statpath

contains


subroutine READ_SETTINGS(tspan,tstep)
! Description:
!    Reads the settings file.
! ==============================================================================

implicit none
real(dk),intent(out) ::  tspan,tstep
integer,parameter    ::  hlines = 12  ! <- Check this when modifying input.txt
integer  :: i
character(len=4096)  ::  dummy
real(dk)  ::  rmxstep

! Open and skip header lines
open(unit=id_set,file='./in/input.txt',status='old',action='read')
read(id_set,'(a)') (dummy, i = 1,hlines)

read(id_set,'(a11,i3)') dummy, model
read(id_set,'(a11,i3)') dummy, insgrav
read(id_set,'(a11,i3)') dummy, isun
read(id_set,'(a11,i3)') dummy, imoon
read(id_set,'(a11,i3)') dummy, idrag
read(id_set,'(a11,i3)') dummy, iephem
read(id_set,'(a11,i3)') dummy, gdeg
read(id_set,'(a11,i3,8(/))') dummy, gord
read(id_set,'(a11,i3)') dummy, integ
read(id_set,'(a11,i3)') dummy, ntol
read(id_set,'(a11,2(e22.15,2x))') dummy, tolarr
read(id_set,'(a11,e22.15)') dummy, tspan
read(id_set,'(a11,e22.15)') dummy, tstep
read(id_set,'(a11,e10.1,4(/))') dummy, rmxstep
read(id_set,'(a11,i3,2(/))') dummy, eqs
read(id_set,'(a4,a)') dummy,outpath
read(id_set,'(a5,a)') dummy,statpath

close(id_set)

mxstep = int(rmxstep)

end subroutine READ_SETTINGS

end module SETTINGS
