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
! Settings file id and path
integer,parameter      ::  id_set = 13
character(len=4096)    ::  input_path
! Physical model
integer  ::  insgrav            ! Non-spherical gravity field flag.
integer  ::  isun               ! 0 = no Sun perturbation, 1 = otherwise.
integer  ::  imoon              ! 0 = no Moon perturbation, 1 = otherwise.
integer  ::  idrag              ! 0 = no atmospheric drag, 1 = Vallado model, 2 = US76 (Carmichael)
integer  ::  iSRP               ! 0 = no SRP, 1 = otherwise.
integer  ::  iephem             ! Ephemerides source. 1 = DE431 ephemerides. 2 = Simpl. Meeus & Brown
integer  ::  gdeg,gord          ! Gravitational potential - maximum degree and order
! Integrator settings
integer  ::  mxstep             ! Max. number of integration/output steps.
real(dk) ::  tol                ! Integrator tolerance.
! Equations of motion settings
integer  ::  eqs                ! Equations of motion type. 1 = Cowell,
!                                 2 = EDromo(t), 3 = EDromo(c), 4 = EDromo(l)
! Output settings
character(len=512)  ::  outpath
integer  ::  verb
integer  ::  iEQUAT,iECLIP

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
open(unit=id_set,file=adjustl(trim(input_path)),status='old',action='read')
read(id_set,'(a)') (dummy, i = 1,hlines)

read(id_set,'(a11,i3)') dummy, insgrav
read(id_set,'(a11,i3)') dummy, isun
read(id_set,'(a11,i3)') dummy, imoon
read(id_set,'(a11,i3)') dummy, idrag
read(id_set,'(a11,i3)') dummy, iSRP
read(id_set,'(a11,i3)') dummy, iephem
read(id_set,'(a11,i3)') dummy, gdeg
read(id_set,'(a11,i3,7(/))') dummy, gord
read(id_set,'(a11,e22.15)') dummy, tol
read(id_set,'(a11,e22.15)') dummy, tspan
read(id_set,'(a11,e22.15)') dummy, tstep
read(id_set,'(a11,e10.1,5(/))') dummy, rmxstep
read(id_set,'(a11,i3,4(/))') dummy, eqs
read(id_set,'(a11,i3)') dummy, verb
read(id_set,'(a11,a)') dummy, outpath
read(id_set,'(a11,i3)') dummy, iEQUAT
read(id_set,'(a11,i3)') dummy, iECLIP

close(id_set)

mxstep = int(rmxstep)

end subroutine READ_SETTINGS

end module SETTINGS
