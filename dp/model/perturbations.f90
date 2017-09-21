module PERTURBATIONS
! Description:
!    Contains wrapper subroutines to include perturbations in the equations of
!    motion.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    d.amato@upm.es
!
! ==============================================================================

use KINDS,       only: dk
use NSGRAV
use SUN_MOON
use DRAG_EXPONENTIAL
use Atmosphere1976, only: ATMOS76 => Atmosphere,dens_SL => RHOZERO
use SRP
use PHYS_CONST,  only: Clm,Slm
use AUXILIARIES, only: DU,TU
implicit none


contains




function PPOTENTIAL(insgrav,GM,RE,r,rm,t)
! Description:
!    Computes the perturbing potential due to a non-spherical gravity field.
!    Units are the same as those of the inputs.
!
! ==============================================================================

! VARIABLES
implicit none
! Arguments
integer,intent(in)    ::  insgrav   ! Non-sph. gravity flag
real(dk),intent(in)   ::  GM				! Gravitational parameter
real(dk),intent(in)   ::  RE				! Equatorial radius
real(dk),intent(in)   ::  r(1:3)    ! Position in inertial coordinates
real(dk),intent(in)   ::  rm        ! Magnitude of position vector
real(dk),intent(in)   ::  t         ! Time (dimensionless)
real(dk)              ::  PPOTENTIAL

! ==============================================================================

PPOTENTIAL = 0._dk
if (insgrav /= 0) then
   PPOTENTIAL = POT_NSG(GM,RE,Clm,Slm,r,rm,t)

end if

end function




function PACC_EJ2K(insgrav,isun,imoon,idrag,iSRP,r,v,rm,t,gradU_sph_out)
! Description:
!    Computes the perturbing acceleration in the EMEJ2000 frame due to a non-sph
!    gravity field, Sun and Moon, drag and solar radiation pressure.
!    Units are DIMENSIONLESS.
!
! ==============================================================================

! MODULES
use PHYS_CONST,  only: GE,GE_nd,GS,GM,RE_nd,ERR_constant,secsPerDay,twopi,RE,&
&CD,A2M_Drag,pSRP,CR,A2M_SRP
use AUXILIARIES, only: DU,TU

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)  ::  r(1:3),rm           ! Radius vector and its magnitude
real(dk),intent(in)  ::  v(1:3)              ! Velocity vector
real(dk),intent(in)  ::  t                   ! Time (dimensionless)
integer,intent(in)   ::  insgrav,isun        ! Perturbation flags
integer,intent(in)   ::  imoon,idrag,iSRP    ! More perturbation flags
! OPTIONAL OUTPUT: gradient of potential in sph. coordinates, used by
! EDromo right-hand-side.
real(dk),optional,intent(out)  ::  gradU_sph_out(1:3)
! Function definition
real(dk)            ::  PACC_EJ2K(1:3)
! LOCALS
! Non-spherical gravity
real(dk)  ::  gradU_sph(1:3)
real(dk)  ::  p_nsg(1:3)
! Third bodies
real(dk)  ::  r_sun(1:3),v_sun(1:3),p_sun(1:3)
real(dk)  ::  r_moon(1:3),v_moon(1:3),p_moon(1:3)
! Drag, density (US76)
real(dk)  ::  p_drag(1:3)
real(dk)  ::  drag_term
real(dk)  ::  density,temp1,temp2
! Velocity wrt atmosphere (km/s)
real(dk)  ::  v_rel(1:3),v_atm(1:3),v_relNorm
! Dimensional quantities for drag computation
real(dk)  ::  r_D(1:3),v_D(1:3),h_D
real(dk)  ::  wE_D
! Solar radiation pressure
real(dk)  ::  p_SRP(1:3)


! ==============================================================================

PACC_EJ2K = 0._dk

! ==============================================================================
! 01. NON-SPHERICAL GRAVITY
! ==============================================================================

gradU_sph = 0._dk; p_nsg = 0._dk
if (insgrav /= 0) then
  ! Compute the potential and its derivatives in non-dimensional units.
  ! EARTH
  gradU_sph = POTPAR(GE_nd,RE_nd,Clm,Slm,r,rm,t)
  p_nsg     = ACC_NSG(r,rm,gradU_sph)

end if

PACC_EJ2K = p_nsg + PACC_EJ2K
if (present(gradU_sph_out)) then
  ! Save gradient of U in spherical coordinates if present
  gradU_sph_out = gradU_sph

end if

! ==============================================================================
! 02. LUNISOLAR PERTURBATIONS
! ==============================================================================

p_sun = 0._dk; p_moon = 0._dk
if (isun /= 0) then
  ! SUN
  call EPHEM(1,DU,TU,t,r_sun,v_sun)
  p_sun = ACC_THBOD_EJ2K_ND(r,r_sun,GE,GS)

end if

if (imoon /= 0 ) then
  ! MOON
  call EPHEM(2,DU,TU,t,r_moon,v_moon)
  p_moon = ACC_THBOD_EJ2K_ND(r,r_moon,GE,GM)

end if

PACC_EJ2K = p_sun + p_moon + PACC_EJ2K

! ==============================================================================
! 03. ATMOSPHERIC DRAG
! ==============================================================================

p_drag = 0._dk
if (idrag == 1) then
  ! Exponential density model (Vallado)
  ! Make quantities dimensional
  r_D = r*DU; v_D = v*DU*TU
  wE_D = ERR_constant*twopi/secsPerDay

  ! Compute drag and non-dimensionalize
  p_drag = DRAG_ACC([r_D,v_D],wE_D,RE,CD,A2M_Drag)
  p_drag = p_drag/(DU*TU**2)

else if (idrag == 2) then
  ! US76 Atmosphere (code by R. Carmichael)
  ! Make quantities dimensional
  r_D = r*DU; v_D = v*DU*TU
  h_D = sqrt(dot_product(r_D,r_D)) - RE

  ! Get density
  call ATMOS76(h_D,density,temp1,temp2)
  density = density*dens_SL                ! Dimensionalize to kg/m^3
  
  ! The following is taken from DRAG_ACC. Relative velocity wrt atmosphere
  wE_D = ERR_constant*twopi/secsPerDay
  v_atm = wE_D*[-r_D(2),r_D(1),0._dk]
  v_rel = v_D - v_atm
  v_relNorm = sqrt(dot_product(v_rel,v_rel))

  ! Acceleration (in km/s^2)
  drag_term = -0.5_dk*1.E3_dk*CD*A2M_Drag*density
  p_drag    = drag_term*v_relNorm*v_rel

  ! Acceleration (non-dimensionalized)
  p_drag = p_drag/(DU*TU**2)

end if

PACC_EJ2K = p_drag + PACC_EJ2K

! ==============================================================================
! 04. SOLAR RADIATION PRESSURE
! ==============================================================================

p_SRP = 0._dk
if (iSRP == 1) then
  ! If the Sun gravitational perturbation is disabled, get its ephemerides
  if (isun == 0) then
    call EPHEM(1,DU,TU,t,r_sun,v_sun)

  end if
  ! Computation is in dimensional units.
  p_SRP = SRP_ACC(pSRP,CR,A2M_SRP,r*DU,r_sun*DU)
  p_SRP = p_SRP/(DU*TU**2)

end if
PACC_EJ2K = p_SRP + PACC_EJ2K

end function PACC_EJ2K




end module PERTURBATIONS
