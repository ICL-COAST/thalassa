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
use US76_PATRIUS
use SRP
use PHYS_CONST,  only: Clm,Slm
use AUXILIARIES, only: DU,TU
implicit none

! Jacchia 77 dynamical atmospheric model
external ISDAMO
! NRLMSISE-00 atmospheric model
external GTD7,GTD7D,METERS




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




function PACC_ICRF(insgrav,isun,imoon,idrag,iSRP,r,v,rm,t,gradU_sph_out)
! Description:
!    Computes the perturbing acceleration in the ICRF frame due to a non-sph
!    gravity field, Sun and Moon, drag and solar radiation pressure.
!    Units are DIMENSIONLESS.
!
! ==============================================================================

! MODULES
use PHYS_CONST,  only: GE,GE_nd,GS,GM,RE_nd,ERR_constant,secsPerDay,twopi,RE,&
&CD,A2M_Drag,pSRP_1au,au,CR,A2M_SRP,MJD_J1950,GMST_UNIFORM,F107,Kp,Ap,JD2CAL,&
&r2d,cutoff_height
use AUXILIARIES, only: DU,TU,MJD0

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
real(dk)            ::  PACC_ICRF(1:3)
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
real(dk)  ::  density,pressure,temperature
! Drag and associated q.ties (J77, NRLMSISE-00)
real(dk)  ::  RA,DEC
real(dk)  ::  RA_sun,DEC_sun,r_sun_m
real(dk)  ::  JD_UT1,MJD_UT1,RJUD,DAFR,GMST,GMST_deg
real(dk)  ::  tempK(1:2),nDens(1:6),wMol
real(dk)  ::  SEC
real(dk)  ::  GLAT_deg,GLONG_deg,RA_deg,STL_hrs
real(dk)  ::  dens_MSIS00(1:6),temp_MSIS00(1:2)
integer   ::  IYD
! Time and date
integer   ::  year,month
real(dk)  ::  dayOfMonth,dayOfYear
! Velocity wrt atmosphere (km/s)
real(dk)  ::  v_rel(1:3),v_atm(1:3),v_relNorm
! Dimensional quantities for drag computation
real(dk)  ::  r_D(1:3),v_D(1:3),h_D
real(dk)  ::  wE_D
! Solar radiation pressure
real(dk)  ::  p_SRP(1:3)


! ==============================================================================

PACC_ICRF = 0._dk

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

PACC_ICRF = p_nsg + PACC_ICRF
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
  call EPHEM_ICRF(1,DU,TU,t,r_sun,v_sun)
  p_sun = ACC_THBOD_EJ2K_ND(r,r_sun,GE,GS)

end if

if (imoon /= 0 ) then
  ! MOON
  call EPHEM_ICRF(2,DU,TU,t,r_moon,v_moon)
  p_moon = ACC_THBOD_EJ2K_ND(r,r_moon,GE,GM)

end if

PACC_ICRF = p_sun + p_moon + PACC_ICRF

! ==============================================================================
! 03. ATMOSPHERIC DRAG
! ==============================================================================
! NOTE:
! Currently, the following approximations are made in the computation of the
! atmospheric drag:
! - Solar flux and geomagnetic activity are both constant (their values are
!   specified in ./in/physical_constants.txt)
! - Geodetic height is assumed = geometric height and geodetic
!   latitude/longitude = geometric, i.e. the ellipticity of the Earth is
!   neglected.

p_drag = 0._dk
! Make quantities dimensional
r_D = r*DU; v_D = v*DU*TU
h_D = sqrt(dot_product(r_D,r_D)) - RE
if (idrag /= 0 .and. h_D <= cutoff_height) then
  select case (idrag)
    case (1)
      ! Piecewise exponential density model (Vallado)
      density = ATMOS_VALLADO(h_D)
    
    case (2)
      ! US76 Atmosphere (Fortran porting of the PATRIUS 3.4.1 code)
      call US76_UPPER(h_D,temperature,pressure,density)
    
    case (3)
      ! Jacchia 77 atmosphere (code by V. Carrara - INPE)
      density = 0._dk
      if (h_D <= 2000._dk) then ! Check on max altitude for J77
        ! Right ascension and declination
        RA  = atan2(r(2),r(1))
        RA  = mod(RA + twopi,twopi)
        DEC = asin(r(3)/rm)
        
        ! Ephemerides of the Sun (if they haven't been computed earlier)
        if (isun == 0) then
          call EPHEM_ICRF(1,DU,TU,t,r_sun,v_sun)
        
        end if
        r_sun_m = sqrt(dot_product(r_sun,r_sun))
        RA_sun  = atan2(r_sun(2),r_sun(1))
        RA_sun  = mod(RA_sun +twopi,twopi)
        DEC_sun = asin(r_sun(3)/r_sun_m)
        MJD_UT1 = MJD0 + t/TU/secsPerDay
        RJUD    = MJD_UT1 - MJD_J1950
        DAFR    = RJUD - int(RJUD)
        GMST    = GMST_UNIFORM(MJD_UT1)
        call ISDAMO([RA,DEC,h_D*1.E3_dk],[RA_sun,DEC_sun],&
        & [F107,F107,Kp],RJUD,DAFR,GMST,tempK,nDens,wMol,density)
      
      end if
    
    case (4)
      ! NRLMSISE-00 atmospheric model
      density = 0._dk
      ! Get date and year
      MJD_UT1 = MJD0 + t/TU/secsPerDay
      JD_UT1  = MJD_UT1 + 2400000.5_dk
      call JD2CAL(JD_UT1,year,month,dayOfMonth,dayOfYear)
      
      ! Note: year number is ignored in NRLMSISE-00.
      IYD       = int(dayOfYear)
      SEC       = (dayOfYear - IYD)*secsPerDay
      GMST_deg  = GMST_UNIFORM(MJD_UT1)*r2d
      GLAT_deg  = asin(r(3)/rm)*r2d
      RA_deg    = mod(atan2(r(2),r(1)) + twopi, twopi)*r2d
      GLONG_deg = mod(RA_deg - GMST_deg + 360._dk,360._dk)
      STL_hrs  = SEC/3600._dk + GLONG_deg/15._dk
      ! Compute density
      call METERS(.true.)
      if (h_D <= 500._dk) then
        call GTD7(IYD,SEC,h_D,GLAT_deg,GLONG_deg,STL_hrs,F107,F107,Ap,48,&
        &dens_MSIS00,temp_MSIS00)
      
      else
        ! Do not neglect contribution from anomalous oxygen
        call GTD7D(IYD,SEC,h_D,GLAT_deg,GLONG_deg,STL_hrs,F107,F107,Ap,48,&
        &dens_MSIS00,temp_MSIS00)
      
      end if
      density = dens_MSIS00(6)

    end select

    ! Velocity wrt atmosphere
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

PACC_ICRF = p_drag + PACC_ICRF

! ==============================================================================
! 04. SOLAR RADIATION PRESSURE
! ==============================================================================

p_SRP = 0._dk
if (iSRP == 1) then
  ! If the Sun gravitational perturbation is disabled, get its ephemerides
  if (isun == 0) then
    call EPHEM_ICRF(1,DU,TU,t,r_sun,v_sun)

  end if
  ! Computation is in dimensional units (m/s^2).
  p_SRP = SRP_ACC(pSRP_1au,au,CR,A2M_SRP,r*DU,r_sun*DU)
  p_SRP = p_SRP/((DU*1.E3_dk)*TU**2)

end if
PACC_ICRF = p_SRP + PACC_ICRF

end function PACC_ICRF




function PACC_MMEIAUE(insgrav,isun,iearth,iSRP,r,v,rm,t,gradU_sph_out)
! Description:
!    Computes the perturbing acceleration in the MMEIAUE frame due to a non-sph
!    gravity field (TO DO), Sun (TO DO) and Earth, and solar radiation pressure
!    (TO DO).
!    Units are DIMENSIONLESS.
!
! ==============================================================================

! MODULES
use PHYS_CONST,  only: GM,GE,GS
use AUXILIARIES, only: DU,TU

! VARIABLES
implicit none
! Arguments IN
real(dk),intent(in)  ::  r(1:3),rm           ! Radius vector and its magnitude
real(dk),intent(in)  ::  v(1:3)              ! Velocity vector
real(dk),intent(in)  ::  t                   ! Time (dimensionless)
integer,intent(in)   ::  insgrav,isun        ! Perturbation flags
integer,intent(in)   ::  iearth,iSRP         ! More perturbation flags
! OPTIONAL OUTPUT: gradient of potential in sph. coordinates, used by
! EDromo right-hand-side.
real(dk),optional,intent(out)  ::  gradU_sph_out(1:3)
! Function definition
real(dk)            ::  PACC_MMEIAUE(1:3)

! LOCALS
! Non-spherical gravity
real(dk)  ::  p_nsg(1:3)
! Third-body perturbations
real(dk)  ::  p_sun(1:3),p_earth(1:3)
real(dk)  ::  r_sun(1:3),v_sun(1:3),r_earth(1:3),v_earth(1:3)
! Solar radiation pressure
real(dk)  ::  p_SRP(1:3)

! ==============================================================================

PACC_MMEIAUE = 0._dk

! ==============================================================================
! 01. NON-SPHERICAL GRAVITY (TO DO)
! ==============================================================================

p_nsg = 0._dk
! Insert computation of non-spherical gravity field of the Moon here
PACC_MMEIAUE = p_nsg + PACC_MMEIAUE

! ==============================================================================
! 02. THIRD-BODY PERTURBATIONS
! ==============================================================================

p_sun = 0._dk; p_earth = 0._dk
if (isun /= 0) then
  ! SUN
  call EPHEM_ICRF(1,DU,TU,t,r_sun,v_sun)

  ! Convert to MMEIAUE CS (TBD)

  p_sun = ACC_THBOD_EJ2K_ND(r,r_sun,GE,GS)

end if

if (iearth /= 0 ) then
  ! EARTH
  ! Position and velocity of the Earth in MMEIAUE are opposite of position and
  ! velocity of the Moon in ICRF.
  call EPHEM_ICRF(2,DU,TU,t,r_earth,v_earth)
  r_earth = -r_earth; v_earth = -v_earth
  p_earth = ACC_THBOD_EJ2K_ND(r,r_earth,GM,GE)

end if

PACC_MMEIAUE = p_sun + p_earth + PACC_MMEIAUE

! ==============================================================================
! 03. SOLAR RADIATION PRESSURE (TBD)
! ==============================================================================

p_SRP = 0._dk

PACC_MMEIAUE = p_SRP + PACC_MMEIAUE


end function PACC_MMEIAUE


end module PERTURBATIONS
