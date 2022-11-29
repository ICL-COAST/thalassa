module CTHALASSA_TYPES
    use, intrinsic :: iso_c_binding
    implicit none

    type, bind(C) :: THALASSA_PHYSICALMODEL_STRUCT
        ! Physical model
        integer(c_int) :: insgrav
        integer(c_int) :: isun
        integer(c_int) :: imoon
        integer(c_int) :: idrag
        integer(c_int) :: iF107
        integer(c_int) :: iSRP
        integer(c_int) :: iephem
        integer(c_int) :: gdeg
        integer(c_int) :: gord        
    end type THALASSA_PHYSICALMODEL_STRUCT

    type, bind(C) :: THALASSA_PROPAGATOR_STRUCT
        ! Integrations
        real(c_double) :: tol
        real(c_double) :: tspan
        real(c_double) :: tstep
        real(c_double) :: mxstep
        integer(c_int) :: imcoll

        ! Equations of motion
        integer(c_int) :: eqs
    end type THALASSA_PROPAGATOR_STRUCT

    type, bind(C) :: THALASSA_OBJECT_STRUCT
        ! Physical characteristics
        real(c_double) :: mass
        real(c_double) :: area_drag
        real(c_double) :: area_srp
        real(c_double) :: cd
        real(c_double) :: cr
    end type THALASSA_OBJECT_STRUCT

    type, bind(C) :: THALASSA_STATE_STRUCT
        ! Epoch
        real(c_double) :: mjd

        ! State
        real(c_double) :: RV(6)
    end type THALASSA_STATE_STRUCT

    type, bind(C) :: THALASSA_PATHS_STRUCT
        character(kind=c_char) :: phys_path(512)
        integer(c_size_t)      :: phys_path_len
        character(kind=c_char) :: earth_path(512)
        integer(c_size_t)      :: earth_path_len
        character(kind=c_char) :: kernel_path(512)
        integer(c_size_t)      :: kernel_path_len
    end type THALASSA_PATHS_STRUCT

end module CTHALASSA_TYPES
