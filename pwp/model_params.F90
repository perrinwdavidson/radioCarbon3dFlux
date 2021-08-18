! ----------------------------------------------------------------------
!                       Price Weller Pinkel Model
! ----------------------------------------------------------------------
!                            Parameter Module
! ----------------------------------------------------------------------
! About:
! Define global parameters.
! ----------------------------------------------------------------------
! Define module --------------------------------------------------------
MODULE model_params

    ! Ensure all variables defined ------------------------
    IMPLICIT NONE

    ! Important model parameters --------------------------
    DOUBLE PRECISION, PARAMETER :: days = 1.0  ! days to run integration, d
    DOUBLE PRECISION, PARAMETER :: dt = 1.0E3  ! time step, s
    DOUBLE PRECISION, PARAMETER :: dz = 1.0  ! depth bin, m

    ! Other model parameters ------------------------------
    DOUBLE PRECISION, PARAMETER :: dtd = dt / 8.64E4 ! time step, d
    INTEGER, PARAMETER :: metu = 1  ! indicates idealized diurnal cycle
    INTEGER, PARAMETER :: mzts = 1  ! indicates idealized initial profiles
    INTEGER, PARAMETER :: nts = 4  ! number of data points for intialize arrays
    INTEGER, PARAMETER :: iwrt = 2  ! writing frequency to terminal
    INTEGER, PARAMETER :: istor = 1  ! indicates writing to output to file
    INTEGER, PARAMETER :: itfreq = 2  ! time writing frequency
    INTEGER, PARAMETER :: izfreq = 5  ! depth writing frequency

    ! Max dimensions --------------------------------------
    INTEGER, PARAMETER :: nzmax = 500  ! depth bin max
    INTEGER, PARAMETER :: nmet = 2500  ! time length max
    INTEGER, PARAMETER :: nnn = 5000  ! time steps max

    ! Numerical constants ---------------------------------
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265359  ! pi

    ! Physical parameters ---------------------------------
    DOUBLE PRECISION, PARAMETER :: rlat = 31.0  ! latitude, degrees north
    DOUBLE PRECISION, PARAMETER :: beta1 = 0.60  ! shortwave extinction coefficient
    DOUBLE PRECISION, PARAMETER :: beta2 = 20.0  ! longwave extinction coefficient
    DOUBLE PRECISION, PARAMETER :: udrag = 9999.0  ! decay rate of velocity from unresolved processes
    DOUBLE PRECISION, PARAMETER :: omega = 7.292E-5  ! frequency of earth rotation, s-1
    DOUBLE PRECISION, PARAMETER :: f = 2.0 * omega &  ! coriolis parameter
                                     * SIN(rlat * (pi / 180.0))
    DOUBLE PRECISION, PARAMETER :: qimax = 978.0  ! amplitude of solar insolation
    DOUBLE PRECISION, PARAMETER :: ql = -126.0  ! heat loss
    DOUBLE PRECISION, PARAMETER :: tx = 0.07  ! eastward wind-stress
    DOUBLE PRECISION, PARAMETER :: emp = 0.0  ! precipitation minus evaporation
    DOUBLE PRECISION, PARAMETER :: pqfaci = 10.0  ! duration of daylight in hours
    DOUBLE PRECISION, PARAMETER :: rkz = 0.0  ! diffusion coefficient
    DOUBLE PRECISION, PARAMETER :: rkzmax = (dz ** 2) / (2.0 * dt)  ! max diffusion coefficient
    DOUBLE PRECISION, PARAMETER :: dstab = (dt * rkz) / (dz ** 2)  ! CFL condition
    DOUBLE PRECISION, PARAMETER :: rb = 0.65  ! initial critical bulk richardson number
    DOUBLE PRECISION, PARAMETER :: rg = 0.25  ! initial critical gradient richardson number
    DOUBLE PRECISION, PARAMETER :: em1 = 0.0  ! coefficient of buoyancy production, 0
    DOUBLE PRECISION, PARAMETER :: em2 = 0.0  ! coefficient of shear production, 0
    DOUBLE PRECISION, PARAMETER :: em3 = 0.0  ! coefficient of ustar**3, 0
    DOUBLE PRECISION, PARAMETER :: g = 9.807  ! gravity, m s-2
    DOUBLE PRECISION, PARAMETER :: ro = 1.024E3  ! reference density
    DOUBLE PRECISION, PARAMETER :: cpw = 4183.3  ! heat capacity of water
    DOUBLE PRECISION, PARAMETER :: hcon = dz * ro * cpw  ! heat inventory

    ! Save output -----------------------------------------
    SAVE

END MODULE model_params

! ----------------------------------------------------------------------
