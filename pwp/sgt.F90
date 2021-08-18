! ----------------------------------------------------------------------
!                       Price Weller Pinkel Model
! ----------------------------------------------------------------------
!                           Density Routine
! ----------------------------------------------------------------------
! About:
! A sigma-t subroutine taken from seaprop; sigma-t knudsen to calculate
! density for the depth bin.
! ----------------------------------------------------------------------
! Define function ------------------------------------------------------
FUNCTION SGT(t, s)

    ! Define variables ------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Input scalars:
    DOUBLE PRECISION :: t
    DOUBLE PRECISION :: s

    ! Functions:
    DOUBLE PRECISION :: SG0
    EXTERNAL :: SG0

    ! Dummy scalars:
    DOUBLE PRECISION :: sg

    ! Output scalars:
    DOUBLE PRECISION :: SGT

    ! Calculate density -----------------------------------
    ! Calculate initial density:
    sg = SG0(s)

    ! Updated density:
    SGT = (((((-1.43803061E-7*t - 1.98248399E-3)*t - 0.545939111)*t &
          + 4.53168426)*t) / (t + 67.26)) &
          + ((((1.667E-8*t - 8.164E-7)*t + 1.803E-5)*t)*sg &
          + ((-1.0843E-6*t + 9.8185E-5)*t - 4.7867E-3)*t + 1.0)*sg

    ! Return values ---------------------------------------
    RETURN

ENDFUNCTION SGT

! ----------------------------------------------------------------------
