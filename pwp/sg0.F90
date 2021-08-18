! ----------------------------------------------------------------------
!                       Price Weller Pinkel Model
! ----------------------------------------------------------------------
!                        Initial Density Routine
! ----------------------------------------------------------------------
! About:
! A sigma-0 subroutine neede by the sigma-t subroutine; sigma-0 knudsen
! to calculate density at depth.
! ----------------------------------------------------------------------
! Define function ------------------------------------------------------
FUNCTION SG0(s)

    ! Define Variables ------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Input value:
    DOUBLE PRECISION :: s

    ! Output value:
    DOUBLE PRECISION :: SG0

    ! Calculate density -----------------------------------
    SG0 = (((6.76786136E-6*s - 4.8249614E-4)*s + 0.814876577)*s) &
          - 0.0934458632

    ! Return ----------------------------------------------
    RETURN

ENDFUNCTION SG0

! ----------------------------------------------------------------------
