! ----------------------------------------------------------------------
!                       Price Weller Pinkel Model
! ----------------------------------------------------------------------
!                          Absorption Routine
! ----------------------------------------------------------------------
! About:
! Compute solar radiation absorption profile. This subroutine assumes
! two wavelengths, and a double exponential depth dependence for
! absorption. Subscript 1 is for red, non-penetrating light, and 2 is
! for blue, penetrating light. rs1 is the fraction assumed to be red.
! ----------------------------------------------------------------------
! Define subroutine ----------------------------------------------------
SUBROUTINE ABSORB(absrb, nz)

    ! Set module use:
    USE model_params

    ! Define variables -------------------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Input dimensions:
    DOUBLE PRECISION, INTENT(IN) :: nz

    ! Input arrays:
    DOUBLE PRECISION, DIMENSION(nz), INTENT(INOUT) :: absrb

    ! Physical parameters:
    DOUBLE PRECISION, PARAMETER :: abss = 0.0
    DOUBLE PRECISION, PARAMETER :: rs1 = 0.6
    DOUBLE PRECISION, PARAMETER :: rs2 = 1.0 - rs1

    ! Indices:
    INTEGER :: j

    ! Initialize arrays ------------------------------------------------
    DO j = 1, nz, 1

        absrb(j) = 0.0

    ENDDO

    ! Calculate adsorbtion ---------------------------------------------
    DO j = 1, nz, 1

        ! Bin edges:
        z1 = FLOAT(j - 1) * dz
        z2 = z1 + dz

        ! Amount of light penetration at each depth:
        z1b1 = z1 / beta1
        z2b1 = z2 / beta1
        z1b2 = z1 / beta2
        z2b2 = z2 / beta2

        ! Check for math overflows:
        IF (z2b1 .LT. 70.0) THEN

            absrb(j) = absrb(j) + (rs1 * (exp(-z1b1) - exp(-z2b1)))

        ENDIF
        IF (z2b2 .LT. 70.0) THEN

            absrb(j) = absrb(j) + (rs2 * (exp(-z1b2) - exp(-z2b2)))

        ENDIF

        ! Make sure that absrb(z) integrates to 1 at large depth:
        abss = abss + absrb(j)

    ENDDO

    ! Return values ----------------------------------------------------
    RETURN

ENDSUBROUTINE ABSORB

! ----------------------------------------------------------------------
