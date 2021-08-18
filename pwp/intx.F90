! ----------------------------------------------------------------------
!                       Price Weller Pinkel Model
! ----------------------------------------------------------------------
!                        Inteprolation Routine
! ----------------------------------------------------------------------
! About:
! This subroutine interploates the array s to find the value at the
! point ri, where r is the independent variable. The interpolation is
! linear, and assumes that r increases montonically, and that ri is
! within the range of r. If ri is not within the range of r, si is set
! equal to the appropriate endpoint of s, and a warning is printed.
! ----------------------------------------------------------------------
! Define routine -------------------------------------------------------
SUBROUTINE INTX(r, s, n, ri, si)

    ! Define variables -------------------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Dimensions:
    INTEGER, INTENT(IN) :: n

    ! Input arrays:
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: r
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: s

    ! Input scalars:
    DOUBLE PRECISION, INTENT(IN) :: ri

    ! Dummy indices:
    INTEGER :: k

    ! Output scalars:
    DOUBLE PRECISION, INTENT(OUT) :: si

    ! Quality Control --------------------------------------------------
    ! Check for out of range ri - less than ---------------
    IF (ri .LT. r(1)) THEN

        ! Set appropriate boundaries:
        si = s(1)

        ! Print warning:
        PRINT *, 'WARNING: ri was below r(1) in call to INTX.'

    ENDIF

    ! Check for out of range ri - greater than ------------
    IF (ri .GT. r(n)) THEN

        ! Set appropriate boundaries:
        si = s(n)

        ! Print warning:
        PRINT *, 'Warning: ri was above r(n) in call to INTX.'

    ENDIF

    ! Find approriate bounding sample values --------------
    ! Loop through all values:
    DO k = 1, (n - 1), 1

        ! Test if bracketed:
        IF (ri .GE. r(k) .AND. ri .LE. r(k + 1)) THEN

            GOTO 1

        ENDIF

    ENDDO

    ! Either exit, if didn't bracket:
    GOTO 9

    ! Or interpolate, if did bracket:
1   CONTINUE

    ! Interpolate ------------------------------------------------------
    si = ((s(k) * (r(k + 1) - ri)) + (s(k + 1) * (ri - r(k)))) &
        / (r(k + 1) - r(k))

    ! End --------------------------------------------------------------
    ! Move on to ending routine:
9   CONTINUE

    ! Return values:
    RETURN

ENDSUBROUTINE INTX

! ----------------------------------------------------------------------
