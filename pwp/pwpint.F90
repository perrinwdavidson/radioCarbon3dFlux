! ----------------------------------------------------------------------
!                       Price Weller Pinkel Model
! ----------------------------------------------------------------------
!                          PWP Implementation
! ----------------------------------------------------------------------
! About:
! This model also implements an energy budget form of an entrainment
! parameterization that is very similar to that described in Price,
! Mooers, Van Leer, JPO 8, 4, 582-599 (and references therein). This part
! of the model should be treated as developmental only, as there are
! several features that are arbitrary to this model (i.e., the depth of
! the ml during times of heating can be the grid interval, dz). To use
! this parameterization set the bulk Richardson number to zero, and set
! em1, or em2, or em3 to non-zero. The gradient Richardson number can be
! zero or not.
!   qi and ql are the insolation and heat loss (qi > 0, ql < 0 almost
! always) in W/m**2. emp is evaporation minus precipitation (m/sec) and
! (emp > if evaporation exceeds precipitation). tx and ty are the east
! and north stress components (Pa).
! ----------------------------------------------------------------------
! Define subroutine ----------------------------------------------------
SUBROUTINE PWPINT(t, s, absrb, nz)

    ! Set module use:
    USE model_params

    ! Define variables -------------------------------------------------
    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Input dimensions:
    DOUBLE PRECISION, INTENT(IN) :: nz

    ! Input arrays:
    DOUBLE PRECISION, DIMENSION(nz), INTENT(IN) :: t
    DOUBLE PRECISION, DIMENSION(nz), INTENT(IN) :: s
    DOUBLE PRECISION, DIMENSION(nz), INTENT(INOUT) :: absrb

    ! Define physical variables:
    DOUBLE PRECISION :: ucon
    DOUBLE PRECISION :: tr  ! reference temp
    DOUBLE PRECISION :: sr  ! reference sal
    DOUBLE PRECISION :: dr  ! reference density
    DOUBLE PRECISION :: alpha  ! thermal expansion coeff
    DOUBLE PRECISION :: beta  ! haline expansion coeff
    DOUBLE PRECISION :: energy  ! energy

    ! Functions:
    DOUBLE PRECISION :: SGT
    EXTERNAL :: SGT

    ! Save some variables to subroutine:
    SAVE ucon, tr, sr, dr, alpha, beta, energy

    ! Initialize variables ---------------------------------------------
    ! Velocity processes:
    ucon = 0.0
    IF (udrag .LT. 100.0) THEN

        ucon = 1.0 / (8.64E4 * udrag)

    ENDIF

    ! References:
    tr = t(1)
    sr = s(1)
    dr = SGT(tr, sr)

    ! Expansion coefficients:
    alpha = SGT(tr + 0.5, sr) - SGT(tr - 0.5, sr)
    beta = SGT(tr, sr + 0.5) - SGT(tr, sr - 0.5)

    ! Compute the absorption profile:
    CALL ABSORB(absrb, nz)

    ! Finished with initialization:
    RETURN

    ! <- I am here -> !

    ! Step forward one time step ---------------------------------------
    ENTRY PWPGO(qi, ql, emp, tx, ty, dml, dtl)

    ! Apply heat and fresh water fluxes to the top most grid cell:
    q1 = qi*absrb(1) + ql
    t(1) = t(1) + dt*q1/(dz*ro*cpw)

    ! Absorb solar radiation at depth:
      do 70 j=2,nz
      t(j) = t(j) + qi*absrb(j)*dt/(dz*ro*cpw)
   70 continue
      s(1) = s(1) + s(1)*emp*dt/dz
!
!  Compute the density (sigma units), and relieve static
!  instability, if it occurs.
!
      do 71 j=1,nz
      d(j) = dr + (t(j) - tr)*alpha + (s(j) - sr)*beta
   71 continue
!
      call mldep(ds1)
      call si(ml,sipe)
      if(ml.gt.1) call mix5(ml)
      call mldep(ds2)
!
!  At this point the density proifile should be statically stable.
!
!  Time step the momentum equation.
!
!  Rotate the current throughout the water column through an
!  angle equal to inertial rotation for half a time step.
!
      ang = f*dt/2.
      sa = sin(ang)
      ca = cos(ang)
!
      do 32 j=1,nz
      call rot(u(j),v(j),sa,ca)
   32 continue
!
!  Apply the wind stress to the mixed layer as it now exists.
!
!  Find the ml depth.
!
      call mldep(dml)
      ml = int(dml/dz)
!
      du = (tx/(dml*ro))*dt
      dv = (ty/(dml*ro))*dt
      do 33 j=1,ml
      u(j) = u(j) + du
      v(j) = v(j) + dv
   33 continue
!
!  Apply drag to the current (this is a horrible parameterization of
!  inertial-internal wave dispersion).
!
      if(ucon.gt.1.e-10) then
      do 37 j = 1,nz
      u(j) = u(j)*(1. - dt*ucon)
      v(j) = v(j)*(1. - dt*ucon)
   37 continue
      end if
!
!  Rotate another half time step.
!
      do 34 j = 1,nz
      call rot(u(j),v(j),sa,ca)
   34 continue
!
!  Finished with the momentum equation for this time step.
!
!  Cause the mixed-layer to deepen.
!
!
      if(rb.le.0.00001) go to 61
!
!  Come here to do the bulk Richardson number instability
!  form of mixing (as in PWP).
!
      nzm = nz - 1
!
      rvc = rb
      do 80 j = ml,nzm
      delr = (d(j+1) - d(1))/ro
      ds2 = (u(j+1) - u(1))**2 + (v(j+1) - v(1))**2
      ds2 = ds2 + 1.e-8
      h1 = float(j)*dz
      rv = g*delr*h1/ds2
      if(rv.gt.rvc) go to 81
!
      jp = j + 1
      call mix5(jp)
!
   80 continue
   81 continue
      go to 65
!
!
   61 continue
!
!
      if(em1.le.0..and.em2.le.0..and.em3.le.0.) go to 65
!
!  Come here to mix according to an energy balance closure;
!  this implementation follows Price, Mooers, Van Leer, JPO,
!  8, 582-599, July 1978. This is not part of the PWP model
!  per se, and has some rather unusual properties, so watch out !
!
!
      call mldep(dml)
      ml = dml/dz
      taumag = sqrt(tx**2 + ty**2)
      ustar = sqrt(taumag/ro)
      ustarc = taumag*ustar*dt
!
!  Compute the energy available to do entrainment.
!
      energy = energy + em3*ustarc + em1*sipe
!
      delkes = 0.
      penes = 0.
!
      nzm = nz - 1
      do 60 l=ml,nzm
!
!  Determine how much energy is required to entrain the next
!  level below the base of the mixed-layer.
!
      call pe(l+1,pene)
      if(energy.lt.pene) go to 65
      if(energy.gt.pene) then
!
!  Come here if sufficient energy is avaliable to entrain at
!  least one level into the mixed layer.
!
!  Compute the change in kinetic energy due to entrainment.
!
      eke1 = eke(l+1)
      call mix5(l+1)
      eke2 = eke(l+1)
      delke = eke1 - eke2
!
      penes = penes + pene
      delkes = delkes + delke
!
!  Add and subtract energy losses and gains due to entrainment.
!
      energy = energy + em2*delke - pene
!
!
      end if
   60 continue
!
!  End of the energy budget method for estimation entrainment.
!
   65 continue
!
!  Relieve gradient Richardson number instability if it occurs.
!
      if(rg.gt.0.) call rimix(nmix,kcd)
!
      call mldep(dml)
      dtl = dz*float(kcd)
!
      return
      end
!
      subroutine pe(m,pene)
!
!  Compute the energy required to entrain density level
!  d(m) into a mixed-layer above.
!
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, rg, rb, em1, em2, em3
!
      dr = d(m) - d(m-1)
      ho2 = dz*float(m)/2.
      pene = g*dz*ho2*dr
      return
      end
