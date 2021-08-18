! ----------------------------------------------------------------------
!                       Price Weller Pinkel Model
! ----------------------------------------------------------------------
!                            Main Program
! ----------------------------------------------------------------------
! Start program --------------------------------------------------------
PROGRAM pwp

    ! Use module:
    USE model_params

    ! Define variables -------------------------------------------------
    ! Ensure all variables defined ------------------------
    IMPLICIT NONE

    ! Arrays ----------------------------------------------
    DOUBLE PRECISION, DIMENSION(nzmax) :: t  ! temperature, C
    DOUBLE PRECISION, DIMENSION(nzmax) :: s  ! salinity, PSU
    DOUBLE PRECISION, DIMENSION(nzmax) :: u  ! eastward current, m s-1
    DOUBLE PRECISION, DIMENSION(nzmax) :: v  ! northward current, m s-1
    DOUBLE PRECISION, DIMENSION(nzmax) :: d  ! density, sigma
    DOUBLE PRECISION, DIMENSION(nzmax) :: adsrb  ! fraction of the solar insolation that will be absorbed within the grid cell j.
    DOUBLE PRECISION, DIMENSION(nzmax) :: z  ! depth bin array
    DOUBLE PRECISION, DIMENSION(nzmax) :: zo  ! initial depth array
    DOUBLE PRECISION, DIMENSION(nzmax) :: to  ! initial temperature array
    DOUBLE PRECISION, DIMENSION(nzmax) :: so  ! initial salinity array
    DOUBLE PRECISION, DIMENSION(nzmax) :: do  ! initial density array
    DOUBLE PRECISION, DIMENSION(nmet) :: daya  ! sample time in days
    DOUBLE PRECISION, DIMENSION(nmet) :: qia  ! shortwave radiation
    DOUBLE PRECISION, DIMENSION(nmet) :: qla  ! outgoing radition
    DOUBLE PRECISION, DIMENSION(nmet) :: txa  ! eastward wind stress
    DOUBLE PRECISION, DIMENSION(nmet) :: tya  ! nortward wind stress
    DOUBLE PRECISION, DIMENSION(nmet) :: empa  ! precpiration minus evaporation

    ! Physical variables ----------------------------------
    DOUBLE PRECISION :: pqfac  ! non-dimensionalized duration of daylight per diurnal cucle.

    ! Other model variables -------------------------------
    CHARACTER(16) :: ifile  ! output file name
    INTEGER :: nzo ! number of initial depth bins
    INTEGER :: nz  ! number of depth bins

    ! Define indices --------------------------------------
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k

    ! Functions -------------------------------------------
    ! Kind:
    DOUBLE PRECISION :: SGT

    ! Specify dependence:
    EXTERNAL :: SGT

    ! Initialize arrays ------------------------------------------------
    ! Depth:
    zo(1) = 0.0
    zo(2) = 10.0
    zo(3) = 11.0
    zo(4) = 250.0

    ! Temperature:
    to(1) = 20.0
    to(2) = 20.0
    to(3) = 10.0
    to(4) = 10.0

    ! Salinity:
    so(1) = 36.0
    so(2) = 36.0
    so(3) = 36.0
    so(4) = 36.0

    ! Determine number of depth bins:
    nzo = INT(zo(nts) / dz)
    nz = nzo

    ! Interpolate arrays -----------------------------------------------
    ! Loop through all depth bins -------------------------
    DO j = 1, nz, 1

        ! Calculate depth:
        z(j) = FLOAT(j - 1) * dz

        ! Interpolate temperature:
        CALL INTX(zo, to, nts, z(j), t(j))

        ! Interpolate salinity:
        CALL INTX(zo, so, nts, z(j), s(j))

        ! Calculate density:
        d(j) = SGT(t(j), s(j))

    ENDDO

    ! Quality Control parameters ---------------------------------------
    ! Ensure daylight fraction a minimum value ------------
    ! Nondimensionalize:
    pqfac = pqfaci / 12.0

    ! Check to make sure acceptable minimum:
    IF (pqfac .LT. 0.00001) THEN

        pqfac = 0.00001

    ENDIF

    ! Ensure max depth bin is actually max ----------------
    IF (nz .GT. nzmax) THEN

        nz = nzmax

    ENDIF

    ! Alert to unstable diffusion -------------------------
    IF (dstab .GT. 0.5) THEN

        PRINT *, 'Warning, this value of kz will be unstable.'

    ENDIF

    ! Data storage -----------------------------------------------------
    ifile = 'outputs/pwp.out'
    OPEN(UNIT=7, FILE=ifile)

    ! Initialize model -------------------------------------------------
    CALL PWPINT(t, s, absrb, nz)

    ! Begin main loop --------------------------------------------------
    DO m = 1, nnn, 1

    ! <- I am here -> !

    ! Compute the solar insolation as a function of the time, where
    ! timed is time since start in days:
        timed = float(m-1)*dtd
        if (timed.gt.days) then
            go to 499
        endif
!
        ! initialize values:
      p2 = 3.141*2.
      tims = amod(timed,1.)
      tims = tims - 0.50
      qi = 0.
      if(cos(p2*tims).gt.0.) then
      qi = qimax*cos(tims*p2/pqfac)
      if(qi.lt.0.) qi = 0.
      end if
      emp = 0.
      end if
!
!  nn is the number of time steps actually taken.
!
      nn = nn + 1
!
!  Time step the model by an increment of time, dt.
!
      call pwpgo(qi,ql,emp,tx,ty,dml,dtl)
!
!
!  Apply a "background" diffusion if rkz is non-zero.
!
      if(rkz.gt.0) then
!
      call diffus (rkz,nz,dz,dt,d)
      call diffus (rkz,nz,dz,dt,s)
      call diffus (rkz,nz,dz,dt,t)
      call diffus (rkz,nz,dz,dt,u)
      call diffus (rkz,nz,dz,dt,v)
!
      end if
!
!  Store the density profile if this is the first step.
!
      if(m.eq.1) then
      do 2278 j=1,nz
      do(j) = d(j)
 2278 continue
      end if
!
!
!  This is the place to compute diagnostics.
!
!  Find max/min values on the first day.
!
      if(m.eq.1) dmlmin = 100.
      if(m.eq.1) tmin = t(1)
      if(timed.lt.1.) then
      if(t(1).lt.tmin) tmin = t(1)
      if(t(1).gt.tmax) tmax = t(1)
      spd = sqrt(u(1)**2 + v(1)**2)
      if(spd.gt.spdmax) spdmax = spd
      if(dml.lt.dmlmin) dmlmin = dml
      end if
!
!
!  Write out a few things on every iwrtth time step.
!
      if(iwrt.lt.1) go to 1339
!
      if(m.eq.1) write (6,444)
  444 format (1x,/,1x,'model solution variables are:',/,1x,
     x '   m timed    qi     ql      tx   t(1)   s(1)   u(1) ',
     x '  v(1)  dml   dtl',/,1x)
!
      if(mod(m,iwrt).ne.0) go to 1339
      write (6,133) m,timed,qi,ql,tx,t(1),s(1),
     ! u(1),v(1),dml,dtl
  133 format (1x,i4,f6.2,2f7.0,f7.2,4f7.2,2f6.0)
!
 1339 continue
!
!
!  Compute the net heat flux as a check on heat conservation.
!
      qnet = qnet + dt*(qi + ql)
!
!
!  Store data on disk, logical unit 7, if you told it to.
!
      if(istor.eq.1) then
      if(mod(m,itfreq).eq.0) then
      write (7,3188) timed,qi,ql,tx,t(1),s(1),u(1),v(1),
     x dml
 3188 format (1x,20f10.3)
      end if
      end if
!
    ENDDO

  499 continue
!
!
      write (6,870) tmin,tmax,dmlmin,spdmax
  870 format (1x,/,1x,'tmin, tmax, dmlmin, spdmax are',4f8.2)
!
!
!  Compute heat and potential energy diagnostics.
!  dpe will be the change in potential energy over the
!  integration period, and heat is the change in heat
!  content. heat should nearly equal the time-integrated
!  surface heat flux, qnet, if heat has been conserved.
!
      dpe = 0.
      heat = 0.
!
      do 8376 j=1,nz
      call intx(zo,to,nzs,nts,z(j),ti)
      tanom = t(j) - ti
      heat = heat + tanom*dz
      if(j.eq.1) tanoms = tanom
!
      danom = d(j) - do(j)
      dpe = dpe - g*danom*z(j)*dz

!     write (6,4465) j,do(j),d(j),danom,dpe
 4465 format (1x,'do,d,danom,dpe',i4,3f10.3,e12.2)
 8376 continue
!
      qnet = qnet/(ro*cpw)
      write (6,3877) heat, qnet, dpe
 3877 format (1x,/,1x,'heat, qnet, dpe are', 3f9.2)
!
!  End of the diagnostics.
!
      end
!

!
!
      subroutine stir (rc,r,a,j)
!
!
!  This subroutine mixes cells j and j+1 just enough so that
!  the Richardson number after the mixing is brought up to
!  the value rnew. In order to have this mixing process
!  converge, rnew must exceed the critical value of the
!  richardson number where mixing is presumed to start. if
!  r critical = rc = 0.25 (the nominal value), and r = 0.20, then
!  rnew = 0.3 would be reasonable. If r were smaller, then a
!  larger value of rnew - rc is used to hasten convergence.
!
!  This subroutine was modified by JFP in Sep 93 to allow for an
!  aribtrary rc and to achieve faster convergence.
!
!  rc is the value of the critical gradient Richardson number,
!  r is the Richardson number before this mixing,
!  a is the array to be mixed,
!  j and j+1 define the cell (grid level) where r was found to
!   be critical.
!
      dimension a(500)
!
!  Set the convergence parameter.
!
      rcon = 0.02 + (rc - r)/2.
      rnew = rc + rcon/5.               !  an arbitrary change
!
      f = 1. - r/rnew
      da = (a(j+1) - a(j))*f/2.
      a(j+1) = a(j+1) - da
      a(j) = a(j) + da
      return

END

! ----------------------------------------------------------------------
