c
       program pwpmain
c
c
c  This main program is used to drive the pwp subroutine with
c  a daily cycle of heating, or, an arbitrary (observed) time
c  series of surface flux data.
c
c  The program assumes the following logical units:
c   5 is the terminal for input,
c   6 is the terminal for output,
c   7 may be opened for data storage (optional),
c   8 may be opened to read air/sea flux data (optional),
c   9 may be opened to read the initial t, s profile (optional).
c
c
c  Written by Jim Price, april 27, 1989. Please direct any
c  comments or questions to Jim Price, WHOI, Woods Hole, MA
c  02543, USA, tel. 508-289-2526.
c
c  Version 2:  changed the gradient Ri relaxation method.
c  Version 3:  clarified the bulk Ri relaxation.
c  Version 4:  fixed an error that disabled background diffusion.
c
c  Last editted on 20 Jan, 1999 by JFP.
c
c
      character*16 ifile
c
c  The arrays below store the values of:
c  east and north current (u, v),
c  temperature and salinity (t, s),
c  density (d, sigma units), and
c  depth of the grid points, z.
c
c  All units are mks, except density, d, which is in sigma
c  units.
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
      dimension z(500)
c
c  The arrays below are used to store the initial z,t,s.
c
      dimension zo(500), to(500), so(500), do(500)
c
c  The arrays below store times series of surface flux data.
c  Depending upon the time step and length of the record you
c  could easily run out of these.
c
      dimension daya(2500), qia(2500), qla(2500), txa(2500),
     x tya(2500), empa(2500)
c
c  nzmax should match the size of the initial profile arrays
c  dimensioned above, and nmet should match the flux arrays.
c
      nzmax = 500
      nmet = 2500
c
      write (6,800)
  800 format (
     x 1x, 'This program provides a means to run the pwp',/,
     x 1x, 'upper ocean model. To run this model you must',/,
     x 1x, 'prescribe the air/sea fluxes, the initial',/,
     x 1x, 'temperature and salinity profile, and several',/,
     x 1x, 'parameters defining the site and the model.',/,
     x 1x, 'The site and model parameters are listed below,',/,
     x 1x, 'where default values and units given in parentheses.',/,
     x 1x, 'To choose the default enter a comma.',/)
      write (6,802)
  802 format (
     x 1x, 'beta1, longwave extinction coefficient (0.6 m)',/,
     x 1x, 'beta2, shortwave extinction coefficient (20.0 m)',/,
     x 1x, 'rlat,  latitude (31.0)',/,
     x 1x, 'udrag, decay time scale of current (9999. days)',/,
     x 1x,/,
     x 1x, 'dt,    time step (1000. sec)',/,
     x 1x, 'days,  the number of days to run (1.0 day)',/,
     x 1x, 'dz,    grid interval (1.0 m)',/,
     x 1x,/,
     x 1x, 'rb,    critical bulk richardson number (0.65)',/,
     x 1x, 'rg,    critical gradient richardson number (0.25)',/,
     x 1x, 'em1,   coefficient of buoyancy production (0.)',/,
     x 1x, 'em2,   coefficient of shear production (0.)',/,
     x 1x, 'em3,   coefficient of ustar**3 (0.)',/)
c
c
c  Initialize rlat, the latitude.
c
      rlat = 31.
c
c  Initialize the solar radiation absorption;
c  beta1 is the  non-penetrating part of the insolation,
c  beta2 is the penetrating (shortwave) part.
c  values below are appropriate for fairly clear waters.
c
      beta1 = 0.60
      beta2 = 20.
c
c  Initialize udrag, the e-folding time (days) of current
c  due to unresolved processes.
c
      udrag = 9999.
c
c
      write (6,888)
  888 format (1x,'Enter beta1, beta2, rlat, udrag.')
      read (5,*)  beta1, beta2, rlat, udrag
c
c
c
c  dt and dz are the time step (sec) and the vertical grid
c  interval (m), and  days is the number of days to run
c  the integration.
c
      dt = 1.0e3
      dz = 1.0
      days = 1.
c
      write (6,706)
  706 format (1x,/,1x,'Enter dt, days, dz.')
      read (5,*) dt, days, dz
c
      dtd = dt/8.64e4
      f = 2.*7.292e-5*sin(rlat*3.141/180.)
c
c
c  Model constants are rb, the critical bulk richardson number,
c  rg, the critical gradient richardson number, and emn, the
c  coefficients that multiply the energy source terms in
c  the turbulent energy budget. Values set below are the PWP
c  model values.
c
      rb = 0.65
      rg = 0.25
      em1 = 0.
      em2 = 0.
      em3 = 0.
c
      write (6,772)
  772 format (1x,/,1x,'Enter rb, rg, em1, em2, em3.')
      read (5,*) rb,rg,em1,em2,em3
c
c
c
c  Set up to input air/sea flux data.
c
      write (6,457)
  457 format (1x,//,
     x 1x,'You can prescribe the air/sea flux data by:',/,
     x 1x,'specifying the parameters of an idealized diurnal',/,
     x 1x,'cycle from terminal input (enter 1 (default)), or, ',/,
     x 1x,'read the air/sea flux data from a disk file (enter 2).')
c
      metu = 1
      read (5,*) metu
c
      if(metu.eq.1) then
c
      write (6,2790)
 2790 format (1x,/,1x,'In order to define the air-sea fluxes ',/,
     x 1x, 'for a diurnal cycle, the program will now ',/,
     x 1x, 'query for the values of the following variables, all ',/,
     x 1x, 'of which may be defaulted to the values in parens.',/,
     x 1x, '(which are from day 131 of the PWP paper):',/)
      write (6,801)
  801 format (
     x 1x, 'qimax,  noon amplitude of insolation (978. w/m**2)',/,
     x 1x, 'ql,     steady heat loss (-126. w/m**2)',/,
     x 1x, 'tx,     steady east wind stress (0.07 pa)',/,
     x 1x, 'emp,    evaporation minus precipitation (0.0 m/sec)',/,
     x 1x, 'pqfac,  duration of insolation (10.0 hours)',/)
c
c  Initialize the amplitude of solar insolation (qi), the
c  heat loss (ql), and the wind stress (tx, uses east only).
c
      qimax = 978.
      ql = -126.
      tx = 0.07
      emp = 0.
c
c  Initialize pqfac, the duration of daylight in hours.
c
      pqfac = 10.
c
c
      write (6,777)
  777 format (1x,'Enter qimax, ql, tx, emp, pqfac.')
      read (5,*) qimax,ql,tx,emp,pqfac
c
c  Convert pqfac to a non-d form.
c
      pqfac = pqfac/12.
      if(pqfac.lt.0.00001) then
      pqfac = 0.00001
      end if

      end if
c
c  Come here if you intend to read a time series of flux data from
c  a disk file.
c
      if(metu.eq.2) then
c
      write (6,117)
  117 format (1x,/,1x,'Enter the name of the the air/sea flux file.')
      read (5,116) ifile
  116 format (a16)
      open (unit=8,file=ifile)
c
c
      do 233 j=1,nmet
      read (8,*,end=459) daya(j), qia(j), qla(j), txa(j), tya(j),
     x empa(j)
c
c
      if(daya(j).lt.-9.) go to 459
c
      nmeto = nmeto + 1
c
      if(nmeto.lt.25) then
      write (6,8838) nmeto, daya(j), qia(j), qla(j), txa(j), tya(j),
     x empa(j)
 8838 format (1x,'First 25 flux data:', i4, f7.2,2f7.0,2f7.2,e10.3)
      end if
c
c
  233 continue
  459 continue
c
c  Set nmet equal to actual number of observations
c
      nmet = nmeto
c
      write (6,8833)
 8833 format (1x,/,1x,'Enter the starting time (day).')
      read (5,*) day0
c
      end if
c
c  Finished with specifying air/sea flux data.
c
c
c
c  Initialize the z, t, s, d profiles.
c
c  Read in the data used to define the initial t,s profiles
c  from logical unit nunit (read in below). The data
c  is presumed to be free format, with each record being a
c  triple (z, t, s).  The first record should be the surface
c  value, i.e., 0., 20., 36. is a reasonable surface value,
c  and the last true value should be deep enough to avoid
c  entering into the caluclation (unless intentional).
c  The end of file is defined either by the end
c  of the data (if read from a disk file), or by a z value less
c  than -10.  These data are then linearly interpolated to give
c  profiles at dz resolution needed by the model.
c
      write (6,458)
  458 format (1x,//,
     x 1x,'The initial temperature and salinity profile can',/,
     x 1x,'be defined by:',/,
     x 1x,'data hardwired into the program (enter 1, default), or',/,
     x 1x,'enter z, t, s data from the terminal (enter 2), or,'/,
     x 1x,'read the data from a disk file (enter 3).')
c
      mzts = 1
      read (5,*) mzts
c
      if(mzts.eq.1) then
      nts = 4
      zo(1) = 0.
      zo(2) = 10.
      zo(3) = 11.
      zo(4) = 250.
      to(1) = 20.
      to(2) = 20.
      to(3) = 10.
      to(4) = 10.
      so(1) = 36.
      so(2) = 36.
      so(3) = 36.
      so(4) = 36.
      end if
      if(mzts.eq.1) go to 53
c
      if(mzts.eq.3) then
      write (6,3290)
 3290 format (1x,/,1x,'Enter the name of the z,t,s data file.')
      read (5,116) ifile
      open (unit=9,file=ifile)
      end if
c
c
      do 50 j=1,nzmax
      if(mzts.eq.2) then
      write (6,876)
  876 format (1x,'Enter a z,t,s (end = z < -10).')
      read (5,*) zo(j),to(j),so(j)
      end if
c
      if(mzts.eq.3) then
      read(9,*,end=53) zo(j),to(j),so(j)
      end if
c
      if(zo(j).lt.-9.) go to 53
      nts = nts + 1
c     if(mzts.eq.3) write (6,503) zo(j),to(j),so(j)
c  503 format (1x,'Initial z, t, s are:', 3f10.2)
   50 continue
   53 continue
c
c  Compute the size of the arrays actually used.
c
      nzo = int(zo(nts)/dz)
      nz = nzo
      if(nz.gt.nzmax) nz = nzmax
c
c  Interpolate the initial data to get a profile at dz resolution.
c
      do 52 j=1,nz
      z(j) = float(j-1)*dz
      call intx(zo,to,nzs,nts,z(j),t(j))
      call intx(zo,so,nzs,nts,z(j),s(j))
      d(j) = sgt(t(j),s(j),gg)
c      write (6,54) j,z(j),t(j),s(j),d(j)
   54 format (1x,'Density prof.:',i4,2x,6f10.3)
   52 continue
c
c  Initialization of z, t, s, d is finished.
c
c
c
c  Section below allows for application of a simple "background"
c  diffusion to be applied to the profiles. this is not part of
c  the pwp model per se, but is included here for the hell of it.
c
      rkz = 0.
c
      rkzmax = (dz**2)/(2.*dt)
      write (6,3355) rkzmax
 3355 format (1x,//,1x
     x 'If you want to have a background diffusion, then',
     x /,1x,'enter Kz (m**2/sec) (default = 0; max =',e10.3,')')
      read (5,*) rkz
      dstab = dt*rkz/dz**2
      if(dstab.gt.0.5) then
      write (6,8433) dstab
 8433 format (1x,/,1x,'Warning, this value of kz will be unstable:',
     x /,1x,'dstab = ', e10.3)
      end if
c
c
      write (6,8899)
 8899 format (1x,/,1x,'Do you want a screen listing of data ?',/,
     x 1x,'(0 = no, n = every nth step (default = 2)')
      iwrt = 2
      read (5,*) iwrt
c
c  Section below stores some data onto a disk file for later
c  analysis or plotting.
c
      istor = 0
      write (6,3186)
 3186 format (1x,/,1x,'Do you want to store the data on disk ?',/,1x,
     x '(0 = no (default), 1 = yes)')
      read (5,*) istor
c
      if(istor.eq.1) then
      write (6,1188)
 1188 format (1x,/,1x,'Enter the name of the file to store data.')
      read (5,116) ifile
      open (unit=7,file=ifile)
c
      itfreq = 2
      izfreq = 5
c
      write (6,3187)
 3187 format (1x,/,1x,'Enter the decimation rate for time and depth.',
     x /,1x,'(default is 2 and 5)')
      read (5,*) itfreq,izfreq
      end if
c
c  At this point, we are ready to call the initialization
c  entry point in the pwp subroutine.
c
      call pwpint(beta1,beta2,rlat,udrag)
c
c
c  Time step the model for up to nnn steps.
c
      nnn = 5000
c
c  nnn is set to 5000 to avoid running forever
c  in case the limit timed > days is missed.
c
      do 40 m=1,nnn
c
c  Compute the solar insolation as a function of the time,
c  where timed is time since start in days.
c
c
      timed = float(m-1)*dtd
      if(timed.gt.days) go to 499
c
      if(metu.eq.1) then
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
c
      if(metu.eq.2) then
c
      day = timed + day0
      call intx (daya,qia,nzs,nmet,day,qi)
      call intx (daya,qla,nzs,nmet,day,ql)
      call intx (daya,txa,nzs,nmet,day,tx)
      call intx (daya,tya,nzs,nmet,day,ty)
      call intx (daya,empa,nzs,nmet,day,emp)
c
      end if
c
c  nn is the number of time steps actually taken.
c
      nn = nn + 1
c
c  Time step the model by an increment of time, dt.
c
      call pwpgo(qi,ql,emp,tx,ty,dml,dtl)
c
c
c  Apply a "background" diffusion if rkz is non-zero.
c
      if(rkz.gt.0) then
c
      call diffus (rkz,nz,dz,dt,d)
      call diffus (rkz,nz,dz,dt,s)
      call diffus (rkz,nz,dz,dt,t)
      call diffus (rkz,nz,dz,dt,u)
      call diffus (rkz,nz,dz,dt,v)
c
      end if
c
c  Store the density profile if this is the first step.
c
      if(m.eq.1) then
      do 2278 j=1,nz
      do(j) = d(j)
 2278 continue
      end if
c
c
c  This is the place to compute diagnostics.
c
c  Find max/min values on the first day.
c
      if(m.eq.1) dmlmin = 100.
      if(m.eq.1) tmin = t(1)
      if(timed.lt.1.) then
      if(t(1).lt.tmin) tmin = t(1)
      if(t(1).gt.tmax) tmax = t(1)
      spd = sqrt(u(1)**2 + v(1)**2)
      if(spd.gt.spdmax) spdmax = spd
      if(dml.lt.dmlmin) dmlmin = dml
      end if
c
c
c  Write out a few things on every iwrtth time step.
c
      if(iwrt.lt.1) go to 1339
c
      if(m.eq.1) write (6,444)
  444 format (1x,/,1x,'model solution variables are:',/,1x,
     x '   m timed    qi     ql      tx   t(1)   s(1)   u(1) ',
     x '  v(1)  dml   dtl',/,1x)
c
      if(mod(m,iwrt).ne.0) go to 1339
      write (6,133) m,timed,qi,ql,tx,t(1),s(1),
     c u(1),v(1),dml,dtl
  133 format (1x,i4,f6.2,2f7.0,f7.2,4f7.2,2f6.0)
c
 1339 continue
c
c
c  Compute the net heat flux as a check on heat conservation.
c
      qnet = qnet + dt*(qi + ql)
c
c
c  Store data on disk, logical unit 7, if you told it to.
c
      if(istor.eq.1) then
      if(mod(m,itfreq).eq.0) then
      write (7,3188) timed,qi,ql,tx,t(1),s(1),u(1),v(1),
     x dml
 3188 format (1x,20f10.3)
      end if
      end if
c
   40 continue
  499 continue
c
c
      write (6,870) tmin,tmax,dmlmin,spdmax
  870 format (1x,/,1x,'tmin, tmax, dmlmin, spdmax are',4f8.2)
c
c
c  Compute heat and potential energy diagnostics.
c  dpe will be the change in potential energy over the
c  integration period, and heat is the change in heat
c  content. heat should nearly equal the time-integrated
c  surface heat flux, qnet, if heat has been conserved.
c
      dpe = 0.
      heat = 0.
c
      do 8376 j=1,nz
      call intx(zo,to,nzs,nts,z(j),ti)
      tanom = t(j) - ti
      heat = heat + tanom*dz
      if(j.eq.1) tanoms = tanom
c
      danom = d(j) - do(j)
      dpe = dpe - g*danom*z(j)*dz

c     write (6,4465) j,do(j),d(j),danom,dpe
 4465 format (1x,'do,d,danom,dpe',i4,3f10.3,e12.2)
 8376 continue
c
      qnet = qnet/(ro*cpw)
      write (6,3877) heat, qnet, dpe
 3877 format (1x,/,1x,'heat, qnet, dpe are', 3f9.2)
c
c  End of the diagnostics.
c
      end
c
c
      subroutine pwpint (beta1,beta2,rlat,udrag)
c
c
c  This subroutine is an implementation of the Price, Weller,
c  Pinkel upper ocean model described in JGR 91, C7 8411-8427
c  July 15, 1986 (PWP). This version was coded and documented by
c  Jim Price in April 1989.
c
c  Edited on 20 September 1993 by JFP to allow for a critical
c  gradient Richardson number other than 1/4, and to implement a
c  different and a priori better means of achieving convergence
c  of gradient ri mixing (see subroutine stir for details).
c  The major difference is that the revised scheme gives a more
c  smoothly varying mixed layer depth over a diurnal cycle.
c  Edited on 14 December 1998 by JFP to clarify the bulk Ri
c  relaxation.
c
c  This model also implements an energy budget form of an
c  entrainment parameterization that is very similar to that
c  described in Price, Mooers, Van Leer, JPO 8, 4, 582-599 (and
c  references therein). This part of the model should be
c  treated as developmental only, as there are several features
c  that are arbitrary to this model (i.e., the depth of the ml
c  during times of heating can be the grid interval, dz). To
c  use this parameterization set the bulk Richardson number to
c  zero, and set em1, or em2, or em3 to non-zero. The gradient
c  Richardson number can be zero or not.
c
c
c  Please direct questions and comments to Jim Price, WHOI,
c  Woods Hole, MA, 02543, MA, tel. 508-548-1400, x2526.
c
c  All units within the model are mks, except density, d,
c  which is in sigma units.
c
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
c
c
c  nz is the number of grid points in the vertical.
c  dz is the grid interval (meters).
c  beta1 and beta2 are the longwave and shortwave extinction
c   coefficients for solar insolation (meters).
c  rlat is the latitude (deg).
c  tr, sr are the reference values of temperature and
c   salinity used in the density equation (C, and ppt).
c  udrag is the e-folding time scale of current due to
c   unresolved processes (days)  (set to infinite if > 100).
c  rg is the critical gradient richardson number (= 1/4,
c   in the pwp model).
c  rb is the critical bulk richardson number (= 0.65,
c   in the pwp model).
c  em1, em2, em3 are coefficients that multiply energy
c   source terms in the turbulent kinetic energy budget:
c   buoyancy production, shear production, and ustar**3; all
c   are zero in the pwp model, but can be used here by setting
c   rb = 0, and then at least one of the emns to non-zero.
c  There are some peculiar features of any energy budget model
c  that arise when the ocean is heated, so be careful !
c
      save d2r, hcon, ucon, f, tr, sr, dr, alpha, beta, energy
c
c  Set a few constants that will not change during the run.
c
      d2r = 2.*3.14159/360.
      g = 9.8
      ro = 1.024e3
      cpw = 4183.3
      hcon = dz*ro*cpw
      ucon = 0.
      if(udrag.lt.100.) ucon = 1./(8.64e4*udrag)
      f = 2.*7.29e-5*sin(rlat*d2r)
c
c  Evaluate alpha and beta, the thermal and haline expansion
c  coefficients used in a linear equation of state.
c  dr is the reference density.
c
      tr = t(1)
      sr = s(1)
      dr = sgt(tr,sr,gg)
      alpha = sgt(tr+.5,sr,gg) - sgt(tr-.5,sr,gg)
      beta = sgt(tr,sr+.5,gg) - sgt(tr,sr-.5,gg)
c
c  Compute the absorption profile; absrb(j) is the fraction
c  of the solar insolation that will be absorbed within
c  the grid cell j.
c
      call absorb(beta1,beta2)
c
c  Finished with initialization.
c
c
      return
c
c
c  Come here to step ahead by an amount dt.
c
      entry pwpgo (qi,ql,emp,tx,ty,dml,dtl)
c
c  Input to the subroutine are the following:
c
c  qi and ql are the insolation and heat loss (qi > 0,
c  ql < 0 almost always)  (W/m**2).
c  emp is evaporation minus precipitation (m/sec)
c  (emp > if evaporation exceeds precipitation).
c  tx and ty are the east and north stress components (Pa).
c
c  Apply heat and fresh water fluxes to the top most grid cell.
c
      q1 = qi*absrb(1) + ql
      t(1) = t(1) + dt*q1/(dz*ro*cpw)
c
c  Absorb solar radiation at depth.
c
      do 70 j=2,nz
      t(j) = t(j) + qi*absrb(j)*dt/(dz*ro*cpw)
   70 continue
      s(1) = s(1) + s(1)*emp*dt/dz
c
c  Compute the density (sigma units), and relieve static
c  instability, if it occurs.
c
      do 71 j=1,nz
      d(j) = dr + (t(j) - tr)*alpha + (s(j) - sr)*beta
   71 continue
c
      call mldep(ds1)
      call si(ml,sipe)
      if(ml.gt.1) call mix5(ml)
      call mldep(ds2)
c
c  At this point the density proifile should be statically stable.
c
c  Time step the momentum equation.
c
c  Rotate the current throughout the water column through an
c  angle equal to inertial rotation for half a time step.
c
      ang = f*dt/2.
      sa = sin(ang)
      ca = cos(ang)
c
      do 32 j=1,nz
      call rot(u(j),v(j),sa,ca)
   32 continue
c
c  Apply the wind stress to the mixed layer as it now exists.
c
c  Find the ml depth.
c
      call mldep(dml)
      ml = int(dml/dz)
c
      du = (tx/(dml*ro))*dt
      dv = (ty/(dml*ro))*dt
      do 33 j=1,ml
      u(j) = u(j) + du
      v(j) = v(j) + dv
   33 continue
c
c  Apply drag to the current (this is a horrible parameterization of
c  inertial-internal wave dispersion).
c
      if(ucon.gt.1.e-10) then
      do 37 j = 1,nz
      u(j) = u(j)*(1. - dt*ucon)
      v(j) = v(j)*(1. - dt*ucon)
   37 continue
      end if
c
c  Rotate another half time step.
c
      do 34 j = 1,nz
      call rot(u(j),v(j),sa,ca)
   34 continue
c
c  Finished with the momentum equation for this time step.
c
c  Cause the mixed-layer to deepen.
c
c
      if(rb.le.0.00001) go to 61
c
c  Come here to do the bulk Richardson number instability
c  form of mixing (as in PWP).
c
      nzm = nz - 1
c
      rvc = rb
      do 80 j = ml,nzm
      delr = (d(j+1) - d(1))/ro
      ds2 = (u(j+1) - u(1))**2 + (v(j+1) - v(1))**2
      ds2 = ds2 + 1.e-8
      h1 = float(j)*dz
      rv = g*delr*h1/ds2
      if(rv.gt.rvc) go to 81
c
      jp = j + 1
      call mix5(jp)
c
   80 continue
   81 continue
      go to 65
c
c
   61 continue
c
c
      if(em1.le.0..and.em2.le.0..and.em3.le.0.) go to 65
c
c  Come here to mix according to an energy balance closure;
c  this implementation follows Price, Mooers, Van Leer, JPO,
c  8, 582-599, July 1978. This is not part of the PWP model
c  per se, and has some rather unusual properties, so watch out !
c
c
      call mldep(dml)
      ml = dml/dz
      taumag = sqrt(tx**2 + ty**2)
      ustar = sqrt(taumag/ro)
      ustarc = taumag*ustar*dt
c
c  Compute the energy available to do entrainment.
c
      energy = energy + em3*ustarc + em1*sipe
c
      delkes = 0.
      penes = 0.
c
      nzm = nz - 1
      do 60 l=ml,nzm
c
c  Determine how much energy is required to entrain the next
c  level below the base of the mixed-layer.
c
      call pe(l+1,pene)
      if(energy.lt.pene) go to 65
      if(energy.gt.pene) then
c
c  Come here if sufficient energy is avaliable to entrain at
c  least one level into the mixed layer.
c
c  Compute the change in kinetic energy due to entrainment.
c
      eke1 = eke(l+1)
      call mix5(l+1)
      eke2 = eke(l+1)
      delke = eke1 - eke2
c
      penes = penes + pene
      delkes = delkes + delke
c
c  Add and subtract energy losses and gains due to entrainment.
c
      energy = energy + em2*delke - pene
c
c
      end if
   60 continue
c
c  End of the energy budget method for estimation entrainment.
c
   65 continue
c
c  Relieve gradient Richardson number instability if it occurs.
c
      if(rg.gt.0.) call rimix(nmix,kcd)
c
      call mldep(dml)
      dtl = dz*float(kcd)
c
      return
      end
c
      subroutine pe(m,pene)
c
c  Compute the energy required to entrain density level
c  d(m) into a mixed-layer above.
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, rg, rb, em1, em2, em3
c
      dr = d(m) - d(m-1)
      ho2 = dz*float(m)/2.
      pene = g*dz*ho2*dr
      return
      end
c
      subroutine rimix(nmix,kcd)
c
c  This subroutine performs the gradient Richardson number
c  relaxation by mixing adjacent cells just enough to bring
c  them to a new richardson number = rnew (defined in
c  subroutine stir). rg is the critical gradient Richardson number
c  at which mixing is presumed to begin. rg = 0.25 is a nominal
c  value.
c
c  All of u, v, d, t, s, dz, nz are exactly as in the pwpri subroutine.
c  nmix is the number of iterations required to achieve the relaxation,
c  kcd is the deepest grid level reached (the base of the transition
c  layer).
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cwp, rg, rb, em1, em2, em3
c
      dimension r(500)
c
      rc = rg
      nzm = nz - 1
      kcd = 0
      nmix = 0
c
      j1 = 1
      j2 = nzm
   10 continue
c
c  Compute the gradient Richardson number, taking care
c  to avoid dividing by zero in the mixed-layer. The numerical
c  values of the minimum allowable density and velocity differences
c  are entirely arbitary, and should not effect the calculations
c  (except that on some occasions thay evidently have!)
c
      do 1 j = j1,j2
      dd = d(j+1) - d(j)
      if(dd.lt.1.e-3) dd = 1.e-3
      dv = (u(j+1) - u(j))**2 + (v(j+1) - v(j))**2
      if(dv.lt.1.e-6) dv = 1.e-6
      r(j) = g*dz*dd/(dv*ro)
    1 continue
c
c  Find the smallest value of r in profile.
c
      rs = r(1)
      js = 1
      do 2 j=2,nzm
      if(r(j).lt.rs) go to 3
      go to 2
    3 continue
      rs = r(j)
      js = j
    2 continue
c
c  Check to see whether the smallest r is critical or not.
c
      if(rs.gt.rc) go to 99
c
c  Mix the cells js and js+1 that had the smallest Richardson number.
c
      if(js.ge.kcd) kcd = js + 1
      call stir(rc,rs,t,js)
      call stir(rc,rs,s,js)
      call stir(rc,rs,d,js)
      call stir(rc,rs,u,js)
      call stir(rc,rs,v,js)
      nmix = nmix + 1
c
c  Recompute the Richardson number over the part of the profile
c  that has changed.
c
      j1 = js - 2
      if(j1.lt.1) j1 = 1
      j2 = js + 2
      if(j2.gt.nzm) j2 = nzm
      go to 10
c
   99 continue
      return
      end
c
c
      subroutine stir (rc,r,a,j)
c
c
c  This subroutine mixes cells j and j+1 just enough so that
c  the Richardson number after the mixing is brought up to
c  the value rnew. In order to have this mixing process
c  converge, rnew must exceed the critical value of the
c  richardson number where mixing is presumed to start. if
c  r critical = rc = 0.25 (the nominal value), and r = 0.20, then
c  rnew = 0.3 would be reasonable. If r were smaller, then a
c  larger value of rnew - rc is used to hasten convergence.
c
c  This subroutine was modified by JFP in Sep 93 to allow for an
c  aribtrary rc and to achieve faster convergence.
c
c  rc is the value of the critical gradient Richardson number,
c  r is the Richardson number before this mixing,
c  a is the array to be mixed,
c  j and j+1 define the cell (grid level) where r was found to
c   be critical.
c
      dimension a(500)
c
c  Set the convergence parameter.
c
      rcon = 0.02 + (rc - r)/2.
      rnew = rc + rcon/5.               !  an arbitrary change
c
      f = 1. - r/rnew
      da = (a(j+1) - a(j))*f/2.
      a(j+1) = a(j+1) - da
      a(j) = a(j) + da
      return
      end
c
c
c
      subroutine mldep(dml)
c
c  This subroutine scans through the density array d to find
c  the depth of the surface mixed-layer. The degree of density
c  homogeneity, deps, is arbitary and should not effect the
c  results.
c
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
c
      deps = 1.e-4
      nzm = nz - 1
      do 1 j=1,nzm
      dd = abs(d(j+1) - d(j))
      if(dd.gt.deps) go to 2
    1 continue
    2 continue
      dml = float(j)*dz
      return
      end
c
c
      subroutine mix1(b,k)
c
c  This subroutine homogenizes the array b down to level k.
c
      dimension b(500)
      bs = 0.
      do 1 j=1,k
      bs = bs + b(j)
    1 continue
      bs = bs/float(k)
      do 2 j=1,k
      b(j) = bs
    2 continue
      return
      end
c
c
      subroutine mix5(k)
c
c  This subroutine mixes arrays t, s, u, v, d down to level k.
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
c
      ts = 0.
      ss = 0.
      us = 0.
      vs = 0.
      ds = 0.
      do 1 j=1,k
      ts = ts + t(j)
      ss = ss + s(j)
      us = us + u(j)
      vs = vs + v(j)
      ds = ds + d(j)
    1 continue
      x = float(k)
      do 2 j=1,k
      u(j) = us/x
      v(j) = vs/x
      t(j) = ts/x
      s(j) = ss/x
      d(j) = ds/x
    2 continue
      return
      end
c
c
      subroutine rot(u,v,sa,ca)
c
c  This subroutine rotates the vector (u,v) through an
c  angle whose sine and cosine are sa and ca.
c
      u0 = u
      v0 = v
      u = u0*ca + v0*sa
      v = v0*ca - u0*sa
      return
      end
c
c
      subroutine si(ml,sipe)
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
c
      dimension di(500)
c
c  Find and relieve static instability that may occur in the
c  density array d. This simulates free convection.
c  ml is the depth of the surface mixed layer after adjustment,
c  sipe is the change in potential energy due to free convection.
c
      do 10 j=1,nz
      di(j) = d(j)
   10 continue
c
      mml = 1
      do 1 j=2,nz
      if(d(j).gt.d(j-1)) go to 2
      if(d(j).lt.d(j-1)) then
      mml = j
      call mix1(d,mml)
      end if
    1 continue
    2 continue
      ml = mml
c
      nj = mml + 3
      do 30 j=1,nj
   30 continue
c
c
      pesum = 0.
      if(ml.ge.2) then
c
c  Compute the change in potential energy.
c
      do 11 j= ml,1,-1
      z = float(j)*dz - dz/2.
      pesum = pesum + (d(j) - di(j))*g*z*dz
   11 continue
      end if
      sipe = pesum
c
      return
      end
c
c
      function eke(k)
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
c
c  Compute the kinetic energy in the u,v profiles from the surface
c  down to level k.
c
      ekes = 0.
      do 1 j=1,k
      ekes = ekes + u(j)**2 + v(j)**2
    1 continue
      eke = ekes*ro*dz/2.
      return
      end
c
c
      subroutine absorb(beta1,beta2)
c
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
c
c  Compute solar radiation absorption profile. This
c  subroutine assumes two wavelengths, and a double
c  exponential depth dependence for absorption.
c
c  Subscript 1 is for red, non-penetrating light, and
c  2 is for blue, penetrating light. rs1 is the fraction
c  assumed to be red.
c
      abss = 0.
      rs1 = 0.6
      rs2 = 1.0 - rs1
c
      do 3 j=1,nz
      absrb(j) = 0.
    3 continue
c
      do 1 j=1,nz
      z1 = float(j-1)*dz
      z2 = z1 + dz
c
      z1b1 = z1/beta1
      z2b1 = z2/beta1
      z1b2 = z1/beta2
      z2b2 = z2/beta2
c
c  The checks below are to avoid math overflows.
c
      if(z2b1.lt.70.) absrb(j) = absrb(j) +
     x  rs1*(exp(-z1b1) - exp(-z2b1))
c
      if(z2b2.lt.70.) absrb(j) = absrb(j) +
     x  rs2*(exp(-z1b2) - exp(-z2b2))
c
c  Make sure that ab(z) integrates to 1 at large depth.
c
      abss = abss + absrb(j)
c      if(j.lt.20) write (6,2) j,z2,absrb(j),abss
c    2 format (1x,'absorption',i5,3f8.3)
    1 continue
      return
      end
c
c
      subroutine intx(r,s,nzs,n,ri,si)
c
c  This subroutine interploates the array s to find the value
c  at the point ri, where r is the independent variable.
c  The interpolation is linear, and assumes that r increases
c  montonically, and that ri is within the range of r. If
c  ri is not within the range of r, si is set equal to the
c  appropriate endpoint of s, and a warning is printed.

c
      dimension r(1000), s(1000)
c
c  Check for out of range ri.
c
      if(ri.lt.r(1)) then
      si = s(1)
      write (6,20) ri, r(1)
  20  format (1x,/,1x,'warning, ri was below r(1) in call to intx',/,
     c 'ri and r(1) were', 2e12.4)
      end if
c
c
      if(ri.gt.r(n)) then
      si = s(n)
      write (6,21) ri, r(n)
  21  format (1x,/,1x,'warning, ri was above r(n) in call to intx',/,
     c 'ri and r(n) were', 2e12.4)
      end if
c
c  OK, data are within bounds for an interpolation.
c
      nm = n - 1
      do 2 k=1,nm
      if(ri.ge.r(k).and.ri.le.r(k+1)) go to 1
    2 continue
c
      go to 9
c
    1 continue
      si = (s(k)*(r(k+1) - ri) + s(k+1)*(ri - r(k)))/
     x (r(k+1) - r(k))
c
c
    9 continue
      return
      end
c
c
c
      subroutine diffus (rkz,nz,dz,dt,a)
      dimension a(1000), work(1000)
c
c
c  This subroutine applies a simple diffusion
c  operation to the array a. It leaves the endpoints
c  unchanged (assumes nothing about the
c  boundary conditions).
c
c
      dconst = dt*rkz/dz**2
      nzm = nz - 1
c

      do 1 j=2,nzm
      work(j) = dconst*(a(j-1) + a(j+1) - 2.*a(j))
    1 continue
c
      do 2 j=2,nzm
    2 a(j) = a(j) + work(j)
c
      return
      end
c

c
      function sg0(s)
c
c  A sigma-0 subroutine neede by the sigma-t subroutine;
c  taken from seaprop.
c
c  sigma-0 knudsen
c
      sg0 = ((6.76786136e-6*s-4.8249614e-4)*s+0.814876577)*s
     x -0.0934458632
      return
      end
c
c
      function sgt(t,s,sg)
c
c  A sigma-t subroutine taken from seaprop;
c  sigma-t knudsen
c
      sg = sg0(s)
   20 sgt = ((((-1.43803061e-7*t-1.98248399e-3)*t-0.545939111)*t
     x +4.53168426)*t)/(t+67.26)+((((1.667e-8*t-8.164e-7)*t
     x +1.803e-5)*t)*sg+((-1.0843e-6*t+9.8185e-5)*t-4.7867e-3)*t
     x +1.0)*sg
      return
      end
