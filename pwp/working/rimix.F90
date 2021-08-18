      subroutine rimix(nmix,kcd)
!
!  This subroutine performs the gradient Richardson number
!  relaxation by mixing adjacent cells just enough to bring
!  them to a new richardson number = rnew (defined in
!  subroutine stir). rg is the critical gradient Richardson number
!  at which mixing is presumed to begin. rg = 0.25 is a nominal
!  value.
!
!  All of u, v, d, t, s, dz, nz are exactly as in the pwpri subroutine.
!  nmix is the number of iterations required to achieve the relaxation,
!  kcd is the deepest grid level reached (the base of the transition
!  layer).
!
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cwp, rg, rb, em1, em2, em3
!
      dimension r(500)
!
      rc = rg
      nzm = nz - 1
      kcd = 0
      nmix = 0
!
      j1 = 1
      j2 = nzm
   10 continue
!
!  Compute the gradient Richardson number, taking care
!  to avoid dividing by zero in the mixed-layer. The numerical
!  values of the minimum allowable density and velocity differences
!  are entirely arbitary, and should not effect the calculations
!  (except that on some occasions thay evidently have!)
!
      do 1 j = j1,j2
      dd = d(j+1) - d(j)
      if(dd.lt.1.e-3) dd = 1.e-3
      dv = (u(j+1) - u(j))**2 + (v(j+1) - v(j))**2
      if(dv.lt.1.e-6) dv = 1.e-6
      r(j) = g*dz*dd/(dv*ro)
    1 continue
!
!  Find the smallest value of r in profile.
!
      rs = r(1)
      js = 1
      do 2 j=2,nzm
      if(r(j).lt.rs) go to 3
      go to 2
    3 continue
      rs = r(j)
      js = j
    2 continue
!
!  Check to see whether the smallest r is critical or not.
!
      if(rs.gt.rc) go to 99
!
!  Mix the cells js and js+1 that had the smallest Richardson number.
!
      if(js.ge.kcd) kcd = js + 1
      call stir(rc,rs,t,js)
      call stir(rc,rs,s,js)
      call stir(rc,rs,d,js)
      call stir(rc,rs,u,js)
      call stir(rc,rs,v,js)
      nmix = nmix + 1
!
!  Recompute the Richardson number over the part of the profile
!  that has changed.
!
      j1 = js - 2
      if(j1.lt.1) j1 = 1
      j2 = js + 2
      if(j2.gt.nzm) j2 = nzm
      go to 10
!
   99 continue
      return
      end
