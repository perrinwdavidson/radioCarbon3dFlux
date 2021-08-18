      subroutine mldep(dml)
!
!  This subroutine scans through the density array d to find
!  the depth of the surface mixed-layer. The degree of density
!  homogeneity, deps, is arbitary and should not effect the
!  results.
!
!
      common t(500), s(500), u(500), v(500), d(500), absrb(500),
     1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
!
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
