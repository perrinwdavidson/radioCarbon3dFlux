
      subroutine rot(u,v,sa,ca)
!
!  This subroutine rotates the vector (u,v) through an
!  angle whose sine and cosine are sa and ca.
!
      u0 = u
      v0 = v
      u = u0*ca + v0*sa
      v = v0*ca - u0*sa
      return
      end
