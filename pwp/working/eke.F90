function eke(k)
!
common t(500), s(500), u(500), v(500), d(500), absrb(500),
1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
!
!  Compute the kinetic energy in the u,v profiles from the surface
!  down to level k.
!
ekes = 0.
do 1 j=1,k
ekes = ekes + u(j)**2 + v(j)**2
1 continue
eke = ekes*ro*dz/2.
return
end
