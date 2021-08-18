subroutine si(ml,sipe)
!
common t(500), s(500), u(500), v(500), d(500), absrb(500),
1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
!
dimension di(500)
!
!  Find and relieve static instability that may occur in the
!  density array d. This simulates free convection.
!  ml is the depth of the surface mixed layer after adjustment,
!  sipe is the change in potential energy due to free convection.
!
do 10 j=1,nz
di(j) = d(j)
10 continue
!
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
!
nj = mml + 3
do 30 j=1,nj
30 continue
!
!
pesum = 0.
if(ml.ge.2) then
!
!  Compute the change in potential energy.
!
do 11 j= ml,1,-1
z = float(j)*dz - dz/2.
pesum = pesum + (d(j) - di(j))*g*z*dz
11 continue
end if
sipe = pesum
!
return
end
