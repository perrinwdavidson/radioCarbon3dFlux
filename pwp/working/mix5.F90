subroutine mix5(k)
!
!  This subroutine mixes arrays t, s, u, v, d down to level k.
!
common t(500), s(500), u(500), v(500), d(500), absrb(500),
1 dz, dt, nz, g, ro, cpw, rg, rb, em1, em2, em3
!
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
