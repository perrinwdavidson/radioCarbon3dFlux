subroutine diffus (rkz,nz,dz,dt,a)
dimension a(1000), work(1000)
!
!
!  This subroutine applies a simple diffusion
!  operation to the array a. It leaves the endpoints
!  unchanged (assumes nothing about the
!  boundary conditions).
!
!
dconst = dt*rkz/dz**2
nzm = nz - 1
!

do 1 j=2,nzm
work(j) = dconst*(a(j-1) + a(j+1) - 2.*a(j))
1 continue
!
do 2 j=2,nzm
2 a(j) = a(j) + work(j)
!
return
end
