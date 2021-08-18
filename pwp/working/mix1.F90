subroutine mix1(b,k)
!
!  This subroutine homogenizes the array b down to level k.
!
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
