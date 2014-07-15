program main
use, intrinsic :: iso_c_binding
implicit none
real(8) :: aval1, bval1, centre1(3), orient1(3), aval2, bval2, centre2(3), orient2(3), &
            s1, s2, r1, r2, d
integer :: i

interface

	subroutine min_dist(aval1,bval1,centre1,orient1,aval2,bval2,centre2,orient2,s1,s2,r1,r2,d) BIND(C,NAME='min_dist')
	use :: iso_c_binding
	real(c_double),value :: aval1, bval1, aval2, bval2
	real(c_double) :: centre1(*),orient1(*),centre2(*),orient2(*)
	real(c_double) :: s1, s2, r1, r2, d
	end subroutine

end interface

aval1 = 5
bval1 = 3
centre1 = (/ 0, 0, 0 /)
orient1 = (/ 1, 0, 0 /)
aval2 = 5
bval2 = 3
centre2 = (/ 0, 10, 0 /)
orient2 = (/ 0, 1, 0 /)
s1 = 0.5
s2 = 0.5
do i = 1,100000
    s1 = 0.2
    s2 = 0.5
    call min_dist(aval1,bval1,centre1,orient1,aval2,bval2,centre2,orient2,s1,s2,r1,r2,d)
enddo
write(*,*) 's1, s2, r1, r2, d: ',s1,s2,r1,r2,d
end