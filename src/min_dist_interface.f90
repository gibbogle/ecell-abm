! Interface to min_dist() in ellipsoid_dist.dll

interface

	subroutine min_dist(a1,b1,centre1,orient1,a2,b2,centre2,orient2,s1,s2,r1,r2,d) BIND(C,NAME='min_dist')
	use :: iso_c_binding
	real(c_double),value :: a1, b1, a2, b2
	real(c_double) :: centre1(*), orient1(*), centre2(*), orient2(*)
	real(c_double) :: s1, s2, r1, r2, d
	end subroutine

end interface