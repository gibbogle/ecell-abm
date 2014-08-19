! To explore ellipsoid geometry

module ellipsoid
use global
implicit none

type line_type
	real(REAL_KIND) :: p1(3)
	real(REAL_KIND) :: p2(3)
end type

contains

!----------------------------------------------------------------------------------------
! A general point on the line connecting P1(x1,y1,z1) and P2(x2,y2,z2) is parametrised:
! x = x1 + s(x2-x1)
! y = y1 + s(y2-y1)
! z = z1 + s(z2-z1)
! where s = 0 at P1, and s = 1 at P2
! REVISION
! We can associate a radius with a point P(s) = 
!----------------------------------------------------------------------------------------
subroutine NearestPoints(line1, line2, s1, s2)
type(line_type) :: line1, line2
real(REAL_KIND) :: s1, s2
real(REAL_KIND) :: a1, b1, c1, d1, e1, f1
real(REAL_KIND) :: a2, b2, c2, d2, e2, f2
real(REAL_KIND) :: a, b, c, d, e, f

a1 = line1%p1(1)
b1 = line1%p2(1) - line1%p1(1)
c1 = line1%p1(2)
d1 = line1%p2(2) - line1%p1(2)
e1 = line1%p1(3)
f1 = line1%p2(3) - line1%p1(3)
a2 = line2%p1(1)
b2 = line2%p2(1) - line2%p1(1)
c2 = line2%p1(2)
d2 = line2%p2(2) - line2%p1(2)
e2 = line2%p1(3)
f2 = line2%p2(3) - line2%p1(3)
a = b1*b1 + d1*d1 + f1*f1
b = -(b1*b2 + d1*d2 + f1*f2)
c = b1*(a2-a1) + d1*(c2-c1) + f1*(e2-e1)
d = b
e = b2*b2 + d2*d2 + f2*f2
f = b2*(a1-a2) + d2*(c1-c2) + f2*(e1-e2)
s2 = (c*d - a*f)/(b*d - a*e)
s1 = (c - b*s2)/a
s1 = max(s1,0.0)
s1 = min(s1,1.0)
s2 = max(s2,0.0)
s2 = min(s2,1.0)
end subroutine

!----------------------------------------------------------------------------------------
! P1 is the "nearest point" on line1, P2 is the "nearest point" on line2
! d12 is the distance from P1 to P2
! v is the unit vector from P1 to P2
!----------------------------------------------------------------------------------------
subroutine ContactAngles(ell1, ell2, s1, s2, theta1, theta2, v, d12)
type(cell_type) :: ell1, ell2
real(REAL_KIND) :: s1, s2, theta1, theta2, v(3)
real(REAL_KIND) :: p1(3), p2(3), d12, cosa

!v = line2%p1 + s2*(line2%p2 - line2%p1) - line1%p1 - s1*(line1%p2 - line1%p1)
p1 = ell1%centre + (s1-0.5)*ell1%a*ell1%orient
p2 = ell2%centre + (s2-0.5)*ell2%a*ell2%orient
v = p2 - p1
d12 = sqrt(sum(v**2))
v = v/d12
cosa = dot_product(v,ell1%orient)
theta1 = acos(cosa)
cosa = -dot_product(v,ell2%orient)
theta2 = acos(cosa)

end subroutine

!----------------------------------------------------------------------------------------
! For an ellipse with axis parameters(a,b), a point P parametrised by s on the long axis,
! computes the distance from P to the ellipse along a line at angle theta to the long axis.
!
!                                   \
!                                   /\
!                                  /
!                               d /
!                                / theta
!   ----------------------------/-------
!  -a                          P(s)    a
!----------------------------------------------------------------------------------------
subroutine InsideDistance(a, b, s, theta, d)
real(REAL_KIND) :: a, b, s, theta, d
real(REAL_KIND) :: x, y, cosa, sina
real(REAL_KIND) :: aa, bb, cc, dd, d1, d2, xb, yb

x = -a + 2*a*s
cosa = cos(theta)
sina = sin(theta)
aa = (cosa/a)**2 + (sina/b)**2
bb = 2*x*cosa/(a*a)
cc = (x/a)**2 - 1
dd = sqrt(bb*bb - 4*aa*cc)
d1 = (-bb - dd)/(2*aa)
d2 = (-bb + dd)/(2*aa)
if (d1 >= 0) then
	d = d1
else
	d = d2
endif
!write(*,*) 'd1, d2: ',d1,d2
!xb = x + d2*cosa
!yb = d2*sina
!write(*,*) 'xb, yb: ',xb,yb,(xb/a)**2+(yb/b)**2
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine NormaliseOrient(ell)
type(cell_type) :: ell
real(REAL_KIND) :: d

!d = sqrt(ell%orient(1)**2 + ell%orient(2)**2 + ell%orient(3)**2)
d = sqrt(sum(ell%orient(:)**2))
ell%orient = ell%orient/d
end subroutine

!----------------------------------------------------------------------------------------
! The peak attractive force = 1, at delta = e (e.g. = 2).  
! With these parameters:
! a = 3, b = 2, c = 5, e = 3.0, g = e/2.5
! an equal repulsive force occurs at approximately delta = -2
! (See ellipsoid_forces.xlsx)
!----------------------------------------------------------------------------------------
subroutine CellContactForce(delta, s1, s2, F)
real(REAL_KIND) :: delta, s1, s2, F
real(REAL_KIND) :: d, t1, t2, tsum, afactor

if (Fcase == 1) then
    if (delta > 2*FP1%e) then
        F = 0
    elseif (delta > 0) then
	    afactor = 1
	    if (POLARITY) then
		    t1 = abs(s1-0.5)
		    t2 = abs(s2-0.5)
		    tsum = t1 + t2
		    if (tsum < 0.5) then
			    afactor = 2*(1 + cos(tsum*PI))
		    else
			    tsum = 1 - tsum
			    afactor = 2*(1 + cos(tsum*PI))
		    endif
	    endif
        F = afactor*exp(-((delta-FP1%e)/FP1%g)**2)      ! attraction
    else
        d = delta
        if (d < -FP1%a*0.9999) d = -FP1%a*0.9999
        F = FP1%c*(-FP1%b/(FP1%a**2-d**2) + FP1%b/FP1%a**2) ! repulsion
    endif
elseif (Fcase == 2) then
    ! g = amplitude of Fa
    ! e = value of delta at the peak of Fa
    ! a = multiplying factor of Fr: when delta = -b, F = a
    ! b = scaling factor of delta
    if (delta > 2*FP2%e) then
        F = 0
    elseif (delta >= 0) then
        F = (FP2%g/2)*(cos((delta-FP2%e)*PI/FP2%e) + 1)      ! attraction
    else
        F = -FP2%a*(-delta/FP2%b)**1.5
    endif
endif    
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine CellInteraction1(ell1, ell2)
type(cell_type) :: ell1, ell2
real(REAL_KIND) :: tol
real(REAL_KIND) :: s1, s2
real(REAL_KIND) :: delta
logical :: incontact
real(REAL_KIND) :: p1(3), p2(3)		! Nearest points on main axes of ell1 and ell2
real(REAL_KIND) :: r1(3), r2(3)	    ! Displacement of contact points on ell1 and ell2 relative to centres
real(REAL_KIND) :: rad1, rad2, vamp, v(3), famp, F(3), M1(3), M2(3)
integer :: res

!write(nflog,'(8f8.4)') ell1%a,ell1%b,ell1%centre,ell1%orient,ell2%a,ell2%b,ell2%centre,ell2%orient

call min_dist(ell1%a,ell1%b,ell1%centre,ell1%orient,ell2%a,ell2%b,ell2%centre,ell2%orient,tol,s1,s2,rad1,rad2,delta,res)
! delta is the cell separation at the "closest" points
incontact = (delta < 2*FP2%e)	! there is a force to compute in the direction given by v (P1 -> P2)
if (incontact) then
!    write(nflog,'(a,5f6.2)') 's1, s2, rad1, rad2, delta: ',s1,s2,rad1,rad2,delta
    p1 = (2*s1-1)*ell1%a*ell1%orient    ! vector offset of sphere centre from ellipsoid centre
    p2 = (2*s2-1)*ell2%a*ell2%orient
    v = (ell2%centre + p2) - (ell1%centre + p1)
    vamp = sqrt(dot_product(v,v))
    v = v/vamp        ! unit vector in direction P1 -> P2
	call CellContactForce(delta, s1, s2, famp)  ! Note: famp > 0 ==> attraction
	famp = max(famp,-Flimit)
    F = famp*v          ! force in the direction of v, i.e. from P1 to P2 
	if (dbug) write(nflog,'(a,6f8.4)') 's1,s2,delta,F: ',s1,s2,delta,F
	ell1%F = ell1%F + F
!	ell2%F = ell2%F - F		! for now, apply force to the current cell only
	if (simulate_rotation) then
		! M1 = r1 x F, M2 = r2 x F      ! signs???
		r1 = p1 + rad1*v    ! vector offset of contact point from ellipsoid centre
	!	r2 = p2 - rad2*v
		call cross_product(r1,F,M1)
	!	call cross_product(r2,-F,M2)
		ell1%M = ell1%M + M1
	!	ell2%M = ell2%M + M2
	endif
endif
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine CellInteraction(a1,b1,centre1,orient1,a2,b2,centre2,orient2,s1,s2,delta,F,M1,M2,ok)
!type(cell_type) :: ell1, ell2
real(REAL_KIND) :: a1, b1, centre1(3), orient1(3)
real(REAL_KIND) :: a2, b2, centre2(3), orient2(3)
real(REAL_KIND) :: s1, s2
real(REAL_KIND) :: delta
real(REAL_KIND) :: F(3), M1(3), M2(3)
real(REAL_KIND) :: tol = 1.0d-5
logical :: ok
logical :: incontact
real(REAL_KIND) :: p1(3), p2(3)		! Nearest points on main axes of cell1 and cell2 
real(REAL_KIND) :: r1(3), r2(3)	    ! Displacement of contact points on cell1 and cell2 relative to centres
real(REAL_KIND) :: rad1, rad2, vamp, v(3), famp, dlimit
integer :: res

!write(nflog,'(8f8.4)') ell1%a,ell1%b,ell1%centre,ell1%orient,ell2%a,ell2%b,ell2%centre,ell2%orient

call min_dist(a1,b1,centre1,orient1,a2,b2,centre2,orient2,tol,s1,s2,rad1,rad2,delta,res)
!delta = 0
!res = 0
if (res /= 0 &
	.or. isnan(delta) .or. isnan(s1) .or. isnan(s2) .or. isnan(rad1) .or. isnan(rad2)) then
	write(*,'(a,i3,5f9.4)') 'Error: CellInteraction: res: ',res,s1,s2,rad1,rad2,delta
    write(*,*) 'a, b, centre, orient:'
	write(*,'(a,8f8.3)') 'Cell 1: ',a1,b1,centre1,orient1
	write(*,'(a,8f8.3)') 'Cell 2: ',a2,b2,centre2,orient2
	ok = .false.
	return
    endif
! delta is the cell separation at the "closest" points
if (Fcase == 1) then
    dlimit = 2*FP1%e
elseif (Fcase == 2) then
    dlimit = 2*FP2%e
endif
incontact = (delta < dlimit)	! there is a force to compute in the direction given by v (P1 -> P2)
F = 0
M1 = 0
M2 = 0
if (incontact) then
!    write(nflog,'(a,5f6.2)') 's1, s2, rad1, rad2, delta: ',s1,s2,rad1,rad2,delta
    p1 = (2*s1-1)*a1*orient1    ! vector offset of sphere centre from ellipsoid centre
    p2 = (2*s2-1)*a2*orient2
    v = (centre2 + p2) - (centre1 + p1)
    vamp = sqrt(dot_product(v,v))
    v = v/vamp        ! unit vector in direction P1 -> P2
	call CellContactForce(delta, s1, s2, famp)  ! Note: famp > 0 ==> attraction
!	famp = max(famp,-Flimit)
    F = famp*v          ! force in the direction of v, i.e. from P1 to P2
	if (dbug) write(nflog,'(a,6f8.4)') 's1,s2,delta,F: ',s1,s2,delta,F
!	ell1%F = ell1%F + F
!	ell2%F = ell2%F - F		! for now, apply force to the current cell only
	! M1 = r1 x F, M2 = r2 x F      ! signs???
	r1 = p1 + rad1*v    ! vector offset of contact point from ellipsoid centre
	r2 = p2 - rad2*v
	call cross_product(r1,F,M1)
	call cross_product(r2,-F,M2)
!	ell1%M = ell1%M + M1
!	ell2%M = ell2%M + M2
endif
ok = .true.
end subroutine

end module

