! To solve cell motion using RK (r8_rkf45)

module solve
use ellipsoid
use rkf45

implicit none

real(REAL_KIND), allocatable :: F(:,:,:)    ! cell-cell forces
real(REAL_KIND), allocatable :: M(:,:,:)    ! cell-cell moments
!real(REAL_KIND), allocatable :: cp(:,:,:)  ! offset of the contact points from centres
!integer, allocatable :: nbrs(:)             ! number of cell neighbours
!integer, allocatable :: nbrlist(:,:)        ! cell neighbours

contains

!------------------------------------------------------------------------------------------
! To construct the Jacobian matrix:
!
! row 1: df(1)/dx(1) df(1)/dx(2) df(1)/dx(3) ... df(1)/dx(n)
! row 2: df(2)/dx(1) df(2)/dx(2) df(2)/dx(3) ... df(2)/dx(n)
! ...
! row n: df(n)/dx(1) df(n)/dx(2) df(n)/dx(3) ... df(n)/dx(n)
!
! where dx(i)/dt = f(i)
!------------------------------------------------------------------------------------------
subroutine JacobianSlow(Jac,x,n)
integer :: n
real(REAL_KIND) :: x(n), Jac(n,n)
real(REAL_KIND), allocatable :: xx(:), f(:), ff(:)
real(REAL_KIND) :: t, dx
integer :: i, j

allocate(xx(n))
allocate(f(n))
allocate(ff(n))

dx = 0.001
call fderiv(t,x,f)
xx = x
do j = 1,n
	xx(j) = x(j) + dx
	call fderiv(t,xx,ff)
	do i = 1,n
		Jac(i,j) = (ff(i) - f(i))/dx
	enddo
	xx(j) = x(j)
enddo
deallocate(xx)
deallocate(f)
deallocate(ff)

end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
logical function inlist(i,list,n)
integer :: i, n, list(:)
integer :: k

do k = 1,n
	if (list(k) == i) then
		inlist = .true.
		return
	endif
enddo
inlist = .false.
end function

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
subroutine remove(i,list,n)
integer :: i, n, list(:)
integer :: k

do k = 1,n
	if (list(k) == i) then
		list(k) = -i
		return
	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------
! Note: only those forces (F, M) on cells that are neighbours of the cell associated with
! variable j are changed when xx(j) is varied.  Therefore it is necessary to compute
! only a small subset of the cell-cell forces.  The other F, M are zero.
! df(i)/dx(j) is non-zero (variable i belonging to cell i1) only in the case that cell i1
! owns variable j or has a neighbour cell i2 that owns variable j.  This means that when
! looking at the variation of x(j), we need only to compute interactions for cell i1 and
! neighbours cells of i1 that own variable j. 
! We can precompute the list of cells that need to be considered for variable j.  This 
! could perhaps be done in the j=0 pass, but for clarity maybe not.
!------------------------------------------------------------------------------------------
subroutine JacobianFast(Jac, x, n, np, ok)
real(REAL_KIND) :: Jac(n,n), x(n)
integer :: n, np
logical :: ok
real(REAL_KIND), allocatable :: xx(:), f0(:), f1(:), Fbase(:,:,:), Mbase(:,:,:)
real(REAL_KIND) :: t, dx
integer :: i, j, jj
integer :: i1, i2, jn, k1, k2, k, nzero, ntzero, nfeval, kpar=0
real(REAL_KIND) :: a1, b1, centre1(3), orient1(3), a2, b2, centre2(3), orient2(3), s1, s2, delta
real(REAL_KIND) :: FF(3), MM1(3), MM2(3), Fsum(3), Msum(3), R
type(cell_type), pointer :: p1, p2
integer, allocatable :: connected(:,:), nconnected(:)
logical :: hit
logical :: clean_list = .true.

allocate(xx(n))
allocate(f0(n))
allocate(f1(n))
allocate(connected(n,MAX_NBRS+1))
allocate(nconnected(n))
allocate(Fbase(ncells,ncells,3))
if (np == 6) allocate(Fbase(ncells,ncells,3))

Jac(1:n,1:n) = 0

! This can be improved by doing it in the j=0 pass, since that will enable dropping of cells that are 
! neighbours but too far away to interact with cell i1.
nconnected = 0
do i1 = 1,ncells
	p1 => cell_list(i1)
	k1 = (i1-1)*np
	do k = 1,np
		j = k1 + k
		if (.not.inlist(i1,connected(j,:),nconnected(j))) then
			nconnected(j) = nconnected(j) + 1
			connected(j,nconnected(j)) = i1
		endif
	enddo
	do jn = 1,p1%nbrs
		i2 = p1%nbrlist(jn)
		k2 = (i2-1)*np
		do k = 1,np
			j = k1 + k
			if (.not.inlist(i2,connected(j,:),nconnected(j))) then
				nconnected(j) = nconnected(j) + 1
				connected(j,nconnected(j)) = i2
			endif
		enddo
	enddo
enddo

!do j = 1,n
!	write(nflog,'(19i4)') j,nconnected(j),connected(j,	1:nconnected(j))
!enddo
!write(*,*)
xx = x
dx = 0.001
nfeval = 0
ntzero = 0
F = 0
if (np == 6) M = 0

do j = 0,n		! column index, except for j=0 which evaluates the current f(:)
	if (j /= 0) then
		xx(j) = x(j) + dx
		F = Fbase
		if (np == 6) M = Mbase
	endif
	nzero = 0
	do i1 = 1,ncells
		p1 => cell_list(i1)
		a1 = p1%a
		b1 = p1%b
		k1 = (i1-1)*np
		centre1 = xx(k1+1:k1+3)
		if (np == 6) then
			orient1 = xx(k1+4:k1+6)
		else
			orient1 = p1%orient
		endif
		if (np == 6) then
			M(i1,1:ncells,:) = 0
		endif
		do jn = 1,p1%nbrs
			i2 = p1%nbrlist(jn)
			if (j > 0) then
				if (.not.inlist(i2,connected(j,:),nconnected(j))) cycle
			endif
			p2 => cell_list(i2)
			a2 = p2%a
			b2 = p2%b
			k2 = (i2-1)*np
			centre2 = xx(k2+1:k2+3)
			if (np == 6) then
				orient2 = xx(k2+4:k2+6)
			else
				orient2 = p2%orient
			endif
!			s1 = 0.5	! These are the initial guesses
!			s2 = 0.5
			s1 = s1s2(i1,i2,1)
			s2 = s1s2(i1,i2,2)
			call CellInteraction(a1,b1,centre1,orient1,a2,b2,centre2,orient2,s1,s2,delta,FF,MM1,MM2,ok)
			nfeval = nfeval + 1
			if (FF(1)==0 .and. FF(2)==0 .and. FF(3)==0) then
				nzero = nzero + 1
				ntzero = ntzero + 1
!				write(*,*) 'FF=0: ',j,i1,i2,nzero,ntzero,nfeval
				
				if (clean_list) then
					do k = 1,np
						jj = k1 + k
						call remove(i2,connected(jj,:),nconnected(jj))
					enddo
					do k = 1,np
						jj = k2 + k
						call remove(i1,connected(jj,:),nconnected(jj))
					enddo
				endif
	
			endif
			if (.not.ok) then
				write(*,*) 'istep, i1, i2: ',istep,i1,i2
				return
			endif
			s1s2(i1,i2,1) = s1
			s1s2(i1,i2,2) = s2
			s1s2(i2,i1,1) = s1
			s1s2(i2,i1,2) = s2
			F(i1,i2,:) = FF
			F(i2,i1,:) = -F(i1,i2,:)
			if (np == 6) then
		        M(i1,i2,:) = MM1
		        M(i2,i1,:) = MM2
		    endif
		enddo
	enddo
	if (j == 0) then
		Fbase = F
		if (np == 6) Mbase = M
	endif
	
	f1 = 0
	do i1 = 1,ncells
		Fsum = 0
		Msum = 0
		p1 => cell_list(i1)
		do jn = 1,p1%nbrs
			i2 = p1%nbrlist(jn)
			Fsum = Fsum + F(i1,i2,:)
			if (np == 6) then
		        Msum = Msum + M(i1,i2,:)
		    endif
		enddo
		k1 = (i1-1)*np
		f1(k1+1:k1+3) = Fsum/Fdrag
		if (np == 6) then
			call cross_product(Msum/Mdrag,orient1,f1(k1+4:k1+6))
		endif
	enddo
	if (j == 0) then
		f0 = f1
	else
		do i = 1,n
			Jac(i,j) = (f1(i) - f0(i))/dx
		enddo
		xx(j) = x(j)	
	endif
enddo

deallocate(xx)
deallocate(f0)
deallocate(f1)
deallocate(connected)
deallocate(nconnected)
deallocate(Fbase)
if (np == 6) deallocate(Mbase)

end subroutine

!------------------------------------------------------------------------------------------
! x(:) holds the position (centre) and orientation (orient) of each cell.
! The neighbour list nbrlist(:,:) holds each cell's neighbours, assumed not to change
! over the duration of this solving step.
! The force between two cells is stored in F(:,:,:), and the contact points on the two cells
! is in cp(:,:,:).  If the force between cells i1 and i2 has already been computed: F(i1,i2,:),
! then F(i2,i1,:) is set = -F(i1,i2,:) and is not recomputed later when cell i2 is treated.
!
!Note: could take account of growth of a, b with t!
!------------------------------------------------------------------------------------------
subroutine fderiv(t,x,xp)
real(REAL_KIND) :: t, x(*), xp(*)
logical :: ok
integer :: i1, i2, j, k1, k2, kpar=0
real(REAL_KIND) :: a1, b1, centre1(3), orient1(3), a2, b2, centre2(3), orient2(3), s1, s2, delta
real(REAL_KIND) :: FF(3), MM1(3), MM2(3), Fsum(3), Msum(3), R, mamp, vm(3), dangle 
type(cell_type), pointer :: p1, p2
integer :: np = 3	! 6

F = 0
do i1 = 1,ncells
    p1 => cell_list(i1)
    a1 = p1%a
    b1 = p1%b
    k1 = (i1-1)*np
    centre1 = x(k1+1:k1+3)
    if (np == 6) then
		orient1 = x(k1+4:k1+6)
	else
		orient1 = p1%orient
	endif
    do j = 1,p1%nbrs
        i2 = p1%nbrlist(j)
        if (F(i1,i2,1) /= 0) cycle		!??????????
        p2 => cell_list(i2)
        a2 = p2%a
        b2 = p2%b
        k2 = (i2-1)*np
        centre2 = x(k2+1:k2+3)
	    if (np == 6) then
			orient2 = x(k2+4:k2+6)
		else
			orient2 = p2%orient
		endif
!		s1 = 0.5	! initial guesses
!		s2 = 0.5
		s1 = s1s2(i1,i2,1)
		s2 = s1s2(i1,i2,2)
        call CellInteraction(a1,b1,centre1,orient1,a2,b2,centre2,orient2,s1,s2,delta,FF,MM1,MM2,ok)
        if (.not.ok) then
            write(*,*) 'istep, i1, i2: ',istep,i1,i2
            return
        endif
		s1s2(i1,i2,1) = s1
		s1s2(i1,i2,2) = s2
		s1s2(i2,i1,1) = s1
		s1s2(i2,i1,2) = s2
        F(i1,i2,:) = FF
        F(i2,i1,:) = -F(i1,i2,:) + (/1.0e-6,0.0,0.0/)
        if (np == 6) then
			M(i1,i2,:) = MM1
			M(i2,i1,:) = MM2
		endif
    enddo
enddo
Fsum = 0
Msum = 0
do i1 = 1,ncells
    p1 => cell_list(i1)
    do j = 1,p1%nbrs
        i2 = p1%nbrlist(j)
        Fsum = Fsum + F(i1,i2,:)
        if (np == 6) then
	        Msum = Msum + M(i1,i2,:)
	    endif
    enddo
!	if (Fjiggle) then
!		do i = 1,3
!			R = par_rnor(kpar)
!			Fsum = Fsum + R*Fjigglefactor
!		enddo
!	endif
!	if (Mjiggle) then
!		do i = 1,3
!			R = par_rnor(kpar)
!			Msum = Msum + R*Mjigglefactor
!		enddo
!	endif
    k1 = (i1-1)*np
    xp(k1+1:k1+3) = Fsum/Fdrag
!!    mamp = sqrt(dot_product(Msum,Msum))
!!    vm = Msum/mamp     ! unit vector of moment axis
!!    dangle = mamp/Mdrag
!!    orient1 = x(k1+4:k1+6)
!!    call rotate(orient1,vm,dangle,xp(k1+4:k1+6))
	if (np == 6) then
		call cross_product(Msum/Mdrag,orient1,xp(k1+4:k1+6))
	endif
enddo
!write(*,*) 'xp:'
!write(*,'(3f7.3,4x,3f7.3)') xp(1:6*ncells)
ok = .true.
end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
subroutine solver(dt, nt,ok)
real(REAL_KIND) :: dt
integer :: nt
logical :: ok
integer :: kcell, k, j, nvars, flag, res
type(cell_type), pointer :: p
real(REAL_KIND), allocatable :: x(:)        ! the cell position and orientation variables
real(REAL_KIND), allocatable :: xp(:)       ! derivs of cell position and orientation variables
real(REAL_KIND) :: tstart, tend, relerr, abserr
real(REAL_KIND) :: amp
real(REAL_KIND), allocatable :: Jac(:,:)
integer :: np 

if (simulate_rotation) then
	np = 6
else
	np = 3
endif

!do kcell = 1,ncells
!    p =>cell_list(kcell)
!	if (p%a < p%b) then
!		write(*,*) 'Error: Solver (1): a < b: cell: ',kcell,p%a,p%b
!		stop
!	endif
!enddo
nvars = np*ncells
allocate(x(nvars),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: x'
	stop
endif
allocate(xp(nvars),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: xp'
	stop
endif
allocate(F(ncells,ncells,3),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: F'
	stop
endif
allocate(Jac(nvars,nvars),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: Jac'
	stop
endif
allocate(M(ncells,ncells,3),stat=res)
if (res /= 0) then
	write(*,*) 'Allocate failed for: M'
	stop
endif
!write(*,*) 'Memory required: ',8*(2*nvars + 2*3*ncells*ncells)
!allocate(cp(ncells,ncells,3))
!allocate(nbrs(ncells))
!allocate(nbrlist(ncells,MAX_NBRS))

do kcell = 1,ncells
    p =>cell_list(kcell)
    k = (kcell-1)*np
    x(k+1:k+3) = p%centre
    if (np == 6) then
		x(k+4:k+6) = p%orient
	endif
enddo

if (isolver == EULER_SOLVER) then
	do k = 1,2
		call JacobianFast(Jac, x, nvars, np, ok)
	!	write(nflog,'(9f8.3)') Jac
		write(*,*) 'got Jac'
	enddo
	stop
else

	tstart = 0
	xp = 0
	F = 0
	M = 0
	call fderiv(tstart,x,xp)

	abserr = sqrt ( epsilon ( abserr ) )
	relerr = sqrt ( epsilon ( relerr ) )

	flag = 1
	do j = 1,nt
		tstart = (j-1)*dt
		tend = tstart + dt
		call r8_rkf45 ( fderiv, nvars, x, xp, tstart, tend, relerr, abserr, flag )
		if (flag == 4) then
			call r8_rkf45 ( fderiv, nvars, x, xp, tstart, tend, relerr, abserr, flag )
		endif
		if (flag /= 2) then
			write(logmsg,*) 'Bad flag: ',flag
			call logger(logmsg)
			deallocate(x)
			deallocate(xp)
			deallocate(F)
			deallocate(M)
			ok = .false.
			return
		endif
		flag = 2
	enddo
endif
do kcell = 1,ncells
    p => cell_list(kcell)
    k = (kcell-1)*np
    p%centre = x(k+1:k+3)
    if (np == 6) then
		p%orient = x(k+4:k+6)
		amp = sqrt(dot_product(p%orient,p%orient))
		p%orient = p%orient/amp
	endif
enddo
deallocate(x)
deallocate(xp)
deallocate(F)
deallocate(M)
deallocate(Jac)
ok = .true.
end subroutine

end module

