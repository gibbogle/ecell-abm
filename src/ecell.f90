! ABM for ellipsoidal cells, off-lattice 
module ecell

use ellipsoid
use solve

implicit none

contains

!----------------------------------------------------------------------------------------
! Note: ellipsoid volume = 4/3.PI.a.b^2
! Say the "normal" cell volume is that of a sphere with radius Rn = 10 um
! Then normal volume Vn = 4/3.PI.Rn^3
! Further, assume that a dividing cell splits when its volume reaches 1.6Vn, and
! the progeny each have volume 0.8Vn.  This implies that cells in a random population
! have volumes in the range 0.8Vn - 1.6Vn.
! How to choose a and b?
! Case 1: (Simplest model)
! Assume that b does not vary as the cell grows, i.e. the cell just elongates.  This is
! an extreme case. For a normal cell, with volume = Vn, let a/b = alpha, the aspect ratio.
! 4/3.PI.alpha.b^3 = Vn = 4/3.PI.Rn^3
! therefore constant b = Rn.alpha^(-1/3) and for the normal cell a = a_n = Rn.alpha^(2/3)
! In this case the cell grows by increasing a until a = 1.6*a_n, then it is replaced 
! by two cells, each with a = 0.8*a_n (i.e. total length is unchanged)
! Case 2:
! Assume that the aspect ratio is always preserved.  In this case both a and b vary.
! The "normal" cell dimensions are a_n = Rn.alpha^(2/3) and b_n = Rn.alpha^(-1/3)
! and as the cell grows both a and b increase.  
! When volume = 1.6*Vn, a and b have increased by a factor beta over a_n and b_n, where
! beta^3 = 1.6, i.e. beta = 1.17
! For a progeny cell, the factor beta is: beta^3 = 0.8, i.e. beta = 0.928.
! In other words, the total length increases from 1.17*a_n to 1.86*a_n, a factor = 1.587)
! This increase in length is a potential problem, if it occurs instantaneously on
! cell division.  Large repulsion forces will immediately be created, risking instability.
! Two possible approaches:
! (a) the two progeny cells could be positioned to lie, as much as possble, side by side
! within the space occupied by the parent.  This is tricky.
! (b) over some time interval before division, the cell aspect ratio changes such that
! at the point of division the b value is equal to that of the progeny cell, 0.928*b_n
! instead of 1.17*b_n.  This implies that a must increase by a factor of (1.17/0.928)^2
! = 1.587 (presumably).  In other words, at the time of cell division the parent cell
! is twice as long as the progeny cell, as in Case 1.
!
! We can think of a growth model with two parameters.
! alpha = aspect ratio of "normal" cell at this location
! beta = level of constancy of the aspect ratio.
!    beta = 0 gives Case 1, in which b = b_n is constant
!    beta = 1 gives Case 2, in which a/b = alpha is constant
!    for intermediate values of beta: 
!
!    a = a_n(0.8(g+1))^(1 - (2/3)beta)
!    b = b_n(0.8(g+1))^(beta/3)
!
! Before the instant of cell division the cell shape must transition to a1 = 2.a0, such that the
! total length after division matches the parent cell length.  
! At g=1, a1 = 2.a0, where a0 = a_n.0.8^(1 - (2/3)beta)
! This transition should occur over some fixed time interval (e.g. 30 min).
!----------------------------------------------------------------------------------------
subroutine test
type(cell_type) :: ell1, ell2
ell1%a = 1.0
ell1%b = 0.50
ell1%centre = (/0.5, 0.0, 0.0/)
ell1%orient = (/1.0, 0.0, 0.0/)
ell1%F = 0
ell1%M = 0
ell1%growthrate = 0
ell1%nbrs = 0
allocate(ell1%nbrlist(20))
ell2%a = 2.0
ell2%b = 0.75
ell2%centre = (/0.5, 0.0, 1.0/)
ell2%orient = (/2.0, 1.0, -0.2/)
ell2%F = 0
ell2%M = 0
ell2%growthrate = 0
ell2%nbrs = 0
allocate(ell2%nbrlist(20))

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
real(REAL_KIND) :: t_div_median, t_div_shape, hours, deltat
integer :: nb0, nt_anim, igrow
character*(256) :: cellmlfile

open(nfin,file=inputfile)
read(nfin,*) nb0
read(nfin,*) igrow
read(nfin,*) t_div_median
read(nfin,*) t_div_shape
read(nfin,*) hours
read(nfin,*) deltat
read(nfin,*) isolver
read(nfin,*) seed(1)
read(nfin,*) seed(2)
read(nfin,*) ncpu_input
read(nfin,*) nt_anim
read(nfin,*) Fdrag
read(nfin,*) Mdrag
read(nfin,*) Falpha
read(nfin,*) Malpha
read(nfin,*) Fjigglefactor
read(nfin,*) Mjigglefactor
read(nfin,*) cellmlfile
close(nfin)

NX = nb0
NY = nb0
NZ = nb0
simulate_rotation = (Mdrag > 0)
simulate_growth = (igrow == 1)
DELTA_T = deltat/60 ! sec -> min
nsteps = hours*60./DELTA_T
write(logmsg,*) 'isolver: ',isolver
call logger(logmsg)
write(logmsg,*) 'hours,DELTA_T, nsteps: ',hours,DELTA_T,nsteps
call logger(logmsg)
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: kcell, j, error

ok = .true.
initialized = .false.
par_zig_init = .false.
istep = 0

inputfile = infile
outputfile = outfile
call logger("ReadCellParams")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")

if (ncpu == 0) then
	ncpu = ncpu_input
endif
Mnodes = ncpu
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call RngInitialisation

call PlaceCells
do kcell = 1,ncells
	call SetNeighbours(kcell)
    write(*,*) 'cell: ',kcell, cell_list(kcell)%nbrs
    write(*,'(20i4)') (cell_list(kcell)%nbrlist(j),j=1,cell_list(kcell)%nbrs)
enddo
cell_list0 = cell_list
write(nflog,*) 'ncells: ',ncells


FP1%a = 2
FP1%b = 2
FP1%c = 10
FP1%e = 1.5
FP1%g = FP1%e/2.5

FP2%a = 10   ! Fr = a when d/b = 1
FP2%b = 2
FP2%e = 1.5
FP2%g = 1   ! Fa = g when d = e

run_kmax = 0
run_dmax = 0

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine PlaceCells
integer :: ix, iy, iz, kcell, i, kpar=0
real(REAL_KIND) :: dx, dy, dz, x, y, z, a, b
real(REAL_KIND) :: centre(3), orient(3), orient0(3), Vn, aspect, beta
real(REAL_KIND) :: cycletime, age
type(cell_type), pointer :: p

if (allocated(cell_list)) then
	deallocate(cell_list)
	deallocate(cell_list0)
	deallocate(s1s2)
endif
allocate(cell_list(MAX_CELLS))
allocate(cell_list0(MAX_CELLS))
allocate(s1s2(MAX_CELLS,MAX_CELLS,2))
s1s2 = 0.5

a = 5
b = 3
Vn = (4./3.)*PI*a*b**2
aspect = a/b
beta = 1.0
dx = 2*b
dy = dx
dz = 2.0*a
orient0 = (/ 0., 0., 1. /)
kcell = 0
do ix = 1,nx
    do iy = 1,ny
        do iz = 1,nz
            x = ix*dx
            y = iy*dy
            z = iz*dz
            centre = (/ x, y, z /)
            do i = 1,3
                orient(i) = orient0(i) + 0.05*(par_uni(kpar) -0.5)
            enddo
            cycletime = CYCLETIME0*(1 + 0.2*(par_uni(kpar)-0.5))
            age = 0.8*cycletime*par_uni(kpar)
            kcell = kcell + 1
!            call CreateCell(kcell,centre,orient,a,g,alpha_n,beta)
			call CreateCell(kcell,centre,orient,Vn,aspect,beta,cycletime,age)
        enddo
    enddo
enddo
ncells = kcell
!do kcell = 1,ncells
!	p =>cell_list(kcell)
!	if (p%a < p%b) then
!		write(*,*) 'Error: PlaceCells: a < b: cell: ',kcell,p%a,p%b
!		stop
!	endif
!enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Should specify:
!    Vn, aspect, beta, cycletime, age
!    centre, orient
!
!    a = a_n(0.8(g+1))^(1 - (2/3)beta)
!    b = b_n(0.8(g+1))^(beta/3) 
!-----------------------------------------------------------------------------------------
subroutine CreateCell(kcell,centre,orient,Vn,aspect,beta,cycletime,age)
integer :: kcell
real(REAL_KIND) :: centre(3),orient(3),Vn,aspect,beta,cycletime,age
real(REAL_KIND) :: g, tnow
type(cell_type), pointer :: p

tnow = max(istep*DELTA_T,0.0)
p => cell_list(kcell)
p%centre = centre
p%orient = orient
p%cycletime = cycletime
p%birthtime = tnow - age
g = GrowthFunction(age,cycletime)
p%aspect_n = aspect
p%beta = beta
p%V_n = Vn
p%b_n = (3*Vn/(4*PI*aspect))**(1./3.)
p%a_n = aspect*p%b_n
p%a = p%a_n*(0.8*(g+1))**(1-2*beta/3)
p%b = p%b_n*(0.8*(g+1))**(beta/3)
!if (p%a < p%b) then
!	write(*,*) 'Error: CreateCell: a < b: cell: ',kcell,p%a,p%b
!	stop
!endif
p%Fprev = 0
p%Mprev = 0
allocate(p%nbrlist(MAX_NBRS))
call NormaliseOrient(p)
if (dbug) write(nflog,'(a,i6,2f8.2)') 'a, b: ',kcell,p%a,p%b
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetNeighbours(kcell)
integer :: kcell
integer :: n, i, jcell, nbrlist(1000)
real(REAL_KIND) :: c0(3), c1(3), r(3), d2, d2lim0, d2lim, d2list(1000)
integer, allocatable :: indx(:)
type(cell_type), pointer :: p
logical, parameter :: sort = .true.
real(REAL_KIND) :: nbfactor = 3.0

p => cell_list(kcell)
d2lim0 = (nbfactor*p%a)**2
c0 = p%centre
n = 0
do jcell = 1,ncells
    if (jcell == kcell) cycle
    c1 = cell_list(jcell)%centre
    r = c1 - c0
    d2lim = max(d2lim0,(nbfactor*cell_list(jcell)%a)**2)
    d2 = dot_product(r,r)
    !if (istep == 400 .and. (kcell==8 .or. kcell==9)) then
    !    write(*,'(2i4,3f8.3,2x,L3)') kcell,jcell,d2,d2lim0,d2lim,(d2 < d2lim)
    !endif
    if (d2 < d2lim) then
		if (.not.sort .and. n == MAX_NBRS) then
			write(*,*) 'Error: SetNeighbours: n > MAX_NBRS'
			write(nflog,*) 'Error: SetNeighbours: n > MAX_NBRS'
			stop
		endif
        n = n+1
        if (sort) then
	        nbrlist(n) = jcell
		    d2list(n) = d2
		else
	        p%nbrlist(n) = jcell
	    endif
    endif
enddo
if (sort) then
	! Now we need to sort by d2
	allocate(indx(n))
	do i = 1,n
		indx(i) = i
	enddo
	call qsort(d2list,n,indx)     ! sorts in increasing order 
	n = min(n,MAX_NBRS)
endif
!write(nflog,*) 'SetNeighbours: ',kcell, n
p%nbrs = n
if (sort) then
	do i = 1,n
		jcell = indx(i)
		p%nbrlist(i) = nbrlist(jcell)
	enddo
	deallocate(indx)
endif
!if (kcell < 20) then
!	write(nflog,'(20i4)') p%nbrlist(1:n)
!endif
!    write(nflog,'(a,i6,a,i2)') 'cell: ',kcell,' nbrs: ',n
!    write(nflog,'(15i5)') p%nbrlist(1:n)
end subroutine

logical function nancomponent(v)
real(REAL_KIND) :: v(3)

nancomponent =  (isnan(v(1)) .or. isnan(v(2)) .or. isnan(v(3)))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SumContacts(ok)
integer :: kcell, k, jcell
real(REAL_KIND) :: delta, s1, s2, famp, mamp, fsum, msum
real(REAL_KIND) :: F(3), M1(3), M2(3)
logical :: incontact, ok
type(cell_type), pointer :: p, p1

do kcell = 1,ncells
    p => cell_list(kcell)
    p%F = 0
    p%M = 0
enddo
fsum = 0
msum = 0
do kcell = 1,ncells
    p => cell_list(kcell)
    p%F = 0
    p%M = 0
    do k = 1,p%nbrs
        jcell = p%nbrlist(k)
        if (jcell == kcell) then
            write(*,*) 'jcell = kcell: ',jcell
            stop
        endif
	    if (dbug) write(nflog,*) 'Cell: ',kcell,' nbr: ',p%nbrs,k,jcell
        p1 => cell_list(jcell)
        if (nancomponent(p%centre)) then
			write(*,*) 'Bad p%centre: ',kcell,p%centre
			ok = .false.
			return
		endif
        if (nancomponent(p%orient)) then
			write(*,*) 'Bad p%orient: ',kcell,p%orient
			ok = .false.
			return
		endif
        if (nancomponent(p1%centre)) then
			write(*,*) 'Bad p1%centre: ',kcell,p1%centre
			ok = .false.
			return
		endif
        if (nancomponent(p1%orient)) then
			write(*,*) 'Bad p1%orient: ',kcell,p1%orient
			ok = .false.
			return
		endif
		s1 = 0.5	! initial guesses
		s2 = 0.5
		s1 = s1s2(kcell,jcell,1)
		s2 = s1s2(kcell,jcell,2)
        call CellInteraction(p%a,p%b,p%centre,p%orient,p1%a,p1%b,p1%centre,p1%orient,s1,s2,delta,F,M1,M2,ok)
        if (.not.ok) return
		s1s2(kcell,jcell,1) = s1
		s1s2(kcell,jcell,2) = s2
		s1s2(jcell,kcell,1) = s1
		s1s2(jcell,kcell,2) = s2
	    p%F = p%F + F
	    p%M = p%M + M1
	    if (nancomponent(F)) then
			write(*,'(a,2i6,3f8.4)') 'Bad F: ',kcell,k,F
			stop
		endif
	    if (nancomponent(M1)) then
			write(*,'(a,2i6,3f8.4)') 'Bad M1: ',kcell,k,M1
			stop
		endif
!		if (kcell == 14) then
!			write(nflog,'(2i4,f8.4,2x,3f8.4,2x,3f8.4)') kcell,jcell,delta,F,M1
!		endif
    enddo
    famp = sqrt(dot_product(p%F,p%F))
    if (dbug) write(nflog,*) 'famp: ',famp
    fsum = fsum + famp
    mamp = sqrt(dot_product(p%M,p%M))
    msum = msum + mamp
enddo
do kcell = 1,ncells
    p => cell_list(kcell)
    p%F = (1-Falpha)*p%F + Falpha*p%Fprev
    p%M = (1-Malpha)*p%M + Malpha*p%Mprev
    p%Fprev = p%F
    p%Mprev = p%M
enddo
!write(logmsg,'(a,2e12.4)') 'SumContacts: average F, M: ',fsum/ncells,msum/ncells
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Mover(dt)
real(REAL_KIND) :: dt
integer :: k, kcell, i, kpar=0
integer, allocatable :: permc(:)
real(REAL_KIND) :: mamp, dangle, vm(3), v(3), R, del(3), rot(3), amp
type(cell_type), pointer :: p
logical :: Fjiggle, Mjiggle

!write(nflog,*) 'mover: ',ncells
Fjiggle = (Fjigglefactor > 0)
Mjiggle = (Mjigglefactor > 0)
allocate(permc(ncells))
do k = 1,ncells
    permc(k) = k
enddo
call permute(permc,ncells,kpar)

do k = 1,ncells
    kcell = permc(k)
    p => cell_list(kcell)
!    write(nflog,'(2i6,6e12.4)') k,kcell,p%F,p%M
	if (Fjiggle) then
		do i = 1,3
			R = par_rnor(kpar)
			p%F(i) = p%F(i) + R*Fjigglefactor
		enddo
	endif
	if (simulate_rotation .and. Mjiggle) then
		do i = 1,3
			R = par_rnor(kpar)
			p%M(i) = p%M(i) + R*Mjigglefactor
		enddo
	endif
	del = (dt/Fdrag)*p%F
    p%centre = p%centre + del
    if (dbug) write(nflog,'(i6,i3,6f8.4)') istep,k,del,p%centre
	if (simulate_rotation) then
		mamp = sqrt(dot_product(p%M,p%M))
		if (mamp == 0) cycle
		vm = p%M/mamp     ! unit vector of moment axis
		dangle = (dt/Mdrag)*mamp
		! Two possible ways to apply the moment
		! 1.
!		call rotate(p%orient,vm,dangle,v)
!		p%orient = v
		! 2.
		call cross_product(p%M/Mdrag,p%orient,rot)
		p%orient = p%orient + dt*rot
		amp = sqrt(dot_product(p%orient,p%orient))
		p%orient = p%orient/amp
		
!		if (kcell == 14) then
!			write(nflog,'(8f9.4)') mamp,dangle,vm,p%orient
!		endif
!		if (nancomponent(v)) then
!			write(*,'(a,i6,3f8.4)') 'Bad rotation: ',k,p%orient
!			write(*,'(a,3f8.3)') 'p%M: ',p%M
!			write(*,'(a,5f8.3)') 'mamp,vm,dangle: ',mamp,vm,dangle
!			write(*,'(a,3f8.4)') 'to: ',v
!			stop
!		endif
	endif
enddo

deallocate(permc)
end subroutine

!-----------------------------------------------------------------------------------------
! This is the simplest form of growth function - linear with time.
!-----------------------------------------------------------------------------------------
function GrowthFunction(t,tcycle) result(g)
real(REAL_KIND) :: t, tcycle, g

g = t/tcycle
g = min(g,1.0)
end function

!-----------------------------------------------------------------------------------------
! We can assume that for a given tissue region/stage there is an idea of a "normal" cell.
! This is the cell in its resting state.  
! The "normal" volume Vn = (4/3).Pi.r_n^3 where r_n is the radius of the equivalent sphere.
! If we assume that a cell divides when its volume = 1.6Vn, we can define a growth parameter
! g: 0 -> 1 (g is stage of growth) and V/Vn = 0.8(g+1).
! Note that the "normal" volume Vn is reached when g = 0.25.
!
! In the general case, g = g(t) where t = time since birth (not necessarily linear)
!
! For an ellipsoid V = (4/3).Pi.a.b^2, Vn = (4/3).Pi.a_n.b_n^2
! We can define the"normal" aspect ratio aspect_n = a_n/b_n
! There are two corner cases for cell growth:
! Case 1: b is constant, b = b_n (beta = 0)
! Case 2: aspect ratio = a/b is constant: a/b = aspect_n (beta = 1)
!
! Case 1.  (a.b^2)/(a_n.b_n^2) = V/Vn = 0.8(g+1)
!          a = a_n.0.8(g+1)
!          b = b_n
! Case 2.  (a.b^2)/(a_n.b_n^2) = 0.8(g+1)
!          a = aspect_n.b
!          a^3 = (aspect_n^2.a_n.b_n^2).0.8(g+1) = a_n^3.0.8(g+1)
!          a = a_n.[0.8(g+1)]^(1/3)
!          b = b_n.[0.8(g+1)]^(1/3)
!
! In the general case, 0 < beta < 1:
!          a = a_n.[0.8(g+1)]^(1-(2/3)beta)
!          b = b_n.[0.8(g+1)]^(beta/3)
! and
!          a/b = (a_n/b_n).[0.8(g+1)]^(1-beta)
!
! Around cell division the shape must be progressively modified to ensure that
! at the instant of cell division the cell's a and b are such that:
! with g=1, a1 = the twice the length of a progeny cell at birth, 
! i.e. (with g=0) = a0 = (a_n/b_n).[0.8]^(1-(2/3)beta), we want a1 = 2.a0.
! Starting at g = g_c, apply this adjustment to a:
! a = ((1-g)/(1-g_c)).a + ((g-g_c)/(1-g_c)).2.a0
! writing alpha = (g-g_c)/(1-g_c)
! a = (1-alpha).a_n.[0.8(g+1)]^(1-(2/3)beta) + alpha.2.a_n.(0.8)^(1-(2/3)beta)
!-----------------------------------------------------------------------------------------
subroutine Grower(ok)
logical :: ok
integer :: kcell, ncells0, nbmax
real(REAL_KIND) :: a, b, t, g, V, e, alfa, tnow, a0
type(cell_type), pointer :: p
logical :: divided

divided = .false.
tnow = istep*DELTA_T
ncells0 = ncells
do kcell = 1,ncells0
	p => cell_list(kcell)
	t = tnow - p%birthtime
	g = GrowthFunction(t,p%cycletime)
	V = p%V_n*0.8*(g+1)
	e = 1 - (2./3.)*p%beta
	a = p%a_n*(0.8*(g+1))**e
	if (g >= 1) then
		call Divider(kcell,ok)
        divided = .true.
		if (.not.ok) return
		cycle
	endif
	if (g > g_threshold) then
		alfa = (g-g_threshold)/(1-g_threshold)
		a0 = p%a_n*0.8**e
		a = (1-alfa)*a + 2*alfa*a0
	endif
	b = sqrt(3.*V/(4.*PI*a))
	p%a = a
	p%b = b
!	if (p%a < p%b) then
!		write(*,*) 'Error: Grower: a < b: cell: ',kcell,p%a,p%b
!		stop
!	endif
enddo
if (divided) then
    nbmax = 0
    do kcell = 1,ncells
        call SetNeighbours(kcell)
        if (cell_list(kcell)%nbrs > nbmax) then
            nbmax = cell_list(kcell)%nbrs
        endif
    enddo
    !write(*,*) 
    !write(*,*) 'nbmax: ',nbmax
    !write(*,*)
endif

ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Divider(kcell0,ok)
logical :: ok
integer :: kcell0
integer :: kcell1, kpar=0
real(REAL_KIND) :: tnow, a, c(3)
type(cell_type), pointer :: p0, p1

!write(nflog,*) 'Divider: ',kcell0
tnow = istep*DELTA_T
ncells = ncells + 1
if (ncells == MAX_CELLS) then
	write(logmsg,*) 'Maximum number of cells reached: ',MAX_CELLS
	call logger(logmsg)
	ok = .false.
	return
endif
kcell1 = ncells
p0 => cell_list(kcell0)
c = p0%centre
p0%a = p0%a/2
if (p0%a < p0%b) then
	write(*,*) 'Error: Divider: a < b: cell: ',kcell1,p0%a,p0%b
	write(*,*) 'aspect_n, a_n, b_n, beta: ',p0%aspect_n, p0%a_n, p0%b_n, p0%beta
	write(logmsg,*) 'Error: Divider: a < b: cell: ',kcell1,p0%a,p0%b
	call logger(logmsg)
	ok = .false.
	return
endif

p0%birthtime = tnow
p0%cycletime = CYCLETIME0*(1 + 0.2*(par_uni(kpar)-0.5))
p0%Fprev = 0
p0%Mprev = 0
p0%centre = c - p0%a*p0%orient
cell_list(kcell1) = cell_list(kcell0)
p1 => cell_list(kcell1)
p1%centre = c + p0%a*p0%orient
p1%cycletime = CYCLETIME0*(1 + 0.2*(par_uni(kpar)-0.5))
if (allocated(p1%nbrlist)) then
	deallocate(p1%nbrlist)
endif
allocate(p1%nbrlist(MAX_NBRS))
!call SetNeighbours(kcell0)
!call SetNeighbours(kcell1)
ok = .true.
end subroutine

!--------------------------------------------------------------------------------
! The unit vector that is the axis of rotation to take the x axis (1 0 0) to v
! is (0 -vz/s vy/s) where s = sin(theta) = sqrt(vy^2 + vz^2), and c = cos(theta) = vx
!--------------------------------------------------------------------------------
subroutine get_scene(nEC_list,EC_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nEC_list 
real(c_double) :: EC_list(*)
type(cell_type), pointer :: p
real(REAL_KIND) :: a, b, vx, vy, vz, s, c, s2, cmax
integer :: kcell, j

cmax = 0
nEC_list = ncells
do kcell = 1,nEC_list
    p => cell_list(kcell)
    j = (kcell-1)*12
    a = p%a
    b = p%b
    vx = p%orient(1)
    vy = p%orient(2)
    vz = p%orient(3)
    s2 = vy*vy+vz*vz
    s = sqrt(s2);	! sin(theta)
    c = vx;			! cos(theta)
    EC_list(j+1:j+3) = p%centre
    EC_list(j+4:j+12) = (/ a*c, a*vy, a*vz, &
                          -b*vy, b*(c+(1-c)*vz*vz/s2), -b*(1-c)*vy*vz/s2, &
                          -b*vz, -b*(1-c)*vy*vz/s2, b*(c+(1-c)*vy*vy/s2) /)
    cmax = max(cmax,c)
enddo
!write(logmsg,*) 'get_scene: ',istep,cmax
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: it, nt
real(REAL_KIND) :: dt
integer :: kcell
type(cell_type), pointer :: p
logical :: ok

istep = istep + 1
if (mod(istep,1) == 0) then
	write(logmsg,*) 'step: ',istep,ncells
	call logger(logmsg)
endif
ok = .true.
if (simulate_growth) then
	call Grower(ok)
endif
if (.not.ok) then
	res = 1
	return
endif
if (isolver == EXPLICIT_SOLVER) then
	nt = 10
	dt = DELTA_T/nt
	do it = 1,nt
		call SumContacts(ok)
		if (.not.ok) then
			res = 1
			return
		endif
		call Mover(dt)
	enddo
	ok = .true.
else
    if (isolver == SLOW_EULER_SOLVER .or. isolver == FAST_EULER_SOLVER) then
	    nt = 1
        dt = DELTA_T/nt
	    call solver(dt,nt,ok)
    elseif (isolver == RKF45_SOLVER) then
	    nt = 10
	    dt = DELTA_T/nt
	    call solver(dt,nt,ok)
    endif
    call FixCentre
endif
if (ok) then
    res = 0
else
    res = 1
endif
end subroutine

subroutine FixCentre
real(REAL_KIND) :: origin(3)
integer :: kcell

origin = 0
do kcell = 1,ncells
    origin = origin + cell_list(kcell)%centre
enddo
origin = origin/ncells
do kcell = 1,ncells
    cell_list(kcell)%centre = cell_list(kcell)%centre - origin
enddo
end subroutine
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(nsteps_dim, deltat) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: nsteps_dim
real(c_double) :: deltat

nsteps_dim = nsteps
deltat = DELTA_T

end subroutine


!void get_dimensions(int *,int *,int *, int *, double *, int *, bool *, double *);
!void get_scene(int *, int *);
!void get_summary(int *, int *, int *);


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C)
!subroutine Execute() BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
integer :: i, res
logical :: ok, isopen

infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

awp_0%is_open = .false.
awp_1%is_open = .false.
par_zig_init = .false.
logfile = 'ecell.log'
inquire(unit=nflog,OPENED=isopen)
if (.not.isopen) then
    open(nflog,file=logfile,status='replace')
endif
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile 
call logger(logmsg)
if (use_TCP) then
	write(nflog,*) 'call connecter'
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
call Setup(ncpu,infile,outfile,ok)
if (ok) then
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
	res = 0
else
	call logger('=== Setup failed ===')
	res = 1
	stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine DisableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from spheroid_main()	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
else
    write(nflog,*) 'connection: awp open: ',port, error
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call Connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CheckAngles
real(REAL_KIND) :: cosa, cosamin
integer :: k, kmin
cosamin = 1.0e10
do k = 1,ncells
	cosa = dot_product(cell_list(k)%orient,cell_list0(k)%orient)
	if (cosa < cosamin) then
		cosamin = cosa
		kmin = k
	endif
enddo
write(nflog,*)
write(nflog,*) 'Minimum cosine of rotation angle: ',kmin,cosamin
write(nflog,*) 'angle: ',acos(cosamin)*180/PI
write(nflog,*)
write(nflog,'(a,3f8.4)') 'Initial orient: ',cell_list0(kmin)%orient
write(nflog,'(a,3f8.4)') 'Final orient:   ',cell_list(kmin)%orient
write(nflog,*)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
call CheckAngles
call Wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == -1) then
	call logger(' Execution stopped')
else
	call logger('  === Execution failed ===')
endif
!write(logmsg,'(a,f10.2)') 'Execution time (min): ',(wtime() - execute_t1)/60
call logger(logmsg)

!close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
		endif
	endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
integer :: ierr, ichemo
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
!call logger('deallocated all arrays')

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
end subroutine


end module
