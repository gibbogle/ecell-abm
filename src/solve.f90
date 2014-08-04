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
! x(:) holds the position (centre) and orientation (orient) of each cell.
! The neighbour list nbrlist(:,:) holds each cell's neighbours, assumed not to change
! over the duration of this solving step.
! The force between two cells is stored in F(:,:,:), and the contact points on the two cells
! is in cp(:,:,:).  If the force between cells i1 and i2 has already been computed: F(i1,i2,:),
! then F(i2,i1,:) is set = -F(i1,i2,:) and is not recomputed later when cell i2 is treated.
!------------------------------------------------------------------------------------------
subroutine fderiv(t,x,xp)
real(REAL_KIND) :: t, x(:), xp(:)
integer :: i1, i2, j, k1, k2, kpar=0
real(REAL_KIND) :: a1, b1, centre1(3), orient1(3), a2, b2, centre2(3), orient2(3) 
real(REAL_KIND) :: FF(3), MM1(3), MM2(3), Fsum(3), Msum(3), R, mamp, vm(3), dangle
type(cell_type), pointer :: p1, p2

F = 0
do i1 = 1,ncells
    p1 => cell_list(i1)
    a1 = p1%a
    b1 = p1%b
    k1 = (i1-1)*6
    centre1 = x(k1+1:k1+3)
    orient1 = x(k1+4:k1+6)
    do j = 1,p1%nbrs
        i2 = p1%nbrlist(j)
        if (F(i1,i2,1) /= 0) cycle
        p2 => cell_list(i2)
        a2 = p2%a
        b2 = p2%b
        k2 = (i2-1)*6
        centre2 = x(k2+1:k2+3)
        orient2 = x(k2+4:k2+6)
        call CellInteraction(a1,b1,centre1,orient1,p2%a,p2%b,centre2,orient2,FF,MM1,MM2)
        F(i1,i2,:) = FF
        F(i2,i1,:) = -F(i1,i2,:) + (/1.0e-10,0.0,0.0/)
        M(i1,i2,:) = MM1
        M(i2,i1,:) = MM2
    enddo
enddo
Fsum = 0
do i1 = 1,ncells
    p1 => cell_list(i1)
    do j = 1,p1%nbrs
        i2 = p1%nbrlist(j)
        Fsum = Fsum + F(i1,i2,:)
        Msum = Msum + M(i1,i2,:)
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
    k1 = (i1-1)*6
    xp(k1+1:k1+3) = Fsum/Fdrag
    mamp = sqrt(dot_product(Msum,Msum))
    vm = Msum/mamp     ! unit vector of moment axis
    dangle = mamp/Mdrag
    call rotate(x(k1+4:k1+6),vm,dangle,xp(k1+4:k1+6))
enddo
end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
subroutine solver(dt, nt)
real(REAL_KIND) :: dt
integer :: nt
integer :: kcell, k, nvars, flag
type(cell_type), pointer :: p
real(REAL_KIND), allocatable :: x(:)        ! the cell position and orientation variables
real(REAL_KIND), allocatable :: xp(:)       ! derivs of cell position and orientation variables
real(REAL_KIND) :: tstart, tend, relerr, abserr

nvars = 6*ncells
allocate(x(nvars))
allocate(xp(nvars))
allocate(F(ncells,ncells,3))
allocate(M(ncells,ncells,3))
!allocate(cp(ncells,ncells,3))
!allocate(nbrs(ncells))
!allocate(nbrlist(ncells,MAX_NBRS))

do kcell = 1,ncells
    p => cell_list(kcell)
    x((kcell-1)*6+1:(kcell-1)*6+3) = p%centre
    x((kcell-1)*6+4:(kcell-1)*6+6) = p%orient
!    nbrs(kcell) = p%nbrs
!    nbrlist(kcell,:) = p%nbrlist
enddo

tstart = 0
xp = 0
call fderiv(tstart,x,xp)

abserr = sqrt ( epsilon ( abserr ) )
relerr = sqrt ( epsilon ( relerr ) )

flag = 1
do k = 1,nt
	tstart = (k-1)*dt
	tend = tstart + dt
    call r8_rkf45 ( fderiv, nvars, x, xp, tstart, tend, relerr, abserr, flag )
	if (flag /= 2) then
		write(logmsg,*) 'Bad flag: ',flag
		call logger(logmsg)
	endif
	flag = 2
enddo

end subroutine

end module
