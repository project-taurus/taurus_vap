!==============================================================================!
! MODULE Wavefunctions                                                         !
!                                                                              !
! This module contains the variables and routines related to the wave funct-   !
! ions (including their densities).                                            !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_wavefunctions                                               !
! - subroutine generate_wavefunction                                           !
! - subroutine set_random_generation                                           !
! - subroutine generate_wavefunction_BCS                                       !
! - subroutine generate_wavefunction_slater                                    !
! - subroutine generate_unitary_matrix                                         !
! - subroutine block_quasiparticle                                             !
! - subroutine print_wavefunction                                              !
! - subroutine construct_canonical_basis                                       !
! - subroutine read_wavefunction                                               !
! - subroutine write_wavefunction                                              !
! - subroutine calculate_densities                                             !
! - subroutine calculate_densities_real                                        !
! - subroutine calculate_spatial_density                                       !
!==============================================================================!
MODULE Wavefunctions 

use Basis

implicit none
public

!!! Label to identify the wave function
integer(i64) :: bogo_label

!!! Parameters that determines the seed wave function
integer :: seed_type,   & ! type of seed
           seed_rand,   & ! seed to initialize the random generation
           seed_text,   & ! format of seed 
           seed_symm,   & ! option to check the symmetries
           seed_allemp, & ! option to include the empty states (overlap)
           blocking_dim   ! number of quasiparticle to block
integer, dimension(:), allocatable :: blocking_id ! indices to block 
real(r64) :: seed_occeps ! cutoff for occupied single-particle states
                
!!! Bogoliubov matrices
integer :: bogo_nparity ! Number parity of the state
real(r64), dimension(:,:), allocatable :: bogo_U0, & ! U present iteration
                                          bogo_V0, & ! V    "        "
                                          bogo_U1, & ! U evolved
                                          bogo_V1    ! V evolved
complex(r64), dimension(:,:), allocatable :: bogo_zU0,  & ! U present iteration
                                             bogo_zV0,  & ! V    "       "
                                             bogo_zU0c, & ! U canonical basis
                                             bogo_zV0c, & ! V    "        "        
                                             bogo_zD0     ! transf. to can. bas.    
                
!!! Density matrices
real(r64), dimension(:,:), allocatable :: dens_rhoRR, &      ! <R|a+a|R> (real)
                                          dens_kappaRR       ! <R|aa|R>  (real)
complex(r64), dimension(:,:), allocatable :: dens_rhoLR,   & ! <L|a+a|R> 
                                             dens_kappaLR, & ! <L|aa|R>
                                             dens_kappaRL    ! <L|a+a+|R>^*
integer :: dens_spatial, & ! option to compute the spatial one-body density
           dens_nr(3)   ! number of points in the discretization
real(r64) :: dens_dr(3) ! spacing between the points 

!!! Information on the good quantum number of the wave function
logical :: is_good_Z=.false.,      & ! Slater determinant for protons
           is_good_N=.false.,      & ! Slater determinant for neutrons
           is_separate_NZ=.false., & ! no p-n pairing
           is_good_P=.false.,      & ! good parity
           is_good_K=.false.         ! good angular-momentum third component
        
CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_wavefunctions                                                 !
!                                                                              !
! Initializes the arrays related to the wave functions                         !
!------------------------------------------------------------------------------!
subroutine set_wavefunctions

integer :: ialloc=0

!!! Wavefunctions
allocate( bogo_U0(HOsp_dim,HOsp_dim),   bogo_V0(HOsp_dim,HOsp_dim),   &
          bogo_U1(HOsp_dim,HOsp_dim),   bogo_V1(HOsp_dim,HOsp_dim),   &
          bogo_zU0(HOsp_dim,HOsp_dim),  bogo_zV0(HOsp_dim,HOsp_dim),  &
          bogo_zU0c(HOsp_dim,HOsp_dim), bogo_zV0c(HOsp_dim,HOsp_dim), & 
          bogo_zD0(HOsp_dim,HOsp_dim), stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of wave functions'

bogo_U0 = zero
bogo_V0 = zero
bogo_U1 = zero
bogo_V1 = zero
bogo_zU0  = zzero
bogo_zV0  = zzero
bogo_zD0  = zzero
bogo_zU0c = zzero
bogo_zV0c = zzero

!!! Densities
allocate( dens_rhoRR(HOsp_dim,HOsp_dim), dens_kappaRR(HOsp_dim,HOsp_dim), &
          dens_rhoLR(HOsp_dim,HOsp_dim), dens_kappaLR(HOsp_dim,HOsp_dim), &
          dens_kappaRL(HOsp_dim,HOsp_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of densities'

dens_rhoRR   = zero
dens_kappaRR = zero
dens_rhoLR   = zzero
dens_kappaLR = zzero
dens_kappaRL = zzero

end subroutine set_wavefunctions

!------------------------------------------------------------------------------!
! subroutine generate_wavefunction                                             !
!                                                                              !
! Generates the seed wave function or reads it from file. Different type of    !
! symmetries are possible for the seed wave function.                          !
!                                                                              !
! seed_type = 0 fully arbitrary wavefunction                                   !
!           = 1 read from file "initial_wf"                                    !
!           = 2 BCS state spherical (includes good parity + no p-n mixing)     !
!           = 3 BCS state axial (includes good parity + no p-n mixing)         !
!           = 4 good parity                                                    !
!           = 5 no p-n mixing                                                  !
!           = 6 good parity + no p-n mixing                                    !
!           = 7 Slater determinant                                             !
!           = 8 Slater + good parity                                           !
!           = 9 Slater HO (includes good parity + good Jz)                     !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!------------------------------------------------------------------------------!
subroutine generate_wavefunction(ndim)

integer, intent(in) :: ndim
integer :: i, nocc0, nemp0
real(r64) :: ovac0                                                              
real(r64), dimension(ndim,ndim) :: A1, A2, A3
!cmpi integer :: ierr=0

!!! Sets the pseudo random number generation
call set_random_generation

!cmpi if ( paral_myrank == 0 ) then        
select case (seed_type)

  !!! Reads the wave function from a file
  case (1)
    call read_wavefunction

  !!! Generates U,V for the HFB states. Starts by generating U,V in the canon-
  !!! ical basis (BCS) and then performs the Bloch-Messiash transformations C,D 
  case (0,2:6)
    call generate_wavefunction_BCS(valence_Z,valence_N,bogo_U0,bogo_V0,ndim)
    if ((seed_type /= 2) .and. (seed_type /= 3) ) then 
      ! Transformation D (real for now)
      call generate_unitary_matrix(seed_type,A3,ndim)
      call dgemm('n','n',ndim,ndim,ndim,one,A3,ndim,bogo_U0,ndim,zero,A1,ndim)
      call dgemm('n','n',ndim,ndim,ndim,one,A3,ndim,bogo_V0,ndim,zero,A2,ndim)
      ! Transformation C (real for now)
      call generate_unitary_matrix(seed_type,A3,ndim)
      call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,A3,ndim,zero,bogo_U0,ndim)
      call dgemm('n','n',ndim,ndim,ndim,one,A2,ndim,A3,ndim,zero,bogo_V0,ndim)
    endif
  
  !!! Generates U,V for Slater determinants.
  case (7:9)
    call generate_wavefunction_slater(int(valence_Z),int(valence_N),bogo_U0, &
                                      bogo_V0,ndim)
    if ( seed_type /= 9 ) then
      ! Bloch-Messiah: transformation D (real for now)
      call generate_unitary_matrix(seed_type,A3,ndim)
      call dgemm('n','n',ndim,ndim,ndim,one,A3,ndim,bogo_U0,ndim,zero,A1,ndim)
      call dgemm('n','n',ndim,ndim,ndim,one,A3,ndim,bogo_V0,ndim,zero,A2,ndim)
      bogo_U0 = A1
      bogo_V0 = A2
    endif

  case default
    stop 'Wrong argument for seed_type.'

end select

!!! Performs the blocking of quasiparticles
do i = 1, blocking_dim
  call block_quasiparticle(blocking_id(i),bogo_U0,bogo_V0,ndim)
enddo
!cmpi endif

!!! Broadcasts the wave function from rank = 0 to all MPI processes
!cmpi call mpi_bcast(bogo_U0,ndim**2,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(bogo_V0,ndim**2,mpi_double_precision,0,mpi_comm_world,ierr)

!!! Builds the canonical basis and prints the important informations
call construct_canonical_basis(bogo_U0,bogo_V0,bogo_zU0c,bogo_zV0c,bogo_zD0, &
                               ovac0,nocc0,nemp0,ndim)
bogo_nparity = (-1)**nocc0

!cmpi if ( paral_myrank == 0 ) then        
call print_wavefunction(nocc0)
!cmpi endif

end subroutine generate_wavefunction

!------------------------------------------------------------------------------!
! subroutine generate_wavefunction_BCS                                         !
!                                                                              ! 
! This subroutine sets the seed for the pseudo random generation using the     ! 
! input parameter seed_rand.                                                   ! 
!                                                                              ! 
! seed_rand = 0 uses the state of the processor (changes with each run)        ! 
!           > 0 uses the values given in input + a small algorithm             ! 
!------------------------------------------------------------------------------!
subroutine set_random_generation

integer :: i, mdim
integer, dimension(:), allocatable :: mentry 

!!! Initialization using the state of the processor
call random_seed()

!!! Sets a seed with a little bit of variability using seed_rand
if ( seed_rand > 0 ) then
  call random_seed(size=mdim)
  allocate (mentry(mdim))

  do i = 1, mdim
    mentry(i) = seed_rand + (i-1) * (2020 + mod(i,2)*10 + mod(i,3)*27)
  enddo

  call random_seed(put=mentry(1:mdim))
endif

end subroutine set_random_generation

!------------------------------------------------------------------------------!
! subroutine generate_wavefunction_BCS                                         !
!                                                                              ! 
! This subroutine generates random Bogoliubov matrices U and V corresponding   ! 
! to a BCS product state with the correct expectation values for the number of ! 
! neutrons and protons.                                                        ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      ! 
!        nprot = number of protons                                             ! 
!        npneut= number of neutrons                                            ! 
!                                                                              ! 
! Output: U,V = Bogoliubov matrices                                            ! 
!------------------------------------------------------------------------------!
subroutine generate_wavefunction_BCS(nprot,nneut,U,V,ndim)

integer, intent(in) :: ndim
real(r64), intent(in) :: nprot, nneut
real(r64), dimension(ndim,ndim), intent(out) :: U, V
integer :: i, j, k, m, mmax, ipos, ipos1, ipos2, ipos3, ipos4
real(r64) :: ranval, fac_v2c, num_v2c, factor
real(r64), dimension(ndim/2) :: v2c

U = zero
V = zero
v2c = zero

!!! Generates random occupancies v^2 for protons. The code first increases
!!! the occupancies until \sum_p v^2 > nprot but the check is done only after
!!! increasing v^2 for all single-particle states (I don't want to favor the
!!! first states in the basis). Then the v^2 are rescaled to obtain the 
!!! correct values \sum_p v^2 = nprot
num_v2c = -1.d0
fac_v2c =  1.d0
  
do while ( num_v2c < 0.98*nprot )
  num_v2c = 0.d0
  do i = 1, ndim/4
    call random_number(ranval)

    !!! Equal filling of an orbit for spherical seed
    if ( (seed_type == 2) .and. (i > 1) ) then
      if ( HOsp_sh(2*i-1) == HOsp_sh(2*i-3) ) then
        ranval = 0.d0
        v2c(i) = v2c(i-1)
      endif
    endif

    v2c(i) = min(0.99,v2c(i)+ranval)
    num_v2c = num_v2c + 2.d0 * v2c(i)
  enddo
enddo

if ( (num_v2c > nprot) .and. (nprot > 0.d0) ) then
  fac_v2c = num_v2c / nprot
  v2c(1:ndim/4) = v2c(1:ndim/4) / fac_v2c
endif

!!! Generates random occupancies v^2 for neutrons. The same strategy is applied
!!! but now to build \sum_n v^2 = nneut.
num_v2c = -1.d0
fac_v2c =  1.d0

do while ( num_v2c < 0.98*nneut )
  num_v2c = 0.d0
  do i = 1+ndim/4, ndim/2
    call random_number(ranval)
    
    if ( (seed_type == 2) .and. (i > 1) ) then
      if ( HOsp_sh(2*i-1) == HOsp_sh(2*i-3) ) then
        ranval = 0.d0
        v2c(i) = v2c(i-1)
      endif
    endif

    v2c(i) = min(0.99,v2c(i)+ranval)
    num_v2c = num_v2c + 2.d0 * v2c(i)
  enddo
enddo

if ( (num_v2c > nneut) .and. (nneut > 0.d0) ) then
  fac_v2c = num_v2c / nneut  
  v2c(1+ndim/4:ndim/2) = v2c(1+ndim/4:ndim/2) / fac_v2c
endif

!!! Builds U and V in the canonical basis where u^2 + v^2 = 1. The pairing is 
!!! built between time-reversal partners: |jm> and |j-m>
ipos = 0
k = 0
do i = 1, HOsh_dim
  j = HOsh_2j(i)
  mmax = (j+1)/2
  do m = 1, mmax
    k = k + 1
    ipos1 = ipos + m
    ipos2 = ipos + j + 2 - m
    ipos3 = ipos1 + ndim/2
    ipos4 = ipos2 + ndim/2

    if ( seed_type == 2 ) then
      factor = (-one)**((HOsp_2j(ipos1) - HOsp_2mj(ipos1))/2)
    else
      factor = one
    endif
 
    U(ipos1,ipos1) =  sqrt(1.d0 - v2c(k))
    U(ipos2,ipos2) =  U(ipos1,ipos1)
    U(ipos3,ipos3) =  sqrt(1.d0 - v2c(k+ndim/4))
    U(ipos4,ipos4) =  U(ipos3,ipos3)
  
    V(ipos1,ipos2) =  factor * sqrt(v2c(k))
    V(ipos2,ipos1) = -V(ipos1,ipos2)
    V(ipos3,ipos4) =  factor * sqrt(v2c(k+ndim/4))
    V(ipos4,ipos3) = -V(ipos3,ipos4)
  enddo
  ipos = ipos + j + 1
enddo

end subroutine generate_wavefunction_BCS    

!------------------------------------------------------------------------------!
! subroutine generate_wavefunction_slater                                      !
!                                                                              !
! This subroutine generates Bogoliubov matrices U and V corresponding to a     ! 
! Slater determinant. This is achieved by occupying the first (nprot,nneut)    !
! single-particle states by storing order. A subsequent random unitary         !
! transformation can then be performed just after to generate a random Slater  !
! determinant.                                                                 !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        nprot = number of protons                                             !
!        npneut= number of neutrons                                            !
!                                                                              !
! Output: U,V = Bogoliubov matrices                                            !
!------------------------------------------------------------------------------!
subroutine generate_wavefunction_slater(nprot,nneut,U,V,ndim)

integer, intent(in) :: ndim, nprot, nneut
real(r64), dimension(ndim,ndim), intent(out) :: U, V                
integer, dimension(ndim) :: idx    
integer :: i, k, l
real(r64) :: ranval

!!! Builds a fully empty state 
U = zero
V = zero

do i = 1, ndim
  U(i,i) = one
enddo

!!! Occupies the sp states for protons. The srategy is to generate a random
!!! integer value between 1 and (ndim/2 - iteration +1) (=> k), then look in an 
!!! array the state to occupy (=> l), then occupy the state, then update the
!!! array to remove the possibility to find the same index at next iteration
do i = 1, ndim/2
  idx(i) = i
enddo

do i = 1, nprot  
  call random_number(ranval)
  k = 1 + int((ndim/2-i+1) * ranval)
  l = idx(k)

  U(l,l) = zero
  V(l,l) = one 

  do l = k+1, ndim/2
    idx(l-1) = idx(l)
  enddo
enddo


!!! Occupies the sp states for neutrons. Same strategy.
do i = 1, ndim/2
  idx(i) = i
enddo

do i = 1, nneut  
  call random_number(ranval)
  k = 1 + int((ndim/2-i+1) * ranval)
  l = idx(k)

  U(l+ndim/2,l+ndim/2) = zero
  V(l+ndim/2,l+ndim/2) = one 

  do l = k+1, ndim/2
    idx(l-1) = idx(l)
  enddo
enddo

end subroutine generate_wavefunction_slater

!------------------------------------------------------------------------------!
! subroutine generate_unitary_matrix                                           !
!                                                                              !
! Generates a unitary matrix by diagonalizing a randomly generated symmetric   !
! matrix. The matrix returned is the matrix made of the eigenvectors of the    !
! diagonalized matrix. The matrix is generated with the appropriate block      !
! structure to preserve the symmetries of the seed wave function (seed_type).  !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        opt_bloc = option to choose the block structure of the matrix.        !
!                                                                              !
! Output: A = random unitary matrix with the appropriate block structure.      !
!------------------------------------------------------------------------------!
subroutine generate_unitary_matrix(opt_bloc,A,ndim)

integer, intent(in) :: ndim, opt_bloc
real(r64), dimension(ndim,ndim), intent(out) :: A
integer :: i, j, k, m, n, nbloc, info
integer, dimension(4) ::  tdim(4), tpar(4), tmax(4), tmin(4)
real(r64) :: ranval
real(r64), dimension(:), allocatable :: wev, work
real(r64), dimension(:,:), allocatable :: B

A = zero

!!! Determines the number of block to use and their dimensions. 
tdim = 0
tpar = 0
tmin = 1
tmax = ndim

select case (opt_bloc)
 case(0)
   nbloc = 1
   tdim(1) = ndim

 case(4)
   do i = 1, ndim
     if ( (-1)**HOsp_l(i) == 1 ) tdim(1) = tdim(1) + 1
   enddo
   tdim(2) = ndim - tdim(1)
   if ( (tdim(1) /= 0) .and. (tdim(2) /= 0) ) then
     nbloc = 2
     tpar(1) =  1    
     tpar(2) = -1    
   else
     nbloc = 1
     tdim(1) = ndim
   endif 

 case(5,7)
   nbloc = 2
   tdim(1) = ndim/2
   tdim(2) = ndim/2
   tmax(1) = ndim/2
   tmin(2) = 1 + ndim/2

 case(6,8)
   do i = 1, ndim/2
     if ( (-1)**HOsp_l(i) == 1 ) tdim(1) = tdim(1) + 1
   enddo
   tdim(2) = ndim/2 - tdim(1)
   if ( (tdim(1) /= 0) .and. (tdim(2) /= 0) ) then
     nbloc = 4
     tdim(3) = tdim(1)
     tdim(4) = tdim(2)
     tpar(1) =  1    
     tpar(2) = -1    
     tpar(3) =  1    
     tpar(4) = -1   
     tmax(1) = ndim/2
     tmax(2) = ndim/2
     tmin(3) = 1 + ndim/2
     tmin(4) = 1 + ndim/2
   else 
     nbloc = 2
     tdim(1) = ndim/2 
     tdim(2) = ndim/2 
     tmax(1) = ndim/2
     tmin(2) = 1 + ndim/2
   endif
end select

!!! Computes each bloch by diagonalizing a symmetric matrix of the same dimen-
!!! sion. As the single-particle states are not in any particular order regard-
!!! ing the parity, one has to keep track of it. 
do k = 1, nbloc

  allocate ( B(tdim(k),tdim(k)), wev(tdim(k)), work(3*tdim(k)-1) ) 
  do j = 1, tdim(k)
    do i = 1, tdim(k)
      call random_number(ranval)
      B(i,j) = ranval
      B(j,i) = ranval
    enddo 
  enddo 
  call dsyev('v','L',tdim(k),B,tdim(k),wev,work,3*tdim(k)-1,info)

  m = 1
  n = 1
  do j = 1, ndim
    if ( (tpar(k)*(-1)**HOsp_l(j) >= 0) .and. (j >= tmin(k)) .and. &
         (j <= tmax(k)) ) then
      do i = 1, ndim
        if ( (tpar(k)*(-1)**HOsp_l(i) >= 0) .and. (i >= tmin(k)) .and. &
           (i <= tmax(k)) ) then
          A(i,j) = B(m,n)
          m = m + 1
        endif
      enddo
      m = 1
      n = n + 1
    endif
  enddo
  deallocate(B, wev, work)

enddo

end subroutine generate_unitary_matrix 

!------------------------------------------------------------------------------!
! subroutine read_wavefunction                                                 !
!                                                                              !
! Reads the wave function from "initial_wf.bin/txt" if the option has been     !
! chosen in the input parameters (seed_type = 1). The wave function is either  !
! read as binary (.bin) or text (.txt) file.                                   !
!                                                                              !
! seed_text = 0 or 2 reads a binary file                                       !
!           = 1 or 3   "   " text    "                                         !
!                                                                              !
! The format of the binary file is                                             !
! Hosh_dim                         [number of shells in the model space]       !
! (HOsh_na(i), i=1,HOsh_dim)       [list of shells in the model space]         !
! bogo_label                       [label identifying the wave function]       !
! ((bogo_U0(j,i), j=1,HOsp_dim), i=1,HOsp_dim)     [Bogoliubov matrix U]       !
! ((bogo_V0(j,i), j=1,HOsp_dim), i=1,HOsp_dim)     [Bogoliubov matrix V]       !
!                                                                              !
! In the case of a text file, the file is read as a single column.             !
!------------------------------------------------------------------------------!
subroutine read_wavefunction            

integer :: i, j, icheck, HOsh_dim0       
integer, dimension(:), allocatable :: HOsh_na0
character(4) :: filetype  
character(11) :: fileform  
character(14) :: filename  
logical :: is_exist, is_binary

!!! Determines the name of the file (binary or text file)
select case (seed_text)
  case(0,2)
    is_binary = .true.
  case(1,3)
    is_binary = .false.
end select

if ( is_binary ) then
  fileform = 'unformatted'
  filetype = '.bin'
else
  fileform = 'formatted'
  filetype = '.txt'
endif

filename = 'initial_wf' // filetype

!!! Checks if the initial wave function exists, otherwise stops the run 
inquire (file=filename, exist=is_exist)
if ( is_exist .eqv. .false. ) then
  print '(/,"The file containing the seed wave function is missing.")'
  stop 
endif 

!!! Opens the file
open(utw, file=filename, status='old', action='read', form=fileform)

!!! Reads the model space 
if ( is_binary ) then
  read(utw) HOsh_dim0
  allocate(HOsh_na0(HOsh_dim0))
  read(utw) (HOsh_na0(i), i=1,HOsh_dim0)
else
  read(utw,*) HOsh_dim0
  allocate(HOsh_na0(HOsh_dim0))
  do i = 1, HOsh_dim0
    read(utw,*) HOsh_na0(i)
  enddo
endif

!!! Stops the run if the model spaces of the wave func. and interaction differ    
icheck = 0
if ( HOsh_dim0 /= HOsh_dim ) icheck = icheck + 1
do i = 1, min(HOsh_dim,HOsh_dim0)
  if ( HOsh_na0(i) /= HOsh_na(i) ) icheck = icheck + 1
enddo
if ( icheck /= 0 ) then
  print '(/,"The model space of the seed wave function is not consistent", & 
        & " with the one of the interaction.")'
  print*, 'Inter:', HOsh_dim, (HOsh_na(i), i=1,HOsh_dim)
  print*, 'State:', HOsh_dim0, (HOsh_na0(i), i=1,HOsh_dim0)
  stop 
endif

!!! Reads the wave function                                                  
if ( is_binary ) then
  read(utw) bogo_label   
  read(utw) ((bogo_U0(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
  read(utw) ((bogo_V0(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
else
  read(utw,*) bogo_label   
  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      read(utw,*) bogo_U0(j,i)
    enddo
  enddo
  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      read(utw,*) bogo_V0(j,i)
    enddo
  enddo
endif

close (utw, status='keep')

deallocate(HOsh_na0)

end subroutine read_wavefunction    

!------------------------------------------------------------------------------!
! subroutine write_wavefunction                                                !
!                                                                              !
! Writes the wave function in a file:                                          !
! - "final_wf.bin/txt"        after  the iterative procedure                   !
! - "intermediate_wf.bin/txt" during  "      "        "                        !
! The wave function is either written in a binary (.bin) or text (.txt) file.  !
! Before being written, the wavefunctions is being assigned a label, which is  !
! an integer between 0 and 10^18 to identify uniquely the state (used in       !
! the GCM).                                                                    !
!                                                                              !
! seed_text = 0 or 3 writes a binary file                                      !
!           = 1 or 2   "    " text    "                                        !
!                                                                              !
! Input: iopt = 0 writes "final_wf.bin/txt"                                    !
!             = 1   "    "intermediate_wf.bin/txt"                             !
!        ener = particle-number projected energy (optional)                    !
!------------------------------------------------------------------------------!
subroutine write_wavefunction(iopt,ener)           

integer, intent(in) :: iopt
real(r64), intent(in) :: ener 
integer :: i, j
real(r64) :: absener
character(4) :: filetype  
character(11) :: fileform  
character(19) :: filename  
logical :: is_binary

!!! Determines the label of the state based on a random number and the 
!!! particle-number projected energy.
absener = abs(ener)

do while (absener >= one) 
  absener = absener / 10.d0
enddo

bogo_label = int(absener*(1.0d18),i64)

!!! Opens the file depending if intermediate or final writing and if binary
!!! or text writing
select case (seed_text)
  case(0,3)
    is_binary = .true.
  case(1,2)
    is_binary = .false.
end select

if ( is_binary ) then
  fileform = 'unformatted'
  filetype = '.bin'
else
  fileform = 'formatted'
  filetype = '.txt'
endif

if ( iopt == 0 ) then
  filename = 'final_wf' // filetype
else
  filename = 'intermediate_wf' // filetype
endif

open(utw, file=filename, status='replace', action='write',form=fileform)

!!! Writes the wave function and model space
if ( is_binary ) then
  write(utw) HOsh_dim
  write(utw) (HOsh_na(i), i=1,HOsh_dim)
  write(utw) bogo_label
  write(utw) ((bogo_U0(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
  write(utw) ((bogo_V0(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
else
  write(utw,*) HOsh_dim
  do i = 1, HOsh_dim
    write(utw,*) HOsh_na(i)
  enddo
  write(utw,*) bogo_label
  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      write(utw,*) bogo_U0(j,i)
    enddo
  enddo
  do i = 1, HOsp_dim
    do j = 1, HOsp_dim
      write(utw,*) bogo_V0(j,i)
    enddo
  enddo
endif

close(utw, status='keep')

end subroutine write_wavefunction    

!------------------------------------------------------------------------------!
! subroutine construct_canonical_basis                                         !
!                                                                              !
! Constructs the canonical basis of the Bogoliubov quasiparticle states with   !
! matrices U,V. Also, computes the number of fully occupied and empty states   !
! in this basis, which is needed to compute the overlap later.                 !
! The algorithm used is the following:                                         !
! 1) We first diagonalize rho                                                  !
! 2) We transforms kappa in the basis that diagonalizes rho                    !
! 3) If they are non-canonical block in kappa (due to degenerate eigenvalues   !
!    of rho), we perform a Schur decomposition on these blocks Given kappa is  !
!    real, it will automatically put it in its canonical form.                 !
!                                                                              !
! Remark: the determination of non-canoncial block can only be numerical and   !
!         therefore may fail in some cases.                                    !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        U,V  Bovoliubov matrices                                              !
!                                                                              !
! Output: zUc,zVc = U and V in the canonical basis                             !
!         zDc = unitary transformation to change the basis                     !
!         ovacc = overlap with the bare vacuum of the fully-paired part of wf  !
!         nocc,nemp = number of fully occupied/empty states in canon. basis    !
!------------------------------------------------------------------------------!
subroutine construct_canonical_basis (U,V,zUc,zVc,zDc,ovacc,nocc,nemp,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: U, V
integer, intent(out) :: nocc, nemp
real(r64), intent(out) :: ovacc       
complex(r64), dimension(ndim,ndim), intent(out) :: zUc, zVc, zDc
integer :: i, j, k, m, info, sdim, ialloc=0
integer, dimension(ndim) :: submax, subdim
real(r64) :: occu_u, occu_v, occu_v2, eps
real(r64), dimension(ndim) :: eigen_rho
real(r64), dimension(3*ndim-1) :: work1   
real(r64), dimension(:), allocatable :: wr, wi, work2
real(r64), dimension(ndim,ndim) :: rho, kappa, Dc, rhoc, kappac, A1, A2 
real(r64), dimension(:,:), allocatable :: B1, vs 
logical, dimension(:), allocatable :: bwork

!!! Cutoff for occupied single-particle states
if ( seed_occeps < tiny(zero) ) then
  eps = 1.0d-8
else
  eps = seed_occeps
endif

!!! Builds the density rho, then diagonalizes it, then transforms rho in the  
!!! basis where it is diagonal.
call dgemm('n','t',ndim,ndim,ndim,one,V,ndim,V,ndim,zero,rho,ndim)

Dc = rho
call dsyev('V','U',ndim,Dc,ndim,eigen_rho,work1,3*ndim-1,info)

call dgemm('t','n',ndim,ndim,ndim,one,Dc,ndim,rho,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,Dc,ndim,zero,rhoc,ndim)

!!! Counts the number of occupied/empty single-particle states. To determine if 
!!! a state is occupied/empty, we use a cutoff: eps.
!!! Then computes the overlap of the fully paired part with the bare vacuum. As
!!! the quantity occu_u can be quite small, we build \sum log(occu_u) rather 
!!! than \product occu_u.
nocc = 0
nemp = 0
ovacc = 0.0d0

do i = 1, ndim
  occu_v2 = abs(rhoc(i,i))
  occu_u  = sqrt(abs(1.d0 - occu_v2))  
  if ( occu_v2 >= 1.0d0 - eps ) then
    nocc = nocc + 1    
  elseif ( (occu_v2 <= eps) .and. (seed_allemp == 0) ) then
    nemp = nemp + 1
  else
    ovacc = ovacc + log(occu_u) 
  endif
enddo

!!! Builds the density kappa, then builds it in the basis that diagonlizes rho.
!!! At this point, there is no guaranty that kappa will be in its canonical form
!!! yet. There could be degeneracies in the v2,u2 that would prevent a full     
!!! canonical structure.
call dgemm('n','t',ndim,ndim,ndim,one,V,ndim,U,ndim,zero,kappa,ndim) 
call dgemm('t','n',ndim,ndim,ndim,one,Dc,ndim,kappa,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,Dc,ndim,zero,kappac,ndim)

!!! Checks if the kappa is already in its canonical form by counting the number
!!! of non-zero marix in a column of kappac. If it is greater than one (i.e. 
!!! the dimension of the subspace subdim > 2), we need to further reduce kappac.
submax = 0

do i = 1+nemp, ndim-nocc
  submax(i) = i 
  do j = i+1, ndim-nocc
    if ( abs(kappac(j,i)) > sqrt(eps)*1.0d-2 ) submax(i) = j
  enddo
enddo

!!! Tries to determine the dimension of each subspace
k = 1 + nemp

do while ( k < ndim - nocc )
  j = submax(k)
  subdim = 0
  subdim(k:j) = submax(k:j)
  m = maxval(subdim)
  submax(k:j) = m
  if ( m == j ) k = j + 1
enddo

k = 1 + nemp

do while ( k < ndim - nocc )
  subdim(k) = submax(k) - k + 1
  if ( (subdim(k) > 2) .and. (mod(subdim(k),2) /= 0) ) then
    !print*, k, subdim(k)
    subdim(k) = 2
    print*,"Warning: subspace of odd dimension when building kappa canonical. &
          & Intenting the calculation assuming dim = 2."
  endif
  k = k + subdim(k)
enddo

!!! Transforms the degenerates subspaces in their Schur form, which automat-
!!! ically puts them in their canonical form (kappa being real)
A1 = zero
do i = 1, ndim
  A1(i,i) = one
enddo

k = 1 + nemp

do while ( k < ndim - nocc )
  m = subdim(k)
  !!! Schur factorization of this block 
  if ( m > 2 ) then
    !print*,"Diangostic purpose (Schur), dim = :", m
    allocate(B1(m,m), wr(m), wi(m), work2(3*m), bwork(m), vs(m,m), &
             stat=ialloc)
    if ( ialloc /= 0 ) stop 'Error during allocation for Schur factorization'
    
    B1(1:m,1:m) = kappac(k:k+m-1,k:k+m-1)
    call dgees('V','N',lapack_sel,m,B1,m,sdim,wr,wi,vs,m,work2,3*m,bwork,info)
  
    ! store the transformation for this block   
    A1(k:k+m-1,k:k+m-1) = vs(1:m,1:m)  
    kappac(k:k+m-1,k:k+m-1) = B1(1:m,1:m)  
    deallocate(B1,wr,wi,work2,bwork,vs)
  endif
  
  k = k + m
enddo    

!!! Update the transformation D and recomptues kappa in the canonical basis
!!! if necessary
if ( maxval(subdim) > 2 ) then
  call dgemm('n','n',ndim,ndim,ndim,one,Dc,ndim,A1,ndim,zero,A2,ndim)
  Dc = A2
endif

zDc = zone * Dc

!!! Construct U and V in the canonical basis. Note that as we know rho/kappa in
!!! the canonical basis, we do not need to use D or construct C (à la BMZ).
!!! Nevertheless, we have to be careful about the sign of the matrix elements
!!! to reproduce kappac.
zUc = zzero
zVc = zzero

k = 0

do i = 1, ndim
  occu_v2 = abs(rhoc(i,i))
  occu_v = sqrt(occu_v2)

  zUc(i,i) = sqrt(abs(1.d0 - occu_v2))

  if ( occu_v2 > 1.0d0 - eps ) then   ! occupied states
    zVc(i,i) = occu_v * zone          
  elseif ( occu_v2 <= eps ) then      ! empty states
    zVc(i,i) = occu_v * zone          
  else                                ! paired states
    if ( mod(k,2) == 0 ) then
      zVc(i,i+1) = sign(occu_v,kappac(i,i+1)) * zone 
    else
      zVc(i,i-1) = sign(occu_v,kappac(i,i-1)) * zone
    endif
    k = k + 1
  endif
enddo

end subroutine construct_canonical_basis

!------------------------------------------------------------------------------!
! subroutine block_quasiparticle                                               !
!                                                                              ! 
! Performs the blocking of quasiparticle states by exchanging the colmuns of   ! 
! Bogoliubov matrices U and V.                                                 ! 
! For index 1 we have : (Um1,Vm1) => (Vm1^*,Um1^*) as described in equation    ! 
! (7.21) of the book by Ring and Schuck (ISBN: 978-3-540-21206-5).             ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        nidx = the index of the qp to block                                   !
!        UR,VR = Bovoliubov matrices                                           !
!------------------------------------------------------------------------------!
subroutine block_quasiparticle (nidx,U,V,ndim)

integer, intent(in) :: nidx, ndim
real(r64), dimension(ndim,ndim), intent(inout) :: U, V
integer :: i
real(r64), dimension(ndim,ndim) :: Ub, Vb

!!! Exits the routine if no blocking is required
if ( nidx == 0 ) return

!!! Checks if the blocking index is consistent with the model space 
if ( nidx > ndim ) then
 print '(1a,1x,1i5,1a,1i5)', 'The blocking index = ',nidx,' is greater than &
       & the dimension of the single-particle space = ',ndim  
 stop 
endif 

!!! Exchanges the nidx-th column of U and V                                
Ub = U
Vb = V
do i = 1, ndim
  U(i,nidx) = Vb(i,nidx)
  V(i,nidx) = Ub(i,nidx)
enddo

end subroutine block_quasiparticle

!------------------------------------------------------------------------------!
! subroutine print_wavefunction                                                !
!                                                                              !
! Prints the informations about the seed wave function.                        !
!                                                                              !
! Input: nocc = number of single-particle states fully occupied (canon. basis) !
!------------------------------------------------------------------------------!
subroutine print_wavefunction(nocc)

integer, intent(in) :: nocc
integer :: i, j, nline
character(len=:), allocatable :: info_type
character(len=*), parameter :: format1 = "(1a20,1x,1i5,3x,1a)", &
                               format2 = "(1a20,5i6)", &
                               format3 = "(1a20,1x,1i5)"

nline = blocking_dim/5 + 1
if ( mod(blocking_dim,5) == 0 ) nline = nline - 1

!!! Additional information on the type of seed
select case (seed_type)
  case (0)
    info_type = '(general, unrestricted)'
  case (1)
    info_type = '(read from file)'
  case (2)
    info_type = '(spherical BCS, good P + J, separate N/Z)'
  case (3)
    info_type = '(axial BCS, good P + Jz, separate N/Z)'
  case (4)
    info_type = '(general, good P)'
  case (5)
    info_type = '(general, separate N/Z)'
  case (6)
    info_type = '(general, good P, separate N/Z)'
  case (7)
    info_type = '(Slater, unrestricted)'
  case (8)
    info_type = '(Slater, good P)'
  case (9)
    info_type = '(Slater, good P + Jz, HO)'
  case default
    info_type = '(not recognized)'
end select 

print '(/,60("%"),/,24x,"WAVE FUNCTION",23x,/,60("%"),//, &
      & 3x,"Description",9x,"Value",/,28("-"))'
print format1, 'Initial type of seed', seed_type, info_type  
print format2, 'Number of qp blocked', blocking_dim
if ( blocking_dim /= 0 ) then 
  print format2, 'Indices for blocking',(blocking_id(i), i=1, &
                                         min(5,blocking_dim))
  do i = 1, nline-1
    print format3, '                    ',(blocking_id(j),j=1+i*5, &
                                           min((i+1)*5,blocking_dim))
  enddo
endif
print format3, 'Total number parity ', bogo_nparity
print format3, 'No. of fully occ. sp', nocc

end subroutine print_wavefunction  

!------------------------------------------------------------------------------!
! subroutine calculate_densities                                               !
!                                                                              !
! Calculates the (possibly non-diagonal) densities according to                !
!     rhoLR =   V_R^* V_L^T                                                    !
!   kappaLR =   V_R^* U_L^T                                                    !
! kappaRL^* = - U_R^* V_L^T = V_L U_R^dagger                                   !
! Be careful that what is called kappaRL in the code is actually kappaRL^*     !
!                                                                              !
! Actually, the projection routine will inject the complex conjugate of Utilde !
! and Vtilde because we use in that case the formulae (that can be found in    !
! T.R. Rodriguez master thesis)                                                !
!     rhoLR =   Vtilde V_L^T                                                   !
!   kappaLR =   Vtilde U_L^T                                                   !
! kappaRL^* = - Utilde V_L^T = V_L Utilde^T                                    !
! where                                                                        !
! Utilde =  U_L^* + V_L Aphi                                                   !
! Vtilde =  V_L^* + U_L Aphi                                                   !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        UL,VL = Left Bogoliubov matrices                                      !
!        UR,VR = Right Bogoliubov matrices                                     !
!                                                                              !
! Output: rhoLR,kappaLR,kappaRL = densities                                    !
!------------------------------------------------------------------------------!
subroutine calculate_densities(UL,VL,UR,VR,rhoLR,kappaLR,kappaRL,ndim)
                                                                                 
integer, intent(in) :: ndim                                                                  
complex(r64), dimension(ndim,ndim), intent(in) :: UL, VL, UR, VR
complex(r64), dimension(ndim,ndim), intent(out) :: rhoLR, kappaLR, kappaRL
complex(r64), dimension(ndim,ndim) :: URc, VRc

URc = conjg(UR)                                                                 
VRc = conjg(VR)                                                                 
                                                                                 
call zgemm('n','t',ndim,ndim,ndim,zone,VRc,ndim,VL,ndim,zzero,rhoLR,ndim)             
call zgemm('n','t',ndim,ndim,ndim,zone,VRc,ndim,UL,ndim,zzero,kappaLR,ndim)             
call zgemm('n','t',ndim,ndim,ndim,zone,VL,ndim,URc,ndim,zzero,kappaRL,ndim)             
                                                                                 
end subroutine calculate_densities

!------------------------------------------------------------------------------!
! subroutine calculate_densities_real                                          !
!                                                                              !
! Calculates the real diagonal densities according to                          !
!     rhoRR = V_R V_R^T                                                        !
!   kappaRR = V_R U_R^T                                                        !
! Here all quantities are real so no need of complex conjugate.                !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        UR,VR = Bovoliubov matrices                                           !
!                                                                              !
! Output: rhoRR,kappaRR = densities                                            !
!------------------------------------------------------------------------------!
subroutine calculate_densities_real(UR,VR,rhoRR,kappaRR,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: UR, VR
real(r64), dimension(ndim,ndim), intent(out) :: rhoRR, kappaRR

call dgemm('n','t',ndim,ndim,ndim,one,VR,ndim,VR,ndim,zero,rhoRR,ndim)
call dgemm('n','t',ndim,ndim,ndim,one,VR,ndim,UR,ndim,zero,kappaRR,ndim)

end subroutine calculate_densities_real

!------------------------------------------------------------------------------!
! subroutine calculate_spatial_density                                         !
!                                                                              !
! Calculates the spatial one-body densities for proton and neutrons separately !
! and writes it in a file.                                                     !
! The spatial density reads                                                    !
!                                                                              !
! The spatial density integrated over the angles reads                         !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        rhoRR = one-body density                                              !
!        rhoNZ = one-body density projected on NZ                              !
!------------------------------------------------------------------------------!
subroutine calculate_spatial_density(rhoRR,rhoNZ,ndim)

integer, intent(in) ::  ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoRR, rhoNZ
integer :: a, b, na, nb, la, lb, mla, mlb, ja, jb, mja, mjb, ms, hdim, &
           r1, r2, r3, r1min, r1max, r2min, r2max, r3min, r3max
real(r64) :: x, y, z, r, theta, phi, dr1, dr2, dr3, Rnla, Rnlb, cga, cgb, &   
             fac, rho_p, rho_n, rhoNZ_p, rhoNZ_n
complex(r64) :: Ylma, Ylmb, angular
character(50) :: filename1, filename2, filename3
character(len=*), parameter :: format1 = "(1f7.3,4f10.5)", &
                               format2 = "(3f8.4,4f10.5)"

hdim = HOsp_dim/2

!!! Name for the file where the densities are written
filename1 = adjustl('spatial_density_R.dat')
filename2 = adjustl('spatial_density_RThetaPhi.dat')
filename3 = adjustl('spatial_density_XYZ.dat')

!!! If dens_spatial = 1: r1 = r, r2 = theta, r3 = phi
!!! If dens_spatial = 2: r1 = x, r2 = y    , r3 = z  
if ( dens_spatial <= 2 ) then
  r1min = 0 
  r1max = dens_nr(1) - 1
  r2min = 0 
  r2max = dens_nr(2) - 1
  r3min = 0 
  r3max = dens_nr(3) - 1
  dr1 = dens_dr(1)
  dr2 =   pi / (r2max+1)
  dr3 = 2*pi / (r3max+1)
else 
  r1min = -dens_nr(1) + 1
  r1max =  dens_nr(1)
  r2min = -dens_nr(2) + 1
  r2max =  dens_nr(2)
  r3min = -dens_nr(3) + 1
  r3max =  dens_nr(3)
  dr1 = dens_dr(1)
  dr2 = dens_dr(2)
  dr3 = dens_dr(3)
endif

!!! Computes the spatial one-body density integrated over the angular part
if ( dens_spatial <= 2 ) then  
  open(utd, file=filename1, status='replace', action='write', form='formatted')        
  write(utd,*) "            unprojected          projected   "
  write(utd,*) "   r      rho_p     rho_n     rho_p     rho_n"

  do r1 = r1min, r1max
    r = r1 * dr1 

    rho_p = zero
    rho_n = zero
    rhoNZ_p = zero
    rhoNZ_n = zero

    do a = 1, hdim
      na = HOsp_n(a)
      la = HOsp_l(a)
      ja = HOsp_2j(a)
      mja= HOsp_2mj(a)
      Rnla = radial_function(na,la,r)

      do b = a, hdim
        nb = HOsp_n(b)
        lb = HOsp_l(b)
        jb = HOsp_2j(b)
        mjb= HOsp_2mj(b)
        Rnlb = radial_function(nb,lb,r)
     
        if ( la /= lb ) cycle   
        if ( ja /= jb ) cycle   
        if ( mja /= mjb ) cycle   
 
        fac = (2.0d0 - kdelta(a,b)) * Rnla * Rnlb
 
        rho_p = rho_p + fac * real(rhoRR(b,a))
        rho_n = rho_n + fac * real(rhoRR(b+hdim,a+hdim))
        rhoNZ_p = rhoNZ_p + fac * real(rhoNZ(b,a))
        rhoNZ_n = rhoNZ_n + fac * real(rhoNZ(b+hdim,a+hdim))
 
      enddo
    enddo

    write(utd,format1) r, rho_p/(4*pi), rho_n/(4*pi), rhoNZ_p/(4*pi), & 
                          rhoNZ_n/(4*pi) 
  enddo

  close(utd)
endif

if ( dens_spatial == 1 ) return

!!! Computes the spatial one-body density 
if ( dens_spatial == 2 ) then
  open(utd, file=filename2, status='replace', action='write', form='formatted')        
  write(utd,*) "                             unprojected          projected   "
  write(utd,*) "   r      theta   phi      rho_p     rho_n     rho_p     rho_n"
else
  open(utd, file=filename3, status='replace', action='write', form='formatted')        
  write(utd,*) "                             unprojected          projected   "
  write(utd,*) "   r       y       z       rho_p     rho_n     rho_p     rho_n"
endif

do r1 = r1min, r1max
  do r2 = r2min, r2max
    do r3 = r3min, r3max

      if ( dens_spatial == 2 ) then
        r = r1 * dr1 
        theta = r2 * dr2 + dr2/2
        phi = r3 * dr3 + dr3/2
      else
        x = r1 * dr1 + dr1/2
        y = r2 * dr2 + dr2/2
        z = r3 * dr3 + dr3/2

        r = sqrt( x**2 + y**2 + z**2 )
        theta = acos( z / r )
        phi = acos( x / sqrt(x**2 + y**2) )
      endif

      rho_p = zero
      rho_n = zero
      rhoNZ_p = zero
      rhoNZ_n = zero

      do a = 1, hdim
        na = HOsp_n(a)
        la = HOsp_l(a)
        ja = HOsp_2j(a)
        mja= HOsp_2mj(a)
        Rnla = radial_function(na,la,r)

        do b = a, hdim
          nb = HOsp_n(b)
          lb = HOsp_l(b)
          jb = HOsp_2j(b)
          mjb= HOsp_2mj(b)
          Rnlb = radial_function(nb,lb,r)
     
          fac = (2.0d0 - kdelta(a,b)) * Rnla * Rnlb
     
          angular = zzero
          do mla = -la, la
            ms = mja - 2*mla
            mlb = (mjb - mja + 2*mla)/2
            if ( abs(mlb) > lb ) cycle

            Ylma = conjg(spherharmonic(la,mla,theta,phi))
            Ylmb = spherharmonic(lb,mlb,theta,phi)
            call ClebschGordan(1,2*la,ja,ms,2*mla,mja,cga)
            call ClebschGordan(1,2*lb,jb,ms,2*mlb,mjb,cgb)
            angular = angular + Ylma * Ylmb * cga * cgb 
          enddo
         
          rho_p = rho_p + fac * real(rhoRR(b,a) * angular)
          rho_n = rho_n + fac * real(rhoRR(b+hdim,a+hdim) * angular)
          rhoNZ_p = rhoNZ_p + fac * real(rhoNZ(b,a) * angular)
          rhoNZ_n = rhoNZ_n + fac * real(rhoNZ(b+hdim,a+hdim) * angular)
        enddo
      enddo

      if ( dens_spatial == 2 ) then 
        write(utd,format2) r, theta, phi, rho_p, rho_n, rhoNZ_p, rhoNZ_n
      else 
        write(utd,format2) x, y, z, rho_p, rho_n, rhoNZ_p, rhoNZ_n
      endif
    enddo
  enddo
enddo

close(utd)

end subroutine calculate_spatial_density

END MODULE Wavefunctions 
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
