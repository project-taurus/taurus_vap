!==============================================================================!
! MODULE Constraints                                                           !
!                                                                              !
! This module contains the variables and routines related to the constraints   !
! applied during the minimization.                                             !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_constraints                                                 !
! - subroutine adjust_constraints                                              !
! - subroutine update_constraints_qpbasis                                      !
! - subroutine calculate_lagrange_multipliers                                  !
! - subroutine cholesky_factorization                                          !
! - subroutine evolve_wavefunction                                             !
! - subroutine check_symmetries                                                !
!==============================================================================!
MODULE Constraints

use Wavefunctions
use Operators
use Fields, only: factor_delta
use Projection, only: proj_Mphip, proj_Mphin
                     

implicit none

!!! Number of constraints and options for specific cases
integer, parameter :: constraint_types=27, & ! Number of constraint types     
                      constraint_max=42      ! Maximum number of constraints
integer :: constraint_dim=0,  & ! Number of constraints to be applied
           constraint_pair=0, & ! Number of pairing constraints
           constraint_switch(constraint_types), & ! Activate constraints
           constraint_id(constraint_types), & ! Remember order of constraints
           enforce_NZ, & ! Apply constraint on N and Z even if PNVAP
           opt_betalm    ! Constraint on (beta,gamma) instead of (Q20,Q22)

!!! Matrix elements and values of constraints
real(r64), dimension(:), allocatable :: constraint_HO,  & ! Mat. elem. (HO)
                                        constraint_val, & ! Values 
                                        constraint_20,  & ! Mat. elem. 20 (QP) 
                                        constraint_11     ! Mat. elem. 11 (QP) 
real(r64) :: constraint_eps, & ! Tolerance for constraints convergence 
             constraint_read(constraint_types,2) ! Values read as inputs

!!! Lagrange multipliers                        
real(r64), dimension(:), allocatable :: lagrange_lambda0, & ! Lag. mult. before
                                        lagrange_lambda1    ! Lag. mult. after
!!! Others
real(r64), dimension(:,:), allocatable, private :: Lchol ! Choleksy L matrix

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_constraints                                                   !
!                                                                              !
! Constructs the matrix containing all the constraints. This is based on the   !
! input parameters but also on the good quantum numbers of the seed wave func- !
! tion.                                                                        !
!                                                                              !
! For betal_lm constraints (opt_blm = 1), the definition is                    !
!   Qlm = beta_lm * coef                                                       !
! whereas for beta/gamma constraints (opt_blm = 2), it is                      !
!   Q20 = beta * cos(gamma) * coef                                             !
!   Q22 = beta * sin(gamma) * coef * sqrt(2)                                   !
! where                                                                        !
!   coef = 3/(4*pi) * r0**2 * A**(5/3)                                         !
!                                                                              !
! For constraint_switch = 3, the values read in input correspond to the        !
! isoscalar/isovector constraints and are transformed into the proton/neutron  !
! values as                                                                    !
! Q_p = (Q_is - Q_iv) / 2                                                      !
! Q_n = (Q_is + Q_iv) / 2                                                      !
!                                                                              !
! For the radius, the values read are squared as we compute <r^2> and we note  !
! that it is also necessary to divide by the expectation values <Z>, <N> or    !
! <A> (done when calculating the expectation values).                          !
!------------------------------------------------------------------------------!
subroutine set_constraints

integer :: i, j, m, n, p, ialloc=0
integer(i64) :: k
real(r64) :: constraint_beta(2), constraint_gamma(2), value_is, value_iv, sig
real(r64), dimension(:), allocatable :: Q

!!! Lifts the constraint on the number of protons if PNVAP and no forcing, or
!!! if the seed wave function is a Slater determinant without pairing constr.
if ( ((proj_Mphip > 1) .and. (enforce_NZ == 0)) .or. is_good_Z ) then
  constraint_switch(1) = 0
  constraint_read(1,1) = 0.0d0
endif

if ( constraint_switch(1) == 0 ) constraint_dim = constraint_dim - 1

!!! Lifts the constraint on the number of neutrons if PNVAP and no forcing, or
!!! if the seed wave function is a Slater determinant without pairing constr.
if ( ((proj_Mphin > 1) .and. (enforce_NZ == 0)) .or. is_good_N ) then
  constraint_switch(2) = 0
  constraint_read(2,1) = 0.0d0
endif

if ( constraint_switch(2) == 0 ) constraint_dim = constraint_dim - 1

!!! Initializes the arrays knowing the final number of constraints
allocate ( Q(HOsp_dim2), Lchol(HOsp_dim,HOsp_dim),  &
           constraint_HO(constraint_dim*HOsp_dim2), &
           constraint_20(constraint_dim*HOsp_dim2), &
           constraint_11(constraint_dim*HOsp_dim2), &
           constraint_val(constraint_dim),   &  
           lagrange_lambda0(constraint_dim), & 
           lagrange_lambda1(constraint_dim), &
           stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of constraints'

Q = zero
constraint_HO = zero
constraint_20 = zero
constraint_11 = zero
constraint_val = zero
lagrange_lambda0 = zero
lagrange_lambda1 = zero
Lchol = zero

!!! If opt_betalm > 0, calculates the corresponding value of Qlm
constraint_beta(:)  = constraint_read(5,:)
constraint_gamma(:) = constraint_read(7,:) * pi / 180.d0
                             
if ( opt_betalm > 0 ) then                             
  m = 2
  do i = 1, 4
    do j = 0, i 
      m = m + 1      
      p = 2 * (kdelta(constraint_switch(m),1) + kdelta(constraint_switch(m),3))
      do n = 1, min(constraint_switch(m),2)
        constraint_read(m,n) = constraint_read(m,n) / coeff_betalm(n+p,i)
      enddo
    enddo           
  enddo           
endif           
                             
if ( opt_betalm > 1 ) then                             
  !!! beta triaxial          
  p = 2 * (kdelta(constraint_switch(5),1) + kdelta(constraint_switch(5),3))
                             
  do n = 1, min(constraint_switch(5),2)
    constraint_read(5,n) = constraint_beta(n) * cos(constraint_gamma(n)) &
                           / coeff_betalm(n+p,2)
  enddo
                             
  !!! gamma triaxial         
  p = 2 * (kdelta(constraint_switch(7),1) + kdelta(constraint_switch(7),3))
                             
  do n = 1, min(constraint_switch(7),2)
    constraint_read(7,n) = constraint_beta(n) * sin(constraint_gamma(n)) & 
                           / (sqrt(2.d0) * coeff_betalm(n+p,2))
  enddo
endif

!!! Calculates the square of the constraints on the radius
p = 2 * (kdelta(constraint_switch(17),1) + kdelta(constraint_switch(17),3))
sig = sign(one,constraint_read(17,2))

constraint_read(17,1) = constraint_read(17,1)**2 / coeff_r2(1+p)
constraint_read(17,2) =  sig * constraint_read(17,2)**2 / coeff_r2(2+p)
                             
!!! If constraint_switch = 3, calculates the proton/neutron constraints
!!! conrresponding to the isoscalar/isovector constraints
do i = 1, constraint_types
  if ( constraint_switch(i) == 3 ) then 
    value_is = constraint_read(i,1)
    value_iv = constraint_read(i,2)

    constraint_read(i,1) = (value_is - value_iv) / 2
    constraint_read(i,2) = (value_is + value_iv) / 2
  endif 
enddo

!!! Collects all single-particle matrix elements from the different constraints
j = 0
constraint_id = 0 

do i = 1, constraint_types
  do m = 1, min(constraint_switch(i),2)
    j = j + 1
    k = (j - 1) * HOsp_dim2 
    if ( (i >= 21) .and. (i < 27) ) constraint_pair = constraint_pair + 1
 
    constraint_id(j) = i 
 
    select case (i)
      case (1)
        Q = partnumb_Z
      case (2)
        Q = partnumb_N
      case (3,4)
        n = i - 3
        if ( constraint_switch(i) == 1 ) then
          Q(:) = multipole_Q1m(:,1,n) + multipole_Q1m(:,2,n)
        else
          Q(:) = multipole_Q1m(:,m,n)
        endif
      case (5,6,7)
        n = i - 5
        if ( constraint_switch(i) == 1 ) then
          Q(:) = multipole_Q2m(:,1,n) + multipole_Q2m(:,2,n)
        else
          Q(:) = multipole_Q2m(:,m,n)
        endif
      case (8,9,10,11)
        n = i - 8
        if ( constraint_switch(i) == 1 ) then
          Q(:) = multipole_Q3m(:,1,n) + multipole_Q3m(:,2,n)
        else
          Q(:) = multipole_Q3m(:,m,n)
        endif
      case (12,13,14,15,16)
        n = i - 12
        if ( constraint_switch(i) == 1 ) then
          Q(:) = multipole_Q4m(:,1,n) + multipole_Q4m(:,2,n)
        else
          Q(:) = multipole_Q4m(:,m,n)
        endif
      case (17)
        if ( constraint_switch(i) == 1 ) then
          Q(:) = radius_r2(:,1) + radius_r2(:,2)
        else
          Q(:) = radius_r2(:,m)
        endif
      case (18)
        Q = angumome_Jx
      case (19)
        Q = angumome_Jy
      case (20)
        Q = angumome_Jz
      case (21)
        Q = pairs_T00_J10
      case (22)
        Q = pairs_T00_J1m1
      case (23)
        Q = pairs_T00_J1p1
      case (24)
        Q = pairs_T10_J00
      case (25)
        Q = pairs_T1m1_J00
      case (26)
        Q = pairs_T1p1_J00
      case (27)
        factor_delta = constraint_read(i,1)
    end select
   
    if ( i < 27 ) then
      constraint_HO(k+1:k+HOsp_dim2) = Q(1:HOsp_dim2)
      constraint_val(j) = constraint_read(i,m)
    endif
  enddo
enddo

deallocate(Q)

end subroutine set_constraints

!------------------------------------------------------------------------------!
! subroutine adjust_constraints                                                !
!                                                                              ! 
! Adjusts the constraints of the wavefunction until they satisfy the values    ! 
! of the constraints read in the input file. At the end of the process, the    ! 
! Bogoliubov matrices U0,V0 are updated                                        ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        Zi = gradient                                                         !
!------------------------------------------------------------------------------!
subroutine adjust_constraints(Zi,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim**2), intent(inout) :: Zi
integer, parameter :: iter_constr=100 
integer :: i, j, k, iter, m, ia, idig, info
integer(i64), parameter :: d1=1
integer(i64) :: ndim2, indi, indj, indk
real(r64) :: expvalp, expvaln
real(r64), dimension(constraint_dim) :: B, expecval, ipiv, work
real(r64), dimension(constraint_dim,constraint_dim) :: A
real(r64), dimension(((constraint_dim)*(constraint_dim+1)/2)*d1*ndim**2) :: A1

ndim2 = ndim**2

do iter = 1, iter_constr

  call cholesky_factorization(Zi,Lchol,ndim)
  call evolve_wavefunction(Zi,ndim)
  call calculate_densities_real(bogo_U1,bogo_V1,dens_rhoRR,dens_kappaRR,ndim)

  !!! Computes the expectation values of the operators
  do m = 1, constraint_dim
    if ( m <= (constraint_dim - constraint_pair) ) then
      call calculate_expectval_obo(dens_rhoRR,constraint_HO(1+(m-1)*ndim2), &
                                   expvalp,expvaln,ndim)
      expecval(m) =  expvalp + expvaln
    else 
      call calculate_expectval_pair(dens_kappaRR,constraint_HO(1+(m-1)*ndim2), &
                                    expecval(m),ndim)
    endif  
  enddo        

  !!! Updates the matrix elements of the constraints in the QP basis.
  !!! Then, exits the loop if all constraints are satisfied (up to tolerance)
  ia = 0
  do i = 1, constraint_dim
    indi = (i-1)*ndim2 + 1
    B(i) = constraint_val(i) - expecval(i)
    if ( abs(B(i)) > constraint_eps ) ia = ia + 1
    if ( i <= constraint_dim-constraint_pair ) then
      call set_operator_qpbasis('f20',bogo_U1,bogo_V1,constraint_HO(indi), & 
                                constraint_20(indi),ndim)
    else 
      call set_operator_qpbasis('g20',bogo_U1,bogo_V1,constraint_HO(indi), &
                                constraint_20(indi),ndim)
    endif 
  enddo

  if ( ia == 0 ) exit

  !!! Product of constraints C20^T * C20
  k = 0
  do i = 1, constraint_dim
    indi = (i-1)*ndim2 + 1
    do j = i, constraint_dim
      indj = (j-1)*ndim2 + 1
      k  = k + 1
      indk = (k-1)*ndim2 + 1
      call dgemm('t','n',ndim,ndim,ndim,one,constraint_20(indi),ndim, &
                 constraint_20(indj),ndim,zero,A1(indk),ndim)
    enddo
  enddo
        
  k = 0
  do i = 1, constraint_dim
    do j = i, constraint_dim
      k = k + 1
      indk = (k-1)*(ndim2) + 1
      A(j,i) = zero   
      do m = 1, ndim
        idig = (m-1)*ndim + m - 1
        A(j,i) = A(j,i) + A1(indk+idig)   
      enddo 
    enddo 
  enddo
        
  do i = 1, constraint_dim
    do j = i+1, constraint_dim
      A(i,j) = A(j,i)   
    enddo
  enddo

  call dsysv('U',constraint_dim,1,A,constraint_dim,ipiv,B,constraint_dim,work, &
             constraint_dim,info)
       
  !!! Update the gradient 
  k = 0
  do i = 1, ndim
    do j = 1, ndim
      k = k + 1
      do m = 1, constraint_dim
        Zi(k) = Zi(k) + B(m) * constraint_20(k+(m-1)*ndim2)
      enddo
    enddo
  enddo

enddo

!!! Update the Bogoliubov matrices
bogo_U0 = bogo_U1
bogo_V0 = bogo_V1

end subroutine adjust_constraints

!------------------------------------------------------------------------------!
! subroutine update_constraints_qpbasis                                        !
!                                                                              ! 
! Updates the matrix elements of the constraints in the quasiparticle basis    ! 
! of the new wave function.                                                    ! 
!------------------------------------------------------------------------------!
subroutine update_constraints_qpbasis

integer :: i
integer(i64) :: j

do i = 1, constraint_dim                 
  j = (i-1)*HOsp_dim2 + 1
  if ( i <= (constraint_dim - constraint_pair) ) then
    call set_operator_qpbasis('f20',bogo_U0,bogo_V0,constraint_HO(j), &
                              constraint_20(j),HOsp_dim)
  else
    call set_operator_qpbasis('g20',bogo_U0,bogo_V0,constraint_HO(j), &
                              constraint_20(j),HOsp_dim)
  endif
enddo

end subroutine update_constraints_qpbasis

!------------------------------------------------------------------------------!
! subroutine calculate_lagrange_multipliers                                    !
!                                                                              ! 
! Computes the Lagrange multipliers lambdas solving a system of equations.     ! 
!   A x = b                                                                    ! 
! where                                                                        ! 
!   A_ij = Tr(Q20(i) Q20^T(j))                                                 ! 
!   b_i  = Tr(H20 Q20^T(i))                                                    ! 
!                                                                              ! 
! This is taken from appendix A of Egido.1980.NuclPhysA.334.1                  !
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        H20 = 20 part of the Hamiltonian in the quasi-particle basis          ! 
!------------------------------------------------------------------------------!
subroutine calculate_lagrange_multipliers(H20,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim*ndim), intent(in) :: H20
integer :: i, j, l, m, indi, indj, indl, idig, info
real(r64), dimension(constraint_dim) :: ipiv, work, B
real(r64), dimension(constraint_dim,constraint_dim) :: A
real(r64), dimension(((constraint_dim+1)*constraint_dim/2)*ndim*ndim) :: A1
real(r64), dimension(constraint_dim*ndim*ndim) :: A2

!!! Exists the routine if a fully unconstrained calculation
if ( constraint_dim == 0 ) return

l = 0
do i = 1, constraint_dim
  indi = (i-1) * ndim * ndim + 1
  do j = i, constraint_dim
    indj = (j-1) * ndim * ndim + 1
    l = l + 1
    indl = (l-1) * ndim * ndim + 1
    call dgemm('n','t',ndim,ndim,ndim,one,constraint_20(indi),ndim, & 
               constraint_20(indj),ndim,zero,A1(indl),ndim)
  enddo
  call dgemm('n','t',ndim,ndim,ndim,one,H20(1),ndim,constraint_20(indi), &
             ndim,zero,A2(indi),ndim)
enddo

l = 0
do i = 1, constraint_dim
  indi = (i-1) * ndim * ndim + 1
  do j = i, constraint_dim
    l = l + 1
    indl = (l-1) * ndim * ndim + 1
    A(j,i) = zero
    do m = 1, ndim
      idig= (m-1) * ndim + m - 1
      A(j,i) = A(j,i) + A1(indl+idig)
    enddo
    A(i,j) = A(j,i)
  enddo
  
  B(i) = zero
  do m = 1, ndim
    idig = (m-1) * ndim + m - 1
    B(i) = B(i) + A2(indi+idig)
  enddo
enddo

call dsysv('U',constraint_dim,1,A,constraint_dim,ipiv,B,constraint_dim,work, &
           constraint_dim,info)

do i = 1, constraint_dim
  lagrange_lambda1(i) = B(i)
enddo

!!! BB: to avoid the singular system when constraining things that are exactly
!!! zero, I guess it should be necessary to remove the constraint for the     
!!! current iteration when it is too small.
if ( info /= 0 ) then
  print '("Warning: impossible to obtain the new lagrange multipliers. Using &
         &the one of the previous iteration instead.")' 
  do i = 1, constraint_dim
    lagrange_lambda1(i) = lagrange_lambda0(i)
  enddo
endif

!!! Update the constraints
do i = 1, constraint_dim
  lagrange_lambda0(i) = lagrange_lambda1(i)
enddo

end subroutine calculate_lagrange_multipliers

!------------------------------------------------------------------------------!
! subroutine cholesky_factorization                                            !
!                                                                              !
! Computes the cholesky factorization                                          !
!   LL^+ = 1 + Z^t Z^*   where Z is the Thouless matrix                        !  
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!         Z = (real) Thouless matrix                                           !
!                                                                              !
! Output: L = matrix L of the Cholesky factorization                           !
!------------------------------------------------------------------------------!
subroutine cholesky_factorization(Z,L,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: Z
real(r64), dimension(ndim,ndim), intent(out) :: L
integer :: i, j, k, m, info
real(r64), dimension(ndim,ndim) :: A1, A2
real(r64), dimension(ndim+1,ndim) :: A3

!!! 1 + Z^t Z^* 
call dgemm('t','n',ndim,ndim,ndim,one,Z,ndim,Z,ndim,zero,A1,ndim)

do i = 1, ndim
  do j = i+1, ndim
    A2(i,j) = A1(i,j)
    A2(j,i) = A2(i,j)
  enddo
  A2(i,i) = A1(i,i) + one 
enddo

! Switch to dpbtrf format
do i = 1, ndim
  do j = 1, i   
    k = 1 + i - j
    A3(k,j) = A2(i,j)
  enddo
enddo

!!! Cholesky                             
call dpbtrf('L',ndim,ndim,A3,ndim+1,info)

if ( info /= 0 ) then
  print*, 'Cholesky factorization has failed, info = ', info
  stop 
endif

! Switch back to our format
L = zero
do i = 1, ndim
  do j = 1, ndim+1-i
    k = i - 1 + j
    m = j
    L(k,m) = A3(i,j)
  enddo
enddo

end subroutine cholesky_factorization

!------------------------------------------------------------------------------!
! subroutine evolve_wavefunction                                               !
!                                                                              !
! Computes the new matrices U1,V1 for the next itration using the equations    !
! that can be found in NuclPhysA, 594, 70, 1995                                !
!   U1 = (U0 + V0 * Z) * Lchol^{-t}                                            !
!   V1 = (V0 + U0 * Z) * Lchol^{-t}                                            !
! where                                                                        !
!   Z is the gradient (for the constraint)                                     !
!   Lchol is the cholesky matrix (see "cholesky_factorization")                !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        Z = gradient                                                          !
!------------------------------------------------------------------------------!
subroutine evolve_wavefunction(Z,ndim)

integer, intent(in) :: ndim 
real(r64), dimension(ndim,ndim), intent(in) :: Z
integer :: i, j, info
real(r64), dimension(ndim,ndim) :: A1, A2

!!! Inverses Cholesky L matrix
call dtrtri('L','N',ndim,Lchol,ndim,info)
if ( info /= 0 ) then
  print*, 'Error in the inverse of L (Cholesky), info = ', info
  stop 
endif

!!! A1 = U0 + V0 * Z
!!! A2 = V0 + U0 * Z
call dgemm('n','n',ndim,ndim,ndim,one,bogo_V0,ndim,Z,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,bogo_U0,ndim,Z,ndim,zero,A2,ndim)        

do j = 1, ndim
  do i = 1, ndim
    A1(i,j) = A1(i,j) + bogo_U0(i,j)  
    A2(i,j) = A2(i,j) + bogo_V0(i,j)  
  enddo
enddo
        
!!! U1 = A1 * Lchol^{-t}   
!!! V1 = A2 * Lchol^{-t}
call dgemm('n','t',ndim,ndim,ndim,one,A1,ndim,Lchol,ndim,zero,bogo_U1,ndim)    
call dgemm('n','t',ndim,ndim,ndim,one,A2,ndim,Lchol,ndim,zero,bogo_V1,ndim)  
 
end subroutine evolve_wavefunction

!------------------------------------------------------------------------------!
! subroutine check_symmetries                                                  !
!                                                                              !
! Checks the symmetries of the seed wave function and if the constraint cons-  !
! erves them. Then simplifies the projection (if seed_symm = 0).               !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!------------------------------------------------------------------------------!
subroutine check_symmetries(ndim)

use Nucleus, only: valence_A
use Projection, only: proj_Mphip, proj_Mphin, proj_facpi, calculate_overlap, &
                      generate_rotation_parity, ROT

integer, intent(in) :: ndim
integer :: i, j, li, lj, hdim, nocc0, nemp0, nparity
real(r64) :: eps=1.0d-16, sum_pp, sum_pn, sum_np, sum_nn, sum_p, sum_pair, &
             value_Z, value_N, value_P, value_K, value_Kp, value_Kn, value_K2, &
             ovac0
real(r64), dimension(:), allocatable :: voveru0
complex(r64) :: pari
complex(r64), dimension(ndim,ndim) :: ROTP
logical :: is_good_N0, is_good_Z0, is_separate_NZ0
character(len=:), allocatable :: action_N, action_Z, action_NZ
character(3) :: ch_good_N=" no", ch_good_Z=" no", ch_good_P=" no", &
                ch_separate_NZ=" no", ch_good_K=" no"
character(len=*), parameter :: format1 = "(1a12,4x,1a3,1x,1f10.5,3x,1a)", &
                               format2 = "(1a12,4x,1a3,14x,1a)"

hdim = ndim/2
is_separate_NZ0 = .false.
is_good_Z0 = .false.
is_good_N0 = .false.

!!! Checks if the number parity of the state is compatible with the projection
if ( max(proj_Mphip,proj_Mphin) > 1 ) then
  nparity = (-1)**(int(valence_A))
  if ( nparity /= bogo_nparity ) then
    !cmpi if ( paral_myrank == 0 ) then        
    print '(/,1a)', "Critical warning: the number parity of the seed wave &
                  & function is not consistent with the particle number in &
                  & the projection operator."
    !cmpi endif
    stop
  endif
endif

!!! Calculates N and Z
call calculate_densities_real(bogo_U0,bogo_V0,dens_rhoRR,dens_kappaRR,ndim)
call calculate_expectval_obo(dens_rhoRR,partnumb_A,value_Z,value_N,ndim)

!!! Symmetries for N and Z
sum_pp = zero
sum_pn = zero
sum_np = zero
sum_nn = zero

do j = 1, hdim
  do i = 1, hdim
    sum_pp = sum_pp + abs(dens_kappaRR(i,j))
    sum_pn = sum_pn + abs(dens_kappaRR(i,j+hdim))
    sum_np = sum_np + abs(dens_kappaRR(i+hdim,j))
    sum_nn = sum_nn + abs(dens_kappaRR(i+hdim,j+hdim))
  enddo
enddo

if ( (sum_np+sum_pn) < eps ) then
  is_separate_NZ  = .true.
  is_separate_NZ0 = .true.
  ch_separate_NZ = "yes"

  if ( (sum_np+sum_pn+sum_pp) < eps ) then 
    is_good_Z  = .true.
    is_good_Z0 = .true.
    ch_good_Z = "yes"
  endif

  if ( (sum_np+sum_pn+sum_nn) < eps ) then
    is_good_N  = .true.
    is_good_N0 = .true.    
    ch_good_N = "yes"
  endif
endif

!!! Calculates parity  
call construct_canonical_basis(bogo_U0,bogo_V0,bogo_zU0c,bogo_zV0c,bogo_zD0, &
                               ovac0,nocc0,nemp0,ndim)
if ( ndim-nocc0-nemp0 == 0 ) ovac0 = 0.0d0

allocate(voveru0(ndim-nocc0))
voveru0 = zero

j = 1
do i = 1+nemp0, ndim-nocc0
  if ( (-1)**j == -1 ) then 
    voveru0(i) = real( bogo_zV0c(i,i+1) / bogo_zU0c(i,i) )
  else
    voveru0(i) = real( bogo_zV0c(i,i-1) / bogo_zU0c(i,i) )
  endif
  j = j + 1
enddo

call generate_rotation_parity(1,ROT,ROTP,ndim)
call calculate_overlap(nocc0,nemp0,ovac0,voveru0,bogo_zD0,nocc0,nemp0,ovac0, &
                       voveru0,bogo_zD0,ROTP,pari,ndim)

value_P = real(pari)

!!! Symmetries for P
sum_p = zero

do j = 1, ndim
  lj = HOsp_l(j)
  do i = 1, ndim
    li = HOsp_l(i)
    if ( (-1)**li /= (-1)**lj ) sum_p = sum_p + abs(dens_rhoRR(i,j)) & 
                                              + abs(dens_kappaRR(i,j))
  enddo
enddo

if ( sum_p < eps ) then
  is_good_P = .true.
  ch_good_P = "yes"
endif

!!! Calculates <Jz> and <Jz^2>
call calculate_expectval_obo(dens_rhoRR,angumome_Jz,value_Kp,value_Kn,ndim)
call calculate_expectval_obos(dens_rhoRR,dens_kappaRR,angumome_Jz,angumome_Jz2, &
                              value_K2,ndim)
value_K = value_Kp + value_Kn

!!! Symmetries for K
if ( abs(value_K**2 -value_K2) < sqrt(eps) ) then
  is_good_K = .true.
  ch_good_K = "yes"
endif

!!! Takes into account the constraints that could break the symmetries
sum_pair = 0.0d0
do i = 13, 16
  sum_pair = sum_pair + abs(constraint_read(i,1))
enddo
if ( sum_pair > eps ) is_separate_NZ = .false.

sum_pair = 0.0d0
do i = 13, 17
  sum_pair = sum_pair + abs(constraint_read(i,1))
enddo
if ( sum_pair > eps ) is_good_Z = .false.

sum_pair = 0.0d0
do i = 13, 18
  if ( i == 17 ) cycle
  sum_pair = sum_pair + abs(constraint_read(i,1))
enddo
if ( sum_pair > eps ) is_good_N = .false.

!!! Simplification separate N/Z
if ( is_separate_NZ ) then
  if ( seed_symm == 0 ) then
    proj_facpi = pi
    action_NZ = "     integral [1,pi]     "
  else
    action_NZ = "     none (seed_symm)    "
  endif
else
  if ( is_separate_NZ0 ) then
    action_NZ = "    none (constraints)   "
  else
    action_NZ = " "
  endif
endif

!!! Simplification good Z      
if ( is_good_Z ) then
  if ( seed_symm == 0 ) then
    proj_Mphip = 1
    action_Z = "no projection (Mphip = 1)"
  else
    action_Z = "     none (seed_symm)    "
  endif
else
  if ( is_good_Z0 ) then
    action_Z = "    none (constraints)   "
  else
    action_Z = " "
  endif
endif

!!! Simplification good N
if ( is_good_N ) then
  if ( seed_symm == 0 ) then
    proj_Mphin = 1
    action_N = "no projection (Mphin = 1)"
  else
    action_N = "     none (seed_symm)    "
  endif
else
  if ( is_good_N0 ) then
    action_Z = "    none (constraints)   "
  else
    action_N = " "
  endif
endif

!!! Prints the summary
!cmpi if ( paral_myrank == 0 ) then        
print '(/,2x,"Symmetry",7x,"?",6x,"Mean",14x,"Action",/,58("-"))'
print format2, 'Separate N/Z',ch_separate_NZ, action_NZ
print format1, 'Good Z      ',ch_good_Z, value_Z, action_Z
print format1, 'Good N      ',ch_good_N, value_N, action_N
print format1, 'Good P      ',ch_good_P, value_P
print format1, 'Good Jz     ',ch_good_K, value_K
!cmpi endif

deallocate(voveru0)

end subroutine check_symmetries

END MODULE Constraints
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
