!==============================================================================!
! MODULE Operators                                                             !
!                                                                              !
! This module contains the variables and routines related to the general man-  !
! ipulation of operators.                                                      !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_operator_qpbasis                                            !
! - subroutine calculate_expectval_obo                                         !
! - subroutine calculate_expectval_obo_cplx                                    !
! - subroutine calculate_expectval_obos                                        !
! - subroutine calculate_expectval_obos_cplx                                   !
! - subroutine calculate_expectval_pair                                        !
!==============================================================================!
MODULE Operators 

use ParticleNumber 
use AngularMomentum
use Multipoles      
use Pairs
use Radius

implicit none
public 

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_operator_qpbasis                                              !
!                                                                              !
! Computes the matrix elements of operator A in the quasi-particle basis       !
! according to the equations (E.12) and (E.13) of the book by Ring and Schuck  !
! (ISBN: 978-3-540-21206-5).                                                   !
!                                                                              !
! Considering an operator                                                      !
!    A = \sum_ij f_ij c^+_i c_j + 1/2 ( g_ij c^+_i c^_j + hc )                 !
! the matrix elements are written                                              !
!    A11 = U^T f U - V^T f^T V + U^T g V - V^T g U                             !
!    A20 = U^T f V - V^T f^T U + U^T g U - V^T g V                             !
! where U, V, f and g have been taken real.                                    !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        A = matrix of the one-body operator in the HO basis                   !
!        U,V = Bogoliubov matrices                                             !
!        part = part to condiser: "f11", "f20", "g11", "g20"                   !
!                                                                              !
! Output: B = matrix of the one-body operator in the QP basis                  !
!------------------------------------------------------------------------------!
subroutine set_operator_qpbasis(part,U,V,A,B,ndim)

integer, intent(in) :: ndim
character(3), intent(in) :: part
real(r64), dimension(ndim,ndim), intent(in) :: U, V, A  
real(r64), dimension(ndim,ndim), intent(out) :: B             
integer :: i, j
real(r64) :: factor
real(r64), dimension(ndim,ndim) :: A1, A2, A3

if ( part == 'f11' ) then 
  factor = one
  call dgemm('t','n',ndim,ndim,ndim,one, U,ndim,A,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,U,ndim,zero,A2,ndim)
  call dgemm('t','t',ndim,ndim,ndim,one, V,ndim,A,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,V,ndim,zero,A3,ndim)
elseif ( part == 'g11' ) then 
  factor = one
  call dgemm('t','n',ndim,ndim,ndim,one, U,ndim,A,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,V,ndim,zero,A2,ndim)
  call dgemm('t','n',ndim,ndim,ndim,one, V,ndim,A,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,U,ndim,zero,A3,ndim)
elseif ( part == 'f20' ) then 
  factor = -one
  call dgemm('t','n',ndim,ndim,ndim,one, U,ndim,A,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,V,ndim,zero,A2,ndim)
  call dgemm('t','t',ndim,ndim,ndim,one, V,ndim,A,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,U,ndim,zero,A3,ndim)
elseif ( part == 'g20' ) then 
  factor = -one
  call dgemm('t','n',ndim,ndim,ndim,one, U,ndim,A,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,U,ndim,zero,A2,ndim)
  call dgemm('t','n',ndim,ndim,ndim,one, V,ndim,A,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,V,ndim,zero,A3,ndim)
endif

do j = 1, ndim
  do i = 1, ndim
    A1(i,j) = A2(i,j) - A3(i,j)
  enddo
enddo

do j = 1, ndim
  do i = 1, ndim
    B(i,j) = 0.5d0 * (A1(i,j) + factor * A1(j,i))
  enddo
enddo

end subroutine set_operator_qpbasis

!------------------------------------------------------------------------------!
! subroutine calculate_expectval_obo                                           !
!                                                                              !
! Computes the expectation value for a one-body operator B of the form         !
!    B = \sum_ij b_ij c^+_i c_j                                                !
! which is simply                                                              !
!    < B > = Tr(b rhoRR) = Tr(b rhoRR_p) + Tr(b rhoRR_n)                       !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        B = matrix of the one-body operator in the HO basis                   !
!        rhoRR = density                                                       !
!                                                                              !
! Output: expval_p = expectation value for the proton  sp states               !
!         expval_n =      "        "    "   "  neutron  "   "                  !
!------------------------------------------------------------------------------!
subroutine calculate_expectval_obo(rhoRR,B,expval_p,expval_n,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: rhoRR, B
real(r64), intent(out) :: expval_p, expval_n
integer :: i
real(r64), dimension(ndim,ndim) :: A1

call dgemm('n','n',ndim,ndim,ndim,one,B,ndim,rhoRR,ndim,zero,A1,ndim)

expval_p = zero
expval_n = zero

do i = 1, ndim/2
  expval_p = expval_p + A1(i,i)
  expval_n = expval_n + A1(i+ndim/2,i+ndim/2)
enddo

end subroutine calculate_expectval_obo

!------------------------------------------------------------------------------!
! subroutine calculate_expectval_obo_cplx                                      !
!                                                                              !
! Same as calculate_expectval_obo but for complex quantities.                  !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        B = matrix of the one-body operator in the HO basis                   !
!        rhoRR = density                                                       !
!                                                                              !
! Output: expval_p = expectation value for the proton  sp states               !
!         expval_n =      "        "    "   "  neutron  "   "                  !
!------------------------------------------------------------------------------!
subroutine calculate_expectval_obo_cplx(rhoLR,B,expval_p,expval_n,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, B
complex(r64), intent(out) :: expval_p, expval_n
integer :: i
complex(r64), dimension(ndim,ndim) :: A1

call zgemm('n','n',ndim,ndim,ndim,zone,B,ndim,rhoLR,ndim,zzero,A1,ndim)

expval_p = zzero
expval_n = zzero

do i = 1, ndim/2
  expval_p = expval_p + A1(i,i)
  expval_n = expval_n + A1(i+ndim/2,i+ndim/2)
enddo

end subroutine calculate_expectval_obo_cplx

!------------------------------------------------------------------------------!
! subroutine calculate_expectval_obos                                          !
!                                                                              !
! Computes the expect value for a two-body operator that is obtained as the    !
! square of a one-body opeartor, for example Jx^2, in the diagonal case with   !
! real densities.                                                              !
! Considering an operator                                                      !
!    O = \sum_ij o_ij c^+_c a_j                                                !
! the expectation value can be written                                         !
!    < O^2 > = Tr(o rhoRR)^2 + Tr(o^2 rhoRR) - Tr([O rhoRR]^2)                 !
!              + Tr(o^T kappaRR^T o kappaRR)                                   !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        rhoRR,kappaRR = densities                                             !
!        O,O2 = matrix elements of O and O^2 in the HO basis                   !
!                                                                              !
! Output: expval = expectation value of O^2                                    !
!------------------------------------------------------------------------------!
subroutine calculate_expectval_obos(rhoRR,kappaRR,O,O2,expval,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: rhoRR, kappaRR, O, O2
real(r64), intent(out) :: expval
integer :: i
real(r64) :: tr1, tr2, tr5, tr6
real(r64), dimension(ndim,ndim) :: A1, A2, A3, A4, A5, A6

call dgemm('n','n',ndim,ndim,ndim,one, O,ndim,  rhoRR,ndim,zero,A1,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,O2,ndim,  rhoRR,ndim,zero,A2,ndim)
call dgemm('t','t',ndim,ndim,ndim,one, O,ndim,kappaRR,ndim,zero,A3,ndim)   
call dgemm('n','n',ndim,ndim,ndim,one, O,ndim,kappaRR,ndim,zero,A4,ndim)   
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,     A1,ndim,zero,A5,ndim)   
call dgemm('n','n',ndim,ndim,ndim,one,A3,ndim,     A4,ndim,zero,A6,ndim)

tr1 = zero
tr2 = zero
tr5 = zero
tr6 = zero

do i = 1, ndim
  tr1 = tr1 + A1(i,i)
  tr2 = tr2 + A2(i,i)
  tr5 = tr5 + A5(i,i)
  tr6 = tr6 + A6(i,i)
enddo

expval = tr1**2 + tr2 - tr5 + tr6

end subroutine calculate_expectval_obos

!------------------------------------------------------------------------------!
! subroutine calculate_expectval_obos_cplx                                     !
!                                                                              !
! Same as calcualte_expectval_obos but for complex quantities.                 !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        rhoRR,kappaRR = densities                                             !
!        O,O2 = matrix elements of O and O^2 in the HO basis                   !
!                                                                              !
! Output: expval = expectation value of O^2                                    !
!------------------------------------------------------------------------------!
subroutine calculate_expectval_obos_cplx(rhoLR,kappaLR,kappaRL,O,O2,expval,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL, O, O2
complex(r64), intent(out) :: expval
integer :: i
complex(r64) :: tr1, tr2, tr5, tr6
complex(r64), dimension(ndim,ndim) :: A1, A2, A3, A4, A5, A6

call zgemm('n','n',ndim,ndim,ndim,zone, O,ndim,  rhoLR,ndim,zzero,A1,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,O2,ndim,  rhoLR,ndim,zzero,A2,ndim)
call zgemm('t','t',ndim,ndim,ndim,zone, O,ndim,kappaRL,ndim,zzero,A3,ndim)   
call zgemm('n','n',ndim,ndim,ndim,zone, O,ndim,kappaLR,ndim,zzero,A4,ndim)   
call zgemm('n','n',ndim,ndim,ndim,zone,A1,ndim,     A1,ndim,zzero,A5,ndim)   
call zgemm('n','n',ndim,ndim,ndim,zone,A3,ndim,     A4,ndim,zzero,A6,ndim)

tr1 = zzero
tr2 = zzero
tr5 = zzero
tr6 = zzero

do i = 1, ndim
  tr1 = tr1 + A1(i,i)
  tr2 = tr2 + A2(i,i)
  tr5 = tr5 + A5(i,i)
  tr6 = tr6 + A6(i,i)
enddo

expval = tr1**2 + tr2 - tr5 + tr6

end subroutine calculate_expectval_obos_cplx

!------------------------------------------------------------------------------!
! subroutine calculate_expectval_pair                                          !
!                                                                              ! 
! Computes the expectation value for a pair operator of the form               ! 
!    B = \sum_ij b_ij c^+_i c^+_j                                              !
! which is simply                                                              !
!    < B > = Tr(b kappaRR^T)                                                   !
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        kappaRR = density                                                     !
!        B = matrix elements in the HO basis                                   !
!                                                                              !
! Output: expval = expectation value of B                                      !
!------------------------------------------------------------------------------!
subroutine calculate_expectval_pair(kappaRR,B,expval,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: B, kappaRR
real(r64), intent(out) :: expval
integer :: i              
real(r64), dimension(ndim,ndim) :: A

call dgemm('n','t',ndim,ndim,ndim,one,B,ndim,kappaRR,ndim,zero,A,ndim)

expval = zero

do i = 1, ndim
  expval = expval + A(i,i)
enddo

end subroutine calculate_expectval_pair

END MODULE Operators 
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
