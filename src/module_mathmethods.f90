!==============================================================================!
! MODULE MathMethods                                                           !
!                                                                              !
! This module contains the variables and routines related to mathematical fun- !
! ctions and numerical methods.                                                !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine GaussLaguerre                                                   !
! - subroutine ClebschGordan                                                   !
! - function factorial                                                         !
! - function dfactorial                                                        !
! - function kdelta                                                            !
! - function laguerre                                                          !
! - function lapack_sel                                                        !
!==============================================================================!
MODULE MathMethods 

use Constants

implicit none
public

CONTAINS 

!------------------------------------------------------------------------------!
! function factorial                                                           !
!                                                                              !
! Computes the factorial: n! = n * (n-1) * ... * 1                             !
!------------------------------------------------------------------------------!
recursive function factorial(n) result(facto)

integer, intent(in) :: n
real(r64) :: facto 

if ( n <= 0 ) then 
  facto = one  
else
  facto = n * factorial(n-1)
endif

end function

!------------------------------------------------------------------------------!
! function dfactorial                                                          !
!                                                                              !
! Computes the double factorial: n!! = n * (n-2) * ... * 1                     !
!------------------------------------------------------------------------------!
recursive function dfactorial(n) result(facto)

integer, intent(in) :: n
real(r64) :: facto 

if ( n <= 0 ) then 
  facto = one  
else
  facto = n * dfactorial(n-2)
endif

end function

!------------------------------------------------------------------------------!
! function kdelta                                                              !
!                                                                              !
! Computes the Kronecker delta: \delta_ij = 1 if i = j                         !
!                                         = 0 otherwise                        !
!------------------------------------------------------------------------------!
function kdelta(i,j) result(delta)

integer, intent(in) :: i, j
integer :: delta

if ( i == j ) then 
  delta = 1    
else
  delta = 0    
endif

end function kdelta

!------------------------------------------------------------------------------!
! function laguerre                                                            !
!                                                                              !
! Calculates the value of the generalized Laguerre polynomial L^a_n(x)         !
! using the reccurence relation:                                               !
! L^a_n(x) = ( (2n-1+a-x)*L^a_n-1(x) - (n-1+a)*L^a_n-2(x) ) / n                !
!------------------------------------------------------------------------------!
function laguerre(x,n,a) result(l1)

integer, intent(in) :: n
real(r64), intent(in) :: x, a
integer :: i
real(r64) :: l1, l2, l3 

l1 = 1.d0
l2 = 0.d0

do i = 1, n  
  l3 = l2
  l2 = l1
  l1 = ( (2*i-1+a-x)*l2 - (i-1+a)*l3 ) / i
enddo

end function laguerre

!------------------------------------------------------------------------------!
! subroutine GaussLaguerre                                                     !
!                                                                              !
! Computes the abscissas and weight for the Gauss-Laguerre quadrature.         !
! The algorithm is taken from the book: Numerical Recipe in Fortran 77         !
! (ISBN: 978-0521430647).                                                      !
!------------------------------------------------------------------------------!
subroutine GaussLaguerre(x,w,n,alf)

integer, intent(in) :: n
real(r64), intent(in) :: alf
real(r64), intent(out) :: w(n), x(n)
integer, parameter :: maxit=20
real(r64), parameter :: eps=1.0d-16
integer :: i, its, j
real(r64) :: ai, p1, p2, p3, pp, z, z1

do i = 1, n
  if ( i == 1 ) then
    z = (1.0 + alf) * (3.0 + 0.92*alf) / (1.0 + 2.4*n + 1.8*alf)
  elseif ( i == 2 ) then
    z = z + (15.0 + 6.25*alf) / (1.0 + 0.9*alf + 2.5*n)
  else
    ai = i - 2
    z = z + ( ((1.0 + 2.55*ai) / (1.9*ai)) + (1.26*ai*alf / (1.0 + 3.5*ai)) ) &
            * (z - x(i-2)) / (1.0 + 0.3*alf)
  endif

  do its = 1, maxit
    p1 = 1.d0
    p2 = 0.d0
    do j = 1, n
      p3 = p2
      p2 = p1
      p1 = ( (2*j-1+alf-z)*p2 - (j-1+alf)*p3 ) / j
    enddo
    pp = (n*p1 - (n+alf)*p2) / z
    z1 = z
    z = z1 - p1/pp
    if ( abs(z - z1) <= eps) exit  
  enddo

  x(i) = z
  w(i) = -exp(log_gamma(alf+n) - log_gamma(n+0.0d0)) / (pp*n*p2)
enddo

end subroutine GaussLaguerre

!------------------------------------------------------------------------------!
! subroutine ClebschGordan                                                     !
!                                                                              !
! Computes the ClebschGordan for the group SU(2). The algorithm used was taken !
! from technical notes from NASA written by W. F. Ford and R. C. Bruley.       !
! Ref: NASA TN D-6173                                                          !
!                                                                              !
! (j1,m1,j2,m2|j3,m3) = c * g                                                  !
! with                                                                         !
! c = D(j1j2j3) * [(j1-m1)!(j2+m2)!(j1+m1)!(j2-m2)!(j3+m3)!(j3-m3)!]**1/2      !
! g = sqrt(2*j3+1) sum_l (-1)**l [(j1+j2-j3-l)!(j1-m1-l)!(j2+m2-l)!            !
!                                 (j3-j2+m1+l)!(j3-j1-m1+l)!l!]**-1            !
! D(j1j2j3) = [(j1+j2-j3)!(j2+j3-j1)!(j3+j1-j2)!]**1/2 / [(j1+j2+j3+1)!]**1/2  !
!                                                                              !
! Be aware that all input values of j and m are supposed to be twice their     !
! values (such that we can treat easily half-integer angular momenta).         !
!------------------------------------------------------------------------------!
subroutine ClebschGordan (j1,j2,j3,m1,m2,m3,cg)

integer, intent(in) :: j1, j2, j3, m1, m2, m3
real(r64), intent(out) :: cg 
integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, k1, k2, k3, k4, k5, k6, &
           l, l1, l2
real(r64) :: c, g, p, q, h, hl, hlm1

cg = zero

!!! Computes the factor c
n1 = 1 + (j1 + j2 - j3)/2
n2 = 1 + (j2 + j3 - j1)/2
n3 = 1 + (j3 + j1 - j2)/2
n4 = 1 + (j1 - m1)/2
n5 = 1 + (j2 + m2)/2
n6 = 1 + (j1 + m1)/2
n7 = 1 + (j2 - m2)/2
n8 = 1 + (j3 + m3)/2
n9 = 1 + (j3 - m3)/2
n10= n1 + n2 + n3 - 1

if ( (min(n1,n2,n3,n4,n5,n6,n7,n8,n9) < 1) .or. (m1+m2 /= m3) ) return

p =  log_gamma(n1+zero) + log_gamma(n2+zero) + log_gamma(n3+zero) &
   + log_gamma(n4+zero) + log_gamma(n5+zero) + log_gamma(n6+zero) &
   + log_gamma(n7+zero) + log_gamma(n8+zero) + log_gamma(n9+zero) &
   - log_gamma(n10+zero)

c = exp(0.5d0*p)

!!! Computes the factor g
k1 = n1
k2 = n4
k3 = n5
k4 = n4 - n3
k5 = n5 - n2

l1 = max(0,k4,k5)
l2 = min(k1,k2,k3)

h  = one 
hl = one

do l = l1+1, l2
  hlm1 = hl
  hl = hlm1 * (l - k1) * (l - k2) * (l - k3) / ((l - k4) * (l - k5) * l)
  h = h + hl
enddo

k1 = k1 - l1
k2 = k2 - l1
k3 = k3 - l1
k4 = l1 + 1 - k4 
k5 = l1 + 1 - k5 
k6 = l1 + 1 

q =  log_gamma(k1+zero) + log_gamma(k2+zero) + log_gamma(k3+zero) &
   + log_gamma(k4+zero) + log_gamma(k5+zero) + log_gamma(k6+zero)

g = sqrt(j3 + one) * (-1)**l1 * exp(-q) * h

!!! Computes the final value combining the two parts.
cg = c * g

end subroutine ClebschGordan

!------------------------------------------------------------------------------!
! function lapack_sel                                                          !
!                                                                              !
! Dummy logical function "select" for LAPACK (used for Schur decomposition).   !
! See LAPACK documentation for dgees.                                          !
!------------------------------------------------------------------------------!
function lapack_sel(wr,wi) result(selec)

logical :: selec
real(r64), intent(in) :: wr, wi

!!! Set to false because the function is only used as dummy argument
selec=.false.

end function lapack_sel

END MODULE MathMethods 
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
