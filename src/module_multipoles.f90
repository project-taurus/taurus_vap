!==============================================================================!
! MODULE Multipoles                                                            !
!                                                                              !
! This module contains the variables and routines related to deformation mult- !
! ipole operators (up to Q44 for now).                                         !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_multipoles                                                  !
!==============================================================================!
MODULE Multipoles

use Basis

implicit none 

!!! Qlm(i,j,k) with i=number matrix elem, j=isospin (1=p,2=n), k=0...m
real(r64), dimension(:,:,:), allocatable :: multipole_Q1m, & ! Dipoles
                                            multipole_Q2m, & ! Quadrupoles
                                            multipole_Q3m, & ! Octupoles  
                                            multipole_Q4m    ! Hexadecapoles

real(r64) :: coeff_betalm(4,4) ! Factor to convert from beta_lm to Qlm

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_multipoles                                                    !
!                                                                              !
! Defines the HO matrix elements of the multipoles operators. The following    !
! formula is used:                                                             !
! <a|Q_lambda,m|b> = b^lambda * R^lambda_ab * sqrt((2lb+1)*(2lambda+1)/(4pi *  !
!                    (2la+1))) * CG(lb,0,lambda,0,la,0) *                      !
!                    sum_ms,mla,mlb CG(la,mla,1/2,ms,ja,ma) *                  !
!                       CG(lb,mlb,1/2,ms,jb,mb) * CG(lb,mlb,lambda,m,la,mla)   !
! where b is the oscillator length, R is the radial part, and CGs are Clebsch- !
! Gordan coefficients.                                                         !
!                                                                              !
! To obtain an hermitian operator, we constrain  the average                   !
! 0.5 * (Qlm + Qlm^\dagger) = 0.5 * (Qlm + (-)^m Ql-m)                         !
!                                                                              !
! The operators are vectorized (for the calculation of constraints).           !
!------------------------------------------------------------------------------!
subroutine set_multipoles    

integer :: i, j, incr, ia, ib, la, lb, ja, jb, mja, mla, mlb, mjb, mta, mtb, &
           ms, lambda, mu, mt, mq, ialloc=0
real(r64) :: sumcb, fac, cb1, cb2, cb3, coeff_1, coeff_2, numbpart

allocate( multipole_Q1m(HOsp_dim2,2,0:1), multipole_Q2m(HOsp_dim2,2,0:2), &
          multipole_Q3m(HOsp_dim2,2,0:3), multipole_Q4m(HOsp_dim2,2,0:4), &
          source=zero, stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of multipole operators'

!!! Factors to convert from beta_lm to Qlm
do j = 1, 4
  do i = 1, 4
    if ( i == 1 ) then
      numbpart = nucleus_Z
    elseif ( i == 2 ) then
      numbpart = nucleus_N
    else 
      numbpart = nucleus_A
    endif
    
    coeff_1 = (4.d0 * pi) / (3.d0 * numbpart)
    coeff_2 = 1.d0 / (radius_r0 * numbpart**(1.d0/3.d0))

    coeff_betalm(i,j) = coeff_1 * coeff_2**j
  enddo
enddo

!!! Matrix elements                             
incr = 0

do ia = 1, HOsp_dim
  la = HOsp_l(ia) * 2
  ja = HOsp_2j(ia)
  mja = HOsp_2mj(ia)
  mta = HOsp_2mt(ia)
  do ib = 1, HOsp_dim
    lb = HOsp_l(ib) * 2
    jb = HOsp_2j(ib)
    mjb = HOsp_2mj(ib)
    mtb = HOsp_2mt(ib)
    incr = incr + 1
    if ( mta /= mtb ) cycle

    do lambda = 2, 8, 2
      if ( mod(la+lb+lambda,4) /= 0 ) cycle
      call ClebschGordan (lb,lambda,la,0,0,0,cb1)
      fac = cb1 * radial_integral_even(ia,ib,lambda/2) * &
            sqrt( ((lb+1)*(lambda+1)) / (4*pi*(la+1)) )
      do mu = -lambda, lambda, 2
        sumcb = zzero
        do ms = -1, 1, 2    
          do mla = -la, la, 2    
            do mlb = -lb, lb, 2    
              call ClebschGordan (la,1,ja,mla,ms,mja,cb1)
              call ClebschGordan (lb,1,jb,mlb,ms,mjb,cb2)
              call ClebschGordan (lb,lambda,la,mlb,mu,mla,cb3)
              sumcb = sumcb + cb1 * cb2 * cb3
            enddo
          enddo
        enddo

        sumcb = sumcb * fac
        if ( mu < 0 ) sumcb = (-1)**(mu/2) * sumcb
        if ( mu /= 0 ) sumcb = sumcb * 0.5d0
       
        mt = (mta + 3)/2
        mq = abs(mu/2) 
       
        if ( lambda == 2 ) then
          multipole_Q1m(incr,mt,mq) = multipole_Q1m(incr,mt,mq) + sumcb 
        elseif ( lambda == 4 ) then                        
          multipole_Q2m(incr,mt,mq) = multipole_Q2m(incr,mt,mq) + sumcb 
        elseif ( lambda == 6 ) then                        
          multipole_Q3m(incr,mt,mq) = multipole_Q3m(incr,mt,mq) + sumcb 
        else                                               
          multipole_Q4m(incr,mt,mq) = multipole_Q4m(incr,mt,mq) + sumcb 
        endif
      enddo
    enddo 

  enddo
enddo

end subroutine set_multipoles 

END MODULE Multipoles
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
