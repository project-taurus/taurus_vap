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

real(r64), dimension(:), allocatable :: multipole_Q10, & ! Dipoles    
                                        multipole_Q11, &
                                        multipole_Q20, & ! Quadrupoles
                                        multipole_Q21, &
                                        multipole_Q22, & 
                                        multipole_Q30, & ! Octupoles
                                        multipole_Q31, &
                                        multipole_Q32, &
                                        multipole_Q33, &
                                        multipole_Q40, & ! Hexadecapoles
                                        multipole_Q41, &
                                        multipole_Q42, &
                                        multipole_Q43, &
                                        multipole_Q44

real(r64) :: coeff_betalm(4) ! Factor to convert from beta_lm to Qlm

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

integer :: i, incr, ia, ib, la, lb, ja, jb, mja, mla, mlb, mjb, mta, mtb, ms, &
           lambda, mu, ialloc=0
real(r64) :: sumcb, fac, cb1, cb2, cb3, coeff_1, coeff_2

allocate( multipole_Q10(HOsp_dim2), multipole_Q11(HOsp_dim2), &
          multipole_Q20(HOsp_dim2), multipole_Q21(HOsp_dim2), &
          multipole_Q22(HOsp_dim2), multipole_Q30(HOsp_dim2), &
          multipole_Q31(HOsp_dim2), multipole_Q32(HOsp_dim2), &
          multipole_Q33(HOsp_dim2), multipole_Q40(HOsp_dim2), &
          multipole_Q41(HOsp_dim2), multipole_Q42(HOsp_dim2), &
          multipole_Q43(HOsp_dim2), multipole_Q44(HOsp_dim2), &
          source=zero, stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of multipole operators'

!!! Factors to convert from beta_lm to Qlm
coeff_1 = (4.d0 * pi) / (3.d0 * nucleus_A)
coeff_2 = 1.d0 / (radius_r0 * nucleus_A**(1.d0/3.d0))

do i = 1, 4
  coeff_betalm(i) = coeff_1 * coeff_2**i
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
      fac = cb1 * radial_even(ia,ib,lambda/2) * &
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
        if ( mu /= 0) sumcb = sumcb * 0.5d0

        if ( lambda == 2 ) then
          if ( mu == -2 ) multipole_Q11(incr) = multipole_Q11(incr) - sumcb 
          if ( mu ==  0 ) multipole_Q10(incr) = multipole_Q10(incr) + sumcb 
          if ( mu ==  2 ) multipole_Q11(incr) = multipole_Q11(incr) + sumcb 
        elseif ( lambda == 4 ) then
          if ( mu == -4 ) multipole_Q22(incr) = multipole_Q22(incr) + sumcb 
          if ( mu == -2 ) multipole_Q21(incr) = multipole_Q21(incr) - sumcb 
          if ( mu ==  0 ) multipole_Q20(incr) = multipole_Q20(incr) + sumcb 
          if ( mu ==  2 ) multipole_Q21(incr) = multipole_Q21(incr) + sumcb 
          if ( mu ==  4 ) multipole_Q22(incr) = multipole_Q22(incr) + sumcb 
        elseif ( lambda == 6 ) then
          if ( mu == -6 )  multipole_Q33(incr) = multipole_Q33(incr) - sumcb
          if ( mu == -4 )  multipole_Q32(incr) = multipole_Q32(incr) + sumcb
          if ( mu == -2 )  multipole_Q31(incr) = multipole_Q31(incr) - sumcb
          if ( mu ==  0 )  multipole_Q30(incr) = multipole_Q30(incr) + sumcb
          if ( mu ==  2 )  multipole_Q31(incr) = multipole_Q31(incr) + sumcb
          if ( mu ==  4 )  multipole_Q32(incr) = multipole_Q32(incr) + sumcb
          if ( mu ==  6 )  multipole_Q33(incr) = multipole_Q33(incr) + sumcb
        else 
          if ( mu == -8 )  multipole_Q44(incr) = multipole_Q44(incr) + sumcb
          if ( mu == -6 )  multipole_Q43(incr) = multipole_Q43(incr) - sumcb
          if ( mu == -4 )  multipole_Q42(incr) = multipole_Q42(incr) + sumcb
          if ( mu == -2 )  multipole_Q41(incr) = multipole_Q41(incr) - sumcb
          if ( mu ==  0 )  multipole_Q40(incr) = multipole_Q40(incr) + sumcb
          if ( mu ==  2 )  multipole_Q41(incr) = multipole_Q41(incr) + sumcb
          if ( mu ==  4 )  multipole_Q42(incr) = multipole_Q42(incr) + sumcb
          if ( mu ==  6 )  multipole_Q43(incr) = multipole_Q43(incr) + sumcb
          if ( mu ==  8 )  multipole_Q44(incr) = multipole_Q44(incr) + sumcb
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
