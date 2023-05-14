!==============================================================================!
! MODULE AngularMomentum                                                       !
!                                                                              !
! This module contains the variables and routines related to angular momentum  !
! operators.                                                                   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_angular_momentum                                            !
!==============================================================================!
MODULE AngularMomentum

use Basis

implicit none 
public

real(r64), dimension(:), allocatable :: angumome_Jx,  & ! Jx
                                        angumome_Jy,  & ! Jy
                                        angumome_Jz,  & ! Jz
                                        angumome_Jx2, & ! Jx^2
                                        angumome_Jy2, & ! Jy^2
                                        angumome_Jz2, & ! Jz^2
                                        angumome_so     ! spin-orbit

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_angular_momentum                                              !
!                                                                              ! 
! Defines the HO matrix elements of the angular-momentum operators.            ! 
! The operators are vectorized (for the calculation of constraints).           ! 
!------------------------------------------------------------------------------!
subroutine set_angular_momentum

integer :: ia, ina, ila, ija, ima, ita, ib, inb, ilb, ijb, imb, itb, &
           ja, la, incr, ialloc=0
real(r64) :: rja, rma, rjb, rmb, xja, &
             fp1, fm1, fp1p2, fm1m2, f02, &
             delta_ma_mbp1, delta_ma_mbm1, delta_ma_mb, &
             delta_ma_mbp2, delta_ma_mbm2 

allocate( angumome_Jx(HOsp_dim2), angumome_Jx2(HOsp_dim2), &
          angumome_Jy(HOsp_dim2), angumome_Jy2(HOsp_dim2), &
          angumome_Jz(HOsp_dim2), angumome_Jz2(HOsp_dim2), &
          stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of angular-momentum operators'

angumome_Jx  = zero
angumome_Jy  = zero
angumome_Jz  = zero
angumome_Jx2 = zero
angumome_Jy2 = zero
angumome_Jz2 = zero

incr = 0
do ia = 1, HOsp_dim
  ina = HOsp_n(ia)
  ila = HOsp_l(ia)
  ija = HOsp_2j(ia)
  ima = HOsp_2mj(ia)
  ita = HOsp_2mt(ia)
  rja = ija/2.d0
  rma = ima/2.d0
  do ib = 1, HOsp_dim
    incr = incr + 1
    inb = HOsp_n(ib)
    ilb = HOsp_l(ib)
    ijb = HOsp_2j(ib)
    imb = HOsp_2mj(ib)
    itb = HOsp_2mt(ib)
    rjb = ijb/2.d0
    rmb = imb/2.d0

    fp1 = rjb*(rjb+1.d0) - rmb*(rmb+1.d0)
    fm1 = rjb*(rjb+1.d0) - rmb*(rmb-1.d0)
    fp1p2 = rjb*(rjb+1.d0) - (rmb+1.d0)*(rmb+2.d0)
    fm1m2 = rjb*(rjb+1.d0) - (rmb-1.d0)*(rmb-2.d0)
    f02 = rjb*(rjb+1.d0) - rmb**2

    if ( (ina == inb) .and. (ila == ilb) .and. (ija == ijb) .and. &
         (ita == itb) ) then
      delta_ma_mbp2 = one * kdelta(ima,imb+4)
      delta_ma_mbp1 = one * kdelta(ima,imb+2)
      delta_ma_mb   = one * kdelta(ima,imb)
      delta_ma_mbm1 = one * kdelta(ima,imb-2)
      delta_ma_mbm2 = one * kdelta(ima,imb-4)
      
      angumome_Jx(incr) = 0.5d0 * ( sqrt(fp1)*delta_ma_mbp1 & 
                                  + sqrt(fm1)*delta_ma_mbm1 )
      angumome_Jy(incr) = 0.5d0 * ( sqrt(fp1)*delta_ma_mbp1 &
                                  - sqrt(fm1)*delta_ma_mbm1 )
      angumome_Jz(incr) = rmb * delta_ma_mb
      
      angumome_Jx2(incr) = 0.25d0 * ( sqrt(fp1*fp1p2)*delta_ma_mbp2 &
                          + sqrt(fm1*fm1m2)*delta_ma_mbm2 + 2*f02*delta_ma_mb )
      angumome_Jy2(incr) = 0.25d0 * ( sqrt(fp1*fp1p2)*delta_ma_mbp2 &
                          + sqrt(fm1*fm1m2)*delta_ma_mbm2 - 2*f02*delta_ma_mb )
      angumome_Jz2(incr) = (rmb * delta_ma_mb)**2
    endif

  enddo
enddo

!!! Spin-orbit
allocate( angumome_so(HOsp_dim2), source=zero, stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of spin-orbit operator'

incr = 0

do ia = 1, HOsp_dim
  la = HOsp_l(ia)
  ja = HOsp_2j(ia)
  xja = ja / 2.d0
  do ib = 1, HOsp_dim
    incr = incr + 1
    if ( ia /= ib ) cycle
    angumome_so(incr) = xja*(xja+1) - la*(la+1) - 0.75d0 
  enddo
enddo

end subroutine set_angular_momentum

END MODULE AngularMomentum
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
