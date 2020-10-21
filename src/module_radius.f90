!==============================================================================!
! MODULE Radius                                                                !
!                                                                              !
! This module contains the variables and routines related to radius squared    !
! one-body operator.                                                           !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_radius                                                      !
!==============================================================================!
MODULE Radius 

use Basis

implicit none 
public

real(r64), dimension(:), allocatable :: radius_r2 ! squared radius

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_radius                                                        !
!                                                                              ! 
! Defines the HO matrix elements of the radius square one-body operator.       ! 
! The operator is vectorized (for the calculation of constraints).             ! 
!------------------------------------------------------------------------------!
subroutine set_radius        

integer :: incr, ia, ib, ja, jb, mja, mjb, la, lb, mta, mtb, ialloc=0

allocate( radius_r2(HOsp_dim2), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of square radius operator'

radius_r2 = zero

incr = 0

do ia = 1, HOsp_dim
  la = HOsp_l(ia)
  ja = HOsp_2j(ia)
  mja = HOsp_2mj(ia)
  mta = HOsp_2mt(ia)
  do ib = 1, HOsp_dim
    lb = HOsp_l(ib)
    jb = HOsp_2j(ib)
    mjb = HOsp_2mj(ib)
    mtb = HOsp_2mt(ib)
    incr = incr + 1
    if ( (ja /= jb) .or. (mja /= mjb) .or. (mta /= mtb) .or. (la /= lb) ) cycle
    radius_r2(incr) = radial(ib,ia,2)
  enddo
enddo

end subroutine set_radius     

END MODULE Radius    
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
