!==============================================================================!
! MODULE Pairs                                                                 !
!                                                                              !
! This module contains the variables and routines related to pair operators.   !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_pairs                                                       !
!==============================================================================!
MODULE Pairs     

use Constants
use MathMethods 
use Basis

implicit none 
integer :: pairs_scheme
real(r64), dimension(:), allocatable :: pairs_T00_J1p1, & !T=0 MT= 0, J=1 MJ=+1
                                        pairs_T00_J10,  & !T=0 MT= 0, J=1 MJ= 0
                                        pairs_T00_J1m1, & !T=0 MT= 0, J=1 MJ=-1
                                        pairs_T1p1_J00, & !T=1 MT=+1, J=0 MJ= 0
                                        pairs_T10_J00,  & !T=1 MT= 0, J=0 MJ= 0
                                        pairs_T1m1_J00    !T=1 MT=-1, J=0 MJ= 0

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_pairs                                                         !
!                                                                              ! 
! Defines the HO matrix elements of the pair operators.                        ! 
! The operators are vectorized (for the calculation of constraints).           ! 
!                                                                              ! 
! pairs_scheme = 0 Agnostic (no special factor or selection rule)              ! 
!               = 1 Seniority (only pairs among particles of the same shell)   ! 
!                 see Eq(19.1) of "Simple Models of Complex Nuclei" by Talmi   ! 
!                 (ISBN: 978-3718605507)                                       ! 
!               = 2 LST coupling (among particles with same n and l)           ! 
!                 see Hinohara.2014.PhysRevC.90.031301                         ! 
!                 BB: it has never been benchmarked                            ! 
!------------------------------------------------------------------------------!
subroutine set_pairs           

integer :: ia, ib, na, la, ja, ma, ta, sha, nb, lb, jb, mb, tb, shb, incr, &
           ml, ms, ialloc=0
real(r64) :: factor, facl, facs, cg1, cg2, cg3, cg4, cg5, cg6

allocate( pairs_T00_J1p1(HOsp_dim2), pairs_T1p1_J00(HOsp_dim2), &
          pairs_T00_J10(HOsp_dim2),  pairs_T10_J00(HOsp_dim2),  &
          pairs_T00_J1m1(HOsp_dim2), pairs_T1m1_J00(HOsp_dim2), &
          stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of pair operators.'

pairs_T00_J1p1 = zero
pairs_T00_J10  = zero
pairs_T00_J1m1 = zero
pairs_T1p1_J00 = zero
pairs_T10_J00  = zero
pairs_T1m1_J00 = zero

incr = 0
do ia = 1, HOsp_dim
  na  = HOsp_n(ia)
  la  = 2*HOsp_l(ia)
  ja  = HOsp_2j(ia)
  ma  = HOsp_2mj(ia)
  ta  = HOsp_2mt(ia)
  sha = HOsp_sh(ia)
  do ib = 1, HOsp_dim
    nb  = HOsp_n(ib)
    lb  = 2*HOsp_l(ib)
    jb  = HOsp_2j(ib)
    mb  = HOsp_2mj(ib)
    tb  = HOsp_2mt(ib)
    shb = HOsp_sh(ib)
    incr = incr + 1

    !!! Conditions on quantum numbers
    if ( pairs_scheme == 1 ) then
      if ( sha /= shb ) cycle
    elseif ( pairs_scheme == 2 ) then 
      if ( (na /= nb) .or. (la /= lb) ) cycle 
    endif

    !!! Factors for different options 
    select case (pairs_scheme)
      case (0)
        factor = one
        if ( sha == shb ) factor =  1.0d0 / sqrt(2.0d0)
      case (1)
        factor = sqrt(ja + 1.0d0) / 2.0d0
      case (2) 
        factor =  1.0d0 / sqrt(2.0d0)
    end select

    !!! T=0 (isoscalar) 
    call ClebschGordan(1,1,0,ta,tb,0,cg1)
    if ( pairs_scheme < 2 ) then
      call ClebschGordan(ja,jb,2,ma,mb, 2,cg2)
      call ClebschGordan(ja,jb,2,ma,mb, 0,cg3)
      call ClebschGordan(ja,jb,2,ma,mb,-2,cg4)
    else  
      cg2 = zero
      cg3 = zero
      cg4 = zero
      do ml = -la, la, 2
        facl = (-1)**((la-ml)/2) 
        call ClebschGordan(la,1,ja, ml,1,ma,cg5)
        call ClebschGordan(lb,1,jb,-ml,1,mb,cg6)
        cg2 = cg2 + facl * cg5 * cg6
        call ClebschGordan(la,1,ja, ml,-1,ma,cg5)
        call ClebschGordan(lb,1,jb,-ml,-1,mb,cg6)
        cg4 = cg4 + facl * cg5 * cg6
        do ms = -1, 1, 2
          facs = 1.0d0 / sqrt(2.0d0)
          call ClebschGordan(la,1,ja, ml, ms,ma,cg5)
          call ClebschGordan(lb,1,jb,-ml,-ms,mb,cg6)
          cg3 = cg3 + facl * facs * cg5 * cg6 
        enddo
      enddo
    endif
    pairs_T00_J1p1(incr) = factor * cg1 * cg2
    pairs_T00_J10(incr)  = factor * cg1 * cg3
    pairs_T00_J1m1(incr) = factor * cg1 * cg4

    !!! T=1 (isovector)
    call ClebschGordan(1,1,2,ta,tb, 2,cg2)
    call ClebschGordan(1,1,2,ta,tb, 0,cg3)
    call ClebschGordan(1,1,2,ta,tb,-2,cg4)
    if ( pairs_scheme < 2 ) then
      call ClebschGordan(ja,jb,0,ma,mb,0,cg1)
    else
      cg1 = zero
      do ml = -la, la, 2
        facl = (-1)**((la-ml)/2) 
        do ms = -1, 1, 2
          facs = (-1)**((1-ms)/2) / sqrt(2.0d0)
          call ClebschGordan(la,1,ja, ml, ms,ma,cg5)
          call ClebschGordan(lb,1,jb,-ml,-ms,mb,cg6)
          cg1 = cg1 + facl * facs * cg5 * cg6
        enddo
      enddo
    endif
    pairs_T1p1_J00(incr) = factor * cg1 * cg2
    pairs_T10_J00(incr)  = factor * cg1 * cg3
    pairs_T1m1_J00(incr) = factor * cg1 * cg4

  enddo
enddo

end subroutine set_pairs           

END MODULE Pairs
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
