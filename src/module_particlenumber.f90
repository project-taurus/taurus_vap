!==============================================================================!
! MODULE ParticleNumber                                                        !
!                                                                              !
! This module contains the variables and routines related to particle-number   !
! operators (including their square and occupation number).                    !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_particle_number                                             !
! - subroutine calculate_particle_number                                       !
! - subroutine calculate_occupation_number                                     !
!==============================================================================!
MODULE ParticleNumber

use Basis

implicit none 

real(r64), dimension(:), allocatable :: partnumb_Z, & ! Proton  number operator
                                        partnumb_N, & ! Neutron    "      "    
                                        partnumb_A    ! Nucleon    "      "    
CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_particle_number                                               !
!                                                                              !
! Defines the HO matrix elements of the proton, neutron and nucleon number op- !
! erators. They are simply the identity matrices in their relative subspace.   !
! The operators are vectorized (for the calculation of constraints).           !
!------------------------------------------------------------------------------!
subroutine set_particle_number

integer :: i, j, k, ialloc=0

allocate( partnumb_Z(HOsp_dim2), partnumb_N(HOsp_dim2), partnumb_A(HOsp_dim2), &
          stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of particle-number operator'

partnumb_Z = zero 
partnumb_N = zero 
partnumb_A = zero 

k = 0
do i = 1, HOsp_dim
  do j = 1, HOsp_dim
    k = k + 1
    if ( i /= j ) cycle
    if ( i <= HOsp_dim/2 ) then
      partnumb_Z(k) = one
    else 
      partnumb_N(k) = one
    endif
  enddo
enddo

partnumb_A = partnumb_Z + partnumb_N

end subroutine set_particle_number

!------------------------------------------------------------------------------!
! subroutine calculate_particle_number                                         !
!                                                                              !
! Computes the expectation values of the proton and neutron number operators.  !
! < Z > = Tr(rhoLR_p) = Tr(rhoLR) from 1            to HOsp_dim/2              !
! < N > = Tr(rhoLR_n) = Tr(rhoLR) from 1+HOsp_dim/2 to HOsp_dim                !
! < Z^2 > = Tr(rhoLR_p) + Tr(rhoLR_p)^2 - Tr(rhoLR_p^2)                        !
!           + Tr(kappaRL^T_p kappaLR_p)                                        !
! < N^2 > = Tr(rhoLR_n) + Tr(rhoLR_n)^2 - Tr(rhoLR_n^2)                        !
!           + Tr(kappaRL^T_n kappaLR_n)                                        !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        io = 0 computes < Z >, < N >, < Z^2 >, < N^2 >                        !
!           = 1 computes < Z >, < N >                                          !
!        rhoLR,kappaLR,kappaRL = transition densities                          !
!                                                                              !
! Output: prot, prot2 = < Z > and < Z^2 >                                      !
!         neut, neut2 = < N > and < N^2 >                                      !
!------------------------------------------------------------------------------!
subroutine calculate_particle_number(io,rhoLR,kappaLR,kappaRL,prot,neut, & 
                                     prot2,neut2,ndim) 

integer, intent(in) :: ndim, io
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), intent(out) :: prot, neut, prot2, neut2
integer :: hdim, j, i 
complex(r64) :: tr1, tr2, tr3, tr4
complex(r64), dimension(ndim/2,ndim/2) :: rhoLRp, kapLRp, kapRLp, A1, A2, &
                                          rhoLRn, kapLRn, kapRLn, A3, A4

if ( (io /= 0) .and. (io /= 1) ) then 
  print*,'Wrong argument in calculate_particle_number: io = ', io 
  stop 
endif

hdim = ndim/2

prot  = zzero
neut  = zzero
prot2 = zzero
neut2 = zzero

!!! N and Z
do i = 1, hdim  
  prot = prot + rhoLR(i,i)
  neut = neut + rhoLR(i+hdim,i+hdim)
enddo

!!! N^2 and Z^2
if ( io == 0 ) then 

  do j = 1, hdim
    do i = 1, hdim 
      rhoLRp(i,j) = rhoLR(i,j)
      rhoLRn(i,j) = rhoLR(i+hdim,j+hdim)
      kapLRp(i,j) = kappaLR(i,j)
      kapLRn(i,j) = kappaLR(i+hdim,j+hdim)
      kapRLp(i,j) = kappaRL(i,j)
      kapRLn(i,j) = kappaRL(i+hdim,j+hdim)
    enddo
  enddo

  call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRp,hdim,rhoLRp,hdim,zzero,A1,hdim)
  call zgemm('t','n',hdim,hdim,hdim,zone,kapRLp,hdim,kapLRp,hdim,zzero,A2,hdim)
  call zgemm('n','n',hdim,hdim,hdim,zone,rhoLRn,hdim,rhoLRn,hdim,zzero,A3,hdim)
  call zgemm('t','n',hdim,hdim,hdim,zone,kapRLn,hdim,kapLRn,hdim,zzero,A4,hdim)

  tr1 = zzero
  tr2 = zzero
  tr3 = zzero
  tr4 = zzero

  do i = 1, hdim  
    tr1 = tr1 + A1(i,i)
    tr2 = tr2 + A2(i,i)
    tr3 = tr3 + A3(i,i)
    tr4 = tr4 + A4(i,i)
  enddo

  prot2 = prot + prot**2 - tr1 + tr2
  neut2 = neut + neut**2 - tr3 + tr4

endif

end subroutine calculate_particle_number

!------------------------------------------------------------------------------!
! subroutine calculate_occupation_number                                       !
!                                                                              !
! Calculates the occupation numbers of the spherical shell, i.e. the number    !
! of particcle in all the 2j+1 multiplet forming a shell.                      !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        sdim = dimension of the shells                                        !
!        rhoLR = transition density                                            !
!                                                                              !
! Output: occnum = occupation numbers                                          !
!------------------------------------------------------------------------------!
subroutine calculate_occupation_number(rhoLR,occnum,ndim,sdim)

integer, intent(in) :: ndim, sdim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR
complex(r64), dimension(sdim,2), intent(out) :: occnum
integer :: i, j, k, jmax

k = 0
occnum = zzero

do i = 1, sdim
  jmax = HOsh_2j(i) + 1
  do j = 1, jmax
    k = k + 1
    occnum(i,1) = occnum(i,1) + rhoLR(k,k)
    occnum(i,2) = occnum(i,2) + rhoLR(k+ndim/2,k+ndim/2)
  enddo
enddo

end subroutine calculate_occupation_number 

END MODULE ParticleNumber
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
