!==============================================================================!
! MODULE Fields                                                                !
!                                                                              !
! This module contains the variables and routines related to the fields and    !
! the matrix elements of H in the QP basis.                                    !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_fields                                                      !
! - subroutine calculate_fields                                                !
! - subroutine calculate_fields_diag                                           !
! - subroutine calculate_H11_real                                              !
! - subroutine calculate_H20                                                   !
! - subroutine calculate_H20_real                                              !
!==============================================================================!
MODULE Fields           

use Wavefunctions  
use Hamiltonian    

implicit none
public

!!! Fields in the single-particle basis
real(r64), dimension(:,:), allocatable :: field_hspRR,   & ! h^RR 
                                          field_gammaRR, & ! Gamma^RR
                                          field_deltaRR    ! Delta^RR 

complex(r64), dimension(:,:), allocatable :: field_hspLR,   & ! h^LR         
                                             field_gammaLR, & ! Gamma^LR
                                             field_deltaLR, & ! Delta^LR 
                                             field_deltaRL    ! Delta^RL

!!! Matrix elements of H in QP basis
real(r64), dimension(:,:), allocatable :: field_H11, & ! 11 part
                                          field_H20    ! 20 part 
real(r64), dimension(:), allocatable :: field_H20v     ! 20 part (vector)

!!! Additional quantities                                                        
real(r64) :: factor_delta = one ! factor to constraint delta  

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_fields                                                        !
!                                                                              !
! Initializes the fields.                                                      !
!------------------------------------------------------------------------------!
subroutine set_fields

integer :: ialloc=0

allocate( field_hspRR(HOsp_dim,HOsp_dim), field_gammaRR(HOsp_dim,HOsp_dim),   &
          field_deltaRR(HOsp_dim,HOsp_dim), field_hspLR(HOsp_dim,HOsp_dim),   &
          field_gammaLR(HOsp_dim,HOsp_dim), field_deltaLR(HOsp_dim,HOsp_dim), &
          field_deltaRL(HOsp_dim,HOsp_dim), field_H11(HOsp_dim,HOsp_dim), &
          field_H20v(HOsp_dim*HOsp_dim), field_H20(HOsp_dim,HOsp_dim), &
          stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of fields'

field_hspRR   = zero
field_gammaRR = zero
field_deltaRR = zero
field_hspLR   = zzero
field_gammaLR = zzero
field_deltaLR = zzero
field_deltaRL = zzero
field_H11 = zero
field_H20 = zero
field_H20v= zero
   
end subroutine set_fields

!------------------------------------------------------------------------------!
! subroutine calculate_fields                                                  !
!                                                                              !
! Calculates the HFB fields h, Gamma an Delta which are then used to compute   !
! other quantities of interest (in particular the energy).                     !
!                                                                              !
! What is calculated:                                                          !
!   h^LR_{ac} = Gamma^LR_{ac} + t_{ac}                                         !
!   Gamma^LR_{ac} = sum_{bd} V_{abcd} rho^LR_{db}                              !
!   Delta^LR_{ab} = 1/2 sum_{cd}  V_{abcd} kappa^LR_{cd}                       !
!                 =     sum_{c<d} V_{abcd} kappa^LR_{cd}                       !
!   Delta^RL_{ab} = 1/2 sum_{cd}  V_{abcd} kappa^RL_{cd}                       !
!                 =     sum_{c<d} V_{abcd} kappa^RL_{cd}                       !
!                                                                              !
! The fields have the symmetry properties:                                     !
!   Delta^LR_{ab} = - Delta^LR_{ba}    (Skew symmetry)                         !
!   Delta^RL_{ab} = - Delta^RL_{ba}    (Skew symmetry)                         !
!                                                                              !
! The actual calculation below uses these symmetries and the fact that only    !
! a subset of the matrix elements of the interaction are stored.               !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        rhoLR,kappaLR,kappaRL = transition densities                          !
!                                                                              !
! Output: gammaLR,deltaLR,deltaLR = transisition fields                        !
!------------------------------------------------------------------------------!
subroutine calculate_fields(rhoLR,kappaLR,kappaRL,gammaLR,hspLR,deltaLR, &
                            deltaRL,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR, kappaRL
complex(r64), dimension(ndim,ndim), intent(out) :: gammaLR, hspLR, deltaLR, & 
                                                   deltaRL
integer :: i, j, ia, ib, ic, id, it, perm
integer(i64) :: kk
real(r64) :: h2b, f2b
!cmpi integer :: ierr=0
!cmpi complex(r64), dimension(ndim,ndim) :: gammaLR_red, deltaLR_red, & 
!cmpi                                       deltaRL_red

gammaLR = zzero
deltaLR = zzero
deltaRL = zzero

!$OMP PARALLEL DO FIRSTPRIVATE(rhoLR,kappaLR,kappaRL) &
!$OMP             PRIVATE(ia,ib,ic,id,h2b,f2b,perm,it) &
!$OMP             REDUCTION(+:gammaLR,deltaLR,deltaRL)
do kk = 1, hamil_H2dim  
  ia = hamil_abcd(1+4*(kk-1))
  ib = hamil_abcd(2+4*(kk-1))
  ic = hamil_abcd(3+4*(kk-1))
  id = hamil_abcd(4+4*(kk-1))
  h2b = hamil_H2(kk)
  perm = hamil_trperm(kk)

  !!! Loop on time reversal
  do it = 1, 2 
    if ( it == 2 ) then
      if ( HOsp_2mj(ia) + HOsp_2mj(ib) == 0 ) cycle
      call find_timerev(perm,ia,ib,ic,id)
      h2b = sign(one,perm*one) * h2b
    endif

    !!! Faster than using if ((a /= c).or.(b /= d))
    f2b = h2b * (1 - kdelta(ia,ic) * kdelta(ib,id))

    !!! Calculation of Gamma
    gammaLR(ia,ic) = gammaLR(ia,ic) + h2b * rhoLR(id,ib)
    gammaLR(ia,id) = gammaLR(ia,id) - h2b * rhoLR(ic,ib)
    gammaLR(ib,ic) = gammaLR(ib,ic) - h2b * rhoLR(id,ia)
    gammaLR(ib,id) = gammaLR(ib,id) + h2b * rhoLR(ic,ia)

    gammaLR(ic,ia) = gammaLR(ic,ia) + f2b * rhoLR(ib,id)
    gammaLR(ic,ib) = gammaLR(ic,ib) - f2b * rhoLR(ia,id)
    gammaLR(id,ia) = gammaLR(id,ia) - f2b * rhoLR(ib,ic)
    gammaLR(id,ib) = gammaLR(id,ib) + f2b * rhoLR(ia,ic)
  
    !!! Calculation of Delta^10 and Delta^01  
    deltaLR(ib,ia) = deltaLR(ib,ia) + h2b * kappaLR(id,ic)
    deltaRL(ib,ia) = deltaRL(ib,ia) + h2b * kappaRL(id,ic)
  
    deltaLR(id,ic) = deltaLR(id,ic) + f2b * kappaLR(ib,ia)
    deltaRL(id,ic) = deltaRL(id,ic) + f2b * kappaRL(ib,ia)

  enddo
enddo
!$OMP END PARALLEL DO
       
!!! Reduces the values for the processes in the same team
!cmpi if ( paral_myteamsize > 1 ) then 
!cmpi   call mpi_reduce(gammaLR,gammaLR_red,ndim**2,mpi_double_complex, &
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)
!cmpi   call mpi_reduce(deltaLR,deltaLR_red,ndim**2,mpi_double_complex, &
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)
!cmpi   call mpi_reduce(deltaRL,deltaRL_red,ndim**2,mpi_double_complex, &
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)
!cmpi   gammaLR = gammaLR_red
!cmpi   deltaLR = deltaLR_red
!cmpi   deltaRL = deltaRL_red
!cmpi endif
       
!!! h = Gamma + 1body
do j = 1, HOsp_dim
  do i = 1, HOsp_dim
    hspLR(i,j) = gammaLR(i,j) + hamil_H1(i,j)
  enddo
enddo

!!! Skew symmetry 
do j = 1, HOsp_dim
  do i = 1, j-1      
    deltaLR(i,j) = -1.0d0 * deltaLR(j,i) 
    deltaRL(i,j) = -1.0d0 * deltaRL(j,i) 
  enddo
enddo

end subroutine calculate_fields

!------------------------------------------------------------------------------!
! subroutine calculate_fields_diag                                             !
!                                                                              !
! Same as calculate_fields but for same left and right states. In that case,   !
! one can use additional symmetries to simplify the calculation:               !
!   Gamma^LR_{ac} = (Gamma^LR_{ca})^*  (Hermicity)                             !
!   h^LR_{ac} = (h^LR_{ca})^*          (Hermicity)                             !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        rhoLR,kappLR = transition densities                                   !
!                                                                              !
! Output: gammaLR,deltaLR,deltaLR = transisition fields                        !
!------------------------------------------------------------------------------!
subroutine calculate_fields_diag(rhoLR,kappaLR,gammaLR,hspLR,deltaLR,deltaRL, &
                                 ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaLR
complex(r64), dimension(ndim,ndim), intent(out) :: gammaLR, hspLR, deltaLR
complex(r64), dimension(ndim,ndim), intent(out), optional :: deltaRL
integer :: i, j, ia, ib, ic, id, it, perm
integer(i64) :: kk
complex(r64) :: h2b
complex(r64), dimension(ndim,ndim) :: delta_tmp
!cmpi integer :: ierr=0
!cmpi complex(r64), dimension(ndim,ndim) :: gammaLR_red, deltaLR_red

gammaLR = zzero
deltaLR = zzero

!$OMP PARALLEL DO FIRSTPRIVATE(rhoLR,kappaLR) &
!$OMP             PRIVATE(ia,ib,ic,id,h2b,perm,it) &
!$OMP             REDUCTION(+:gammaLR,deltaLR)
do kk = 1, hamil_H2dim  
  ia = hamil_abcd(1+4*(kk-1))
  ib = hamil_abcd(2+4*(kk-1))
  ic = hamil_abcd(3+4*(kk-1))
  id = hamil_abcd(4+4*(kk-1))
  h2b = hamil_H2(kk)
  perm = hamil_trperm(kk)
 
  !!! Loop on time reversal
  do it = 1, 2 
    if ( it == 2 ) then
      if ( HOsp_2mj(ia) + HOsp_2mj(ib) == 0 ) cycle
      call find_timerev(perm,ia,ib,ic,id)
      h2b = sign(one,perm*one) * h2b
    endif

    !!! Calculation of Gamma
    gammaLR(ia,ic) = gammaLR(ia,ic) + h2b * rhoLR(id,ib)
    gammaLR(ia,id) = gammaLR(ia,id) - h2b * rhoLR(ic,ib)
 
    if ( (ic == ia) .and. (ib /= id) ) then 
      gammaLR(ic,ia) = gammaLR(ic,ia) + h2b * rhoLR(ib,id)
    endif                
    if ( ib <= id ) then 
      gammaLR(ib,id) = gammaLR(ib,id) + h2b * rhoLR(ic,ia)
    endif                
    if ( ib <= ic ) then 
      gammaLR(ib,ic) = gammaLR(ib,ic) - h2b * rhoLR(id,ia)
    endif                
 
    if ( (ia /= ic) .or. (ib /= id) ) then 
      if ( id <= ib ) then 
        gammaLR(id,ib) = gammaLR(id,ib) + h2b * rhoLR(ia,ic)
      endif                
      if ( ic <= ib ) then 
        gammaLR(ic,ib) = gammaLR(ic,ib) - h2b * rhoLR(ia,id)
      endif                
    endif                
 
    !!! Calculation of Delta^10 and Delta^01  
    deltaLR(ib,ia) = deltaLR(ib,ia) + h2b * kappaLR(id,ic)
 
    if ( (ia /= ic) .or. (ib /= id) ) then 
      deltaLR(id,ic) = deltaLR(id,ic) + h2b * kappaLR(ib,ia)
    endif                

  enddo
enddo
!$OMP END PARALLEL DO

!!! Reduces the values for the processes in the same team
!cmpi if ( paral_myteamsize > 1 ) then 
!cmpi   call mpi_reduce(gammaLR,gammaLR_red,ndim**2,mpi_double_complex, &
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)
!cmpi   call mpi_reduce(deltaLR,deltaLR_red,ndim**2,mpi_double_complex, &
!cmpi                   mpi_sum,0,mpi_comm_team,ierr)
!cmpi   gammaLR = gammaLR_red
!cmpi   deltaLR = deltaLR_red
!cmpi endif

!!! Hermicity 
do j = 1, ndim
  do i = j+1, ndim  
    gammaLR(i,j) = conjg(gammaLR(j,i)) 
  enddo
enddo

!!! h = Gamma + 1-body
do j = 1, ndim
  do i = 1, ndim
    hspLR(i,j) = gammaLR(i,j) + hamil_H1(i,j)
  enddo
enddo

!!! Skew symmetry 
do j = 1, ndim
  do i = 1, j-1      
    deltaLR(i,j) = -1.d0 * deltaLR(j,i) 
  enddo
enddo

!!! Complex conjuate. Although this looks complicated, I do it this way to 
!!! avoid a bug when compiling with option -O3.
do j = 1, ndim
  do i = 1, ndim     
    delta_tmp(i,j) = conjg(deltaLR(i,j))
  enddo
enddo

if ( present(deltaRL) ) then
  deltaRL = delta_tmp
endif

end subroutine calculate_fields_diag

!------------------------------------------------------------------------------!
! subroutine calculate_H11_real                                                !
!                                                                              !
! Computes the H11 part of the Hamiltonian in the QP basis, i.e.               !
!   H11 = U^t h U - V^t h^t V + U^t delta V - V^t delta U                      !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!------------------------------------------------------------------------------!
subroutine calculate_H11_real(ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim) :: A1, A2, A3, A4, A5

!!! U^t h U
call dgemm('t','n',ndim,ndim,ndim,one,bogo_U0,ndim,field_hspRR,ndim,zero,A1, &
           ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,bogo_U0,ndim,zero,A2,ndim)

!!! V^t h^t V 
call dgemm('t','t',ndim,ndim,ndim,one,bogo_V0,ndim,field_hspRR,ndim,zero,A1, &
           ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,bogo_V0,ndim,zero,A3,ndim)

!!! U^t delta V
call dgemm('t','n',ndim,ndim,ndim,one,bogo_U0,ndim,field_deltaRR,ndim,zero,A1, &
           ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,bogo_V0,ndim,zero,A4,ndim)

!!! V^t delta U
A5 = -transpose(A4)

!!! H11
field_H11 = A2 - A3 + (A4 - A5) * factor_delta

end subroutine calculate_H11_real

!------------------------------------------------------------------------------!
! subroutine calculate_H20                                                     !
!                                                                              !
! Calculates the 20 part of the Hamiltonian in the QP basis                    !
!   H20 = U^t h V - V^t h^t U + U^t deltaLR U - V^t deltaRL V                  !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        U,V = Bogoliubov matrices                                             !
!                                                                              !
! Output: H20 = 20 part of the Hamiltonian in the qp basis                     !
!------------------------------------------------------------------------------!
subroutine calculate_H20(U,V,H20,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: U, V
complex(r64), dimension(ndim,ndim), intent(out) :: H20
complex(r64), dimension(ndim,ndim) :: A1, A2, A3, A4, A5

!!! U^t * hLR * V
call zgemm('t','n',ndim,ndim,ndim,zone,U,ndim,field_hspLR,ndim,zzero,A1,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,A1,ndim,V,ndim,zzero,A2,ndim)

!!! V^t * hLR^t * U 
A3 = transpose(A2)

!!! U^t * deltaLR * U
call zgemm('t','n',ndim,ndim,ndim,zone,U,ndim,field_deltaLR,ndim,zzero,A1,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,A1,ndim,U,ndim,zzero,A4,ndim)

!!! V^t * deltaRL * V
call zgemm('t','n',ndim,ndim,ndim,zone,V,ndim,field_deltaRL,ndim,zzero,A1,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,A1,ndim,V,ndim,zzero,A5,ndim)

!!! H20
H20 = A2 - A3 + (A4 - A5) * factor_delta

end subroutine calculate_H20 

!------------------------------------------------------------------------------!
! subroutine calculate_H20_real                                                !
!                                                                              !
! Computes the H20 part of the Hamiltonian in the QP basis, i.e.               !
!   H20 = U^t h V - V^t h^t U + U^t delta U - V^t delta V                      !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!------------------------------------------------------------------------------!
subroutine calculate_H20_real(ndim)

integer, intent(in) :: ndim
integer :: i, j, k
real(r64), dimension(ndim,ndim) :: A1, A2, A3, A4, A5

!!! U^t hRR V
call dgemm('t','n',ndim,ndim,ndim,one,bogo_U0,ndim,field_hspRR,ndim,zero,A1, &
           ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,bogo_V0,ndim,zero,A2,ndim)

!!! V^t hRR^t U   
A3 = transpose(A2)

!!! U^t deltaRR U
call dgemm('t','n',ndim,ndim,ndim,one,bogo_U0,ndim,field_deltaRR,ndim,zero,A1, &
           ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,bogo_U0,ndim,zero,A4,ndim)

!!! V^t deltaRR V
call dgemm('t','n',ndim,ndim,ndim,one,bogo_V0,ndim,field_deltaRR,ndim,zero,A1, & 
           ndim)
call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,bogo_V0,ndim,zero,A5,ndim)

!!! H20
field_H20 = A2 - A3 + (A4 - A5) * factor_delta

k = 0
do j = 1, ndim
  do i = 1, ndim
    k = k + 1
    field_H20v(k) = field_H20(i,j)
  enddo
enddo 

end subroutine calculate_H20_real

!------------------------------------------------------------------------------!
! subroutine broadcast_densities                                               !
!                                                                              !
! Broadcasts the Bogoliubov matrices and diagonal densities to all processes.  !
!------------------------------------------------------------------------------!
!cmpi subroutine broadcast_densities

!cmpi integer :: ndim2, ierr=0

!cmpi ndim2 = HOsp_dim**2

!cmpi call mpi_bcast(bogo_U0,ndim2,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(bogo_V0,ndim2,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(dens_rhoRR,ndim2,mpi_double_precision,0,mpi_comm_world, &
!cmpi                ierr)
!cmpi call mpi_bcast(dens_kappaRR,ndim2,mpi_double_precision,0,mpi_comm_world, &
!cmpi                ierr)

!cmpi end subroutine broadcast_densities

END MODULE Fields        
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
