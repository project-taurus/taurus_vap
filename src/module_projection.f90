!==============================================================================!
! MODULE Projection                                                            !
!                                                                              !
! This module contains the variables and routines related to the projection    !
! onto good particle number (and more in the future)                           !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_projection                                                  !
! - subroutine project_wavefunction                                            !
! - subroutine reset_pnp                                                       !
! - subroutine sum_gauge                                                       !
! - subroutine generate_rotation_gauge                                         !
! - subroutine generate_rotation_parity                                        !
! - subroutine rotate_wavefunction                                             !
! - subroutine calculate_thouless                                              !
! - subroutine calculate_projected_H20                                         !
! - subroutine calculate_overlap                                               !
! - subroutine calculate_norm                                                  !
! - subroutine print_results                                                   !
!==============================================================================!
MODULE Projection    

use Fields
use Operators      

implicit none
public

!!! Number of points in the particle-number projections 
integer :: proj_Mphip, & ! number of angles in the projection for protons
           proj_Mphin, & !   "    "    "    "   "      "       "  neutrons
           proj_myangles !   "    "    "    "   "  loop (for routine/process)
real(r64), dimension(:), allocatable :: proj_phip, &    ! proton  gauge angles
                                        proj_phin       ! neutron   "     "   
complex(r64), dimension(:), allocatable :: proj_weip, & ! proton    "   weights
                                           proj_wein    ! neutron   "      "   
real(r64) :: proj_facpi=2*pi ! factor for integral reduction      

!!! Rotated quantities
complex(r64) :: rot_over,    & ! overlap
                rot_ener,    & ! energy
                rot_pari,    & ! parity
                rot_rad2(2), & ! radius^2 (1=p,2=n)
                rot_prot,    & ! proton number 
                rot_neut,    & ! neutron number
                rot_prot2,   & ! proton variance
                rot_neut2,   & ! neutron variance
                rot_amj(3),  & ! angular momentum J_i (1=x, 2=y, 3=z)
                rot_amj2(3), & ! angular momentum J_i^2 (1=x, 2=y, 3=z) 
                rot_spor(2), & ! spin-orbit (1=p,2=n)
                rot_Qlm(2,0:4,4) ! multipole Qlm (1=p,2=n), m, l
complex(r64), dimension(12) :: rot_ecomp ! "components" of the rotated energy    
complex(r64), dimension(:,:), allocatable :: rot_occnum, & ! occupation numbers
                                             rot_H20, & ! Quantities for the
                                             rot_A,   & ! projected gradient
                                             ROT ! Rotation matrix
!!! Projected quantities 
complex(r64) :: pnp_over,    & ! overlap
                pnp_ener,    & ! energy
                pnp_pari,    & ! parity
                pnp_rad2(2), & ! radius^2 (1=p,2=n) 
                pnp_prot,    & ! proton number 
                pnp_neut,    & ! neutron number
                pnp_prot2,   & ! proton variance
                pnp_neut2,   & ! neutron variance
                pnp_amj(3),  & ! angular momentum J_i (1=x, 2=y, 3=z)
                pnp_amj2(3), & ! angular momentum J_i^2 (1=x, 2=y, 3=z) 
                pnp_spor(2), & ! spin-orbit (1=p,2=n) 
                pnp_Qlm(2,0:4,4) ! multipole Qlm (1=p,2=n), m, l
complex(r64), dimension(12) :: pnp_ecomp ! "components" of the projected energy  
complex(r64), dimension(:,:), allocatable :: pnp_occnum, & ! occupation numbers
                                             pnp_H20,  & ! quantities for the
                                             pnp_A,    & ! projected gradient
                                             pnp_AH20, &
                                             pnp_rho     ! projected density

CONTAINS 

!------------------------------------------------------------------------------!
! subroutine set_projection                                                    !
!                                                                              ! 
! Initializes the arrays and weights for the projection on particle number.    ! 
!------------------------------------------------------------------------------!
subroutine set_projection

integer :: i, ialloc=0

allocate( rot_occnum(HOsh_dim,2), pnp_occnum(HOsh_dim,2), &
          rot_H20(HOsp_dim,HOsp_dim), pnp_H20(HOsp_dim,HOsp_dim), & 
          rot_A(HOsp_dim,HOsp_dim), pnp_A(HOsp_dim,HOsp_dim), & 
          pnp_AH20(HOsp_dim,HOsp_dim), ROT(HOsp_dim,HOsp_dim), &
          pnp_rho(HOsp_dim,HOsp_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of arrays for the projection'

ROT = zzero
do i = 1, HOsp_dim
 ROT(i,i) = zone
enddo

end subroutine set_projection

!------------------------------------------------------------------------------!
! subroutine assign_angles                                                     !
!                                                                              ! 
! Determines the angles (with their weights) to be calculated. In case of MPI, ! 
! distributes the angles between the different teams.                          ! 
!                                                                              ! 
! Input: Mphip = number of angles for the discretization of proton   PNR       ! 
!        Mphin =   "    "    "     "   "        "        "  neutrons  "        ! 
!------------------------------------------------------------------------------!

subroutine assign_angles(Mphip,Mphin)

integer, intent(in) :: Mphip, Mphin
integer :: i, j, ialloc=0, jp, jn, myangles, myoffset
real(r64) :: phip, phin
!cmpi integer :: divide, rest

proj_myangles = Mphip * Mphin

myangles = proj_myangles
myoffset = 0

!!! MPI distribution of angles
!cmpi divide = proj_myangles / paral_teams
!cmpi rest = modulo(proj_myangles,paral_teams)

!cmpi paral_myangles = divide
!cmpi paral_myoffset = paral_myangles * paral_myteam

!cmpi if ( paral_myteam < rest ) then
!cmpi   paral_myangles = paral_myangles + 1
!cmpi   paral_myoffset = paral_myoffset + paral_myteam
!cmpi else
!cmpi   paral_myoffset = paral_myoffset + rest
!cmpi endif

!cmpi proj_myangles = paral_myangles
!cmpi myangles = paral_myangles
!cmpi myoffset = paral_myoffset

!if ( paral_myrank == 0 ) then
!print*, "divide,rest", divide, rest
!endif
!print*, "myrank,myoffeset",paral_myrank, myoffset

!!! Store the angles and weights in arrays
allocate( proj_phip(0:myangles), proj_weip(0:myangles), &
          proj_phin(0:myangles), proj_wein(0:myangles), &
          stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of arrays for the angles'

proj_phip = zero
proj_phin = zero
proj_weip = zzero
proj_weip = zzero

i = 0
j = 0
do jp = 1, Mphip
  do jn = 1, Mphin
    i = i + 1
    if ( (i < myoffset+1) .or. (i > myoffset+myangles) ) cycle
      j = j + 1

      if ( modulo(Mphip,2) == 1 ) then
        phip = proj_facpi * (jp-one) / Mphip
      else
        phip = proj_facpi * (jp-one/2) / Mphip
      endif

      if ( modulo(Mphin,2) == 1 ) then
        phin = proj_facpi * (jn-one) / Mphin
      else
        phin = proj_facpi * (jn-one/2) / Mphin
      endif

      proj_phip(j) = phip
      proj_phin(j) = phin
      proj_weip(j) = exp(-zimag * phip * valence_Z) / Mphip
      proj_wein(j) = exp(-zimag * phin * valence_N) / Mphin
  enddo
enddo

end subroutine assign_angles 

!------------------------------------------------------------------------------!
! subroutine project_wavefunction                                              !
!                                                                              ! 
! Performs the particle-number projection of the Bogoliubov state at the curr- ! 
! ent iteration. The routine depends on the option iopt that controls what is  ! 
! calculated.                                                                  ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        iopt = 0 computes quantities for the gradient but a few expect. val.  ! 
!             = 1 computes many expectation values for the final printing.     ! 
!        iprint = option to control the printing at each iteration.            ! 
!        Mphip = number of angles for the discretization of proton   PNR       ! 
!        Mphin =   "    "    "     "   "        "        "  neutrons  "        ! 
!------------------------------------------------------------------------------!
subroutine project_wavefunction(iopt,iprint,Mphip,Mphin,ndim) 

integer, intent(in) :: ndim, iopt, iprint, Mphip, Mphin
integer :: i, j, nocc0, nemp0, nangle, nangle_min, nangle_max, ialloc=0
real(r64) :: phip, phin, ovac0
real(r64), dimension(:), allocatable :: voveru0
complex(r64) :: weip, wein, amjx_p, amjx_n, amjy_p, amjy_n, amjz_p, amjz_n
complex(r64), dimension(ndim,ndim) :: bogo_zU0bar, bogo_zV0bar, & 
                                      bogo_zU0tilde, bogo_zV0tilde, ROTG
complex(r64), dimension(ndim**2) :: Qlm, rad2

!!! Initialization: sets most gauge-dependent quantity to zero
call reset_pnp(iopt)
bogo_zU0 = zone * bogo_U0
bogo_zV0 = zone * bogo_V0

!!! Prepares the calculation of the overlap: builds the canonical basis, stores
!!! the occupancies and the number of occupied/empty states
call construct_canonical_basis(bogo_U0,bogo_V0,bogo_zU0c,bogo_zV0c,bogo_zD0, &
                               ovac0,nocc0,nemp0,ndim)
if ( ndim-nocc0-nemp0 == 0 ) ovac0 = 0.0d0

allocate(voveru0(ndim-nocc0), source=zero, stat=ialloc)
if ( ialloc /= 0 ) stop 'Error during allocation of array for the overlap'

j = 1
do i = 1+nemp0, ndim-nocc0
  if ( (-1)**j == -1 ) then 
    voveru0(i) = real( bogo_zV0c(i,i+1) / bogo_zU0c(i,i) )
  else
    voveru0(i) = real( bogo_zV0c(i,i-1) / bogo_zU0c(i,i) )
  endif
  j = j + 1
enddo

!!! Loop over gauge angles                                    
call assign_angles(Mphip,Mphin)
  
nangle_min = min(1,proj_myangles)
nangle_max = proj_myangles
  
do nangle = nangle_min, nangle_max
  phip = proj_phip(nangle)
  phin = proj_phin(nangle)
  weip = proj_weip(nangle)
  wein = proj_wein(nangle)

  !!! Rotates the wave function 
  call generate_rotation_gauge(phip,phin,ROT,ROTG,ndim)
  call rotate_wavefunction(ROTG,bogo_zU0,bogo_zV0,bogo_zU0bar,bogo_zV0bar,& 
                           ndim)

  !!! Computes all the important quantites (overlap, densities, fields)
  !!! between for different left and right states 
  if ( max(Mphip,Mphin) > 1 ) then 
    call calculate_overlap(nocc0,nemp0,ovac0,voveru0,bogo_zD0, &
                           nocc0,nemp0,ovac0,voveru0,bogo_zD0, &
                           ROTG,rot_over,ndim)
  else
    rot_over = zone
  endif
  call calculate_thouless(bogo_zU0,bogo_zV0,bogo_zU0bar,bogo_zV0bar,rot_A,&
                          bogo_zU0tilde,bogo_zV0tilde,ndim)
  call calculate_densities(bogo_zU0,bogo_zV0,conjg(bogo_zU0tilde), &
                           conjg(bogo_zV0tilde),dens_rhoLR,dens_kappaLR, &
                           dens_kappaRL,ndim)
  if ( max(Mphip,Mphin) > 1 ) then
    call calculate_fields(dens_rhoLR,dens_kappaLR,dens_kappaRL, & 
                          field_gammaLR,field_hspLR,field_deltaLR, & 
                          field_deltaRL,ndim)
  else 
    call calculate_fields_diag(dens_rhoLR,dens_kappaLR, & 
                               field_gammaLR,field_hspLR,field_deltaLR, & 
                               field_deltaRL,ndim)
  endif 
!cmpi  if ( nangle_max == 0 ) exit
!cmpi  if ( paral_myteamrank > 0 ) cycle 

  if ( iopt+iprint > 0 ) then
    call generate_rotation_parity(1,zone*ROTG,ROTG,ndim)
    call calculate_overlap(nocc0,nemp0,ovac0,voveru0,bogo_zD0, &
                           nocc0,nemp0,ovac0,voveru0,bogo_zD0, &
                           ROTG,rot_pari,ndim)
  endif

  !!! Particle number (always computed)
  call calculate_particle_number(0,dens_rhoLR,dens_kappaLR,dens_kappaRL, & 
                                 rot_prot,rot_neut,rot_prot2,rot_neut2,ndim)

  !!! Computes all other observables 
  !!! BB: this should be simplified/compacted at some point 
  if ( iopt == 0 ) then 
    !!! Energy and H20 for gradient
    call calculate_expectval_energy(dens_rhoLR,dens_kappaRL,field_gammaLR, &
                                    field_deltaLR,rot_ener,ndim)
    call calculate_H20(bogo_zU0tilde,bogo_zV0tilde,rot_H20,ndim)
  else 
    !!! Energy (decomposition)
    call calculate_decompo_energy(dens_rhoLR,dens_kappaRL,field_gammaLR, &
                                  field_deltaLR,rot_ecomp,ndim)
    
    !!! Occupation numbers
    call calculate_occupation_number(dens_rhoLR,rot_occnum,ndim,HOsh_dim)
    
    !!! Jx, Jy, Jz
    call calculate_expectval_obo_cplx(dens_rhoLR,zone*angumome_Jx,amjx_p, &
                                      amjx_n,ndim)
    call calculate_expectval_obo_cplx(dens_rhoLR,zimag*angumome_Jy,amjy_p, &
                                      amjy_n,ndim) 
    call calculate_expectval_obo_cplx(dens_rhoLR,zone*angumome_Jz,amjz_p, &
                                      amjz_n,ndim)
    rot_amj(1) = amjx_p + amjx_n
    rot_amj(2) = amjy_p + amjy_n
    rot_amj(3) = amjz_p + amjz_n
    
    !!! Jx^2, Jy^2, Jz^2
    call calculate_expectval_obos_cplx(dens_rhoLR,dens_kappaLR,dens_kappaRL, &
                                       zone*angumome_Jx,zone*angumome_Jx2, &
                                       rot_amj2(1),ndim)
    call calculate_expectval_obos_cplx(dens_rhoLR,dens_kappaLR,dens_kappaRL, &
                                       -zimag*angumome_Jy,-zone*angumome_Jy2,&
                                       rot_amj2(2),ndim)
    call calculate_expectval_obos_cplx(dens_rhoLR,dens_kappaLR,dens_kappaRL, &
                                       zone*angumome_Jz,zone*angumome_Jz2, &
                                       rot_amj2(3),ndim)

    !!! spin orbit
    call calculate_expectval_obo_cplx(dens_rhoLR,zone*angumome_so,rot_spor(1), &
                                      rot_spor(2),ndim)
    
    !!! Multipole Qlm
    do i = 1, 4
      do j = 0, i
    
        if ( i == 1 ) then 
          Qlm(:) = zone * (multipole_Q1m(:,1,j) + multipole_Q1m(:,2,j))
        elseif ( i == 2 ) then 
          Qlm(:) = zone * (multipole_Q2m(:,1,j) + multipole_Q2m(:,2,j))
        elseif ( i == 3 ) then 
          Qlm(:) = zone * (multipole_Q3m(:,1,j) + multipole_Q3m(:,2,j))
        else
          Qlm(:) = zone * (multipole_Q4m(:,1,j) + multipole_Q4m(:,2,j))
        endif 
    
        call calculate_expectval_obo_cplx(dens_rhoLR,Qlm,rot_Qlm(1,j,i), &
                                          rot_Qlm(2,j,i),ndim)
      enddo
    enddo
  
    !!! radius square      
    rad2(:) = zone * (radius_r2(:,1) + radius_r2(:,2))
    call calculate_expectval_obo_cplx(dens_rhoLR,rad2,rot_rad2(1),rot_rad2(2), &
                                      ndim)
  endif
  
  !!! Sums this gauge angle constribution to the integral
  call sum_gauge(iopt,weip,wein)

enddo  !end loop nangle

!!! In MPI runs: communicates all quantites through allreduce
!cmpi if ( (paral_teams > 1) .and. (paral_myteamrank == 0) ) then        
!cmpi   call reduce_projected_quantities(iopt)
!cmpi endif

!cmpi if ( paral_myrank == 0 ) then        
if ( iopt == 0 ) call calculate_projected_H20(ndim)
!cmpi endif

deallocate(voveru0)
deallocate(proj_phip,proj_weip,proj_phin,proj_wein)

end subroutine project_wavefunction

!------------------------------------------------------------------------------!
! subroutine reset_pnp                                                         !
!                                                                              ! 
! Sets the values of the particle-number projected quantities before a new     ! 
! call of the projection routine                                               ! 
!                                                                              ! 
! Input: iopt = 0 sets to zero only quantities used during the iterative proc. ! 
!             = 1  "   "   "   all       "                                     ! 
!------------------------------------------------------------------------------!
subroutine reset_pnp(iopt)

integer, intent(in) :: iopt

pnp_over = zzero  
pnp_ener = zzero  
pnp_prot = zzero  
pnp_neut = zzero  
pnp_pari = zzero  
pnp_prot2 = zzero
pnp_neut2 = zzero

if ( iopt == 0 ) then
  pnp_A = zzero
  pnp_H20 = zzero
  pnp_AH20 = zzero
else
  pnp_rad2 = zzero  
  pnp_amj = zzero
  pnp_amj2 = zzero
  pnp_spor = zzero  
  pnp_Qlm = zzero
  pnp_occnum = zzero
  pnp_ecomp = zzero
  pnp_rho = zzero
endif

end subroutine reset_pnp  

!------------------------------------------------------------------------------!
! subroutine sum_gauge                                                         !
!                                                                              ! 
! Adds the weighted gauge rotated quantities to the global sum that will rep-  !
! reseent the (discretized) integral.                                          !
!                                                                              ! 
! Input: iopt = option to consider more or less quantitites                    ! 
!             = 0 only the most important                                      ! 
!             = 1 more quantities (used when printing results)                 ! 
!        weip = weight in the Fomenko discretization for PNP for protons       ! 
!        wein =   "    "   "     "           "        "   "   "  neutrons      ! 
!------------------------------------------------------------------------------!
subroutine sum_gauge(iopt,weip,wein)

integer, intent(in) :: iopt
complex(r64), intent(in) :: weip, wein
complex(r64) :: weight, factor1, factor2

weight = weip * wein
factor1 = weight * rot_over
factor2 = factor1 * rot_ener

pnp_over = pnp_over + factor1
pnp_ener = pnp_ener + rot_ener * factor1
pnp_prot = pnp_prot + rot_prot * factor1
pnp_neut = pnp_neut + rot_neut * factor1
pnp_pari = pnp_pari + rot_pari * weight
pnp_prot2 = pnp_prot2 + rot_prot2 * factor1
pnp_neut2 = pnp_neut2 + rot_neut2 * factor1

if ( iopt == 0 ) then
  pnp_A = pnp_A + rot_A * factor1     
  pnp_H20 = pnp_H20 + rot_H20 * factor1
  pnp_AH20 = pnp_AH20 + rot_A * factor2
else
  pnp_rad2 = pnp_rad2 + rot_rad2 * factor1
  pnp_amj = pnp_amj + rot_amj * factor1
  pnp_amj2 = pnp_amj2 + rot_amj2 * factor1
  pnp_spor = pnp_spor + rot_spor * factor1
  pnp_Qlm = pnp_Qlm + rot_Qlm * factor1
  pnp_occnum = pnp_occnum + rot_occnum * factor1
  pnp_ecomp = pnp_ecomp + rot_ecomp * factor1
  pnp_rho = pnp_rho + dens_rhoLR * factor1
endif

end subroutine sum_gauge                  

!------------------------------------------------------------------------------!
! subroutine reduce_projected_quantities                                       !
!                                                                              !
! Reduces all important projected quantities using ALLREDUCE such that in the  !
! end all processes have the same values (particularly important because they  !
! need to have the same gradient/evolution).                                   !
!                                                                              !
! Input: iopt = option to consider more or less quantitites                    ! 
!------------------------------------------------------------------------------!
!cmpi subroutine reduce_projected_quantities(iopt)

!cmpi integer, intent(in) :: iopt
!cmpi integer :: ndim2, ndim3, ierr=0
!cmpi complex(r64) :: pnp_over_red, pnp_ener_red, pnp_pari_red, & 
!cmpi                 pnp_rad2_red(2), pnp_prot_red, pnp_neut_red, &
!cmpi                 pnp_prot2_red, pnp_neut2_red, &
!cmpi                 pnp_amj_red(3), pnp_amj2_red(3), pnp_spor_red(2), &
!cmpi                 pnp_Qlm_red(2,0:4,4) 
!cmpi complex(r64), dimension(12) :: pnp_ecomp_red
!cmpi complex(r64), dimension(:,:), allocatable :: pnp_occnum_red, pnp_H20_red, &
!cmpi                                              pnp_A_red, pnp_AH20_red, &
!cmpi                                              pnp_rho_red
!cmpi
!cmpi  call mpi_reduce(pnp_over,pnp_over_red,1,mpi_double_complex, &
!cmpi                  mpi_sum,0,mpi_comm_peers,ierr)
!cmpi  call mpi_reduce(pnp_ener,pnp_ener_red,1,mpi_double_complex, &
!cmpi                  mpi_sum,0,mpi_comm_peers,ierr)
!cmpi  call mpi_reduce(pnp_prot,pnp_prot_red,1,mpi_double_complex, &
!cmpi                  mpi_sum,0,mpi_comm_peers,ierr)
!cmpi  call mpi_reduce(pnp_neut,pnp_neut_red,1,mpi_double_complex, &
!cmpi                  mpi_sum,0,mpi_comm_peers,ierr)
!cmpi  call mpi_reduce(pnp_pari,pnp_pari_red,1,mpi_double_complex, &
!cmpi                  mpi_sum,0,mpi_comm_peers,ierr)
!cmpi  call mpi_reduce(pnp_prot2,pnp_prot2_red,1,mpi_double_complex, &
!cmpi                  mpi_sum,0,mpi_comm_peers,ierr)
!cmpi  call mpi_reduce(pnp_neut2,pnp_neut2_red,1,mpi_double_complex, &
!cmpi                  mpi_sum,0,mpi_comm_peers,ierr)
!cmpi  pnp_over = pnp_over_red
!cmpi  pnp_ener = pnp_ener_red
!cmpi  pnp_prot = pnp_prot_red
!cmpi  pnp_neut = pnp_neut_red
!cmpi  pnp_pari = pnp_pari_red
!cmpi  pnp_prot2 = pnp_prot2_red
!cmpi  pnp_neut2 = pnp_neut2_red
!cmpi 
!cmpi  if ( iopt == 0 ) then
!cmpi    ndim2 = HOsp_dim**2
!cmpi    pnp_A_red = pnp_A
!cmpi    pnp_H20_red = pnp_H20
!cmpi    pnp_AH20_red = pnp_AH20
!cmpi    call mpi_reduce(pnp_A,pnp_A_red,ndim2,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_H20,pnp_H20_red,ndim2,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_AH20,pnp_AH20_red,ndim2,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    pnp_A = pnp_A_red
!cmpi    pnp_H20 = pnp_H20_red
!cmpi    pnp_AH20 = pnp_AH20_red
!cmpi  else
!cmpi    ndim2 = 2 * HOsh_dim
!cmpi    ndim3 = HOsh_dim**2
!cmpi    pnp_occnum_red = zzero * pnp_occnum
!cmpi    pnp_rho_red = zzero * pnp_rho   
!cmpi    call mpi_reduce(pnp_rad2,pnp_rad2_red,2,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_amj,pnp_amj_red,3,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_amj2,pnp_amj2_red,3,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_spor,pnp_spor_red,2,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_Qlm,pnp_Qlm_red,40,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_occnum,pnp_occnum_red,ndim2,mpi_double_complex, & 
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_rho,pnp_rho_red,ndim3,mpi_double_complex, & 
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    call mpi_reduce(pnp_ecomp,pnp_ecomp_red,12,mpi_double_complex, &
!cmpi                    mpi_sum,0,mpi_comm_peers,ierr)
!cmpi    pnp_rad2 = pnp_rad2_red
!cmpi    pnp_amj = pnp_amj_red
!cmpi    pnp_amj2 = pnp_amj2_red
!cmpi    pnp_spor = pnp_spor_red
!cmpi    pnp_Qlm = pnp_Qlm_red
!cmpi    pnp_occnum = pnp_occnum_red
!cmpi    pnp_rho = pnp_rho_red
!cmpi    pnp_ecomp = pnp_ecomp_red
!cmpi  endif

!cmpi end subroutine reduce_projected_quantities

!------------------------------------------------------------------------------!
! subroutine generate_rotation_gauge                                           !
!                                                                              !
! Generates the transformation matrix for gauge rotations. In the HO basis,    !
! the gauge rotation reads                                                     !
!   < a | R(phi_p, phi_n) | b > =  e^(i phi_{mt_a}) delta_ab                   !
! where                                                                        !
!   | a > = | n_a l_a j_a m_a mt_a >                                           !
!   | b > = | n_b l_b j_b m_b mt_b >                                           !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        phip = gauge angle for protons                                        !
!        phin = gauge angle for neutrons                                       !
!        RI = initial rotation matrix                                          !
!                                                                              !
! Output: RF = final rotation matrix                                           !
!------------------------------------------------------------------------------!
subroutine generate_rotation_gauge(phi_p,phi_n,RI,RF,ndim)

integer, intent(in) :: ndim
real(r64), intent(in) :: phi_p, phi_n
complex(r64), dimension(ndim,ndim), intent(in) :: RI 
complex(r64), dimension(ndim,ndim), intent(out) :: RF
integer :: i, j
complex(r64) :: phase_p, phase_n 

phase_p = exp(zimag*phi_p) 
phase_n = exp(zimag*phi_n) 

do j = 1, ndim
  do i = 1, ndim/2
    RF(i,j) = phase_p * RI(i,j)
    RF(i+ndim/2,j) = phase_n * RI(i+ndim/2,j)
  enddo
enddo

end subroutine generate_rotation_gauge

!------------------------------------------------------------------------------!
! subroutine generate_rotation_parity                                          !
!                                                                              !
! Generates the transformation matrix for parity operation. In the HO basis,   !
! the parity operation reads                                                   !
!   < a | P | b > = (-1)**l_a delta_ab                                         !
! where                                                                        !
!   | a > = | n_a l_a j_a m_a mt_a >                                           !
!   | b > = | n_b l_b j_b m_b mt_b >                                           !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        ip = 0 for identity                                                   !
!           = 1 for P                                                          !
!        RI = initial rotation matrix                                          !
!                                                                              !
! Output: RF = final rotation matrix                                           !
!------------------------------------------------------------------------------!
subroutine generate_rotation_parity(ip,RI,RF,ndim)

integer, intent(in) :: ip, ndim
complex(r64), dimension(ndim,ndim), intent(in) :: RI 
complex(r64), dimension(ndim,ndim), intent(out) :: RF
integer :: i, j

!!! If no parity (ip=0), returns the same matrix
if ( ip == 0 ) then 
  RF = RI
else                
  do j = 1, ndim
    do i = 1, ndim
      RF(i,j) = (-1)**(HOsp_l(i)) * RI(i,j)
    enddo
  enddo
endif

end subroutine generate_rotation_parity

!------------------------------------------------------------------------------!
! subroutine rotate_wavefunction                                               !
!                                                                              !
! Computes the Bogoliubov matrices of the rotated state                        !
!   UF = R UI                                                                  !
!   VF = R^* VI                                                                !
! Note that here R is an arbitrary rotation matrix (e.g. Euler+Parity+Gauge)   !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        UI,VI = initial Bogoliubov matrices                                   !
!        R = rotation matrix                                                   !
!                                                                              !
! Output: UF,VF = final Bogoliubov matrices                                    ! 
!------------------------------------------------------------------------------!
subroutine rotate_wavefunction(R,UI,VI,UF,VF,ndim)        

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: R, UI, VI
complex(r64), dimension(ndim,ndim), intent(out) :: UF, VF
complex(r64), dimension(ndim,ndim) :: Rc

Rc = conjg(R)

call zgemm('n','n',ndim,ndim,ndim,zone, R,ndim,UI,ndim,zzero,UF,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,Rc,ndim,VI,ndim,zzero,VF,ndim)

end subroutine rotate_wavefunction 

!------------------------------------------------------------------------------!
! subroutine calculate_thouless                                                !
!                                                                              ! 
! Calculates the Thouless matrix relating the rotated state to the unrotated   ! 
! one following:                                                               ! 
!   Aphi = (Vphi Uphi^-1)^*                                                    ! 
! where                                                                        ! 
!   Uphi = U^dagger Ubar + V^dagger Vbar                                       !
!   Vphi = U^T Vbar + V^T Ubar                                                 !
!                                                                              ! 
! Also calculates the matrices                                                 ! 
!   Utilde = U^* + V Aphi                                                      ! 
!   Vtilde = V^* + U Aphi                                                      ! 
! used to comute the non-diagonal density matrices as defined in T. Rodriguez  ! 
! master tesis:                                                                ! 
!     rhoLR =   Vtilde V^T                                                     ! 
!   kappaLR =   Vtilde U^T                                                     ! 
!   kappaRL = - Utilde V^T                                                     ! 
!                                                                              ! 
! Note that as U and V are real, a lot of complex conjugates disappear in the  !
! computations below.                                                          !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        U,V = Bogoliubov matrices of the original state                       !
!        Ubar,Vbar = Bogoliubov matrices of the rotated state                  !
!                                                                              !
! Output: Aphi = Thouless matrices relating the two states                     !
!         Utilde,Vtilde = matrices used to compute non-diagonal densities      !
!------------------------------------------------------------------------------!
subroutine calculate_thouless(U,V,Ubar,Vbar,Aphi,Utilde,Vtilde,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: U, V, Ubar, Vbar
complex(r64), dimension(ndim,ndim), intent(out) :: Utilde, Vtilde, Aphi 
integer :: info1, info2
integer, dimension(ndim) :: ipiv
complex(r64), dimension(ndim) :: work
complex(r64), dimension(ndim,ndim) :: Uphi, Vphi, Uinv, A1, A2
  
!!! Uphi
call zgemm('t','n',ndim,ndim,ndim,zone,V,ndim,Vbar,ndim,zzero,Uphi,ndim)
call zgemm('t','n',ndim,ndim,ndim,zone,U,ndim,Ubar,ndim, zone,Uphi,ndim)

!!! Vphi
call zgemm('t','n',ndim,ndim,ndim,zone,V,ndim,Ubar,ndim,zzero,Vphi,ndim)
call zgemm('t','n',ndim,ndim,ndim,zone,U,ndim,Vbar,ndim, zone,Vphi,ndim)

!!! Uphi^-1
Uinv = Uphi
call zgetrf(ndim,ndim,Uinv,ndim,ipiv,info1)
call zgetri(ndim,Uinv,ndim,ipiv,work,ndim,info2)

!!! Aphi  
call zgemm('n','n',ndim,ndim,ndim,zone,Vphi,ndim,Uinv,ndim,zzero,Aphi,ndim)
Aphi = conjg(Aphi)

!!! Utilde and Vtilde     
call zgemm('n','n',ndim,ndim,ndim,zone,V,ndim,Aphi,ndim,zzero,A1,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,U,ndim,Aphi,ndim,zzero,A2,ndim)

Utilde = U + A1
Vtilde = V + A2
 
end subroutine calculate_thouless

!------------------------------------------------------------------------------!
! subroutine calculate_projected_H20                                           !
!                                                                              ! 
! Computes the projected H20 (vectorized) that is needed to compute the grad-  ! 
! ient. The formulae are taken from PNVAP equations in the master thesis of    ! 
! Tomás R. Rodriguez.                                                          ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!------------------------------------------------------------------------------!
subroutine calculate_projected_H20(ndim)

integer, intent(in) :: ndim
integer :: i, j, l
real(r64) :: eps=1.0d-8, imag_grad
complex(r64) :: over, ener
complex(r64), dimension(ndim,ndim) :: A, H20, AH20

over = pnp_over                                                            
ener = pnp_ener / pnp_over        

!!! Computes projected intermediate quantities
A    = 0.5d0 * (pnp_A    - transpose(pnp_A))    / over     
H20  = 0.5d0 * (pnp_H20  - transpose(pnp_H20))  / over     
AH20 = 0.5d0 * (pnp_AH20 - transpose(pnp_AH20)) / over     

!!! Computes H20
l = 0
do j = 1, ndim
  do i = 1, ndim
    l = l + 1
    field_H20v(l) = real(AH20(i,j) - ener*A(i,j) + H20(i,j))
    imag_grad = aimag((pnp_AH20(i,j) - ener*pnp_A(i,j) + pnp_H20(i,j)) / over)   
    if (imag_grad > eps) then
      print*, 'The gradient has non-zero imaginary part:', imag_grad
      stop 
    endif  
  enddo
enddo

end subroutine calculate_projected_H20

!------------------------------------------------------------------------------!
! subroutine calculate_overlap                                                 !
!                                                                              !
! Computes the overlap <0|R|1> for general Bogoliubov quasiparticle states     !
! |0> and |1> and where R is a general rotation matrix (e.g. spatial + gauge   !
! + parity).                                                                   !
! The overlap is computed through the pfaffian formula found in the reference  !
! Avez.2012.PhysRevC.85.034325.                                                !
! The routine to computes the pfaffian of a skew-symmetric matrix is taken     !
! from the reference Wimmer.2012.ACM.TransMathSoftware.38.30.                  !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        nocc0 = number of occupied states in the canonical basis of |0>       !
!        nemp0 = number of empty    states in the canonical basis of |0>       !
!        ovac0 = overlap of |0> with the single-particule vacuum               !
!        voveru0 = values of v/u in the canonical basis of |0>                 !
!        D0 = matrix D of the BMZ transformation for |0> (canon. basis)        !
!        nocc1, nemp1, ovac1, voveru1, D1 = same but for the state |1>         !
!                                                                              !
! Ouput: overlap = complex overlap <0|R|1>                                     !
!------------------------------------------------------------------------------!
subroutine calculate_overlap(nocc0,nemp0,ovac0,voveru0,D0, &
                             nocc1,nemp1,ovac1,voveru1,D1,ROT,overlap,ndim)

integer, intent(in) :: ndim, nocc0, nemp0, nocc1, nemp1 
real(r64), intent(in) :: ovac0, ovac1
real(r64), dimension(ndim-nocc0), intent(in) :: voveru0                 
real(r64), dimension(ndim-nocc1), intent(in) :: voveru1                 
complex(r64), dimension(ndim,ndim), intent(in) :: D0, D1, ROT
complex(r64), intent(out) :: overlap
integer :: i, j, nreg, nempm, nbloc0, nbloc1, n0, n1, n2, info, isla0, isla1, &
           ialloc=0
complex(r64) :: detR, normvac, over_module, over_phase
real(r64) :: sgn, fac1, fac2, fac3, signfac
complex(r64), dimension(ndim,ndim) :: DRc, Rn
complex(r64), dimension(:,:), allocatable :: Mreg, Rn2, Rinv
integer   , dimension(:), allocatable :: ipiv, iwork        
real(r64), dimension(:), allocatable :: rwork        
complex(r64), dimension(:), allocatable :: zwork, zwork2

!!! Computes the dimensions for the regularized matrix M, and computes the 
!!! phase factor accordingly 
nempm = min(nemp1,nemp0) 
nbloc0 = ndim - nocc0 - nempm
nbloc1 = ndim - nocc1 - nempm
nreg = nbloc1 + nbloc0 ! = 2*(n-nempm) - (nocc0+nocc1)

n0 = ndim - nocc0
n1 = ndim - nocc1
n2 = ndim - nempm
fac1 = (-1)**(n2*(n2+1)/2)
fac2 = (-1)**(nocc0*(nocc0-1)/2)
fac3 = (-1)**(nocc0*(nocc0+1)/2 + nocc1*(nocc1+1)/2 + nocc0*(n0+n2) + nocc1*n1)
signfac = fac1 * fac2 * fac3

!!! When constructing the canonical basis, we store ovac = log(prod u) such 
!!! that we use sqrt(a*b) = exp(1/2 [log(a) + log(b)])
normvac = zone * (0.5d0 * (ovac0 + ovac1))

!!! Check if we are dealing with Slater determinants
isla0 = 0
if ( n0-nemp0 == 0) isla0 = 1
isla1 = 0
if ( n1-nemp1 == 0) isla1 = 2

select case ( isla0+isla1 )
  case (1)
    if ( nocc0 > ndim-nemp1 ) then
      overlap = zzero
      return
    endif
  case (2)
    if ( nocc1 > ndim-nemp0 ) then
      overlap = zzero
      return
    endif
  case (3)
    if ( nocc0 == nocc1 ) then
      overlap = zone  
    else
      overlap = zzero
      return
    endif
  case default
    continue
end select

!!! Determine the overlap matrix R between the canonical bases of the original
!!! left state and the rotated right state. The inversion and determinant 
!!! calculations of R are done using LU factorization
allocate( Rn2(n2,n2), Rinv(n2,n2), ipiv(n2), zwork(n2), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation (1) in calculate overlap'

call zgemm('n','n',ndim,ndim,ndim,zone,ROT,ndim,D1,ndim,zzero,DRc,ndim)
call zgemm('c','n',ndim,ndim,ndim,zone,D0,ndim,DRc,ndim,zzero,Rn,ndim)

Rn2(1:n2,1:n2) = Rn(1+nempm:ndim,1+nempm:ndim) 
Rinv = Rn2

call zgetrf(n2,n2,Rinv,n2,ipiv,info)
if (info /= 0) then
  print*,'In calculate_overlap got info = ',info,' from zgetrf'
  stop 
endif

detR = zone
do i = 1, n2  
  detR = detR * Rinv(i,i)
enddo
 
sgn = one   
do i= 1, n2
  if (ipiv(i) /= i) then
    sgn = -sgn
  endif
enddo
detR = sgn * detR
         
call zgetri(n2,Rinv,n2,ipiv,zwork,n2,info)
if (info /= 0) then
  print*,'In calculate_overlap got info = ',info,' from zgetri'
  stop 
endif

!!! Construct the M matrix 
if ( isla0+isla1 /= 3 ) then
  allocate( Mreg(nreg,nreg), zwork2(nreg**2), iwork(nreg), rwork(nreg), &
           stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation (2) in calculate overlap'

  Mreg = zzero

  ! Upper left corner
  do i = 1, nbloc1
    if ( (-1)**i == -1 ) then 
      Mreg(i,i+1) = voveru1(i+nempm) 
    else
      Mreg(i,i-1) = voveru1(i+nempm) 
    endif
  enddo

  ! Lower right corner
  do j = 1, nbloc0
    i = j + nbloc1
    if ( (-1)**i == -1 ) then 
      Mreg(i,i+1) = -voveru0(j+nempm) 
    else
      Mreg(i,i-1) = -voveru0(j+nempm) 
    endif
  enddo

  ! Upper right and bottom left corners
  do j = 1, nbloc0
    do i = 1, nbloc1
      Mreg(i,j+nbloc1) = -Rinv(i,j) 
      Mreg(j+nbloc1,i) =  Rinv(i,j) 
    enddo
  enddo

  ! Computes the pfaffian using routines from M. Wimmer
  call zskpfa('u','p',nreg,Mreg,nreg,overlap,iwork,zwork2,nreg**2,rwork,info)

  deallocate(Mreg,zwork2,iwork,rwork)
endif 


!!! Final value of the overlap
if ( abs(overlap) > 0.0d0 ) then
  over_phase = zimag * atan2(aimag(overlap),real(overlap))
  over_module = zone * log(abs(overlap)) 

  overlap = detR * signfac * exp(over_phase + over_module + normvac)
endif

deallocate(Rinv,Rn2,ipiv,zwork)

end subroutine calculate_overlap

!------------------------------------------------------------------------------!
! subroutine calculate_norm                                                    !
!                                                                              ! 
! Computes the module of the overlap using the Onishi formula                  ! 
!   |<L|R>| = sqrt(|det(U)|) = sqrt(|det(UR^\dagger UL + VR^\dagger VL)|)      ! 
!                                                                              ! 
! Input: ndim = dimension of the sp basis                                      !
!        UL,VL = Bogoliubov matrices of the left  state                        ! 
!        UR,VR =      "         "    "   "  right   "                          ! 
!                                                                              ! 
! Output: norm = norm between the two states                                   ! 
!------------------------------------------------------------------------------!
subroutine calculate_norm(UL,VL,UR,VR,norm,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: UL, VL, UR, VR
real(r64), intent(out) :: norm 
integer :: i, info
integer, dimension(ndim) :: ipiv
complex(r64) :: detU 
complex(r64), dimension(ndim,ndim) :: UI       

!!! U
call zgemm('c','n',ndim,ndim,ndim,zone,UR,ndim,UL,ndim,zzero,UI,ndim)
call zgemm('c','n',ndim,ndim,ndim,zone,VR,ndim,VL,ndim, zone,UI,ndim)

call zgetrf(ndim,ndim,UI,ndim,ipiv,info)
if (info /= 0) then
  print*,' In calculate_norm got info = ',info,' from zgetrf'
  !stop 
endif

detU = zone
do i = 1, ndim
  detU = detU * UI(i,i)
enddo

norm = sqrt(abs(detU))

end subroutine calculate_norm

!------------------------------------------------------------------------------!
! subroutine print_results                                                     !
!                                                                              ! 
! Performs a last projection with the final state, computing more obsevables,  ! 
! and prints the results in the standard output. Also, occupation numbers are  ! 
! written in a file.                                                           ! 
!                                                                              ! 
! Input: Mphip = number of angles for the discretization of proton   PNR       ! 
!        Mphin =   "    "    "     "   "        "        "  neutrons  "        ! 
!------------------------------------------------------------------------------!
subroutine print_results(Mphip,Mphin)

integer, intent(in) :: Mphip, Mphin
integer :: i, j, ialloc=0
real(r64) :: over, pari, prot, neut, prot2, neut2,& 
             rad2_p, rad2_n, rad2_ch, rad2_is, rad2_iv, &
             amj(3), amj2(3), spor_p, spor_n, & 
             Qlm(4,0:4,4), betalm(4,0:4,4), betaT(4), &
             gammT(4), P_T00_J1m1, P_T00_J10, P_T00_J1p1, &
             P_T1m1_J00, P_T10_J00, P_T1p1_J00, ener_0b, ener_1b_p, ener_1b_n, & 
             ener_2bPH_pp, ener_2bPH_pn, ener_2bPH_np, ener_2bPH_nn,&
             ener_2bPP_pp, ener_2bPP_pn, ener_2bPP_np, ener_2bPP_nn, &
             energy_p, energy_n
complex(r64) :: fac
real(r64), dimension(:,:), allocatable, save :: occnum_pnp, occnum
complex(r64), dimension(:,:), allocatable, save :: rho_pnp
character(1) :: i_ch, j_ch
character(4) :: Qlm_ch
character(7) :: betalm_ch
character(18) :: label_ch  
character(19) :: filename  
character(len=*), parameter :: format1 = "(1a20,2f12.7)", &                              
                               format2 = "(1a9,39x,1f13.6)", &
                               format3 = "(1a9,2f13.6,13x,1f13.6)", &
                               format4 = "(1a9,4f13.6)", &
                               format5 = "(1a20,1a18,/)", &                              
                               format6 = "(1a4,4f12.6)", &
                               format7 = "(1a7,5f12.6)", &
                               format8 = "(1a5,3f12.6)", &
                               format9 = "(1a5,12x,1f12.6)", &
                               format10 = "(1a13,3f12.6)", &
                               format11 = "(1i3,1x,4i6,1i9,2f13.6)", &
                               format12 = "(1a5,2f12.6)", &
                               format13 = "(1a61)", &
                               format14 = "(1a64)", &
                               format15 = "(1a3,34x,2f13.6)"
!!!
!!! Initialization
!!!
if ( allocated(rho_pnp) .eqv. .false. ) then 
  allocate( rho_pnp(HOsp_dim,HOsp_dim), stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of projected density'
  rho_pnp = zzero
endif

if ( allocated(occnum_pnp) .eqv. .false. ) then 
  allocate( occnum_pnp(HOsp_dim,HOsp_dim), stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of occupation numbers'
  occnum_pnp = zero
endif

if ( allocated(occnum) .eqv. .false. ) then 
  allocate( occnum(HOsh_dim,2), stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of occupation numbers'
  occnum = zero
endif

!cmpi if ( paral_myrank == 0 ) then
if ( max(Mphip,Mphin) > 1 ) then
  print '(/,60("%"),/,17x,"PROJECTED STATE PROPERTIES",17x,/,60("%"),/)'
else
!!!
  print '(/,60("%"),/,15x,"QUASIPARTICLE STATE PROPERTIES",15x,/,60("%"),/)'
endif
!cmpi endif

!!!
!!! Projects the wave function
!!!
call project_wavefunction(1,0,Mphip,Mphin,HOsp_dim)

if ( max(Mphip,Mphin) > 1 ) then
  rho_pnp = pnp_rho / pnp_over
endif

!cmpi if ( paral_myrank == 0 ) then

!!!
!!! Basic properties 
!!!
over  = real( pnp_over )  ! If needed, the norm can be calculated via Onishi
pari  = real( pnp_pari / pnp_over )
prot  = real( pnp_prot / pnp_over )
neut  = real( pnp_neut / pnp_over )
prot2 = real( pnp_prot2 / pnp_over - (pnp_prot / pnp_over)**2 )
neut2 = real( pnp_neut2 / pnp_over - (pnp_neut / pnp_over)**2 )

if ( max(Mphip, Mphin) == 1 ) then
  write(label_ch,'(1i18)') bogo_label                                              
  label_ch = adjustl(label_ch)                                  
  write(uto,format5) 'Label of the state: ',label_ch                        
endif

print '(6x,"Quantity",12x,"Mean",6x,"Variance",/,44("-"))'                      
if ( max(Mphip, Mphin) > 1 ) then
  write(uto,format1) 'Projected overlap   ', over
else
  write(uto,format1) 'Norm                ', over
endif
write(uto,format1) 'Number of protons   ', prot, prot2
write(uto,format1) 'Number of neutrons  ', neut, neut2
write(uto,format1) 'Parity              ', pari

!!!
!!!  Energy
!!!
ener_0b      = real(hamil_H0)
ener_1b_p    = real( pnp_ecomp(1) / pnp_over ) 
ener_1b_n    = real( pnp_ecomp(2) / pnp_over ) 
ener_2bPH_pp = real( pnp_ecomp(3) / pnp_over ) 
ener_2bPH_pn = real( pnp_ecomp(4) / pnp_over ) 
ener_2bPH_np = real( pnp_ecomp(5) / pnp_over ) 
ener_2bPH_nn = real( pnp_ecomp(6) / pnp_over ) 
ener_2bPP_pp = real( pnp_ecomp(7) / pnp_over ) 
ener_2bPP_pn = real( pnp_ecomp(8) / pnp_over ) 
ener_2bPP_np = real( pnp_ecomp(9) / pnp_over ) 
ener_2bPP_nn = real( pnp_ecomp(10)/ pnp_over ) 
energy_p     = real( pnp_ecomp(11)/ pnp_over ) 
energy_n     = real( pnp_ecomp(12)/ pnp_over ) 

print '(/,"ENERGY DECOMPOSITION",/,20("="),//, &
        & "Part \ Iso",6x,"p-p",10x,"n-n",10x,"p-n",10x,"Total",/,61("-"))'
write(uto,format2) 'Zero-body', ener_0b
write(uto,format3) 'One-body ', ener_1b_p, ener_1b_n, (ener_1b_p + ener_1b_n)
write(uto,format4) ' ph part ', ener_2bPH_pp, ener_2bPH_nn, &
                     (ener_2bPH_pn + ener_2bPH_np), &
                     (ener_2bPH_pp + ener_2bPH_nn + ener_2bPH_pn + ener_2bPH_np)
write(uto,format4)  ' pp part ', ener_2bPP_pp, ener_2bPP_nn, & 
                     (ener_2bPP_pn + ener_2bPP_np), &
                     (ener_2bPP_pp + ener_2bPP_nn + ener_2bPP_pn + ener_2bPP_np)
write(uto,format4)  'Two-body ', (ener_2bPH_pp + ener_2bPP_pp), &
                  (ener_2bPH_nn + ener_2bPP_nn), &
                  (ener_2bPH_pn + ener_2bPP_pn + ener_2bPH_np + ener_2bPP_np), &
                  (ener_2bPH_pp + ener_2bPH_nn + ener_2bPH_pn + ener_2bPH_np + & 
                   ener_2bPP_pp + ener_2bPP_nn + ener_2bPP_pn + ener_2bPP_np)
write(uto,format4)  'Full H   ', (ener_1b_p + ener_2bPH_pp + ener_2bPP_pp), &
                        (ener_1b_n + ener_2bPH_nn + ener_2bPP_nn), &
                  (ener_2bPH_pn + ener_2bPP_pn + ener_2bPH_np + ener_2bPP_np), &
                  (ener_0b + ener_1b_p + ener_1b_n & 
                 + ener_2bPH_pp + ener_2bPH_nn + ener_2bPH_pn + ener_2bPH_np &
                 + ener_2bPP_pp + ener_2bPP_nn + ener_2bPP_pn + ener_2bPP_np)

!!!
!!! Multipole deformation
!!!
print '(/,"MULTIPOLE DEFORMATIONS",/,22("="),/, &
       & 37x,"Nucleons",/, &
       &"Q_lm",5x,"Protons",4x,"Neutrons",3x,"Isoscalar",3x,"Isovector", & 
       & /,52("-"))'

Qlm(1:2,:,:) = real( pnp_Qlm(1:2,:,:) / pnp_over )
Qlm(3,:,:) =  Qlm(1,:,:) + Qlm(2,:,:)
Qlm(4,:,:) = -Qlm(1,:,:) + Qlm(2,:,:)

do i = 1, 4
  write(i_ch,'(1i1)') i
  do j = 0, i
    write(j_ch,'(1i1)') j
    Qlm_ch = "Q_" // i_ch // j_ch
    write(uto,format6) Qlm_ch, Qlm(1,j,i), Qlm(2,j,i), Qlm(3,j,i), Qlm(4,j,i)
  enddo
enddo

!!! Beta parameters
do i = 1, 4
  do j = 0, i
    betalm(:,j,i) = Qlm(:,j,i) * coeff_betalm(:,i)
  enddo
enddo

print '(/,40x,"Nucleons",/, & 
       & "Beta_lm",5x,"Protons",4x,"Neutrons",3x,"Isoscalar",3x,"Isovector", &
       & /,55("-"))'
do i = 1, 4
  write(i_ch,'(1i1)') i
  do j = 0, i
    write(j_ch,'(1i1)') j
    betalm_ch = "Beta_" // i_ch // j_ch 
    write(uto,format7) betalm_ch, betalm(1,j,i), betalm(2,j,i), betalm(3,j,i), &
                       betalm(4,j,i)
  enddo
enddo

!!! Triaxial parameters
betaT(:) = sqrt( Qlm(:,0,2)**2 + 2.0d0*(Qlm(:,2,2)**2) ) * coeff_betalm(:,2)
gammT(:) = atan( sqrt(2.0d0) * abs( Qlm(:,2,2) / Qlm(:,0,2) ) )

where ( abs(Qlm(:,0,2)) <= epsilon0 ) gammT = 0.d0
where ( (Qlm(:,0,2) > 0.d0) .and. (Qlm(:,2,2) <  0.d0) ) gammT = 2.d0*pi - gammT
where ( (Qlm(:,0,2) < 0.d0) .and. (Qlm(:,2,2) >= 0.d0) ) gammT = pi - gammT
where ( (Qlm(:,0,2) < 0.d0) .and. (Qlm(:,2,2) <  0.d0) ) gammT = pi + gammT
gammT = gammT * 180.0/pi

print '(/,40x,"Nucleons",/,"Triaxial",4x,"Protons",4x,"Neutrons",3x, &
       & "Isoscalar",3x,"Isovector",/,55("-"))'
write(uto,format7) 'Beta   ', betaT(1), betaT(2), betaT(3), betaT(4)
write(uto,format7) 'Gamma  ', gammT(1), gammT(2), gammT(3), gammT(4)

!!!
!!! Radius RMS
!!!
rad2_p = real( pnp_rad2(1) / pnp_prot )
rad2_n = real( pnp_rad2(2) / pnp_neut )

rad2_is = real( ( pnp_rad2(1) + pnp_rad2(2)) / (pnp_prot + pnp_neut) )
rad2_iv = real( (-pnp_rad2(1) + pnp_rad2(2)) / (pnp_prot + pnp_neut) )

spor_p = real( pnp_spor(1) / pnp_over ) 
spor_n = real( pnp_spor(2) / pnp_over ) 

! RMS charge radius with corrections
rad2_ch = rad2_p + radius_r2p + (neut / prot) * radius_r2n & 
          + 0.75d0 * (hbarc / mass_ma)**2 &
          + (1 / prot) * ((hbarc / mass_ma)**2) * ( magmome_mup * spor_p + &
                                                    magmome_mun * spor_n ) 

print '(/,"RADIUS",/,6("="),/, &
       & 40x,"Nucleons",/, &
       & "Quantity",4x,"Protons",4x,"Neutrons",3x,"Isoscalar",3x,"Isovector", &
       & 5x,"Charge"/,67("-"))'
write(uto,format7) '  r    ', sqrt(rad2_p), sqrt(rad2_n), sqrt(rad2_is), &
                              sign(sqrt(abs(rad2_iv)),rad2_iv), sqrt(rad2_ch)
write(uto,format7) '  r^2  ', rad2_p, rad2_n, rad2_is, rad2_iv, rad2_ch

print '(/,"Warning: no center-of-mass correction.")'

!!!
!!! Angular momentum
!!!
amj  = real( pnp_amj  / pnp_over )
amj2 = real( pnp_amj2 / pnp_over )

print '(/,"ANGULAR MOMENTUM",/,16("="),//, &
      & 2x,"i",8x,"J_i",9x,"J_i^2",5x,"Var(J_i)",/,41("-"))'
write(uto,format8) '  X  ', amj(1), amj2(1), amj2(1) - amj(1)**2
write(uto,format8) '  Y  ', amj(2), amj2(2), amj2(2) - amj(2)**2
write(uto,format8) '  Z  ', amj(3), amj2(3), amj2(3) - amj(3)**2
write(uto,format9) 'Total', amj2(1) + amj2(2) + amj2(3)

!!!
!!! Pair coupling
!!!
if ( max(Mphip,Mphin) == 1 ) then
  call calculate_expectval_pair(dens_kappaRR,pairs_T00_J1m1,P_T00_J1m1,HOsp_dim)
  call calculate_expectval_pair(dens_kappaRR,pairs_T00_J10, P_T00_J10, HOsp_dim)
  call calculate_expectval_pair(dens_kappaRR,pairs_T00_J1p1,P_T00_J1p1,HOsp_dim)
  call calculate_expectval_pair(dens_kappaRR,pairs_T1m1_J00,P_T1m1_J00,HOsp_dim)
  call calculate_expectval_pair(dens_kappaRR,pairs_T10_J00, P_T10_J00, HOsp_dim)
  call calculate_expectval_pair(dens_kappaRR,pairs_T1p1_J00,P_T1p1_J00,HOsp_dim)

  print '(/,"PAIR COUPLING",/,13("="),//, & 
        & 3x,"MJ or MT =",7x,"-1",10x," 0",10x,"+1",/,49("-"))'
  write(uto,format10) 'T = 0 ; J = 1', abs(P_T00_J1m1), abs(P_T00_J10), & 
                                       abs(P_T00_J1p1)
  write(uto,format10) 'T = 1 ; J = 0', abs(P_T1m1_J00), abs(P_T10_J00), &
                                       abs(P_T1p1_J00)
endif

!!!
!!! Occupation numbers            
!!!
if ( max(Mphip,Mphin) > 1 ) then 
  occnum_pnp = real( pnp_occnum / pnp_over )
else 
  occnum = real( pnp_occnum / pnp_over )
endif

if ( max(Mphip,Mphin) == 1 ) then
  open(ute, file='occupation_numbers.dat', status='replace', action='write', &       
            form='formatted')                                                     
  write(ute,format13) "Occupation numbers"        
  write(ute,format14) "  #    2*mt    n     l    2*j   label   unprojected &
                    &   projected"

  do i = 1, 2      
    write(ute,'(64("-"))')
    do j = 1, HOsh_dim
      write(ute,format11) j, (-1)**i, HOsh_n(j), HOsh_l(j), HOsh_2j(j), & 
                          HOsh_na(j), occnum(j,i), occnum_pnp(j,i) 
    enddo
    write(ute,format15) "sum", sum(occnum(:,i)), sum(occnum_pnp(:,i))
  enddo

  close(ute, status='keep')                                                       
endif

!!!
!!! Spatial one-body density      
!!!
if ( max(proj_Mphip,proj_Mphin) > 1 ) then
  fac = zone 
else
  fac = zzero
endif

if ( (dens_spatial > 0) .and. (max(Mphip,Mphin) == 1) ) then
  call calculate_spatial_density(zone*dens_rhoRR,fac*rho_pnp,HOsp_dim)
endif

!!!
!!! Complementary informations    
!!!
if ( max(Mphip,Mphin) > 1 ) return

select case (seed_text)
  case(0,3)
    filename = 'final_wf.bin'
  case(1,2)
    filename = 'final_wf.txt'
end select

print '(/,60("%"),/,20x,"COMPLEMENTARY FILES",21x,/,60("%"),/)'
print '(5x,"Description",17x,"File",/,44("-"))'
if ( dens_spatial == 1 ) then
  print*,"Spatial density int: spatial_density_R.dat"
elseif ( dens_spatial == 2 ) then
  print*,"Spatial density int: spatial_density_R.dat"
  print*,"Spatial density 3D : spatial_density_RThetaPhi.dat"
elseif ( dens_spatial == 3 ) then
  print*,"Spatial density 3D : spatial_density_XYZ.dat"
endif
print*,"Occupation numbers : occupation_numbers.dat"
print*,"Canonical basis    : canonicalbasis.dat"
print*,"Eigenbasis h       : eigenbasis_h.dat"
print*,"Eigenbasis H11     : eigenbasis_H11.dat"
print*,"Final wave function: ", trim(adjustl(filename))
if ( hamil_read == 0 ) then
  print*,"Reduced hamiltonian: ", hamil_fred
endif

!cmpi endif

end subroutine print_results

END MODULE Projection     
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
