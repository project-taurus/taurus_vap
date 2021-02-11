!==============================================================================!
! MODULE Gradient                                                              !
!                                                                              !
! This module contains the variables and routines related to the gradient des- !
! cent (also the diagonalization of hsp and H11)                               !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_gradient                                                    !
! - subroutine calculate_gradient                                              !
! - subroutine diagonalize_hsp_and_H11                                         !
! - subroutine print_iteration                                                 !
!==============================================================================!
MODULE Gradient   

use Fields
use Constraints

implicit none
public

!!! Parameters controlling the gradient descent
integer :: gradient_type     ! type of gradient descent
real(r64) :: gradient_eta, & ! step parameter for gradient descent
             gradient_mu,  & ! weight of the momentum 
             gradient_eps    ! tolerance for gradient convergence

!!! Value of the gradient
real(r64):: gradient_norm ! frobenius norm of the gradient
real(r64), dimension(:), allocatable :: gradient_Zi, &  ! gradient iter.  i
                                        gradient_Zim1   !    "      "     i-1  

!!! Others
integer, private :: info_H11, & ! check if problem during diag(field_H11)
                    info_hsp    !   "   "     "      "    diag(field_hspRR)
real(r64), dimension(:), allocatable, private :: eigen_hsp, & ! sp energies 
                                                 eigen_H11    ! qp    "      
integer, dimension(:), private :: control_NZ(5)=0

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_gradient                                                      !
!                                                                              !
! Initializes the arrays related to the gradient                               !
!------------------------------------------------------------------------------!
subroutine set_gradient     

integer :: ialloc=0

allocate( gradient_Zi(HOsp_dim*HOsp_dim), gradient_Zim1(HOsp_dim*HOsp_dim), & 
          eigen_hsp(HOsp_dim), eigen_H11(HOsp_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of gradient'
   
gradient_Zi = zero
gradient_Zim1 = zero
gradient_norm = zero
eigen_hsp = zero
eigen_H11 = zero

end subroutine set_gradient

!------------------------------------------------------------------------------!
! subroutine calculate_gradient                                                !
!                                                                              !
! Calculates the gradient to construct the next step in the iterative minimiz- !
! ation process. This depdends on the option for the type of algorithm used:   !
! Standard Gradient or Heavy Ball (= gradient with momentum).                  !
!                                                                              !
! The direction in the descent is determined according to                      !
! Z(i) = - eta * (field_H20-constraints) + mu * Z(i-1)                         !
!                                                                              !
! The formulae for eta and mu are taken from Ryssens.2019.EPJA.55.9 and the    !
! approximation to the eigenvalues of the second derivative (stability matrix) !
! is taken from Robledo.2011.PhysRevC.84.014312. Be aware that all these reci- !
! pes are empirical.                                                           !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!------------------------------------------------------------------------------!
subroutine calculate_gradient(ndim)

integer, intent(in) :: ndim
integer :: i, j, k, opt_grad
real(r64) :: sum_grad, cond_numb, eigen_max, eigen_min, eta, mu 
real(r64), dimension(ndim,ndim) :: Zi2

opt_grad = gradient_type

!!! If something wrong happened during the diagonalization of H11 or hsp, the
!!! code uses eta and mu read from the input parameters                                 
if ( (info_H11 /= 0) .and. (gradient_type == 1) ) opt_grad = 0
if ( (info_hsp /= 0) .and. (gradient_type == 2) ) opt_grad = 0

!!! Computes the parameters for the evolution at this iteration
select case (opt_grad)

  !!! Parameters read from input file 
  case (0)
    eta = gradient_eta
    mu = gradient_mu

  !!! Heavy Ball with approximate optimal parameters + approximating the Hessian
  !!! à la Robledo                                                              
  case (1)
    eigen_max = 4.0d0 * maxval(abs(eigen_H11))   
    eigen_min = 2.0d0 * minval(abs(eigen_H11))
    cond_numb = eigen_max / eigen_min 
    eta = ( 2.0d0 / (sqrt(eigen_max) + sqrt(eigen_min)) )**2
    mu = ( (sqrt(cond_numb)-1) / (sqrt(cond_numb)+1) )**2

  !!! Heavy Ball with approximate optimal parameters + approximating the Hessian
  !!! eigenvalues by twice the single-particles energies (seems to work better                              
  !!! for odd-mass nuclei)                                                                              
  case (2)
    eigen_max = 2.0d0 * maxval(abs(eigen_hsp))
    eigen_min = 2.0d0 * minval(abs(eigen_hsp))
    cond_numb = eigen_max / eigen_min 
    eta = ( 2.0d0 / (sqrt(eigen_max) + sqrt(eigen_min)) )**2
    mu = ( (sqrt(cond_numb)-1) / (sqrt(cond_numb)+1) )**2

end select

!!! Computes the new gradient matrix Z
do i = 1, ndim**2
  sum_grad = zero 
  do j = 1, constraint_dim
    k = (j-1) * (ndim**2)
    sum_grad = sum_grad + lagrange_lambda1(j) * constraint_20(k+i)
  enddo
  gradient_Zi(i) = -eta * (field_H20v(i) - sum_grad) + mu * gradient_Zim1(i)
enddo

!!! Computes the Frobenius norm of Z to evaluate if the state is converged
call dgemm('t','n',ndim,ndim,ndim,one,gradient_Zi,ndim,gradient_Zi,ndim,zero, &
           Zi2,ndim)

sum_grad = zero
do i = 1, ndim
  sum_grad = sum_grad + Zi2(i,i)
enddo

gradient_norm = sqrt(sum_grad) / eta

end subroutine calculate_gradient

!------------------------------------------------------------------------------!
! subroutine diagonalize_hsp_and_H11                                           !
!                                                                              !
! Constructs the eigenvectors and eigenvalues of h and H11 in order to compute !
! the optimal minimization parameters (or so we hope) and writes in a file     !
! after the final iteration.                                                   !
!                                                                              !
! opt = 0 diagonalize what is needed for the gradient                          !
!     = 1 diagonlize but also prints the SP and QP basis in files              !
!                                                                              !
! Input: ndim = dimension of the sp basis                                      !
!        opt = option controling what is calculated                            !
!------------------------------------------------------------------------------!
subroutine diagonalize_hsp_and_H11(opt,ndim)

integer, intent(in) :: opt, ndim
integer :: i, j, k, l, m, nocc0, nemp0
integer, dimension(1) :: tabmin
integer, dimension(ndim) :: eigenh_order
real(r64), dimension(ndim) :: eigenh_tmp
real(r64), dimension(3*ndim-1) :: work
real(r64), dimension(ndim,ndim) :: D0, rhoc, hspc, A1, A2
real(r64) :: q11_aux, q00_aux, xn, xl2, xl, xneut, xprot, xpar, xjz, xj2, xj, &
             fermi_p, fermi_n, ovac0
real(r64), dimension(:,:), allocatable :: hspr
real(r64), dimension(:), allocatable :: workr, eigenr
complex(r64), dimension(ndim,ndim) :: hspRR, gammaRR, deltaRR
character(len=*), parameter :: format1 = "(1i4,7f9.3,1x,2f12.6)", &
                               format2 = "(1a77,/,80('-'))", &
                               format3 = "(1a89,/,92('-'))"

!!! Computes the fields
if ( (opt == 0) .and. (max(proj_Mphip,proj_Mphip) == 1) ) then 
   field_hspRR = real(field_hspLR)
   field_deltaRR = real(field_deltaRL)
else
  call calculate_fields_diag(zone*dens_rhoRR,zone*dens_kappaRR,gammaRR,hspRR, &
                             deltaRR,ndim=ndim)
  field_hspRR = real(hspRR)
  field_deltaRR = real(deltaRR)
endif

!cmpi if ( paral_myrank /= 0 ) return

call calculate_H11_real(ndim)

!!! Takes into account the constraints
if ( opt == 0 ) then 
  do m = 1, constraint_dim-constraint_pair
  k = (m-1)*(ndim**2) + 1
  if ( m <= constraint_dim-constraint_pair ) then
    call set_operator_qpbasis('f11',bogo_U0,bogo_V0,constraint_HO(k), &
                              constraint_11(k),ndim) 
  else
    call set_operator_qpbasis('g11',bogo_U0,bogo_V0,constraint_HO(k), &
                              constraint_11(k),ndim)
    endif
  enddo

  l = 0
  do j = 1, ndim
   do i = 1, ndim
      l = l + 1 
      q11_aux = zero
      q00_aux = zero 
      do m = 1, constraint_dim
        k = (m-1)*(ndim**2)
        q11_aux = q11_aux + lagrange_lambda1(m) * constraint_11(k+l)
        q00_aux = q00_aux + lagrange_lambda1(m) * constraint_HO(k+l)
      enddo
      field_H11(i,j) = field_H11(i,j) - q11_aux
      field_hspRR(i,j) = field_hspRR (i,j) - q00_aux
   enddo 
  enddo 
endif

!!! hsp in canonical basis
if ( opt == 1 ) then
  call construct_canonical_basis(bogo_U0,bogo_V0,bogo_zU0c,bogo_zV0c,bogo_zD0, &
                                 ovac0,nocc0,nemp0,ndim)
  D0 = real(bogo_zD0)

  call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,dens_rhoRR,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,rhoc,ndim)

  call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,field_hspRR,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,hspc,ndim)

  !!! Further reduces h in case of fully empty/occupides states
  if ( nemp0 > 0 ) then
    allocate (hspr(nemp0,nemp0), eigenr(nemp0),workr(3*nemp0-1))
    hspr(1:nemp0,1:nemp0) = hspc(1:nemp0,1:nemp0)
    call dsyev('v','u',nemp0,hspr,nemp0,eigenr,workr,3*nemp0-1,info_hsp)
    A1 = zero
    A2 = D0
    do i = 1, ndim
      A1(i,i) = one
    enddo
    A1(1:nemp0,1:nemp0) = hspr(1:nemp0,1:nemp0)
    call dgemm('n','n',ndim,ndim,ndim,one,A2,ndim,A1,ndim,zero,D0,ndim)
    deallocate(hspr, eigenr, workr)
  endif

  if ( nocc0 > 0 ) then
    allocate (hspr(nocc0,nocc0), eigenr(nocc0),workr(3*nocc0-1))
    hspr(1:nocc0,1:nocc0) = hspc(ndim-nocc0+1:ndim,ndim-nocc0+1:ndim)
    call dsyev('v','u',nocc0,hspr,nocc0,eigenr,workr,3*nocc0-1,info_hsp)
    A1 = zero
    A2 = D0
    do i = 1, ndim
      A1(i,i) = one
    enddo
    A1(ndim-nocc0+1:ndim,ndim-nocc0+1:ndim) = hspr(1:nocc0,1:nocc0)
    call dgemm('n','n',ndim,ndim,ndim,one,A2,ndim,A1,ndim,zero,D0,ndim)
    deallocate(hspr, eigenr, workr)
  endif

  call dgemm('t','n',ndim,ndim,ndim,one,D0,ndim,field_hspRR,ndim,zero,A1,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A1,ndim,D0,ndim,zero,hspc,ndim)

  !!! Ordering of energies
  l = 0
  eigenh_order = 0
  eigenh_tmp = 999

  do i = 1, ndim
    if ( abs(rhoc(i,i)) > 1.d-7 ) then
      l = l + 1
      eigenh_tmp(i) = hspc(i,i)
    endif
  enddo

  do i = 1, l
    tabmin = minloc(eigenh_tmp)
    eigenh_order(i) = tabmin(1)
    eigenh_tmp(tabmin(1)) = 1000
  enddo

  eigenh_tmp = 999

  do i = 1, ndim
    if ( abs(rhoc(i,i)) <= 1.d-7 ) then
      eigenh_tmp(i) = hspc(i,i)
    endif
  enddo

  do i = l+1, ndim
    tabmin = minloc(eigenh_tmp)
    eigenh_order(i) = tabmin(1)
    eigenh_tmp(tabmin(1)) = 1000
  enddo
endif

!!! Diagonalizes hsp                               
call dsyev('v','u',ndim,field_hspRR,ndim,eigen_hsp,work,3*ndim-1,info_hsp)

!!! Writes the properties of the single-particle states in a file
if ( opt == 1 ) then                                                             

  fermi_p = 0.d0
  fermi_n = 0.d0
 
  if ( constraint_switch(1) == 1 ) then 
    fermi_p = lagrange_lambda1(1)
  endif
  if ( constraint_switch(2) == 1 ) then 
    fermi_n = lagrange_lambda1(1 + constraint_switch(1))
  endif

  !!! Basis that diagonalizes h
  open(ute, file='eigenbasis_h.dat', status='replace', action='write', &       
           form='formatted')                                                     
  write(ute,"(1a,1f12.6)")   "Proton  fermi energy = ",fermi_p
  write(ute,"(1a,1f12.6,/)") "Neutron fermi energy = ",fermi_n
  write(ute,format2) "   #      Z        N        n        l        p &
                     &      jz        j         h  " 
  do i = 1, ndim                                                                 
    xneut = zero                                                                 
    xprot = zero                                                                 
    xpar  = zero                                                                 
    xjz   = zero                                                                 
    xj2   = zero                                                                 
    xn    = zero                                                                 
    xl2   = zero                                                                 
    do j = 1, ndim                                                               
      xprot = xprot + field_hspRR(j,i)**2 * (-HOsp_2mt(j) + 1)/2.0d0                
      xneut = xneut + field_hspRR(j,i)**2 * ( HOsp_2mt(j) + 1)/2.0d0                
      xpar  = xpar  + field_hspRR(j,i)**2 * (-1.d0)**HOsp_l(j)                     
      xn    = xn    + field_hspRR(j,i)**2 * HOsp_n(j)                              
      xjz   = xjz   + field_hspRR(j,i)**2 * HOsp_2mj(j)/2.0d0                       
      xj2   = xj2   + field_hspRR(j,i)**2 * (HOsp_2j(j)*(HOsp_2j(j)+2))/4.0d0       
      xl2   = xl2   + field_hspRR(j,i)**2 * (HOsp_l(j)*(HOsp_l(j)+1))             
    enddo                                                                        
    xj = 0.5d0 * (-1.d0 + sqrt(1+4*abs(xj2)))                                       
    xl = 0.5d0 * (-1.d0 + sqrt(1+4*abs(xl2)))                                       
    write(ute,format1) i, xprot, xneut, xn, xl, xpar, xjz, xj, eigen_hsp(i)           
  enddo                                                                          
  close(ute, status='keep')                                                       

  !!! Canonical basis
  open(ute, file='canonicalbasis.dat', status='replace', action='write', &       
           form='formatted')                                                     
  write(ute,"(1a,1f12.6)")   "Proton  fermi energy = ",fermi_p
  write(ute,"(1a,1f12.6,/)") "Neutron fermi energy = ",fermi_n
  write(ute,format3) "   #      Z        N        n        l        p &
                     &      jz        j         v2           h " 
  do i = 1, ndim                                                                 
    xneut = zero                                                                 
    xprot = zero                                                                 
    xpar  = zero                                                                 
    xjz   = zero                                                                 
    xj2   = zero                                                                 
    xn    = zero                                                                 
    xl2   = zero                                                                 
    k = eigenh_order(i)
    do j = 1, ndim                                                               
      xprot = xprot + D0(j,k)**2 * (-HOsp_2mt(j) + 1)/2.0d0                
      xneut = xneut + D0(j,k)**2 * ( HOsp_2mt(j) + 1)/2.0d0                
      xpar  = xpar  + D0(j,k)**2 * (-1.d0)**HOsp_l(j)                     
      xn    = xn    + D0(j,k)**2 * HOsp_n(j)                              
      xjz   = xjz   + D0(j,k)**2 * HOsp_2mj(j)/2.0d0                       
      xj2   = xj2   + D0(j,k)**2 * (HOsp_2j(j)*(HOsp_2j(j)+2))/4.0d0       
      xl2   = xl2   + D0(j,k)**2 * (HOsp_l(j)*(HOsp_l(j)+1))             
    enddo                                                                        
    xj = 0.5d0 * (-1.d0 + sqrt(1+4*abs(xj2)))                                       
    xl = 0.5d0 * (-1.d0 + sqrt(1+4*abs(xl2)))                                       
    write(ute,format1) i, xprot, xneut, xn, xl, xpar, xjz, xj, rhoc(k,k), & 
                       hspc(k,k)
  enddo                                                                          
  close(ute, status='keep')                                                       
endif                      

!!! Diaonalizes H11
call dsyev('v','u',ndim,field_H11,ndim,eigen_H11,work,3*ndim-1,info_H11)

!!! Writes the properties of the quasi-particle states in a file                                        
if ( opt == 1 ) then                                                             
  open(ute, file='eigenbasis_H11.dat', status='replace', action='write', &       
           form='formatted')                                                     
  write(ute,format2) "   #      Z        N        n        l        p & 
                     &      jz        j         H11"                                              
  do i = 1, ndim                                                                 
    xneut = zero                                                                 
    xprot = zero                                                                 
    xpar  = zero                                                                 
    xjz   = zero                                                                 
    xj2   = zero                                                                 
    xn    = zero                                                                 
    xl2   = zero                                                                 
    do j = 1, ndim                                                               
      xprot = xprot + field_H11(j,i)**2 * (-HOsp_2mt(j) + 1)/2.0d0                
      xneut = xneut + field_H11(j,i)**2 * ( HOsp_2mt(j) + 1)/2.0d0                
      xpar  = xpar  + field_H11(j,i)**2 * (-1.d0)**HOsp_l(j)                     
      xn    = xn    + field_H11(j,i)**2 * HOsp_n(j)                              
      xjz   = xjz   + field_H11(j,i)**2 * HOsp_2mj(j)/2.0d0                       
      xj2   = xj2   + field_H11(j,i)**2 * (HOsp_2j(j)*(HOsp_2j(j)+2))/4.0d0       
      xl2   = xl2   + field_H11(j,i)**2 * (HOsp_l(j)*(HOsp_l(j)+1))             
    enddo                                                                        
    xj = 0.5d0 * (-1.d0 + sqrt(1+4*abs(xj2)))                                       
    xl = 0.5d0 * (-1.d0 + sqrt(1+4*abs(xl2)))                                       
    write(ute,format1) i, xprot, xneut, xn, xl, xpar, xjz, xj, eigen_H11(i)           
  enddo                                                                          
  close(ute, status='keep')                                                       
endif                      

if ( (info_hsp+info_H11 /= 0) .and. (opt == 0) ) then
 print*,'Error in diagonalize_hsp_and_H11', info_hsp, info_H11
 !stop
endif 
 
end subroutine diagonalize_hsp_and_H11

!------------------------------------------------------------------------------!
! subroutine print_iteration                                                   !
!                                                                              !
! Prints selected projected quantities at the present iteration.               !
!                                                                              !
! Input: iter_print = option for printing at each iteration                    !
!        iter = iteration number                                               !
!------------------------------------------------------------------------------!
subroutine print_iteration(iter_print,iter)

use Projection, only : pnp_over, pnp_ener, pnp_pari, pnp_prot, pnp_neut, &
                       pnp_prot2, pnp_neut2 

integer, intent(in) :: iter_print, iter
integer :: i
real(r64) :: ener, pari, prot, neut, prot2, neut2, beta, gamm, &
             q20_p, q20_n, q20_a, q22_p, q22_n, q22_a
character(len=*), parameter :: format1 = "(1i6,5x,1es12.5,2x,5f12.6)", &
                               format2 = "(1i6,5x,1es12.5,2x,5f12.6,1f11.6, &  
                                         & 1f8.3,1f8.2)"

!!! Prints the "caption"
if ( iter == 1 ) then
  if ( iter_print == 0 ) then
    print '(" ")'
    print '("Iteration",5x,"Gradient",7x,"Energy",6x,"Protons",4x,"Var(Prot)", &
          & 4x,"Neutrons",3x,"Var(Neut)",/,85("-"))'
  elseif ( iter_print == 1 ) then
    print '(99x,"(unprojected)")'
    print '("Iteration",5x,"Gradient",7x,"Energy",6x,"Protons",4x,"Var(Prot)", &
          & 4x,"Neutrons",3x,"Var(Neut)",4x,"Parity",5x,"Beta",3x,"Gamma", &
          & /,112("-"))'
  endif
endif

!!! Normalization of projected value
ener = real(pnp_ener / pnp_over)
prot = real(pnp_prot / pnp_over)
neut = real(pnp_neut / pnp_over)
pari = real(pnp_pari / pnp_over)
prot2 = real(pnp_prot2 / pnp_over)
neut2 = real(pnp_neut2 / pnp_over)

prot2 = abs(prot2 - prot**2)
neut2 = abs(neut2 - neut**2)


if ( iter_print == 1 ) then
  call calculate_expectval_obo(dens_rhoRR,multipole_Q2m(1:HOsp_dim**2,1,0) + &
                          multipole_Q2m(1:HOsp_dim**2,2,0),q20_p,q20_n,HOsp_dim)
  call calculate_expectval_obo(dens_rhoRR,multipole_Q2m(1:HOsp_dim**2,1,2) + &
                          multipole_Q2m(1:HOsp_dim**2,2,2),q22_p,q22_n,HOsp_dim)
  q20_a = q20_p + q20_n
  q22_a = q22_p + q22_n

  beta = sqrt( q20_a**2 + 2.0d0*(q22_a**2) ) * coeff_betalm(2,3)
  gamm = atan( sqrt(2.d0) * abs(q22_a) / abs(q20_a) )

  if (  abs(q20_a) <= epsilon0 ) gamm = 0.d0
  if ( (q20_a > 0.d0) .and. (q22_a <  0.d0) )  gamm = 2.d0*pi - gamm
  if ( (q20_a < 0.d0) .and. (q22_a >= 0.d0) )  gamm = pi - gamm
  if ( (q20_a < 0.d0) .and. (q22_a <  0.d0) )  gamm = pi + gamm
  gamm = gamm * 180.0/pi
endif 

!!! Printing
if ( iter_print == 0 ) then
  print format1, iter, gradient_norm, ener, prot, prot2, neut, neut2
else
  print format2, iter, gradient_norm, ener, prot, prot2, neut, neut2, pari, &
        beta, gamm
endif

!!! Test on the average number of particles. Can stop the run. 
do i = 4, 1, -1
  control_NZ(i+1) = control_NZ(i)
enddo

control_NZ(1) = 0

if ( (abs(prot - valence_Z) > 0.5d0) .or. (abs(neut - valence_N) > 0.5d0) ) then
  control_NZ(1) = 1
  print '("Warning: the numbers of particles are far from the values set in &
         &the input parameters.")'
endif

if ( sum(control_NZ) >= 3 ) then
  print '("Critical error: the numbers of particles were wrong three times in &
         &in the last five iterations. The code will stop.")'
 stop
endif  

end subroutine print_iteration

END MODULE Gradient    
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
