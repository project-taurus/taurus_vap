!==============================================================================!
! MODULE Basis                                                                 !
!                                                                              !
! This module contains the variables and routines related to the Harmonic Osc- !
! illator model space.                                                         !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_basis                                                       !
! - subroutine print_basis                                                     !
! - function radial_function                                                   !
! - function radial_integral                                                   !
! - function radial_integral_even                                              !
! - function spherharmonic                                                     !
!==============================================================================!
MODULE Basis      

!cmpi use MPI            
!cmpi use Parallelization
use Constants
use MathMethods
use Nucleus, only: valence_Z, valence_N, valence_A, & 
                   nucleus_Z, nucleus_N, nucleus_A

implicit none
public

!!! HO parameters
real(r64) :: HO_hw, & ! hbar*omega 
             HO_b     ! b = sqrt(hbar/(m*omega)) 
                      !   = hbar*c / sqrt((m*c**2)*(hbar*omega)) 

!!! HO single-particle states
integer :: HOsp_dim       ! dimension of the basis
integer(i64) :: HOsp_dim2 ! dimension squared (64bits)
integer, dimension(:), allocatable :: HOsp_n,   & ! quantum number    n 
                                      HOsp_l,   & !     "      "      l 
                                      HOsp_2j,  & !     "      "    2*j 
                                      HOsp_2mj, & !     "      "   2*mj 
                                      HOsp_2mt, & !     "      "   2*mt 
                                      HOsp_sh, &  ! shell where the sp is  
                                      HOsp_tr     ! time-reversal of indices

!!! HO shells
integer :: HOsh_dim ! dimension of the basis
integer, dimension(:), allocatable :: HOsh_n,  & ! quantum number  n
                                      HOsh_l,  & !     "      "    l
                                      HOsh_2j, & !     "      "   2j
                                      HOsh_na    ! label of the shell

!!! Other HO quantities
integer :: HO_Nmax, & ! maximum value of N = 2n + l + 1 
           HO_lmax, & !    "      "   "  l                 
           HO_2jmax   !    "      "   "  j1+j2 

!!! Private routines
private :: print_basis

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_basis                                                         !
!                                                                              !
! Sets HO basis: number of levels, quantum numbers, oscillator parameters, ... !
! To determine the basis, the codes reads the shells (same number for protons  !
! and neutrons) that defines the model space in the main hamiltonian file.     !
!                                                                              !
! If hamil_type = 1,2  ANTOINE                                                 !
!                      format of shells: HOsh_na = 1000*n + 100*l + 2*j        !
!               >= 3   General                                                 !
!                      format:of shells: HOsh_na = 10000*n + 100*l + 2*j       !
!                      (factor 10000 because the values of l can be > 10)      !
!------------------------------------------------------------------------------!
subroutine set_basis    

integer, parameter :: max_columns=100, max_length=5000
integer :: i, j, k, htype, opt_hw=0, facn, shinc, mjinc, jmax, ialloc=0, itmp
real(r64) :: xtmp
real(r64), dimension(1:max_columns) :: columns
character(len=max_length) :: line

!!! Recovers the basis information from hamiltonian file 
rewind(uth)
read(uth,*) 
read(uth,*) htype

if ( (htype == 1) .or. (htype == 2) ) then
  backspace uth
  read(uth,*) htype, HOsh_dim
else
  read(uth,*) HOsh_dim
endif

allocate( HOsh_n(HOsh_dim), HOsh_l(HOsh_dim), HOsh_2j(HOsh_dim), &
          HOsh_na(HOsh_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of shells'

backspace uth

if ( (htype == 1) .or. (htype == 2) ) then
  read(uth,*) itmp, HOsh_dim, (HOsh_na(i),i=1,HOsh_dim)

  !!! Small algorithm that counts the number of columns to see if there is
  !!! an optional value of HO_hw in ANTOINE files (not native)
  do i = 1, htype
    read(uth,*) 
  enddo

  read(uth,'(a)') line
  do i = 1, max_columns
    read(line,*,iostat=j) columns(1:i)
    if ( j == -1 ) exit
  enddo
  backspace uth
 
  if ( i-1 == 4 ) then
    read(uth,*)
    opt_hw = 0  
  else
    read(uth,*) itmp, itmp, itmp, xtmp, opt_hw, HO_hw
  endif 
else
  read(uth,*) HOsh_dim, (HOsh_na(i),i=1,HOsh_dim)
  read(uth,*) 
  read(uth,*) opt_hw, HO_hw
endif

!!! Computes oscillator parameters
select case (opt_hw)
  !!! Formula from Blomqvist et al., NPA 106, 545 (1968). 
  !!! See also formula (3.45) in the book by J. Suhonen
  case (0)
    HO_hw = 45.0d0 * nucleus_A**(-1.0d0/3.0d0) &
            - 25.0d0 * nucleus_A**(-2.0d0/3.0d0)
  
  !!! Formula from Bohr and Mottelson. 
  !!! See also (3.44) in the book by J. Suhonen
  case (1)
    HO_hw = 41.0d0 * nucleus_A**(-1.0d0/3.0d0) 

  !!! Takes the value read in the .sho file
  case (2) 
    continue

  !!! Default case: stops the code
  case default
    print*, "The option for the oscillator frequency, opt_hw =", opt_hw, &
            " should be equal to 0, 1 or 2."
    stop
end select
 
if ( HO_hw < epsilon0 ) then
  print "(a,1es12.5,a)", "The oscillator frequency (HO_hw) = ",HO_hw, &
        " should positive."
  stop
endif
 
HO_b = hbarc / sqrt(mass_ma * HO_hw)

!!! Determines the quantum numbers of the shells
if ( (htype == 1) .or. (htype == 2) ) then
  facn = 1000
else
  facn = 10000
endif

do i = 1, HOsh_dim
  HOsh_n(i)  =  HOsh_na(i) / facn 
  HOsh_l(i)  = (HOsh_na(i) - HOsh_n(i)*facn) / 100
  HOsh_2j(i) =  HOsh_na(i) - HOsh_n(i)*facn - HOsh_l(i)*100
enddo

!!! Computes the maximum values reachable in the basis
HO_2jmax = maxval(HOsh_2j)
HO_lmax  = maxval(HOsh_l)

HO_Nmax = 0
do i = 1, HOsh_dim
  j = 2*HOsh_n(i) + HOsh_l(i)
  HO_Nmax = max(j,HO_Nmax)
enddo

!!! Determines the dimension of the basis and check against the particle numbers
HOsp_dim = 0
do i = 1, HOsh_dim
  HOsp_dim = HOsp_dim + HOsh_2j(i) + 1
enddo
HOsp_dim = 2 * HOsp_dim 
HOsp_dim2 = HOsp_dim**2

if ( valence_Z > HOsp_dim/2 ) then 
  print*, 'The number of valence protons = ',valence_Z,' should be less', &
          ' than or equal to the number of s.p. states in the basis =', &
          HOsp_dim/2
  stop 
endif
if ( valence_N > HOsp_dim/2 ) then 
  print*, 'The number of valence neutrons = ',valence_N,' should be less', &
          ' than or equal to the number of s.p. states ls in the basis =', &
          HOsp_dim/2
  stop 
endif

!!! Determines the quantum numbers of the single-particle states
allocate( HOsp_n(HOsp_dim), HOsp_l(HOsp_dim), HOsp_2j(HOsp_dim), &
          HOsp_2mj(HOsp_dim), HOsp_2mt(HOsp_dim), HOsp_sh(HOsp_dim), &
          HOsp_tr(HOsp_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of single-particle states'

k = 0
shinc = 1
do i = 1, HOsh_dim
  jmax = HOsh_2j(i) + 1
  mjinc = 0
  do j = 1, jmax
    k = k + 1

    !!! protons
    HOsp_n(k)   = HOsh_n(i)
    HOsp_l(k)   = HOsh_l(i)
    HOsp_2j(k)  = HOsh_2j(i)
    HOsp_2mj(k) = HOsh_2j(i) - mjinc
    HOsp_2mt(k) = -1
    HOsp_sh(k)  = shinc   
    HOsp_tr(k)  = k + HOsp_2mj(k)

    !!! neutrons
    HOsp_n(k+HOsp_dim/2)   = HOsh_n(i)
    HOsp_l(k+HOsp_dim/2)   = HOsh_l(i)
    HOsp_2j(k+HOsp_dim/2)  = HOsh_2j(i)
    HOsp_2mj(k+HOsp_dim/2) = HOsh_2j(i) - mjinc
    HOsp_2mt(k+HOsp_dim/2) = 1
    HOsp_sh(k+HOsp_dim/2)  = shinc   
    HOsp_tr(k+HOsp_dim/2)  = k + HOsp_2mj(k) + HOsp_dim/2

    mjinc = mjinc + 2
  enddo
  shinc = shinc + 1
enddo

!!! Prints the basis informations in the standard output
!cmpi if ( paral_myrank == 0 ) then        
call print_basis
!cmpi endif

end subroutine set_basis    

!------------------------------------------------------------------------------!
! subroutine print_basis                                                       !
!                                                                              !
! Prints the characteristics of the basis.                                     !
!------------------------------------------------------------------------------!
subroutine print_basis       

integer :: i, j, nline
character(len=*), parameter :: format1 = "(1a18,2x,1i6)",   &
                               format2 = "(1a18,1x,1f7.3)", &
                               format3 = "(1a18,1x,8i7)"

nline = HOsh_dim/5 + 1
if ( mod(HOsh_dim,5) == 0 ) nline = nline - 1

print '(/,60("%"),/,26x,"HO BASIS",26x,/,60("%"),//, & 
      & 4x,"Quantity",10x,"Value",/,27("-"))'
print format1, 'No. of sp states  ', HOsp_dim
print format1, 'Max. value of N   ', HO_Nmax 
print format1, 'Max. value of l   ', HO_lmax  
print format1, 'Max. value of 2j  ', HO_2jmax
print format2, 'hbar*omega (MeV)  ', HO_hw  
print format2, 'Osc. length b (fm)', HO_b   
print format1, 'No. of shells     ', HOsh_dim
print format3, 'List of shells    ',(HOsh_na(i), i=1, min(5,HOsh_dim))
do i = 1, nline-1
  print format3, '                  ',(HOsh_na(j),j=1+i*5,min((i+1)*5,HOsh_dim))
enddo

end subroutine print_basis       

!------------------------------------------------------------------------------!
! function radial_function                                                     !
!                                                                              ! 
! Computes the radial function for the spherical HO single-particle states.    !
! R_nl(r) = Anl * (r/b)^2 * exp(-0.5*(r/b)**2) * L^(l+0.5)_n((r/b)**2)         !
!                                                                              !
! where                                                                        ! 
! A_nl = normalization factor                                                  ! 
! b = oscillator length                                                        ! 
! L^i_j = Laguerrre polynomial                                                 ! 
!------------------------------------------------------------------------------!
function radial_function(n,l,r) result(Rnl)

integer, intent(in) :: n, l
real(r64), intent(in) :: r
real(r64) :: Rnl, Anl, rob

!!! Checks the validity of the arguments
if ( (l < 0) .or. (n < 0) ) then 
  print "(a)", "Wrong argument(s) in function radial_function"
endif

!!! Computes the radial function
rob = r / HO_b
Anl = sqrt( (2**(n+l+2) * factorial(n)) / &
             (sqrt(pi) * dfactorial(2*n+2*l+1)) )

Rnl = Anl * (rob**l) * exp(-0.5d0 * rob**2) * genlaguerre(l+0.5d0,n,rob**2) / & 
      sqrt(HO_b)**3

end function radial_function

!------------------------------------------------------------------------------!
! function radial_integral                                                     !
!                                                                              ! 
! Computes the radial integral of r^lambda in the HO basis in the general case !
! (still la=lb) using numerical Gauss-Laguerre integration.                    !
!                                                                              ! 
! <a|r^lambda|b> = 0.5 * A_nala * A_nblb * b^lambda * \int du u^alpha e^-u *   ! 
!                  u^{lambda/2} L^alpha_na(u) L^alpha_nb(u)                    ! 
!                 |________________________________________| = f(u)            ! 
!                                                                              ! 
!                = 0.5 * A_nala * A_nblb * b^lambda * \sum_i w_i f(x_i)        ! 
!                                                                              ! 
! where                                                                        ! 
! A_nala, A_nblb = normalization factor                                        ! 
! b = oscillator length                                                        ! 
! L^i_j = Laguerrre polynomial                                                 ! 
! w_i, x_i = weights and abcissas of the Gauss-Laguerre quadrature.            ! 
! u = (r/b)^2                                                                  ! 
!------------------------------------------------------------------------------!
function radial_integral(a,b,lambda) result(integral)

integer, intent(in) :: a, b, lambda
integer :: na, nb, la, lb, nmax, np, i
real(r64) :: Anla, Anlb, alpha, integral
real(r64), dimension(:), allocatable :: xLag, wLag

na = HOsp_n(a)
nb = HOsp_n(b)
la = HOsp_l(a)
lb = HOsp_l(b) 

if ( la /= lb ) then
  integral = 0.d0
  return
endif

!!! Normalization factors 
Anla = sqrt( (2**(na+la+2) * factorial(na)) / &
              (sqrt(pi) * dfactorial(2*na+2*la+1)) )
Anlb = sqrt( (2**(nb+lb+2) * factorial(nb)) / & 
              (sqrt(pi) * dfactorial(2*nb+2*lb+1)) )

!!! Abscissas and weights of Gauss Laguerre
alpha = la + 0.5d0
nmax = max(na,nb)
np = nmax + lambda/2 + 10 !+10 just for safety
allocate(xLag(np))  
allocate(wLag(np))  
call GaussLaguerre(xLag,wLag,np,alpha)

!!! Integration
integral = 0.0d0
do i = 1, np
  integral = integral + wLag(i) * xLag(i)**(lambda/2) * &
                        genlaguerre(alpha,na,xLag(i)) * &
                        genlaguerre(alpha,nb,xLag(i))
enddo

integral = integral * Anla * Anlb * HO_b**lambda / 2.0d0

deallocate(xLag, wLag)

end function radial_integral

!------------------------------------------------------------------------------!
! function radial_integral_even                                                !
!                                                                              ! 
! Computes the radial integral of r^lambda in the HO basis in the case where   !
! lambda = even. The formula is the equation (6.41) in the book                !
! "From nucleons to nucleus" by J. Suhonen (ISBN:978-3-540-48859-0)            !
!------------------------------------------------------------------------------!
function radial_integral_even(a,b,lambda) result(integral)

integer, intent(in) :: a, b, lambda
integer :: na, nb, la, lb, sigma, sigma_max, sigma_min, taua, taub
real(r64) :: integral, prefactor, sum_sigma

na = HOsp_n(a) 
nb = HOsp_n(b) 
la = HOsp_l(a)
lb = HOsp_l(b)

integral = 0.0d0

if ( mod(la+lb+lambda,2) == 0 ) then
  taua = (lb-la+lambda)/2
  taub = (la-lb+lambda)/2
  if ( (taua >= 0) .and. (taub >= 0) ) then 
    prefactor = ((-1.d0)**(na+nb)) * sqrt( (factorial(na)*factorial(nb)) &
                  / (gamma(na+la+1.5d0)*gamma(nb+lb+1.5d0)) ) &
                * factorial(taua) * factorial(taub)

    sigma_min = max(0,na-taua,nb-taub)
    sigma_max = min(na,nb)
    sum_sigma = 0.d0
    do sigma = sigma_min, sigma_max
      sum_sigma = sum_sigma + gamma((la+lb+lambda)/2.d0+sigma+1.5d0) / &
                  ( factorial(sigma)*factorial(na-sigma)*factorial(nb-sigma) &
                    *factorial(sigma+taua-na)*factorial(sigma+taub-nb) )
    enddo
    integral = sum_sigma * prefactor * HO_b**lambda
  endif
endif

end function radial_integral_even

!------------------------------------------------------------------------------!
! function spherharmonic                                                       !
!                                                                              !
! Calculates the spherical harmonic Y_lm(theta,phi) using its relation to the  !
! associated Legendre polynomials:                                             !
! Y^m_l(theta,phi) = sqrt((2l+1) (l-m)! / 4pi (l+m)!) e^imphi P^_l(cos(theta)) !
!------------------------------------------------------------------------------!
function spherharmonic(l,m,theta,phi) result(Ylm)

integer, intent(in) :: l, m
real(r64), intent(in) :: theta, phi
real(r64) :: factor1, factor2, factor3
complex(r64) :: factor, Ylm

!!! Checks the validity of the arguments
if ( (l < 0) .or. (abs(m) > l) ) then 
  print "(a)", "Wrong argument(s) in function spherharmonic"
endif

!!! Computes the value of the spherical harmonic
factor1 = sqrt( (2*l+1) / (4*pi) )
factor2 = sqrt( one * factorial(l-m) )
factor3 = sqrt( one * factorial(l+m) )

factor = (factor1 * factor2 / factor3) * exp(zimag * m * phi)

Ylm = factor * assolegendre(l,m,cos(theta))

end function spherharmonic

END MODULE Basis      
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
