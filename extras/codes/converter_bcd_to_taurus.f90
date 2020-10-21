!==============================================================================!
! PROGRAM converter_bcd_to_taurus                                              !
!                                                                              !
! This program transforms an interaction written in the format of the code     !
! BoccaDorata to an interaction written in the format of TAURUS.               !
! Both interactions are in J-scheme.                                           !
!                                                                              !
! To compile: gfortran/ifort -o bcd2taurus.exe converter_bcd_to_taurus.f90     !
! To use: ./bcd2taurus.exe < input.txt                                         !
!                                                                              !
! Formpat input file:     (Example)                                            !
! a11, a100                BCD file = my_bcd_file                              !
! a11, a100                TAU name = my_mainname_for_taurus_files             !
! a11, i2                  emax     = 6                                        !
! a11, f10.4               hw       = 20.0                                     !
! a11, i1                  Only 2b  = 1                                        !
! a11, i1                  Only com = 0                                        !
!                                                                              !
! opt_only2b = 0 reads an extended file conaitining 0- and 1-body pieces       !
!            = 1 reads the only 2-body piece (see format in the BoccaDorata    !
!                manual)                                                       !
! opt_onlycom = 0 for the main part of the interaction                         !
!             = 1 for the 2-body COM correction                                !
!                                                                              !
! The "generalized" format including the 0- and 1-body has the following       !
! format:                                                                      !
! L00                                                                          !
! 32.782169                                                                    !
! L11 format : n nn twoj pi ch                                                 !
! 144                                                                          !
! 0 0 1 0 0  4.682396                                                          !
! ... (144 matrix elements in total)                                           !
! L22 format is: na 2ja ipa cha nb... chd J Vpp                                !
! 482383                                                                       !
! 0 1 0 0 0 1 0 0  0 1 0 0 0 1 0 0  0 -7.346050                                !
! ... (482383 matrix elements in total)                                        !
!==============================================================================!
PROGRAM converter_bcd_to_taurus

use Iso_fortran_env

implicit none

!!! Definitions for portability  
integer, parameter :: uti = input_unit,  &
                      uto = output_unit, &
                      uth = uti + uto, &
                      r64 = real64 

!!! Basis
real(r64) :: HO_hw
integer :: HOsh_dim, HOsp_dim, HOsh_maxjj
integer, dimension(:), allocatable :: HOsh_n, HOsh_l, HOsh_2j, HOsh_na, &
                                      HOsp_n, HOsp_l, HOsp_2j, HOsp_2mj, &
                                      HOsp_2mt, HOsp_sh
integer, dimension(:,:,:), allocatable :: nlj_tab

!!! Hamiltonian
integer :: hamil_type, opt_only2b, opt_onlycom
real(r64) :: hamil_H0
real(r64), dimension(:,:,:), allocatable :: hamil_H1
real(r64), dimension(:,:,:,:,:,:), allocatable :: hamil_H2
character(len=100) :: hamil_fbcd, hamil_name
character(len=:), allocatable :: hamil_fsho, hamil_f01b, hamil_f2b, hamil_fcom
character(4) :: appsho='.sho', app01b='.01b', appcom='.com'
character(3) ::  app2b='.2b'

!!! Various   
integer ::  i, j, k, l, m, n, s, t, emax, Nmax, lmin, & 
            t1, t2, t3, t4, n1, n2, n3, n4, l1, l2, l3, l4, & 
            j1, j2, j3, j4, k1, k2, k3, k4, p1, p2, p3, p4, &
            ch1, ch2, ch3, ch4, &
            jj, jjmin, jjmax, j12min, j12max, j34min, j34max, &
            find_iso
real(r64) :: xj1, xj2, xj3, xj4, delta_12, delta_34, phas12, phas34, &
             NJ_12, NJ_34, phasJ, xj, Htmp
character(11), dimension(6) :: inp_ch
logical :: lexists 
character(2) :: emax_ch
character(10) :: HO_hw_ch 
character(len=*), parameter :: format1 = "(1a)", &
                               format2 = "(1i4,500i7)", &
                               format3 = "(1es20.12)", &
                               format4 = "(2i4,1x,4i7,1x,2i4)", &
                               format5 = "(6es20.12)"

!------------------------------------------------------------------------------!
! Reads the inputs                                                             !
!------------------------------------------------------------------------------!

read(uti, '(a11, a100)')   inp_ch(1), hamil_fbcd 
read(uti, '(a11, a100)')   inp_ch(2), hamil_name
read(uti, '(a11, i2)')     inp_ch(3), emax   
read(uti, '(a11, f10.4)')  inp_ch(4), HO_hw  
read(uti, '(a11, i1)')     inp_ch(5), opt_only2b
read(uti, '(a11, i1)')     inp_ch(6), opt_onlycom

!!! Name of hamiltonian files 
j = len_trim(adjustl(hamil_name))
allocate (character(j+4) :: hamil_fsho)
allocate (character(j+4) :: hamil_f01b)
allocate (character(j+3) :: hamil_f2b)
allocate (character(j+4) :: hamil_fcom)

hamil_fsho = trim(hamil_name) // appsho
hamil_f01b = trim(hamil_name) // app01b
hamil_f2b  = trim(hamil_name) // app2b
hamil_fcom = trim(hamil_name) // appcom

if ( opt_only2b == 0 ) then
  hamil_type = 3
else
  hamil_type = 4
endif

!!! Checks
if ( emax < 0 ) then
  print "(a,1i3,a)","The value for the number of oscillator shells emax = ", &  
         emax," should be positive!"
  stop 
endif 

if ( HO_hw < 0 ) then
  print "(a,1f10.4,a)","The value for the oscillator frequency hw = ",HO_hw, &
         " should be positive!"
  stop 
endif 

if ( opt_only2b > 1 ) then
  print "(a,1i1,a)","The value for the opt_only2b = ",opt_only2b, &
         " should be 0 or 1!"
  stop 
endif 

if ( opt_onlycom > 1 ) then
  print "(a,1i1,a)","The value for the opt_onlycom = ",opt_onlycom, &
         " should be 0 or 1!"
  stop 
endif 

inquire (file = hamil_fbcd, exist = lexists)
if ( lexists .eqv. .false. ) then
  print "(a,a,a)","The interaction file hamil_fbcd = ",hamil_fbcd, &
        " does not exists!"
  stop
endif 

!!! Prints the inputs
write(emax_ch,'(1i2)') emax
emax_ch = adjustl(emax_ch)
write(HO_hw_ch,'(1f10.4)') HO_hw
HO_hw_ch = adjustl(HO_hw_ch)

write(uto, '(a11, a100)')   inp_ch(1), hamil_fbcd 
write(uto, '(a11, a100)')   inp_ch(2), hamil_name
write(uto, '(a11, a2)')     inp_ch(3), emax_ch
write(uto, '(a11, a10)')    inp_ch(4), HO_hw_ch
write(uto, '(a11, i1)')     inp_ch(5), opt_only2b
write(uto, '(a11, i1)')     inp_ch(6), opt_onlycom

!!! Open the units
open(uth+1, file=hamil_fbcd, status='old', action='read')

if ( opt_only2b == 0 ) then
  open(uth+3, file=hamil_f01b, status='replace', action='write', &
       form='formatted')
endif

if ( opt_onlycom == 0 ) then
  open(uth+2, file=hamil_fsho, status='replace',action='write',form='formatted')
  open(uth+4, file=hamil_f2b,  status='replace',action='write',form='formatted')
else
  open(uth+4, file=hamil_fcom, status='replace',action='write',form='formatted')
endif

!------------------------------------------------------------------------------!
! Builds the basis                                                             !
!------------------------------------------------------------------------------!

Nmax = emax 
HOsh_dim = (Nmax*(Nmax+3))/2 + 1
HOsh_maxjj = 2*Nmax + 1 

allocate( HOsh_n(HOsh_dim), HOsh_l(HOsh_dim), HOsh_2j(HOsh_dim), &
          HOsh_na(HOsh_dim) )
allocate( nlj_tab(0:Nmax/2,0:Nmax,HOsh_maxjj) , source=0 )

i = 1

do N = 0, Nmax
  lmin = 0 + (1-(-1)**N)/2
  do l = lmin, N, 2
    do s = 1, -1, -2
      if ( 2*l+s < 1 ) cycle
      HOsh_n(i)  = (N - l)/2
      HOsh_l(i)  = l
      HOsh_2j(i) = 2*l + s
      HOsh_na(i) = 10000*HOsh_n(i) + 100*HOsh_l(i) + HOsh_2j(i)
      nlj_tab(HOsh_n(i),HOsh_l(i),HOsh_2j(i)) = i
      i = i + 1
    enddo
  enddo
enddo

!!! Writes the file .sho 
if ( opt_onlycom == 0 ) then
  write(uth+2,format1) hamil_name
  write(uth+2,format2) hamil_type
  write(uth+2,format2) HOsh_dim, (HOsh_na(i), i=1, HOsh_dim)
  write(uth+2,format3) HO_hw

  close(uth+2, status='keep')
endif

!------------------------------------------------------------------------------!
! Reads/Writes the 0- and 1- body part of the Hamiltonian                      !
!------------------------------------------------------------------------------!

if ( (opt_onlycom == 0) .and. (opt_only2b == 0) ) then
  allocate ( hamil_H1(HOsh_dim,HOsh_dim,-1:1) )
  hamil_H1 = 0.d0

  !!! Reads from the bcd file
  read(uth+1,*)
  read(uth+1,*) hamil_H0
  read(uth+1,*)
  read(uth+1,*) n

  do i = 1, n
    read(uth+1,*) n1, n2, j1, p1, ch1, Htmp

    l1 = (j1+1)/2 - (1+(-1)**(p1+j1/2))/2
    t1 = 1 - 2*ch1

    k1 = nlj_tab(n1,l1,j1)
    k2 = nlj_tab(n2,l1,j1)

    hamil_H1(k1,k2,t1) = Htmp
    hamil_H1(k2,k1,t1) = Htmp
  enddo

  !!! Writes the file .01b
  write(uth+3,format1) hamil_name
  write(uth+3,*) hamil_H0

  do k1 = 1, HOsh_dim
    do k2 = k1, HOsh_dim
      if ( HOsh_2j(k1) /= HOsh_2j(k2) ) cycle
      if ( (-1)**HOsh_l(k1) /= (-1)**HOsh_l(k2) ) cycle
      write(uth+3,*) HOsh_na(k1), HOsh_na(k2), hamil_H1(k2,k1,-1), & 
                 hamil_H1(k2,k1,+1)
    enddo 
  enddo 

  close(uth+3, status='keep')
endif

!------------------------------------------------------------------------------!
! Reads the 2- body part of the Hamiltonian                                    !
!                                                                              !
! The different "isospin" components are:                                      !
! T = 0   pp pp                                                                !
! T = 1   pn pn                                                                !
! T = 2   pn np                                                                !
! T = 3   np pn                                                                !
! T = 4   np np                                                                !
! T = 5   nn nn                                                                !
!------------------------------------------------------------------------------!

allocate ( hamil_H2(HOsh_dim,HOsh_dim,HOsh_dim,HOsh_dim,0:HOsh_maxjj,0:5) )
hamil_H2 = 0.d0

if ( (opt_onlycom == 0) .and. (opt_only2b == 0) ) then
  read(uth+1,*) 
endif
read(uth+1,*) 

DO
  read(uth+1,*,end=999) n1, j1, p1, ch1, n2, j2, p2, ch2, n3, j3, p3, ch3, &
                        n4, j4, p4, ch4, J, Htmp

  l1 = (j1+1)/2 - (1+(-1)**(p1+j1/2))/2
  t1 = 1 - 2*ch1

  l2 = (j2+1)/2 - (1+(-1)**(p2+j2/2))/2
  t2 = 1 - 2*ch2

  l3 = (j3+1)/2 - (1+(-1)**(p3+j3/2))/2
  t3 = 1 - 2*ch3

  l4 = (j4+1)/2 - (1+(-1)**(p4+j4/2))/2
  t4 = 1 - 2*ch4

  k1 = nlj_tab(n1,l1,j1)
  k2 = nlj_tab(n2,l2,j2)
  k3 = nlj_tab(n3,l3,j3)
  k4 = nlj_tab(n4,l4,j4)

  xj  = J * 1.d0                                                              
  xj1 = j1 / 2.d0
  xj2 = j2 / 2.d0
  xj3 = j3 / 2.d0
  xj4 = j4 / 2.d0
  phas12 = (-1.d0)**( xj1 + xj2 + xj + 1.d0 )
  phas34 = (-1.d0)**( xj3 + xj4 + xj + 1.d0 )
  phasJ  = (-1.d0)**( xj + 1.d0 )                                              

  delta_12 = 0.d0                                                          
  if ( (k1 == k2) .and. (t1 == t2) ) delta_12 = 1.d0 
  delta_34 = 0.d0                                                          
  if ( (k3 == k4) .and. (t3 == t4) ) delta_34 = 1.d0

  NJ_12 = sqrt(1.d0 - delta_12*phasJ) / (1.d0 + delta_12)               
  NJ_34 = sqrt(1.d0 - delta_34*phasJ) / (1.d0 + delta_34)               

  if (     t1 == -1 .and. t2 == -1 .and. t3 == -1 .and. t4 == -1 ) then
    T = 0
  elseif ( t1 == -1 .and. t2 == +1 .and. t3 == -1 .and. t4 == +1 ) then
    T = 1
  elseif ( t1 == -1 .and. t2 == +1 .and. t3 == +1 .and. t4 == -1 ) then
    T = 2
  elseif ( t1 == +1 .and. t2 == -1 .and. t3 == -1 .and. t4 == +1 ) then
    T = 3
  elseif ( t1 == +1 .and. t2 == -1 .and. t3 == +1 .and. t4 == -1 ) then
    T = 4
  elseif ( t1 == +1 .and. t2 == +1 .and. t3 == +1 .and. t4 == +1 ) then
    T = 5
  endif

  hamil_H2(k1,k2,k3,k4,J,T) = Htmp 
  hamil_H2(k1,k2,k4,k3,J,find_iso(T,1,1)) = Htmp * phas34
  hamil_H2(k2,k1,k3,k4,J,find_iso(T,2,1)) = Htmp * phas12
  hamil_H2(k2,k1,k4,k3,J,find_iso(T,3,1)) = Htmp * phas12 * phas34
  hamil_H2(k3,k4,k1,k2,J,find_iso(T,4,1)) = Htmp
  hamil_H2(k3,k4,k2,k1,J,find_iso(T,5,1)) = Htmp * phas12
  hamil_H2(k4,k3,k1,k2,J,find_iso(T,6,1)) = Htmp * phas34
  hamil_H2(k4,k3,k2,k1,J,find_iso(T,7,1)) = Htmp * phas12 * phas34

ENDDO
999 continue 

close(uth+1, status='keep')

!------------------------------------------------------------------------------!
! Writes the 2-body matrix elements à la Antoine.                              !
!                                                                              !
! 0 5 shell1 shell2 shell3 shell4 jmin jmax                                    !
! V(jmin,T=0) ...  V(jmin,T=51)                                                !
! ...                                                                          !
! V(jmax,T=0) ...  V(jmax,T=51)                                                !
!------------------------------------------------------------------------------!

write(uth+4,format1) hamil_name

do i = 1, HOsh_dim
  n1 = HOsh_n(i)
  l1 = HOsh_l(i)
  j1 = HOsh_2j(i)
  do k = i, HOsh_dim
    n3 = HOsh_n(k)
    l3 = HOsh_l(k)
    j3 = HOsh_2j(k)
    do l = k, HOsh_dim 
      n4 = HOsh_n(l)
      l4 = HOsh_l(l)
      j4 = HOsh_2j(l)
      m = HOsh_dim 
      if ( k == i ) m = l
      do j = i, m    
        n2 = HOsh_n(j)
        l2 = HOsh_l(j)
        j2 = HOsh_2j(j)
    
        if ( (-1)**(l1+l2) /= (-1)**(l3+l4) ) cycle 
 
        j12min = abs(j1-j2)/2
        j12max = (j1+j2)/2
        j34min = abs(j3-j4)/2
        j34max = (j3+j4)/2

        jjmin = max(j12min,j34min)
        jjmax = min(j12max,j34max)
        
        if ( jjmin <= jjmax ) then
          write(uth+4,format4) 0, 5, HOsh_na(i), HOsh_na(j), HOsh_na(k), &
                           HOsh_na(l), jjmin, jjmax
          do jj = jjmin, jjmax
            write(uth+4,format5) hamil_H2(i,j,k,l,jj,0), &
                                 hamil_H2(i,j,k,l,jj,1), &
                                 hamil_H2(i,j,k,l,jj,2), & 
                                 hamil_H2(i,j,k,l,jj,3), & 
                                 hamil_H2(i,j,k,l,jj,4), &
                                 hamil_H2(i,j,k,l,jj,5)
          enddo
        endif
      enddo
    enddo
  enddo
enddo

close(uth+4, status='keep')

write(uto,'("This is the end, my only friend, the end.")') 

END PROGRAM converter_bcd_to_taurus

!------------------------------------------------------------------------------!
! function find_iso                                                            ! 
!                                                                              !
! Determines the correct isospin after symmetry reconstruction. This function  !
! is needed because symmetries are used to reduced the number of matrix eleme- !
! nts stored.                                                                  !
! The definition of the permutations is:                                       !
! 1 = abdc                                                                     !
! 2 = bacd                                                                     !
! 3 = badc                                                                     !
! 4 = cdab                                                                     !
! 5 = cdba                                                                     !
! 6 = dcab                                                                     !
! 7 = dcba                                                                     !
!                                                                              !
! Input                                                                        !
! Tini = initial value of T                                                    !
! perm = permutation to consider                                               !
! inter = type of the interaction                                              !
!                                                                              !
! Output                                                                       !
! Tfin = final value of T                                                      !
!------------------------------------------------------------------------------!
integer function find_iso(Tini,perm,inter) result(Tfin)

implicit none
integer, intent(in) :: Tini, perm, inter

if ( inter == 0 ) then
  Tfin = Tini
elseif ( inter == 1 ) then
  if ( Tini == 0 .or. Tini == 5 ) then
    Tfin = Tini
  elseif ( Tini == 1 ) then
    if ( perm == 1 ) Tfin = 2
    if ( perm == 2 ) Tfin = 3        
    if ( perm == 3 ) Tfin = 4         
    if ( perm == 4 ) Tfin = 1         
    if ( perm == 5 ) Tfin = 2         
    if ( perm == 6 ) Tfin = 3         
    if ( perm == 7 ) Tfin = 4        
  elseif ( Tini == 2 ) then
    if ( perm == 1 ) Tfin = 1
    if ( perm == 2 ) Tfin = 4        
    if ( perm == 3 ) Tfin = 3         
    if ( perm == 4 ) Tfin = 3         
    if ( perm == 5 ) Tfin = 4         
    if ( perm == 6 ) Tfin = 1         
    if ( perm == 7 ) Tfin = 2        
  elseif ( Tini == 3 ) then
    if ( perm == 1 ) Tfin = 4
    if ( perm == 2 ) Tfin = 1        
    if ( perm == 3 ) Tfin = 2         
    if ( perm == 4 ) Tfin = 2         
    if ( perm == 5 ) Tfin = 1         
    if ( perm == 6 ) Tfin = 4         
    if ( perm == 7 ) Tfin = 3        
  elseif ( Tini == 4 ) then
    if ( perm == 1 ) Tfin = 3
    if ( perm == 2 ) Tfin = 2        
    if ( perm == 3 ) Tfin = 1         
    if ( perm == 4 ) Tfin = 4         
    if ( perm == 5 ) Tfin = 3         
    if ( perm == 6 ) Tfin = 2         
    if ( perm == 7 ) Tfin = 1        
  endif
endif

end function find_iso      

!==============================================================================!
! End of file                                                                  !
!==============================================================================!
