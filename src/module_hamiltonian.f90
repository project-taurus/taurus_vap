!==============================================================================!
! MODULE Hamiltonian                                                           !
!                                                                              !
! This module contains the variables and routines related to the Hamiltonian.  !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_hamiltonian                                                 !
! - subroutine set_hamiltonian_1body                                           !
! - subroutine read_hamiltonian_1body_antoine                                  !
! - subroutine read_hamiltonian_1body_general                                  !
! - subroutine set_hamiltonian_2body                                           !
! - subroutine read_hamiltonian_2body_antoine                                  !
! - subroutine read_hamiltonian_2body_Jscheme                                  !
! - subroutine read_hamiltonian_2body_com                                      !
! - subroutine read_hamiltonian_reduced                                        !
! - subroutine write_hamiltonian_reduced                                       !
! - subroutine print_hamiltonian                                               !
! - subroutine set_kinetic_bare                                                ! 
! - subroutine uncouple_J                                                      !
! - subroutine uncouple_JT                                                     !
! - subroutine reconstruct_2body_timerev                                       !
! - subroutine find_timerev                                                    !
! - subroutine calculate_expectval_energy                                      !
! - subroutine calculate_decompo_energy                                        !
! - function find_iso                                                          !
!==============================================================================!
MODULE Hamiltonian

use Basis

implicit none
public       

!!! Basic quantities 
character(100) :: hamil_name ! name
integer :: hamil_read, &     ! option to read the reduced file
           hamil_type, &     ! type of interaction (ChiEFT, ANTOINE, ...)
           hamil_com         ! option for center-of-mass correction

!!! Name of files 
character(len=:), allocatable :: hamil_file, & ! main name            
                                 hamil_fsho, & ! name file model space
                                 hamil_f01b, & !  "    "   0+1-body
                                 hamil_f2b,  & !  "    "   2-body
                                 hamil_fred, & !  "    "   reduced 
                                 hamil_fcom    !  "    "   center of mass

!!! Matrix elements in coupled scheme: JT or J
!real(r64), dimension(:,:,:), allocatable :: hamil_H1cpd ! 1-body
real(r64), dimension(:,:,:,:,:,:), allocatable, private :: hamil_H2cpd ! 2-body 

!!! Matrix elements in m-scheme (HO basis)  
complex(r64) :: hamil_H0                               ! 0-body part 
complex(r64), dimension(:,:), allocatable :: hamil_H1  ! 1-body part 
real(rH2), dimension(:), allocatable :: hamil_H2       ! 2-body part 
integer(i64) :: hamil_H2dim, hamil_H2dim_all           ! number of 2BME stored 
integer(i16), dimension(:), allocatable :: hamil_abcd  ! indices of 2BME   
integer(i8), dimension(:), allocatable :: hamil_trperm ! time reversal permut.

!!! Other quantities
real(r64), private :: scaling_1b, & ! scaling factor 1-body mat. elem.
                      scaling_2b    !    "       "   2-body  "     " 

!!! Private routines
private :: set_hamiltonian_01body, read_hamiltonian_1body_antoine, &
           read_hamiltonian_1body_general, print_hamiltonian, &
           set_hamiltonian_2body, read_hamiltonian_2body_antoine, & 
           read_hamiltonian_2body_Jscheme, read_hamiltonian_reduced, &
           write_hamiltonian_reduced, uncouple_J, uncouple_JT, find_iso 

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_hamiltonian                                                   !
!                                                                              ! 
! Sets the matrices of the Hamiltonian and reads their matrix elements from    ! 
! the interaction file specified in the input file.                            ! 
!------------------------------------------------------------------------------!
subroutine set_hamiltonian       

integer :: ialloc=0

!!! Initiliazes arrays
allocate(hamil_H1(HOsp_dim,HOsp_dim), stat=ialloc)
if ( ialloc /= 0 ) stop 'Error during allocation of Hamiltonian'

hamil_H0 = zzero
hamil_H1 = zzero

!!! Reads and/or computes the matrix elements (in HO basis)
if ( hamil_read == 0 ) then 
  call set_hamiltonian_01body
  call set_hamiltonian_2body
  call write_hamiltonian_reduced
else 
  call read_hamiltonian_reduced
endif

!!! Prints the basic information on the Hamiltonian
!cmpi if ( paral_myrank == 0 ) then        
call print_hamiltonian
!cmpi endif

close(uth, status='keep')

end subroutine set_hamiltonian       

!------------------------------------------------------------------------------!
! subroutine set_hamiltonian_1body                                             !
!                                                                              ! 
! Sets the one-body part of the Hamiltonian: either reads from file in case    ! 
! of shell-model and normal-order interactions or constructs the bare kinetic  !
! operator (possibly with 1b part of the com correction)                       ! 
!------------------------------------------------------------------------------!
subroutine set_hamiltonian_01body

scaling_1b = 1.0d0
scaling_2b = 1.0d0

rewind(uth)
read(uth,'(a)') hamil_name
read(uth,*) hamil_type 

select case (hamil_type)

  !!! single-particle energies à la ANTOINE 
  case (1:2)
    call read_hamiltonian_1body_antoine

  !!! 1-body + (com)
  case (3)
    call read_hamiltonian_1body_general
    if ( hamil_com == 1 ) call set_kinetic_bare

  !!! Kinetic + com
  case(4) 
    call set_kinetic_bare

end select 
  
end subroutine set_hamiltonian_01body

!------------------------------------------------------------------------------!
! subroutine read_hamiltonian_1body_antoine                                    !
!                                                                              !
! Reads the one-body matrix elements of the hamiltonian in antoine format.     !
! I tried to respect as much as possible their variables names (based on the   !
! what can be read on their website).                                          !
!------------------------------------------------------------------------------!
subroutine read_hamiltonian_1body_antoine
  
use Nucleus, only: core_A

integer(i32) :: i, j, k, qmax, jmax, idens, core(2), ialloc=0
real(r64) :: x1, x2, zmass
real(r64), dimension(:,:), allocatable :: zepsi      

qmax = hamil_type

allocate( zepsi(HOsh_dim,qmax), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of s.p. energies'

do i = 1, qmax
  read(uth,*) (zepsi(k,i), k=1, HOsh_dim)
enddo
read(uth,*) idens, core(1), core(2), zmass

! Mass scaling for the shell-model matrix elements
x1 = nucleus_A   
x2 = core_A + 2

if ( idens == 1 ) then
  scaling_1b = 1.d0
  scaling_2b = (x2/x1)**(zmass)
elseif ( idens == 2 ) then
  scaling_1b = (x2/x1)**(zmass)
  scaling_2b = (x2/x1)**(zmass)
endif

! 1-body part
k = 0
do i = 1, HOsh_dim
  jmax = HOsh_2j(i) + 1
  do j = 1, jmax
    k = k + 1
    hamil_H1(k,k) = zone * zepsi(i,1) * scaling_1b
    hamil_H1(k+HOsp_dim/2,k+HOsp_dim/2) = zone * zepsi(i,qmax) * scaling_1b
  enddo
enddo

deallocate (zepsi)
  
end subroutine read_hamiltonian_1body_antoine

!------------------------------------------------------------------------------!
! subroutine read_hamiltonian_1body_general                                    !
!                                                                              !
! Reads the one-body matrix elements of the hamiltonian in the general format. !
!------------------------------------------------------------------------------!
subroutine read_hamiltonian_1body_general
  
integer(i32) :: k, l, mjk, mjl, a, b, iexit !, ialloc=0
real(r64) :: H0, hamil_H1cpd_p, hamil_H1cpd_n

!!! Allocate hamil_H1cpd (not raelly needed but prettier)
!allocate( hamil_H1cpd(0:1,HOsh_dim,HOsh_dim), stat=ialloc ) 
!if ( ialloc /= 0 ) stop 'Error during allocation of coupled hamiltonian'
!hamil_H1cpd = zero

!!! Reads the file
read(uth1,*) 
read(uth1,*) H0
hamil_H0 = zone * H0

do 
  read(uth1,*,iostat=iexit) a, b, hamil_H1cpd_p, hamil_H1cpd_n
  if ( iexit /= 0 ) exit

  do k = 1, HOsp_dim/2
    mjk = HOsp_2mj(k)
    do l = 1, HOsp_dim/2
      mjl = HOsp_2mj(l)
      if ( mjk /= mjl ) cycle
      if ( (HOsh_na(HOsp_sh(k)) == a) .and. &
           (HOsh_na(HOsp_sh(l)) == b) ) then
         hamil_H1(k,l) = zone * hamil_H1cpd_p
         hamil_H1(l,k) = zone * hamil_H1cpd_p
         hamil_H1(k+HOsp_dim/2,l+HOsp_dim/2) = zone * hamil_H1cpd_n
         hamil_H1(l+HOsp_dim/2,k+HOsp_dim/2) = zone * hamil_H1cpd_n
      endif
    enddo
  enddo
enddo 

close (uth1, status='keep')
 
!deallocate (hamil_H1cpd)
 
end subroutine read_hamiltonian_1body_general

!------------------------------------------------------------------------------!
! subroutine set_hamiltonian_2body                                             !
!                                                                              ! 
! Reads the 2-body part of the interaction from file "hamil_f2bb" and uncouples ! 
! it before storing only the non-zero matrix elements to be used later.        ! 
!------------------------------------------------------------------------------!
subroutine set_hamiltonian_2body

integer(i16) :: ared, bred, cred, dred
integer(i32) :: ht, j, t, tmax, uth6=uth+6, uth7=uth+7, fac_ht, ialloc=0, &
                a, ma, la, ta, b, mb, lb, tb, bmax, &
                c, mc, lc, tc, d, md, ld, td
integer(i64) :: kk
real(r64) :: xja, xjb, xjc, xjd, xjtot, xttot, phasab, phascd, Vtmp, &
             Vcut, Vdec
real(rH2) :: Vred

!!! Allocates the arrays and reads the matrix elements
select case (hamil_type)
  case (1:2)
    allocate( hamil_H2cpd(0:3,0:HO_2jmax,HOsh_dim,HOsh_dim,HOsh_dim,HOsh_dim), &
              stat=ialloc ) 
    if ( ialloc /= 0 ) stop 'Error during allocation of coupled hamiltonian'
    hamil_H2cpd = zero
 
    call read_hamiltonian_2body_antoine
 
  case (3:4)
    allocate( hamil_H2cpd(0:5,0:HO_2jmax,HOsh_dim,HOsh_dim,HOsh_dim,HOsh_dim), &
              stat=ialloc ) 
    if ( ialloc /= 0 ) stop 'Error during allocation of coupled hamiltonian'
    hamil_H2cpd = zero
 
    call read_hamiltonian_2body_Jscheme
    if ( hamil_com == 1 ) call read_hamiltonian_2body_com    
 
  case default
    stop 'Format of the hamiltonian not recognized!'
end select

!!! Symmetry reconstruction in case only part of the matrix elements are
!!! written in the file (usually the case)
!call reconstruct_hamiltonian_2body (hamil_type)
ht = hamil_type
fac_ht = kdelta(ht,1) + kdelta(ht,2)
tmax = 5 
if ( (hamil_type == 1) .or. (hamil_type == 2) ) tmax = 3

do a = 1, HOsh_dim 
  do b = 1, HOsh_dim
    do c = 1, HOsh_dim
      do d = 1, HOsh_dim
        do J = 0, HO_2jmax
          do T = 0, tmax 
            xja = HOsh_2j(a) / 2.d0
            xjb = HOsh_2j(b) / 2.d0
            xjc = HOsh_2j(c) / 2.d0
            xjd = HOsh_2j(d) / 2.d0
            xjtot = 1.d0 * J      
            xttot = 1.d0 - 1.0d0 * kdelta(T,0) * fac_ht
            phasab = (-1.d0)**(xja + xjb + xjtot + xttot)
            phascd = (-1.d0)**(xjc + xjd + xjtot + xttot)
            
            Vtmp = hamil_H2cpd(T,J,a,b,c,d)
            if ( abs(Vtmp) > 1.0d-12 ) then
              hamil_H2cpd(find_iso(T,1,ht),J,a,b,d,c) = Vtmp * phascd
              hamil_H2cpd(find_iso(T,2,ht),J,b,a,c,d) = Vtmp * phasab
              hamil_H2cpd(find_iso(T,3,ht),J,b,a,d,c) = Vtmp * phasab * phascd
              hamil_H2cpd(find_iso(T,4,ht),J,c,d,a,b) = Vtmp
              hamil_H2cpd(find_iso(T,5,ht),J,c,d,b,a) = Vtmp * phasab
              hamil_H2cpd(find_iso(T,6,ht),J,d,c,a,b) = Vtmp * phascd
              hamil_H2cpd(find_iso(T,7,ht),J,d,c,b,a) = Vtmp * phasab * phascd
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
enddo

!!! Computes the two-body matrix elements in m-scheme 
open (uth6, status='scratch', action='readwrite', access='stream', & 
           form='unformatted')              
open (uth7, status='scratch', action='readwrite', access='stream', & 
           form='unformatted')              

Vcut = 1.0d-16
hamil_H2dim = 0

do a = 1, HOsp_dim
  la = HOsp_l(a)
  ma = HOsp_2mj(a)
  ta = HOsp_2mt(a)
  do c = a, HOsp_dim
    lc = HOsp_l(c)
    mc = HOsp_2mj(c)
    tc = HOsp_2mt(c)
    do d = c+1, HOsp_dim
      ld = HOsp_l(d)
      md = HOsp_2mj(d)
      td = HOsp_2mt(d)
      bmax = HOsp_dim
      if ( c == a ) bmax = d
      do b = a+1, bmax
        lb = HOsp_l(b)
        mb = HOsp_2mj(b)
        tb = HOsp_2mt(b)

        if ( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
        if ( ta + tb /= tc + td ) cycle
        if ( ma + mb /= mc + md ) cycle
        if ( ma + mb < 0 ) cycle
 
        if ( (hamil_type == 1) .or. (hamil_type == 2) ) then
          call uncouple_JT(a,b,c,d,Vdec) 
        else
          call uncouple_J(a,b,c,d,Vdec) 
        endif

        Vdec = Vdec * scaling_2b 
    
        !!! Select only matrix elements above a given cutoff to reduce the
        !!! CPU time and storage
        if ( abs(Vdec) > Vcut ) then 
          hamil_H2dim = hamil_H2dim + 1
          ared = int(a,i16)
          bred = int(b,i16)
          cred = int(c,i16)
          dred = int(d,i16)
          Vred = real(Vdec,r32)
          write(uth6) ared, bred, cred, dred      
          write(uth7) Vred      
        endif  
         
      enddo  !end loop b
    enddo  !end loop d
  enddo  !end loop c
enddo  !end loop a
  
deallocate (hamil_H2cpd)

hamil_H2dim_all = hamil_H2dim

!!! Final allocation and reading of two-body matrix elements
allocate ( hamil_H2(hamil_H2dim), hamil_abcd(4*hamil_H2dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of array of indices'

rewind(uth6) 
rewind(uth7) 

read(uth6) (hamil_abcd(kk), kk=1, 4*hamil_H2dim)
read(uth7) (hamil_H2(kk), kk=1, hamil_H2dim)

close(uth6)
close(uth7)

!!! Determines the permutation needed to obtain the time-reversed two-body
!!! matrix elements
call reconstruct_2body_timerev

end subroutine set_hamiltonian_2body

!------------------------------------------------------------------------------!
! subroutine read_hamiltonian_2body_antoine                                    !
!                                                                              !
! Reads the two-body matrix elements of the hamiltonian in antoine format      !
! (with JT-scheme).                                                            !
!------------------------------------------------------------------------------!
subroutine read_hamiltonian_2body_antoine
  
integer(i32) :: i, j, k, T, tmin, tmax, jmin, jmax, ialloc=0, iexit,  &
                a, b, c, d, nlj_tmp(4)
real(r64), dimension(:), allocatable :: H2_tmp   
  
allocate( H2_tmp(0:HO_2jmax), stat=ialloc  ) 
if ( ialloc /= 0 ) stop 'Error during allocation of coupled hamiltonian'

H2_tmp = zero
  
do 
  read(uth,*,iostat=iexit) tmin, tmax, (nlj_tmp(k),k=1,4), jmin, jmax
  if ( iexit /= 0 ) exit

  do i = 1, HOsh_dim
    if ( nlj_tmp(1) == HOsh_na(i) ) a = i
    if ( nlj_tmp(2) == HOsh_na(i) ) b = i
    if ( nlj_tmp(3) == HOsh_na(i) ) c = i
    if ( nlj_tmp(4) == HOsh_na(i) ) d = i
  enddo

  do T = tmin, tmax
    H2_tmp = zero    
    read(uth,*) (H2_tmp(J), J=jmin,jmax)
  
    select case (T)
      case (0)
        hamil_H2cpd(0,:,a,b,c,d) = H2_tmp(:) 
      case (1)
        do i = 2, 5-tmax
          j = i
          if ( i == 4 ) j = 1
          hamil_H2cpd(j,:,a,b,c,d) = H2_tmp(:) 
        enddo
      case (2)
        hamil_H2cpd(T+(-1)**(tmax+1),:,a,b,c,d) = H2_tmp(:) 
      case (3)
        hamil_H2cpd(1,:,a,b,c,d) = H2_tmp(:) 
      case default
        stop 'Wrong interaction (Antoine) format: isospin!'
    end select
  enddo 
enddo
 
deallocate (H2_tmp)
  
end subroutine read_hamiltonian_2body_antoine

!------------------------------------------------------------------------------!
! subroutine read_hamiltonian_2body_Jscheme                                    !
!                                                                              !
! Reads the two-body matrix elements of the hamiltonian in J-scheme.           !
!------------------------------------------------------------------------------!
subroutine read_hamiltonian_2body_Jscheme
  
integer(i32) :: i, j, T, tmin, tmax, jmin, jmax, iexit, a, b, c, d, nlj_tmp(4)
  
!Skipping the name
read(uth2,*)    
  
do 
  read(uth2,*,iostat=iexit) tmin, tmax, (nlj_tmp(i),i=1,4), jmin, jmax
  if ( iexit /= 0 ) exit

  do i = 1, HOsh_dim
    if ( nlj_tmp(1) == HOsh_na(i) ) a = i
    if ( nlj_tmp(2) == HOsh_na(i) ) b = i
    if ( nlj_tmp(3) == HOsh_na(i) ) c = i
    if ( nlj_tmp(4) == HOsh_na(i) ) d = i
  enddo

  do j = jmin, jmax
    read(uth2,*) (hamil_H2cpd(T,j,a,b,c,d), T=tmin,tmax)
  enddo 
enddo

close (uth2, status='keep')
  
end subroutine read_hamiltonian_2body_Jscheme
  
!------------------------------------------------------------------------------!
! subroutine read_hamiltonian_2body_com                                        !
!                                                                              !
! Reads the two-body matrix elements of the center-of-mass correction and      !
! adds them to the rest of the Hamiltonian (J-scheme).                         !
!------------------------------------------------------------------------------!
subroutine read_hamiltonian_2body_com    
  
integer(i32) :: i, j, T, tmin, tmax, jmin, jmax, iexit, a, b, c, d, nlj_tmp(4)
real(r64) :: factor, H2com(0:5)
 
factor = HO_hw / (nucleus_A * hbarmass)
 
!Skipping the name
read(uthc,*)    
  
do 
  read(uthc,*,iostat=iexit) tmin, tmax, (nlj_tmp(i),i=1,4), jmin, jmax
  if ( iexit /= 0 ) exit

  do i = 1, HOsh_dim
    if ( nlj_tmp(1) == HOsh_na(i) ) a = i
    if ( nlj_tmp(2) == HOsh_na(i) ) b = i
    if ( nlj_tmp(3) == HOsh_na(i) ) c = i
    if ( nlj_tmp(4) == HOsh_na(i) ) d = i
  enddo

  do j = jmin, jmax
    read(uthc,*) (H2com(T), T=tmin,tmax)
    hamil_H2cpd(:,j,a,b,c,d) = hamil_H2cpd(:,j,a,b,c,d) + factor * H2com(:)
  enddo 
enddo

close (uthc, status='keep')
  
end subroutine read_hamiltonian_2body_com
  
!------------------------------------------------------------------------------!
! subroutine read_hamiltonian_reduced                                          !
!                                                                              !
! Reads the reduced Hamiltonian (for the TBME, only the one above a cutoff and !
! only a minimal symmetric subset of them)                                     !
!------------------------------------------------------------------------------!
subroutine read_hamiltonian_reduced
  
integer(i32) :: i, j, bytes_H2, ialloc=0
integer(i64) :: kk, pos_ini, pos_abcd, pos_H2, pos_perm
!cmpi integer(i64) :: divide, rest
  
read(uthr, pos=1) hamil_name
read(uthr) hamil_type     
read(uthr) hamil_H0      
read(uthr) ((hamil_H1(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
read(uthr) hamil_H2dim_all

hamil_H2dim = hamil_H2dim_all

!!! Determines the position of the arrays in the file. This is not needed in a 
!!! sequential run but I prefer to have the same code structure as in the 
!!! MPI version.
inquire(uthr, pos=pos_ini)

if ( rH2 == r32 ) then 
  bytes_H2 = 4
else
  bytes_H2 = 8
endif

pos_abcd = pos_ini
pos_H2 = pos_ini + hamil_H2dim_all * 8 
pos_perm = pos_ini + hamil_H2dim_all * (8 + bytes_H2)

!!! MPI: distribution of matrix elements among the membrs of the teams. 
!!! First computes the number of matrix elements and then their position in the
!!! files. 
!cmpi divide = hamil_H2dim_all / paral_myteamsize
!cmpi rest = modulo(hamil_H2dim_all,paral_myteamsize)

!cmpi hamil_H2dim = divide
!cmpi if ( paral_myteamrank < rest ) then
!cmpi   hamil_H2dim = hamil_H2dim + 1
!cmpi endif

!print*,"myrank,H2dim", paral_myrank, hamil_H2dim

!cmpi pos_abcd = pos_ini + paral_myteamrank * hamil_H2dim * 8 
!cmpi pos_H2 = pos_ini + hamil_H2dim_all * 8 + paral_myteamrank * hamil_H2dim & 
!cmpi                                           * bytes_H2
!cmpi pos_perm = pos_ini + hamil_H2dim_all * (8 + bytes_H2)  & 
!cmpi            + paral_myteamrank * hamil_H2dim * 1

!print*,"myrank,abcd,H2,perm", paral_myrank, pos_abcd, pos_H2, pos_perm

!cmpi if ( paral_myteamrank >= rest ) then
!cmpi   pos_abcd = pos_abcd + rest * 8
!cmpi   pos_H2 = pos_H2 + rest * bytes_H2
!cmpi   pos_perm = pos_perm + rest 
!cmpi endif

!!! Reads the arrays related to the two-body matrix elements
allocate ( hamil_H2(hamil_H2dim), hamil_abcd(4*hamil_H2dim), & 
           hamil_trperm(hamil_H2dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of array of indices'

read(uthr, pos=pos_abcd) (hamil_abcd(kk), kk=1, 4*hamil_H2dim)
read(uthr, pos=pos_H2) (hamil_H2(kk), kk=1, hamil_H2dim)
read(uthr, pos=pos_perm) (hamil_trperm(kk), kk=1, hamil_H2dim)
  
close(uthr, status='keep')

end subroutine read_hamiltonian_reduced

!------------------------------------------------------------------------------!
! subroutine write_hamiltonian_reduced                                         !
!                                                                              ! 
! Writes the reduced Hamiltonian (for the TBME, only the one above a cutoff    ! 
! and only a minimal symmetric subset of them)                                 ! 
!------------------------------------------------------------------------------!
subroutine write_hamiltonian_reduced
  
integer(i32) :: i, j
integer(i64) :: kk
  
open(uthr, file=hamil_fred, status='replace', action='write', access='stream', &
          form='unformatted')              
  
write(uthr) hamil_name 
write(uthr) hamil_type 
write(uthr) hamil_H0      
write(uthr) ((hamil_H1(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
write(uthr) hamil_H2dim
write(uthr) (hamil_abcd(kk), kk=1, 4*hamil_H2dim)
write(uthr) (hamil_H2(kk), kk=1, hamil_H2dim)
write(uthr) (hamil_trperm(kk), kk=1, hamil_H2dim)
  
close(uthr, status='keep')

end subroutine write_hamiltonian_reduced

!------------------------------------------------------------------------------!
! print_hamiltonian                                                            !
!                                                                              ! 
! Prints the properties of the Hamiltonian file defined in the inputs.         !
!------------------------------------------------------------------------------!
subroutine print_hamiltonian

character(12) :: hamil_H2dim_ch
character(len=:), allocatable :: info_type, info_read
character(len=*), parameter :: format1 = "(1a19,3x,1a)", &
                               format2 = "(1a19,3x,1i3,3x,1a)", &
                               format3 = "(1a19,3x,1a12)"

write(hamil_H2dim_ch,'(1i12)') hamil_H2dim_all
hamil_H2dim_ch = adjustl(hamil_H2dim_ch)

!!! Additional information on the type of interaction
select case (hamil_type)
  case (1)
    info_type = '(ANTOINE, spe_p==spe_n)'
  case (2)
    info_type = '(ANTOINE, spe_p/=spe_n)'
  case (3)
    info_type = '(J-scheme, NO2B)'
  case (4)
    info_type = '(J-scheme, K + 2B)'
  case default
    info_type = '(not recognized)'
end select 

!!! Additional information on the file read to obtain the matrix elements
select case (hamil_read)
  case (0)
    info_read = '(normal)'
  case (1:2)
    info_read = '(reduced)'
end select 

print '(/,60("%"),/,25x,"HAMILTONIAN",24x,/,60("%"),//, &
      & 3x,"Description",8x,"Value",/,27("-"))'
print format1, 'Main name of files ', hamil_file
print format1, 'Name of hamiltonian', trim(adjustl(hamil_name))
print format2, 'Type of hamiltonian', hamil_type, info_type
print format2, 'Read from file     ', hamil_read, info_read
print format3, 'Number of 2BME     ', hamil_H2dim_ch 

end subroutine print_hamiltonian

!------------------------------------------------------------------------------! 
! subroutine set_kinetic_bare                                                  ! 
!                                                                              ! 
! Calculates the matrix elements of the kinetic one-body operator using the    ! 
! formula:                                                                     ! 
! <na la ja ma ta|T|nb lb jb mb tb> = delta_la,lb * delta_ja,jb * delta_ma,mb  ! 
!                                    * delta ta,tb *                           ! 
!                                     (delta_na,nb * (2*nb + lb + 3/2)         ! 
!                                    + delta_na,nb-1 * sqrt(nb*(nb + lb + 1/2) ! 
!                                    + delta_na,nb+1 * sqrt(na*(na + lb + 1/2))! 
!                                                                              ! 
! If hamil_com = 1, also includes the one-body part of the center-of-mass      ! 
! correction using:                                                            !
! T + com(1b) = (1 - 1/A) * T                                                  ! 
!------------------------------------------------------------------------------! 
subroutine set_kinetic_bare                                                      
                                                                                 
integer(i32) :: a, na, la, ja, ma, ta, b, nb, lb, jb, mb, tb                     
real(r64) :: fac_na_nb, fac_na_nbm1, fac_na_nbp1, & 
             delta_na_nb, delta_na_nbm1, delta_na_nbp1 
complex(r64) :: fachw, faccom, fackin, factor                                              
                                                                            
fackin = kdelta(hamil_type,4) * zone
faccom = hamil_com * zone / nucleus_A
fachw = 0.5d0 * HO_hw * zone 

factor = fachw * (fackin - faccom)
                                                                                 
do a = 1, HOsp_dim                                                               
  na = HOsp_n(a)                                                                 
  la = HOsp_l(a)                                                                 
  ja = HOsp_2j(a)                                                                
  ma = HOsp_2mj(a)                                                               
  ta = HOsp_2mt(a)                                                               
  do b = 1, HOsp_dim                                                             
    nb = HOsp_n(b)                                                               
    lb = HOsp_l(b)                                                               
    jb = HOsp_2j(b)                                                              
    mb = HOsp_2mj(b)                                                             
    tb = HOsp_2mt(b)                                                             
    fac_na_nb = 2.d0 * nb + lb + 1.5d0                                           
    fac_na_nbm1 = sqrt(nb * (nb + lb + 0.5d0))                                   
    fac_na_nbp1 = sqrt(na * (na + la + 0.5d0))                                   
                                                                                 
   if ( (ta == tb) .and. (la == lb) .and. (ja == jb) .and. (ma == mb) ) then    
      delta_na_nb = one * kdelta(na,nb)                                         
      delta_na_nbm1 = one * kdelta(na,nb-1)                                     
      delta_na_nbp1 = one * kdelta(na,nb+1)                                     
      hamil_H1(a,b) =  hamil_H1(a,b) + factor * (delta_na_nb   * fac_na_nb    &                   
                                               + delta_na_nbm1 * fac_na_nbm1  &                   
                                               + delta_na_nbp1 * fac_na_nbp1 )                    
    endif                                                                       
                                                                                
  enddo                                                                         
enddo            
                                                                                 
end subroutine set_kinetic_bare                                                  

!------------------------------------------------------------------------------!
! subroutine uncouple_J                                                        !
!                                                                              ! 
! Decouples J-scheme two-body matrix elements to m-scheme: ab|V|cd>            ! 
! The different "isospin" components are:                                      ! 
! T = 0   pp pp                                                                ! 
! T = 1   pn pn                                                                ! 
! T = 2   pn np                                                                ! 
! T = 3   np pn                                                                ! 
! T = 4   np np                                                                ! 
! T = 5   nn nn                                                                ! 
!------------------------------------------------------------------------------!
subroutine uncouple_J (a,b,c,d,sumando)

integer(i32), intent(in) :: a, b, c, d
real(r64), intent(out) :: sumando
integer(i32) :: ja, ma, ta, sha, jb, mb, tb, shb, & 
                jc, mc, tc, shc, jd, md, td, shd, J, MJ, T
real(r64) :: cg1, cg2, N_AB, N_CD, N_INV, delta_ab, delta_cd, phasJ

ja = HOsp_2j(a)
jb = HOsp_2j(b)
jc = HOsp_2j(c)
jd = HOsp_2j(d)

ma = HOsp_2mj(a)
mb = HOsp_2mj(b)
mc = HOsp_2mj(c)
md = HOsp_2mj(d)
MJ = ma + mb

ta = HOsp_2mt(a)
tb = HOsp_2mt(b)
tc = HOsp_2mt(c)
td = HOsp_2mt(d)

sha = HOsp_sh(a)
shb = HOsp_sh(b)
shc = HOsp_sh(c)
shd = HOsp_sh(d)

if (     ta == -1 .and. tb == -1 .and. tc == -1 .and. td == -1 ) then
  T = 0
elseif ( ta == -1 .and. tb == +1 .and. tc == -1 .and. td == +1 ) then
  T = 1
elseif ( ta == -1 .and. tb == +1 .and. tc == +1 .and. td == -1 ) then
  T = 2
elseif ( ta == +1 .and. tb == -1 .and. tc == -1 .and. td == +1 ) then
  T = 3
elseif ( ta == +1 .and. tb == -1 .and. tc == +1 .and. td == -1 ) then
  T = 4
elseif ( ta == +1 .and. tb == +1 .and. tc == +1 .and. td == +1 ) then
  T = 5
endif

delta_ab = 0.d0
delta_cd = 0.d0
if ( T == 0 .or. T == 5 ) then
  delta_ab = one * kdelta(sha,shb)
  delta_cd = one * kdelta(shc,shd)
endif

sumando = 0.d0                                                            

do J = abs(ma+mb)/2, min((ja+jb)/2,(jc+jd)/2)                    
  if ( abs(hamil_H2cpd(T,J,sha,shb,shc,shd)) < 1.d-12 ) cycle  

  ! Normalization factor
  phasJ = (-1.d0)**(J + 1)                                          
  N_AB = sqrt(1.d0 - delta_ab * phasJ) / (1.d0 + delta_ab)          
  N_CD = sqrt(1.d0 - delta_cd * phasJ) / (1.d0 + delta_cd)          
  if ( (N_AB <= 1.d-12) .or. (N_CD <= 1.d-12) ) cycle                                       
  N_INV = 1.d0 / (N_AB * N_CD)                                      

  call ClebschGordan(ja,jb,2*J,ma,mb,MJ,cg1)                                 
  call ClebschGordan(jc,jd,2*J,mc,md,MJ,cg2)                                 
  sumando = sumando + N_INV * cg1 * cg2 * hamil_H2cpd(T,J,sha,shb,shc,shd)                   
enddo                                                                  

end subroutine uncouple_J

!------------------------------------------------------------------------------!
! subroutine uncouple_JT                                                       !
!                                                                              ! 
! Decouples JT-scheme two-body matrix elements to m-scheme: ab|V|cd>           ! 
! The different isospin components are:                                        ! 
! MMT = 0   T = 0   MT =  0                                                    ! 
! MMT = 1   T = 1   MT = -1                                                    ! 
! MMT = 2   T = 1   MT =  0                                                    ! 
! MMT = 3   T = 1   MT = +1                                                    ! 
! Note that, for now, it is asusmed that there is no-mixing between T=1 MT=0   ! 
! and T=0 MT=0 states.                                                         ! 
!------------------------------------------------------------------------------!
subroutine uncouple_JT (a,b,c,d,sumando)

integer(i32), intent(in) :: a, b, c, d
real(r64), intent(out) :: sumando
integer(i32) :: ja, ma, ta, sha, jb, mb, tb, shb, & 
                jc, mc, tc, shc, jd, md, td, shd, J, MJ, T, MT, MMT
real(r64) :: cg1, cg2, cg3, cg4, N_AB, N_CD, N_INV, delta_ab, delta_cd, phasJT

ja = HOsp_2j(a)
jb = HOsp_2j(b)
jc = HOsp_2j(c)
jd = HOsp_2j(d)

ma = HOsp_2mj(a)
mb = HOsp_2mj(b)
mc = HOsp_2mj(c)
md = HOsp_2mj(d)
MJ = ma + mb

ta = HOsp_2mt(a)
tb = HOsp_2mt(b)
tc = HOsp_2mt(c)
td = HOsp_2mt(d)
MT = ta + tb

sha = HOsp_sh(a)
shb = HOsp_sh(b)
shc = HOsp_sh(c)
shd = HOsp_sh(d)

delta_ab = one * kdelta(sha,shb)
delta_cd = one * kdelta(shc,shd)

sumando = 0.d0                                                            

do J = abs(ma+mb)/2, min((ja+jb)/2,(jc+jd)/2)
  call ClebschGordan(ja,jb,2*J,ma,mb,MJ,cg1)
  call ClebschGordan(jc,jd,2*J,mc,md,MJ,cg2)

  do T = max(abs(ta+tb)/2,0), 1
    MMT = MT/2 + 2*T
    if ( abs(hamil_H2cpd(MMT,J,sha,shb,shc,shd)) < 1.d-12 ) cycle

    ! Normalization factor
    phasJT = (-1.d0)**(J + T)
    N_AB = sqrt(1.d0 - delta_ab * phasJT) / (1.d0 + delta_ab)
    N_CD = sqrt(1.d0 - delta_cd * phasJT) / (1.d0 + delta_cd)
    if ( (N_AB <= 1.d-12) .or. (N_CD <= 1.d-12) ) cycle                                       
    N_INV = 1.d0 / (N_AB * N_CD)

    call ClebschGordan(1,1,2*T,ta,tb,MT,cg3)
    call ClebschGordan(1,1,2*T,tc,td,MT,cg4)
    sumando = sumando + N_INV * cg1 * cg2 * cg3 * cg4 &
                        * hamil_H2cpd(MMT,J,sha,shb,shc,shd)
  enddo 
enddo 

end subroutine uncouple_JT

!------------------------------------------------------------------------------!
! function find_iso                                                            ! 
!                                                                              !
! Determines the correct isospin after symmetry reconstruction. This function  !
! is needed because symmetries are used to reduced the number of matrix eleme- !
! nts stored.                                                                  !
!                                                                              !
! The definition of the permutations is:                                       !
! 1 = abdc                                                                     !
! 2 = bacd                                                                     !
! 3 = badc                                                                     !
! 4 = cdab                                                                     !
! 5 = cdba                                                                     !
! 6 = dcab                                                                     !
! 7 = dcba                                                                     !
!                                                                              !
! The different "isospin" components are:                                      ! 
! T = 0   pp pp                                                                ! 
! T = 1   pn pn                                                                ! 
! T = 2   pn np                                                                ! 
! T = 3   np pn                                                                ! 
! T = 4   np np                                                                ! 
! T = 5   nn nn                                                                ! 
!------------------------------------------------------------------------------!
integer(i32) function find_iso(Tini,perm,htype) result(Tfin)

integer(i32), intent(in) :: Tini, perm, htype

!!! In case of JT-scheme, returns the initial value
if ( htype < 3 ) then
  Tfin = Tini
  return
endif

!!! In case of J-scheme, looks for the correct isospin after permutation
select case (Tini)
  case (0,5)
    Tfin = Tini
  case (1)
    if ( perm == 1 ) Tfin = 2
    if ( perm == 2 ) Tfin = 3        
    if ( perm == 3 ) Tfin = 4         
    if ( perm == 4 ) Tfin = 1         
    if ( perm == 5 ) Tfin = 2         
    if ( perm == 6 ) Tfin = 3         
    if ( perm == 7 ) Tfin = 4        
  case (2)
    if ( perm == 1 ) Tfin = 1
    if ( perm == 2 ) Tfin = 4        
    if ( perm == 3 ) Tfin = 3         
    if ( perm == 4 ) Tfin = 3         
    if ( perm == 5 ) Tfin = 4         
    if ( perm == 6 ) Tfin = 1         
    if ( perm == 7 ) Tfin = 2        
  case (3)
    if ( perm == 1 ) Tfin = 4
    if ( perm == 2 ) Tfin = 1        
    if ( perm == 3 ) Tfin = 2         
    if ( perm == 4 ) Tfin = 2         
    if ( perm == 5 ) Tfin = 1         
    if ( perm == 6 ) Tfin = 4         
    if ( perm == 7 ) Tfin = 3        
  case (4)
    if ( perm == 1 ) Tfin = 3
    if ( perm == 2 ) Tfin = 2        
    if ( perm == 3 ) Tfin = 1         
    if ( perm == 4 ) Tfin = 4         
    if ( perm == 5 ) Tfin = 3         
    if ( perm == 6 ) Tfin = 2         
    if ( perm == 7 ) Tfin = 1        
end select

end function find_iso         

!------------------------------------------------------------------------------!
! subroutine reconstruct_2body_timerev                                         !
!                                                                              !
! Determines the correct order of indices (a,b,c,d) to have the time-reversed  !
! two-body matrix elements in the correct interval (required because of the    !
! symmetry reconstruction) and assigns a permumtation index (the phase is      !
! stored as the sign of the index).                                            !
! To reduce the size of the array, 1-bit integers are used.                    !
!------------------------------------------------------------------------------!
subroutine reconstruct_2body_timerev

integer :: ia, ib, ic, id, ta, tb, tc, td, tmp, ialloc=0
integer(i64) :: kk
real(r64) :: phase

allocate( hamil_trperm(hamil_H2dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of time reversal'

hamil_trperm = 0

do kk = 1, hamil_H2dim
  ia = hamil_abcd(1+4*(kk-1)) 
  ib = hamil_abcd(2+4*(kk-1))
  ic = hamil_abcd(3+4*(kk-1))
  id = hamil_abcd(4+4*(kk-1))
  
  phase = (-one)**((HOsp_2j(ia) + HOsp_2j(ib) + HOsp_2j(ic) + HOsp_2j(id))/2)
  
  ta = HOsp_tr(ia)
  tb = HOsp_tr(ib)
  tc = HOsp_tr(ic)
  td = HOsp_tr(id)
  
  if ( ta > tb ) then
    tmp = ta
    ta = tb
    tb = tmp
    phase = -phase
    hamil_trperm(kk) = hamil_trperm(kk) + int(1,i8)
  endif
    
  if ( tc > td ) then
    tmp = tc
    tc = td
    td = tmp
    phase = -phase
    hamil_trperm(kk) = hamil_trperm(kk) + int(2,i8)
  endif
    
  if ( (ta > tc) .or. ((ta == tc) .and. (tb > td)) ) then 
    hamil_trperm(kk) = hamil_trperm(kk) + int(4,i8)
  endif

  if ( phase < 0.d0 ) hamil_trperm(kk) = hamil_trperm(kk) - int(8,i8)

enddo

end subroutine reconstruct_2body_timerev

!------------------------------------------------------------------------------!
! subroutine find_timerev                                                      !
!                                                                              !
! Given the indices (a,b,c,d) and the predetermined permutation (p), finds     !
! the correct order of the indices for the time-reversed two-body matrix       !
! element.                                                                     !
! This routine depends on the routine reconstruct_2body_timerev.               !
!------------------------------------------------------------------------------!
subroutine find_timerev(p,a,b,c,d)

integer, intent(in) :: p
integer, intent(inout) :: a, b, c, d
integer :: ta, tb, tc, td

ta = HOsp_tr(a)
tb = HOsp_tr(b)
tc = HOsp_tr(c)
td = HOsp_tr(d)

select case (p)
  case(-8,0)
    a = ta 
    b = tb
    c = tc
    d = td

  case(-7,1)
    a = tb 
    b = ta
    c = tc
    d = td

  case(-6,2)
    a = ta 
    b = tb
    c = td
    d = tc

  case(-5,3)
    a = tb 
    b = ta
    c = td
    d = tc

  case(-4,4)
    a = tc 
    b = td
    c = ta
    d = tb

  case(-3,5)
    a = tc 
    b = td
    c = tb
    d = ta

  case(-2,6)
    a = td 
    b = tc
    c = ta
    d = tb

  case(-1,7)
    a = td 
    b = tc
    c = tb
    d = ta
end select

end subroutine find_timerev

!------------------------------------------------------------------------------!
! subroutine calculate_expectval_energy                                        !
!                                                                              !
! Computes the expectation value for the total energy                          !
!   < H > = H_0 + Tr(H_1 * rhoRL) + 1/2 Tr(gammaLR * rhoRL)                    !
!           + 1/2 Tr(deltaLR * kappaRL)                                        !
!------------------------------------------------------------------------------!
subroutine calculate_expectval_energy(rhoLR,kappaRL,gammaLR,deltaLR,energy, &
                                      ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaRL, gammaLR, &
                                                  deltaLR
complex(r64), intent(out) :: energy
integer :: i
complex(r64) :: ener_1b, ener_2b_ph, ener_2b_pp
complex(r64), dimension(ndim,ndim) :: A1, A2, A3

call zgemm('n','n',ndim,ndim,ndim,zone,hamil_H1,ndim,rhoLR,ndim,zzero,A1,ndim)
call zgemm('n','n',ndim,ndim,ndim,zone,gammaLR,ndim,rhoLR,ndim,zzero,A2,ndim)
call zgemm('n','t',ndim,ndim,ndim,zone,deltaLR,ndim,kappaRL,ndim,zzero,A3,ndim)

energy  = zzero
ener_1b = zzero
ener_2b_ph = zzero
ener_2b_pp = zzero

do i = 1, ndim
  ener_1b    = ener_1b    + A1(i,i)
  ener_2b_ph = ener_2b_ph + A2(i,i)
  ener_2b_pp = ener_2b_pp + A3(i,i)
enddo

energy = hamil_H0 + ener_1b + 0.5d0 * (ener_2b_ph + ener_2b_pp)

end subroutine calculate_expectval_energy

!------------------------------------------------------------------------------!
! subroutine calculate_decompo_energy                                          !
!                                                                              !
! Computes the decomposition of the expectation value of the energy, i.e.      !
! the 0-body, 1-body, 2-body (ph and pp).                                      !
!------------------------------------------------------------------------------!
subroutine calculate_decompo_energy(rhoLR,kappaRL,gammaLR,deltaLR,ecomp,ndim)

integer, intent(in) :: ndim
complex(r64), dimension(ndim,ndim), intent(in) :: rhoLR, kappaRL,gammaLR, &
                                                  deltaLR
complex(r64), dimension(12), intent(out) :: ecomp
integer :: i, j, hdim
complex(r64) :: energy_p, energy_n, ener_1b_p, ener_1b_n, &
                ener_2bPH_pp, ener_2bPH_pn, ener_2bPH_np, ener_2bPH_nn, &
                ener_2bPP_pp, ener_2bPP_pn, ener_2bPP_np, ener_2bPP_nn, &
                ener_2bPH_p, ener_2bPH_n, ener_2bPP_p, ener_2bPP_n
complex(r64), dimension(ndim/2,ndim/2) :: rhoLR_pp, rhoLR_pn, &
                                          rhoLR_np, rhoLR_nn, &
                               gammaLR_pp, gammaLR_pn, gammaLR_np, gammaLR_nn, &
                               kappaRL_pp, kappaRL_pn, kappaRL_np, kappaRL_nn, &
                               deltaLR_pp, deltaLR_pn, deltaLR_np, deltaLR_nn, &
                               APH_pp, APH_pn, APH_np, APH_nn, & 
                               APP_pp, APP_pn, APP_np, APP_nn
complex(r64), dimension(ndim,ndim) :: A1, A2, A3

!!! Densities in each block ij (i,j = N or Z)
hdim = ndim/2

do j = 1, hdim
  do i = 1, hdim
    gammaLR_pp(i,j) = gammaLR(i,j)
    gammaLR_pn(i,j) = gammaLR(i,j+hdim)
    gammaLR_np(i,j) = gammaLR(i+hdim,j)
    gammaLR_nn(i,j) = gammaLR(i+hdim,j+hdim)
    
    rhoLR_pp(i,j) = rhoLR(i,j)
    rhoLR_pn(i,j) = rhoLR(i,j+hdim)
    rhoLR_np(i,j) = rhoLR(i+hdim,j)
    rhoLR_nn(i,j) = rhoLR(i+hdim,j+hdim)
    
    deltaLR_pp(i,j) = deltaLR(i,j)
    deltaLR_pn(i,j) = deltaLR(i,j+hdim)
    deltaLR_np(i,j) = deltaLR(i+hdim,j)
    deltaLR_nn(i,j) = deltaLR(i+hdim,j+hdim)
    
    kappaRL_pp(i,j) = kappaRL(i,j)
    kappaRL_pn(i,j) = kappaRL(i,j+hdim)
    kappaRL_np(i,j) = kappaRL(i+hdim,j)
    kappaRL_nn(i,j) = kappaRL(i+hdim,j+hdim)
  enddo 
enddo

!!! Energies in each block ij (i,j = N or Z)
call zgemm('n','n',hdim,hdim,hdim,zone,gammaLR_pp,hdim,rhoLR_pp,hdim,zzero, &
           APH_pp,hdim)  
call zgemm('n','n',hdim,hdim,hdim,zone,gammaLR_pn,hdim,rhoLR_np,hdim,zzero, &
           APH_pn,hdim)  
call zgemm('n','n',hdim,hdim,hdim,zone,gammaLR_np,hdim,rhoLR_pn,hdim,zzero, &
           APH_np,hdim)  
call zgemm('n','n',hdim,hdim,hdim,zone,gammaLR_nn,hdim,rhoLR_nn,hdim,zzero, &
           APH_nn,hdim)  

call zgemm('n','t',hdim,hdim,hdim,zone,deltaLR_pp,hdim,kappaRL_pp,hdim,zzero, &
           APP_pp,hdim)  
call zgemm('n','t',hdim,hdim,hdim,zone,deltaLR_pn,hdim,kappaRL_pn,hdim,zzero, &
           APP_pn,hdim)  
call zgemm('n','t',hdim,hdim,hdim,zone,deltaLR_np,hdim,kappaRL_np,hdim,zzero, &
           APP_np,hdim)  
call zgemm('n','t',hdim,hdim,hdim,zone,deltaLR_nn,hdim,kappaRL_nn,hdim,zzero, &
           APP_nn,hdim)  

ener_2bPH_pp = zzero
ener_2bPH_pn = zzero
ener_2bPH_np = zzero
ener_2bPH_nn = zzero
ener_2bPP_pp = zzero
ener_2bPP_pn = zzero
ener_2bPP_np = zzero
ener_2bPP_nn = zzero

do i = 1, hdim
  ener_2bPH_pp = ener_2bPH_pp + 0.5d0 * APH_pp(i,i)
  ener_2bPH_pn = ener_2bPH_pn + 0.5d0 * APH_pn(i,i)
  ener_2bPH_np = ener_2bPH_np + 0.5d0 * APH_np(i,i)
  ener_2bPH_nn = ener_2bPH_nn + 0.5d0 * APH_nn(i,i)
  ener_2bPP_pp = ener_2bPP_pp + 0.5d0 * APP_pp(i,i)
  ener_2bPP_pn = ener_2bPP_pn + 0.5d0 * APP_pn(i,i)
  ener_2bPP_np = ener_2bPP_np + 0.5d0 * APP_np(i,i)
  ener_2bPP_nn = ener_2bPP_nn + 0.5d0 * APP_nn(i,i)
enddo

!!! Energies in each block i (i = N or Z)       
call zgemm('n','n',ndim,ndim,ndim,zone,hamil_H1,ndim,rhoLR,ndim,zzero,A1,ndim)  
call zgemm('n','n',ndim,ndim,ndim,zone,gammaLR,ndim,rhoLR,ndim,zzero,A2,ndim)  
call zgemm('n','t',ndim,ndim,ndim,zone,deltaLR,ndim,kappaRL,ndim,zzero,A3,ndim) 

energy_p = zzero
energy_n = zzero
ener_1b_p = zzero
ener_1b_n = zzero
ener_2bPH_p = zzero
ener_2bPH_n = zzero
ener_2bPP_p = zzero
ener_2bPP_n = zzero

do i = 1, hdim
  ener_1b_p = ener_1b_p + A1(i,i)
  ener_1b_n = ener_1b_n + A1(i+hdim,i+hdim)
  ener_2bPH_p = ener_2bPH_p  + 0.5d0 * A2(i,i)
  ener_2bPH_n = ener_2bPH_n  + 0.5d0 * A2(i+hdim,i+hdim)
  ener_2bPP_p = ener_2bPP_p  + 0.5d0 * A3(i,i)
  ener_2bPP_n = ener_2bPP_n  + 0.5d0 * A3(i+hdim,i+hdim)
enddo 

energy_p = ener_1b_p + ener_2bPH_p + ener_2bPP_p
energy_n = ener_1b_n + ener_2bPH_n + ener_2bPP_n

!!! Assigns the energies to ouput array       
ecomp(1)  = ener_1b_p   
ecomp(2)  = ener_1b_n   
ecomp(3)  = ener_2bPH_pp
ecomp(4)  = ener_2bPH_pn
ecomp(5)  = ener_2bPH_np
ecomp(6)  = ener_2bPH_nn
ecomp(7)  = ener_2bPP_pp
ecomp(8)  = ener_2bPP_pn
ecomp(9)  = ener_2bPP_np
ecomp(10) = ener_2bPP_nn
ecomp(11) = energy_p   
ecomp(12) = energy_n   

end subroutine calculate_decompo_energy

END MODULE Hamiltonian
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
