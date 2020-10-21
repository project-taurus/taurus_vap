!==============================================================================!
! PROGRAM create_reduced_hamiltonian                                           !
!                                                                              !
! This small program can be used to create a reduced hamiltonian file that     !
! can be used in TAURUS codes. This is particularly useful for large model     !
! spaces. The code includes a simple MPI parallelization to make the           !
! calculation faster (it can be very long).                                    !
! The code works only for J-scheme Hamiltonians.                               !
!                                                                              !
! Note that in the MPI case, the code writes the interaction in many different !
! files that have to be eventuall combined. This can be easily done by using   !
! the commands:  "cat fort.100 >> name_of_hamiltonian.red"                     !
!                "cat fort.101 >> name_of_hamiltonian.red"                     !
!                "cat fort.102 >> name_of_hamiltonian.red"                     !
!                          etc.                                                !
! This is much faster than making the code read/write all the files after the  !
! uncoupling.                                                                  !
!                                                                              !
! To compile: gfortran/ifort -o uncouple.exe create_reduced_hamiltonian.f90    !
!             For mpi compilation, it is necessary to remove the flags "!cmpi" !
!             for example using sed:                                           !
!             sed -i "s/\!cmpi//g" create_reduced_hamiltonian.f90              !
!             mpiifort -o uncouple.exe create_reduced_hamiltonian.f90          !
! To use: ./uncouple.exe < input.txt (sequential)                              !
!         mpirun -np X ./uncouple.exe < input.txt (MPI)                        !
!                                                                              !
! Formpat input file:     (Example)                                            !
! a19, a100                Hamiltonian name = my_mainname_for_taurus_files     !
!==============================================================================!
PROGRAM create_reduced_hamiltonian

use Iso_fortran_env
!cmpi use MPI            

implicit none

!!! Input/output units (for portability)
integer, parameter :: uti = input_unit,  &
                      uto = output_unit, &
                      uth = 95, &       ! hamiltonian (main)
                      uth1 = uth + 1, & !     "       (1b)
                      uth2 = uth + 2, & !     "       (2b)
                      uthc = uth + 4, & !     "       (2b com)
                      uthr = uth + 5    !     "       (red)

!!! Kind parameters (for portability)
integer, parameter :: i8  = int8,   & ! integer  8 bits
                      i16 = int16,  & !    "    16  "   
                      i32 = int32,  & !    "    32  "   
                      i64 = int64,  & !    "    64  "   
                      r32 = real32, & ! real 32 bits (single precision)
                      r64 = real64, & !  "   64  "   (double     "    )
                      rH2 = real32    ! for 2BME of the hamiltonian

!!! Simple names for numerical values
real(r64), parameter :: one  = 1.0d0, &
                        zero = 0.0d0
complex(r64), parameter :: zone  = (1.0d0,0.0d0), &
                           zzero = (0.0d0,0.0d0), &
                           zimag = (0.0d0,1.0d0)

!!! Physical constants
real(r64), parameter :: hbarc = 197.327053d0,     & ! hbar*c in Mev.fm
                        mass_mp = 938.27208816d0, & ! proton mass
                        mass_mn = 939.56542052d0, & ! neutron mass
                        mass_ma = (mass_mp + mass_mn)/2, & ! nucleon mass
                        hbarmass = hbarc**2 / (2*mass_ma)  ! factor kin. energy

!!! Basis
real(r64) :: HO_hw
integer :: HOsh_dim, HOsp_dim, HO_2jmax
integer, dimension(:), allocatable :: HOsh_n, HOsh_l, HOsh_2j, HOsh_na, &
                                      HOsp_n, HOsp_l, HOsp_2j, HOsp_2mj, &
                                      HOsp_2mt, HOsp_sh, HOsp_tr

!!! Hamiltonian
integer :: hamil_type
integer(i64) :: hamil_H2dim, hamil_H2red, hamil_H2cpd_dim
complex(r64) :: hamil_H0
complex(r64), dimension(:,:), allocatable :: hamil_H1
integer(i64), dimension(:,:,:,:), allocatable :: shells_idx
real(r64), dimension(:), allocatable :: hamil_H2cpd
character(len=100) :: hamil_name, hamil_dummy
character(len=:), allocatable :: hamil_fsho, hamil_f01b, hamil_f2b, &
                                 hamil_fcom, hamil_fred, hamil_file
character(4) :: appsho='.sho', app01b='.01b', appcom='.com', appred='.red'
character(3) ::  app2b='.2b'

!!! Various
integer(i8) :: trperm
integer(i16) :: ared, bred, cred, dred
integer :: ht, i, j, t, fac_ht, ialloc=0, &
           a, ja, ma, la, ta, b, jb, mb, lb, tb, bmax, tmax, &
           c, jc, mc, lc, tc, d, jd, md, ld, td, &
           uth6, uth7, uth8
real(r64) :: Vtmp, Vcut, Vdec
real(rH2) :: Vred
character(19), dimension(3) :: inp_ch
character(11) :: file_uth6, file_uth7, file_uth8
character(6) ::  uth6_ch, uth7_ch, uth8_ch

!!! Parallel variables
integer :: worldsize=1, myrank=0, myjump, ierr

!------------------------------------------------------------------------------!
! Initialization                                                               !
!------------------------------------------------------------------------------!

!cmpi call mpi_init(ierr)
!cmpi call mpi_comm_size(MPI_COMM_WORLD,worldsize,ierr)
!cmpi call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)

if ( myrank == 0 ) then
  read(uti, '(a19, a100)') inp_ch(1), hamil_dummy
endif

!cmpi call mpi_bcast(hamil_dummy,100,mpi_character,0,mpi_comm_world,ierr)

!!! Assigns the names for the hamiltonian files and opens them
j = len_trim(adjustl(hamil_dummy))

allocate (character(j  ) :: hamil_file)
allocate (character(j+4) :: hamil_fsho)
allocate (character(j+4) :: hamil_f01b)
allocate (character(j+3) :: hamil_f2b )
allocate (character(j+4) :: hamil_fred)
allocate (character(j+4) :: hamil_fcom)

hamil_file = trim(adjustl(hamil_dummy))
hamil_fsho = hamil_file // appsho
hamil_f01b = hamil_file // app01b
hamil_f2b  = hamil_file // app2b
hamil_fred = hamil_file // appred
hamil_fcom = hamil_file // appcom

!!! No check is performed on the files, so it better be correct
open(uth,  file=hamil_fsho, status='old', action='read', form='formatted')
open(uth1, file=hamil_f01b, status='old', action='read', form='formatted')
open(uth2, file=hamil_f2b,  status='old', action='read', form='formatted')

open(uthr, file=hamil_fred, status='replace', access='stream', action='write', &
     form='unformatted')

!------------------------------------------------------------------------------!
! Defines the basis                                                            !
!------------------------------------------------------------------------------!

call set_basis

!------------------------------------------------------------------------------!
! Reads/writes the zero-one body parts                                         !
!------------------------------------------------------------------------------!

call read_hamiltonian_1body_general

if ( myrank == 0 ) then
  write(uthr) hamil_name 
  write(uthr) hamil_type 
  write(uthr) hamil_H0      
  write(uthr) ((hamil_H1(j,i), j=1,HOsp_dim), i=1,HOsp_dim)
endif 

!------------------------------------------------------------------------------!
! Reads/uncouples/writes the two-body parts                                    !
!------------------------------------------------------------------------------!
    
allocate( shells_idx(HOsh_dim,HOsh_dim,HOsh_dim,HOsh_dim), stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of shell index'
shells_idx = 0

hamil_H2cpd_dim = 0
call read_hamiltonian_2body_Jscheme(0)

allocate( hamil_H2cpd(hamil_H2cpd_dim), stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of hamil_H2cpd'
hamil_H2cpd = zero

rewind(uth2)
call read_hamiltonian_2body_Jscheme(1)

close (uth2, status='keep')

!------------------------------------------------------------------------------!
! Uncouples/writes the two-body parts                                          !
!------------------------------------------------------------------------------!

uth6 = uthr + 1 + myrank
uth7 = uth6 + worldsize
uth8 = uth6 + 2*worldsize

write(uth6_ch,'(1i6)') uth6
write(uth7_ch,'(1i6)') uth7
write(uth8_ch,'(1i6)') uth8

uth6_ch = adjustl(uth6_ch)
uth7_ch = adjustl(uth7_ch)
uth8_ch = adjustl(uth8_ch)

file_uth6 = 'fort.' // uth6_ch
file_uth7 = 'fort.' // uth7_ch
file_uth8 = 'fort.' // uth8_ch

open(uth6, file=file_uth6, status='replace', action='write', access='stream', & 
     form='unformatted')  
open(uth7, file=file_uth7, status='replace', action='write', access='stream', & 
     form='unformatted')  
open(uth8, file=file_uth8, status='replace', action='write', access='stream', & 
     form='unformatted')  
  
Vcut = 1.0d-16
hamil_H2dim = 0
hamil_H2red = 0

!cmpi myjump = 0          
do a = 1, HOsp_dim
  !cmpi if ( a /= 1 + myrank + myjump ) cycle
  !cmpi myjump = myjump + worldsize
  la = HOsp_l(a)
  ma = HOsp_2mj(a)
  ta = HOsp_2mt(a)
  do c = a, HOsp_dim
    lc = HOsp_l(c)
    mc = HOsp_2mj(c)
    tc = HOsp_2mt(c)
    do d = c, HOsp_dim
      ld = HOsp_l(d)
      md = HOsp_2mj(d)
      td = HOsp_2mt(d)
      bmax = HOsp_dim
      if ( c == a ) bmax = d
      do b = a, bmax
        lb = HOsp_l(b)
        mb = HOsp_2mj(b)
        tb = HOsp_2mt(b)

        if ( (-1)**(la+lb) /= (-1)**(lc+ld) ) cycle
        if ( ta + tb /= tc + td ) cycle
        if ( ma + mb /= mc + md ) cycle
        if ( ma + mb < 0 ) cycle
 
        call uncouple_J(a,b,c,d,Vdec) 
    
        !!! Select only matrix elements above a given cutoff to reduce the
        !!! CPU time and storage
        if ( abs(Vdec) > Vcut ) then 
          call reconstruct_2body_timerev(a,b,c,d,trperm)
          hamil_H2dim = hamil_H2dim + 1
          ared = a
          bred = b
          cred = c
          dred = d
          Vred = Vdec
          write(uth6) ared, bred, cred, dred      
          write(uth7) Vred      
          write(uth8) trperm    
        endif  
         
      enddo  !end loop b
    enddo  !end loop d
  enddo  !end loop c
enddo  !end loop a
  
deallocate (hamil_H2cpd)

close(uth6, status='keep')
close(uth7, status='keep')
close(uth8, status='keep')

!cmpi call mpi_reduce(hamil_H2dim,hamil_H2red,1,mpi_integer8,mpi_sum,0,&
!cmpi                 mpi_comm_world,ierr)
!cmpi hamil_H2dim = hamil_H2red

if ( myrank == 0 ) then
! print*, hamil_H2dim
  write(uthr) hamil_H2dim   
endif

!cmpi call mpi_finalize(ierr)

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_basis                                                         !
!                                                                              !
! Sets HO basis: number of levels, quantum numbers, oscillator parameters, ... !
! To determine the basis, the codes reads the shells (same number for protons  !
! and neutrons) that defines the model space in the main hamiltonian file.     !
!------------------------------------------------------------------------------!
subroutine set_basis    

integer :: i, j, k, facn, shinc, mjinc, jmax, ialloc=0

!!! Recovers the basis information from hamiltonian file 
rewind(uth)
read(uth,*) hamil_name
read(uth,*) hamil_type
read(uth,*) HOsh_dim

allocate( HOsh_n(HOsh_dim), HOsh_l(HOsh_dim), HOsh_2j(HOsh_dim), &
          HOsh_na(HOsh_dim), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of shells'

backspace uth

read(uth,*) HOsh_dim, (HOsh_na(i),i=1,HOsh_dim)

!!! Determines the quantum numbers of the shells
facn = 10000

do i = 1, HOsh_dim
  HOsh_n(i)  =  HOsh_na(i) / facn 
  HOsh_l(i)  = (HOsh_na(i) - HOsh_n(i)*facn) / 100
  HOsh_2j(i) =  HOsh_na(i) - HOsh_n(i)*facn - HOsh_l(i)*100
enddo

!!! Computes the maximum values reachable in the basis
HO_2jmax = maxval(HOsh_2j)

!!! Determines the dimension of the basis and check against the particle numbers
HOsp_dim = 0
do i = 1, HOsh_dim
  HOsp_dim = HOsp_dim + HOsh_2j(i) + 1
enddo
HOsp_dim = 2 * HOsp_dim 

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

end subroutine set_basis    

!------------------------------------------------------------------------------!
! subroutine read_hamiltonian_1body_general                                    !
!                                                                              !
! Reads the one-body matrix elements of the hamiltonian in the general format. !
!------------------------------------------------------------------------------!
subroutine read_hamiltonian_1body_general
  
integer(i32) :: k, l, mjk, mjl, a, b, iexit, ialloc=0
real(r64) :: H0, hamil_H1cpd_p, hamil_H1cpd_n

allocate( hamil_H1(HOsp_dim,HOsp_dim), stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of coupled hamiltonian'

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
 
end subroutine read_hamiltonian_1body_general

!------------------------------------------------------------------------------!
! subroutine read_hamiltonian_2body_Jscheme                                    !
!                                                                              !
! Reads the two-body matrix elements of the hamiltonian in J-scheme.           !
!                                                                              !
! opt = 0 : first reading to determine the number of matrix elements           !
!     = 1 : second reading to actually the matrix elements in a vector         !
!------------------------------------------------------------------------------!
subroutine read_hamiltonian_2body_Jscheme(opt)
 
integer, intent(in) :: opt 
integer(i32) :: i, j, T, tmin, tmax, jmin, jmax, iexit, a, b, c, d, nlj_tmp(4)
integer(i64) :: k
  
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

  if ( opt == 0 ) then 

    shells_idx(a,b,c,d) = 1 + hamil_H2cpd_dim
    hamil_H2cpd_dim = hamil_H2cpd_dim + (jmax-jmin+1)*6
    do j = jmin, jmax
      read(uth2,*) 
    enddo 

  else

    do j = 0, jmax-jmin
      k = shells_idx(a,b,c,d)
      read(uth2,*) (hamil_H2cpd(k + T + 6*j), T=tmin,tmax)
    enddo 

  endif
enddo
  
end subroutine read_hamiltonian_2body_Jscheme
  
!------------------------------------------------------------------------------!
! subroutine reconstruct_2body_timerev                                         !
!                                                                              !
! Determines the correct order of indices (a,b,c,d) to have the time-reversed  !
! two-body matrix elements in the correct interval (required because of the    !
! symmetry reconstruction) and assigns a permumtation index (the phase is      !
! stored as the sign of the index).                                            !
! To reduce the size of the array, 1-bit integers are used.                    !
!------------------------------------------------------------------------------!
subroutine reconstruct_2body_timerev(ia,ib,ic,id,trperm)

integer, intent(in) :: ia, ib, ic, id
integer(i8), intent(out) :: trperm
integer :: ta, tb, tc, td, tmp
real(r64) :: phase

trperm = 0

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
  trperm = trperm + 1
endif
  
if ( tc > td ) then
  tmp = tc
  tc = td
  td = tmp
  phase = -phase
  trperm = trperm + 2
endif
  
if ( (ta > tc) .or. ((ta == tc) .and. (tb > td)) ) then 
  trperm = trperm + 4
endif

if ( phase < 0.d0 ) trperm = trperm - 8

end subroutine reconstruct_2body_timerev
  
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
                jc, mc, tc, shc, jd, md, td, shd, J, MJ, T, &
                jinc, perm, tsym
integer(i64) :: ksym
real(r64) :: cg1, cg2, N_AB, N_CD, N_INV, delta_ab, delta_cd, phasJ, HJ, &
             xja, xjb, xjc, xjd, xj, xt, phasab, phascd, phas

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

xja = ja / 2.d0
xjb = jb / 2.d0
xjc = jc / 2.d0
xjd = jd / 2.d0

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

!!! Needed to determine what is needed.
if ( shells_idx(sha,shb,shc,shd) /= 0 ) then
  ksym = shells_idx(sha,shb,shc,shd)
  perm = 0
  tsym = T
elseif (shells_idx(sha,shb,shd,shc) /= 0 ) then
  ksym = shells_idx(sha,shb,shd,shc)
  perm = 1
  tsym = find_osi(T,1,hamil_type)
elseif (shells_idx(shb,sha,shc,shd) /= 0 ) then
  ksym = shells_idx(shb,sha,shc,shd)
  perm = 2
  tsym = find_osi(T,2,hamil_type)
elseif (shells_idx(shb,sha,shd,shc) /= 0 ) then
  ksym = shells_idx(shb,sha,shd,shc)
  perm = 3
  tsym = find_osi(T,3,hamil_type)
elseif (shells_idx(shc,shd,sha,shb) /= 0 ) then
  ksym = shells_idx(shc,shd,sha,shb)
  perm = 4
  tsym = find_osi(T,4,hamil_type)
elseif (shells_idx(shc,shd,shb,sha) /= 0 ) then
  ksym = shells_idx(shc,shd,shb,sha)
  perm = 6
  tsym = find_osi(T,6,hamil_type)
elseif (shells_idx(shd,shc,sha,shb) /= 0 ) then
  ksym = shells_idx(shd,shc,sha,shb)
  perm = 5
  tsym = find_osi(T,5,hamil_type)
elseif (shells_idx(shd,shc,shb,sha) /= 0 ) then
  ksym = shells_idx(shd,shc,shb,sha)
  perm = 7
  tsym = find_osi(T,7,hamil_type)
endif

!!! The sum over the J values
sumando = 0.d0                                                            
jinc = 0 + abs(ma+mb)/2 - max(abs(ja-jb)/2,abs(jc-jd)/2) - 1

do J = abs(ma+mb)/2, min((ja+jb)/2,(jc+jd)/2)                    
  jinc = jinc + 1

  xj = 1.d0 * J      
  xt = 1.d0
  phasab = (-1.d0)**(xja + xjb + xj + xt)
  phascd = (-1.d0)**(xjc + xjd + xj + xt)
  if ( perm == 0 ) phas = 1.0d0 
  if ( perm == 1 ) phas = 1.0d0 * phascd
  if ( perm == 2 ) phas = 1.0d0 * phasab
  if ( perm == 3 ) phas = 1.0d0 * phasab * phascd
  if ( perm == 4 ) phas = 1.0d0
  if ( perm == 6 ) phas = 1.0d0 * phasab
  if ( perm == 5 ) phas = 1.0d0 * phascd
  if ( perm == 7 ) phas = 1.0d0 * phasab * phascd

  HJ = phas * hamil_H2cpd( ksym + tsym + 6*jinc )
  if ( abs(HJ) < 1.d-12 ) cycle  

  ! Normalization factor
  phasJ = (-1.d0)**(J + 1)                                          
  N_AB = sqrt(1.d0 - delta_ab * phasJ) / (1.d0 + delta_ab)          
  N_CD = sqrt(1.d0 - delta_cd * phasJ) / (1.d0 + delta_cd)          
  if ( (N_AB <= 1.d-12) .or. (N_CD <= 1.d-12) ) cycle                                       
  N_INV = 1.d0 / (N_AB * N_CD)                                      

  call ClebschGordan(ja,jb,2*J,ma,mb,MJ,cg1)                                 
  call ClebschGordan(jc,jd,2*J,mc,md,MJ,cg2)                                 
  sumando = sumando + N_INV * cg1 * cg2 * HJ
enddo                                                                  

end subroutine uncouple_J

!------------------------------------------------------------------------------!
! subroutine ClebschGordan                                                     !
!                                                                              !
! Computes the ClebschGordan for the group SU(2). The algorithm used was taken !
! from technical notes from NASA written by W. F. Ford and R. C. Bruley.       !
! Ref: NASA TN D-6173                                                          !
!                                                                              !
! (j1,m1,j2,m2|j3,m3) = c * g                                                  !
! with                                                                         !
! c = D(j1j2j3) * [(j1-m1)!(j2+m2)!(j1+m1)!(j2-m2)!(j3+m3)!(j3-m3)!]**1/2      !
! g = sqrt(2*j3+1) sum_l (-1)**l [(j1+j2-j3-l)!(j1-m1-l)!(j2+m2-l)!            !
!                                 (j3-j2+m1+l)!(j3-j1-m1+l)!l!]**-1            !
! D(j1j2j3) = [(j1+j2-j3)!(j2+j3-j1)!(j3+j1-j2)!]**1/2 / [(j1+j2+j3+1)!]**1/2  !
!                                                                              !
! Be aware that all input values of j and m are supposed to be twice their     !
! values (such that we can treat easily half-integer angular momenta).         !
!------------------------------------------------------------------------------!
subroutine ClebschGordan (j1,j2,j3,m1,m2,m3,cg)

integer, intent(in) :: j1, j2, j3, m1, m2, m3
real(r64), intent(out) :: cg 
integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, k1, k2, k3, k4, k5, k6, &
           l, l1, l2
real(r64) :: c, g, p, q, h, hl, hlm1

cg = zero

!!! Computes the factor c
n1 = 1 + (j1 + j2 - j3)/2
n2 = 1 + (j2 + j3 - j1)/2
n3 = 1 + (j3 + j1 - j2)/2
n4 = 1 + (j1 - m1)/2
n5 = 1 + (j2 + m2)/2
n6 = 1 + (j1 + m1)/2
n7 = 1 + (j2 - m2)/2
n8 = 1 + (j3 + m3)/2
n9 = 1 + (j3 - m3)/2
n10= n1 + n2 + n3 - 1

if ( (min(n1,n2,n3,n4,n5,n6,n7,n8,n9) < 1) .or. (m1+m2 /= m3) ) return

p =  log_gamma(n1+zero) + log_gamma(n2+zero) + log_gamma(n3+zero) &
   + log_gamma(n4+zero) + log_gamma(n5+zero) + log_gamma(n6+zero) &
   + log_gamma(n7+zero) + log_gamma(n8+zero) + log_gamma(n9+zero) &
   - log_gamma(n10+zero)

c = exp(0.5d0*p)

!!! Computes the factor g
k1 = n1
k2 = n4
k3 = n5
k4 = n4 - n3
k5 = n5 - n2

l1 = max(0,k4,k5)
l2 = min(k1,k2,k3)

h  = one 
hl = one

do l = l1+1, l2
  hlm1 = hl
  hl = hlm1 * (l - k1) * (l - k2) * (l - k3) / ((l - k4) * (l - k5) * l)
  h = h + hl
enddo

k1 = k1 - l1
k2 = k2 - l1
k3 = k3 - l1
k4 = l1 + 1 - k4 
k5 = l1 + 1 - k5 
k6 = l1 + 1 

q =  log_gamma(k1+zero) + log_gamma(k2+zero) + log_gamma(k3+zero) &
   + log_gamma(k4+zero) + log_gamma(k5+zero) + log_gamma(k6+zero)

g = sqrt(j3 + one) * (-1)**l1 * exp(-q) * h

!!! Computes the final value combining the two parts.
cg = c * g

end subroutine ClebschGordan

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
! function find_iso                                                            ! 
!                                                                              !
! Reverse of find_iso (the two function could actually be combined at some     !
! point).                                                                      !
!------------------------------------------------------------------------------!
integer(i32) function find_osi(Tini,perm,htype) result(Tfin)

implicit none
integer(i32) :: Tini, perm, htype

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
    if ( perm == 5 ) Tfin = 3         
    if ( perm == 6 ) Tfin = 2         
    if ( perm == 7 ) Tfin = 4        
  case (2)
    if ( perm == 1 ) Tfin = 1
    if ( perm == 2 ) Tfin = 4        
    if ( perm == 3 ) Tfin = 3         
    if ( perm == 4 ) Tfin = 3         
    if ( perm == 5 ) Tfin = 1         
    if ( perm == 6 ) Tfin = 4         
    if ( perm == 7 ) Tfin = 2        
  case (3)
    if ( perm == 1 ) Tfin = 4
    if ( perm == 2 ) Tfin = 1        
    if ( perm == 3 ) Tfin = 2         
    if ( perm == 4 ) Tfin = 2         
    if ( perm == 5 ) Tfin = 4         
    if ( perm == 6 ) Tfin = 1         
    if ( perm == 7 ) Tfin = 3        
  case (4)
    if ( perm == 1 ) Tfin = 3
    if ( perm == 2 ) Tfin = 2        
    if ( perm == 3 ) Tfin = 1         
    if ( perm == 4 ) Tfin = 4         
    if ( perm == 5 ) Tfin = 2         
    if ( perm == 6 ) Tfin = 3         
    if ( perm == 7 ) Tfin = 1        
end select

end function find_osi         

!------------------------------------------------------------------------------!
! function kdelta                                                              !
!                                                                              !
! Computes the Kronecker delta: \delta_ij = 1 if i = j                         !
!                                         = 0 otherwise                        !
!------------------------------------------------------------------------------!
function kdelta(i,j) result(delta)

integer, intent(in) :: i, j
integer :: delta

if ( i == j ) then 
  delta = 1    
else
  delta = 0    
endif

end function kdelta

END PROGRAM create_reduced_hamiltonian
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
