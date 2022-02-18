!==============================================================================!
! MODULE Initialization                                                        !
!                                                                              !
! This module contains the variables and routines related to the reading of    !
! the input parameters and files.                                              !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine print_version                                                   !
! - subroutine read_input                                                      !
! - subroutine print_input                                                     !
! - subroutine check_input                                                     !
! - subroutine broadcast_inputs (MPI)                                          !
! - subroutine open_files_hamiltonian                                          !
!==============================================================================!
MODULE Initialization

use Constants
!cmpi use MPI            
!cmpi use Parallelization
use Nucleus, only: valence_Z, valence_N
use Hamiltonian, only: hamil_file, hamil_fsho, hamil_f01b, hamil_f2b, &
                       hamil_fred, hamil_fcom, hamil_com, hamil_type, &
                       hamil_read 
use WaveFunctions, only: seed_type, blocking_dim, blocking_id, seed_text, &
                         seed_rand, seed_symm, seed_occeps, seed_allemp, &
                         dens_spatial, dens_nr, dens_dr
use Pairs, only: pairs_scheme
use Projection, only: proj_Mphip, proj_Mphin
use Constraints, only: constraint_eps, constraint_max, constraint_dim, &
                       constraint_read, constraint_switch, opt_betalm, & 
                       constraint_types, enforce_NZ
use Gradient, only: gradient_type, gradient_eta, gradient_mu, gradient_eps

implicit none
private 

character(30), dimension(69) :: input_names ! Name of inputs
character(30), dimension(:), allocatable :: input_block ! Name for blocking

!!! Public routines
public :: print_version, read_input

CONTAINS 

!------------------------------------------------------------------------------!
! subroutine print_version                                                     !
!                                                                              !
! Prints the panel with the version, the description and the logo of the code. !
! This is the first action of the code.                                        !
!------------------------------------------------------------------------------!
subroutine print_version

!cmpi if ( paral_myrank == 0 ) then        
print '("  _________________________________________________________ ",/, &
      & " |                                                         |",/, &
      & " |  (______)  TAURUS_vap                Benjamin Bally     |",/, &
      & " |  <(0  0)>  2022.02.18                Tomás R. Rodríguez |",/, &
      & " |    (°°)                              Adrián Sánchez-F.  |",/, &
      & " |                                                         |",/, &
      & " | This code performs the particle-number variation after  |",/, &
      & " | projection of real general Bogoliubov quasi-particle    |",/, &
      & " | states represented in a spherical harmonic oscillator   |",/, &
      & " | basis.                                                  |",/, &
      & " |                                                         |",/, &
      & " | Licence: GNU General Public License version 3 or later  |",/, &
      & " | DOI: https://doi.org/10.5281/zenodo.4130680             |",/, &
      & " | Git: https://github.com/project-taurus/taurus_vap.git   |",/, &
      & " |_________________________________________________________|",/)' 
!cmpi endif

end subroutine print_version

!------------------------------------------------------------------------------!
! subroutine read_input                                                        !
!                                                                              !
! Reads the input parameters from the input file and perform somes checks to   !
! see if they have approriate values. Also, opens the files containing the     !
! information on the Hamiltonian and model space.                              !
!                                                                              !
! Input: iter_max = maximum number of iterations                               !
!        iter_write = step for the intermediate writing of wave functions      !
!        iter_print = option to calculate and print more at each iteration     !
!------------------------------------------------------------------------------!
subroutine read_input(iter_max,iter_write,iter_print)

integer, intent(out) :: iter_max, iter_write, iter_print
integer :: i, j, dummy_teamssize
character(3) :: app2b='.2b'
character(4) :: appsho='.sho', app01b='.01b', appred='.red', appcom='.com'
character(100) :: hamil_dummy
character(len=*), parameter :: format1 = "(1a)", &
                               format2 = "(1a30,1a100)", &
                               format3 = "(1a30,1i1)", &
                               format4 = "(1a30,1i5)", &
                               format5 = "(1a30,1f7.2)", &
                               format6 = "(1a30,1es10.3)", &
                               format7 = "(1a30,1i1,1x,1f7.3,1x,1f7.3)", &
                               format8 = "(1a30,1i3,1x,1f5.2)"

!!! Reads the input parameters
!cmpi if ( paral_myrank == 0 ) then        
read(uti,format1) input_names(1)
read(uti,format1) input_names(2)
read(uti,format2) input_names(3),  hamil_dummy
read(uti,format3) input_names(4),  hamil_com  
read(uti,format3) input_names(5),  hamil_read
read(uti,format4) input_names(6),  dummy_teamssize
read(uti,format1) input_names(7) 
read(uti,format1) input_names(8) 
read(uti,format1) input_names(9) 
read(uti,format5) input_names(10), valence_Z
read(uti,format5) input_names(11), valence_N
read(uti,format4) input_names(12), proj_Mphip
read(uti,format4) input_names(13), proj_Mphin
read(uti,format1) input_names(14)
read(uti,format1) input_names(15)
read(uti,format1) input_names(16)
read(uti,format3) input_names(17), seed_type
read(uti,format4) input_names(18), blocking_dim
allocate( blocking_id(blocking_dim), input_block(blocking_dim) )
do i = 1, blocking_dim
  read(uti,format4) input_block(i), blocking_id(i)
enddo
read(uti,format3) input_names(19), seed_symm    
read(uti,format4) input_names(20), seed_rand
read(uti,format3) input_names(21), seed_text    
read(uti,format6) input_names(22), seed_occeps   
read(uti,format3) input_names(23), seed_allemp   
read(uti,format3) input_names(24), dens_spatial  
read(uti,format8) input_names(25), dens_nr(1), dens_dr(1) 
read(uti,format8) input_names(26), dens_nr(2), dens_dr(2) 
read(uti,format8) input_names(27), dens_nr(3), dens_dr(3) 
read(uti,format1) input_names(28)
read(uti,format1) input_names(29)
read(uti,format1) input_names(30)
read(uti,format4) input_names(31), iter_max
read(uti,format4) input_names(32), iter_write
read(uti,format3) input_names(33), iter_print
read(uti,format3) input_names(34), gradient_type
read(uti,format6) input_names(35), gradient_eta
read(uti,format6) input_names(36), gradient_mu
read(uti,format6) input_names(37), gradient_eps
read(uti,format1) input_names(38)
read(uti,format1) input_names(39)
read(uti,format1) input_names(40)
read(uti,format3) input_names(41), enforce_NZ
read(uti,format3) input_names(42), opt_betalm
read(uti,format3) input_names(43), pairs_scheme  
read(uti,format6) input_names(44), constraint_eps
do i = 3, constraint_types  
  read(uti,format7) input_names(42+i), constraint_switch(i), &
                    constraint_read(i,1), constraint_read(i,2)
enddo
!cmpi endif

!!! MPI adjustment and broadcast
!cmpi paral_teamssize = max(dummy_teamssize,1) 
!cmpi call broadcast_inputs(hamil_dummy,iter_max,iter_write,iter_print)

!!! Counts the number of constraints. By default, switches on the constraint for
!!! the particle numbers. The constraint on delta is treated separately.
constraint_switch(1) = 1
constraint_read(1,1) = valence_Z
constraint_switch(2) = 1
constraint_read(2,1) = valence_N

constraint_dim = 0

do i = 1, constraint_types-1
  constraint_dim = constraint_dim + min(constraint_switch(i),2)
enddo

where (constraint_switch == 0) constraint_read(:,1) = 0.d0
where (constraint_switch <= 1) constraint_read(:,2) = 0.d0

!!! Assigns the value 1 to the number of gauge angles if 0 is read
if ( proj_Mphip > -1 ) proj_Mphip = max(1,proj_Mphip)
if ( proj_Mphin > -1 ) proj_Mphin = max(1,proj_Mphin)

!!! Determines the spacing for the spatial density in spherical coordinates
if ( dens_spatial == 0 ) then
  dens_nr = 0   
  dens_dr = zero
elseif ( dens_spatial == 1 ) then
  dens_nr(2:3) = 0   
  dens_dr(2:3) = zero
endif

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

!!! Performs some tests on the value of the inputs and link the units for
!!! the Hamiltonian files
!cmpi if ( paral_myrank == 0 ) then        
call print_input(iter_max,iter_write,iter_print)
!cmpi endif
call check_input(iter_max,iter_write,iter_print)
call open_files_hamiltonian

end subroutine read_input

!------------------------------------------------------------------------------!
! subroutine print_input                                                       !
!                                                                              !
! Prints the input parameters at the beginning of the calculation in the same  !
! format such that it can be copied and reused in an input file.               !
!                                                                              !
! Input: iter_max = maximum number of iterations                               !
!        iter_write = step for the intermediate writing of wave functions      !
!        iter_print = option to calculate and print more at each iteration     !
!------------------------------------------------------------------------------!
subroutine print_input(iter_max,iter_write,iter_print)

integer, intent(in) :: iter_max, iter_write, iter_print
integer :: i, dummy_teamssize=0
character(3) :: dens_nr_ch(3)
character(5) :: dens_dr_ch(3)
character(5) :: blocking_dim_ch, iter_max_ch, iter_write_ch, proj_Mphip_ch, &
                proj_Mphin_ch, paral_teamssize_ch, seed_rand_ch
character(10) :: gradient_eta_ch, gradient_mu_ch, gradient_eps_ch, &
                 constraint_eps_ch, valence_N_ch, valence_Z_ch, &
                 seed_occeps_ch
character(7) :: constraint_read_ch(constraint_types,2)
character(5), dimension(:), allocatable :: blocking_id_ch
character(len=*), parameter :: format1 = "(1a)", &
                               format2 = "(1a30,1a)", &
                               format3 = "(1a30,1i1)", &
                               format4 = "(1a30,1a5)", &
                               format5 = "(1a30,1a10)", &
                               format6 = "(1a30,1i1,1x,1a7,1x,1a7)", &
                               format7 = "(1a30,1a3,1x,1a5)"

allocate(blocking_id_ch(blocking_dim))

!!! Formats the variable to eliminate the unpleasant blanck spaces
write(valence_Z_ch,'(1f7.2)') valence_Z
write(valence_N_ch,'(1f7.2)') valence_N
valence_Z_ch = adjustl(valence_Z_ch)
valence_N_ch = adjustl(valence_N_ch)

!cmpi dummy_teamssize = paral_teamssize 
write(paral_teamssize_ch,'(1i5)') dummy_teamssize
paral_teamssize_ch = adjustl(paral_teamssize_ch)

write(proj_Mphip_ch,'(1i5)') proj_Mphip
write(proj_Mphin_ch,'(1i5)') proj_Mphin
proj_Mphip_ch = adjustl(proj_Mphip_ch)
proj_Mphin_ch = adjustl(proj_Mphin_ch)

write(blocking_dim_ch,'(1i5)') blocking_dim
blocking_dim_ch = adjustl(blocking_dim_ch)
do i = 1, blocking_dim
  write(blocking_id_ch(i),'(1i5)') blocking_id(i)
  blocking_id_ch(i) = adjustl(blocking_id_ch(i))
enddo

write(seed_rand_ch,'(1i5)') seed_rand
write(seed_occeps_ch,'(1es10.3)') seed_occeps
seed_rand_ch = adjustl(seed_rand_ch)
seed_occeps_ch = adjustl(seed_occeps_ch)

write(dens_nr_ch(1),'(1i3)') dens_nr(1)
write(dens_nr_ch(2),'(1i3)') dens_nr(2)
write(dens_nr_ch(3),'(1i3)') dens_nr(3)
dens_nr_ch(1) = adjustl(dens_nr_ch(1))
dens_nr_ch(2) = adjustl(dens_nr_ch(2))
dens_nr_ch(3) = adjustl(dens_nr_ch(3))

write(dens_dr_ch(1),'(1f5.2)') dens_dr(1)
write(dens_dr_ch(2),'(1f5.2)') dens_dr(2)
write(dens_dr_ch(3),'(1f5.2)') dens_dr(3)
dens_dr_ch(1) = adjustl(dens_dr_ch(1))
dens_dr_ch(2) = adjustl(dens_dr_ch(2))
dens_dr_ch(3) = adjustl(dens_dr_ch(3))

write(iter_max_ch,'(1i5)') iter_max
write(iter_write_ch,'(1i5)') iter_write
iter_max_ch = adjustl(iter_max_ch)
iter_write_ch = adjustl(iter_write_ch)

write(gradient_eta_ch,'(1es10.3)') gradient_eta
write(gradient_mu_ch,'(1es10.3)') gradient_mu
write(gradient_eps_ch,'(1es10.3)') gradient_eps
gradient_eta_ch = adjustl(gradient_eta_ch)
gradient_mu_ch = adjustl(gradient_mu_ch)
gradient_eps_ch = adjustl(gradient_eps_ch)

write(constraint_eps_ch,'(1es10.3)') constraint_eps
constraint_eps_ch = adjustl(constraint_eps_ch)

do i = 3, constraint_types
  write(constraint_read_ch(i,1),'(1f7.3)') constraint_read(i,1)
  write(constraint_read_ch(i,2),'(1f7.3)') constraint_read(i,2)
  constraint_read_ch(i,1) = adjustr(constraint_read_ch(i,1))
  constraint_read_ch(i,2) = adjustr(constraint_read_ch(i,2))
enddo


!!! Prints the input parameters
print '(60("%"),/,22x,"INPUT PARAMETERS",22x,/,60("%"),/)'
write(uto,format1) input_names(1)
write(uto,format1) input_names(2)
write(uto,format2) input_names(3),  trim(adjustl(hamil_file))
write(uto,format3) input_names(4),  hamil_com  
write(uto,format3) input_names(5),  hamil_read
write(uto,format4) input_names(6),  paral_teamssize_ch
write(uto,format1) input_names(7) 
write(uto,format1) input_names(8) 
write(uto,format1) input_names(9) 
write(uto,format5) input_names(10), valence_Z_ch
write(uto,format5) input_names(11), valence_N_ch
write(uto,format4) input_names(12), proj_Mphip_ch
write(uto,format4) input_names(13), proj_Mphin_ch
write(uto,format1) input_names(14)
write(uto,format1) input_names(15)
write(uto,format1) input_names(16)
write(uto,format3) input_names(17), seed_type
write(uto,format4) input_names(18), blocking_dim_ch
do i = 1, blocking_dim
  write(uto,format4) input_block(i), blocking_id_ch(i)
enddo
write(uto,format3) input_names(19), seed_symm    
write(uto,format4) input_names(20), seed_rand_ch
write(uto,format3) input_names(21), seed_text    
write(uto,format5) input_names(22), seed_occeps_ch
write(uto,format3) input_names(23), seed_allemp  
write(uto,format3) input_names(24), dens_spatial 
write(uto,format7) input_names(25), dens_nr_ch(1), dens_dr_ch(1)
write(uto,format7) input_names(26), dens_nr_ch(2), dens_dr_ch(2)
write(uto,format7) input_names(27), dens_nr_ch(3), dens_dr_ch(3)
write(uto,format1) input_names(28)
write(uto,format1) input_names(29)
write(uto,format1) input_names(30)
write(uto,format4) input_names(31), iter_max_ch
write(uto,format4) input_names(32), iter_write_ch
write(uto,format3) input_names(33), iter_print
write(uto,format3) input_names(34), gradient_type
write(uto,format5) input_names(35), gradient_eta_ch
write(uto,format5) input_names(36), gradient_mu_ch
write(uto,format5) input_names(37), gradient_eps_ch
write(uto,format1) input_names(38)
write(uto,format1) input_names(39)
write(uto,format1) input_names(40)
write(uto,format3) input_names(41), enforce_NZ
write(uto,format3) input_names(42), opt_betalm
write(uto,format3) input_names(43), pairs_scheme  
write(uto,format5) input_names(44), constraint_eps_ch
do i = 3, constraint_types
  if ( constraint_switch(i) < 2 ) then 
    write(uto,format6) input_names(42+i), constraint_switch(i), &
                       constraint_read_ch(i,1)
  else
    write(uto,format6) input_names(38+i), constraint_switch(i), &
                       constraint_read_ch(i,1), constraint_read_ch(i,2)
  endif
enddo
print*,' '

deallocate(input_block,blocking_id_ch)
 
end subroutine print_input

!------------------------------------------------------------------------------!
! subroutine broadcast_inputs                                                  !
!                                                                              !
! Broadcasts the inputs read by process rank = 0 to all others. This is safer  !
! I think than making all the processes read the input file.                   !
!------------------------------------------------------------------------------!
!cmpi subroutine broadcast_inputs(hamil_dummy,iter_max,iter_write,iter_print)

!cmpi integer, intent(inout) :: iter_max, iter_write, iter_print
!cmpi character(100), intent(in) :: hamil_dummy
!cmpi integer :: ierr=0

!!! Hamiltonian
!cmpi call mpi_bcast(hamil_dummy,100,mpi_character,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(hamil_com,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(hamil_read,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(paral_teamssize,1,mpi_integer,0,mpi_comm_world,ierr)

!!! Particle Number
!cmpi call mpi_bcast(valence_Z,1,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(valence_N,1,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(proj_Mphip,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(proj_Mphin,1,mpi_integer,0,mpi_comm_world,ierr)

!!! Wave function
!cmpi call mpi_bcast(seed_type,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(blocking_dim,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(blocking_id,blocking_dim,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(seed_text,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(seed_symm,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(seed_occeps,1,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(seed_allemp,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(dens_spatial,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(dens_nr,3,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(dens_dr,3,mpi_double_precision,0,mpi_comm_world,ierr)

!!! Iterative procedure
!cmpi call mpi_bcast(iter_max,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(iter_write,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(iter_print,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(gradient_type,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(gradient_eta,1,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(gradient_mu,1,mpi_double_precision,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(gradient_eps,1,mpi_double_precision,0,mpi_comm_world,ierr)

!!! Constraint
!cmpi call mpi_bcast(enforce_NZ,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(opt_betalm,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(pairs_scheme,1,mpi_integer,0,mpi_comm_world,ierr)
!cmpi call mpi_bcast(constraint_eps,1,mpi_double_precision,0, &
!cmpi                 mpi_comm_world,ierr)
!cmpi call mpi_bcast (constraint_switch,constraint_max,mpi_integer,0, & 
!cmpi                 mpi_comm_world,ierr)
!cmpi call mpi_bcast (constraint_read,constraint_max,mpi_double_precision,0, & 
!cmpi                 mpi_comm_world,ierr)

!cmpi end subroutine broadcast_inputs

!------------------------------------------------------------------------------!
! subroutine check_input                                                       !
!                                                                              !
! Checks the values of the input parameters to see if they are appropriate.    !
! If not, the code will stop.                                                  !
!                                                                              !
! Input: iter_max = maximum number of iterations                               !
!        iter_write = step for the intermediate writing of wave functions      !
!        iter_print = option to calculate and print more at each iteration     !
!------------------------------------------------------------------------------!
subroutine check_input(iter_max,iter_write,iter_print)

integer, intent(in) :: iter_max, iter_write, iter_print
integer :: i, isum, ierror
integer, dimension(constraint_types) :: switch_check=0
real(r64) :: integerness
!cmpi integer :: ierr=0

!!! Counter for the number of errors (should be 0 at the end)
ierror = 0

!!!
!!! Hamiltonian    
!!!

!cmpi if ( paral_myrank == 0 ) then
if ( (hamil_com < 0) .or. (hamil_com > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for the COM correction (hamil_com) = ", & 
         hamil_com," should be 0 or 1."
endif 

if ( (hamil_read < 0) .or. (hamil_read > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to read an uncoupled file (hamil_read) = ", & 
         hamil_read," should be 0 or 1."
endif 

!cmpi if ( hamil_read /= 1 ) then
!cmpi   ierror = ierror + 1
!cmpi   print "(a,1i1,a)","The option to read an uncoupled file (hamil_read) &  
!cmpi   = ",hamil_read," should be 1 when doing MPI calculations."
!cmpi endif

!cmpi if ( (paral_teamssize < 0) .or. (paral_teamssize > paral_worldsize) ) then
!cmpi   ierror = ierror + 1
!cmpi   print "(a,1i1,a)","The numer of processes per team (MPI) = ", & 
!cmpi          paral_teamssize," should be positive and smaller than the", &
!cmpi         " total number of processes."
!cmpi endif 

!!!
!!! Particle number
!!!

if ( valence_Z < 0 ) then
  ierror = ierror + 1
  print "(a,1f10.3,a)","The number of active protons (valence_Z) = ", & 
         valence_Z," should be positive."
endif 

if ( valence_N < 0 ) then
  ierror = ierror + 1
  print "(a,1f10.3,a)","The number of active neutrons (valence_N) = ", & 
         valence_N," should be positive."
endif 

if ( (abs(valence_Z) + abs(valence_N)) <= epsilon0 ) then
  ierror = ierror + 1
  print "(a,1f10.3,a)","The number of active nucleons (valence_Z + valence_N) &
        &= ",valence_Z+valence_N," should be strictly positive."
endif 

if ( proj_Mphip < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of gauge angles for protons (proj_Mphip) = ", & 
         proj_Mphip," should be positive."
endif 

if ( proj_Mphin < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of gauge angles for neutrons (proj_Mphin) = ", & 
         proj_Mphin," should be positive."
endif 

integerness = valence_Z - int(valence_Z)  
if ( ((proj_Mphip > 1) .or. (seed_type > 6)) .and. & 
     (abs(integerness) > epsilon0) ) then
  ierror = ierror + 1
  print "(a,1f10.3,a)","The number of active protons (valence_Z) = ", & 
         valence_Z," is not an integer (needed for PNR or Slater)."
endif 

integerness = valence_N - int(valence_N)  
if ( ((proj_Mphip > 1) .or. (seed_type > 6)) .and. &
     (abs(integerness) > epsilon0) ) then
  ierror = ierror + 1
  print "(a,1f10.3,a)","The number of active neutrons (valence_N) = ", & 
         valence_N," is not an integer (needed for PNR or Slater)."
endif 

!!!
!!! Wave function
!!!

if ( (seed_type < 0) .or. (seed_type > 9) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for the type of seed (seed_type) = ", & 
         seed_type," should be between 0 and 9."
endif 

if ( blocking_dim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of blocked quasiparticles (blocking_dim) = " & 
         ,blocking_dim," should be positive."
endif

do i = 1, blocking_dim
  if ( blocking_id(i) < 0 ) then
    ierror = ierror + 1
    print "(a,1i5,a,1i5)","The index of the blocked quasiparticle &
          &(blocking_id(i)) = ", blocking_id(i), " for i = ", i
  endif
enddo

if ( (seed_symm < 0) .or. (seed_symm > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for checking the symmetries (seed_symm) = ", & 
         seed_symm," should be 0 or 1."
endif 

if ( seed_rand < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The seed for the random number generation (seed_rand) = ",& 
         seed_rand," should be positive or null."
endif 

if ( (seed_text < 0) .or. (seed_text > 3) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to write the wf in binary/text (seed_text) = " & 
         ,seed_text," should be between 0 or 3."
endif 

if ( seed_occeps < 0.0d0 ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The cutoff for the occupied single-particle states &
        &(seed_occeps) = ", seed_occeps," should be positive."
endif 

if ( (seed_allemp < 0) .or. (seed_allemp > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for including all the empty states &
        &(seed_allemp) = ", seed_allemp," should be 0 or 1."
endif 

if ( (dens_spatial < 0) .or. (dens_spatial > 3) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for including all the empty states &
        &(seed_allemp) = ", seed_allemp," should be 0 or 1."
endif 

do i = 1, 3
  if ( (dens_nr(i) < 0.d0) .or. (dens_dr(i) < 0.d0) ) then
    ierror = ierror + 1
    print "(a,1i1,a)","The values for the discretization of the spatial &
          &one-body density should be positive."
  endif
enddo

!!!
!!! Iterations
!!!

if ( iter_max < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The maximum number of iterations (iter_max) = ", & 
         iter_max," should be positive."
endif

if ( iter_write < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The step for intermediate wf writing (iter_write) = ", & 
         iter_write," should be positive."
endif

if ( (iter_print < 0) .or. (iter_print > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for extensive printing (iter_print) = ", & 
         iter_print," should be 0 or 1."
endif 

if ( (gradient_type < 0) .or. (gradient_type > 2) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option for the type of gradient (gradient_type) = ", & 
         gradient_type," should be beteween 0 and 2."
endif 

if ( gradient_eta < 0.0d0 ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The gradient parameter eta (gradient_eta) = ", & 
         gradient_eta," should be positive."
endif 

if ( (gradient_mu < 0.0d0) .or. (gradient_mu > 1.0d0) ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The gradient parameter mu (gradient_mu) = ", & 
         gradient_mu," should be between 0.0 and 1.0."
endif 

if ( gradient_eps <= 0.0d0 ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The tolerance for the gradient (gradient_eps) = ", & 
         gradient_eps," should be strictly positive."
endif 

!!!
!!! Constraints
!!!

where ((constraint_switch < 0) .and. (constraint_switch > 3)) switch_check = 1

do i = 18, constraint_types
 if ( constraint_switch(i) > 1 ) switch_check(i) = 1
enddo

isum = sum(switch_check)

if ( (enforce_NZ < 0) .or. (enforce_NZ > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to force the constraint on N/Z (enforce_NZ) = "& 
         ,enforce_NZ," should be 0 or 1."
endif 

if ( (opt_betalm < 0) .or. (opt_betalm > 2) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to constraint beta_lm (opt_betalm) = ", & 
         opt_betalm," should be between 0 and 2."
endif 

if ( (opt_betalm == 2) .and. ((constraint_read(5,1) < 0.0d0) .or. &
     (constraint_read(5,2) < 0.0d0)) ) then
  ierror = ierror + 1
  print "(a,2es11.3,a)","The values of beta (constraint_read(5,1:2)) = ", & 
         constraint_read(5,1),constraint_read(5,2)," should be positive."
endif 

if ( (constraint_switch(17) == 3) .and. (constraint_read(17,1) - &
     constraint_read(17,2) < 0.0d0) ) then
  ierror = ierror + 1
  print "(a,1es10.3,a,1es10.3,a)","The value of the isoscalar radius &
        &constraint (constraint_read(17,1))",constraint_read(17,1)," should &
        &be larger than the value of the isovector radius constraint &
        &(constraint_read(17,2))",constraint_read(17,2)
endif

if ( (pairs_scheme < 0) .or. (pairs_scheme > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The scheme for the coupling of pairs (pairs_scheme) = ", & 
         pairs_scheme," should be 0 or 1."
endif 

if ( constraint_switch(19) == 1 ) then
  ierror = ierror + 1
  print "(a)", "The constraint on <Jy> has to be switched off for now as the &
        &wave functions are real."
endif 

if ( constraint_eps <= 0.0d0 ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The tolerance for the constraints (constraint_eps) = "& 
         ,constraint_eps," should be strictly positive."
endif 

if ( isum /= 0 ) then
  ierror = ierror + 1
  print "(a,1i1,a)", "The flags to switch on/off the constraints & 
        &(constraint_switch) accept the values: 0, 1 (all), 2, 3 (radius and &
        &multipoles)."
endif
!cmpi endif 

!!!
!!! Stops the code if an error has been found in the input file
!!!

!cmpi call mpi_bcast(ierror,1,mpi_integer,0,mpi_comm_world,ierr)

if ( ierror /= 0 ) then
!cmpi   if ( paral_myrank == 0 ) then
  print "(a,1i2,a)", "The code has dectected ",ierror," problem(s) with the &
        &input parameters and will stop. Please check the manual."
!cmpi   endif 
  stop 
endif

end subroutine check_input

!------------------------------------------------------------------------------!
! subroutine open_files_hamiltonian                                            !
!                                                                              !
! Opens the hamiltonian files to link the units if they exist. If not, the     !
! code will stop.                                                              !
! BB: contrarily to subroutine check_input, when performing MPI calculations   !
! all the processes will print the error messages. I might want to change that !
! latter but it might be of advantage to keep it in case only part of the      !
! processes cannot access the files.                                           !
!------------------------------------------------------------------------------!
subroutine open_files_hamiltonian

integer :: htype, ierror, iwarn
character(100) :: hname, hname1, hname2, hnamer
logical :: is_exist            

!!! Counter for the number of errors (should be 0 at the end)
ierror = 0
iwarn  = 0

!!! Main file 
inquire (file=hamil_fsho, exist=is_exist)

if ( is_exist ) then
  open(uth, file=hamil_fsho, status='old', action='read', form='formatted')
  read(uth,'(a)') hname
  read(uth,*) htype
  rewind(uth)
else
  ierror = ierror + 1
  print "(a,a,a)", "The main file (hamil_fsho) = ", hamil_fsho, &
        " can not be found."
endif 

!!! Reduced file
if ( hamil_read == 1 ) then
  inquire (file=hamil_fred, exist=is_exist)
  if ( is_exist ) then
    open(uthr, file=hamil_fred, status='old', action='read', access='stream',&
         form='unformatted')              
    read(uthr) hnamer
    if ( hnamer /= hname ) then
      iwarn = iwarn + 1
      print "(a,a,a,a)", "Warning: the name in the reduced file = ", &
             trim(adjustl(hnamer)), &
            " does not correspond with the one of the main file = ", &
             trim(adjustl(hname))     
    endif
  else  
    ierror = ierror + 1
    print "(a,a,a)", "The uncoupled binary file (hamil_fred) = ",hamil_fred, &
          " can not be found."
  endif
endif

!!! 0+1 and 2-body files
select case ( htype*(1-hamil_read) )
  case (0:2)
    continue

  case (3:4)
    !!! 1-body
    if ( htype == 3 ) then
      inquire (file=hamil_f01b, exist=is_exist)
      if ( is_exist ) then
        open(uth1, file=hamil_f01b, status='old', action='read', & 
             form='formatted')
        read(uth1,'(a)') hname1
        rewind(uth1)

        if ( hname1 /= hname ) then
          iwarn = iwarn + 1
          print "(a,a,a,a)", "Warning: the name in the 0+1-body file = ", &
                 trim(adjustl(hname1)), &
                " does not correspond with the one of the main file = ",  &
                 trim(adjustl(hname))     
        endif
      else  
        ierror = ierror + 1
        print "(a,a,a)", "The 0+1-body file (hamil_f01b) = ",hamil_f01b, &
              " can not be found."
      endif
    endif

    !!! 2-body
    inquire (file=hamil_f2b, exist=is_exist)
    if ( is_exist ) then
      open(uth2, file=hamil_f2b, status='old', action='read', form='formatted')
      read(uth2,'(a)') hname2
      rewind(uth2)
      if ( hname2 /= hname ) then
        iwarn = iwarn + 1
        print "(a,a,a,a)", "Warning: the name in the 2-body file = ",  &  
               trim(adjustl(hname2)), &
              " does not correspond with the one of the main file = ", &
               trim(adjustl(hname))     
      endif
    else 
      ierror = ierror + 1
      print "(a,a,a)", "The 2-body file (hamil_f1) = ",hamil_f2b, &
            " can not be found."
    endif

    !!! Center of mass (2-body)
    if ( hamil_com == 1 ) then
      inquire (file=hamil_fcom, exist=is_exist)
      if ( is_exist ) then
        open(uthc, file=hamil_fcom, status='old', action='read', & 
             form='formatted')
      else  
        ierror = ierror + 1
        print "(a,a,a)", "The center-of-mass file (hamil_fcom) = ",hamil_fcom, &
              " can not be found."
      endif
    endif

  case default
    ierror = ierror + 1
    print "(a,1i1,a)", "The hamiltonian format (hamil_type) = ",htype, & 
          " should be between 1 and 4."

end select

!!!
!!! Stops the code if an error has been found with the hamiltonian files
!!!

if ( ierror /= 0 ) then
  print "(a,1i1,a)", "The code has dectected ",ierror," problem(s) with the &
        &hamiltonian files and will stop. Please check the files."
  stop 
endif

if ( iwarn /= 0 ) print*,' '

end subroutine open_files_hamiltonian

END MODULE Initialization
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
