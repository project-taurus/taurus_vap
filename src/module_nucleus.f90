!==============================================================================!
! MODULE Nucleus                                                               !
!                                                                              !
! This module contains the variables and routines related to the target nucl-  !
! eus.                                                                         !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_nucleus                                                     !
! - subroutine print_nucleus                                                   !
!==============================================================================!
MODULE Nucleus

!cmpi use MPI            
!cmpi use Parallelization
use Constants

implicit none
public

!!! Number of valence particles 
real(r64) :: valence_Z, & ! protons
             valence_N, & ! neutrons
             valence_A    ! nucleons
                        
!!! Number of core particles 
real(r64) :: core_Z, & ! protons
             core_N, & ! neutrons
             core_A    ! nucleons

!!! Total number of particles
real(r64) :: nucleus_Z, & ! protons
             nucleus_N, & ! neutrons
             nucleus_A    ! nucleons

!!! Private routines
private :: print_nucleus

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_nucleus                                                       !
!                                                                              !
! Computes the different number of particles. The values of valence_Z and      !
! valence_N are read in the input file.                                        !
!------------------------------------------------------------------------------!
subroutine set_nucleus       

integer :: i, ityp, itmp1

!!! Reads the number of core particle in the interaction file
rewind(uth)
read(uth,*)
read(uth,*) ityp
if ( (ityp == 1) .or. (ityp == 2) ) then
  do i = 1, ityp
    read(uth,*) 
  enddo
  read(uth,*) itmp1, core_Z, core_N
else 
  read(uth,*) 
  read(uth,*) core_Z, core_N
endif

!!! Computes all the remaining number of particles 
valence_A = valence_Z + valence_N
core_A    = core_Z    + core_N

nucleus_Z = core_Z + valence_Z
nucleus_N = core_N + valence_N
nucleus_A = core_A + valence_A

!!! Print the nucleus informations in the standard output
!cmpi if ( paral_myrank == 0 ) then        
call print_nucleus
!cmpi endif

end subroutine set_nucleus       

!------------------------------------------------------------------------------!
! subroutine print_nucleus                                                     !
!                                                                              !
! Prints the different number of particles.                                    !
!------------------------------------------------------------------------------!
subroutine print_nucleus       

character(len=*), parameter :: format1 = "(1a7,1x,1f10.2,1x,1f10.2,1x,1f10.2)"

print '(60("%"),/,27x,"NUCLEUS",26x,/,60("%"),//, &
      & "Part \ No.",5x,"Z",10x,"N",10x,"A",/,40("-"))'
print format1, 'Valence', valence_Z, valence_N, valence_A
print format1, 'Core   ',    core_Z,    core_N,    core_A
print format1, 'Nucleus', nucleus_Z, nucleus_N, nucleus_A

end subroutine print_nucleus       

END MODULE Nucleus
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
