!==============================================================================!
! MODULE Constants                                                             !
!                                                                              !
! This module contains the variables related to the physical, mathematical and !
! numerical constants.                                                         !
!==============================================================================!
MODULE Constants

use Iso_fortran_env

implicit none
public

!!! Definition of input/output units (for portability)
integer, parameter :: uti = input_unit,  &
                      uto = output_unit, &
                      utw = uti + uto + 1, &  ! wave functions
                      ute = uti + uto + 2, &  ! sp/qp energies, occ. numb.
                      uth = uti + uto + 3, &  ! hamiltonian (main)
                      uth1 = uth + 1, &       !     "       (1b)
                      uth2 = uth + 2, &       !     "       (2b)
                      uthc = uth + 4, &       !     "       (2b com)
                      uthr = uth + 5          !     "       (red)

!!! Definition of kind parameters (for portability)
integer, parameter :: i8  = int8,   & ! integer  8 bits
                      i16 = int16,  & !    "    16  "   
                      i32 = int32,  & !    "    32  "   
                      i64 = int64,  & !    "    64  "   
                      r32 = real32, & ! real 32 bits (single precision)
                      r64 = real64, & !  "   64  "   (double     "    )
                      rH2 = real32    ! for 2BME of the hamiltonian

!!! Definition of simple names for numerical values
real(r64), parameter :: one  = 1.0d0, &
                        zero = 0.0d0, &
                        epsilon0 = 1.0d-16
complex(r64), parameter :: zone  = (1.0d0,0.0d0), &
                           zzero = (0.0d0,0.0d0), &
                           zimag = (0.0d0,1.0d0)

!!! Definition of physical constants
real(r64), parameter :: pi = 4.0d0 * atan(1.0d0), & ! 3.14159265...
                        hbarc = 197.327053d0,     & ! hbar*c in Mev.fm
                        radius_r0 = 1.2d0,        & ! radius factor 
                        mass_mp = 938.27208816d0, & ! proton mass
                        mass_mn = 939.56542052d0, & ! neutron mass
                        mass_ma = (mass_mp + mass_mn)/2, & ! nucleon mass
                        hbarmass = hbarc**2 / (2*mass_ma)  ! factor kin. energy

END MODULE Constants
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
