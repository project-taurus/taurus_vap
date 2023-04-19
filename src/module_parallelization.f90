!==============================================================================!
! MODULE Parallelization                                                       !
!                                                                              !
! This module contains the variables and routines related to the initialization!
! of MPI parallel runs.                                                        !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_parallel_teams                                              !
!==============================================================================!
MODULE Parallelization   

use MPI
use MathMethods

implicit none
public

!!! Informations on the divide of tasks (ranks, team, job, ...)
integer :: paral_worldsize,  & ! Total number of processes
           paral_myrank,     & ! Rank of a given process in world
           mpi_comm_team,    & ! Internal team communicator
           mpi_comm_peers,   & ! Internal team communicator
           paral_teams,      & ! Number of teams in world
           paral_teamssize,  & ! Number of desired processes per team 
           paral_myteam,     & ! Rank of a given team in world
           paral_myteamrank, & ! Rank of a given process in its team
           paral_myteamsize, & ! Number of processes in a given team 
           paral_myangles,   & ! Number of gauge angles for a given team
           paral_myoffset      ! Offset for gauge angles for a given team

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_parallel_teams                                                !
!                                                                              !
! Sets the characteristics of the parallel teams (number, communicator, assign-!
! ment, size, rank of processes) which compute different gauge angles.         !
!------------------------------------------------------------------------------!
subroutine set_parallel_teams            

integer :: divide, rest, ierr=0

!!! Determines the number of teams and creates a new communicator
divide = paral_worldsize / paral_teamssize
rest = modulo(paral_worldsize,paral_teamssize)

paral_teams = divide + 1 - kdelta(rest,0)
paral_myteam = paral_myrank / paral_teamssize

!!! Creates the different teams
call mpi_comm_split(mpi_comm_world,paral_myteam,paral_myrank,mpi_comm_team,ierr)
call mpi_comm_size(mpi_comm_team,paral_myteamsize,ierr)
call mpi_comm_rank(mpi_comm_team,paral_myteamrank,ierr)

!!! Creates the different peers (same teamrank across teams)       
call mpi_comm_split(mpi_comm_world,paral_myteamrank,paral_myteam, & 
                    mpi_comm_peers,ierr)

end subroutine set_parallel_teams                

END MODULE Parallelization 
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
