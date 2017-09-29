!##########################################################################
!        FILE: incompact3dlmn.f90
!      AUTHOR: Paul Bartholomew <ptb08@imperial.ac.uk>
!     CREATED: 29#SEP#2017
! DESCRIPTION: This is the main file for incompact3dlmn.
!##########################################################################

PROGRAM incompact3dlmn

  USE MPI
  USE IBM

  IMPLICIT NONE

  INTEGER :: ierr
  INTEGER :: rank, size
  
  !--------------------------------------------
  ! Initialisation
  !--------------------------------------------

  CALL MPI_INIT(ierr)
  print *, "Initialising incompact3dlmn"

  !--------------------------------------------
  ! Time loop
  !--------------------------------------------
  print *, "Starting time loop"

  !--------------------------------------------
  ! Finalise
  !--------------------------------------------
  print *, "Finalising"

  CALL MPI_FINALIZE(ierr)

ENDPROGRAM incompact3dlmn

