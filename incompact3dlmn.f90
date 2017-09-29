!##########################################################################
!        FILE: incompact3dlmn.f90
!      AUTHOR: Paul Bartholomew <ptb08@imperial.ac.uk>
!     CREATED: 29#SEP#2017
! DESCRIPTION: This is the main file for incompact3dlmn.
!##########################################################################

PROGRAM incompact3dlmn

  USE MPI
  !USE IBM

  IMPLICIT NONE

  INTEGER :: ierr
  INTEGER :: rank, size

  ! INTEGER :: itime, ifirst, ilast
  ! REAL :: flowtime, dt
  
  !--------------------------------------------
  ! Initialisation
  !--------------------------------------------

  ! MPI init
  CALL MPI_INIT(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  PRINT *, "Initialising incompact3dlmn"

  !--------------------------------------------
  ! Time loop
  !--------------------------------------------

  PRINT *, "Starting time loop"
!   DO itime = ifirst, ilast

!      flowtime = (itime - 1) * dt
!      WRITE(*, 1001) itime, flowtime
! 1001 FORMAT("Timestep = ", i7, ", Time unit = ", F9.3)

!      ! Convection-Diffusion
!      ! Scalar Convection-Diffusion?
!      ! Velocity divergence
!      ! Pressure-Poisson
!      ! Pressure gradient
!      ! Velocity correction

!   ENDDO

  !--------------------------------------------
  ! Finalise
  !--------------------------------------------
  PRINT *, "Finalising"

  CALL MPI_FINALIZE(ierr)

ENDPROGRAM incompact3dlmn

