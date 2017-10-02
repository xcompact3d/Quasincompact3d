!##########################################################################
!        FILE: incompact3dlmn.f90
!      AUTHOR: Paul Bartholomew <ptb08@imperial.ac.uk>
!     CREATED: 29-SEP-2017
! DESCRIPTION: This is the main file for incompact3dlmn.
!##########################################################################

PROGRAM incompact3dlmn

  USE MPI

  USE decomp_2d
  !USE IBM

  IMPLICIT NONE

  INTEGER :: ierr
  INTEGER :: rank, size

  INTEGER :: nx, ny, nz
  REAL :: lx, ly, lz

  nx = 0
  ny = 0
  nz = 0
  lx = 0.0
  ly = 0.0
  lz = 0.0

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

  CALL decomp_2d_init(nx, ny, nz, p_row, p_col)

  IF(rank == 0) THEN
    WRITE(*, 1001) lx, ly, lz
1001 FORMAT("Initialised domain with dimensions (", F9.3, ", ", F9.3, ", ", F9.3, ")")
    WRITE(*, 1002) nx, ny, nz
1002 FORMAT("Initialised mesh with              (", i9, ", ", i9, ", ", i9, ") cells")
  ENDIF
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

