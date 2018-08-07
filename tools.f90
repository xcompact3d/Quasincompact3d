!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

!*******************************************************************
!
!
!*******************************************************************
subroutine test_scalar_min_max(phi)

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI


  implicit none

  integer :: i,j,k
  real(mytype) :: phimax,phimin,cfl
  real(mytype) :: phimax1,phimin1

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi
  integer :: code

  phimax=-1609._mytype
  phimin=1609._mytype
  ! do k=1,xsize(3)
  !   do j=1,xsize(2)
  !     do i=1,xsize(1)
  !       if (phi(i,j,k).gt.phimax) phimax=phi(i,j,k)
  !       if (phi(i,j,k).lt.phimin) phimin=phi(i,j,k)
  !     enddo
  !   enddo
  ! enddo
  phimax = MAXVAL(phi)
  phimin = MINVAL(phi)

  call MPI_REDUCE(phimax,phimax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(phimin,phimin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
  if (nrank==0) then
    print *,'PHI max=',phimax1
    print *,'PHI min=',phimin1
  endif

  return
end subroutine test_scalar_min_max

!*******************************************************************
!
!
!*******************************************************************
subroutine test_density_min_max(rho)

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  integer :: i,j,k
  real(mytype) :: rhomax,rhomin,cfl
  real(mytype) :: rhomax1,rhomin1

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: rho
  integer :: code

  ! rhomax=-1609._mytype
  ! rhomin=1609._mytype
  ! do k=1,xsize(3)
  !   do j=1,xsize(2)
  !     do i=1,xsize(1)
  !       if (rho(i,j,k).gt.rhomax) rhomax=rho(i,j,k)
  !       if (rho(i,j,k).lt.rhomin) rhomin=rho(i,j,k)
  !     enddo
  !   enddo
  ! enddo
  rhomax = MAXVAL(rho)
  rhomin = MINVAL(rho)

  call MPI_REDUCE(rhomax,rhomax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(rhomin,rhomin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
  if (nrank==0) then
    print *,'RHO max=',rhomax1
    print *,'RHO min=',rhomin1
  endif

  return
end subroutine test_density_min_max

subroutine test_temperature_min_max(temperature)

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  integer :: i,j,k
  real(mytype) :: temperaturemax,temperaturemin,cfl
  real(mytype) :: temperaturemax1,temperaturemin1

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temperature
  integer :: code

  ! temperaturemax=-1609._mytype
  ! temperaturemin=1609._mytype
  ! do k=1,xsize(3)
  !   do j=1,xsize(2)
  !     do i=1,xsize(1)
  !       if (temperature(i,j,k).gt.temperaturemax) temperaturemax=temperature(i,j,k)
  !       if (temperature(i,j,k).lt.temperaturemin) temperaturemin=temperature(i,j,k)
  !     enddo
  !   enddo
  ! enddo
  temperaturemax = MAXVAL(temperature)
  temperaturemin = MINVAL(temperature)

  call MPI_REDUCE(temperaturemax,temperaturemax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(temperaturemin,temperaturemin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
  if (nrank==0) then
    print *,'TEMPERATURE max=',temperaturemax1
    print *,'TEMPERATURE min=',temperaturemin1
  endif

  return
end subroutine test_temperature_min_max

!*******************************************************************
!
!
!*******************************************************************
subroutine test_speed_min_max(ux,uy,uz)

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI


  implicit none

  integer :: i,j,k
  real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin,cfl
  real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
  integer :: code

  ! uxmax=-1609._mytype
  ! uymax=-1609._mytype
  ! uzmax=-1609._mytype
  ! uxmin=1609._mytype
  ! uymin=1609._mytype
  ! uzmin=1609._mytype
  ! do k=1,xsize(3)
  !   do j=1,xsize(2)
  !     do i=1,xsize(1)
  !       if (ux(i,j,k).gt.uxmax) uxmax=ux(i,j,k)
  !       if (uy(i,j,k).gt.uymax) uymax=uy(i,j,k)
  !       if (uz(i,j,k).gt.uzmax) uzmax=uz(i,j,k)
  !       if (ux(i,j,k).lt.uxmin) uxmin=ux(i,j,k)
  !       if (uy(i,j,k).lt.uymin) uymin=uy(i,j,k)
  !       if (uz(i,j,k).lt.uzmin) uzmin=uz(i,j,k)
  !     enddo
  !   enddo
  ! enddo

  uxmax = MAXVAL(ux)
  uxmin = MINVAL(ux)
  uymax = MAXVAL(uy)
  uymin = MINVAL(uy)
  uzmax = MAXVAL(uz)
  uzmin = MINVAL(uz)

  call MPI_REDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(uymax,uymax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(uzmax,uzmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(uymin,uymin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(uzmin,uzmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
  if (nrank==0) then
    print *,'U,V,W max=',uxmax1,uymax1,uzmax1
    print *,'U,V,W min=',uxmin1,uymin1,uzmin1
  endif

  return
end subroutine test_speed_min_max

!*******************************************************************
!
!
!*******************************************************************
subroutine restart(ux1,uy1,uz1,rho1,temperature1,ep1,pp3,phi1,&
     gx1,gy1,gz1,rhos1,px1,py1,pz1,phis1,&
     hx1,hy1,hz1,rhoss1,phiss1,&
     pressure0,phG,irestart)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  TYPE(DECOMP_INFO) :: phG
  integer :: i,j,k,ijk,irestart,nzmsize,fh,ierror,code
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: rho1,temperature1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: rhos1,rhoss1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: phi1, phis1, phiss1
  real(mytype), dimension(phG%zst(1):phG%zen(1),phG%zst(2):phG%zen(2),phG%zst(3):phG%zen(3)) :: pp3
  real(mytype) :: pressure0
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  real(mytype) :: xdt
  integer, dimension(2) :: dims, dummy_coords
  logical, dimension(2) :: dummy_periods

  if (iscalar==0) then
    if (nscheme.ne.4) then
      if (irestart==1) then !write
        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
             MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
             fh, ierror)
        filesize = 0_MPI_OFFSET_KIND
        call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_write_var(fh,disp,1,ux1)
        call decomp_2d_write_var(fh,disp,1,uy1)
        call decomp_2d_write_var(fh,disp,1,uz1)
        call decomp_2d_write_var(fh,disp,1,rho1)
        call decomp_2d_write_var(fh,disp,1,temperature1)
        call decomp_2d_write_var(fh,disp,1,gx1)
        call decomp_2d_write_var(fh,disp,1,gy1)
        call decomp_2d_write_var(fh,disp,1,gz1)
        call decomp_2d_write_var(fh,disp,1,rhos1)
        call decomp_2d_write_var(fh,disp,1,px1)
        call decomp_2d_write_var(fh,disp,1,py1)
        call decomp_2d_write_var(fh,disp,1,pz1)
        call decomp_2d_write_var(fh,disp,1,pp3,phG)
        call MPI_FILE_CLOSE(fh,ierror)
      else
        if (nrank==0) then
          print *,'RESTART'
        endif
        
        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
             MPI_MODE_RDONLY, MPI_INFO_NULL, &
             fh, ierror)
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_read_var(fh,disp,1,ux1)   
        call decomp_2d_read_var(fh,disp,1,uy1)
        call decomp_2d_read_var(fh,disp,1,uz1)
        call decomp_2d_read_var(fh,disp,1,rho1)
        call decomp_2d_read_var(fh,disp,1,temperature1)
        call decomp_2d_read_var(fh,disp,1,gx1)
        call decomp_2d_read_var(fh,disp,1,gy1)
        call decomp_2d_read_var(fh,disp,1,gz1)
        call decomp_2d_read_var(fh,disp,1,rhos1)
        call decomp_2d_read_var(fh,disp,1,px1)
        call decomp_2d_read_var(fh,disp,1,py1)
        call decomp_2d_read_var(fh,disp,1,pz1)
        call decomp_2d_read_var(fh,disp,1,pp3,phG)
        call MPI_FILE_CLOSE(fh,ierror)
      endif
    else !AB3
      if (irestart==1) then !write
        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
             MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
             fh, ierror)
        filesize = 0_MPI_OFFSET_KIND
        call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_write_var(fh,disp,1,ux1)
        call decomp_2d_write_var(fh,disp,1,uy1)
        call decomp_2d_write_var(fh,disp,1,uz1)
        call decomp_2d_write_var(fh,disp,1,rho1)
        call decomp_2d_write_var(fh,disp,1,temperature1)
        call decomp_2d_write_var(fh,disp,1,gx1)
        call decomp_2d_write_var(fh,disp,1,gy1)
        call decomp_2d_write_var(fh,disp,1,gz1)
        call decomp_2d_write_var(fh,disp,1,rhos1)
        call decomp_2d_write_var(fh,disp,1,hx1)
        call decomp_2d_write_var(fh,disp,1,hy1)
        call decomp_2d_write_var(fh,disp,1,hz1)
        call decomp_2d_write_var(fh,disp,1,rhoss1)
        call decomp_2d_write_var(fh,disp,1,px1)
        call decomp_2d_write_var(fh,disp,1,py1)
        call decomp_2d_write_var(fh,disp,1,pz1)
        call decomp_2d_write_var(fh,disp,1,pp3,phG)
        call MPI_FILE_CLOSE(fh,ierror)
      else
        if (nrank==0) print *,'RESTART'
        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
             MPI_MODE_RDONLY, MPI_INFO_NULL, &
             fh, ierror)
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_read_var(fh,disp,1,ux1)   
        call decomp_2d_read_var(fh,disp,1,uy1)
        call decomp_2d_read_var(fh,disp,1,uz1)
        call decomp_2d_read_var(fh,disp,1,rho1)
        call decomp_2d_read_var(fh,disp,1,temperature1)
        call decomp_2d_read_var(fh,disp,1,gx1)
        call decomp_2d_read_var(fh,disp,1,gy1)
        call decomp_2d_read_var(fh,disp,1,gz1)
        call decomp_2d_read_var(fh,disp,1,rhos1)
        call decomp_2d_read_var(fh,disp,1,hx1)
        call decomp_2d_read_var(fh,disp,1,hy1)
        call decomp_2d_read_var(fh,disp,1,hz1)
        call decomp_2d_read_var(fh,disp,1,rhoss1)
        call decomp_2d_read_var(fh,disp,1,px1)
        call decomp_2d_read_var(fh,disp,1,py1)
        call decomp_2d_read_var(fh,disp,1,pz1)
        call decomp_2d_read_var(fh,disp,1,pp3,phG)
        call MPI_FILE_CLOSE(fh,ierror)
      endif
    endif
  else !SCALAR 
    if (nscheme.ne.4) then
      if (irestart==1) then !write
        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
             MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
             fh, ierror)
        filesize = 0_MPI_OFFSET_KIND
        call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_write_var(fh,disp,1,ux1)
        call decomp_2d_write_var(fh,disp,1,uy1)
        call decomp_2d_write_var(fh,disp,1,uz1)
        call decomp_2d_write_var(fh,disp,1,rho1)
        call decomp_2d_write_var(fh,disp,1,temperature1)
        call decomp_2d_write_var(fh,disp,1,gx1)
        call decomp_2d_write_var(fh,disp,1,gy1)
        call decomp_2d_write_var(fh,disp,1,gz1)
        call decomp_2d_write_var(fh,disp,1,rhos1)
        call decomp_2d_write_var(fh,disp,1,px1)
        call decomp_2d_write_var(fh,disp,1,py1)
        call decomp_2d_write_var(fh,disp,1,pz1)
        call decomp_2d_write_var(fh,disp,1,pp3,phG)
        call decomp_2d_write_var(fh,disp,1,phi1)
        call decomp_2d_write_var(fh,disp,1,phis1)
        call MPI_FILE_CLOSE(fh,ierror)
      else
        if (nrank==0) then
          print *,'RESTART'
        endif
        
        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
             MPI_MODE_RDONLY, MPI_INFO_NULL, &
             fh, ierror)
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_read_var(fh,disp,1,ux1)   
        call decomp_2d_read_var(fh,disp,1,uy1)
        call decomp_2d_read_var(fh,disp,1,uz1)
        call decomp_2d_read_var(fh,disp,1,rho1)
        call decomp_2d_read_var(fh,disp,1,temperature1)
        call decomp_2d_read_var(fh,disp,1,gx1)
        call decomp_2d_read_var(fh,disp,1,gy1)
        call decomp_2d_read_var(fh,disp,1,gz1)
        call decomp_2d_read_var(fh,disp,1,rhos1)
        call decomp_2d_read_var(fh,disp,1,px1)
        call decomp_2d_read_var(fh,disp,1,py1)
        call decomp_2d_read_var(fh,disp,1,pz1)
        call decomp_2d_read_var(fh,disp,1,pp3,phG)
        call decomp_2d_read_var(fh,disp,1,phi1)
        call decomp_2d_read_var(fh,disp,1,phis1)
        call MPI_FILE_CLOSE(fh,ierror)
      endif
    else !SCALAR + AB3
      if (irestart==1) then !write
        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
             MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
             fh, ierror)
        filesize = 0_MPI_OFFSET_KIND
        call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_write_var(fh,disp,1,ux1)
        call decomp_2d_write_var(fh,disp,1,uy1)
        call decomp_2d_write_var(fh,disp,1,uz1)
        call decomp_2d_write_var(fh,disp,1,rho1)
        call decomp_2d_write_var(fh,disp,1,temperature1)
        call decomp_2d_write_var(fh,disp,1,gx1)
        call decomp_2d_write_var(fh,disp,1,gy1)
        call decomp_2d_write_var(fh,disp,1,gz1)
        call decomp_2d_write_var(fh,disp,1,rhos1)
        call decomp_2d_write_var(fh,disp,1,hx1)
        call decomp_2d_write_var(fh,disp,1,hy1)
        call decomp_2d_write_var(fh,disp,1,hz1)
        call decomp_2d_write_var(fh,disp,1,rhoss1)
        call decomp_2d_write_var(fh,disp,1,px1)
        call decomp_2d_write_var(fh,disp,1,py1)
        call decomp_2d_write_var(fh,disp,1,pz1)
        call decomp_2d_write_var(fh,disp,1,pp3,phG)
        call decomp_2d_write_var(fh,disp,1,phi1)
        call decomp_2d_write_var(fh,disp,1,phis1)
        call decomp_2d_write_var(fh,disp,1,phiss1)
        call MPI_FILE_CLOSE(fh,ierror)
      else
        if (nrank==0) then
          print *,'RESTART'
        endif

        
        call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
             MPI_MODE_RDONLY, MPI_INFO_NULL, &
             fh, ierror)
        disp = 0_MPI_OFFSET_KIND
        call decomp_2d_read_var(fh,disp,1,ux1)   
        call decomp_2d_read_var(fh,disp,1,uy1)
        call decomp_2d_read_var(fh,disp,1,uz1)
        call decomp_2d_read_var(fh,disp,1,rho1)
        call decomp_2d_read_var(fh,disp,1,temperature1)
        call decomp_2d_read_var(fh,disp,1,gx1)
        call decomp_2d_read_var(fh,disp,1,gy1)
        call decomp_2d_read_var(fh,disp,1,gz1)
        call decomp_2d_read_var(fh,disp,1,rhos1)
        call decomp_2d_read_var(fh,disp,1,hx1)
        call decomp_2d_read_var(fh,disp,1,hy1)
        call decomp_2d_read_var(fh,disp,1,hz1)
        call decomp_2d_read_var(fh,disp,1,hz1)
        call decomp_2d_read_var(fh,disp,1,rhoss1)
        call decomp_2d_read_var(fh,disp,1,px1)
        call decomp_2d_read_var(fh,disp,1,py1)
        call decomp_2d_read_var(fh,disp,1,pz1)
        call decomp_2d_read_var(fh,disp,1,pp3,phG)
        call decomp_2d_read_var(fh,disp,1,phi1)
        call decomp_2d_read_var(fh,disp,1,phis1)
        call decomp_2d_read_var(fh,disp,1,phiss1)
        call MPI_FILE_CLOSE(fh,ierror)
      endif
    endif
  endif

  if (irestart==0) then ! We are READING data
    ! reconstruction of the dp/dx, dp/dy and dp/dz from px1,py1 and pz1
    ! Temporal scheme (1:AB2, 2: RK3, 3:RK4, 4:AB3)
    if (nscheme==1) then
      xdt=adt(1)+bdt(1)
    else if (nscheme==2) then
      xdt=(3._mytype/4._mytype)*dt + (-5._mytype/12._mytype)*dt
    else if (nscheme==3) then
      xdt=0.041717869325_mytype*dt
    else if (nscheme==4) then
      xdt=adt(1)+bdt(1)+cdt(1)
    else
      if (nrank.eq.0) then
        PRINT *, "Unknown nscheme=", nscheme
      endif
    endif

    do k=1,xsize(3)
      do j=1,xsize(2)
        dpdyx1(j,k)=py1(1,j,k)/xdt
        dpdzx1(j,k)=pz1(1,j,k)/xdt
        dpdyxn(j,k)=py1(nx,j,k)/xdt
        dpdzxn(j,k)=pz1(nx,j,k)/xdt
      enddo
    enddo

    if (xsize(3)==1) then
      do j=1,xsize(2)
        do i=1,xsize(1)
          dpdxz1(i,j)=px1(i,j,1)/xdt
          dpdyz1(i,j)=py1(i,j,1)/xdt
        enddo
      enddo
    endif
    if (xsize(3)==nz) then
      do j=1,xsize(2)
        do i=1,xsize(1)
          dpdxzn(i,j)=px1(i,j,nz)/xdt
          dpdyzn(i,j)=py1(i,j,nz)/xdt
        enddo
      enddo
    endif

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, code)

    if (dims(1)==1) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxy1(i,k)=px1(i,1,k)/xdt
          dpdzy1(i,k)=pz1(i,1,k)/xdt
        enddo
      enddo
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxyn(i,k)=px1(i,xsize(2),k)/xdt
          dpdzyn(i,k)=pz1(i,xsize(2),k)/xdt
        enddo
      enddo
    else
      !find j=1 and j=ny
      if (xstart(2)==1) then
        do k=1,xsize(3)
          do i=1,xsize(1)
            dpdxy1(i,k)=px1(i,1,k)/xdt
            dpdzy1(i,k)=pz1(i,1,k)/xdt
          enddo
        enddo
      endif
      !      print *,nrank,xstart(2),ny-(nym/p_row)
      if (ny-(nym/dims(1))==xstart(2)) then
        do k=1,xsize(3)
          do i=1,xsize(1)
            dpdxyn(i,k)=px1(i,xsize(2),k)/xdt
            dpdzyn(i,k)=pz1(i,xsize(2),k)/xdt
          enddo
        enddo
      endif

    endif


    if (nrank==0) print *,'reconstruction pressure gradients done!'

    !! Reconstruct thermodynamic pressure from density and
    !  temperature fields.
    pressure0 = 0._mytype
    do ijk = 1, xsize(1) * xsize(2) * xsize(3)
      pressure0 = pressure0 + rho1(ijk, 1, 1) * temperature1(ijk, 1, 1)
    enddo
    pressure0 = pressure0 / float(xsize(1) * xsize(2) * xsize(3))
    call MPI_Allreduce(MPI_IN_PLACE, pressure0, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
    pressure0 = pressure0 / float(nproc)
  endif

  if (irestart==0) then
    if (ivirt==2) then
      call MPI_FILE_OPEN(MPI_COMM_WORLD, 'epsilon.dat', &
           MPI_MODE_RDONLY, MPI_INFO_NULL, &
           fh, ierror)
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,1,ep1) 
      call MPI_FILE_CLOSE(fh,ierror)
      if (nrank==0) print *,'read epsilon file done from restart'
    endif
  endif

end subroutine restart

!*******************************************************************
!
!
!*******************************************************************
subroutine stretching()
  
  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  real(mytype) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
  integer :: j

  yinf=-yly/2._mytype
  den=2._mytype*beta*yinf
  xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
  alpha=abs(xnum/den)
  xcx=1._mytype/beta/alpha
  if (alpha.ne.0._mytype) then
    if (istret.eq.1) yp(1)=0._mytype
    if (istret.eq.2) yp(1)=0._mytype
    if (istret.eq.1) yeta(1)=0._mytype
    if (istret.eq.2) yeta(1)=-1._mytype/2._mytype
    if (istret.eq.3) yp(1)=0._mytype
    if (istret.eq.3) yeta(1)=-1._mytype/2._mytype
    do j=2,ny!-1
      if (istret==1) then
        if (ncly.eq.0) yeta(j)=(j-1._mytype)*(1._mytype/ny)
        if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1._mytype)*(1._mytype/(ny-1._mytype))
      endif
      if (istret==2) then
        if (ncly.eq.0) yeta(j)=(j-1._mytype)*(1._mytype/ny)-0.5_mytype
        if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1._mytype)*(1._mytype/(ny-1._mytype))-0.5_mytype
      endif
      if (istret==3) then
        if (ncly.eq.0) yeta(j)=((j-1._mytype)*(1._mytype/2._mytype/ny)-0.5_mytype)
        if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=((j-1._mytype)*(1._mytype/2._mytype/(ny-1._mytype))-0.5_mytype)
      endif
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2._mytype*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
      den4=2._mytype*alpha*beta-cos(2._mytype*pi*yeta(j))+1._mytype
      if ((ncly.ne.0).and.(j==ny).and.(istret==1)) then
        xnum1=0._mytype
      else
        xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
      endif
      cst=sqrt(beta)*pi/(2._mytype*sqrt(alpha)*sqrt(alpha*beta+1._mytype))
      if (istret==1) then
        if (yeta(j).lt.0.5_mytype) yp(j)=xnum1-cst-yinf
        if (yeta(j).eq.0.5_mytype) yp(j)=0._mytype-yinf
        if (yeta(j).gt.0.5_mytype) yp(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
        if (yeta(j).lt.0.5_mytype) yp(j)=xnum1-cst+yly
        if (yeta(j).eq.0.5_mytype) yp(j)=0._mytype+yly
        if (yeta(j).gt.0.5_mytype) yp(j)=xnum1+cst+yly
      endif
      if (istret==3) then
        if (yeta(j).lt.0.5_mytype) yp(j)=(xnum1-cst+yly)*2._mytype
        if (yeta(j).eq.0.5_mytype) yp(j)=(0._mytype+yly)*2._mytype
        if (yeta(j).gt.0.5_mytype) yp(j)=(xnum1+cst+yly)*2._mytype
      endif
    enddo

  endif
  !if (nrank==0) then
  !do j=1,ny
  !print *,j,yp(j),yeta(j)
  !enddo
  !endif
  !stop
  if (alpha.eq.0.) then
    yp(1)=-1.e10_mytype
    do j=2,ny
      yeta(j)=(j-1._mytype)*(1._mytype/ny)
      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
    enddo
  endif
  if (alpha.ne.0._mytype) then
    do j=1,ny
      if (istret==1) then
        if (ncly.eq.0) then
          yetai(j)=(j-0.5_mytype)*(1._mytype/ny)
        else if ((ncly.eq.1).or.(ncly.eq.2)) then
          yetai(j)=(j-0.5_mytype)*(1._mytype/(ny-1._mytype))
        endif
      else if (istret==2) then
        if (ncly.eq.0) then
          yetai(j)=(j-0.5_mytype)*(1._mytype/ny)-0.5_mytype
        else if ((ncly.eq.1).or.(ncly.eq.2)) then
          yetai(j)=(j-0.5_mytype)*(1._mytype/(ny-1._mytype))-0.5_mytype
        endif
      else if (istret==3) then
        if (ncly.eq.0) then
          yetai(j)=(j-0.5)*(1./2./ny)-0.5_mytype
        else if ((ncly.eq.1).or.(ncly.eq.2)) then
          yetai(j)=(j-0.5_mytype)*(1._mytype/2._mytype/(ny-1._mytype))-0.5_mytype
        endif
      endif

      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2._mytype*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
      den4=2._mytype*alpha*beta-cos(2._mytype*pi*yetai(j))+1.
      xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2._mytype*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
        if (yetai(j).lt.0.5_mytype) then
          ypi(j)=xnum1-cst-yinf
        else if (yetai(j).gt.0.5_mytype) then
          ypi(j)=xnum1+cst-yinf
        else
          ypi(j)=0._mytype-yinf
        endif
      else if (istret==2) then
        if (yetai(j).lt.0.5_mytype) then
          ypi(j)=xnum1-cst+yly
        else if (yetai(j).gt.0.5_mytype) then
          ypi(j)=xnum1+cst+yly
        else
          ypi(j)=0._mytype+yly
        endif
      endif
      if (istret==3) then
        if (yetai(j).lt.0.5_mytype) then
          ypi(j)=(xnum1-cst+yly)*2._mytype
        else if (yetai(j).gt.0.5_mytype) then
          ypi(j)=(xnum1+cst+yly)*2._mytype
        else
          ypi(j)=(0._mytype+yly)*2._mytype
        endif
      endif
    enddo
  endif
  if (alpha.eq.0._mytype) then
    ypi(1)=-1.e10_mytype
    do j=2,ny
      yetai(j)=(j-1._mytype)*(1._mytype/ny)
      ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
    enddo
  endif

  do j=1,ny
    ppy(j)=yly*(alpha/pi+(1./pi/beta)*sin(pi*yeta(j))* &
         sin(pi*yeta(j)))
    pp2y(j)=ppy(j)*ppy(j)
    pp4y(j)=(-2._mytype/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
  enddo
  do j=1,ny
    ppyi(j)=yly*(alpha/pi+(1./pi/beta)*sin(pi*yetai(j))* &
         sin(pi*yetai(j)))
    pp2yi(j)=ppyi(j)*ppyi(j)
    pp4yi(j)=(-2._mytype/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
  enddo

end subroutine stretching

!*****************************************************************
!
!
!*****************************************************************
subroutine inversion5_v1(aaa,eee,spI)

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16_mytype
#else
  real(mytype), parameter :: epsilon = 1.e-8_mytype
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  do i=1,2
    ja(i)=4-i
    jb(i)=5-i
  enddo
  do m=1,ny/2-2
    do i=1,2
      mi=m+i
      do k=spI%yst(3),spI%yen(3)
        do j=spI%yst(1),spI%yen(1) 
          if (real(aaa(j,m,k,3), kind=mytype).ne.0._mytype) then
            tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
          endif
          if (aimag(aaa(j,m,k,3)).ne.0._mytype) then
            tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
          endif
          sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
          eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
               aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
        enddo
      enddo
      do jc=ja(i),jb(i)
        do k=spI%yst(3),spI%yen(3)
          do j=spI%yst(1),spI%yen(1) 
            aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)&
                 * real(aaa(j,m,k,jc+i), kind=mytype), aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))&
                 * aimag(aaa(j,m,k,jc+i)), kind=mytype)
          enddo
        enddo
      enddo
    enddo
  enddo


  do k=spI%yst(3),spI%yen(3)
    do j=spI%yst(1),spI%yen(1) 
      if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
        tmp1=real(aaa(j,ny/2,k,2), kind=mytype)/real(aaa(j,ny/2-1,k,3), kind=mytype)
      else
        tmp1=0._mytype
      endif
      if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
        tmp2=aimag(aaa(j,ny/2,k,2))/aimag(aaa(j,ny/2-1,k,3))
      else
        tmp2=0._mytype
      endif
      sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
      b1(j,k)=cmplx(real(aaa(j,ny/2,k,3), kind=mytype)-tmp1*real(aaa(j,ny/2-1,k,4), kind=mytype),&
           aimag(aaa(j,ny/2,k,3))-tmp2*aimag(aaa(j,ny/2-1,k,4)), kind=mytype)

      if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
        tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
        tmp3=real(eee(j,ny/2,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,ny/2-1,k), kind=mytype)
      else
        tmp1=0._mytype
        tmp3=0._mytype
      endif
      if (abs(aimag(b1(j,k))).gt.epsilon) then
        tmp2=aimag(sr(j,k))/aimag(b1(j,k))
        tmp4=aimag(eee(j,ny/2,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,ny/2-1,k))
      else
        tmp2=0._mytype
        tmp4=0._mytype
      endif
      a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
      eee(j,ny/2,k)=cmplx(tmp3,tmp4, kind=mytype)

      if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
        tmp1=1._mytype/real(aaa(j,ny/2-1,k,3), kind=mytype)
      else
        tmp1=0.
      endif
      if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
        tmp2=1._mytype/aimag(aaa(j,ny/2-1,k,3))
      else
        tmp2=0._mytype
      endif
      b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
      a1(j,k)=cmplx(real(aaa(j,ny/2-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
           aimag(aaa(j,ny/2-1,k,4))*aimag(b1(j,k)), kind=mytype)
      eee(j,ny/2-1,k)=cmplx(real(eee(j,ny/2-1,k))*real(b1(j,k))-real(a1(j,k))*real(eee(j,ny/2,k)),&
           aimag(eee(j,ny/2-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,ny/2,k)), kind=mytype)
    enddo
  enddo

  do i=ny/2-2,1,-1
    do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1) 
        if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
          tmp1=1._mytype/real(aaa(j,i,k,3), kind=mytype)
        else
          tmp1=0._mytype
        endif
        if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
          tmp2=1._mytype/aimag(aaa(j,i,k,3))
        else
          tmp2=0._mytype
        endif
        sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
        a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
             aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
        b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
             aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
        eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
             real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
             real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
             aimag(eee(j,i,k))*aimag(sr(j,k))-&
             aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
      enddo
    enddo
  enddo

  return

end subroutine inversion5_v1

!*****************************************************************
!
!
!*****************************************************************
subroutine inversion5_v2(aaa,eee,spI)

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16_mytype
#else
  real(mytype), parameter :: epsilon = 1.e-8_mytype
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3),5) :: aaa
  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  do i=1,2
    ja(i)=4-i
    jb(i)=5-i
  enddo
  do m=1,nym-2
    do i=1,2
      mi=m+i
      do k=spI%yst(3),spI%yen(3)
        do j=spI%yst(1),spI%yen(1) 
          if (real(aaa(j,m,k,3), kind=mytype).ne.0._mytype) then
            tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
          endif
          if (aimag(aaa(j,m,k,3)).ne.0._mytype) then
            tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
          endif
          sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
          eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
               aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
        enddo
      enddo
      do jc=ja(i),jb(i)
        do k=spI%yst(3),spI%yen(3)
          do j=spI%yst(1),spI%yen(1) 
            aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)&
                 * real(aaa(j,m,k,jc+i), kind=mytype), aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))&
                 * aimag(aaa(j,m,k,jc+i)), kind=mytype)
          enddo
        enddo
      enddo
    enddo
  enddo
  do k=spI%yst(3),spI%yen(3)
    do j=spI%yst(1),spI%yen(1) 
      if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
        tmp1=real(aaa(j,nym,k,2), kind=mytype)/real(aaa(j,nym-1,k,3), kind=mytype)
      else
        tmp1=0._mytype
      endif
      if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
        tmp2=aimag(aaa(j,nym,k,2))/aimag(aaa(j,nym-1,k,3))
      else
        tmp2=0._mytype
      endif
      sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
      b1(j,k)=cmplx(real(aaa(j,nym,k,3), kind=mytype)-tmp1*real(aaa(j,nym-1,k,4), kind=mytype),&
           aimag(aaa(j,nym,k,3))-tmp2*aimag(aaa(j,nym-1,k,4)), kind=mytype)
      if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
        tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
        tmp3=real(eee(j,nym,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,nym-1,k), kind=mytype)
      else
        tmp1=0._mytype
        tmp3=0._mytype
      endif
      if (abs(aimag(b1(j,k))).gt.epsilon) then
        tmp2=aimag(sr(j,k))/aimag(b1(j,k))
        tmp4=aimag(eee(j,nym,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,nym-1,k))
      else
        tmp2=0._mytype
        tmp4=0._mytype
      endif
      a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
      eee(j,nym,k)=cmplx(tmp3,tmp4, kind=mytype)

      if (abs(real(aaa(j,nym-1,k,3), kind=mytype)).gt.epsilon) then
        tmp1=1._mytype/real(aaa(j,nym-1,k,3), kind=mytype)
      else
        tmp1=0._mytype
      endif
      if (abs(aimag(aaa(j,nym-1,k,3))).gt.epsilon) then
        tmp2=1._mytype/aimag(aaa(j,nym-1,k,3))
      else
        tmp2=0._mytype
      endif
      b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
      a1(j,k)=cmplx(real(aaa(j,nym-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
           aimag(aaa(j,nym-1,k,4))*aimag(b1(j,k)), kind=mytype)
      eee(j,nym-1,k)=cmplx(real(eee(j,nym-1,k), kind=mytype)*real(b1(j,k), kind=mytype)-&
           real(a1(j,k), kind=mytype)*real(eee(j,nym,k), kind=mytype),&
           aimag(eee(j,nym-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,nym,k)), kind=mytype)
    enddo
  enddo

  do i=nym-2,1,-1
    do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1) 
        if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
          tmp1=1._mytype/real(aaa(j,i,k,3), kind=mytype)
        else
          tmp1=0._mytype
        endif
        if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
          tmp2=1._mytype/aimag(aaa(j,i,k,3))
        else
          tmp2=0._mytype
        endif
        sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
        a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
             aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
        b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
             aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
        eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
             real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
             real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
             aimag(eee(j,i,k))*aimag(sr(j,k))-&
             aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
      enddo
    enddo
  enddo

  return

end subroutine inversion5_v2

!********************************************************************
!
!
!********************************************************************
subroutine channel (ux)

  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux

  integer :: j,i,k,code
  real(mytype) :: can,ut3,ut,ut4

  ut3=0.
  do k=1,ysize(3)
    do i=1,ysize(1)
      ut=0.
      do j=1,ny-1
        if (istret.ne.0) then
          ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-0.5_mytype*(ux(i,j+1,k)-ux(i,j,k)))
        else
          ut=ut+(yly/(ny-1))*(ux(i,j+1,k)-0.5_mytype*(ux(i,j+1,k)-ux(i,j,k)))
        endif
      enddo
      ut=ut/yly
      ut3=ut3+ut
    enddo
  enddo
  ut3=ut3/ysize(1)/ysize(3)

  call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
  ut4=ut4/nproc

  !can=-2.*xnu*gdt(itr) ! Poisseuille    
  can=-(2._mytype/3._mytype-ut4) ! constant flow rate

  if (nrank==0) then
    print *,nrank,'UT',ut4,can
  endif

  do k=1,ysize(3)
    do i=1,ysize(1)
      do j=2,ny-1
        ux(i,j,k)=-can+ux(i,j,k)
      enddo
    enddo
  enddo

  return
end subroutine channel

!*****************************************************************
!
!
!*****************************************************************
subroutine collect_data()
  
  USE decomp_2d
  USE decomp_2d_poisson
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  integer :: i,j,imin,ii,code
  real(mytype), dimension(200) :: ta
  integer,dimension(200) :: tai,tbi

  !TOTAL TIME FOR COLLECTION T=50000*0.01=500
  !I want 100 files
  !100 numbers from 0 to 500

  if (nrank==0) then
    call random_number(ta)
    do i=1,200
      tai(i)=int(ta(i)*2000./dt)
    enddo
    do j=1,200
      imin=999998
      ii=0
      do i=1,200
        if (tai(i).le.imin) then
          imin=tai(i)
          ii=i
        endif
      enddo
      tbi(j)=imin
      tai(ii)=999999
    enddo
    idata=tbi
    do i=1,200
      print *,i,idata(i),'VALUE COLLECTION DATA'
    enddo
  endif

  call MPI_BCAST(idata,200,MPI_INTEGER,0,MPI_COMM_WORLD,code)

end subroutine collect_data

!*****************************************************************
!  SUBROUTINE: eval_error
! DESCRIPTION: Given a (known) solution, compare with numerical
!              value and compute (normalised) l2 norm.
!*****************************************************************
SUBROUTINE eval_error(sol_num, sol_exact, name)

  USE var
  USE MPI
  USE decomp_2d
  
  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1),xsize(2),xsize(3)), INTENT(IN) :: sol_num, sol_exact
  CHARACTER(LEN=*), INTENT(IN) :: name

  REAL(mytype) :: err, errlocal
  INTEGER :: ijk
  INTEGER :: nvect1
  INTEGER :: ierr

  nvect1 = xsize(1) * xsize(2) * xsize(3)
  errlocal = 0._mytype
  DO ijk = 1,nvect1
    errlocal = errlocal + (sol_num(ijk,1,1) - sol_exact(ijk,1,1))**2
  ENDDO

  CALL MPI_ALLREDUCE(errlocal, err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
  err = err / float(nxm * nym * nzm)
  err = SQRT(err)

  IF (nrank.eq.0) THEN
    PRINT *, "Error in ", name, " = ", err
  ENDIF
  
ENDSUBROUTINE eval_error

!*****************************************************************
!  SUBROUTINE: eval_error_vel
! DESCRIPTION: Compute an exact solution for rho and compare with
!              the numerically obtained solution.
!*****************************************************************
SUBROUTINE eval_error_vel(ux1_num, uy1_num, uz1_num)

  USE var
  USE param
  USE decomp_2d
  USE MPI

  IMPLICIT NONE

  REAL(mytype), INTENT(IN), DIMENSION(xsize(1),xsize(2),xsize(3)) :: ux1_num, uy1_num, uz1_num

  REAL(mytype), DIMENSION(xsize(1),xsize(2),xsize(3)) :: ux1_exact, uy1_exact, uz1_exact
  REAL(mytype) :: x,y,z
  REAL(mytype) :: xspec,yspec,zspec
  INTEGER :: i,j,k

  ! Compute the exact solution
  DO k = 1, xsize(3)
    z = float(k + xstart(3) - 2) * dz
    zspec = (2._mytype * PI) * (z / zlz)
    DO j = 1, xsize(2)
      y = float(j + xstart(2) - 2) * dy
      yspec = (2._mytype * PI) * (y / yly)
      DO i = 1, xsize(1)
        x = float(i + xstart(1) - 2) * dx
        xspec = (2._mytype * PI) * (x / xlx)

        ux1_exact(i,j,k) = (xlx / (2._mytype * PI)) * SIN(xspec) * COS(yspec) * COS(zspec)
        uy1_exact(i,j,k) = (yly / (2._mytype * PI)) * COS(xspec) * SIN(yspec) * COS(zspec)
        uz1_exact(i,j,k) = -2._mytype * (zlz / (2._mytype * PI)) * COS(xspec) * COS(yspec) * SIN(zspec)
      ENDDO
    ENDDO
  ENDDO

  ! Compare against the numerical solution
  CALL eval_error(ux1_num, ux1_exact, "UVEL")
  CALL eval_error(uy1_num, uy1_exact, "VVEL")
  CALL eval_error(uz1_num, uz1_exact, "WVEL")
  
ENDSUBROUTINE eval_error_vel

!*****************************************************************
!  SUBROUTINE: eval_error_rho
! DESCRIPTION: Compute an exact solution for rho and compare with
!              the numerically obtained solution.
!*****************************************************************
SUBROUTINE eval_error_rho(rho_num)

  USE var
  USE param
  USE decomp_2d
  USE MPI

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1),xsize(2),xsize(3)) :: rho_num

  REAL(mytype), DIMENSION(xsize(1),xsize(2),xsize(3)) :: rho_exact
  REAL(mytype) :: x,y,z
  REAL(mytype) :: xspec,yspec,zspec
  INTEGER :: i,j,k

  ! Compute the exact solution
  DO k = 1, xsize(3)
    z = float(k + xstart(3) - 2) * dz
    zspec = (2._mytype * PI) * (z / zlz)
    DO j = 1, xsize(2)
      y = float(j + xstart(2) - 2) * dy
      yspec = (2._mytype * PI) * (y / yly)
      DO i = 1, xsize(1)
        x = float(i + xstart(1) - 2) * dx
        xspec = (2._mytype * PI) * (x / xlx)
        rho_exact(i, j, k) = 2._mytype + SIN(xspec) * SIN(yspec) * SIN(zspec)
      ENDDO
    ENDDO
  ENDDO

  ! Compare against the numerical solution
  CALL eval_error(rho_num, rho_exact, "RHO")
  
ENDSUBROUTINE eval_error_rho

!*****************************************************************
!  SUBROUTINE: updateXYZ
! DESCRIPTION: Given a set of new values in pencil 'X', update
!              by transposing to 'Y' then 'Z'.
!       INPUT: var1, the variable in 'X' stencil
!              size1, the size array of var1
!              size2, the size array of var2
!              size3, the size array of var3
!              idx1, the actual pencil of var1: 1 = X, 2 = Y, 3 = Z
!              idx2, the actual pencil of var2: 1 = X, 2 = Y, 3 = Z
!              idx3, the actual pencil of var3: 1 = X, 2 = Y, 3 = Z
!      OUTPUT: var2, var3 the variable in 'Y' and 'Z' stencils
!              with updated values
!*****************************************************************
SUBROUTINE updateXYZ(var1, var2, var3, size1, size2, size3, idx1, idx2, idx3)

  USE var
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: idx1, idx2, idx3
  INTEGER, DIMENSION(3), INTENT(IN) :: size1, size2, size3

  REAL(mytype), DIMENSION(size1(1), size1(2), size1(3)), INTENT(IN) :: var1
  REAL(mytype), DIMENSION(size2(1), size2(2), size2(3)), INTENT(OUT) :: var2
  REAL(mytype), DIMENSION(size3(1), size3(2), size3(3)), INTENT(OUT) :: var3

  ! IF((idx1.EQ.idx2).OR.(idx1.EQ.idx3).OR.(idx2.EQ.idx3)) THEN
  !   PRINT *, "Transpose: ", idx1, "->", idx2, "->" idx3, " is invalid!"
  !   STOP
  ! ENDIF

  ! idx1->idx2
  IF(idx1.EQ.1) THEN
    IF(idx2.EQ.2) THEN
      CALL transpose_x_to_y(var1, var2)
    ELSE
      PRINT *, "X->Z not yet implemented"
      STOP
      ! CALL transpose_x_to_z(var1, var2)
    ENDIF
  ELSEIF(idx1.EQ.2) THEN
    IF(idx2.EQ.1) THEN
      CALL transpose_y_to_x(var1, var2)
    ELSE
      CALL transpose_y_to_z(var1, var2)
    ENDIF
  ELSE
    IF(idx2.EQ.1) THEN
      PRINT *, "Z->X not yet implemented"
      STOP
      ! CALL transpose_z_to_x(var1, var2)
    ELSE
      CALL transpose_z_to_y(var1, var2)
    ENDIF
  ENDIF

  ! idx2->idx3
  IF(idx2.EQ.1) THEN
    IF(idx3.EQ.2) THEN
      CALL transpose_x_to_y(var2, var3)
    ELSE
      PRINT *, "X->Z not yet implemented"
      STOP
      ! CALL transpose_x_to_z(var2, var3)
    ENDIF
  ELSEIF(idx2.EQ.2) THEN
    IF(idx3.EQ.1) THEN
      CALL transpose_y_to_x(var2, var3)
    ELSE
      CALL transpose_y_to_z(var2, var3)
    ENDIF
  ELSE
    IF(idx3.EQ.1) THEN
      PRINT *, "Z->X not yet implemented"
      STOP
      ! CALL transpose_z_to_x(var2, var3)
    ELSE
      CALL transpose_z_to_y(var2, var3)
    ENDIF
  ENDIF
    
ENDSUBROUTINE updateXYZ

SUBROUTINE track_front(ux1, rho1)

  USE param
  USE variables
  USE decomp_2d
  USE MPI
  
  IMPLICIT NONE

  INTEGER :: i, j, k, r
  INTEGER :: code
  LOGICAL :: file_exists

  INTEGER status(MPI_STATUS_SIZE)
  INTEGER :: ierr
  
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, rho1
  REAL(mytype) :: x, x_left, x_right, u_left, u_right, rho_left, rho_right, rho_mid
  REAL(mytype) :: x_leftmin, u_leftmin, x_rightmax, u_rightmax

  x_left = 10._mytype * xlx
  x_right = -10._mytype * xlx

  rho_left = 0._mytype
  rho_right = 0._mytype
  DO k = 1, xsize(3)
     DO j = 1, xsize(2)
        rho_left = rho_left + rho1(1, j, k)
        rho_right = rho_right + rho1(xsize(1), j, k)
     ENDDO
  ENDDO
  rho_left = rho_left / xsize(2) / xsize(3)
  rho_right = rho_right / xsize(2) / xsize(3)
  rho_mid = 0.5_mytype * (rho_left + rho_right)

  !! Find the fronts
  DO k = 1, xsize(3)
     DO j = 1, xsize(2)
        DO i = 1, xsize(1) - 1
           x = ((i + 0.5_mytype) + xstart(1) - 2) * dx

           IF (((rho1(i, j, k).GT.rho_mid).AND.(rho1(i + 1, j, k).LT.rho_mid)) &
                .OR.((rho1(i, j, k).LT.rho_mid).AND.(rho1(i + 1, j, k).GT.rho_mid))) THEN
              IF (x.GT.x_right) THEN
                 !! Found the right front
                 x_right = x
                 u_right = 0.5_mytype * (ux1(i, j, k) + ux1(i + 1, j, k))
              ELSEIF (x.LT.x_left) THEN
                 !! Found the left front
                 x_left = x
                 u_left = 0.5_mytype * (ux1(i, j, k) + ux1(i + 1, j, k))
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  IF (nrank.EQ.0) THEN
     x_leftmin = x_left
     u_leftmin = u_left
     x_rightmax = x_right
     u_rightmax = u_right

     DO r = 1, nproc - 1
        CALL MPI_RECV(x_left, 1, real_type, r, 3000 + r, MPI_COMM_WORLD, status, ierr)
        CALL MPI_RECV(u_left, 1, real_type, r, 3000 + nproc + r, MPI_COMM_WORLD, status, ierr)
        CALL MPI_RECV(x_right, 1, real_type, r, 4000 + r, MPI_COMM_WORLD, status, ierr)
        CALL MPI_RECV(u_right, 1, real_type, r, 4000 + nproc + r, MPI_COMM_WORLD, status, ierr)
        IF (x_left.LT.x_leftmin) THEN
           x_leftmin = x_left
           u_leftmin = u_left
        ENDIF
        IF (x_right.GT.x_rightmax) THEN
           x_rightmax = x_right
           u_rightmax = u_right
        ENDIF
     ENDDO
  ELSE
     CALL MPI_SEND(x_left, 1, real_type, 0, 3000 + nrank, MPI_COMM_WORLD, ierr)
     CALL MPI_SEND(u_left, 1, real_type, 0, 3000 + nproc + nrank, MPI_COMM_WORLD, ierr)
     CALL MPI_SEND(x_right, 1, real_type, 0, 4000 + nrank, MPI_COMM_WORLD, ierr)
     CALL MPI_SEND(u_right, 1, real_type, 0, 4000 + nproc + nrank, MPI_COMM_WORLD, ierr)
  ENDIF

  IF (nrank.EQ.0) THEN
     INQUIRE(FILE="FRONTLOC.log", EXIST=file_exists)
     IF (file_exists.EQV..TRUE.) THEN
        OPEN(11, FILE="FRONTLOC.log", STATUS="old", ACTION="write", POSITION="append")
     ELSE
        OPEN(11, FILE="FRONTLOC.log", STATUS="new", ACTION="write")
        WRITE(11, *) "TIME XL UL XR UR"
     ENDIF
     WRITE(11, *) t, x_leftmin, u_leftmin, x_rightmax, u_rightmax
     CLOSE(11)
  ENDIF
  
ENDSUBROUTINE track_front

SUBROUTINE track_front_height(rho1, rho2, rho3)

  USE param
  USE variables
  USE decomp_2d
  USE MPI
  
  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: rho1
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: rho2
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: rho3

  REAL(mytype) :: face_area

  INTEGER :: i, j, k, r

  INTEGER :: N, i0, rem, ctr
  REAL(mytype), DIMENSION(nx) :: x, h
  REAL(mytype) :: densr

  INTEGER status(MPI_STATUS_SIZE)
  INTEGER :: ierr

  LOGICAL :: file_exists

  !! We want to work in y, but first need to average in z
  CALL transpose_x_to_y(rho1, rho2)
  CALL transpose_y_to_z(rho2, rho3)

  DO k = 2, zsize(3)
     DO j = 1, zsize(2)
        DO i = 1, zsize(1)
           rho3(i, j, 1) = rho3(i, j, 1) + rho3(i, j, k)
        ENDDO
     ENDDO
  ENDDO
  DO j = 1, zsize(2)
     DO i = 1, zsize(1)
        rho3(i, j, 1) = rho3(i, j, 1) / float(zsize(3))
        rho3(i, j, :) = rho3(i, j, 1)
     ENDDO
  ENDDO

  CALL transpose_z_to_y(rho3, rho2)

  !! Now compute moving average of front height
  densr = MIN(dens1, dens2) / MAX(dens1, dens2)
  x(:) = -1
  h(:) = 0._mytype
  N = 16 ! Averaging interval
  DO j = 1, ysize(2)
     i0 = 1
     ctr = 1

     IF ((j.EQ.1).OR.(j.EQ.ysize(2))) THEN
        face_area = 0.5_mytype * dy
     ELSE
        face_area = dy
     ENDIF
     
     DO WHILE ((i0 + (N - 1)).LE.ysize(1))
        IF (j.EQ.1) THEN
           x(ctr) = 0._mytype
           h(ctr) = -N * densr
           DO i = i0, i0 + (N - 1)
              x(ctr) = x(ctr) + (i + ystart(1) - 2) * dx
           ENDDO
           x(ctr) = x(ctr) / N
        ENDIF
        DO i = i0, i0 + (N - 1)
           h(ctr) = h(ctr) + rho2(i, j, 1) * face_area
        ENDDO
        i0 = i0 + N
        ctr = ctr + 1
     ENDDO

     IF (i0.LE.ysize(1)) THEN !! We have more work to do!
        rem = ysize(1) - i0
        IF (j.EQ.1) THEN
           x(ctr) = 0._mytype
           h(ctr) = -(rem + 1) * densr
           DO i = i0, ysize(1)
              x(ctr) = x(ctr) + (i + ystart(1) - 2) * dx
           ENDDO
           x(ctr) = x(ctr) / (rem + 1)
        ENDIF
        DO i = i0, ysize(1)
           h(ctr) = h(ctr) + rho2(i, j, 1) * face_area
        ENDDO
        ctr = ctr + 1
     ENDIF
  ENDDO
  ctr = ctr - 1 !! This is the number of points averaged over

  DO i = 1, ctr - 1
     h(i) = h(i) / N
  ENDDO
  rem = ysize(1) - N * (ctr - 1)
  h(ctr) = h(ctr) / rem

  h(:) = h(:) / (1._mytype - densr)

  !! Cleanup
  DO i = 1, ctr
     h(i) = MIN(h(i), yly)
     h(i) = MAX(h(i), 0._mytype)
  ENDDO

  IF (nrank.EQ.0) THEN
     !! Open datafile
     INQUIRE(FILE="FRONTHEIGHT.log", EXIST=file_exists)
     IF (file_exists) THEN
        OPEN(12, FILE="FRONTHEIGHT.log", STATUS="old", ACTION="write", POSITION="append")
     ELSE
        OPEN(12, FILE="FRONTHEIGHT.log", STATUS="new", ACTION="write")
        WRITE(12, *) "TIME X H"
     ENDIF
     
     !! First print out your own data
     DO i = 1, ctr
        WRITE(12, *) t, x(i), h(i)
     ENDDO
     
     !! Listen for data from others and print it out too
     DO r = 1, nproc - 1
        CALL MPI_RECV(x, nx, real_type, r, 1000 + r, MPI_COMM_WORLD, status, ierr)
        CALL MPI_RECV(h, nx, real_type, r, 1000 + nproc + r, MPI_COMM_WORLD, status, ierr)

        DO i = 1, nx
           IF (x(i).LT.0._mytype) THEN
              EXIT
           ENDIF
           WRITE(12, *) t, x(i), h(i)
        ENDDO
     ENDDO
     
     CLOSE(12)
  ELSE
     !! Send data back to rank 0
     CALL MPI_SEND(x, nx, real_type, 0, 1000 + nrank, MPI_COMM_WORLD, ierr)
     CALL MPI_SEND(h, nx, real_type, 0, 1000 + nproc + nrank, MPI_COMM_WORLD, ierr)
  ENDIF
  
ENDSUBROUTINE track_front_height

SUBROUTINE calc_energy_budgets(rho1, ux1, uy1, uz1, mu1, ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, &
     di1, &
     rho2, ux2, uy2, uz2, ta2, tb2, tc2, td2, te2, tf2, di2, &
     rho3, ux3, uy3, uz3, ta3, tb3, tc3, di3)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

  IMPLICIT NONE

  LOGICAL, SAVE :: firstcall = .TRUE.
  REAL(mytype), SAVE :: tlast
  INTEGER :: i, j, k, N
  INTEGER :: ierr
  INTEGER, DIMENSION(2) :: dims, dummy_coords
  LOGICAL, DIMENSION(2) :: dummy_periods

  INTEGER status(MPI_STATUS_SIZE)
  INTEGER :: i0, ctr, r

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: rho1, ux1, uy1, uz1, mu1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ta1, tb1, tc1, td1, te1, tf1, tg1, th1, &
       ti1, di1
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: rho2, ux2, uy2, uz2, ta2, tb2, tc2, &
       td2, te2, tf2, di2
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: rho3, ux3, uy3, uz3, ta3, tb3, tc3, di3

  REAL(mytype) :: deltax, deltay, deltaz, vol
  
  REAL(mytype) :: densr
  REAL(mytype) :: x, y, z
  REAL(mytype) :: KE, EP, KE1, EP1 ! Kinetic and Potential Energy
  REAL(mytype), SAVE :: ED ! Dissipated Energy
  REAL(mytype) :: ED1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: dissipation1
  REAL(mytype), DIMENSION(nx) :: dissipation_int, xpos

  densr = MIN(dens1, dens2) / MAX(dens1, dens2)

  KE = 0._mytype
  EP = 0._mytype
  
  CALL transpose_x_to_y(ux1, ux2)
  CALL transpose_x_to_y(uy1, uy2)
  CALL transpose_x_to_y(uz1, uz2)
  
  CALL transpose_y_to_z(ux2, ux3)
  CALL transpose_y_to_z(uy2, uy3)
  CALL transpose_y_to_z(uz2, uz3)

  CALL derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  CALL derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  CALL derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
  
  CALL transpose_z_to_y(ta3, td2) ! dudz
  CALL transpose_z_to_y(tb3, te2) ! dvdz
  CALL transpose_z_to_y(tc3, tf2) ! dwdz
  
  CALL dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
  CALL dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  CALL dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  
  CALL transpose_y_to_x(ta2, td1) ! dudy
  CALL transpose_y_to_x(tb2, te1) ! dvdy
  CALL transpose_y_to_x(tc2, tf1) ! dwdy
  CALL transpose_y_to_x(td2, tg1) ! dudz
  CALL transpose_y_to_x(te2, th1) ! dvdz
  CALL transpose_y_to_x(tf2, ti1) ! dwdz
  
  CALL derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  CALL derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  CALL derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  
  CALL MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
       dims, dummy_periods, dummy_coords, ierr)
  
  dissipation1(:,:,:) = 0._mytype
  DO k = 1, xsize(3)
     z = (k + xstart(3) - 2) * dz

     IF (dims(2).EQ.1) THEN
        IF ((k.EQ.1).OR.(k.EQ.xsize(3))) THEN
           deltaz = 0.5_mytype * dz
        ELSE
           deltaz = dz
        ENDIF
     ELSE
        IF ((xstart(3).EQ.1).AND.(k.EQ.1)) THEN
           deltaz = 0.5_mytype * dz
        ELSEIF ((nz - (nzm/dims(2)).EQ.xstart(3)).AND.(k.EQ.xsize(3))) THEN
           deltaz = 0.5_mytype * dz
        ELSE
           deltaz = dz
        ENDIF
     ENDIF

     DO j = 1, xsize(2)
        y = (j + xstart(2) - 2) * dy

        IF (dims(1).EQ.1) THEN
           IF ((j.EQ.1).OR.(j.EQ.xsize(2))) THEN
              deltay = 0.5_mytype * dy
           ELSE
              deltay = dy
           ENDIF
        ELSE
           IF ((xstart(2).EQ.1).AND.(j.EQ.1)) THEN
              deltay = 0.5_mytype * dy
           ELSEIF ((ny - (nym/dims(1)).EQ.xstart(2)).AND.(j.EQ.xsize(2))) THEN
              deltay = 0.5_mytype * dy
           ELSE
              deltay = dy
           ENDIF
        ENDIF

        DO i = 1, xsize(1)
           dissipation1(i, j, k) = ta1(i, j, k) * ta1(i, j, k) + ta1(i, j, k) * ta1(i, j, k) &
                + td1(i, j, k) * td1(i, j, k) + td1(i, j, k) * tb1(i, j, k) &
                + tg1(i, j, k) * tg1(i, j, k) + tg1(i, j, k) * tc1(i, j, k) &
                + tb1(i, j, k) * tb1(i, j, k) + tb1(i, j, k) * td1(i, j, k) &
                + te1(i, j, k) * te1(i, j, k) + te1(i, j, k) * te1(i, j, k) &
                + th1(i, j, k) * th1(i, j, k) + th1(i, j, k) * tf1(i, j, k) &
                + tc1(i, j, k) * tc1(i, j, k) + tc1(i, j, k) * tg1(i, j, k) &
                + tf1(i, j, k) * tf1(i, j, k) + tf1(i, j, k) * th1(i, j, k) &
                + ti1(i, j, k) * ti1(i, j, k) + ti1(i, j, k) * ti1(i, j, k)
           dissipation1(i, j, k) = dissipation1(i, j, k) &
                + ta1(i, j, k) * ta1(i, j, k) + ta1(i, j, k) * ta1(i, j, k) &
                + tb1(i, j, k) * td1(i, j, k) + tb1(i, j, k) * tb1(i, j, k) &
                + tc1(i, j, k) * tg1(i, j, k) + tc1(i, j, k) * tc1(i, j, k) &
                + td1(i, j, k) * tb1(i, j, k) + td1(i, j, k) * td1(i, j, k) &
                + te1(i, j, k) * te1(i, j, k) + te1(i, j, k) * te1(i, j, k) &
                + tf1(i, j, k) * th1(i, j, k) + tf1(i, j, k) * tf1(i, j, k) &
                + tg1(i, j, k) * tc1(i, j, k) + tg1(i, j, k) * tg1(i, j, k) &
                + th1(i, j, k) * tf1(i, j, k) + th1(i, j, k) * th1(i, j, k) &
                + ti1(i, j, k) * ti1(i, j, k) + ti1(i, j, k) * ti1(i, j, k)
           dissipation1(i, j, k) = -(2._mytype * xnu) &
                * (rho1(i, j, k) * (0.5_mytype * dissipation1(i, j, k)))

           IF ((i.EQ.1).OR.(i.EQ.xsize(1))) THEN
              deltax = 0.5_mytype * dx
           ELSE
              deltax = dx
           ENDIF

           vol = deltax * deltay * deltaz

           KE = KE + (ux1(i, j, k)**2 + uy1(i, j, k)**2 + uz1(i, j, k)**2) * vol
           EP = EP + (rho1(i, j, k) * y**2) * vol

           IF (firstcall) THEN
              ED = 0._mytype
           ELSE
              ED = ED + dissipation1(i, j, k) * vol * (0.5_mytype * (t - tlast))
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  KE = 0.5_mytype * KE
  EP = 0.5_mytype * EP / (1._mytype - densr) ! The 1/2 comes from int rho y dy = 1/2 rho y^2|^yly_0

  CALL MPI_ALLREDUCE(KE, KE1, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(EP, EP1, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(ED, ED1, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)

  !! Write out energy budgets
  IF (nrank.EQ.0) THEN
     IF (firstcall) THEN
        OPEN(13, FILE="ENERGYBUDGET.log", STATUS="new", ACTION="write")
        WRITE(13, *) "TIME KE EP ED"
     ELSE
        OPEN(13, FILE="ENERGYBUDGET.log", STATUS="old", ACTION="write", POSITION="append")
     ENDIF
     WRITE(13, *) t, KE1, EP1, ED1
     CLOSE(13)
  ENDIF

  !! Integrate dissipation rate vertically (after depth averaging)
  CALL transpose_x_to_y(dissipation1, ta1)
  CALL transpose_y_to_z(ta2, ta3)

  DO k = 2, zsize(3)
     DO j = 1, zsize(2)
        DO i = 1, zsize(1)
           ta3(i, j, 1) = ta3(i, j, 1) + ta3(i, j, k)
        ENDDO
     ENDDO
  ENDDO
  DO j = 1, zsize(2)
     DO i = 1, zsize(1)
        ta3(i, j, 1) = ta3(i, j, 1) / float(zsize(3))
        ta3(i, j, :) = ta3(i, j, 1)
     ENDDO
  ENDDO

  CALL transpose_z_to_y(ta3, ta2)

  xpos(:) = -1
  dissipation_int(:) = 0._mytype
  N = 16 ! Averaging interval
  DO j = 1, ysize(2)
     i0 = 1
     ctr = 1

     IF ((j.EQ.1).OR.(j.EQ.ysize(2))) THEN
        deltay = 0.5_mytype * dy
     ELSE
        deltay = dy
     ENDIF

     DO WHILE ((i0 + (N - 1)).LE.ysize(1))
        IF (j.EQ.1) THEN
           xpos(ctr) = 0._mytype
           DO i = i0, i0 + (N - 1)
              xpos(ctr) = xpos(ctr) + (i + ystart(1) - 2) * dx
           ENDDO
           xpos(ctr) = xpos(ctr) / N
        ENDIF
        DO i = i0, i0 + (N - 1)
           dissipation_int(ctr) = dissipation_int(ctr) + ta2(i, j, 1) * deltay
        ENDDO
        i0 = i0 + N
        ctr = ctr + 1
     ENDDO

     IF (i0.LE.ysize(1)) THEN
        IF (j.EQ.1) THEN
           xpos(ctr) = 0._mytype
           DO i = i0, ysize(1)
              xpos(ctr) = xpos(ctr) + (i + ystart(1) - 2) * dx
           ENDDO
           xpos(ctr) = xpos(ctr) / (ysize(1) - i0 + 1)
        ENDIF
        DO i = i0, ysize(1)
           dissipation_int(ctr) = dissipation_int(ctr) + ta2(i, j, 1) * deltay
        ENDDO
        ctr = ctr + 1
     ENDIF
  ENDDO
  ctr = ctr - 1

  DO i = 1, ctr - 1
     dissipation_int(i) = dissipation_int(i) / N
  ENDDO
  dissipation_int(ctr) = dissipation_int(ctr) / (ysize(1) - N * (ctr - 1))

  IF (nrank.EQ.0) THEN
     IF (firstcall) THEN
        OPEN(14, FILE="DISSIPATION.log", STATUS="new", ACTION="write")
        WRITE(14, *) "TIME X EPSILON"
     ELSE
        OPEN(14, FILE="DISSIPATION.log", STATUS="old", ACTION="write", POSITION="append")
     ENDIF

     DO i = 1, ctr
        WRITE(14, *) t, xpos(i), dissipation_int(i)
     ENDDO

     DO r = 1, nproc - 1
        CALL MPI_RECV(xpos, nx, real_type, r, 2000 + r, MPI_COMM_WORLD, status, ierr)
        CALL MPI_RECV(dissipation_int, nx, real_type, r, 2000 + nproc + r, MPI_COMM_WORLD, status, &
             ierr)

        DO i = 1, nx
           IF (xpos(i).LT.0._mytype) THEN
              EXIT
           ENDIF
           WRITE(14, *) t, xpos(i), dissipation_int(i)
        ENDDO
     ENDDO

     CLOSE(14)
  ELSE
     CALL MPI_SEND(xpos, nx, real_type, 0, 2000 + nrank, MPI_COMM_WORLD, ierr)
     CALL MPI_SEND(dissipation_int, nx, real_type, 0, 2000 + nproc + nrank, MPI_COMM_WORLD, ierr)
  ENDIF

  firstcall = .FALSE.

  !! Record time function called for later use in determining interval
  tlast = t
  
ENDSUBROUTINE calc_energy_budgets

SUBROUTINE force_variable_2d(ta1, ta2, ta3)

  USE param
  USE variables
  USE decomp_2d

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ta1
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: ta2
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: ta3

  INTEGER :: i, j, k

  CALL transpose_x_to_y(ta1, ta2)
  CALL transpose_y_to_z(ta2, ta3)
  DO k = 2, zsize(3)
     DO j = 1, zsize(2)
        DO i = 1, zsize(1)
           ta3(i, j, k) = ta3(i, j, 1)
        ENDDO
     ENDDO
  ENDDO
  CALL transpose_z_to_y(ta3, ta2)
  CALL transpose_y_to_x(ta2, ta1)
  
ENDSUBROUTINE force_variable_2d

SUBROUTINE calc_sedimentation(rho1, rho2, rho3)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: rho1
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: rho2
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: rho3
  REAL(mytype), DIMENSION(nx) :: D, Dtmp, xarr, xtmp

  REAL(mytype) :: us, c, mdot, mdotglob
  REAL(mytype) :: rhomax, rhomin, gamma
  REAL(mytype) :: x, z

  INTEGER :: i, j, k, r

  INTEGER status(MPI_STATUS_SIZE)
  INTEGER :: ierr
  INTEGER, DIMENSION(2) :: dims, dummy_coords
  LOGICAL, DIMENSION(2) :: dummy_periods

  LOGICAL :: file_exists

  us = 0._mytype ! Sedimenting velocity

  rhomax = MAX(dens1, dens2)
  rhomin = MIN(dens1, dens2)
  gamma = rhomax / rhomin

  CALL transpose_x_to_y(rho1, rho2)
  CALL transpose_y_to_z(rho2, rho3)
  
  CALL MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, ierr)

  mdot = 0._mytype
  D(:) = 0._mytype
  xarr(:) = 0._mytype
  
  j = 1
  IF ((dims(1).EQ.1).OR.(zstart(2).EQ.1)) THEN
     DO k = 1, zsize(3)
        z = float(k + zstart(3) - 2) * dz
        DO i = 1, zsize(1)
           x = float(i + zstart(1) - 2) * dx
           xarr(i) = x
           
           c = ((rho3(i, j, k) / rhomin) - 1._mytype) / (gamma - 1._mytype)

           D(i) = D(i) + c * us * dz
        ENDDO
     ENDDO
  ENDIF

  DO i = 1, zsize(1)
     mdot = mdot + D(i) * dx
  ENDDO

  mdot = mdot / (xlx * zlz)
  D(:) = D(:) / zlz

  CALL MPI_ALLREDUCE(mdot, mdotglob, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)

  IF (nrank.EQ.0) THEN
     IF ((dims(1).EQ.1).OR.(zstart(2).EQ.1)) THEN
        i = zsize(1)
     ELSE
        i = 0
     ENDIF

     DO r = 1, nproc - 1
        !! Does procr have info?
        CALL MPI_RECV(j, 1, MPI_INTEGER, r, 1000 + r, MPI_COMM_WORLD, status, ierr)
        
        !! If so, get that info
        IF (j.GT.0) THEN
           CALL MPI_RECV(Dtmp, j, real_type, r, 2000 + r, MPI_COMM_WORLD, status, ierr)
           CALL MPI_RECV(xtmp, j, real_type, r, 3000 + r, MPI_COMM_WORLD, status, ierr)
           DO k = 1, j
              D(i + k) = Dtmp(k)
              xarr(i + k) = xtmp(k)
           ENDDO
           i = i + j
        ENDIF
     ENDDO
     
     !! Write data to file
     INQUIRE(FILE="ms.log", EXIST=file_exists)
     IF (file_exists.EQV..TRUE.) THEN
        OPEN(20, FILE="ms.log", STATUS="old", ACTION="write", POSITION="append")
     ELSE
        OPEN(20, FILE="ms.log", STATUS="new", ACTION="write")
        WRITE(20, *) "| ITIME | ms |"
     ENDIF
     WRITE(20, *) itime, mdotglob
     CLOSE(20)
     
     INQUIRE(FILE="deposit.log", EXIST=file_exists)
     IF (file_exists.EQV..TRUE.) THEN
        OPEN(21, FILE="deposit.log", STATUS="old", ACTION="write", POSITION="append")
     ELSE
        OPEN(21, FILE="deposit.log", STATUS="new", ACTION="write")
        WRITE(21, *) "| ITIME | x | deposit |"
     ENDIF
     DO i = 1, nx
        WRITE(21, *) itime, xarr(i), D(i)
     ENDDO
     CLOSE(21)
  ELSE
     IF ((dims(1).EQ.1).OR.(zstart(2).EQ.1)) THEN
        !! Tell proc0 we have info
        CALL MPI_SEND(zsize(1), 1, MPI_INTEGER, 0, 1000 + nrank, MPI_COMM_WORLD, ierr)
        
        !! Send D array to proc0
        CALL MPI_SEND(D, zsize(1), real_type, 0, 2000 + nrank, MPI_COMM_WORLD, ierr)
        CALL MPI_SEND(xarr, zsize(1), real_type, 0, 3000 + nrank, MPI_COMM_WORLD, ierr)
     ELSE
        !! Tell proc0 we have no info
        CALL MPI_SEND(0, 1, MPI_INTEGER, 0, 1000 + nrank, MPI_COMM_WORLD, ierr)
     ENDIF
  ENDIF
  
ENDSUBROUTINE calc_sedimentation
