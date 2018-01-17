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

!********************************************************************
!
! 
!********************************************************************
subroutine intt (ux,uy,uz,gx,gy,gz,hx,hy,hz,ta1,tb1,tc1,rho)

  USE param
  USE variables
  USE decomp_2d

  implicit none

  integer :: ijk,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: rho
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx,gy,gz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx,hy,hz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1

  nxyz=xsize(1)*xsize(2)*xsize(3)

  if (ilmn.ne.0) then
    !! First, convert velocity to momentum
    if (iskew.ne.2) then
      !! Rotational form or Quasi skew-symmetric
      do ijk = 1, nxyz
        ux(ijk, 1, 1) = rho(ijk, 1, 1) * ux(ijk, 1, 1)
        uy(ijk, 1, 1) = rho(ijk, 1, 1) * uy(ijk, 1, 1)
        uz(ijk, 1, 1) = rho(ijk, 1, 1) * uz(ijk, 1, 1)
      enddo
    else
      !! Skew-symmetric
      ux(:,:,:) = SQRT(rho(:,:,:)) * ux(:,:,:)
      uy(:,:,:) = SQRT(rho(:,:,:)) * uy(:,:,:)
      uz(:,:,:) = SQRT(rho(:,:,:)) * uz(:,:,:)
    endif
  endif

  if ((nscheme.eq.1).or.(nscheme.eq.2)) then
    !! AB2 or RK3
    
    if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
         (nscheme.eq.2.and.itr.eq.1)) then
      do ijk=1,nxyz
        ux(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+ux(ijk,1,1)
        uy(ijk,1,1)=gdt(itr)*tb1(ijk,1,1)+uy(ijk,1,1) 
        uz(ijk,1,1)=gdt(itr)*tc1(ijk,1,1)+uz(ijk,1,1)
        gx(ijk,1,1)=ta1(ijk,1,1)
        gy(ijk,1,1)=tb1(ijk,1,1)
        gz(ijk,1,1)=tc1(ijk,1,1)            
      enddo
    else
      if (nz.gt.1) then
        do ijk=1,nxyz
          ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
          uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)   
          uz(ijk,1,1)=adt(itr)*tc1(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+uz(ijk,1,1)
          gx(ijk,1,1)=ta1(ijk,1,1)
          gy(ijk,1,1)=tb1(ijk,1,1)
          gz(ijk,1,1)=tc1(ijk,1,1)            
        enddo
      else !! End is 3D
        do ijk=1,nxyz
          ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
          uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)   
          gx(ijk,1,1)=ta1(ijk,1,1)
          gy(ijk,1,1)=tb1(ijk,1,1)
        enddo
      endif !! End is 2D
    endif
  else if (nscheme.eq.3) then 
    if (nz.gt.1) then
      ! if (adt(itr)==0._mytype) then
      if (itr.eq.0) then ! XXX The above double comparison is only true for itr=0
        do ijk=1,nxyz
          gx(ijk,1,1)=dt*ta1(ijk,1,1)
          gy(ijk,1,1)=dt*tb1(ijk,1,1)
          gz(ijk,1,1)=dt*tc1(ijk,1,1)
        enddo
      else
        do ijk=1,nxyz
          gx(ijk,1,1)=adt(itr)*gx(ijk,1,1)+dt*ta1(ijk,1,1)
          gy(ijk,1,1)=adt(itr)*gy(ijk,1,1)+dt*tb1(ijk,1,1)
          gz(ijk,1,1)=adt(itr)*gz(ijk,1,1)+dt*tc1(ijk,1,1)
        enddo
      endif
      do ijk=1,nxyz
        ux(ijk,1,1)=ux(ijk,1,1)+bdt(itr)*gx(ijk,1,1)
        uy(ijk,1,1)=uy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)
        uz(ijk,1,1)=uz(ijk,1,1)+bdt(itr)*gz(ijk,1,1)
      enddo
    else
      ! if (adt(itr)==0._mytype) then
      if (itr.eq.0) then
        do ijk=1,nxyz
          gx(ijk,1,1)=dt*ta1(ijk,1,1)
          gy(ijk,1,1)=dt*tb1(ijk,1,1)
        enddo
      else
        do ijk=1,nxyz
          gx(ijk,1,1)=adt(itr)*gx(ijk,1,1)+dt*ta1(ijk,1,1)
          gy(ijk,1,1)=adt(itr)*gy(ijk,1,1)+dt*tb1(ijk,1,1)
        enddo
      endif
      do ijk=1,nxyz
        ux(ijk,1,1)=ux(ijk,1,1)+bdt(itr)*gx(ijk,1,1)
        uy(ijk,1,1)=uy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)
      enddo
    endif
  else if (nscheme==4) then
    if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) then
        print *,'start with Euler',itime
      endif

      do ijk=1,nxyz !start with Euler
        ux(ijk,1,1)=dt*ta1(ijk,1,1)+ux(ijk,1,1)
        uy(ijk,1,1)=dt*tb1(ijk,1,1)+uy(ijk,1,1) 
        uz(ijk,1,1)=dt*tc1(ijk,1,1)+uz(ijk,1,1)
        gx(ijk,1,1)=ta1(ijk,1,1)
        gy(ijk,1,1)=tb1(ijk,1,1)
        gz(ijk,1,1)=tc1(ijk,1,1)            
      enddo
    else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
        if (nrank==0) then
          print *,'then with AB2',itime
        endif
        
        do ijk=1,nxyz
          ux(ijk,1,1)=1.5_mytype*dt*ta1(ijk,1,1)-0.5_mytype*dt*gx(ijk,1,1)+ux(ijk,1,1)
          uy(ijk,1,1)=1.5_mytype*dt*tb1(ijk,1,1)-0.5_mytype*dt*gy(ijk,1,1)+uy(ijk,1,1)
          uz(ijk,1,1)=1.5_mytype*dt*tc1(ijk,1,1)-0.5_mytype*dt*gz(ijk,1,1)+uz(ijk,1,1)
          hx(ijk,1,1)=gx(ijk,1,1)
          hy(ijk,1,1)=gy(ijk,1,1)
          hz(ijk,1,1)=gz(ijk,1,1)
          gx(ijk,1,1)=ta1(ijk,1,1)
          gy(ijk,1,1)=tb1(ijk,1,1)
          gz(ijk,1,1)=tc1(ijk,1,1)
        enddo
      else
        do ijk=1,nxyz
          ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+&
               cdt(itr)*hx(ijk,1,1)+ux(ijk,1,1)
          uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+&
               cdt(itr)*hy(ijk,1,1)+uy(ijk,1,1)
          uz(ijk,1,1)=adt(itr)*tc1(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+&
               cdt(itr)*hz(ijk,1,1)+uz(ijk,1,1)
          hx(ijk,1,1)=gx(ijk,1,1)
          hy(ijk,1,1)=gy(ijk,1,1)
          hz(ijk,1,1)=gz(ijk,1,1)
          gx(ijk,1,1)=ta1(ijk,1,1)
          gy(ijk,1,1)=tb1(ijk,1,1)
          gz(ijk,1,1)=tc1(ijk,1,1)
        enddo
      endif
    endif
  endif


  return
end subroutine intt

!********************************************************************
!********************************************************************
SUBROUTINE inttdensity(rho1, rhos1, rhoss1, rhos01, tg1, drhodt1)

  USE param
  USE variables
  USE decomp_2d

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: tg1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(OUT) :: rhos01, drhodt1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(INOUT) :: rho1, rhos1, rhoss1

  IF ((nscheme.EQ.1).OR.(nscheme.EQ.2)) THEN
    !! AB2 or RK3

    ! First store -rho1 in drhodt1 incase we use simple extrapolation
    drhodt1(:,:,:) = -rho1(:,:,:)

    IF (nscheme.EQ.1) THEN
      !! AB2
      rhos01(:,:,:) = rhoss1(:,:,:)
      rhoss1(:,:,:) = rho1(:,:,:)
    ENDIF

    IF ((nscheme.EQ.1.AND.itime.EQ.1.AND.ilit.EQ.0).OR.&
         (nscheme.EQ.2.AND.itr.EQ.1)) THEN
      rho1(:,:,:) = rho1(:,:,:) + gdt(itr) * tg1(:,:,:)

      IF (nscheme.EQ.2) THEN
        !! RK3
        rhos01(:,:,:) = rhoss1(:,:,:)
        rhoss1(:,:,:) = tg1(:,:,:)
      ENDIF
    ELSE
      rho1(:,:,:) = rho1(:,:,:) + adt(itr) * tg1(:,:,:) &
           + bdt(itr) * rhos1(:,:,:)
    ENDIF
  ELSE IF (nscheme.EQ.3) THEN
    !! RK4
    !! XXX Not implemented!
    IF (nrank.EQ.0) THEN
      PRINT  *, 'LMN: RK4 not ready'
    ENDIF
    STOP
  ELSE
    !! AB3
    IF ((itime.EQ.1).AND.(ilit.EQ.0)) THEN
      IF (nrank.EQ.0) THEN
        PRINT  *, 'start with Euler', itime
      ENDIF
      rho1(:,:,:) = rho1(:,:,:) + dt * tg1(:,:,:)
    ELSE
      IF  ((itime.EQ.2).AND.(ilit.EQ.0)) THEN
        IF (nrank.EQ.0) THEN
          PRINT *, 'then with AB2', itime
        ENDIF
        rho1(:,:,:) = rho1(:,:,:) - 0.5_mytype * dt * (rhos1(:,:,:) - 3._mytype * tg1(:,:,:))
      ELSE
        rho1(:,:,:) = rho1(:,:,:) + adt(itr) * tg1(:,:,:) + bdt(itr) * rhos1(:,:,:) + cdt(itr) &
             * rhoss1(:,:,:)
      ENDIF

      !! Update oldold stage
      rhoss1(:,:,:) = rhos1(:,:,:)
    ENDIF
  ENDIF

  !! Update old stage
  rhos1(:,:,:) = tg1(:,:,:)

ENDSUBROUTINE inttdensity

!********************************************************************
!
! 
!********************************************************************
subroutine corgp (ux,gx,uy,uz,px,py,pz,rho)

  USE decomp_2d
  USE variables
  USE param
  USE var
  USE MPI

  implicit none

  integer :: ijk,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz,rho
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gx

  real(mytype) :: invrho

  nxyz=xsize(1)*xsize(2)*xsize(3)

  if (ilmn.ne.0) then
    if (nrhoscheme.ne.0) then
      !! We are solving constant-coefficient Poisson equation,
      !! first convert momentum->velocity
      if (iskew.ne.2) then
        !! Rotational or quasi skew-symmetric
        do ijk = 1, nxyz
          invrho = 1._mytype / rho(ijk, 1, 1)
          ux(ijk, 1, 1) = ux(ijk, 1, 1) * invrho
          uy(ijk, 1, 1) = uy(ijk, 1, 1) * invrho
          uz(ijk, 1, 1) = uz(ijk, 1, 1) * invrho
        enddo
      else
        !! Skew-symmetric
        do ijk = 1, nxyz
          invrho = 1._mytype / SQRT(rho(ijk, 1, 1))
          ux(ijk, 1, 1) = ux(ijk, 1, 1) * invrho
          uy(ijk, 1, 1) = uy(ijk, 1, 1) * invrho
          uz(ijk, 1, 1) = uz(ijk, 1, 1) * invrho
        enddo
      endif
    endif
    
    if (iskew.ne.2) then
      !! Rotational or quasi skew-symmetric
      do ijk=1, nxyz
        invrho = 1._mytype / rho(ijk, 1, 1)
        ux(ijk, 1, 1) = ux(ijk, 1, 1) - invrho * px(ijk, 1, 1)
        uy(ijk, 1, 1) = uy(ijk, 1, 1) - invrho * py(ijk, 1, 1)
        uz(ijk, 1, 1) = uz(ijk, 1, 1) - invrho * pz(ijk, 1, 1)
      enddo
    else
      !! Skew-symmetric
      do ijk = 1, nxyz
        invrho = 1._mytype / SQRT(rho(ijk, 1, 1))
        ux(ijk, 1, 1) = ux(ijk, 1, 1) - invrho * px(ijk, 1, 1)
        uy(ijk, 1, 1) = uy(ijk, 1, 1) - invrho * py(ijk, 1, 1)
        uz(ijk, 1, 1) = uz(ijk, 1, 1) - invrho * pz(ijk, 1, 1)
      enddo
    endif
  else
    ux(:,:,:) = -px(:,:,:) + ux(:,:,:)
    uy(:,:,:) = -py(:,:,:) + uy(:,:,:)
    uz(:,:,:) = -pz(:,:,:) + uz(:,:,:)
  endif

  if (itype==2) then !channel flow
    call transpose_x_to_y(ux,gx)
    call channel(gx)
    call transpose_y_to_x(gx,ux)
  endif

  return
end subroutine corgp

!*********************************************************
!
!*********************************************************
subroutine inflow (ux, uy, uz, rho, phi)

  USE param
  USE IBM
  USE variables
  USE decomp_2d

  implicit none

  integer :: k, j
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz, rho, phi
  real(mytype) :: r1, r2, r3, y, um

  call ecoule(ux, uy, uz, rho)

  call random_number(bxo)
  call random_number(byo)
  call random_number(bzo)

  if (iin.eq.1) then  
    do k = 1, xsize(3)
      do j = 1, xsize(2)
        bxx1(j, k) = bxx1(j, k)+bxo(j, k) * noise1
        bxy1(j, k) = bxy1(j, k)+byo(j, k) * noise1
        bxz1(j, k) = bxz1(j, k)+bzo(j, k) * noise1
        rho(1, j, k) = 1._mytype
      enddo
    enddo

    if (iscalar==1) then
      do k = 1, xsize(3)
        do j = 1, xsize(2)
          phi(1, j, k) = 1._mytype
        enddo
      enddo
    endif
  endif

  return
end subroutine inflow

!*********************************************************
!
!*********************************************************
subroutine outflow (ux, uy, uz, rho, phi)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

  implicit none

  integer :: j, k, i,  code
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz, rho, phi
  real(mytype) :: udx, udy, udz, uddx, uddy, uddz, uxmax, &
       uxmin, vphase, cx, coef, uxmax1, uxmin1, volflux
  real(mytype) :: Ay

  !! Compute 'convective velocity' at outlet
  udx = 1._mytype / dx
  udy = 1._mytype / dy
  udz = 1._mytype / dz
  uddx = 0.5_mytype / dx
  uddy = 0.5_mytype / dy
  uddz = 0.5_mytype / dz

  ! If inlet velocity specified in terms of u1 and u2
  cx = 0.5_mytype * (u1 + u2) * gdt(itr) * udx

  uxmax = -1609._mytype
  uxmin = 1609._mytype
  do k = 1, xsize(3)
    do j = 1, xsize(2)
      if (ux(nx - 1, j, k).gt.uxmax) uxmax = ux(nx - 1, j, k)
      if (ux(nx - 1, j, k).lt.uxmin) uxmin = ux(nx - 1, j, k)
    enddo
  enddo
  call MPI_ALLREDUCE(uxmax, uxmax1, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE(uxmin, uxmin1, 1, real_type, MPI_MIN, MPI_COMM_WORLD, code)
  vphase = 0.5_mytype * (uxmax1 + uxmin1)
  cx = vphase * gdt(itr) * udx

  !! Compute inlet volume-flux
  volflux = 0._mytype
  do k = 1, xsize(3)
    do j = 1, xsize(2) - 1
      if (istret.eq.0) then
        Ay = yly / (ny - 1)
      else
        Ay = (yp(j + 1) - yp(j))
      endif
      volflux = volflux + 0.5_mytype * (ux(1, j, k) + ux(1, j + 1, k)) * Ay
    enddo
  enddo
  volflux = volflux / (yly * xsize(3))
  call MPI_ALLREDUCE(MPI_IN_PLACE, volflux, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
  cx = volflux * gdt(itr) * udx

  if (itype.ne.9) then
    !! XXX Update density last as we are imposing
    !
    !        ddt (rho phi) + C ddn (rho phi) = 0
    !      phi = 1, u, etc. so that
    !
    !        phi^{k+1} = (rho^k -dt C ddn (rho phi)^k) / rho^{k+1}
    do k = 1, xsize(3)
      do j = 1, xsize(2)
        bxxn(j, k) = rho(nx, j, k) * ux(nx, j, k) &
             - cx * (rho(nx, j, k) * ux(nx, j, k) - rho(nx - 1, j, k) * ux(nx - 1, j, k))
        bxyn(j, k) = rho(nx, j, k) * uy(nx, j, k) &
             - cx * (rho(nx, j, k) * uy(nx, j, k) - rho(nx - 1, j, k) * uy(nx - 1, j, k))
        bxzn(j, k) = rho(nx, j, k) * uz(nx, j, k) &
             - cx * (rho(nx, j, k) * uz(nx, j, k) - rho(nx - 1, j, k) * uz(nx - 1, j, k))
      enddo
    enddo

    if (iscalar.eq.1) then
      do k = 1, xsize(3)
        do j = 1, xsize(2)
          phi(nx, j, k) = rho(nx, j, k) * phi(nx, j, k) &
               - cx * (rho(nx, j, k) * phi(nx, j, k) - rho(nx - 1, j, k) * phi(nx - 1, j, k))
        enddo
      enddo
    endif

    ! Update density and other variables
    do k = 1, xsize(3)
      do j = 1, xsize(2)
        rho(nx, j, k) = rho(nx, j, k) - cx * (rho(nx, j, k) - rho(nx - 1, j, k))

        bxxn(j, k) = bxxn(j, k) / rho(nx, j, k)
        bxyn(j, k) = bxyn(j, k) / rho(nx, j, k)
        bxzn(j, k) = bxzn(j, k) / rho(nx, j, k)
      enddo
    enddo
    if (iscalar.eq.1) then
      do k = 1, xsize(3)
        do j = 1, xsize(2)
          phi(nx, j, k) = phi(nx, j, k) / rho(nx, j, k)
        enddo
      enddo
    endif
  else
    print *, 'NOT READY'
    stop
  endif

  return
end subroutine outflow

!**********************************************************************
!
!
!**********************************************************************
subroutine ecoule(ux1,uy1,uz1,rho1)

  USE param
  USE IBM
  USE variables
  USE decomp_2d

  implicit none

  integer  :: i,j,k,jj1,jj2 
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,rho1
  real(mytype) :: x,y,z,ym
  real(mytype) :: xspec,yspec,zspec
  real(mytype) :: r1,r2,r3,r
  real(mytype) :: uh,ud,um,xv,bruit1
  real(mytype) :: u_disturb, v_disturb, disturb_decay

  bxx1=0._mytype;bxy1=0._mytype;bxz1=0._mytype
  byx1=0._mytype;byy1=0._mytype;byz1=0._mytype
  bzx1=0._mytype;bzy1=0._mytype;bzz1=0._mytype 

  !ITYPE=1 --> Constant flow field
  !ITYPE=2 --> Channel flow
  !ITYPE=3 --> Wake flow
  !ITYPE=4 --> Mixing layer with splitter plate
  !ITYPE=5 --> Channel flow
  !ITYPE=6 --> Taylor Green vortices
  !ITYPE=7 --> Cavity flow
  !ITYPE=8 --> Flat plate Boundary layer
  !ITYPE=9 --> Tank 

  if (itype.eq.1) then
    do k=1,xsize(3)
      z = float(k + xstart(3) - 2) * dz
      zspec = (2._mytype * PI) * (z / zlz)
      do j=1,xsize(2)
        y = float(j + xstart(2) - 2) * dy
        yspec = (2._mytype * PI) * (y / yly)
        do i=1,xsize(1)
          x = float(i + xstart(1) - 2) * dx
          xspec = (2._mytype * PI) * (x / xlx)

          ux1(i,j,k) = (xlx / (2._mytype * PI)) * SIN(xspec) * COS(yspec) * COS(zspec)
          uy1(i,j,k) = (yly / (2._mytype * PI)) * COS(xspec) * SIN(yspec) * COS(zspec)
          uz1(i,j,k) = -2._mytype * (zlz / (2._mytype * PI)) * COS(xspec) * COS(yspec) * SIN(zspec)
        enddo
      enddo
    enddo
  else if (itype.eq.2) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (istret.eq.0) then
          y=(j+xstart(2)-1-1)*dy-yly/2._mytype
        else
          y=yp(j+xstart(2)-1)-yly/2._mytype
        endif
        do i=1,xsize(1)
          ux1(i,j,k)=ux1(i,j,k)+1._mytype-y*y
        enddo
      enddo
    enddo
  else if (itype.eq.3) then
    if (nrank.eq.0) then
      PRINT *, "itype=", itype, " not implemented!"
      STOP
    endif
  else if (itype.eq.4) then
    ! Mixing layer flow

    ! ! Check the BCs
    ! if(nclx.eq.2) then
    !   print *, "Please set x-BC to 0 or 1 (periodic or free-slip)"
    !   stop
    ! endif
    ! if(ncly.ne.1) then
    !   print *, "Please set y-BC to 1 (free-slip)"
    !   stop
    ! endif
    ! if(nclz.ne.0) then
    !   print *, "Please set z-BC to 0 (periodic)"
    !   stop
    ! endif

    ! #ifndef TWOD
    ! Set the flowfield
    do k = 1, xsize(3)
      do j = 1, xsize(2)
        if (istret.eq.0) then
          y = float((j + xstart(2) - 2)) * dy - yly / 2._mytype
        else
          y = yp(j + xstart(2) - 1) - yly / 2._mytype
        endif
        do i = 1, xsize(1)
          x = float((i + xstart(1) - 2)) * dx

          ! Set mean field
          ux1(i, j, k) = ux1(i, j, k) + (u1 + u2) / 2._mytype &
               + (u1 - u2) * TANH(2._mytype * y) / 2._mytype

          ! Calculate disturbance field (as given in Fortune2004)
          ! NB x and y are swapped relative to Fortune2004
          ! `noise' is used to set the intensity of the disturbance
          disturb_decay = noise * (u1 - u2) * EXP(-0.05_mytype * (y**2))
          u_disturb = disturb_decay * (SIN(8._mytype * PI * x / xlx) &
               + SIN(4._mytype * PI * x / xlx) / 8._mytype &
               + SIN(2._mytype * PI * x / xlx) / 16._mytype)
          u_disturb = (0.05_mytype * y * xlx / PI) * u_disturb
          v_disturb = disturb_decay * (COS(8._mytype * PI * x / xlx) &
               + COS(4._mytype * PI * x / xlx) / 8._mytype &
               + COS(2._mytype * PI * x / xlx) / 16._mytype)

          ux1(i, j, k) = ux1(i, j, k) + u_disturb
          uy1(i, j, k) = uy1(i, j, k) + v_disturb 
          uz1(i, j, k) = uz1(i, j, k) + 0._mytype

          if (y.gt.0._mytype) then
            rho1(i, j, k) = rho1(i, j, k) + dens1
          else if (y.lt.0._mytype) then
            rho1(i, j, k) = rho1(i, j, k) + dens2
          else
            rho1(i, j, k) = rho1(i, j, k) + 0.5_mytype * (dens1 + dens2)
          endif
        enddo
      enddo
    enddo
    ! #else
    !   ! 2D mixing layer
    ! #endif
  else if (itype.eq.5) then
    if (nclx.ne.0) then
      print *,'NOT POSSIBLE'
      stop
    endif
    if (nclz.ne.0) then
      print *,'NOT POSSIBLE'
      stop
    endif
    if (ncly==0) then
      print *,'NOT POSSIBLE'
      stop
    endif
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/2._mytype
        if (istret.ne.0) y=yp(j)-yly/2._mytype
        do i=1,xsize(1)
          ux1(i,j,k)=ux1(i,j,k)+1._mytype-y*y
        enddo
      enddo
    enddo
  else if (itype.eq.6) then
    t=0._mytype
    !xv=1._mytype/100._mytype
    !xxk1=twopi/xlx
    !xxk2=twopi/yly
    do k=1,xsize(3)
      z=float((k+xstart(3)-1-1))*dz
      do j=1,xsize(2)
        y=float((j+xstart(2)-1-1))*dy
        do i=1,xsize(1)
          x=float(i-1)*dx
          ux1(i,j,k)=+sin(x)*cos(y)*cos(z)
          uy1(i,j,k)=-cos(x)*sin(y)*cos(z)
          uz1(i,j,k)=0._mytype
          bxx1(j,k)=0._mytype
          bxy1(j,k)=0._mytype
          bxz1(j,k)=0._mytype
        enddo
      enddo
    enddo
!!! CM    call test_min_max('ux1  ','In intt        ',ux1,size(ux1))
!!! CM    call test_min_max('uy1  ','In intt        ',uy1,size(uy1))
  else if (itype.eq.7) then
    if (nrank.eq.0) then
      PRINT *, "itype=", itype, " not implemented!"
      STOP
    endif
  else if (itype.eq.8) then
    if (nrank.eq.0) then
      PRINT *, "itype=", itype, " not implemented!"
      STOP
    endif
  else if (itype.eq.9) then
    if (nrank.eq.0) then
      PRINT *, "itype=", itype, " not implemented!"
      STOP
    endif
  else if (itype.eq.10) then
    do k=1,xsize(3)
      do j=1,xsize(2)
        bxx1(j,k)=0._mytype
        bxy1(j,k)=0._mytype
        bxz1(j,k)=0._mytype
      enddo
    enddo
  else
  endif
    if (nrank.eq.0) then
      PRINT *, "itype=", itype, " not implemented!"
      STOP
    endif
  return
end subroutine ecoule

!********************************************************************
!
!
!********************************************************************
subroutine init (ux1,uy1,uz1,rho1,ep1,phi1,&
     gx1,gy1,gz1,rhos1,phis1,&
     hx1,hy1,hz1,rhoss1,phiss1,&
     pressure0)

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE MPI

  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,ep1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1,phis1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1,phiss1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: rho1,rhos1,rhoss1
  real(mytype) :: pressure0

  real(mytype) :: x,y,z,r,um,r1,r2,r3
  real(mytype) :: xspec,yspec,zspec
  integer :: k,j,i,fh,ierror,ii
  integer :: code
  integer (kind=MPI_OFFSET_KIND) :: disp

  ! LMN: set thermodynamic pressure
  pressure0 = 1._mytype

  if (iin.eq.1) then !generation of a random noise

    call system_clock(count=code)
    call random_seed(size = ii)
    call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /)) !

    call random_number(ux1)
    call random_number(uy1)
    call random_number(uz1)

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          ux1(i,j,k)=0._mytype !noise*ux1(i,j,k)
          uy1(i,j,k)=0._mytype !noise*uy1(i,j,k)
          uz1(i,j,k)=0._mytype !noise*uz1(i,j,k)
        enddo
      enddo
    enddo

    !modulation of the random noise
    do k=1,xsize(3)
      do j=1,xsize(2)
        if (istret.eq.0) then
          y=(j+xstart(2)-1-1)*dy-yly/2._mytype
        else
          y=yp(j+xstart(2)-1)-yly/2._mytype
        endif
        um=exp(-0.2_mytype*y*y)
        do i=1,xsize(1)
          x = (i + xstart(1) - 2) * dx - xlx / 2._mytype
          um = exp(-0.2_mytype * x**2)
          ux1(i,j,k)=um*ux1(i,j,k)
          uy1(i,j,k)=um*uy1(i,j,k)
          uz1(i,j,k)=um*uz1(i,j,k)
        enddo
      enddo
    enddo

    ! LMN: set density
    do k = 1, xsize(3)
      z = float(k + xstart(3) - 2) * dz
      do j = 1, xsize(2)
        y = float(j + xstart(2) - 2) * dy - yly / 2._mytype
        do i = 1, xsize(1)
          x = float(i + xstart(1) - 2) * dx

          rho1(i, j, k) = 0._mytype
        enddo
      enddo
    enddo
    
    if (iscalar==1) then
      do k=1,xsize(3)
        do j=1,xsize(2)
          do i=1,xsize(1)
            phi1(i,j,k)=0._mytype
            phis1(i,j,k)=phi1(i,j,k)
            phiss1(i,j,k)=phis1(i,j,k)
          enddo
        enddo
      enddo
    endif
  else if (iin.eq.2) then !read a correlated noise generated before
  else !set initial fields to zero
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          ux1(i,j,k)=0._mytype
          uy1(i,j,k)=0._mytype
          uz1(i,j,k)=0._mytype
        enddo
      enddo
    enddo
  endif

  !MEAN FLOW PROFILE
  call ecoule(ux1,uy1,uz1,rho1)
  if (ilmn.eq.0) then
    rho1(:,:,:) = 1._mytype
  endif
  !INIT FOR G AND U=MEAN FLOW + NOISE
  do k=1,xsize(3)
    do j=1,xsize(2)
      do i=1,xsize(1)
        ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
        uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
        uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
        gx1(i,j,k)=ux1(i,j,k)
        gy1(i,j,k)=uy1(i,j,k)
        gz1(i,j,k)=uz1(i,j,k)
        hx1(i,j,k)=gx1(i,j,k)
        hy1(i,j,k)=gy1(i,j,k)
        hz1(i,j,k)=gz1(i,j,k)

        rhos1(i,j,k) = rho1(i,j,k)
        rhoss1(i,j,k) = rhos1(i,j,k)
      enddo
    enddo
  enddo

  if (ivirt==2) then
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'epsilon.dat', &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_read_var(fh,disp,1,ep1) 
    call MPI_FILE_CLOSE(fh,ierror)
    if (nrank==0) print *,'read epsilon file done from init'
    print *,ep1
  endif

  return
end subroutine init

!********************************************************************
!
! 
!********************************************************************
subroutine divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
     td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,divu3,pp3,&
     nxmsize,nymsize,nzmsize,ph1,ph3,ph4,nlock)

  USE param
  USE IBM
  USE decomp_2d
  USE variables
  USE MPI

  implicit none

  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

  integer :: nxmsize,nymsize,nzmsize

  !X PENCILS NX NY NZ  -->NXM NY NZ
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1,ux1,uy1,uz1,ep1
  real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: td1,te1,tf1 
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: divu1
  !Y PENCILS NXM NY NZ  -->NXM NYM NZ
  real(mytype),dimension(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)) :: td2,te2,tf2,di2
  real(mytype),dimension(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)) :: ta2,tb2,tc2
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: divu2
  !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
  real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)) :: ta3,tb3,tc3,di3
  real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: td3,te3,tf3,pp3
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: divu3

  integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nlock
  integer :: code
  real(mytype) :: tmax,tmin,tmoy,tmax1,tmin1,tmoy1

  nvect1=xsize(1)*xsize(2)*xsize(3)
  nvect2=ysize(1)*ysize(2)*ysize(3)
  nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

  if (nlock==1) then
    if (ivirt.eq.0) ep1(:,:,:)=0.
    do ijk=1,nvect1
      ta1(ijk,1,1)=(1._mytype-ep1(ijk,1,1))*ux1(ijk,1,1)
      tb1(ijk,1,1)=(1._mytype-ep1(ijk,1,1))*uy1(ijk,1,1)
      tc1(ijk,1,1)=(1._mytype-ep1(ijk,1,1))*uz1(ijk,1,1)
    enddo
  else
    ta1(:,:,:)=ux1(:,:,:)
    tb1(:,:,:)=uy1(:,:,:)
    tc1(:,:,:)=uz1(:,:,:)
  endif

!!! CM call test_min_max('ta1  ','In divergence  ',ta1,size(ta1))
!!! CM call test_min_max('tb1  ','In divergence  ',tb1,size(tb1))
!!! CM call test_min_max('tc1  ','In divergence  ',tc1,size(tc1))

  !WORK X-PENCILS
  call decx6(td1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)
  if (nrhoscheme.eq.0) then
    !! Solving variable-coefficient Poisson equation
    !  Get divu to x pencils and interpolate to pressure points
    call transpose_z_to_y(divu3, divu2)
    call transpose_y_to_x(divu2, divu1)
    call inter6(te1,divu1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    td1(:,:,:) = td1(:,:,:) - te1(:,:,:)
  endif
  call inter6(te1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
  call inter6(tf1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

!!! CM call test_min_max('di1  ','In divergence  ',di1,size(di1))
!!! CM call test_min_max('td1  ','In divergence  ',td1,size(td1))
!!! CM call test_min_max('te1  ','In divergence  ',te1,size(te1))
!!! CM call test_min_max('tf1  ','In divergence  ',tf1,size(tf1))

  call transpose_x_to_y(td1,td2,ph4)!->NXM NY NZ
  call transpose_x_to_y(te1,te2,ph4)
  call transpose_x_to_y(tf1,tf2,ph4)

  !WORK Y-PENCILS
  call intery6(ta2,td2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
  call decy6(tb2,te2,di2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)
  call intery6(tc2,tf2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

  call transpose_y_to_z(ta2,ta3,ph3)!->NXM NYM NZ
  call transpose_y_to_z(tb2,tb3,ph3)
  call transpose_y_to_z(tc2,tc3,ph3)

  !WORK Z-PENCILS
  call interz6(td3,ta3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
       (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)    
  call interz6(te3,tb3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
       (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
  call decz6(tf3,tc3,di3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
       (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)

  do k=1,nzmsize
    do j=ph1%zst(2),ph1%zen(2)
      do i=ph1%zst(1),ph1%zen(1)
        pp3(i,j,k)=td3(i,j,k)+te3(i,j,k)+tf3(i,j,k)
      enddo
    enddo
  enddo

  if (nlock==2) then
    pp3(:,:,:)=pp3(:,:,:)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
  endif

  tmax=-1.e30_mytype
  tmin=+1.e30_mytype
  tmoy=0._mytype
  do k=1,nzmsize
    do j=ph1%zst(2),ph1%zen(2)
      do i=ph1%zst(1),ph1%zen(1)
        if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
        if (pp3(i,j,k).lt.tmin) tmin=pp3(i,j,k)
        tmoy=tmoy+abs(pp3(i,j,k))
      enddo
    enddo
  enddo
  tmoy=tmoy/nvect3

  call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(tmin,tmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)!

  if (nrank==0) then
    if (nlock==2) then
      print *,'DIV U final Max=',tmax1
      print *,'DIV U final Min=',tmin1
      print *,'DIV U final Moy=',tmoy1/real(nproc)
    else
      print *,'DIV U* Max=',tmax1
      print *,'DIV U* Min=',tmin1
      print *,'DIV U* Moy=',tmoy1/real(nproc)
    endif
  endif

  return
end subroutine divergence

!********************************************************************
!  SUBROUTINE: extrapol_rhotrans
! DESCRIPTION: Extrapolates the transient of density at time k+1.
!       INPUT: rho1,rhos1,rhoss1, current density (rho^{k+1}), old
!              divergence of momentum (-div(rho u)^{k+1}), and ?
!      OUTPUT: drhodt1, the predicted transient of continuity
!              equation at time k+1.
!        NOTE: All input and output in X-pencils.
!********************************************************************
SUBROUTINE extrapol_rhotrans(rho1, rhos1, rhoss1, rhos01, drhodt1)

  USE param
  USE decomp_2d
  USE variables

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rho1, rhos1, rhoss1, rhos01, drhodt1
  INTEGER :: subitr

  INTEGER :: ijk, nxyz

  nxyz = xsize(1) * xsize(2) * xsize(3)

  IF (nrhoscheme.EQ.0) THEN
    IF (nrank.EQ.0) THEN
      PRINT *, "nrhoscheme=0 corresponds to variable-coefficient Poisson equation"
      PRINT *, "Shoul not be extrapolating drhodt!!!"
      STOP
    ENDIF
  ELSE IF (nscheme.EQ.1) THEN
    !! AB2
    IF (itime.EQ.1.AND.ilit.EQ.0) THEN
      drhodt1(:,:,:) = drhodt1(:,:,:) + rho1(:,:,:)
    ELSE
      drhodt1(:,:,:) = 3._mytype * rho1(:,:,:) - 4._mytype * rhoss1(:,:,:) + rhos01(:,:,:)
      drhodt1(:,:,:) = 0.5_mytype * drhodt1(:,:,:)
    ENDIF

    drhodt1(:,:,:) = drhodt1(:,:,:) / dt
  ELSE IF (nscheme.EQ.2) THEN
    !! RK3

    IF (nrhoscheme.EQ.1) THEN
      !! Straightforward approximation:
      !    ddt rho^{k+1} approx -div(rho u)^k = -rho^k div(u^k) - u^k cdot grad(rho^k)
      !                                       = rhos1
      drhodt1(:,:,:) = rhos1(:,:,:)
    ELSE IF (nrhoscheme.EQ.2) THEN
      !! Alternative approximation:
      !    ddt rho^{k+1} approx (rho^{k+1} - rho^k) / (c_k dt)
      drhodt1(:,:,:) = (drhodt1(:,:,:) + rho1(:,:,:)) / gdt(itr)
    ELSE
      !! Golanski
      IF (itime.GT.1) THEN
        drhodt1(:,:,:) = rhoss1(:,:,:)
        DO ijk = 1, nxyz
          DO subitr = 1, itr
            !! TODO Check should it be gdt(itr) or gdt(subitr)?
            !
            ! Based on testing, it appears Golanski2005 made a typo,
            ! their expression would use gdt(itr) not gdt(subitr)
            drhodt1(ijk, 1, 1) = drhodt1(ijk, 1, 1) &
                 + (gdt(subitr) / dt) * (rhoss1(ijk, 1, 1) - rhos01(ijk, 1, 1)) 
          ENDDO ! End loop over subitr
        ENDDO ! End loop over ijk
      ELSE
        ! Need to use first order approximation for first
        ! full timestep

        drhodt1(:,:,:) = rhos1(:,:,:)
        ! drhodt1(:,:,:) = (drhodt1(:,:,:) + rho1(:,:,:)) / gdt(itr)
      ENDIF
    ENDIF
  ELSE
    IF (nrank.EQ.0) THEN
      PRINT *, "Extrapolating drhodt only implemented for AB2 and RK3 (nscheme = 0,1)"
      STOP
    ENDIF
  ENDIF

ENDSUBROUTINE extrapol_rhotrans

!********************************************************************
!  SUBROUTINE: divergence_mom
! DESCRIPTION: In LMN with the constant-coefficient poisson equation
!              we need to approximate the divergence of momentum at
!              the new timestep as:
!                div(rho u)^{k+1} = -ddt rho^{k+1}
!              where some appropriate approximation is used to
!              extrapolate ddt rho^{k+1}.
!       INPUT: drhodt1, the approximation of ddt rho^{k+1} in X
!              stencil.
!      OUTPUT: pp3, the RHS of pressure Poisson equation in Z
!              stencil. Part of this has already been computed so
!              this is an addition.
!********************************************************************
SUBROUTINE divergence_mom(drhodt1, pp3, di1, di2, di3, nxmsize, nymsize, nzmsize, ph1, ph3, ph4)

  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  TYPE(DECOMP_INFO) :: ph1, ph3, ph4
  INTEGER :: nxmsize, nymsize, nzmsize

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: di1
  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), ysize(2), ysize(3)) :: di2
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(1):ph1%zen(2), xsize(3)) :: di3

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1
  REAL(mytype), DIMENSION(nxmsize, xsize(2), xsize(3)) :: divmom1
  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), ysize(2), ysize(3)) :: drhodt2
  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), nymsize, ysize(3)) :: divmom2
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), zsize(3)) :: drhodt3
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: divmom3
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: pp3

  ! Interpolate in x
  CALL inter6(divmom1, -drhodt1, di1, sx, cifxp6, cisxp6, ciwxp6, xsize(1), nxmsize, xsize(2), &
       xsize(3), 1)

  ! Interpolate in y
  CALL transpose_x_to_y(-divmom1, drhodt2, ph4) !->NXM NY NZ
  CALL intery6(divmom2, -drhodt2, di2, sy, cifyp6, cisyp6, ciwyp6, (ph1%yen(1) - ph1%yst(1) + 1), &
       ysize(2), nymsize, ysize(3), 1)

  ! Interpolate in z
  CALL transpose_y_to_z(-divmom2, drhodt3, ph3) !->NXM NYM NZ
  CALL interz6(divmom3, -drhodt3, di3, sz, cifzp6, ciszp6, ciwzp6, (ph1%zen(1) - ph1%zst(1) + 1), &
       (ph1%zen(2) - ph1%zst(2) + 1), zsize(3), nzmsize, 1)

  ! Add new divergence of momentum to RHS of Poisson equation
  pp3(:,:,:) = pp3(:,:,:) - divmom3(:,:,:)

ENDSUBROUTINE divergence_mom

!********************************************************************
!
!
!********************************************************************
subroutine gradp(ta1,tb1,tc1,di1,td2,tf2,ta2,tb2,tc2,di2,&
     ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)

  USE param 
  USE decomp_2d
  USE variables

  implicit none

  TYPE(DECOMP_INFO) :: ph2,ph3
  integer :: i,j,k,ijk,nxmsize,nymsize,nzmsize,code
  integer, dimension(2) :: dims, dummy_coords
  logical, dimension(2) :: dummy_periods


  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
  !Z PENCILS NXM NYM NZM-->NXM NYM NZ
  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,tc3,di3 
  !Y PENCILS NXM NYM NZ -->NXM NY NZ
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2,tc2
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,td2,tf2,di2
  !X PENCILS NXM NY NZ  -->NX NY NZ
  real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: td1,te1,tf1 
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1



  !WORK Z-PENCILS

  call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
       (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
  call deciz6(tc3,pp3,di3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
       (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

  !WORK Y-PENCILS
  call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
  call transpose_z_to_y(tc3,tc2,ph3)

  call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
       (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
  call deciy6(td2,ta2,di2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
       (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
  call interiy6(tf2,tc2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
       (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

  !WORK X-PENCILS

  call transpose_y_to_x(tb2,td1,ph2) !nxm ny nz
  call transpose_y_to_x(td2,te1,ph2)
  call transpose_y_to_x(tf2,tf1,ph2)

  call deci6(ta1,td1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
       nxmsize,xsize(1),xsize(2),xsize(3),1)
  call interi6(tb1,te1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
       nxmsize,xsize(1),xsize(2),xsize(3),1)
  call interi6(tc1,tf1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
       nxmsize,xsize(1),xsize(2),xsize(3),1)


  !we are in X pencils:
  do k=1,xsize(3)
    do j=1,xsize(2)
      dpdyx1(j,k)=tb1(1,j,k)/gdt(itr)
      dpdzx1(j,k)=tc1(1,j,k)/gdt(itr)
      dpdyxn(j,k)=tb1(nx,j,k)/gdt(itr)
      dpdzxn(j,k)=tc1(nx,j,k)/gdt(itr)
    enddo
  enddo


  if (xstart(3)==1) then
    do j=1,xsize(2)
      do i=1,xsize(1)
        dpdxz1(i,j)=ta1(i,j,1)/gdt(itr)
        dpdyz1(i,j)=tb1(i,j,1)/gdt(itr)
      enddo
    enddo
  endif
  if (xend(3)==nz) then
    do j=1,xsize(2)
      do i=1,xsize(1)
        dpdxzn(i,j)=ta1(i,j,nz)/gdt(itr)
        dpdyzn(i,j)=tb1(i,j,nz)/gdt(itr)
      enddo
    enddo
  endif

  ! determine the processor grid in use
  call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
       dims, dummy_periods, dummy_coords, code)

  if (dims(1)==1) then
    do k=1,xsize(3)
      do i=1,xsize(1)
        dpdxy1(i,k)=ta1(i,1,k)/gdt(itr)
        dpdzy1(i,k)=tc1(i,1,k)/gdt(itr)
      enddo
    enddo
    do k=1,xsize(3)
      do i=1,xsize(1)
        dpdxyn(i,k)=ta1(i,xsize(2),k)/gdt(itr)
        dpdzyn(i,k)=tc1(i,xsize(2),k)/gdt(itr)
      enddo
    enddo
  else
    !find j=1 and j=ny
    if (xstart(2)==1) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxy1(i,k)=ta1(i,1,k)/gdt(itr)
          dpdzy1(i,k)=tc1(i,1,k)/gdt(itr)
        enddo
      enddo
    endif
    !      print *,nrank,xstart(2),ny-(nym/p_row)
    if (ny-(nym/dims(1))==xstart(2)) then
      do k=1,xsize(3)
        do i=1,xsize(1)
          dpdxyn(i,k)=ta1(i,xsize(2),k)/gdt(itr)
          dpdzyn(i,k)=tc1(i,xsize(2),k)/gdt(itr)
        enddo
      enddo
    endif

  endif


  return
end subroutine gradp

!********************************************************************
!
! 
!********************************************************************
subroutine corgp_IBM (ux,uy,uz,px,py,pz,nlock)

  USE param 
  USE decomp_2d
  USE variables

  implicit none

  integer :: ijk,nlock,nxyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz

  if (itime==1) then
    px(:,:,:)=0._mytype
    py(:,:,:)=0._mytype
    pz(:,:,:)=0._mytype
  endif

  nxyz=xsize(1)*xsize(2)*xsize(3)

  if (nlock.eq.1) then
    if (nz.gt.1) then
      do ijk=1,nxyz
        uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1) 
        uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1) 
        ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
    else
      do ijk=1,nxyz
        uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1) 
        ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
      enddo
    endif
  endif
  if (nlock.eq.2) then
    if (nz.gt.1) then
      do ijk=1,nxyz
        uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1) 
        uz(ijk,1,1)=pz(ijk,1,1)+uz(ijk,1,1) 
        ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
    else
      do ijk=1,nxyz
        uy(ijk,1,1)=py(ijk,1,1)+uy(ijk,1,1) 
        ux(ijk,1,1)=px(ijk,1,1)+ux(ijk,1,1)
      enddo
    endif
  endif

  return
end subroutine corgp_IBM

!*******************************************************************
!
!
!*******************************************************************
subroutine body(ux,uy,uz,ep1)

  USE param 
  USE decomp_2d
  USE variables
  USE IBM

  implicit none

  real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,ep1
  integer :: i,j,k
  real(mytype) :: xm,ym,r

  ep1=0._mytype
  do k=xstart(3),xend(3)
    do j=xstart(2),xend(2)
      do i=xstart(1),xend(1)
        xm=(i-1)*dx 
        ym=yp(j)!(j-1)*dy
        r=sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)) 
        if (r-ra >= 0._mytype) cycle
        ux(i,j,k)=0._mytype
        uy(i,j,k)=0._mytype 
        uz(i,j,k)=0._mytype
        ep1(i,j,k)=1._mytype 
      enddo
    enddo
  enddo


  return  
end subroutine body

!****************************************************************************
!
!****************************************************************************
subroutine pre_correc(ux,uy,uz,rho)

  USE decomp_2d
  USE variables
  USE param
  USE var
  USE MPI


  implicit none

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: rho
  integer :: i,j,k,code
  real(mytype) :: ut,ut1,utt,ut11
  integer, dimension(2) :: dims, dummy_coords
  logical, dimension(2) :: dummy_periods

  if (itime==1) then
    dpdyx1=0._mytype
    dpdzx1=0._mytype
    dpdyxn=0._mytype
    dpdzxn=0._mytype
  endif

  !we are in X pencils:
  do k=1,xsize(3)
    do j=1,xsize(2)
      dpdyx1(j,k)=dpdyx1(j,k)*gdt(itr)
      dpdzx1(j,k)=dpdzx1(j,k)*gdt(itr)
      dpdyxn(j,k)=dpdyxn(j,k)*gdt(itr)
      dpdzxn(j,k)=dpdzxn(j,k)*gdt(itr)
    enddo
  enddo

  if (xstart(3)==1) then
    do j=1,xsize(2)
      do i=1,xsize(1)
        dpdxz1(i,j)=dpdxz1(i,j)*gdt(itr)
        dpdyz1(i,j)=dpdyz1(i,j)*gdt(itr)
      enddo
    enddo
  endif
  if (xend(3)==nz) then
    do j=1,xsize(2)
      do i=1,xsize(1)
        dpdxzn(i,j)=dpdxzn(i,j)*gdt(itr)
        dpdyzn(i,j)=dpdyzn(i,j)*gdt(itr)
      enddo
    enddo
  endif

  if (xstart(2)==1) then
    do k=1,xsize(3)
      do i=1,xsize(1)
        dpdxy1(i,k)=dpdxy1(i,k)*gdt(itr)
        dpdzy1(i,k)=dpdzy1(i,k)*gdt(itr)
      enddo
    enddo
  endif
  if (xend(2)==ny) then
    do k=1,xsize(3)
      do i=1,xsize(1)
        dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
        dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
      enddo
    enddo
  endif

  ! Computatation of the flow rate Inflow/Outflow
  ! XXX ux, etc contain momentum, bxx contain velocity
  ! we are in X pencils:
  if (nclx==2) then
    ut1 = 0._mytype
    do k = 1, xsize(3)
      do j = 1, xsize(2)
        ut1 = ut1 + bxx1(j, k)
      enddo
    enddo
    ut1 = ut1 / (xsize(2) * xsize(3))
    call MPI_ALLREDUCE(ut1, ut11, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    ut11 = ut11 / nproc
    ut  =  0._mytype
    do k = 1, xsize(3)
      do j = 1, xsize(2)
        ut = ut + bxxn(j, k)
      enddo
    enddo
    ut = ut / (xsize(2) * xsize(3))
    call MPI_ALLREDUCE(ut, utt, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    utt = utt / nproc

    if (nrank.eq.0) then
      print *, 'FLOW RATE I/O [m^3 / s]', ut11, utt
    endif

    do k=1, xsize(3)
      do j=1, xsize(2)
        bxxn(j, k) = bxxn(j, k) - (utt - ut11)
      enddo
    enddo
  endif

  !********NCLX==2*************************************
  !****************************************************
  if (nclx.eq.2) then
    do k=1, xsize(3)
      do j=1, xsize(2)
        ux(1 , j, k)=rho(1 , j, k)*bxx1(j, k)
        uy(1 , j, k)=rho(1 , j, k)*bxy1(j, k)+dpdyx1(j, k)
        uz(1 , j, k)=rho(1 , j, k)*bxz1(j, k)+dpdzx1(j, k)
        ux(nx, j, k)=rho(nx, j, k)*bxxn(j, k)
        uy(nx, j, k)=rho(nx, j, k)*bxyn(j, k)+dpdyxn(j, k)
        uz(nx, j, k)=rho(nx, j, k)*bxzn(j, k)+dpdzxn(j, k)
      enddo
    enddo
  endif
  !****************************************************
  !********NCLY==2*************************************
  !****************************************************
  !WE ARE IN X PENCIL!!!!!!
  if (ncly==2) then
    if (itype.eq.2) then

      ! determine the processor grid in use
      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
           dims, dummy_periods, dummy_coords, code)

      if (dims(1)==1) then
        do k=1,xsize(3)
          do i=1,xsize(1)
            ux(i,1,k)=0._mytype+dpdxy1(i,k)
            uy(i,1,k)=0._mytype
            uz(i,1,k)=0._mytype+dpdzy1(i,k)
          enddo
        enddo
        do k=1,xsize(3)
          do i=1,xsize(1)
            ux(i,xsize(2),k)=0._mytype+dpdxyn(i,k)
            uy(i,xsize(2),k)=0._mytype
            uz(i,xsize(2),k)=0._mytype+dpdzyn(i,k)
          enddo
        enddo
      else
        !find j=1 and j=ny
        if (xstart(2)==1) then
          do k=1,xsize(3)
            do i=1,xsize(1)
              ux(i,1,k)=0._mytype+dpdxy1(i,k)
              uy(i,1,k)=0._mytype
              uz(i,1,k)=0._mytype+dpdzy1(i,k)
            enddo
          enddo
        endif
        !      print *,nrank,xstart(2),ny-(nym/p_row)
        if (ny-(nym/dims(1))==xstart(2)) then
          do k=1,xsize(3)
            do i=1,xsize(1)
              ux(i,xsize(2),k)=0._mytype+dpdxyn(i,k)
              uy(i,xsize(2),k)=0._mytype
              uz(i,xsize(2),k)=0._mytype+dpdzyn(i,k)
            enddo
          enddo
        endif

      endif
    endif
  endif
  !****************************************************
  !********NCLZ==2*************************************
  !****************************************************
  !****************************************************

  !##################################################### 

  return
end subroutine pre_correc

!*******************************************************
!*******************************************************
SUBROUTINE pre_correc_tractionfree(ux1, uy1, uz1, rho1, ta1, pp1, di1,&
     ta2, pp2, di2,&
     ta3, pp3, di3,&
     nxmsize, nymsize, nzmsize, ph2, ph3)

  USE MPI
  USE decomp_2d
  USE param
  USE variables

  IMPLICIT NONE

  INTEGER :: nxmsize, nymsize, nzmsize
  TYPE(DECOMP_INFO) :: ph2, ph3

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rho1

  REAL(mytype), DIMENSION(ph3%zst(1):ph3%zen(1), ph3%zst(2):ph3%zen(2), nzmsize), INTENT(IN) :: pp3
  REAL(mytype), DIMENSION(ph3%zst(1):ph3%zen(1), ph3%zst(2):ph3%zen(2), zsize(3)) :: ta3, di3
  REAL(mytype), DIMENSION(ph3%yst(1):ph3%yen(1), nymsize, ysize(3)) :: ta2
  REAL(mytype), DIMENSION(ph3%yst(1):ph3%yen(1), ysize(2), ysize(3)) :: pp2, di2
  REAL(mytype), DIMENSION(nxmsize, xsize(2), xsize(3)) :: ta1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: pp1, di1

  INTEGER :: i, j, k, i_slipcollar
  REAL(mytype) :: x, y, z
  REAL(mytype) :: p0
  REAL(mytype) :: ue
  REAL(mytype) :: l_slipcollar

  !! Set length of slip collar.
  !  For x >= xlx - l_slipcollar, the boundary is treated
  !  as a free-slip plane
  l_slipcollar = 0.1_mytype * xlx
  DO i = 1, xsize(1)
    x = (i + xstart(1) - 2) * dx
    IF (x.LT.(xlx - l_slipcollar)) THEN
      i_slipcollar = i
    ELSE
      EXIT
    ENDIF
  ENDDO

  !!=====================================================
  ! First interpolate pressure to primary grid
  CALL interiz6(ta3, pp3, di3, sz, cifip6z, cisip6z, ciwip6z, cifz6, cisz6, ciwz6,&
       (ph3%zen(1) - ph3%zst(1) + 1), (ph3%zen(2) - ph3%zst(2) + 1), nzmsize, zsize(3), 1)
  CALL transpose_z_to_y(ta3, ta2, ph3)
  CALL interiy6(pp2, ta2, di2, sy, cifip6y, cisip6y, ciwip6y, cify6, cisy6, ciwy6,&
       (ph3%yen(1) - ph3%yst(1) + 1), nymsize, ysize(2), ysize(3), 1)
  CALL transpose_y_to_x(pp2, ta1, ph2)
  CALL interi6(pp1, ta1, di1, sx, cifip6, cisip6, ciwip6, cifx6, cisx6, ciwx6,&
       nxmsize, xsize(1), xsize(2), xsize(3), 1)
  ! The pressure field on the main mesh is in pp1

  !!=====================================================
  ! Evaluate reference pressure (using x=xlx as reference
  ! plane)
  !!=====================================================
  ! Use Bernoulli to evaluate velocity crossing boundary
  IF (ncly.EQ.2) THEN
    DO k = 1,xsize(3)
      DO i = 1, i_slipcollar
        x = (i + xstart(1) - 2) * dx

        ue = p0
        ! Add gravity

        ue = ue - pp1(i, 1, k)
        ue = ue * (2._mytype / rho1(i, 1, k))
        ue = MAX(ue, 0._mytype)
        ue = SQRT(ue)
        uy1(i, 1, k) = ue

        ue = p0
        ! Add gravity

        ue = ue - pp1(i, xsize(2), k)
        ue = ue * (2._mytype / rho1(i, xsize(2), k))
        ue = MAX(ue, 0._mytype)
        ue = SQRT(ue)
        uy1(i, xsize(2), k) = -ue
      ENDDO
    ENDDO
  ENDIF
  IF (nclz.EQ.2) THEN
    DO j = 1, xsize(2)
      DO i = 1, i_slipcollar
        x = (i + xstart(1) - 2) * dx

        ue = p0
        ! Add gravity

        ue = ue - pp1(i, j, 1)
        ue = ue * (2._mytype / rho1(i, j, 1))
        ue = MAX(ue, 0._mytype)
        ue = SQRT(ue)
        uz1(i, j, 1) = ue

        ue = p0
        ! Add gravity

        ue = ue - pp1(i, j, xsize(3))
        ue = ue * (2._mytype / rho1(i, j, xsize(3)))
        ue = MAX(ue, 0._mytype)
        ue = SQRT(ue)
        uz1(i, j, xsize(3)) = -ue
      ENDDO
    ENDDO
  ENDIF
ENDSUBROUTINE pre_correc_tractionfree

