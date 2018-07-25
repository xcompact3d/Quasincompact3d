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

  REAL(mytype) :: udenslim, ldenslim
  INTEGER :: ijk, nxyz

  nxyz = xsize(1) * xsize(2) * xsize(3)

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

  !! Limiting
  CALL test_density_min_max(rho1)
  udenslim = MAX(dens1, dens2)
  ldenslim = MIN(dens1, dens2)
  DO ijk = 1, nxyz
    rho1(ijk, 1, 1) = MAX(rho1(ijk, 1, 1), ldenslim)
    rho1(ijk, 1, 1) = MIN(rho1(ijk, 1, 1), udenslim)
  ENDDO

ENDSUBROUTINE inttdensity

!********************************************************************
!********************************************************************
SUBROUTINE eval_densitycoeffs(rho1, temperature1, ta1, rhos1, rhoss1, rhos01, drhodt1)

  USE param
  USE variables
  USE decomp_2d

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: rho1, temperature1, ta1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rhos1, rhoss1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rhos01, drhodt1
  
  INTEGER :: ijk, nxyz

  nxyz = xsize(1) * xsize(2) * xsize(3)

  IF ((nscheme.EQ.1).OR.(nscheme.EQ.2)) THEN
    !! AB2 or RK3

    ! First store -rho1 in drhodt1 incase we use simple extrapolation
    drhodt1(:,:,:) = -rho1(:,:,:)

    IF (nscheme.EQ.1) THEN
      !! AB2
      rhos01(:,:,:) = rhoss1(:,:,:)
      rhoss1(:,:,:) = rho1(:,:,:)
    ELSE IF (itr.EQ.1) THEN
      !! RK3, first iteration
      rhos01(:,:,:) = rhoss1(:,:,:)
      rhoss1(:,:,:) = -(temperature1(:,:,:) / rho1(:,:,:)) * ta1(:,:,:)
    ENDIF
  ELSE IF (nscheme.EQ.3) THEN
    !! RK4
    !! XXX Not implemented!
    IF (nrank.EQ.0) THEN
      PRINT  *, 'LMN: RK4 not ready'
      STOP
    ENDIF
  ELSE
    !! AB3
    !! XXX Not implemented
    IF (nrank.EQ.0) THEN
      PRINT  *, 'LMN: AB3 not ready'
      STOP
    ENDIF
  ENDIF

  !! Update old stage
  rhos1(:,:,:) = -(temperature1(:,:,:) / rho1(:,:,:)) * ta1(:,:,:)
   
ENDSUBROUTINE eval_densitycoeffs

!********************************************************************
!********************************************************************
SUBROUTINE intttemperature(temperature1, temperatures1, temperaturess1, tg1)

  USE param
  USE variables
  USE decomp_2d

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: tg1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: temperature1, temperatures1, temperaturess1
  INTEGER :: ijk, nxyz

  nxyz = xsize(1) * xsize(2) * xsize(3)

  IF ((nscheme.EQ.1).OR.(nscheme.EQ.2)) THEN
    !! AB2 or RK3

    IF (((nscheme.EQ.1).AND.(itime.EQ.1).AND.(ilit.EQ.0)).OR.&
         ((nscheme.EQ.2).AND.(itr.EQ.1))) THEN
      temperature1(:,:,:) = temperature1(:,:,:) + gdt(itr) * tg1(:,:,:)
    ELSE
      temperature1(:,:,:) = temperature1(:,:,:) + adt(itr) * tg1(:,:,:) &
           + bdt(itr) * temperatures1(:,:,:)
    ENDIF
  ELSE IF (nscheme.EQ.3) THEN
    !! RK3
    IF (nrank.EQ.0) THEN
      PRINT *, "LMN: RK4 not ready!"
      STOP
    ENDIF
  ELSE
    !! AB3
    IF ((itime.EQ.1).AND.(ilit.EQ.0)) THEN
      IF (nrank.EQ.0) THEN
        PRINT  *, 'start with Euler', itime
      ENDIF
      temperature1(:,:,:) = temperature1(:,:,:) + dt * tg1(:,:,:)
    ELSE
      IF  ((itime.EQ.2).AND.(ilit.EQ.0)) THEN
        IF (nrank.EQ.0) THEN
          PRINT *, 'then with AB2', itime
        ENDIF
        temperature1(:,:,:) = temperature1(:,:,:) - 0.5_mytype * dt &
             * (temperatures1(:,:,:) - 3._mytype * tg1(:,:,:))
      ELSE
        temperature1(:,:,:) = temperature1(:,:,:) + adt(itr) * tg1(:,:,:) &
             + bdt(itr) * temperatures1(:,:,:) + cdt(itr) * temperaturess1(:,:,:)
      ENDIF

      !! Update oldold stage
      temperaturess1(:,:,:) = temperatures1(:,:,:)
    ENDIF
  ENDIF

  !! Update old stage
  temperatures1(:,:,:) = tg1(:,:,:)

  !! Limiting
  CALL test_temperature_min_max(temperature1)
  DO ijk = 1, nxyz
    temperature1(ijk, 1, 1) = MAX(temperature1(ijk, 1, 1), 1._mytype)
    temperature1(ijk, 1, 1) = MIN(temperature1(ijk, 1, 1), 1._mytype)
  ENDDO
  
ENDSUBROUTINE intttemperature

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
    if (ivarcoeff.eq.0) then
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
subroutine inflow (ux, uy, uz, rho, temperature, massfrac, phi)

  USE param
  USE IBM
  USE variables
  USE decomp_2d

  implicit none

  integer :: k, j, n, nmodes
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz, rho, temperature, phi, massfrac
  real(mytype) :: r1, r2, r3, y, z, um, theta, freq, St, mf, s

  call ecoule(ux, uy, uz, rho, temperature, massfrac)

  call random_number(bxo)
  call random_number(byo)
  call random_number(bzo)

  nmodes = 6

  St = 0.3
  freq = St * u1 / (1._mytype)

  if (t.LT.1._mytype) then
     s = SIN(t * (PI / 2._mytype))
  else
     s = 1._mytype
  endif

  if (iin.eq.1) then  
    do k = 1, xsize(3)
      z = (k + xstart(3) - 2) * dz - zlz / 2._mytype
      do j = 1, xsize(2)
        y = (j + xstart(2) - 2) * dy - yly / 2._mytype
        r1 = SQRT(y**2 + z**2)

        if (r1.lt.0.5_mytype) then
          mf = rho(1, j, k) * bxx1(j, k)

          IF (z.GT.0._mytype) THEN
            IF (y.GT.0._mytype) THEN
              theta = ATAN(y / z)
            ELSE IF (y.LT.0._mytype) THEN
              theta = 2._mytype * PI - ATAN(-y / z)
            ELSE
              theta = 0._mytype
            ENDIF
          ELSE IF (z.LT.0._mytype) THEN
            IF (y.GT.0._mytype) THEN
              theta = PI - ATAN(y / (-z))
            ELSE IF (y.LT.0._mytype) THEN
              theta = PI + ATAN(y / z)
            ELSE
              theta = PI
            ENDIF
          ELSE
            IF (y.GT.0._mytype) THEN
              theta = 0.5_mytype * PI
            ELSE IF (y.LT.0._mytype) THEN
              theta = 1.5_mytype * PI
            ELSE
              theta = 0._mytype
            ENDIF
          ENDIF

          !! Additional forcing
          um = 0._mytype
          DO n = 1, nmodes
            um = um + SIN(2._mytype * PI * (t * freq) / (DBLE(n)) + theta)
          ENDDO
          um = 0.2_mytype * um / DBLE(nmodes)
          IF (t.LT.1._mytype) THEN
            um = um * SIN(t * (0.5_mytype * PI))
          ENDIF
          um = s * um
          bxx1(j, k) = (1._mytype + um) * bxx1(j, k)

          bxx1(j, k) = bxx1(j, k) + noise1 * (1._mytype - 2._mytype * bxo(j, k))
          bxy1(j, k) = bxy1(j, k) + noise1 * (1._mytype - 2._mytype * byo(j, k))
          bxz1(j, k) = bxz1(j, k) + noise1 * (1._mytype - 2._mytype * bzo(j, k))
          bxx1(j, k) = MAX(bxx1(j, k), u2) ! Prevent backflow
          
          ! if ((mf * bxx1(j, k)).gt.0._mytype) then
          !   rho(1, j, k) = mf / bxx1(j, k)
          ! endif
        endif ! End within jet
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
subroutine outflow (ux, uy, uz, rho, temperature, massfrac, phi)

  USE param
  USE variables
  USE decomp_2d
  USE MPI

  implicit none

  integer :: j, k, i,  code
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz, rho, temperature, massfrac, phi
  real(mytype) :: udx, udy, udz, uddx, uddy, uddz, uxmax, &
       uxmin, vphase, coef, uxmax1, uxmin1, volflux, volflux_out
  real(mytype), dimension(xsize(2), xsize(3)) :: cx
  real(mytype) :: Ay

  real(mytype) :: y, z, yc, zc
  real(mytype) :: r2

  real(mytype) :: ucf, g_umax, gauss
  real(mytype) :: g_rext, g_rext2
  
  !! Compute 'convective velocity' at outlet
  udx = 1._mytype / dx
  udy = 1._mytype / dy
  udz = 1._mytype / dz
  uddx = 0.5_mytype / dx
  uddy = 0.5_mytype / dy
  uddz = 0.5_mytype / dz

  ! ! If inlet velocity specified in terms of u1 and u2
  ! cx(:,:) = 0.5_mytype * (u1 + u2) * gdt(itr) * udx

  ! uxmax = -1609._mytype
  ! uxmin = 1609._mytype
  ! do k = 1, xsize(3)
  !   do j = 1, xsize(2)
  !     if (ux(nx - 1, j, k).gt.uxmax) uxmax = ux(nx - 1, j, k)
  !     if (ux(nx - 1, j, k).lt.uxmin) uxmin = ux(nx - 1, j, k)
  !   enddo
  ! enddo
  ! call MPI_ALLREDUCE(uxmax, uxmax1, 1, real_type, MPI_MAX, MPI_COMM_WORLD, code)
  ! call MPI_ALLREDUCE(uxmin, uxmin1, 1, real_type, MPI_MIN, MPI_COMM_WORLD, code)
  ! vphase = 0.5_mytype * (uxmax1 + uxmin1)
  ! cx(:,:) = vphase * gdt(itr) * udx

  ! ! Compute mean velocity (inlet)
  ! volflux = 0._mytype
  ! do k = 1, xsize(3)
  !   do j = 1, xsize(2) - 1
  !     if (istret.eq.0) then
  !       Ay = yly / (ny - 1)
  !     else
  !       Ay = (yp(j + 1) - yp(j))
  !     endif
  !     volflux = volflux + 0.5_mytype * (ux(1, j, k) + ux(1, j + 1, k)) * Ay * dz
  !   enddo
  ! enddo
  ! PRINT *, volflux
  ! call MPI_ALLREDUCE(MPI_IN_PLACE, volflux, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
  ! volflux = volflux / (yly * zlz)
  ! cx(:,:) = volflux * gdt(itr) * udx

  !! Gaussian outflow (to balance inlet mass flux)
  g_rext = 1.5_mytype / 2.14_mytype ! 2.14 is the magic number of (r/R) giving e^(-r**2 / R**2) = 0.01
  g_rext2 = g_rext**2
  g_umax = 1._mytype ! We will calculate this later to balance mass flux

  yc = 0.5_mytype * yly
  zc = 0.5_mytype * zlz

  if (ilmn.ne.0) then
     volflux = outflux ! The required outflux, computed by compute_outflux_lmn
  else
     volflux = u2 * (yly * zlz)
     volflux = volflux + (u1 - u2) * (PI * (0.5_mytype**2))
  endif
  ! volflux_out = 0._mytype
  ! ucf = 0.1_mytype * u1 * (PI * 0.5_mytype**2) / (yly * zlz)
  ! DO k = 1, nz
  !   z = DBLE(k - 1) * dz - zc
  !   DO j = 1, ny
  !     y = DBLE(j - 1) * dy - yc

  !     r2 = y**2 + z**2
  !     gauss = g_umax * EXP(-r2 / g_rext2) + ucf

  !     volflux_out = volflux_out + gauss * dy * dz
  !   ENDDO
  ! ENDDO
  ! bxxn_scale = volflux / volflux_out

  ! DO k = 1, xsize(3)
  !   z = DBLE(k + xstart(3) - 2) * dz - zc
  !   DO j = 1, xsize(2)
  !     y = DBLE(j + xstart(2) - 2) * dy - yc

  !     r2 = y**2 + z**2
  !     gauss = bxxn_scale * (g_umax * EXP(-r2 / g_rext2) + ucf)
      
  !     cx(j, k) = gauss * (gdt(itr) * udx)
  !   ENDDO
  ! ENDDO

  !! Set average outflux
  volflux = volflux / (xlx * yly)
  DO k = 1, xsize(3)
    DO j = 1, xsize(2)
      cx(j, k) = volflux * (gdt(itr) * udx)
    ENDDO
  ENDDO

  ! volflux = u1 * (PI * 0.5_mytype**2) / (yly * zlz)
  ! cx(:,:) = volflux * gdt(itr) * udx

  ! !! Volume correction
  ! volflux = 0._mytype
  ! DO k = 1, xsize(3)
  !   DO j = 1, xsize(2) - 1
  !     IF (istret.EQ.0) THEN
  !       Ay = yly / (ny - 1)
  !     ELSE
  !       Ay = (yp(j + 1) - yp(j))
  !     ENDIF
  !     volflux = volflux + 0.5_mytype * (bxx1(j, k) + bxx1(j + 1, k)) * Ay * dz * gdt(itr) * udx
  !     volflux = volflux - 0.5_mytype * (cx(j, k) + cx(j + 1, k)) * Ay * dz
  !   ENDDO
  ! ENDDO

  ! IF (ncly.EQ.2) THEN
  !   DO k = 1, xsize(3)
  !     DO i = 1, xsize(1)
  !       volflux = volflux + (byy1(i, k) - byyn(i, k)) * dx * dz * gdt(itr) * udx
  !     ENDDO
  !   ENDDO
  ! ENDIF

  ! IF (nclz.EQ.2) THEN
  !   DO j = 1, xsize(2) - 1
  !     IF (istret.EQ.0) THEN
  !       Ay = yly / (ny - 1)
  !     ELSE
  !       Ay = (yp(j + 1) - yp(j))
  !     ENDIF
  !     DO i = 1, xsize(1)
  !       volflux = volflux + 0.5_mytype * ((bzz1(i, j) + bzz1(i, j + 1)) &
  !            - (bzzn(i, j) + bzzn(i, j + 1))) * dx * Ay * gdt(itr) * udx
  !     ENDDO
  !   ENDDO
  ! ENDIF
  
  ! CALL MPI_ALLREDUCE(MPI_IN_PLACE, volflux, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
  ! volflux = volflux / (yly * zlz)
  ! cx(:,:) = cx(:,:) + volflux
  ! u2 = 0._mytype
  
  if (itype.ne.9) then
    do k = 1, xsize(3)
      do j = 1, xsize(2)
        bxxn(j, k) = ux(nx, j, k) - cx(j, k) * (ux(nx, j, k) - ux(nx - 1, j, k))
        bxyn(j, k) = uy(nx, j, k) - cx(j, k) * (uy(nx, j, k) - uy(nx - 1, j, k))
        bxzn(j, k) = uz(nx, j, k) - cx(j, k) * (uz(nx, j, k) - uz(nx - 1, j, k))
        ! bxyn(j, k) = 0._mytype 
        ! bxzn(j, k) = 0._mytype

        massfrac(nx, j, k) = massfrac(nx, j, k) - cx(j, k) &
             * (massfrac(nx, j, k) - massfrac(nx - 1, j, k))
      enddo
    enddo

    if (isolvetemp.eq.0) then
      do k = 1, xsize(3)
        do j = 1, xsize(2)
          rho(nx, j, k) = rho(nx, j, k) - cx(j, k) * (rho(nx, j, k) - rho(nx - 1, j, k))
        enddo
      enddo
    else
      do k = 1, xsize(3)
        do j = 1, xsize(2)
          temperature(nx, j, k) = temperature(nx, j, k) - cx(j, k) &
               * (temperature(nx, j, k) - temperature(nx - 1, j, k))
        enddo
      enddo
    endif

    if (iscalar.eq.1) then
      do k = 1, xsize(3)
        do j = 1, xsize(2)
          phi(nx, j, k) = phi(nx, j, k) - cx(j, k) * (phi(nx, j, k) - phi(nx - 1, j, k))
        enddo
      enddo
    endif
  else
    print *, 'NOT READY'
    stop
  endif

  return
end subroutine outflow

SUBROUTINE compute_outflux_lmn(temperature1, gradtempx1, di1,&
     temperature2, gradtempy2, di2,&
     temperature3, gradtempz3, di3)

  USE MPI
  USE decomp_2d
  USE param
  USE variables
  
  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: temperature1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: gradtempx1, di1
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: temperature2, gradtempy2, di2
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: temperature3, gradtempz3, di3

  INTEGER :: i, j, k

  REAL(mytype) :: invpr

  INTEGER :: ierr
  REAL(mytype) :: outflux_local

  outflux_local = 0._mytype

  invpr = 1._mytype / pr

  IF (itype.NE.7) THEN
     IF (nclx.EQ.2) THEN
        CALL derx (gradtempx1,temperature1,di1,sx,ffxp,fsxp,fwxp,&
             xsize(1),xsize(2),xsize(3),1)
        DO k = 1, xsize(3)
           DO j = 1, xsize(2)
              outflux_local = outflux_local + bxx1(j, k) * (dy * dz)
              outflux_local = outflux_local + (xnu * invpr) &
                   * (gradtempx1(nx, j, k) - gradtempx1(1, j, k)) * (dy * dz)
           ENDDO
        ENDDO
     ENDIF

     IF (MAX(ncly, nclz).EQ.2) THEN
        CALL transpose_x_to_y(temperature1, temperature2)

        IF (ncly.EQ.2) THEN
           IF (xstart(2).EQ.1) THEN
              DO k = 1, xsize(3)
                 DO i = 1, xsize(1)
                    outflux_local = outflux_local + byy1(i, k) * (dx * dz)
                 ENDDO
              ENDDO
           ENDIF
           IF (xend(2).EQ.ny) THEN
              DO k = 1, xsize(3)
                 DO i = 1, xsize(1)
                    outflux_local = outflux_local - byyn(i, k) * (dx * dz)
                 ENDDO
              ENDDO
           ENDIF

           CALL dery (gradtempy2,temperature2,di2,sy,ffyp,fsyp,fwyp,ppy,&
                ysize(1),ysize(2),ysize(3),1)
           DO k = 1, ysize(3)
              DO i = 1, ysize(1)
                 outflux_local = outflux_local + (xnu * invpr) &
                      * (gradtempy2(i, ny, k) - gradtempy2(i, 1, k)) * (dx * dz)
              ENDDO
           ENDDO
        ENDIF

        IF (nclz.EQ.2) THEN
           IF (xstart(3).EQ.1) THEN
              DO j = 1, xsize(2)
                 DO i = 1, xsize(1)
                    outflux_local = outflux_local + bzz1(i, j) * (dx * dy)
                 ENDDO
              ENDDO
           ENDIF
           IF (xend(3).EQ.nz) THEN
              DO j = 1, xsize(2)
                 DO i = 1, xsize(1)
                    outflux_local = outflux_local - bzz1(i, j) * (dx * dy)
                 ENDDO
              ENDDO
           ENDIF

           CALL transpose_y_to_z(temperature2, temperature3)
           CALL derz (gradtempz3,temperature3,di3,sz,ffzp,fszp,fwzp,&
                zsize(1),zsize(2),zsize(3),1)
           DO j = 1, zsize(2)
              DO i = 1, zsize(1)
                 outflux_local = outflux_local + (xnu * invpr) &
                      * (gradtempz3(i, j, nz) - gradtempz3(i, j, 1)) * (dx * dy)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  
  CALL MPI_ALLREDUCE(outflux_local, outflux, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
  
ENDSUBROUTINE compute_outflux_lmn

SUBROUTINE set_velocity_entrainment_y(clx1, cly1, clz1)

  USE decomp_2d
  USE variables
  USE param
  
  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: clx1, cly1, clz1
  INTEGER :: i, j, k

  IF (xstart(2).EQ.1) THEN
     j = 1
     DO k = 1, xsize(3)
        DO i = 1, xsize(1)
           byx1(i, k) = clx1(i, j, k)
           byy1(i, k) = cly1(i, j, k)
           byz1(i, k) = clz1(i, j, k)
        ENDDO
     ENDDO
  ENDIF
  
  IF (xend(2).EQ.ny) THEN
     j = xsize(2)
     DO k = 1, xsize(3)
        DO i = 1, xsize(1)
           byxn(i, k) = clx1(i, j, k)
           byyn(i, k) = cly1(i, j, k)
           byzn(i, k) = clz1(i, j, k)
        ENDDO
     ENDDO
  ENDIF
  
ENDSUBROUTINE set_velocity_entrainment_y

SUBROUTINE set_density_entrainment_y(rho1, uy1)

  USE decomp_2d
  USE variables
  USE param
  
  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: uy1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rho1
  INTEGER :: i, j, k

  REAL(mytype) :: x
  REAL(mytype) :: cy

  REAL(mytype) :: l_fringe, xph_fringe
  INTEGER :: iph_fringe

  l_fringe = 0.1_mytype * xlx
  xph_fringe = xlx - l_fringe
  !! Find fringe
  DO i = 1, xsize(1)
     x = (i + xstart(1) - 2) * dx
     IF (x.GT.xph_fringe) THEN
        EXIT
     ELSE
        iph_fringe = i
     ENDIF
  ENDDO

  IF (ilmn.NE.0) THEN
    j = 1
    IF (xstart(2).EQ.1) THEN
      DO k = 1, xsize(3)
        DO i = 1, iph_fringe
          IF (uy1(i, j, k).GT.0._mytype) THEN
            !! INFLOW
            rho1(i, j, k) = dens2
          ELSE
            !! OUTFLOW
            cy = uy1(i, j, k) * gdt(itr) / dy
            rho1(i, j, k) = rho1(i, j, k) - cy * (rho1(i, j + 1, k) - rho1(i, j, k))
          ENDIF
       ENDDO
       DO i = iph_fringe + 1, xsize(1)
          rho1(i, j, k) = rho1(i, j + 1, k)
       ENDDO
      ENDDO
    ENDIF
    IF (xend(2).EQ.ny) THEN
      j = xsize(2)
      DO k = 1, xsize(3)
        DO i = 1, iph_fringe
          IF (uy1(i, j, k).LT.0._mytype) THEN
            !! INFLOW
            rho1(i, j, k) = dens2
          ELSE
            !! OUTFLOW
            cy = uy1(i, j, k) * gdt(itr) / dy
            rho1(i, j, k) = rho1(i, j, k) - cy * (rho1(i, j, k) - rho1(i, j - 1, k))
          ENDIF
       ENDDO
       DO i = iph_fringe + 1, xsize(1)
          rho1(i, j, k) = rho1(i, j - 1, k)
       ENDDO
      ENDDO
    ENDIF
  ENDIF
  
ENDSUBROUTINE set_density_entrainment_y

SUBROUTINE set_velocity_entrainment_z(clx1, cly1, clz1)

  USE decomp_2d
  USE variables
  USE param
  
  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: clx1, cly1, clz1
  INTEGER :: i, j, k

  IF (xstart(3).EQ.1) THEN
     k = 1
     DO j = 1, xsize(2)
        DO i = 1, xsize(1)
           bzx1(i, j) = clx1(i, j, k)
           bzy1(i, j) = cly1(i, j, k)
           bzz1(i, j) = clz1(i, j, k)
        ENDDO
     ENDDO
  ENDIF
  
  IF (xend(3).eq.nz) THEN
     k = xsize(3)
     DO j = 1, xsize(2)
        DO i = 1, xsize(1)
           bzxn(i, j) = clx1(i, j, k)
           bzyn(i, j) = cly1(i, j, k)
           bzzn(i, j) = clz1(i, j, k)
        ENDDO
     ENDDO
  ENDIF
  
ENDSUBROUTINE set_velocity_entrainment_z

SUBROUTINE set_density_entrainment_z(rho1, uz1)

  USE decomp_2d
  USE variables
  USE param
  
  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rho1
  INTEGER :: i, j, k

  REAL(mytype) :: x
  REAL(mytype) :: cz

  REAL(mytype) :: l_fringe, xph_fringe
  INTEGER :: iph_fringe

  l_fringe = 0.1_mytype * xlx
  xph_fringe = xlx - l_fringe
  !! Find fringe
  DO i = 1, xsize(1)
     x = (i + xstart(1) - 2) * dx
     IF (x.GT.xph_fringe) THEN
        EXIT
     ELSE
        iph_fringe = i
     ENDIF
  ENDDO

  IF (ilmn.NE.0) THEN
    IF (xstart(3).EQ.1) THEN
      k = 1
      DO j = 1, xsize(2)
         DO i = 1, iph_fringe
            IF (uz1(i, j, k).GT.0._mytype) THEN
               !! INFLOW
               rho1(i, j, k) = dens2
            ELSE
               !! OUTFLOW
               cz = uz1(i, j, k) * gdt(itr) / dz
               rho1(i, j, k) = rho1(i, j, k) - cz * (rho1(i, j, k + 1) - rho1(i, j, k))
            ENDIF
         ENDDO
         DO i = iph_fringe + 1, xsize(1)
            rho1(i, j, k) = rho1(i, j, k + 1)
         ENDDO
      ENDDO
    ENDIF
    IF (xend(3).EQ.nz) THEN
      k = xsize(3)
      DO j = 1, xsize(2)
         DO i = 1, iph_fringe
            IF (uz1(i, j, k).LT.0._mytype) THEN
               !! INFLOW
               rho1(i, j, k) = dens2
            ELSE
               !! OUTFLOW
               cz = uz1(i, j, k) * gdt(itr) / dz
               rho1(i, j, k) = rho1(i, j, k) - cz * (rho1(i, j, k) - rho1(i, j, k - 1))
            ENDIF
         ENDDO
         DO i = iph_fringe + 1, xsize(1)
            rho1(i, j, k) = rho1(i, j, k - 1)
         ENDDO
      ENDDO
    ENDIF
  ENDIF
  
ENDSUBROUTINE set_density_entrainment_z

SUBROUTINE set_density_bcs(rho1, ux1, uy1, uz1)

  USE decomp_2d
  USE param
  USE variables
  USE MPI

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rho1

  INTEGER :: i, j, k

  INTEGER :: ierr
  INTEGER, DIMENSION(2) :: dims, dummy_coords
  LOGICAL, DIMENSION(2) :: dummy_periods

  IF (nclx.EQ.2) THEN
  ENDIF

  IF ((ncly.EQ.2).OR.(nclz.EQ.2)) THEN
     ! determine the processor grid in use
     call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, ierr)
  ENDIF

  IF (ncly.EQ.2) THEN
     IF (dims(1).EQ.1) THEN
        DO k = 1, xsize(3)
           j = 1
           DO i = 1, xsize(1)
              rho1(i, j, k) = rho1(i, j + 1, k)
           ENDDO
           
           j = xsize(2)
           DO i = 1, xsize(1)
              rho1(i, j, k) = rho1(i, j - 1, k)
           ENDDO
        ENDDO
     ELSE
        IF (xstart(2).EQ.1) THEN
           j = 1
           DO k = 1, xsize(3)
              DO i = 1, xsize(1)
                 rho1(i, j, k) = rho1(i, j + 1, k)
              ENDDO
           ENDDO
        ENDIF

        IF (ny - (nym / dims(1)).EQ.xstart(2)) THEN
           j = xsize(2)
           DO k = 1, xsize(3)
              DO i = 1, xsize(1)
                 rho1(i, j, k) = rho1(i, j - 1, k)
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDIF

  IF (nclz.EQ.2) THEN
  ENDIF
  
ENDSUBROUTINE set_density_bcs

!**********************************************************************
!
!
!**********************************************************************
subroutine ecoule(ux1,uy1,uz1,rho1,temperature1,massfrac1)

  USE param
  USE IBM
  USE variables
  USE decomp_2d

  implicit none

  integer  :: i,j,k,jj1,jj2 
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,rho1,temperature1,massfrac1
  real(mytype) :: x,y,z,ym
  real(mytype) :: xspec,yspec,zspec
  real(mytype) :: r1,r2,r3,r
  real(mytype) :: uh,ud,um,xv,bruit1
  real(mytype) :: u_disturb, v_disturb, disturb_decay

  real(mytype) :: b2, D
  real(mytype) :: gauss, g_umax, ucf
  real(mytype) :: g_rext
  real(mytype) :: s

  bxx1=0._mytype;bxy1=0._mytype;bxz1=0._mytype
  byx1=0._mytype;byy1=0._mytype;byz1=0._mytype
  bzx1=0._mytype;bzy1=0._mytype;bzz1=0._mytype 

  !ITYPE=1 --> Constant flow field
  !ITYPE=2 --> Channel flow
  !ITYPE=3 --> Wake flow
  !ITYPE=4 --> Mixing layer with splitter plate
  !ITYPE=5 --> Jet flow
  !ITYPE=6 --> Taylor Green vortices
  !ITYPE=7 --> Cavity flow / lock-exchange
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
    u1 = SQRT(dens2 / dens1) / (SQRT(dens2 / dens1) + 1._mytype)
    u2 = -SQRT(dens1 / dens2) / (1._mytype + SQRT(dens1 / dens2))
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
          uy1(i, j, k) = 0._mytype
          uz1(i, j, k) = 0._mytype
          rho1(i, j, k) = 0._mytype

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
    b2 = 0.25_mytype * 10._mytype
    D = 1._mytype ! Jet diameter (used as length scale)

    if (itime.eq.0) then
      ucf = 0.1_mytype * (u1 * (PI * 0.5_mytype**2)) / (yly * zlz)
      g_umax = 0._mytype
      g_rext = 1.5_mytype / 2.14_mytype
      do k = 1,nz
        z = (k - 1) * dz - 0.5_mytype * zlz
        do j = 1,ny
          y = (j - 1) * dy - 0.5_mytype * yly
          r = SQRT(y**2 + z**2)
          g_umax = g_umax + EXP(-(r / g_rext)**2) * dy * dz
        enddo
      enddo
      g_umax = (u1 * (PI * 0.5_mytype**2) - ucf * (yly * zlz)) / g_umax
    endif

    if (t.lt.1._mytype) then
      s = SIN(t * (0.5_mytype * PI))
    else
      s = 1._mytype
    endif
    
    do k=1,xsize(3)
      z = (k + xstart(3) - 2) * dz - zlz / 2._mytype
      do j=1,xsize(2)
        if (istret.eq.0) then
          y=(j + xstart(2) - 2) * dy - yly / 2._mytype
        else
          y=yp(j)-yly/2._mytype
        endif
        r = SQRT(y**2 + z**2)
        r = r + 1.0e-32 ! Avoid division by zero

        ! Set the mean profile
        bxx1(j, k) = u2 - (u2 - u1) * 0.5_mytype * (1._mytype &
             - TANH(b2 * (2._mytype * r / D - D / (2._mytype * r))))
        bxx1(j, k) = MAX(bxx1(j, k), u2) ! Ensure a minimum velocity

        ! if (r.lt.0.5_mytype * D) then ! jet
        massfrac1(1, j, k) = 0._mytype - (0._mytype - 1._mytype) &
             * 0.5_mytype * (1._mytype - TANH(b2 * (2._mytype * r / D - D / (2._mytype * r))))
        ! else ! co-flow
        !   massfrac1(1, j, k) = 0._mytype
        ! endif

        if (isolvetemp.eq.0) then
           rho1(1, j, k) = dens2 - (dens2 - dens1) &
                * 0.5_mytype * (1._mytype - TANH(b2 * (2._mytype * r / D - D / (2._mytype * r))))
        else
           temperature1(1, j, k) = 1._mytype
        endif

        !! Smooth inflow in time, s = SIN(t * PI / 2) if t < 1, 1 otherwise
        ! bxx1(j, k) = s * bxx1(j, k) + (1._mytype - s) * u2
        rho1(1, j, k) = s * rho1(1, j, k) + (1._mytype - s) * dens2
        temperature1(1, j, k) = s * temperature1(1, j, k) + (1._mytype - s) * 1._mytype
        ! massfrac1(1, j, k) = s * massfrac1(1, j, k) + (1._mytype - s) * 0._mytype

        ! if (itime.eq.0) then
        !   !! Comment out this loop to make a "tube"
        !   do i = 1, xsize(1)
        !     ux1(i,j,k) = ux1(i,j,k) - bxx1(j,k)

        !     ! !! Uncomment this section to have smooth transition from inflow to outflow condition
        !     ! x = (i + xstart(1) - 2) * dx
        !     ! gauss = g_umax * EXP(-(r / g_rext)**2) + ucf
        !     ! ux1(i,j,k) = ux1(i,j,k) + (1._mytype - SIN((x / xlx) * (0.5_mytype * PI))) * bxx1(j,k)
        !     ! ux1(i,j,k) = ux1(i,j,k) + SIN((x / xlx) * (0.5_mytype * PI)) * gauss
        !   end do
        ! endif

        bxy1(j, k) = 0._mytype
        bxz1(j, k) = 0._mytype

        massfrac1(1, j, k) = (1._mytype - float(imulticomponent)) * 1._mytype &
             + float(imulticomponent) * massfrac1(1, j, k)
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
  else if (itype.eq.7) then
     if (itime.eq.0) then
        do k = 1, xsize(3)
           do j = 1, xsize(2)
              y = float((j + xstart(2) - 2)) * dy - 0.5_mytype * yly
              do i = 1, xsize(1)
                 x = float((i + xstart(1) - 2)) * dx - (14_mytype / 32._mytype) * xlx ! The weird 14/32 is from Birman
                 
                 ux1(i, j, k) = 0._mytype
                 uy1(i, j, k) = 0._mytype
                 uz1(i, j, k) = 0._mytype
                 
                 temperature1(i, j, k) = 1._mytype

                 !! Birman2005
                 rho1(i, j, k) = 0.5_mytype * ((dens2 / dens1 + 1._mytype) &
                      - (1._mytype - dens2 / dens1) * ERF(x * SQRT(sc / xnu)))
                 
                 massfrac1(i, j, k) = (dens1 * dens2 / (rho1(i, j, k) * temperature1(i, j, k)) &
                      - dens1) / (dens2 - dens1)
              enddo
           enddo
        enddo
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
    if (nrank.eq.0) then
      PRINT *, "itype=", itype, " not implemented!"
      STOP
    endif
  endif

  return
end subroutine ecoule

!********************************************************************
!
!
!********************************************************************
subroutine init (ux1,uy1,uz1,rho1,temperature1,massfrac1,ep1,phi1,&
     gx1,gy1,gz1,rhos1,temperatures1,massfracs1,phis1,&
     hx1,hy1,hz1,rhoss1,temperaturess1,massfracss1,phiss1,&
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
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: temperature1,temperatures1,temperaturess1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: massfrac1,massfracs1,massfracss1
  real(mytype) :: pressure0

  real(mytype) :: x,y,z,r,um,r1,r2,r3
  real(mytype) :: xspec,yspec,zspec
  integer :: k,j,i,fh,ierror,ii
  integer :: code
  integer (kind=MPI_OFFSET_KIND) :: disp
  real(mytype) :: rhor, rhol
  real(mytype) :: p_front

  p_front = 1._mytype
  
  ! LMN: set thermodynamic pressure
  pressure0 = 1._mytype

  if (iin.eq.1) then !generation of a random noise

    call system_clock(count=code)
    call random_seed(size = ii)
    call random_seed(put = code+63946*nrank*(/ (i - 1, i = 1, ii) /)) !

    call random_number(ux1)
    call random_number(uy1)
    call random_number(uz1)

    do k=1, xsize(3)
      do j=1, xsize(2)
        do i=1, xsize(1)
          ux1(i, j, k) = noise * (1._mytype - 2._mytype * ux1(i, j, k))
          uy1(i, j, k) = noise * (1._mytype - 2._mytype * uy1(i, j, k))
          uz1(i, j, k) = noise * (1._mytype - 2._mytype * uz1(i, j, k))
        enddo
      enddo
    enddo

    !modulation of the random noise
    rhol = dens2
    rhor = dens1
    do k=1,xsize(3)
       z = float(k + xstart(3) - 2) * dz - zlz / 2._mytype
       do j=1,xsize(2)
          if (istret.eq.0) then
             y=(j+xstart(2)-2)*dy-yly/2._mytype
          else
             y=yp(j+xstart(2)-1)-yly/2._mytype
          endif
          do i=1,xsize(1)
             x = (i + xstart(1) - 2) * dx
             x = x - p_front

             um = -((20._mytype / noise**3) / Fry) &
                  * ((rhol * p_front + rhor * (xlx - p_front)) * yly) &
                  / (rhol * (erf(0._mytype) - erf(5._mytype * sqrt(2._mytype) * p_front)) &
                  + (rhor * (erf(5._mytype * sqrt(2._mytype) * (xlx - p_front)) - erf(0._mytype))))
             um = noise * exp(-25._mytype * x**2)
             
             ux1(i,j,k)=um*ux1(i,j,k)
             uy1(i,j,k)=um*uy1(i,j,k)
             uz1(i,j,k)=um*uz1(i,j,k)
          enddo
       enddo
    enddo
    
    !! Set initial massfrac field:
    ! massfrac = 0 -> rho = dens1
    !          = 1 -> rho = dens2
    massfrac1(:,:,:) = 0._mytype
    ! massfrac1(:,:,:) = 1._mytype

    if (isolvetemp.eq.0) then
      ! LMN: set density
      do k = 1, xsize(3)
        z = float(k + xstart(3) - 2) * dz
        do j = 1, xsize(2)
          y = float(j + xstart(2) - 2) * dy - yly / 2._mytype
          do i = 1, xsize(1)
            x = float(i + xstart(1) - 2) * dx

            if (x.LT.p_front) then
               rho1(i, j, k) = dens2
            else
               rho1(i, j, k) = dens1
            endif
          enddo
        enddo
      enddo
    else
      temperature1(:,:,:) = 1._mytype
    endif
    
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
  call ecoule(ux1,uy1,uz1,rho1,temperature1,massfrac1)
  if (ilmn.eq.0) then
    rho1(:,:,:) = 1._mytype
    temperature1(:,:,:) = 1._mytype
  endif
  if (imulticomponent.eq.0) then
    massfrac1(:,:,:) = 0._mytype
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

        massfrac1(i,j,k)=massfrac1(1,j,k)
        massfracs1(i,j,k)=massfrac1(i,j,k)
        massfracss1(i,j,k)=massfracs1(i,j,k)
      enddo
    enddo
  enddo

  if (isolvetemp.eq.0) then
    call calctemp_eos(temperature1,rho1,massfrac1,pressure0,xsize)
  else
    call calcrho_eos(rho1,temperature1,massfrac1,pressure0,xsize)
  endif
  do k=1,xsize(3)
    do j=1,xsize(2)
      do i=1,xsize(1)
        rhos1(i,j,k) = rho1(i,j,k)
        rhoss1(i,j,k) = rhos1(i,j,k)
        temperatures1(i,j,k) = temperature1(i,j,k)
        temperaturess1(i,j,k) = temperatures1(i,j,k)
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
subroutine divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,drhodt1,&
     td2,te2,tf2,di2,ta2,tb2,tc2,&
     ta3,tb3,tc3,di3,td3,te3,tf3,divu3,pp3,&
     nxmsize,nymsize,nzmsize,ph1,ph3,ph4,nlock,quiet)

  USE param
  USE IBM
  USE decomp_2d
  USE variables
  USE MPI

  implicit none

  TYPE(DECOMP_INFO) :: ph1,ph3,ph4

  integer :: nxmsize,nymsize,nzmsize

  !X PENCILS NX NY NZ  -->NXM NY NZ
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1,ux1,uy1,uz1,ep1,drhodt1
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

  logical :: quiet

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

  !WORK X-PENCILS
  call decx6(td1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)
  if ((ilmn.ne.0) &               ! We are solving LMN equations AND
       .and.(((ivarcoeff.ne.0) &  ! (We are solving variable-coefficient Poisson equation
       .and.(nlock.ne.3)) &       !  but not making an initial guess based on const-coeff
       .or.(nlock.eq.2))) then    ! OR regardless of Poisson type we have solved it and would
                                  ! like to check how well we done!)
    !  Get divu to x pencils and interpolate to pressure points
    call transpose_z_to_y(divu3, divu2)
    call transpose_y_to_x(divu2, divu1)
    call inter6(te1,divu1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    td1(:,:,:) = td1(:,:,:) - te1(:,:,:)
  endif
  
  call inter6(te1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
  call inter6(tf1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

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

  if ((ilmn.ne.0) &                                 ! We are solving LMN equations AND
       .and.(((ivarcoeff.eq.0).and.(nlock.eq.1))) & ! (We are solving constant-coefficient Poisson equation
       .or.(nlock.eq.3)) then                       ! OR special case of predicting first step of var-coeff)
    !! Apply extrapolated divergence of momentum
    call divergence_mom(drhodt1,pp3,di1,di2,di3,nxmsize,nymsize,nzmsize,ph1,ph3,ph4)
  endif

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

  if (quiet.eqv..FALSE.) then
    if (nrank==0) then
      if (nlock==2) then
        print *,'ERR DIV U final Max=',tmax1
        print *,'ERR DIV U final Min=',tmin1
        print *,'ERR DIV U final Moy=',tmoy1/real(nproc)
      else
        if ((ilmn.eq.0).or.((ivarcoeff.ne.0).and.(nlock.ne.3))) then
          print *,'ERR DIV U* Max=',tmax1
          print *,'ERR DIV U* Min=',tmin1
          print *,'ERR DIV U* Moy=',tmoy1/real(nproc)
        else 
          print *,'ERR DIV (RHO U)* Max=',tmax1
          print *,'ERR DIV (RHO U)* Min=',tmin1
          print *,'ERR DIV (RHO U)* Moy=',tmoy1/real(nproc)
        endif
      endif
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

  IF (nscheme.EQ.1) THEN
     !! AB2
     IF (itime.EQ.1.AND.ilit.EQ.0) THEN
        drhodt1(:,:,:) = rho1(:,:,:) + drhodt1(:,:,:) ! drhodt1 stores -rho^k
     ELSE
        ! drhodt = (3 rho^{k+1} - 4 rho^k + rho^{k-1}) / (2 dt)
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
        drhodt1(:,:,:) = (rho1(:,:,:) + drhodt1(:,:,:)) / gdt(itr) ! drhodt1 stores -rho^k
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
           ! drhodt1(:,:,:) = (rho1(:,:,:) + drhodt1(:,:,:) / gdt(itr) ! drhodt1 stores -rho^k
        ENDIF
     ENDIF
  ELSE
     IF (nrank.EQ.0) THEN
        PRINT *, "Extrapolating drhodt only implemented for AB2 and RK3 (nscheme = 0,1)"
        STOP
     ENDIF
  ENDIF

ENDSUBROUTINE extrapol_rhotrans

SUBROUTINE birman_rhotrans_corr(rho1, drhodt1, ta1, tb1, di1, rho2, ta2, tb2, di2, rho3, ta3, di3)

  USE param
  USE variables
  USE decomp_2d

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: rho1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1, ta1, tb1, di1
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: rho2, ta2, tb2, di2
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rho3, ta3, di3

  REAL(mytype) :: invpe

  invpe = xnu / sc

  CALL transpose_x_to_y(rho1, rho2)
  CALL transpose_y_to_z(rho2, rho3)

  CALL derzz (ta3,rho3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)

  CALL transpose_z_to_y(ta3, tb2)

  CALL deryy (ta2,rho2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
  ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)

  CALL transpose_y_to_x(ta2, tb1)
  
  CALL derxx (ta1,rho1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
  ta1(:,:,:) = ta1(:,:,:) + tb1(:,:,:)

  drhodt1(:,:,:) = drhodt1(:,:,:) - invpe * ta1(:,:,:)
  
ENDSUBROUTINE birman_rhotrans_corr

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
!  SUBROUTINE: divergence_corr
! DESCRIPTION: In LMN with the variable-coefficient poisson equation
!              compute the correction term
!              1/rho nabla^2(p) - div( 1/rho grad(p) )
!********************************************************************
SUBROUTINE divergence_corr(rho1, px1, py1, pz1, ta1, tb1, tc1, td1, te1, tf1, di1, &
     py2, pz2, ta2, tb2, tc2, td2, di2, &
     pz3, pp3corr, ta3, tb3, tc3, di3, rho0p3, pp3, tg3, &
     nxmsize, nymsize, nzmsize, ph1, ph2, ph3, ph4, &
     divup3norm, poissiter, converged)

  USE MPI
  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param

  IMPLICIT NONE

  INTEGER :: ierr

  TYPE(DECOMP_INFO) :: ph1, ph2, ph3, ph4
  INTEGER, INTENT(IN) :: nxmsize, nymsize, nzmsize

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: px1, py1, pz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: rho1
  
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: td1, te1, tf1, di1
  REAL(mytype), DIMENSION(nxmsize , xsize(2), xsize(3)) :: ta1, tb1, tc1

  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), ysize(2), ysize(3)) :: tc2, td2, py2, di2
  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), nymsize , ysize(3)) :: ta2, tb2, pz2

  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), zsize(3)) :: tc3, pz3, di3
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: ta3, tb3, &
       pp3corr, rho0p3, tg3
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize), INTENT(IN) :: pp3

  REAL(mytype) :: resnorm, resnormlocal, pressnorm, pressnormlocal, deltadivunormlocal, deltadivunorm
  REAL(mytype), INTENT(IN) :: divup3norm
  INTEGER, INTENT(IN) :: poissiter
  LOGICAL, INTENT(OUT) :: converged

  LOGICAL :: file_exists

  REAL(mytype) :: rho0, rho0local
  REAL(mytype) :: rho1max, rho1maxlocal, rho1min, rho1minlocal

  REAL(mytype), DIMENSION(xszV(1), xszV(2), xszV(3)) :: uvisu
  character(len=20) :: filename

  INTEGER :: nxyzp3, ijk

  LOGICAL :: compute_rho0

  !! Correction = 1/rho0 nabla^2 p - div( 1/rho nabla p )
  resnorm = 0._mytype
  
  !! We should already know grad(p) in X-Pencils

  !!===================================================================================================
  ! Check convergence
  pressnormlocal = SUM((pp3corr(:,:,:) - tg3(:,:,:))**2)
  CALL MPI_ALLREDUCE(pressnormlocal, pressnorm, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
  pressnorm = SQRT(pressnorm)
  IF ((poissiter.NE.0).AND.(pressnorm.LE.tol)) THEN
    converged = .TRUE.  
  ENDIF
  tg3(:,:,:) = pp3corr(:,:,:) ! Store current pressure field
  
  !!===================================================================================================
  ! 1) Compute residual r = -div(1/rho grad p) + (div u^* - div u)

  !!---------------------------------------------------------------------------------------------------
  ! X PENCIL
  td1(:,:,:) = px1(:,:,:) / rho1(:,:,:)
  te1(:,:,:) = py1(:,:,:) / rho1(:,:,:)
  tf1(:,:,:) = pz1(:,:,:) / rho1(:,:,:)
  CALL decx6(ta1, td1, di1, sx, cfx6, csx6, cwx6, xsize(1), nxmsize, xsize(2), xsize(3), 0)
  CALL inter6(tb1, te1, di1, sx, cifxp6, cisxp6, ciwxp6, xsize(1), nxmsize, xsize(2), xsize(3), 1)
  CALL inter6(tc1, tf1, di1, sx, cifxp6, cisxp6, ciwxp6, xsize(1), nxmsize, xsize(2), xsize(3), 1)

  CALL transpose_x_to_y(ta1, tc2, ph4) !->NXM NY NZ
  CALL transpose_x_to_y(tb1, py2, ph4)
  CALL transpose_x_to_y(tc1, td2, ph4)
  
  !!---------------------------------------------------------------------------------------------------
  ! Y PENCIL
  CALL intery6(tb2, tc2, di2, sy, cifyp6, cisyp6, ciwyp6, (ph1%yen(1) - ph1%yst(1) + 1), ysize(2), &
       nymsize, ysize(3), 1)
  CALL decy6(ta2, py2, di2, sy, cfy6, csy6, cwy6, ppyi, (ph1%yen(1) - ph1%yst(1) + 1), ysize(2), &
       nymsize, ysize(3), 0)
  CALL intery6(pz2, td2, di2, sy, cifyp6, cisyp6, ciwyp6, (ph1%yen(1) - ph1%yst(1) + 1), ysize(2), &
       nymsize, ysize(3), 1)

  ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)

  CALL transpose_y_to_z(ta2, tc3, ph3) !->NXM NYM NZ
  CALL transpose_y_to_z(pz2, pz3, ph3)

  !!---------------------------------------------------------------------------------------------------
  ! Z PENCIL
  CALL interz6(tb3, tc3, di3, sz, cifzp6, ciszp6, ciwzp6, (ph1%zen(1) - ph1%zst(1) + 1), &
       (ph1%zen(2) - ph1%zst(2) + 1), zsize(3), nzmsize, 1)
  CALL decz6(ta3, pz3, di3, sz, cfz6, csz6, cwz6, (ph1%zen(1) - ph1%zst(1) + 1), &
       (ph1%zen(2) - ph1%zst(2) + 1), zsize(3), nzmsize, 0)
  
  !!----------------------------------------------------------------------------------------------------
  ! Subtract div(1/rho gradp) from div u^* - div u
  pp3corr(:,:,:) = pp3(:,:,:) - (ta3(:,:,:) + tb3(:,:,:))
  
  !!===================================================================================================
  ! Convergence test
  !  | -(div(u^*) - div(u^{k+1})) + dt div(1/rho grad(p^k)) | <= tol * | div(u^{k+1}) |
  resnormlocal = SUM(pp3corr(:,:,:)**2)
  CALL MPI_ALLREDUCE(resnormlocal, resnorm, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
  resnorm = SQRT(resnorm)
  IF (resnorm.LE.(MAX(tol * divup3norm, tol))) THEN
    converged = .TRUE.
  ENDIF

  deltadivunormlocal = SUM(pp3(:,:,:)**2)
  CALL MPI_ALLREDUCE(deltadivunormlocal, deltadivunorm, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
  deltadivunorm = SQRT(deltadivunorm)

  IF (nrank.EQ.0) THEN
    INQUIRE(FILE="VARPOISSON.log", EXIST=file_exists)
    IF (file_exists.EQV..TRUE.) THEN
      OPEN(10, FILE="VARPOISSON.log", STATUS="old", ACTION="write", POSITION="append")
    ELSE
      OPEN(10, FILE="VARPOISSON.log", STATUS="new", ACTION="write")
      WRITE(10, *) "| ITIME | ITR | POISSITER | POISSON RESIDUAL NORM | DIVU^{n+1} NORM | DIVU^* - DIVU^{n+1} NORM | DELTA P NORM |"
    ENDIF
    WRITE(10, *) itime, itr, poissiter, resnorm, divup3norm, deltadivunorm, pressnorm
    CLOSE(10)
  ENDIF

!   !!===================================================================================================
!   ! Visualise the residual
!   !WORK Z-PENCILS
!   call interiz6(tc3,pp3corr,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
!        (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!   !WORK Y-PENCILS
!   call transpose_z_to_y(tc3,ta2,ph3) !nxm nym nz
!   call interiy6(tc2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
!        (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!   !WORK X-PENCILS
!   call transpose_y_to_x(tc2,ta1,ph2) !nxm ny nz
!   call interi6(td1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
!        nxmsize,xsize(1),xsize(2),xsize(3),1)
!   !The pressure field on the main mesh is in td1
!   !PRESSURE
!   uvisu=0.
!   call fine_to_coarseV(1,td1,uvisu)
! 990 format('res',I1,"-",I1,"-",I3.3)
!   write(filename, 990) itime/imodulo, itr, poissiter
!   call decomp_2d_write_one(1,uvisu,filename,2)

  IF (converged.EQV..TRUE.) THEN
     !! Need to reset pp3corr to the pressure field
     pp3corr(:,:,:) = tg3(:,:,:)
  ELSE !! Not yet converged

    !!===================================================================================================
    ! 2) Scale the residual by density

    IF (poissiter.EQ.0) THEN
      compute_rho0 = .TRUE.
    ELSE
      compute_rho0 = .FALSE.
    ENDIF
    
    IF (compute_rho0.EQV..TRUE.) THEN
      !! Compute (and store) rho0
      
      IF (npoissscheme.EQ.0) THEN
        ! Simple minimum
        rho0local = MINVAL(rho1(:,:,:))
        CALL MPI_ALLREDUCE(rho0local, rho0, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
        rho0p3(:,:,:) = rho0
      ELSE
        ! Harmonic average
        td1(:,:,:) = 1._mytype / rho1(:,:,:)
        CALL inter6(ta1, td1, di1, sx, cifxp6, cisxp6, ciwxp6, xsize(1), nxmsize, xsize(2), xsize(3), 1)
        CALL transpose_x_to_y(ta1, tc2, ph4)
        CALL intery6(ta2, tc2, di2, sy, cifyp6, cisyp6, ciwyp6, (ph1%yen(1) - ph1%yst(1) + 1), ysize(2), &
             nymsize, ysize(3), 1)
        CALL transpose_y_to_z(ta2, tc3, ph3) !->NXM NYM NZ
        CALL interz6(tb3, tc3, di3, sz, cifzp6, ciszp6, ciwzp6, (ph1%zen(1) - ph1%zst(1) + 1), &
             (ph1%zen(2) - ph1%zst(2) + 1), zsize(3), nzmsize, 1)
        rho0p3(:,:,:) = 1._mytype / tb3(:,:,:)

        ! ! Perform some conditioning
        ! rho1maxlocal = MAXVAL(rho1(:,:,:))
        ! rho1minlocal = MINVAL(rho1(:,:,:))

        ! CALL MPI_ALLREDUCE(rho1maxlocal, rho1max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        ! CALL MPI_ALLREDUCE(rho1minlocal, rho1min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

        ! nxyzp3 = nzmsize * (ph1%zen(2) - ph1%zst(2)) * (ph1%zen(1) - ph1%zst(1))
        ! DO ijk = 1, nxyzp3
        !   rho0p3(ijk,1,1) = MAX(rho0p3(ijk,1,1), rho1min)
        !   rho0p3(ijk,1,1) = MIN(rho0p3(ijk,1,1), rho1max)
        ! ENDDO
      ENDIF
    ENDIF

    pp3corr(:,:,:) = rho0p3(:,:,:) * pp3corr(:,:,:)
    
    !!===================================================================================================
    ! 3) Compute div(grad p)

    ! td1(:,:,:) = (1._mytype - rho0 / rho1(:,:,:)) * px1(:,:,:)
    ! te1(:,:,:) = (1._mytype - rho0 / rho1(:,:,:)) * py1(:,:,:)
    ! tf1(:,:,:) = (1._mytype - rho0 / rho1(:,:,:)) * pz1(:,:,:)
    td1(:,:,:) = px1(:,:,:)
    te1(:,:,:) = py1(:,:,:)
    tf1(:,:,:) = pz1(:,:,:)
    
    !!---------------------------------------------------------------------------------------------------
    ! X PENCIL
    CALL decx6(ta1, td1, di1, sx, cfx6, csx6, cwx6, xsize(1), nxmsize, xsize(2), xsize(3), 0)
    CALL inter6(tb1, te1, di1, sx, cifxp6, cisxp6, ciwxp6, xsize(1), nxmsize, xsize(2), xsize(3), 1)
    CALL inter6(tc1, tf1, di1, sx, cifxp6, cisxp6, ciwxp6, xsize(1), nxmsize, xsize(2), xsize(3), 1)
    
    CALL transpose_x_to_y(ta1, tc2, ph4) !->NXM NY NZ
    CALL transpose_x_to_y(tb1, py2, ph4)
    CALL transpose_x_to_y(tc1, td2, ph4)
    
    !!---------------------------------------------------------------------------------------------------
    ! Y PENCIL
    CALL intery6(tb2, tc2, di2, sy, cifyp6, cisyp6, ciwyp6, (ph1%yen(1) - ph1%yst(1) + 1), ysize(2), &
         nymsize, ysize(3), 1)
    CALL decy6(ta2, py2, di2, sy, cfy6, csy6, cwy6, ppyi, (ph1%yen(1) - ph1%yst(1) + 1), ysize(2), &
         nymsize, ysize(3), 0)
    CALL intery6(pz2, td2, di2, sy, cifyp6, cisyp6, ciwyp6, (ph1%yen(1) - ph1%yst(1) + 1), ysize(2), &
         nymsize, ysize(3), 1)
    
    ta2(:,:,:) = ta2(:,:,:) + tb2(:,:,:)
    
    CALL transpose_y_to_z(ta2, tc3, ph3) !->NXM NYM NZ
    CALL transpose_y_to_z(pz2, pz3, ph3)
    
    !!---------------------------------------------------------------------------------------------------
    ! Z PENCIL
    CALL interz6(tb3, tc3, di3, sz, cifzp6, ciszp6, ciwzp6, (ph1%zen(1) - ph1%zst(1) + 1), &
         (ph1%zen(2) - ph1%zst(2) + 1), zsize(3), nzmsize, 1)
    CALL decz6(ta3, pz3, di3, sz, cfz6, csz6, cwz6, (ph1%zen(1) - ph1%zst(1) + 1), &
         (ph1%zen(2) - ph1%zst(2) + 1), zsize(3), nzmsize, 0)

    ta3(:,:,:) = ta3(:,:,:) + tb3(:,:,:)

    !!===================================================================================================
    ! 4) Compute final sum

    pp3corr(:,:,:) = pp3corr(:,:,:) + ta3(:,:,:)
  ENDIF
    
ENDSUBROUTINE divergence_corr

!********************************************************************
!  SUBROUTINE: approx_divergence_corr
! DESCRIPTION: In LMN with the variable-coefficient poisson equation
!              approximate the correction term
!              1/rho nabla^2(p) - div( 1/rho grad(p) )
!********************************************************************
SUBROUTINE approx_divergence_corr(ux1, uy1, uz1, rho1, ta1, tb1, tc1, td1, te1, tf1, ep1, di1, &
     rhos1, rhoss1, rhos01, drhodt1, &
     td2, te2, tf2, di2, ta2, tb2, tc2, &
     ta3, tb3, tc3, di3, td3, te3, tf3, tg3, pp3corr, divu3, &
     nxmsize, nymsize, nzmsize, ph1, ph3, ph4, &
     divup3norm)

  USE MPI
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  INTEGER :: ierr

  INTEGER :: nxmsize, nymsize, nzmsize
  TYPE(DECOMP_INFO) :: ph1, ph3, ph4

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: rho1, rhos1, rhoss1, rhos01
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: momx1, momy1, momz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ta1, tb1, tc1, ep1, di1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: drhodt1
  REAL(mytype), DIMENSION(nxmsize, xsize(2), xsize(3)) :: td1, te1, tf1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: divu1

  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), ysize(2), ysize(3)) :: td2, te2, tf2, di2
  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), nymsize, ysize(3)) :: ta2, tb2, tc2
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: divu2

  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), zsize(3)) :: ta3, tb3, tc3, &
       di3
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: td3, te3, tf3, &
       pp3corr, tg3
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: zero3
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3

  REAL(mytype) :: divup3normlocal
  REAL(mytype), INTENT(OUT) :: divup3norm
  
  !! Approximate 1/rho0 nabla^2 p following constant-coefficient Poisson equation
  !  i.e.
  !  nabla^2 p = 1/dt ( div(rho u)^* - div(rho u)^{k+1} )
  !            = 1/dt ( div(rho u)^* + drhodt^{k+1} )

  zero3(:,:,:) = 0._mytype
  momx1(:,:,:) = rho1(:,:,:) * ux1(:,:,:)
  momy1(:,:,:) = rho1(:,:,:) * uy1(:,:,:)
  momz1(:,:,:) = rho1(:,:,:) * uz1(:,:,:)
  CALL extrapol_rhotrans(rho1,rhos1,rhoss1,rhos01,drhodt1)
  CALL divergence (momx1, momy1, momz1, ep1, ta1, tb1, tc1, di1, td1, te1, tf1, drhodt1, &
       td2, te2, tf2, di2, ta2, tb2, tc2, &
       ta3, tb3, tc3, di3, td3, te3, tf3, divu3, pp3corr, &
       nxmsize, nymsize, nzmsize, ph1, ph3, ph4, 3, .TRUE.)
  
  ! !! Approximate div( 1/rho nabla p ) using variable-coefficient Poisson equation
  ! !  i.e.
  ! !  div( 1/rho nabla p ) = 1/Delta t ( div(u^*) - div(u^{k+1}) )
  ! !
  ! !  which is already known in pp3.
  
  ! pp3corr(:,:,:) = pp3corr(:,:,:) - rho0p3(:,:,:) * pp3(:,:,:)

  ! !! Add original RHS (pp3 = 1/dt ( div(u^*) - div(u^{k+1}) )) to correction
  ! pp3corr(:,:,:) = pp3corr(:,:,:) + rho0p3(:,:,:) * pp3(:,:,:)

  !! Compute | div(u) |

  ! Z->Y->X
  CALL transpose_z_to_y(divu3, divu2)
  CALL transpose_y_to_x(divu2, divu1)
  CALL inter6(td1,divu1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

  ! X->Y
  CALL transpose_x_to_y(td1, td2, ph4)
  CALL intery6(ta2,td2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

  ! Y->Z
  CALL transpose_y_to_z(ta2,ta3,ph3)
  CALL interz6(td3,ta3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
       (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)

  divup3normlocal = SUM(td3(:,:,:)**2)
  CALL MPI_ALLREDUCE(divup3normlocal, divup3norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  divup3norm = SQRT(divup3norm)

  !! Set/store old pressure
  IF (itime.EQ.1) THEN
    tg3(:,:,:) = 0._mytype
  ELSE
    tg3(:,:,:) = pp3corr(:,:,:)
  ENDIF
  
ENDSUBROUTINE approx_divergence_corr

SUBROUTINE calc_divup3norm(divup3norm, divu3, ta1, tb1, di1, ta2, tb2, tc2, di2, ta3, tb3, di3, &
  nxmsize, nymsize, nzmsize, ph1, ph3, ph4)

  USE MPI
  USE decomp_2d
  USE variables
  USE param

  IMPLICIT NONE

  INTEGER :: ierr
  TYPE(DECOMP_INFO) :: ph1, ph3, ph4
  INTEGER, INTENT(IN) :: nxmsize, nymsize, nzmsize

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ta1, di1
  REAL(mytype), DIMENSION(nxmsize, xsize(2), xsize(3)) :: tb1
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: ta2, di2
  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), ysize(2), ysize(3)) :: tb2
  REAL(mytype), DIMENSION(ph1%yst(1):ph1%yen(1), nymsize , ysize(3)) :: tc2
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(IN) :: divu3, di3
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), zsize(3)) :: ta3
  REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize) :: tb3
  REAL(mytype) :: divup3normlocal
  REAL(mytype), INTENT(OUT) :: divup3norm

  ! Z->Y->X
  CALL transpose_z_to_y(divu3, ta2)
  CALL transpose_y_to_x(ta2, ta1)
  CALL inter6(tb1,ta1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

  ! X->Y
  CALL transpose_x_to_y(tb1, tb2, ph4)
  CALL intery6(tc2,tb2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

  ! Y->Z
  CALL transpose_y_to_z(tc2,ta3,ph3)
  CALL interz6(tb3,ta3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
       (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)

  divup3normlocal = SUM(tb3(:,:,:)**2)
  CALL MPI_ALLREDUCE(divup3normlocal, divup3norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  divup3norm = SQRT(divup3norm)
  
ENDSUBROUTINE calc_divup3norm

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

  if (dims(2)==1) then
    k = 1
    do j=1,xsize(2)
      do i=1,xsize(1)
        dpdxz1(i,j)=ta1(i,j,k)/gdt(itr)
        dpdyz1(i,j)=tb1(i,j,k)/gdt(itr)
      enddo
    enddo
    
    k = xsize(3)
    do j=1,xsize(2)
      do i=1,xsize(1)
        dpdxzn(i,j)=ta1(i,j,k)/gdt(itr)
        dpdyzn(i,j)=tb1(i,j,k)/gdt(itr)
      enddo
    enddo
  else
    if (xstart(3)==1) then
      k = 1
      do j=1,xsize(2)
        do i=1,xsize(1)
          dpdxz1(i,j)=ta1(i,j,k)/gdt(itr)
          dpdyz1(i,j)=tb1(i,j,k)/gdt(itr)
        enddo
      enddo
    endif
    if (xend(3)==nz) then
      k = xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          dpdxzn(i,j)=ta1(i,j,k)/gdt(itr)
          dpdyzn(i,j)=tb1(i,j,k)/gdt(itr)
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

SUBROUTINE apply_grav(ta1, tb1, tc1, rho1)

  USE decomp_2d
  USE variables
  USE param
  USE MPI

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ta1, tb1, tc1, rho1
  REAL(mytype) :: invfrx, invfry, invfrz
  REAL(mytype) :: bcx, bcy, bcz

  INTEGER :: i, j, k, ijk, nxyz

  LOGICAL :: zerograv

  INTEGER :: code
  INTEGER, DIMENSION(2) :: dims, dummy_coords
  LOGICAL, DIMENSION(2) :: dummy_periods

  nxyz = xsize(1) * xsize(2) * xsize(3)
  zerograv = .TRUE.

  IF (nclx.EQ.0) THEN
     bcx = 0._mytype
  ELSE
     bcx = 1._mytype
  ENDIF
  IF (ncly.EQ.0) THEN
     bcy = 0._mytype
  ELSE
     bcy = 1._mytype
  ENDIF
  IF (nclz.EQ.0) THEN
     bcz = 0._mytype
  ELSE
     bcz = 1._mytype
  ENDIF

  IF ((frx.GT.0._mytype).OR.(frx.LT.0._mytype)) THEN
    invfrx = 1._mytype / frx
    zerograv = .FALSE.
  ELSE
    invfrx = 0._mytype
  ENDIF
  IF ((fry.GT.0._mytype).OR.(fry.LT.0._mytype)) THEN
    invfry = 1._mytype / fry
    zerograv = .FALSE.
  ELSE
    invfry = 0._mytype
  ENDIF
  IF ((frz.GT.0._mytype).OR.(frz.LT.0._mytype)) THEN
    invfrz = 1._mytype / frz
    zerograv = .FALSE.
  ELSE
    invfrz = 0._mytype
  ENDIF

  IF (zerograv.EQV..FALSE.) THEN
     DO ijk = 1, nxyz
        ta1(ijk, 1, 1) = ta1(ijk, 1, 1) + (rho1(ijk, 1, 1) - 1._mytype) * invfrx
        tb1(ijk, 1, 1) = tb1(ijk, 1, 1) + (rho1(ijk, 1, 1) - 1._mytype) * invfry
        tc1(ijk, 1, 1) = tc1(ijk, 1, 1) + (rho1(ijk, 1, 1) - 1._mytype) * invfrz
     ENDDO

     !! Apply BCs
     !  For periodic boundaries the heavier fluid should be able to 'fall' through boundary.
     
     i = 1
     DO k = 1, xsize(3)
        DO j = 1, xsize(2)
           ta1(i, j, k) = ta1(i, j, k) - (bcx * rho1(i, j, k) + (1._mytype - bcx) - 1._mytype) &
                * invfrx
        ENDDO
     ENDDO
     i = xsize(1)
     DO k = 1, xsize(3)
        DO j = 1, xsize(2)
           ta1(i, j, k) = ta1(i, j, k) - (bcx * rho1(i, j, k) + (1._mytype - bcx) - 1._mytype) &
                * invfrx
        ENDDO
     ENDDO
     
     CALL MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, code)

     IF (dims(1).EQ.1) THEN
        j = 1
        DO k = 1, xsize(3)
           DO i = 1, xsize(1)
              tb1(i, j, k) = tb1(i, j, k) - (bcy * rho1(i, j, k) + (1._mytype - bcy) - 1._mytype) &
                   * invfry
           ENDDO
        ENDDO
        j = xsize(2)
        DO k = 1, xsize(3)
           DO i = 1, xsize(1)
              tb1(i, j, k) = tb1(i, j, k) - (bcy * rho1(i, j, k) + (1._mytype - bcy) - 1._mytype) &
                   * invfry
           ENDDO
        ENDDO
     ELSE
        IF (xstart(2).EQ.1) THEN
           j = 1
           DO k = 1, xsize(3)
              DO i = 1, xsize(1)
                 tb1(i, j, k) = tb1(i, j, k) &
                      - (bcy * rho1(i, j, k) + (1._mytype - bcy) - 1._mytype) * invfry
              ENDDO
           ENDDO
        ENDIF
        IF (ny - (ny / dims(1)).EQ.xstart(2)) THEN
           j = xsize(2)
           DO k = 1, xsize(3)
              DO i = 1, xsize(1)
                 tb1(i, j, k) = tb1(i, j, k) &
                      - (bcy * rho1(i, j, k) + (1._mytype - bcy) - 1._mytype) * invfry
              ENDDO
           ENDDO
        ENDIF
     ENDIF

     IF (dims(2).EQ.1) THEN
        k = 1
        DO k = 1, xsize(2)
           DO i = 1, xsize(1)
              tc1(i, j, k) = tc1(i, j, k) - (bcz * rho1(i, j, k) + (1._mytype - bcz) - 1._mytype) &
                   * invfrz
           ENDDO
        ENDDO
        k = xsize(3)
        DO j = 1, xsize(2)
           DO i = 1, xsize(1)
              tc1(i, j, k) = tc1(i, j, k) - (bcz * rho1(i, j, k) + (1._mytype - bcz) - 1._mytype) &
                   * invfrz
           ENDDO
        ENDDO
     ELSE
        IF (xstart(3).EQ.1) THEN
           k = 1
           DO k = 1, xsize(2)
              DO i = 1, xsize(1)
                 tc1(i, j, k) = tc1(i, j, k) &
                      - (bcz * rho1(i, j, k) + (1._mytype - bcz) - 1._mytype) * invfrz
              ENDDO
           ENDDO
        ENDIF
        IF (nz - (nz / dims(2)).EQ.xstart(3)) THEN
           k = xsize(3)
           DO j = 1, xsize(2)
              DO i = 1, xsize(1)
                 tc1(i, j, k) = tc1(i, j, k) &
                      - (bcz * rho1(i, j, k) + (1._mytype - bcz) - 1._mytype) * invfrz
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDIF
  
ENDSUBROUTINE apply_grav

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

  real(mytype) :: Aout

  integer :: istart, iend, jstart, jend, kstart, kend

  real(mytype) :: uthreshold
  real(mytype) :: outflux_addedlocal, outflux_added, outflux_ratio

  if (itime==1) then
    dpdyx1=0._mytype
    dpdzx1=0._mytype
    dpdyxn=0._mytype
    dpdzxn=0._mytype

    dpdxy1 = 0._mytype
    dpdzy1 = 0._mytype
    dpdxyn = 0._mytype
    dpdzyn = 0._mytype

    dpdxy1 = 0._mytype
    dpdzy1 = 0._mytype
    dpdxyn = 0._mytype
    dpdzyn = 0._mytype
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

  ! Computatation of the flow rate Inflow/Outflow
  ! XXX ux, etc contain momentum, bxx contain velocity
  ! we are in X pencils:
  if (nclx==2) then

    uthreshold = 0.001_mytype * u1

    !! Compute inflow
    if (ilmn.eq.0) then
      ut1 = 0._mytype
      do k = 1, xsize(3)
        do j = 1, xsize(2)
          ut1 = ut1 + bxx1(j, k) * dy * dz
        enddo
      enddo
      
      ! if (ncly.eq.2) then
      !   do k = 1, xsize(3)
      !     do i = 1, xsize(1)
      !       ut1 = ut1 + byy1(i, k) * dx * dz
      !       ut1 = ut1 - byyn(i, k) * dx * dz
      !       Ain = Ain + 2._mytype * dx * dz
      !     enddo
      !   enddo
      ! endif
      
      ! if (nclz.eq.2) then
      !   do j = 1, xsize(2)
      !     do i = 1, xsize(1)
      !       ut1 = ut1 + bzz1(i, j) * dx * dy
      !       ut1 = ut1 - bzzn(i, j) * dx * dy
      !       Ain = Ain + 2._mytype * dx * dy
      !     enddo
      !   enddo
      ! endif
      
      call MPI_ALLREDUCE(ut1, ut11, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
      ! ut11 = ut11 / nproc
    else
      ut11 = outflux
    endif

    !! Compute outflow
    ut  =  0._mytype
    do k = 1, xsize(3)
      do j = 1, xsize(2)
        ut = ut + bxxn(j, k) * dy * dz
      enddo
    enddo
    call MPI_ALLREDUCE(ut, utt, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    Aout = yly * zlz
    ! utt = utt / nproc

    if (nrank.eq.0) then
      print *, 'FLOW RATE I/O [m^3 / s]', ut11, utt
    endif

    outflux_addedlocal = 0._mytype
    do k=1, xsize(3)
      do j=1, xsize(2)
        bxxn(j, k) = bxxn(j, k) - (utt - ut11) / Aout
        
        ! Check for/prevent backflow
        if (bxxn(j, k) < uthreshold) then
          outflux_addedlocal = outflux_addedlocal + (uthreshold - bxxn(j, k)) * (dy * dz)
          bxxn(j, k) = uthreshold
        endif
      enddo
    enddo

    ! Balance fluxes accounting for backflow
    call MPI_ALLREDUCE(outflux_addedlocal, outflux_added, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
    if (ut11.ne.0._mytype) then
       outflux_ratio = ut11 / (ut11 + outflux_added)
       bxxn(:,:) = outflux_ratio * bxxn(:,:)
    endif
  endif

  !********NCLX==2*************************************
  !****************************************************
  
  if (nclx.eq.2) then
     do k=1, xsize(3)
        do j=1, xsize(2)
           ux(1 , j, k)=rho(1 , j, k) * bxx1(j, k)
           uy(1 , j, k)=rho(1 , j, k) * bxy1(j, k)+dpdyx1(j, k)
           uz(1 , j, k)=rho(1 , j, k) * bxz1(j, k)+dpdzx1(j, k)
           ux(nx, j, k)=rho(nx, j, k) * bxxn(j, k)
           uy(nx, j, k)=rho(nx, j, k) * bxyn(j, k)+dpdyxn(j, k)
           uz(nx, j, k)=rho(nx, j, k) * bxzn(j, k)+dpdzxn(j, k)
        enddo
     enddo
  endif

  IF (nclx.EQ.2) THEN
     istart = 2
     iend = xsize(1) - 1
  ELSE
     istart = 1
     iend = xsize(1)
  ENDIF
  
  !****************************************************
  !********NCLY==2*************************************
  !****************************************************
  !WE ARE IN X PENCIL!!!!!!
  if (ncly==2) then
    if ((itype.eq.2).or.(itype.eq.7)) then

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
    else if(itype.eq.5) then
      ! determine the processor grid in use
      call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
           dims, dummy_periods, dummy_coords, code)
      if (dims(1).eq.1) then
        j = 1
        do k = 1, xsize(3)
          do i = istart, iend
            ux(i, j, k) = rho(i, j, k) * byx1(i, k) + dpdxy1(i, k)
            uy(i, j, k) = rho(i, j, k) * byy1(i, k)
            uz(i, j, k) = rho(i, j, k) * byz1(i, k) + dpdzy1(i, k)
          enddo
        enddo

        j = xsize(2)
        do k = 1, xsize(3)
          do i = istart, iend
            ux(i, j, k) = rho(i, j, k) * byxn(i, k) + dpdxyn(i, k)
            uy(i, j, k) = rho(i, j, k) * byyn(i, k)
            uz(i, j, k) = rho(i, j, k) * byzn(i, k) + dpdzyn(i, k)
          enddo
        enddo
      else
        if (xstart(2).eq.1) then
          j = 1
          do k = 1, xsize(3)
            do i = istart, iend
              ux(i, j, k) = rho(i, j, k) * byx1(i, k) + dpdxy1(i, k)
              uy(i, j, k) = rho(i, j, k) * byy1(i, k)
              uz(i, j, k) = rho(i, j, k) * byz1(i, k) + dpdzy1(i, k)
            enddo
          enddo
        endif

        if ((ny - (nym / dims(1))).eq.xstart(2)) then
          j = xsize(2)
          do k = 1, xsize(3)
            do i = istart, iend
              ux(i, j, k) = rho(i, j, k) * byxn(i, k) + dpdxyn(i, k)
              uy(i, j, k) = rho(i, j, k) * byyn(i, k)
              uz(i, j, k) = rho(i, j, k) * byzn(i, k) + dpdzyn(i, k)
            enddo
          enddo
        endif
      endif
    endif
  endif
  
  !****************************************************
  !********NCLZ==2*************************************
  !****************************************************
  !WE ARE IN X PENCIL!!!!!!!
  IF (nclz.EQ.2) THEN
    IF (itype.EQ.5) THEN
      ! determine the processor grid in use
      CALL MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
           dims, dummy_periods, dummy_coords, code)
      
      IF (dims(2).EQ.1) THEN
        k = 1
        DO j = 1, xsize(2)
          DO i = istart, iend
            ux(i, j, k) = rho(i, j, k) * bzx1(i, j) + dpdxz1(i, j)
            uy(i, j, k) = rho(i, j, k) * bzy1(i, j) + dpdyz1(i, j)
            uz(i, j, k) = rho(i, j, k) * bzz1(i, j)
          ENDDO
        ENDDO
        
        k = xsize(3)
        DO j = 1, xsize(2)
          DO i = istart, iend
            ux(i, j, k) = rho(i, j, k) * bzxn(i, j) + dpdxzn(i, j)
            uy(i, j, k) = rho(i, j, k) * bzyn(i, j) + dpdyzn(i, j)
            uz(i, j, k) = rho(i, j, k) * bzzn(i, j)
          ENDDO
        ENDDO
      ELSE
        IF (xstart(3).EQ.1) THEN
           k = 1
           DO j = 1, xsize(2)
              DO i = istart, iend
                 ux(i, j, k) = rho(i, j, k) * bzx1(i, j) + dpdxz1(i, j)
                 uy(i, j, k) = rho(i, j, k) * bzy1(i, j) + dpdyz1(i, j)
                 uz(i, j, k) = rho(i, j, k) * bzz1(i, j)
              ENDDO
           ENDDO
        ENDIF

        IF ((nz - (nzm / dims(2))).EQ.xstart(3)) THEN
           k = xsize(3)
           DO j = 1, xsize(2)
              DO i = istart, iend
                 ux(i, j, k) = rho(i, j, k) * bzxn(i, j) + dpdxzn(i, j)
                 uy(i, j, k) = rho(i, j, k) * bzyn(i, j) + dpdyzn(i, j)
                 uz(i, j, k) = rho(i, j, k) * bzzn(i, j)
              ENDDO
           ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  !##################################################### 

  return
end subroutine pre_correc

