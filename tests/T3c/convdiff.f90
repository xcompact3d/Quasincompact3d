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
subroutine convdiff(ux1,uy1,uz1,rho1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ux2,uy2,uz2,rho2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ux3,uy3,uz3,rho3,divu3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)
! 
!********************************************************************
USE param
USE variables
USE decomp_2d

USE MPI

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2 
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: rho1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: rho2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: rho3

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: divu1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: divu2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: divu3

real(mytype) :: ta1min, ta1min1, ta1max, ta1max1
real(mytype) :: tb1min, tb1min1, tb1max, tb1max1
real(mytype) :: tc1min, tc1min1, tc1max, tc1max1

integer :: ijk,nvect1,nvect2,nvect3,i,j,k
integer :: code
real(mytype) :: x,y,z

real(mytype), parameter :: ONETHIRD = 1._mytype / 3._mytype

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

!!! CM call test_min_max('ux1  ','In convdiff    ',ux1,size(ux1))
!!! CM call test_min_max('uy1  ','In convdiff    ',uy1,size(uy1))
!!! CM call test_min_max('uz1  ','In convdiff    ',uz1,size(uz1))

if (iskew==0) then !UROTU!
!WORK X-PENCILS
   call derx (ta1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (tb1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call transpose_x_to_y(ux1,ux2)
   call transpose_x_to_y(uy1,uy2)
   call transpose_x_to_y(uz1,uz2)
   call transpose_x_to_y(ta1,ta2)
   call transpose_x_to_y(tb1,tb2)
!WORK Y-PENCILS
   call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call dery (td2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call transpose_y_to_z(ux2,ux3)
   call transpose_y_to_z(uy2,uy3)
   call transpose_y_to_z(uz2,uz3)
   call transpose_y_to_z(ta2,ta3)
   call transpose_y_to_z(tb2,tb3)
   call transpose_y_to_z(tc2,tc3)
   call transpose_y_to_z(td2,td3)
!WORK Z-PENCILS
   call derz (te3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (tf3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   do ijk=1,nvect3
      ta3(ijk,1,1)=uz3(ijk,1,1)*(te3(ijk,1,1)-tb3(ijk,1,1))-&
           uy3(ijk,1,1)*(ta3(ijk,1,1)-tc3(ijk,1,1))
      tb3(ijk,1,1)=ux3(ijk,1,1)*(ta3(ijk,1,1)-tc3(ijk,1,1))-&
           uz3(ijk,1,1)*(td3(ijk,1,1)-tf3(ijk,1,1))
      tc3(ijk,1,1)=uy3(ijk,1,1)*(td3(ijk,1,1)-tf3(ijk,1,1))-&
           ux3(ijk,1,1)*(te3(ijk,1,1)-tb3(ijk,1,1))
   enddo
else !SKEW!
!WORK X-PENCILS
   td1 = rho1 * ux1 * ux1
   te1 = rho1 * uy1 * ux1
   tf1 = rho1 * uz1 * ux1
   
   call derx (tg1,td1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (th1,te1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (ti1,tf1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (td1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (te1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (tf1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

   ta1 = 0.5_mytype * (tg1 + rho1 * ux1 * td1)
   tb1 = 0.5_mytype * (th1 + rho1 * ux1 * te1)
   tc1 = 0.5_mytype * (ti1 + rho1 * ux1 * tf1)
   
   ! Quasi-skew symmetric terms
   call derx(tg1,rho1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   ta1 = ta1 + 0.5_mytype * ux1 * ux1 * tg1
   tb1 = tb1 + 0.5_mytype * uy1 * ux1 * tg1
   tc1 = tc1 + 0.5_mytype * uz1 * ux1 * tg1

   call transpose_x_to_y(ux1,ux2)
   call transpose_x_to_y(uy1,uy2)
   call transpose_x_to_y(uz1,uz2)
   call transpose_x_to_y(ta1,ta2)
   call transpose_x_to_y(tb1,tb2)
   call transpose_x_to_y(tc1,tc2)

   call transpose_x_to_y(rho1,rho2)
!WORK Y-PENCILS
   td2 = rho2 * ux2 * uy2
   te2 = rho2 * uy2 * uy2
   tf2 = rho2 * uz2 * uy2
   
   call dery (tg2,td2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) 
   call dery (th2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   call dery (ti2,tf2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) 
   call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   
   ta2 = ta2 + 0.5_mytype * (tg2 + rho2 * uy2 * td2)
   tb2 = tb2 + 0.5_mytype * (th2 + rho2 * uy2 * te2)
   tc2 = tc2 + 0.5_mytype * (ti2 + rho2 * uy2 * tf2)

   ! Quasi-skew symmetric terms
   call dery(th2,rho2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   ta2 = ta2 + 0.5_mytype * ux2 * uy2 * th2
   tb2 = tb2 + 0.5_mytype * uy2 * uy2 * th2
   tc2 = tc2 + 0.5_mytype * uz2 * uy2 * th2

   call transpose_y_to_z(ux2,ux3)
   call transpose_y_to_z(uy2,uy3)
   call transpose_y_to_z(uz2,uz3)
   call transpose_y_to_z(ta2,ta3)
   call transpose_y_to_z(tb2,tb3)
   call transpose_y_to_z(tc2,tc3)

   call transpose_y_to_z(rho2,rho3)
!WORK Z-PENCILS
   td3 = rho3 * ux3 * uz3
   te3 = rho3 * uy3 * uz3
   tf3 = rho3 * uz3 * uz3
   
   call derz (tg3,td3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call derz (th3,te3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call derz (ti3,tf3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

   ta3 = ta3 + 0.5_mytype * (tg3 + rho3 * uz3 * td3)
   tb3 = tb3 + 0.5_mytype * (th3 + rho3 * uz3 * te3)
   tc3 = tc3 + 0.5_mytype * (ti3 + rho3 * uz3 * tf3)

   ! Quasi-skew symmetric terms
   call derz(ti3,rho3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   ta3 = ta3 + 0.5_mytype * ux3 * uz3 * ti3
   tb3 = tb3 + 0.5_mytype * uy3 * uz3 * ti3
   tc3 = tc3 + 0.5_mytype * uz3 * uz3 * ti3

   ! Quasi-skew symmetric contribution from div(u)
   ta3 = ta3 + 0.5_mytype * ux3 * rho3 * divu3
   tb3 = tb3 + 0.5_mytype * uy3 * rho3 * divu3
   tc3 = tc3 + 0.5_mytype * uz3 * rho3 * divu3
endif
!ALL THE CONVECTIVE TERMS ARE IN TA3, TB3 and TC3

!!! CM call test_min_max('td3  ','In convdiff    ',td3,size(td3))
!!! CM call test_min_max('te3  ','In convdiff    ',te3,size(te3))
!!! CM call test_min_max('tf3  ','In convdiff    ',tf3,size(tf3))

tg3 = ta3
th3 = tb3
ti3 = tc3

!DIFFUSIVE TERMS IN Z
call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)

! Compute bulk shear contribution
! tg3, th3, ti3 available as work vectors
! TODO need to check ffzp, and whether last terms should be 1 or 0
call derz(tf3, divu3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0)
tc3 = tc3 - 2._mytype * ONETHIRD * tf3

!!! CM call test_min_max('ta3  ','In convdiff    ',ta3,size(ta3))
!!! CM call test_min_max('tb3  ','In convdiff    ',tb3,size(tb3))
!!! CM call test_min_max('tc3  ','In convdiff    ',tc3,size(tc3))


!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2)
call transpose_z_to_y(tb3,tb2)
call transpose_z_to_y(tc3,tc2)
call transpose_z_to_y(tg3,tg2)
call transpose_z_to_y(th3,th2)
call transpose_z_to_y(ti3,ti2)

call transpose_z_to_y(divu3, divu2)

!!! CM call test_min_max('tg2  ','In convdiff    ',tg2,size(tg2))
!!! CM call test_min_max('th2  ','In convdiff    ',th2,size(th2))
!!! CM call test_min_max('ti2  ','In convdiff    ',ti2,size(ti2))

!DIFFUSIVE TERMS IN Y
!-->for ux
if (istret.ne.0) then 
   call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      td2(i,j,k)=td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
   enddo
   enddo
   enddo
else
!!! CM    call test_min_max('ux2  ','In convdiff    ',ux2,size(ux2))
!!! CM    call test_min_max('di2  ','In convdiff    ',di2,size(di2))
!!! CM    write(*,*) ysize(1),ysize(2),ysize(3)
   call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif

!!! CM call test_min_max('td2  ','In convdiff    ',td2,size(td2))

!-->for uy
if (istret.ne.0) then 
   call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
   call dery (tf2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      te2(i,j,k)=te2(i,j,k)*pp2y(j)-pp4y(j)*tf2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0) 
endif
!-->for uz
if (istret.ne.0) then 
   call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   call dery (tj2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      tf2(i,j,k)=tf2(i,j,k)*pp2y(j)-pp4y(j)*tj2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif

!!! CM call test_min_max('td2  ','In convdiff    ',td2,size(td2))
!!! CM call test_min_max('te2  ','In convdiff    ',te2,size(te2))
!!! CM call test_min_max('tf2  ','In convdiff    ',tf2,size(tf2))

ta2 = ta2 + td2
tb2 = tb2 + te2
tc2 = tc2 + tf2

!!! CM call test_min_max('ta2  ','In convdiff    ',ta2,size(ta2))
!!! CM call test_min_max('tb2  ','In convdiff    ',tb2,size(tb2))
!!! CM call test_min_max('tc2  ','In convdiff    ',tc2,size(tc2))

! Compute bulk shear contribution
! td2, te2, tf2 avaiable as work vectors
call dery(te2, divu2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)
tb2 = tb2 - 2._mytype * ONETHIRD * te2

!WORK X-PENCILS
call transpose_y_to_x(ta2,ta1)
call transpose_y_to_x(tb2,tb1)
call transpose_y_to_x(tc2,tc1) !diff
call transpose_y_to_x(tg2,tg1)
call transpose_y_to_x(th2,th1)
call transpose_y_to_x(ti2,ti1) !conv

call transpose_y_to_x(divu2, divu1)

!DIFFUSIVE TERMS IN X
call derxx (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
call derxx (te1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
call derxx (tf1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)

ta1 = ta1 + td1
tb1 = tb1 + te1
tc1 = tc1 + tf1

! Compute bulk shear contribution
! td1, te1, tf1 available as work vectors
! TODO need to check ffzp, and whether last terms should be 1 or 0
call derx(td1, divu1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0)
ta1 = ta1 - 2._mytype * ONETHIRD * td1

!if (nrank==1) print *,'ATTENTION ATTENTION canal tournant',itime
!tg1(:,:,:)=tg1(:,:,:)-2./18.*uy1(:,:,:)
!th1(:,:,:)=th1(:,:,:)-2./18.*ux1(:,:,:)

!INTERMEDIATE SUM: DIFF TERMS + CONV TERMS
ta1 = xnu * ta1 - tg1
tb1 = xnu * tb1 - th1
tc1 = xnu * tc1 - ti1

!! Compute cross-shear
! NB ta1,tb1,tc1 cannot be touched!
! NB u2,u3 have already been updated, no need to transpose velocities!

! X - accumulate d(v,w)dx terms
call derx(te1, uy1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)
call derx(tf1, uz1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)

call transpose_x_to_y(te1, te2) ! te2 contains dvdx
call transpose_x_to_y(tf1, tf2) ! tf2 contains dwdx

! Y - accumulate dwdy terms
call dery(ti2, uz2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)

call transpose_y_to_z(tf2, tf3) ! tf3 contains dwdx
call transpose_y_to_z(ti2, ti3) ! ti3 contains dwdy

! Z - accumulate ddz terms
call derz(ta3, ux3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)
call derz(tb3, uy3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)
call derz(tc3, uz3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0)

! Z - compute ddz(dwdx, dwdy, dwdz)
call derz(td3, tf3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0)
call derz(te3, ti3, di3, sz, ffz, fsz, fwz, zsize(1), zsize(2), zsize(3), 0)
call derz(tf3, tc3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)

call transpose_z_to_y(td3, tg2) ! tg2 contains d2wdzdx
call transpose_z_to_y(te3, th2) ! th2 contains d2wdzdy
call transpose_z_to_y(tf3, ti2) ! ti2 contains d2wdzdz

call transpose_z_to_y(ta3, ta2) ! ta2 contains dudz
call transpose_z_to_y(tb3, tb2) ! tb2 contains dvdz

! Y - compute ddy(dvdx, dvdy, dvdz)

call dery(td2, te2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)
call dery(tc2, uy2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)
call dery(te2, tc2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)
call dery(tf2, tb2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)

td2 = td2 + tg2 ! td2 contains d2vdydx + d2wdzdx
te2 = te2 + th2 ! te2 contains d2vdydy + d2wdzdy
tf2 = tf2 + ti2 ! tf2 contains d2vdydz + d2wdzdz

call transpose_y_to_x(td2, td1) ! td1 contains d2vdydx + d2wdzdx
call transpose_y_to_x(te2, te1) ! te1 contains d2vdydy + d2wdzdy
call transpose_y_to_x(tf2, tf1) ! tf1 contains d2vdydz + d2wdzdz

call dery(td2, ux2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)

call transpose_y_to_x(td2, th1) ! tg1 contains dudy
call transpose_y_to_x(ta2, ti1) ! ti1 contains dudz

! X - compute ddx(dudx, dudy, dudz)

! First make some room to work!
ta1 = ta1 + xnu * td1
tb1 = tb1 + xnu * te1
tc1 = tc1 + xnu * tf1

call derx(tg1, ux1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0)

call derx(td1, tg1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)
call derx(te1, th1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0)
call derx(tf1, ti1, di1, sx, ffx, fsx, fwx, xsize(1), xsize(2), xsize(3), 0)

!! Finish off adding cross-stresses to shear stress
ta1 = ta1 + xnu * td1
tb1 = tb1 + xnu * te1
tc1 = tc1 + xnu * tf1

!! MMS Source term
!call momentum_source_mms(ta1,tb1,tc1)

ta1max=-1.e30_mytype
ta1min=+1.e30_mytype
tb1max=-1.e30_mytype
tb1min=+1.e30_mytype
tc1max=-1.e30_mytype
tc1min=+1.e30_mytype

ta1max = MAXVAL(ta1)
ta1min = MINVAL(ta1)
tb1max = MAXVAL(tb1)
tb1min = MINVAL(tb1)
tc1max = MAXVAL(tc1)
tc1min = MINVAL(tc1)
! do k=xstart(3),xend(3)
!    do j=xstart(2),xend(2)
!       do i=xstart(1),xend(1)
!          if (ta1(i,j,k).gt.ta1max) ta1max=ta1(i,j,k)
!          if (ta1(i,j,k).lt.ta1min) ta1min=ta1(i,j,k)
!          if (tb1(i,j,k).gt.tb1max) tb1max=tb1(i,j,k)
!          if (tb1(i,j,k).lt.tb1min) tb1min=tb1(i,j,k)
!          if (tc1(i,j,k).gt.tc1max) tc1max=tc1(i,j,k)
!          if (tc1(i,j,k).lt.tc1min) tc1min=tc1(i,j,k)
!       enddo
!    enddo
! enddo

call MPI_REDUCE(ta1max,ta1max1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(ta1min,ta1min1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(tb1max,tb1max1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(tb1min,tb1min1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(tc1max,tc1max1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(tc1min,tc1min1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)


!!! CM if (nrank==0) then
!!! CM    write(*,*) 'In convdiff ta1',ta1max1,ta1min1
!!! CM    write(*,*) 'In convdiff tb1',tb1max1,tb1min1
!!! CM    write(*,*) 'In convdiff tc1',tc1max1,tc1min1
!!! CM endif

end subroutine convdiff


!************************************************************
!
subroutine scalar(ux1,uy1,uz1,phi1,phis1,phiss1,di1,ta1,tb1,tc1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,epsi)
!
!************************************************************

USE param
USE variables
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,phis1,&
                                              phiss1,di1,ta1,tb1,tc1,td1,epsi
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2,di2,ta2,tb2,tc2,td2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,di3,ta3,tb3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz
real(mytype) :: x,y,z

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

!X PENCILS
do ijk=1,nvect1
   ta1(ijk,1,1)=ux1(ijk,1,1)*phi1(ijk,1,1)
enddo
call derx (tb1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derxx (ta1,phi1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)

call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)

!Y PENCILS
do ijk=1,nvect2
   ta2(ijk,1,1)=uy2(ijk,1,1)*phi2(ijk,1,1)
enddo
call dery (tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
if (istret.ne.0) then 
   call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      ta2(i,j,k)=ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
   enddo
   enddo
   enddo
else
   call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
endif

call transpose_y_to_z(phi2,phi3)
call transpose_y_to_z(uz2,uz3)

!Z PENCILS
do ijk=1,nvect3
   ta3(ijk,1,1)=uz3(ijk,1,1)*phi3(ijk,1,1)
enddo
call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
call derzz (ta3,phi3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)

call transpose_z_to_y(ta3,tc2)
call transpose_z_to_y(tb3,td2)

!Y PENCILS ADD TERMS
do ijk=1,nvect2
   tc2(ijk,1,1)=tc2(ijk,1,1)+ta2(ijk,1,1)
   td2(ijk,1,1)=td2(ijk,1,1)+tb2(ijk,1,1)
enddo

call transpose_y_to_x(tc2,tc1)
call transpose_y_to_x(td2,td1)

!X PENCILS ADD TERMS
do ijk=1,nvect1
   ta1(ijk,1,1)=ta1(ijk,1,1)+tc1(ijk,1,1) !SECOND DERIVATIVE
   tb1(ijk,1,1)=tb1(ijk,1,1)+td1(ijk,1,1) !FIRST DERIVATIVE
enddo
 
do ijk=1,nvect1
   ta1(ijk,1,1)=xnu/sc*ta1(ijk,1,1)-tb1(ijk,1,1) 
enddo

!TIME ADVANCEMENT
nxyz=xsize(1)*xsize(2)*xsize(3)  

if ((nscheme.eq.1).or.(nscheme.eq.2)) then
   if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
        (nscheme.eq.2.and.itr.eq.1)) then
      do ijk=1,nxyz
         phi1(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      do ijk=1,nxyz
         phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
endif
endif

if (nscheme.eq.3) then 
if (nrank==0) print *,'Not ready'
stop 
endif

if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'start with Euler',itime
      do ijk=1,nxyz !start with Euler
         phi1(ijk,1,1)=dt*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
         if (nrank==0) print *,'then with AB2',itime
         do ijk=1,nxyz
            phi1(ijk,1,1)=1.5_mytype*dt*ta1(ijk,1,1)-0.5_mytype*dt*phis1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo 
      else
         do ijk=1,nxyz
            phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+&
                 cdt(itr)*phiss1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo
      endif
   endif
endif


end subroutine scalar

!!--------------------------------------------------------------------
!!  Subroutine: density
!!
!! Description: Advances density in time for LMN.
!!
!!       Notes: The "diffusion" term in the continuity equation is
!!              given as:
!!
!!                (1 / (Re Pr T)) div(kappa grad(T))
!!
!!              however, we have p = rho T and grad(p) = 0, where p is
!!              the thermodynamic pressure in the LMN approximation.
!!              Therefore can make the substitution T = p / rho giving
!!
!!                (1 / (Re Pr (p / rho))) div(kappa grad(p / rho))
!!                  => (rho / (Re Pr)) div(kappa grad(1 / rho))
!!
!!              making use of grad(p) = 0. This saves
!!              memory/communication as would need to compute, store
!!              and transpose temperature array when we do this with
!!              density anyway.
!!--------------------------------------------------------------------
SUBROUTINE density(ux1, uy1, uz1, rho1, rhos1, rhoss1, di1, ta1, tb1, tc1, td1, &
     uy2, uz2, rho2, di2, ta2, tb2, tc2, td2, &
     uz3, rho3, divu3, di3, ta3, tb3, &
     epsi)
  
  USE param
  USE variables
  USE decomp_2d
  
  IMPLICIT NONE
  
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: rho1, rhos1, rhoss1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: temperature1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: di1, ta1, tb1, tc1, td1, epsi
  
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: uy2, uz2
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: rho2
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: temperature2
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: di2, ta2, tb2, tc2, td2
  
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: uz3
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: rho3
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: divu3
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: di3, ta3, tb3

INTEGER :: i,j,k
REAL(mytype) :: x,y
REAL(mytype) :: tnext

  !------------------------------------------------------------------------
  ! X PENCILS
  ! ta1 = diffusion
  ! tb1 = advection

  ! Advection term (non-conservative)
  CALL derx (tb1, rho1, di1, sx, ffxp, fsxp, fwxp, xsize(1), xsize(2), xsize(3), 1)
  tb1 = ux1 * tb1

  ! Go to Y
  CALL transpose_x_to_y(rho1, rho2)
  CALL transpose_x_to_y(temperature1, temperature2)
  CALL transpose_x_to_y(uy1, uy2)
  CALL transpose_x_to_y(uz1, uz2)

  !------------------------------------------------------------------------
  !Y PENCILS
  ! ta2 = diffusion
  ! tb2 = advection

  ! Advection term (non-conservative)
  CALL dery (tb2, rho2, di2, sy, ffyp, fsyp, fwyp, ppy, ysize(1), ysize(2), ysize(3), 1)
  tb2 = uy2 * tb2

  ! Go to Z
  CALL transpose_y_to_z(rho2, rho3)
  CALL transpose_y_to_z(uz2, uz3)

  !------------------------------------------------------------------------
  ! Z PENCILS
  ! ta3 = diffusion
  ! tb3 = advection

  ! Advection term (non-conservative)
  CALL derz (tb3, rho3, di3, sz, ffzp, fszp, fwzp, zsize(1), zsize(2), zsize(3), 1)
  tb3 = uz3 * tb3

  ! We can now add div u
  tb3 = tb3 + rho3 * divu3

  ! Get back to Y
  CALL transpose_z_to_y(tb3, td2)

  !------------------------------------------------------------------------
  !Y PENCILS ADD TERMS
  td2 = td2 + tb2

  ! Get back to X
  CALL transpose_y_to_x(td2, td1)

  !------------------------------------------------------------------------
  ! X PENCILS ADD TERMS
  tb1 = tb1 + td1 !FIRST DERIVATIVE (CONV)
  
  ta1 = -tb1

  !! MMS Source term
  CALL density_source_mms(ta1)

  ! ! Now store velocity as momentum
  ! ux1 = rho1 * ux1
  ! uy1 = rho1 * uy1
  ! uz1 = rho1 * uz1
  
  !------------------------------------------------------------------------
  ! TIME ADVANCEMENT

  IF ((nscheme.EQ.1).OR.(nscheme.EQ.2)) THEN
    !! AB2 or RK3
    IF ((nscheme.EQ.1.AND.itime.EQ.1.AND.ilit.EQ.0).OR.&
         (nscheme.EQ.2.AND.itr.EQ.1)) THEN
      rho1 = rho1 + gdt(itr) * ta1
    ELSE
      rho1 = rho1 + adt(itr) * ta1 + bdt(itr) * rhos1
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
      rho1 = rho1 + dt * ta1
    ELSE
      IF  ((itime.EQ.2).AND.(ilit.EQ.0)) THEN
        IF (nrank.EQ.0) THEN
          PRINT *, 'then with AB2', itime
        ENDIF
        rho1 = rho1 - 0.5_mytype * dt * (rhos1 - 3._mytype * ta1)
      ELSE
        rho1 = rho1 + adt(itr) * ta1 + bdt(itr) * rhos1 + cdt(itr) &
             * rhoss1
      ENDIF

      !! Update oldold stage
      rhoss1 = rhos1
    ENDIF
  ENDIF

  !! Update old stage
  rhos1 = ta1

      ! XXX Temporary: move velocity forward in time
      if ((nscheme.eq.1).or.(nscheme.eq.2)) then
        !! AB2 or RK3
        if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
          (nscheme.eq.2.and.itr.eq.1)) then
          tnext = t + gdt(itr)
        else
          tnext = t + adt(itr)
        endif
      endif
      do k = 1,xsize(3)
      do j = 1,xsize(2)
      do i = 1,xsize(1)
        x = float(i + xstart(1) - 2) * dx - xlx / 2._mytype
        y = float(j + xstart(2) - 2) * dy - yly / 2._mytype
        ux1(i,j,k) = 2._mytype * (5._mytype - 1._mytype) * SIN(PI * 2._mytype * y) &
            * SIN(PI * 2._mytype * tnext) * COS(PI * 2._mytype * x) & 
            / (2._mytype * 2._mytype * ((1._mytype + 5._mytype) + (1._mytype - 5._mytype) &
            * SIN(PI * 2._mytype * x) * SIN(PI * 2._mytype * y) * COS(PI * 2._mytype * tnext)))
        uy1(i,j,k) = 2._mytype * (5._mytype - 1._mytype) * SIN(PI * 2._mytype * x) &
            * SIN(PI * 2._mytype * tnext) * COS(PI * 2._mytype * y) &
            / (2._mytype * 2._mytype * ((1._mytype + 5._mytype) + (1._mytype - 5._mytype) &
            * SIN(PI * 2._mytype * x) * SIN(PI * 2._mytype * y) * COS(PI * 2._mytype * tnext)))
      enddo
      enddo
      enddo
  
ENDSUBROUTINE density

!!--------------------------------------------------------------------
!!  SUBROUTINE: calc_divu
!! DESCRIPTION: In LMN the divergence of velocity is given in terms of
!!              the temperature field, ensuring the gradient of
!!              thermodynamic pressure is zero.
!!--------------------------------------------------------------------
SUBROUTINE calc_divu(ta1, rho1, temperature1, di1, &
     ta2, tb2, tc2, rho2, temperature2, di2, &
     divu3, ta3, rho3, temperature3, di3, &
     pressure0)

  USE param
  USE variables
  USE decomp_2d
  
  IMPLICIT NONE

  INTEGER i, j, k

  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ta1, di1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: rho1
  REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(OUT) :: temperature1
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)) :: ta2, tb2, tc2, di2
  REAL(mytype), DIMENSION(ysize(1), ysize(2), ysize(3)), INTENT(OUT) :: rho2, temperature2
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)) :: ta3, di3
  REAL(mytype), DIMENSION(zsize(1), zsize(2), zsize(3)), INTENT(OUT) :: divu3, rho3, temperature3
  REAL(mytype), INTENT(IN) :: pressure0

  REAL(mytype) :: invpressure0

  !-------------------------------------------------------------------
  ! X pencil
  !-------------------------------------------------------------------

  ! Update temperature
  CALL calctemp_eos(temperature1, rho1, pressure0, xsize)

  ! Calculate divergence of velocity
  CALL derxx (ta1, temperature1, di1, sx, sfxp, ssxp, swxp, xsize(1), xsize(2), xsize(3), 1)

  ! Transpose to Y
  CALL transpose_x_to_y(rho1, rho2)
  CALL transpose_x_to_y(ta1, ta2)

  !-------------------------------------------------------------------
  ! Y pencil
  !-------------------------------------------------------------------

  ! Update temperature
  CALL calctemp_eos(temperature2, rho2, pressure0, ysize)

  ! Calculate divergence of velocity
  IF(istret.NE.0) THEN
    CALL deryy (tb2, temperature2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
    CALL dery (tc2, temperature2, di2, sy, ffy, fsy, fwy, ppy, ysize(1), ysize(2), ysize(3), 0)
    DO k = 1, ysize(3)
      DO j = 1, ysize(2)
        DO i = 1, ysize(1)
          tb2(i, j, k) = tb2(i, j, k) * pp2y(j) - pp4y(j) * tc2(i, j, k)
        ENDDO
      ENDDO
    ENDDO
  ELSE
    CALL deryy (tb2, temperature2, di2, sy, sfyp, ssyp, swyp, ysize(1), ysize(2), ysize(3), 1)
  ENDIF
  ta2 = ta2 + tb2

  ! Transpose to Z
  CALL transpose_y_to_z(rho2, rho3)
  CALL transpose_y_to_z(ta2, ta3)

  !-------------------------------------------------------------------
  ! Z pencil
  !-------------------------------------------------------------------

  ! Update temperature
  CALL calctemp_eos(temperature3, rho3, pressure0, zsize)

  ! Calculate divergence of velocity
  CALL derzz (divu3, temperature3, di3, sz, sfzp, sszp, swzp, zsize(1), zsize(2), zsize(3), 1)
  divu3 = divu3 + ta3
  divu3 = (xnu / pr) * divu3

  ! XXX add dpdt and additional source terms

  ! divu3 = divu3 / (rho3 * temperature3)
  invpressure0 = 1._mytype / pressure0
  divu3 = invpressure0 * divu3 ! rho*T = pressure0 = constant (in space)

  !-------------------------------------------------------------------
  ! XXX Density and temperature fields are now up to date in all
  !     pencils.
  !     Divergence of velocity is only known in Z pencil.
  !-------------------------------------------------------------------
  
ENDSUBROUTINE calc_divu

!!--------------------------------------------------------------------
!!  SUBROUTINE: calctemp_eos
!! DESCRIPTION: Given the new density field, calculate temperature
!!              using the equation of state.
!!--------------------------------------------------------------------
SUBROUTINE calctemp_eos(temperature1, rho1, pressure0, arrsize)

  USE variables
  USE decomp_2d

  IMPLICIT NONE

  INTEGER, DIMENSION(3), INTENT(IN) :: arrsize

  REAL(mytype), DIMENSION(arrsize(1), arrsize(2), arrsize(3)), INTENT(OUT) :: temperature1
  REAL(mytype), DIMENSION(arrsize(1), arrsize(2), arrsize(3)), INTENT(IN) :: rho1

  REAL(mytype), INTENT(IN) :: pressure0

  !!------------------------------------------------------------------
  !! Very simple EOS
  !!   p = rho T
  !!------------------------------------------------------------------
  temperature1 = pressure0 / rho1
  
ENDSUBROUTINE calctemp_eos
  
!!--------------------------------------------------------------------
!! SUBROUTINE: density_source_mms
!! DESCIPTION: Computes the source term for the density equation in
!!             Method of Manufactured Solutions test and adds it to
!!             the stress/diffusion term.
SUBROUTINE density_source_mms(mms)

  USE var

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1),xsize(2),xsize(3)) :: mms

  REAL(mytype) :: x,y,z
  INTEGER :: i,j,k

  REAL(mytype) :: rhomms, rhoa, rhob
  REAL(mytype) :: Tmms
  REAL(mytype) :: press0
  REAL(mytype) :: wavnum, omega
  REAL(mytype) :: SrhoX, SrhoY, SrhoT
  REAL(mytype) :: MMSource
  REAL(mytype) :: SINX, SINY, SINT
  REAL(mytype) :: COSX, COSY, COST

  press0 = 1._mytype
  rhoa = 5._mytype
  rhob = 1._mytype
  wavnum = 2._mytype
  omega = 2._mytype

  DO k = 1,xsize(3)
    z = float(k + xstart(3) - 2) * dz - zlz / 2._mytype
    DO j = 1,xsize(2)
      y = float(j + xstart(2) - 2) * dy - yly / 2._mytype
      DO i = 1, xsize(1)
        x = float(i + xstart(1) - 2) * dx - xlx / 2._mytype

        SINX = SIN(PI * wavnum * x)
        SINY = SIN(PI * wavnum * y)
        SINT = SIN(PI * omega * t)
        COSX = COS(PI * wavnum * x)
        COSY = COS(PI * wavnum * y)
        COST = COS(PI * omega * t)

	rhomms = ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST) / 2._mytype
        Tmms = press0 / rhomms
        
        !!
        !! Compute nabla.nabla T
        !!

        ! d/dx( d/dx T )
        SrhoX = -rhoa**2 * (SINX**2 * SINY * COST + SINX &
             - 2._mytype * SINY * COST) &
             + rhoa * rhob * (2._mytype * SINX**2 * SINY * COST &
             - 4._mytype * SINY * COST) &
             - rhob**2 * (SINX**2 * SINY * COST - SINX &
             - 2._mytype * SINY * COST)
        SrhoX = ((2._mytype * (PI * wavnum)**2 * press0 * SINY * COST) &
             / ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**3) &
             * SrhoX

        ! d/dy( d/dy T )
        SrhoY = -rhoa**2 * (SINX * SINY**2 * COST + SINY &
             - 2._mytype * SINX * COST) &
             + rhoa * rhob * (2._mytype * SINX * SINY**2 * COST & 
             - 4._mytype * SINX * COST) &
             - rhob**2 * (SINX * SINY**2 * COST - SINY &
             - 2._mytype * SINX * COST)
        SrhoY = ((2._mytype * (PI * wavnum)**2 * press0 * SINX * COST) &
             / ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**3) &
             * SrhoY

        !!
        !! Compute divu = (1 / (Re Pr T)) * nabla.nabla T / rho
        !!
        MMSource = SrhoX + SrhoY
        MMSource = MMSource * (xnu / (pr * (rhomms * Tmms)))

        !!
        !! Finally: S_rho = rho * divu
        !!
        MMSource = rhomms * MMSource

        !!
        !! Compute u.grad(rho)
        !!
        SrhoX = (PI * omega) * (rhoa**2 - 2._mytype * rhoa * rhob + rhob**2) &
             * COSX**2 * SINY**2 * SINT * COST
	SrhoX = SrhoX / (4._mytype * ((rhoa - rhob) * SINX * SINY * COST &
		- (rhoa + rhob)))

        SrhoY = (PI * OMEGA) * (rhoa**2 - 2._mytype * rhoa * rhob + rhob**2) &
             * SINX**2 * COSY**2 * SINT * COST
	SrhoY = SrhoY / (4._mytype * ((rhoa - rhob) * SINX * SINY * COST &
		- (rhoa + rhob)))

        MMSource = MMSource + (SrhoX + SrhoY)

        !!
        !! Compute drhodt
        !!
        SrhoT = 0.5_mytype * (PI * omega) * (rhoa - rhob) &
             * SINX * SINY * SINT
        MMSource = MMSource + SrhoT

        mms(i,j,k) = mms(i,j,k) + MMSource
      ENDDO ! End loop over i
    ENDDO ! End loop over j
  ENDDO ! End loop over k

ENDSUBROUTINE density_source_mms
  
!!--------------------------------------------------------------------
!! SUBROUTINE: momentum_source_mms
!! DESCIPTION: Computes the source term for the momentum equations in
!!             Method of Manufactured Solutions test and adds it to
!!             the stress/diffusion term
SUBROUTINE momentum_source_mms(mmsx1, mmsy1, mmsz1)

  USE var

  IMPLICIT NONE

  REAL(mytype), DIMENSION(xsize(1),xsize(2),xsize(3)) :: mmsx1, mmsy1, mmsz1

  REAL(mytype) :: x,y,z
  INTEGER :: i,j,k

  REAL(mytype) :: wavnum, omega

  REAL(mytype) :: umms, vmms, wmms
  REAL(mytype) :: rhomms, rhoa, rhob
  REAL(mytype) :: press0
  REAL(mytype) :: Tmms
  REAL(mytype) :: divumms
  REAL(mytype) :: MMSource

  REAL(mytype) :: SINX, SINY, SINT
  REAL(mytype) :: COSX, COSY, COST

  press0 = 1._mytype

  rhoa = 5._mytype
  rhob = 1._mytype
  omega = 2._mytype
  wavnum = 2._mytype
  
  DO k = 1,xsize(3)
    z = float(k + xstart(3) - 2) * dz - zlz / 2._mytype
    DO j = 1,xsize(2)
      y = float(j + xstart(2) - 2) * dy - yly / 2._mytype
      DO i = 1, xsize(1)
        x = float(i + xstart(1) - 2) * dx - xlx / 2._mytype

        SINX = SIN(PI * wavnum * x)
        SINY = SIN(PI * wavnum * y)
        SINT = SIN(PI * omega * t)
        COSX = COS(PI * wavnum * x)
        COSY = COS(PI * wavnum * y)
        COST = COS(PI * omega * t)
        
        !! XMOM

        ! Transient
        MMSource = ((PI * omega**2) / (4._mytype * wavnum)) * (rhoa - rhob) &
             * COSX * SINY * COST
        mmsx1(i, j, k) = mmsx1(i, j, k) + MMSource

        ! Advection
        MMSource = (rhoa - rhob) * COSX**2 * SINY * COST &
             - 3._mytype * ((rhoa + rhob) + (rhob - rhoa) &
             * SINX * SINY * COST) * SINX
        MMSource = ((PI * omega**2 * (rhoa - rhob)**2 * COSX * SINY**2 * SINT**2) &
             / (16._mytype * wavnum * ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST))**2) &
             * MMSource
        mmsx1(i, j, k) = mmsx1(i, j, k) + MMSource

        MMSource = (rhoa + rhob) * COSY**2 &
             + ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST) &
             * COS(2._mytype * PI * wavnum * y)
        MMSource = ((PI * omega**2 * (rhoa - rhob)**2 * SINX * COSX * SINT**2) &
             / (16._mytype * wavnum * ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST))**2) &
             * MMSource
        mmsx1(i, j, k) = mmsx1(i, j, k) + MMSource

        ! The first half of the viscous stress tensor (grad u + grad^T u)
        MMSource = rhoa**3 * (SINX * COSY**2 * COST + SINY**3 * COST**2 - SINY) &
             - rhoa**2 * rhob * (SINX * COSY**2 * COST + 3._mytype * SINY**3 * COST**2 + SINY) &
             - rhoa * rhob**2 * (SINX * COSY**2 * COST - 3._mytype * SINY**3 * COST**2 - SINY) &
             + rhob**3 * (SINX * COSY**2 * COST - SINY**3 * COST**2 + SINY)
        MMSource = -((PI**2 * wavnum * omega * COSX * SINT) * xnu &
             / ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**3) * MMSource
        mmsx1(i, j, k) = mmsx1(i, j, k) + MMSource

        MMSource = rhoa**3 * (SINX * COSY**2 * COST + SINY**3 * COST**2 - SINY) &
             - rhoa**2 * rhob * (SINX * COSY**2 * COST + 3._mytype * SINY**3 * COST**2 + SINY) &
             - rhoa * rhob**2 * (SINX * COSY**2 * COST - 3._mytype * SINY**3 * COST**2 - SINY) &
             + rhob**3 * (SINX * COSY**2 * COST - SINY**3 * COST**2 + SINY)
        MMSource = -((PI**2 * wavnum * omega * COSX * SINT) * xnu &
             / ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**3) * MMSource
        mmsx1(i, j, k) = mmsx1(i, j, k) + MMSource

        ! The bulk component of viscous stress tensor
        MMSource = rhoa**3 * (4._mytype * SINX * COSY**2 * COST - 2._mytype * SINX * COST &
             + COSX**2 * SINY**3 * COST**2 + 2._mytype * SINY**3 * COST**2 &
             - COSX**2 * SINY * COST**2 + SINY * COST**2 - SINY) &
             - rhoa**2 * rhob * (4._mytype * SINX * COSY**2 * COST - 2._mytype * SINX * COST &
             + 3._mytype * COSX**2 * SINY**3 * COST**2 + 6._mytype * SINY**3 * COST**2 &
             - 3._mytype * COSX**2 * SINY * COST**2 + 3._mytype * SINY * COST**2 + SINY) &
             - rhoa * rhob**2 * (4._mytype * SINX * COSY**2 * COST - 2._mytype * SINX * COST &
             - 3._mytype * COSX**2 * SINY**3 * COST**2 - 6._mytype * SINY**3 * COST**2 &
             + 3._mytype * COSX**2 * SINY * COST**2 - 3._mytype * SINY * COST**2 - SINY) &
             + rhob**3 * (4._mytype * SINX * COSY**2 * COST - 2._mytype * SINX * COST &
             - COSX**2 * SINY**3 * COST**2 - 2._mytype * SINY**3 * COST**2 &
             + COSX**2 * SINY * COST**2 - SINY * COST**2 + SINY)
        MMSource = ((8._mytype * (PI * wavnum)**3 * COSX * COST) * xnu**2 &
             / (3._mytype * pr * ((rhoa - rhob) * SINX * SINY * COST - (rhoa + rhob))**4))
        mmsx1(i, j, k) = mmsx1(i, j, k) + MMSource

        !! YMOM

        ! Transient
        MMSource = ((PI * omega**2) / (4._mytype * wavnum)) * (rhoa - rhob) &
             * SINX * COSY * COST
        mmsy1(i, j, k) = mmsy1(i, j, k) + MMSource
        
        ! Advection
        MMSource = (rhoa + rhob) * COSX**2 &
             + ((rhoa + rhob) + (rhob - rhoa) * (SINX * SINY * COST) &
             * COS(2._mytype * PI * wavnum * x))
        MMSource = ((PI * omega**2 * (rhoa - rhob)**2 * SINY * COSY * SINT**2) &
             / (16._mytype * wavnum * ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**2))
        mmsy1(i, j, k) = mmsy1(i, j, k) + MMSource

        MMSource = (rhoa - rhob) * SINX * COSY**2 * COST &
             - 3._mytype * ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST) * SINY
        MMSource = ((PI * omega**2 * (rhoa - rhob)**2 * SINX**2 * COSY * SINT**2) &
             / (16._mytype * wavnum * ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**2))
        mmsy1(i, j, k) = mmsy1(i, j, k) + MMSource

        ! The first half of the viscous stress tensor (grad u + grad^T u)
        MMSource = rhoa**3 * (SINX**3 * COST**2 - SINX + COSX**2 * SINY * COST) &
             - rhoa**2 * rhob * (3._mytype * SINX**3 * COST**2 + SINX + COSX**2 * SINY * COST &
             + rhoa * rhob**2 * (3._mytype * SINX**3 * COST**2 + SINX - COSX**2 * SINY * COST) &
             - rhob**3 * (SINX**3 * COST**2 - SINX - COSX**2 * SINY * COST))
        MMSource = -((PI**2 * wavnum * omega * COSY * SINT) * xnu &
             / ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**3)
        mmsy1(i, j, k) = mmsy1(i, j, k) + MMSource

        MMSource = rhoa**3 * (SINX**3 * COST**2 - SINX + COSX**2 * SINY * COST) &
             - rhoa**2 * rhob * (3._mytype * SINX**3 * COST**2 + SINX + COSX**2 * SINY * COST &
             + rhoa * rhob**2 * (3._mytype * SINX**3 * COST**2 + SINX - COSX**2 * SINY * COST) &
             - rhoa**3 * (SINX**3 * COST**2 - SINX - COSX**2 * SINY * COST))
        MMSource = -((PI**2 * wavnum * omega * COSY * SINT) * xnu &
             / ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**3)
        
        ! The bulk component of viscous stress tensor
        MMSource = rhoa**3 * (SINX**3 * COSY**2 * COST**2 + 2._mytype * SINX**3 * COST**2 &
             - SINX * COSY**2 * COST**2 + SINX * COST**2 - SINX &
             + 4._mytype * COSX**2 * SINY * COST - 2._mytype * SINY * COST) &
             - rhoa**2 * rhob * (3._mytype * SINX**3 * COSY**2 * COST**2 &
             + 6._mytype * SINX**3 * COST**2 - 3._mytype * SINX * COSY**2 * COST**2 &
             + SINX * COST**2 + SINX + 4._mytype * COSX**2 * SINY * COST &
             -2._mytype * SINY * COST) &
             + rhoa * rhob**2 * (3._mytype * SINX**3 * COSY**2 * COST**2 &
             + 6._mytype * SINX**3 * COST**2 - 3._mytype * SINX * COSY**2 * COST**2 &
             + 3._mytype * SINX * COST**2 + SINX - 4._mytype * COSX**2 * SINY * COST &
             + 2._mytype * SINY * COST) &
             - rhob**3 * (SINX**3 * COSY**2 * COST**2 + 2._mytype * SINX**2 * COST**2 &
             - SINX * COSY**2 * COST**2 + SINX * COST**2 - SINX &
             - 4._mytype * COSX**2 * SINY * COST + 2._mytype * SINY * COST)
        MMSource = ((8._mytype * (PI * wavnum)**3 * COSY * COST) * (xnu**2) &
             / (-3._mytype * pr * ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**4))
        mmsy1(i, j, k) = mmsy1(i, j, k) + MMSource

        !! Correction for quasi-skew symmetry
        MMSource = (rhoa - rhob) * (-rhoa**2 * (SINX * SINY + COSX**2 * COSY**2 * COST - COST) &
             + 2._mytype * rhoa * rhob * COSX**2 * COSY**2 * COST &
             - 2._mytype * rhoa * rhob * COST &
             + rhob**2 * (SINX * SINY - COSX**2 * COSY**2 * COST + COST))
        MMSource = ((PI**2 * wavnum * omega * COSX * SINY * SINT * COST) &
             / (2._mytype * (pr / xnu) * ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**3))
        mmsx1(i, j, k) = mmsx1(i, j, k) + MMSource

        MMSource = (rhoa - rhob) * (-rhoa**2 * (SINX * SINY + COSX**2 * COSY**2 * COST - COST) &
             + 2._mytype * rhoa * rhob * COSX**2 * COSY**2 * COST &
             - 2._mytype * rhoa * rhob * COST &
             + rhob**2 * (SINX * SINY - COSX**2 * COSY**2 * COST + COST))
        MMSource = ((PI**2 * wavnum * omega * SINX * COSY * SINT * COST) * xnu &
             / (2._mytype * pr * ((rhoa + rhob) + (rhob - rhoa) * SINX * SINY * COST)**3))
        mmsy1(i, j, k) = mmsy1(i, j, k) + MMSource

      ENDDO ! End loop over i
    ENDDO ! End loop over j
  ENDDO ! End loop over k

ENDSUBROUTINE momentum_source_mms
