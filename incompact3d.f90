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

PROGRAM incompact3d

USE decomp_2d
USE decomp_2d_poisson
use decomp_2d_io
USE variables
USE param
USE var
USE MPI
USE IBM
USE derivX
USE derivZ

implicit none

integer :: code,nlock,i,j,k,ii,bcx,bcy,bcz,fh,ierror
real(mytype) :: x,y,z,tmp1
double precision :: t1,t2
character(len=20) :: filename

TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

CALL MPI_INIT(code)
call decomp_2d_init(nx,ny,nz,p_row,p_col)
!start from 1 == true
call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)
call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)
call parameter()

call init_variables

call schemes()

!!! CM call test_min_max('di2  ','In main        ',di2,size(di2))

if (nclx==0) then
   bcx=0
else
   bcx=1
endif
if (ncly==0) then
   bcy=0
else
   bcy=1
endif
if (nclz==0) then
   bcz=0
else
   bcz=1
endif

call decomp_2d_poisson_init(bcx,bcy,bcz)

call decomp_info_init(nxm,nym,nzm,phG)

!!! CM call test_min_max('di2  ','In main 1      ',di2,size(di2))

!if you want to collect 100 snapshots randomly on 50000 time steps
!call collect_data() !it will generate 100 random time steps

if (ilit==0) call init(ux1,uy1,uz1,rho1,rhos1,rhoss1,ep1,phi1,&
     gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)  
if (ilit==1) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
     px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,0)

call test_speed_min_max(ux1,uy1,uz1)
call test_density_min_max(rho1)
! call test_density_min_max(rho1)
! call test_temperature_min_max(temperature1)
if (iscalar==1) call test_scalar_min_max(phi1)

!array for stat to zero
umean=0._mytype;vmean=0._mytype;wmean=0._mytype
uumean=0._mytype;vvmean=0._mytype;wwmean=0._mytype
uvmean=0._mytype;uwmean=0._mytype;vwmean=0._mytype
phimean=0._mytype;phiphimean=0._mytype

!!! CM call test_min_max('di2  ','In main 2      ',di2,size(di2))

t1 = MPI_WTIME()

!div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
call decomp_info_init(nxm, nym, nzm, ph1)
call decomp_info_init(nxm, ny, nz, ph4)

!gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
call decomp_info_init(nxm, ny, nz, ph2)  
call decomp_info_init(nxm, nym, nz, ph3) 

itime=0
call VISU_INSTA(ux1,uy1,uz1,rho1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)

call VISU_INSTB(ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)

!!! CM call test_min_max('di2  ','In main 3      ',di2,size(di2))

do itime=ifirst,ilast

   t=(itime-1)*dt
   if (nrank==0) then

      write(*,1001) itime,t
1001  format('Time step =',i7,', Time unit =',F9.3)
   endif


!!! CM call test_min_max('di2  ','In main 4      ',di2,size(di2))
   do itr=1,iadvance_time

     ! Update density
     call density(ux1,uy1,uz1,rho1,rhos1,rhoss1,di1,tg1,th1,ti1,td1,&
          uy2,uz2,rho2,di2,ta2,tb2,tc2,td2,&
          uz3,rho3,di3,ta3,tb3,ep1)

      if (nclx.eq.2) then
         call inflow (ux1,uy1,uz1,phi1) !X PENCILS
         call outflow(ux1,uy1,uz1,phi1) !X PENCILS 
      endif

      !X-->Y-->Z-->Y-->X
      call convdiff(ux1,uy1,uz1,rho1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)
           
      if (iscalar==1) then
         call scalar(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
              uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,ep1) 
      endif

      !X PENCILS
      call intt (ux1,uy1,uz1,gx1,gy1,gz1,hx1,hy1,hz1,ta1,tb1,tc1) 

!!! CM call test_min_max('ux1  ','In main intt   ',ux1,size(ux1))
!!! CM call test_min_max('uy1  ','In main intt   ',uy1,size(uy1))
!!! CM call test_min_max('uz1  ','In main intt   ',uz1,size(uz1))

      call pre_correc(ux1,uy1,uz1)

!!! CM call test_min_max('ux1  ','In main pre_   ',ux1,size(ux1))
!!! CM call test_min_max('uy1  ','In main pre_   ',uy1,size(uy1))
!!! CM call test_min_max('uz1  ','In main pre_   ',uz1,size(uz1))

!!$      if (ivirt==1) then !solid body old school
!!$         !we are in X-pencil
!!$         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
!!$         call body(ux1,uy1,uz1,ep1)
!!$         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
!!$      endif

      !X-->Y-->Z
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,1)       

!!! CM call test_min_max('ux1  ','In main dive   ',ux1,size(ux1))
!!! CM call test_min_max('uy1  ','In main dive   ',uy1,size(uy1))
!!! CM call test_min_max('uz1  ','In main dive   ',uz1,size(uz1))

      !POISSON Z-->Z 
      call decomp_2d_poisson_stg(pp3,bcx,bcy,bcz)

!!! CM call test_min_max('pp3  ','In main deco   ',pp3,size(pp3))

      !Z-->Y-->X
      call gradp(px1,py1,pz1,di1,td2,tf2,ta2,tb2,tc2,di2,&
           ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)

      !X PENCILS
      call corgp(ux1,ux2,uy1,uz1,px1,py1,pz1) 

      !does not matter -->output=DIV U=0 (in dv3)
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,dv3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,2)

      call test_speed_min_max(ux1,uy1,uz1)
      call test_density_min_max(rho1)
      if (iscalar==1) call test_scalar_min_max(phi1)

   enddo ! End sub-timesteps

!!$   call STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
!!$        uvmean,uwmean,vwmean,phiphimean,tmean)
   

   if (mod(itime,isave)==0) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
        px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,1)
     
   if (mod(itime,imodulo)==0) then
     call VISU_INSTA(ux1,uy1,uz1,rho1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
          ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
          ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
     call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
          ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
   endif

   if (mod(itime,10)==0) then
     call VISU_INSTB(ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
          ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
          ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
   endif

   ! MMS: compare errors
   CALL eval_error_rho(rho1)
   CALL eval_error_vel(ux1,uy1,uz1)

enddo

t2=MPI_WTIME()-t1
call MPI_ALLREDUCE(t2,t1,1,MPI_REAL8,MPI_SUM, &
                   MPI_COMM_WORLD,code)
if (nrank==0) print *,'time per time_step: ', &
     t1/float(nproc)/(ilast-ifirst+1),' seconds'
if (nrank==0) print *,'simulation with nx*ny*nz=',nx,ny,nz,'mesh nodes'
if (nrank==0) print *,'Mapping p_row*p_col=',p_row,p_col


!call decomp_2d_poisson_finalize
call decomp_2d_finalize
CALL MPI_FINALIZE(code)

end PROGRAM incompact3d
