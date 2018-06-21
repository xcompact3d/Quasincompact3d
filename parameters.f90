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
!
!********************************************************************
!
subroutine parameter()
!
!********************************************************************
  
USE param
USE IBM 
USE variables
USE decomp_2d

implicit none

real(mytype) :: re, theta, cfl,cf2 
integer :: longueur ,impi,j
character :: a*80

#ifdef DOUBLE_PREC 
pi=dacos(-1.d0) 
#else
pi=acos(-1.)
#endif

twopi=2.*pi

1000 format(a,80x) 
1003 format(a,80x)
open(10,file='incompact3d.prm',status='unknown',form='formatted') 
read (10,1000) a 
read (10,1000) a 
read (10,1000) a 
read (10,*) xlx
read (10,*) yly 
read (10,*) zlz 
read (10,*) re
read (10,*) pr
read (10,*) sc
read (10,*) frx
read (10,*) fry
read (10,*) frz
read (10,*) u1 
read (10,*) u2
read (10,*) dens1
read (10,*) dens2
read (10,*) noise 
read (10,*) noise1
read (10,*) dt
read (10,1000) a 
read (10,1000) a 
read (10,1000) a 
read (10,*) nclx 
read (10,*) ncly 
read (10,*) nclz 
read (10,*) itype 
read (10,*) iin
read (10,*) ifirst
read (10,*) ilast
read (10,*) nscheme
read (10,*) istret
read (10,*) beta
read (10,*) ilmn
read (10,*) isolvetemp
read (10,*) imulticomponent
read (10,*) nrhoscheme
read (10,*) ivarcoeff
read (10,*) npoissscheme
read (10,*) tol
read (10,*) iskew
read (10,*) iprops
read (10,*) iscalar
read (10,1000) a 
read (10,1000) a 
read (10,1000) a 
read (10,*) ilit 
read (10,*) isave
read (10,*) imodulo
read (10,1000) a 
read (10,1000) a 
read (10,1000) a 
read (10,*) ivirt
read (10,*) cex 
read (10,*) cey 
read (10,*) cez 
read (10,*) ra 
read (10,1000) a 
close(10) 
if (nrank==0) then
print *,'==========================================================='
print *,'==========================================================='
print *,'==========================================================='
print *,'======================Incompact3d=========================='
print *,'===Copyright (c) 2012 Eric Lamballais and Sylvain Laizet==='
print *,'eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com'
print *,'==========================================================='
print *,'==========================================================='
print *,'==========================================================='
print *,''
print *,''
print *,''
if (itype.eq.1) then
  print *,'Constant flow field'
else if (itype.eq.2) then
  print *,'Channel flow'
else if (itype.eq.3) then
  print *,'Wake flow'
else if (itype.eq.4) then
  print *,'Mixing layer with splitter plate'
else if (itype.eq.5) then
  print *,'Channel flow'
else if (itype.eq.6) then
  print *,'Taylor Green vortices'
else if (itype.eq.7) then
  print *,'Cavity flow'
else if (itype.eq.8) then
  print *,'Flat plate Boundary layer'
else if (itype.eq.9) then
  print *,'Water tank'
else
  print *, "Unknown itype=", itype
  STOP
endif
write(*,1101) nx,ny,nz
write(*,1103) xlx,yly,zlz 
write(*,1102) nclx,ncly,nclz 
write(*,1104) u1,u2 
write(*,1105) dens1,dens2 
write(*,1106) re
write(*,1107) dt
write(*,1112) pr
write(*,1114) frx,fry,frz

if (nscheme.eq.1) then
  print *,'Temporal scheme   : Adams-bashforth 2'
else if (nscheme.eq.2) then
  print *,'Temporal scheme   : Runge-Kutta 3'
else if (nscheme.eq.3) then
  print *,'Temporal scheme   : Runge-Kutta 4'
else if (nscheme.eq.4) then
  print *,"Temporal scheme   : Adams-bashforth 4?"
else
  print *,"Unknown nscheme=", nscheme
  STOP
endif

if (iprops.eq.0) then
  print *, "Variable props    : off"
else
  print *, "Variable props    : on"
endif

if (iscalar.eq.0) then
  print *, 'Passive scalar    : off'
else if (iscalar.eq.1) then
  print *, 'Passive scalar : on'
  write (*,1113) sc
else
  print *, "Unknown iscalar=", iscalar
endif

if (ivirt.eq.0) then
  print *,'Immersed boundary : off'
else if (ivirt.eq.1) then
  print *,'Immersed boundary : on old school'
  write(*,1108) cex,cey,cez
  write(*,1110) ra
else if (ivirt.eq.2) then
  print *,'Immersed boundary : on with Lagrangian Poly'
else
  print *, "Unknown ivirt=", ivirt
  STOP
endif

if (ilmn.ne.0) then
  print *, "Low Mach Number: Enabled"

  if (isolvetemp.eq.0) then
    print *, "Solving for density"
  else
    print *, "Solving for temperature"
  endif
  if (imulticomponent.eq.0) then
    print *, "Multicomponent: DISABLED"
  else
    print *, "Multicomponent: ENABLED"
  endif

  if (ivarcoeff.ne.0) then
    print *, "Var-coeff Poisson: ENABLED"
    print *, "Poisson tolerance: ", tol

    if (npoissscheme.eq.0) then
      print *, "Using rho0 = MIN(rho) in var-coeff Poisson equation"
    else
      print *, "Using rho0 = harmonic_avg(rho) in var-coeff Poisson equation"
    endif
  endif
endif

if (iskew.eq.0) then
  print *, "Adv-MOM: Rotational"
else if (iskew.eq.1) then
  print *, "Adv-MOM: Quasi-skew symmetric"
else if (iskew.eq.2) then
  print *, "Adv-MOM: Skew symmetric"
else
  print *, "Unknown iskew=", iskew
  STOP
endif

1101 format(' Spatial Resolution: (nx,ny,nz)=(',I4,',',I4,',',I4,')')
1102 format(' Boundary condition: (nclx,ncly,nclz)=(',I1,',',I1,',',I1,')')
1103 format(' Domain dimension  : (lx,ly,lz)=(',F6.1,',',F6.1,',',F6.1,')')
1104 format(' High and low speed: u1=',F6.2,' and u2=',F6.2)
1105 format(' High and low density: dens1=',F6.2,' and dens2=',F6.2)
1106 format(' Reynolds number Re: ',F15.8)
1107 format(' Time step dt      : ',F15.8)
1108 format(' Object centred at : (',F6.2,',',F6.2,',',F6.2,')')
1110 format(' Object length     : ',F6.2)
1112 format(' Prandtl number    : ',F6.2)
1113 format(' Schmidt number    : ',F6.2)
1114 format(' Froude number     : (frx,fry,frz)=(',F6.2,',',F6.2,',',F6.2,')')
endif
xnu=1._mytype/re 
   
if (nclx==0) dx=xlx/nx 
if (nclx==1 .or. nclx==2) dx=xlx/(nx-1._mytype) 
if (ncly==0) dy=yly/ny 
if (ncly==1.or.ncly==2) dy   =yly/(ny-1._mytype) 
dx2=dx*dx
dy2=dy*dy
#ifndef TWOD
   if (nclz==0) dz=zlz/nz 
   if (nclz==1.or.nclz==2) dz=zlz/(nz-1._mytype) 
   dz2=dz*dz
#endif

if (istret.eq.0) then
   do j=1,ny
      yp(j)=(j-1._mytype)*dy
      ypi(j)=(j-0.5_mytype)*dy
   enddo
else
   call stretching()
endif


!******************************************************************
!
!**TIME ADVANCE***1=AB2***2=RK3***3=RK4C&K****************** 
!
!******************************************************************

adt(:)=0._mytype ; bdt(:)=0._mytype ; cdt(:)=0._mytype ; gdt(:)=0._mytype
if (nscheme==1) then!AB2
   iadvance_time=1 
   adt(1)=1.5_mytype*dt
   bdt(1)=-0.5_mytype*dt
   gdt(1)=adt(1)+bdt(1)
   gdt(3)=gdt(1)
endif
if (nscheme==2) then !RK3
   iadvance_time=3 
   adt(1)=(8._mytype/15._mytype)*dt
   bdt(1)=0._mytype
   gdt(1)=adt(1)
   adt(2)=(5._mytype/12._mytype)*dt
   bdt(2)=(-17._mytype/60._mytype)*dt
   gdt(2)=adt(2)+bdt(2)
   adt(3)=(3._mytype/4._mytype)*dt
   bdt(3)=(-5._mytype/12._mytype)*dt
   gdt(3)=adt(3)+bdt(3)
endif
if (nscheme==3) then !RK4 Carpenter and Kennedy  
   iadvance_time=5 
   adt(1)=0._mytype
   adt(2)=-0.4178904745_mytype
   adt(3)=-1.192151694643_mytype
   adt(4)=-1.697784692471_mytype
   adt(5)=-1.514183444257_mytype
   bdt(1)=0.1496590219993_mytype
   bdt(2)=0.3792103129999_mytype
   bdt(3)=0.8229550293869_mytype
   bdt(4)=0.6994504559488_mytype
   bdt(5)=0.1530572479681_mytype
   gdt(1)=0.1496590219993_mytype*dt
   gdt(2)=0.220741935365_mytype*dt
   gdt(3)=0.25185480577_mytype*dt
   gdt(4)=0.33602636754_mytype*dt
   gdt(5)=0.041717869325_mytype*dt
endif

if (nscheme==4) then!AB3
   iadvance_time=1
   adt(1)= (23._mytype/12._mytype)*dt
   bdt(1)=-(16._mytype/12._mytype)*dt
   cdt(1)= ( 5._mytype/12._mytype)*dt
   gdt(1)=adt(1)+bdt(1)+cdt(1)
   gdt(3)=gdt(1)
endif

return  
end subroutine parameter



