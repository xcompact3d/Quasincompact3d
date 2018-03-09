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
subroutine schemes()
!
!********************************************************************

USE param
USE derivX
USE derivY
USE derivZ
USE variables
USE var

implicit none

integer  :: i,j,k
real(mytype) :: fpi2

alfa1x= 2._mytype
af1x  =-(5._mytype/2._mytype  )/dx
bf1x  = (   2._mytype  )/dx
cf1x  = (1._mytype/2._mytype  )/dx
df1x  = 0._mytype
alfa2x= 1._mytype/4._mytype
af2x  = (3._mytype/4._mytype  )/dx
alfanx= 2._mytype
afnx  =-(5._mytype/2._mytype  )/dx
bfnx  = (   2._mytype  )/dx
cfnx  = (1._mytype/2._mytype  )/dx
dfnx  = 0._mytype
alfamx= 1._mytype/4._mytype
afmx  = (3._mytype/4._mytype  )/dx
alfaix= 1._mytype/3._mytype
afix  = (7._mytype/9._mytype  )/dx
bfix  = (1._mytype/36._mytype )/dx
alsa1x= 11._mytype
as1x  = (13._mytype    )/dx2
bs1x  =-(27._mytype    )/dx2
cs1x  = (15._mytype    )/dx2
ds1x  =-(1._mytype     )/dx2
alsa2x= 1._mytype/10._mytype
as2x  = (6._mytype/5._mytype  )/dx2
alsa3x= 2._mytype/11._mytype
as3x  = (12._mytype/11._mytype)/dx2
bs3x  = (3._mytype/44._mytype )/dx2
alsanx= 11._mytype
asnx  = (13._mytype    )/dx2
bsnx  =-(27._mytype    )/dx2
csnx  = (15._mytype    )/dx2
dsnx  =-(1._mytype     )/dx2
alsamx= 1._mytype/10._mytype
asmx  = (6._mytype/5._mytype  )/dx2
alsatx= 2._mytype/11._mytype
astx  = (12._mytype/11._mytype)/dx2
bstx  = (3._mytype/44._mytype )/dx2

!alsaix= 2./11.
!asix  = (12./11.)/dx2
!bsix  = (3./44. )/dx2
!csix  = 0.
!NUMERICAL DISSIPATION (see publications for help)
fpi2=4._mytype        ! Uncomment this line to add some small scale numerical dissipation
!fpi2=(48./7)/(pi*pi) !
alsaix=(45._mytype*fpi2*pi*pi-272._mytype)/(2._mytype*(45._mytype*fpi2*pi*pi-208._mytype))
asix  =((6._mytype-9._mytype*alsaix)/4._mytype)/dx2
bsix  =((-3._mytype+24._mytype*alsaix)/5._mytype)/(4._mytype*dx2)
csix  =((2._mytype-11._mytype*alsaix)/20._mytype)/(9._mytype*dx2)
!      stop
!if (nrank==0) then
!      write(*,*) '=== derxx ==='
!      write(*,*) alsaix
!      write(*,*) asix*dx2
!      write(*,*) bsix*4*dx2
!      write(*,*) csix*9*dx2
!      write(*,*) '============='

alfa1y= 2._mytype
af1y  =-(5._mytype/2._mytype  )/dy
bf1y  = (   2._mytype  )/dy
cf1y  = (1._mytype/2._mytype  )/dy
df1y  = 0._mytype
alfa2y= 1._mytype/4._mytype
af2y  = (3._mytype/4._mytype  )/dy
alfany= 2._mytype
afny  =-(5._mytype/2._mytype  )/dy
bfny  = (   2._mytype  )/dy
cfny  = (1._mytype/2._mytype  )/dy
dfny  = 0._mytype
alfamy= 1._mytype/4._mytype
afmy  = (3._mytype/4._mytype  )/dy
alfajy= 1._mytype/3._mytype
afjy  = (7._mytype/9._mytype  )/dy
bfjy  = (1._mytype/36._mytype )/dy
alsa1y= 11._mytype
as1y  = (13._mytype    )/dy2
bs1y  =-(27._mytype    )/dy2
cs1y  = (15._mytype    )/dy2
ds1y  =-(1._mytype     )/dy2
alsa2y= 1._mytype/10._mytype
as2y  = (6._mytype/5._mytype  )/dy2
alsa3y= 2._mytype/11._mytype
as3y  = (12._mytype/11._mytype)/dy2
bs3y  = (3._mytype/44._mytype )/dy2
alsany= 11._mytype
asny  = (13._mytype    )/dy2
bsny  =-(27._mytype    )/dy2
csny  = (15._mytype    )/dy2
dsny  =-(1._mytype     )/dy2
alsamy= 1._mytype/10._mytype
asmy  = (6._mytype/5._mytype  )/dy2
alsaty= 2._mytype/11._mytype
asty  = (12._mytype/11._mytype)/dy2
bsty  = (3._mytype/44._mytype )/dy2

!alsajy= 2./11.
!asjy  = (12./11.)/dy2
!bsjy  = (3./44. )/dy2
!csjy   = 0.

alsajy=(45._mytype*fpi2*pi*pi-272._mytype)/(2._mytype*(45._mytype*fpi2*pi*pi-208._mytype))
asjy  =((6._mytype-9._mytype*alsajy)/4._mytype)/dy2
bsjy  =((-3._mytype+24._mytype*alsajy)/5._mytype)/(4._mytype*dy2)
csjy  =((2._mytype-11._mytype*alsajy)/20._mytype)/(9._mytype*dy2)
!if (nrank==0) then
!      write(*,*) '=== deryy ==='
!      write(*,*) alsajy
!      write(*,*) asjy*dy2
!      write(*,*) bsjy*4*dy2
!      write(*,*) csjy*9*dy2
!      write(*,*) '============='
!endif
alcaix6=9._mytype/62._mytype 
acix6=(63._mytype/62._mytype)/dx
bcix6=(17._mytype/62._mytype)/3._mytype/dx 
cfx6(1)=alcaix6 
cfx6(2)=alcaix6 
cfx6(nxm-2)=alcaix6 
cfx6(nxm-1)=alcaix6 
cfx6(nxm)=0._mytype 
if (nclx==0) ccx6(1)=2._mytype  
if (nclx==1) ccx6(1)=1._mytype + alcaix6 
if (nclx==2) ccx6(1)=1._mytype + alcaix6 
ccx6(2)=1._mytype 
ccx6(nxm-2)=1._mytype 
ccx6(nxm-1)=1._mytype 
if (nclx==0) ccx6(nxm)=1._mytype + alcaix6*alcaix6  
if (nclx==1) ccx6(nxm)=1._mytype + alcaix6 
if (nclx==2) ccx6(nxm)=1._mytype + alcaix6
cbx6(1)=alcaix6 
cbx6(2)=alcaix6 
cbx6(nxm-2)=alcaix6 
cbx6(nxm-1)=alcaix6 
cbx6(nxm)=0._mytype 
do i=3,nxm-3 
   cfx6(i)=alcaix6 
   ccx6(i)=1._mytype 
   cbx6(i)=alcaix6 
enddo

cfi6(1)=alcaix6 + alcaix6 
cfi6(2)=alcaix6 
cfi6(nx-2)=alcaix6 
cfi6(nx-1)=alcaix6 
cfi6(nx)=0._mytype 
cci6(1)=1._mytype 
cci6(2)=1._mytype 
cci6(nx-2)=1._mytype 
cci6(nx-1)=1._mytype 
cci6(nx)=1._mytype 
cbi6(1)=alcaix6 
cbi6(2)=alcaix6 
cbi6(nx-2)=alcaix6 
cbi6(nx-1)=alcaix6 + alcaix6 
cbi6(nx)=0._mytype 
do i=3,nx-3 
   cfi6(i)=alcaix6 
   cci6(i)=1._mytype 
   cbi6(i)=alcaix6 
enddo

ailcaix6=0.49_mytype!3._mytype/10._mytype 
aicix6=1._mytype/128._mytype*(75._mytype+70._mytype*ailcaix6)
bicix6=1._mytype/256._mytype*(126._mytype*ailcaix6-25._mytype)
cicix6=1._mytype/256._mytype*(-10._mytype*ailcaix6+3._mytype)
if (nrank==0) print *,'New coef Inter X',aicix6,bicix6,cicix6
cifx6(1)=ailcaix6 
cifx6(2)=ailcaix6 
cifx6(nxm-2)=ailcaix6 
cifx6(nxm-1)=ailcaix6 
cifx6(nxm)=0. 
if (nclx==0) cicx6(1)=2._mytype  
if (nclx==1) cicx6(1)=1._mytype + ailcaix6 
if (nclx==2) cicx6(1)=1._mytype + ailcaix6 
cicx6(2)=1._mytype 
cicx6(nxm-2)=1._mytype 
cicx6(nxm-1)=1._mytype 
if (nclx==0) cicx6(nxm)=1._mytype + ailcaix6*ailcaix6
if (nclx==1) cicx6(nxm)=1._mytype + ailcaix6 
if (nclx==2) cicx6(nxm)=1._mytype + ailcaix6
cibx6(1)=ailcaix6 
cibx6(2)=ailcaix6 
cibx6(nxm-2)=ailcaix6 
cibx6(nxm-1)=ailcaix6 
cibx6(nxm)=0._mytype 
do i=3,nxm-3 
   cifx6(i)=ailcaix6 
   cicx6(i)=1._mytype 
   cibx6(i)=ailcaix6 
enddo
cifi6(1)=ailcaix6 + ailcaix6 
cifi6(2)=ailcaix6 
cifi6(nx-2)=ailcaix6 
cifi6(nx-1)=ailcaix6 
cifi6(nx)=0._mytype 
cici6(1)=1._mytype 
cici6(2)=1._mytype 
cici6(nx-2)=1._mytype 
cici6(nx-1)=1._mytype 
cici6(nx)=1._mytype 
cibi6(1)=ailcaix6 
cibi6(2)=ailcaix6 
cibi6(nx-2)=ailcaix6 
cibi6(nx-1)=ailcaix6 + ailcaix6 
cibi6(nx)=0._mytype 
do i=3,nx-3 
   cifi6(i)=ailcaix6 
   cici6(i)=1._mytype 
   cibi6(i)=ailcaix6 
enddo

alcaiy6=9._mytype/62._mytype 
aciy6=(63._mytype/62._mytype)/dy 
bciy6=(17._mytype/62._mytype)/3._mytype/dy 
cfy6(1)=alcaiy6 
cfy6(2)=alcaiy6 
cfy6(nym-2)=alcaiy6 
cfy6(nym-1)=alcaiy6 
cfy6(nym)=0._mytype 
if (ncly==0) ccy6(1)=2._mytype  
if (ncly==1) ccy6(1)=1._mytype + alcaiy6 
if (ncly==2) ccy6(1)=1._mytype + alcaiy6
ccy6(2)=1._mytype 
ccy6(nym-2)=1._mytype 
ccy6(nym-1)=1._mytype 
if (ncly==0) ccy6(nym)=1._mytype + alcaiy6*alcaiy6  
if (ncly==1) ccy6(nym)=1._mytype + alcaiy6 
if (ncly==2) ccy6(nym)=1._mytype + alcaiy6
cby6(1)=alcaiy6 
cby6(2)=alcaiy6 
cby6(nym-2)=alcaiy6 
cby6(nym-1)=alcaiy6 
cby6(nym)=0._mytype 
do j=3,nym-3 
   cfy6(j)=alcaiy6 
   ccy6(j)=1._mytype 
   cby6(j)=alcaiy6 
enddo
cfi6y(1)=alcaiy6 + alcaiy6 
cfi6y(2)=alcaiy6 
cfi6y(ny-2)=alcaiy6 
cfi6y(ny-1)=alcaiy6 
cfi6y(ny)=0._mytype 
cci6y(1)=1._mytype 
cci6y(2)=1._mytype 
cci6y(ny-2)=1._mytype 
cci6y(ny-1)=1._mytype 
cci6y(ny)=1._mytype 
cbi6y(1)=alcaiy6 
cbi6y(2)=alcaiy6 
cbi6y(ny-2)=alcaiy6 
cbi6y(ny-1)=alcaiy6 + alcaiy6 
cbi6y(ny)=0._mytype 
do j=3,ny-3 
   cfi6y(j)=alcaiy6 
   cci6y(j)=1._mytype 
   cbi6y(j)=alcaiy6 
enddo
ailcaiy6=0.49_mytype!3._mytype/10._mytype
aiciy6=1._mytype/128._mytype*(75._mytype+70._mytype*ailcaiy6)
biciy6=1._mytype/256._mytype*(126._mytype*ailcaiy6-25._mytype)
ciciy6=1._mytype/256._mytype*(-10._mytype*ailcaiy6+3._mytype)
if (nrank==0) print *,'New coef Inter Y',aiciy6,biciy6,ciciy6
cify6(1)=ailcaiy6 
cify6(2)=ailcaiy6 
cify6(nym-2)=ailcaiy6 
cify6(nym-1)=ailcaiy6 
cify6(nym)=0._mytype 
if (ncly==0) cicy6(1)=2._mytype  
if (ncly==1) cicy6(1)=1._mytype + ailcaiy6 
if (ncly==2) cicy6(1)=1._mytype + ailcaiy6 
cicy6(2)=1._mytype 
cicy6(nym-2)=1._mytype 
cicy6(nym-1)=1._mytype 
if (ncly==0) cicy6(nym)=1._mytype + ailcaiy6*ailcaiy6
if (ncly==1) cicy6(nym)=1._mytype + ailcaiy6 
if (ncly==2) cicy6(nym)=1._mytype + ailcaiy6
ciby6(1)=ailcaiy6 
ciby6(2)=ailcaiy6 
ciby6(nym-2)=ailcaiy6 
ciby6(nym-1)=ailcaiy6 
ciby6(nym)=0._mytype 
do j=3,nym-3 
   cify6(j)=ailcaiy6 
   cicy6(j)=1._mytype 
   ciby6(j)=ailcaiy6 
enddo
cifi6y(1)=ailcaiy6 + ailcaiy6 
cifi6y(2)=ailcaiy6 
cifi6y(ny-2)=ailcaiy6 
cifi6y(ny-1)=ailcaiy6 
cifi6y(ny)=0._mytype 
cici6y(1)=1._mytype 
cici6y(2)=1._mytype 
cici6y(ny-2)=1._mytype 
cici6y(ny-1)=1._mytype 
cici6y(ny)=1._mytype 
cibi6y(1)=ailcaiy6 
cibi6y(2)=ailcaiy6 
cibi6y(ny-2)=ailcaiy6 
cibi6y(ny-1)=ailcaiy6 + ailcaiy6 
cibi6y(ny)=0._mytype 
do j=3,ny-3 
   cifi6y(j)=ailcaiy6 
   cici6y(j)=1._mytype 
   cibi6y(j)=ailcaiy6 
enddo

#ifndef TWOD
   alcaiz6=9._mytype/62._mytype 
   aciz6=(63._mytype/62._mytype)/dz 
   bciz6=(17._mytype/62._mytype)/3._mytype/dz 
   cfz6(1)=alcaiz6 
   cfz6(2)=alcaiz6 
   cfz6(nzm-2)=alcaiz6 
   cfz6(nzm-1)=alcaiz6 
   cfz6(nzm)=0._mytype 
   if (nclz==0) ccz6(1)=2._mytype  
   if (nclz==1) ccz6(1)=1._mytype + alcaiz6 
   if (nclz==2) ccz6(1)=1._mytype + alcaiz6
   ccz6(2)=1._mytype 
   ccz6(nzm-2)=1._mytype 
   ccz6(nzm-1)=1._mytype 
   if (nclz==0) ccz6(nzm)=1._mytype + alcaiz6*alcaiz6  
   if (nclz==1) ccz6(nzm)=1._mytype + alcaiz6 
   if (nclz==2) ccz6(nzm)=1._mytype + alcaiz6
   cbz6(1)=alcaiz6 
   cbz6(2)=alcaiz6 
   cbz6(nzm-2)=alcaiz6 
   cbz6(nzm-1)=alcaiz6 
   cbz6(nzm)=0._mytype 
   do k=3,nzm-3 
      cfz6(k)=alcaiz6 
      ccz6(k)=1. 
      cbz6(k)=alcaiz6 
   enddo
   cfi6z(1)=alcaiz6 + alcaiz6 
   cfi6z(2)=alcaiz6 
   cfi6z(nz-2)=alcaiz6 
   cfi6z(nz-1)=alcaiz6 
   cfi6z(nz)=0._mytype 
   cci6z(1)=1._mytype 
   cci6z(2)=1._mytype 
   cci6z(nz-2)=1._mytype 
   cci6z(nz-1)=1._mytype 
   cci6z(nz)=1._mytype 
   cbi6z(1)=alcaiz6 
   cbi6z(2)=alcaiz6 
   cbi6z(nz-2)=alcaiz6 
   cbi6z(nz-1)=alcaiz6 + alcaiz6 
   cbi6z(nz)=0._mytype 
   do k=3,nz-3 
      cfi6z(k)=alcaiz6 
      cci6z(k)=1._mytype 
      cbi6z(k)=alcaiz6 
   enddo
   ailcaiz6=0.49_mytype!3._mytype/10._mytype
   aiciz6=1._mytype/128._mytype*(75._mytype+70._mytype*ailcaiz6)
   biciz6=1._mytype/256._mytype*(126._mytype*ailcaiz6-25._mytype)
   ciciz6=1._mytype/256._mytype*(-10._mytype*ailcaiz6+3._mytype)
   if (nrank==0) print *,'New coef Inter Z',aiciz6,biciz6,ciciz6
   cifz6(1)=ailcaiz6 
   cifz6(2)=ailcaiz6 
   cifz6(nzm-2)=ailcaiz6 
   cifz6(nzm-1)=ailcaiz6 
   cifz6(nzm)=0._mytype 
   if (nclz==0) cicz6(1)=2._mytype  
   if (nclz==1) cicz6(1)=1._mytype + ailcaiz6 
   if (nclz==2) cicz6(1)=1._mytype + ailcaiz6 
   cicz6(2)=1._mytype 
   cicz6(nzm-2)=1._mytype 
   cicz6(nzm-1)=1._mytype 
   if (nclz==0) cicz6(nzm)=1._mytype + ailcaiz6*ailcaiz6
   if (nclz==1) cicz6(nzm)=1._mytype + ailcaiz6 
   if (nclz==2) cicz6(nzm)=1._mytype + ailcaiz6
   cibz6(1)=ailcaiz6 
   cibz6(2)=ailcaiz6 
   cibz6(nzm-2)=ailcaiz6 
   cibz6(nzm-1)=ailcaiz6 
   cibz6(nzm)=0._mytype 
   do k=3,nzm-3 
      cifz6(k)=ailcaiz6 
      cicz6(k)=1._mytype 
      cibz6(k)=ailcaiz6 
   enddo
   cifi6z(1)=ailcaiz6 + ailcaiz6 
   cifi6z(2)=ailcaiz6 
   cifi6z(nz-2)=ailcaiz6 
   cifi6z(nz-1)=ailcaiz6 
   cifi6z(nz)=0._mytype 
   cici6z(1)=1._mytype 
   cici6z(2)=1._mytype 
   cici6z(nz-2)=1._mytype 
   cici6z(nz-1)=1._mytype 
   cici6z(nz)=1._mytype 
   cibi6z(1)=ailcaiz6 
   cibi6z(2)=ailcaiz6 
   cibi6z(nz-2)=ailcaiz6 
   cibi6z(nz-1)=ailcaiz6 + ailcaiz6 
   cibi6z(nz)=0._mytype 
   do k=3,nz-3 
      cifi6z(k)=ailcaiz6 
      cici6z(k)=1._mytype 
      cibi6z(k)=ailcaiz6 
   enddo
   alfa1z= 2._mytype
   af1z  =-(5._mytype/2._mytype  )/dz
   bf1z  = (   2._mytype  )/dz
   cf1z  = (1._mytype/2._mytype  )/dz
   df1z  = 0._mytype
   alfa2z= 1._mytype/4._mytype
   af2z  = (3._mytype/4._mytype  )/dz
   alfanz= 2._mytype
   afnz  =-(5._mytype/2._mytype  )/dz
   bfnz  = (   2._mytype  )/dz
   cfnz  = (1._mytype/2._mytype  )/dz
   dfnz  = 0._mytype
   alfamz= 1._mytype/4._mytype
   afmz  = (3._mytype/4._mytype  )/dz
   alfakz= 1._mytype/3._mytype
   afkz  = (7._mytype/9._mytype  )/dz
   bfkz  = (1._mytype/36._mytype )/dz
   alsa1z= 11._mytype
   as1z  = (13._mytype    )/dz2
   bs1z  =-(27._mytype    )/dz2
   cs1z  = (15._mytype    )/dz2
   ds1z  =-(1._mytype     )/dz2
   alsa2z= 1._mytype/10._mytype
   as2z  = (6._mytype/5._mytype  )/dz2
   alsa3z= 2._mytype/11._mytype
   as3z  = (12._mytype/11._mytype)/dz2
   bs3z  = (3._mytype/44._mytype )/dz2
   alsanz= 11._mytype
   asnz  = (13._mytype    )/dz2
   bsnz  =-(27._mytype    )/dz2
   csnz  = (15._mytype    )/dz2
   dsnz  =-(1._mytype     )/dz2
   alsamz= 1._mytype/10._mytype
   asmz  = (6._mytype/5._mytype  )/dz2
   alsatz= 2._mytype/11._mytype
   astz  = (12._mytype/11._mytype)/dz2
   bstz  = (3._mytype/44._mytype )/dz2

!   alsakz= 2./11.
!   askz  = (12./11.)/dz2
!   bskz  = (3./44. )/dz2
!   cskz  = 0.
         alsakz=(45._mytype*fpi2*pi*pi-272._mytype)/(2._mytype*(45._mytype*fpi2*pi*pi-208._mytype))
         askz  =((6._mytype-9._mytype*alsakz)/4._mytype)/dz2
         bskz  =((-3._mytype+24._mytype*alsakz)/5._mytype)/(4._mytype*dz2)
         cskz  =((2._mytype-11._mytype*alsakz)/20._mytype)/(9._mytype*dz2)
!if (nrank==0) then
!         write(*,*) '=== derzz ==='
!         write(*,*) alsakz
!         write(*,*) askz*dz2
!         write(*,*) bskz*4*dz2
!         write(*,*) cskz*9*dz2
!         write(*,*) '============='
!endif
#endif

if (nclx.eq.0) then
   ffx(1)   =alfaix
   ffx(2)   =alfaix
   ffx(nx-2)=alfaix
   ffx(nx-1)=alfaix
   ffx(nx)  =0._mytype
   fcx(1)   =2._mytype
   fcx(2)   =1._mytype
   fcx(nx-2)=1._mytype
   fcx(nx-1)=1._mytype
   fcx(nx  )=1._mytype+alfaix*alfaix
   fbx(1)   =alfaix
   fbx(2)   =alfaix
   fbx(nx-2)=alfaix
   fbx(nx-1)=alfaix
   fbx(nx  )=0._mytype
   do i=3,nx-3
      ffx(i)=alfaix
      fcx(i)=1._mytype
      fbx(i)=alfaix
   enddo
endif

if (nclx.eq.1) then
   ffx(1)   =alfaix+alfaix
   ffx(2)   =alfaix
   ffx(nx-2)=alfaix
   ffx(nx-1)=alfaix
   ffx(nx)  =0._mytype
   fcx(1)   =1._mytype
   fcx(2)   =1._mytype
   fcx(nx-2)=1._mytype
   fcx(nx-1)=1._mytype
   fcx(nx  )=1._mytype
   fbx(1)   =alfaix 
   fbx(2)   =alfaix
   fbx(nx-2)=alfaix
   fbx(nx-1)=alfaix+alfaix
   fbx(nx  )=0._mytype
   do i=3,nx-3
      ffx(i)=alfaix
      fcx(i)=1._mytype
      fbx(i)=alfaix
   enddo
endif

if (nclx.eq.2) then
   ffx(1)   =alfa1x
   ffx(2)   =alfa2x
   ffx(nx-2)=alfaix
   ffx(nx-1)=alfamx
   ffx(nx)  =0._mytype
   fcx(1)   =1._mytype
   fcx(2)   =1._mytype
   fcx(nx-2)=1._mytype
   fcx(nx-1)=1._mytype
   fcx(nx  )=1._mytype
   fbx(1)   =alfa2x 
   fbx(2)   =alfaix
   fbx(nx-2)=alfamx
   fbx(nx-1)=alfanx
   fbx(nx  )=0._mytype
   do i=3,nx-3
      ffx(i)=alfaix
      fcx(i)=1._mytype
      fbx(i)=alfaix
   enddo
endif

if (ncly.eq.0) then
   ffy(1)   =alfajy
   ffy(2)   =alfajy
   ffy(ny-2)=alfajy
   ffy(ny-1)=alfajy
   ffy(ny)  =0._mytype
   fcy(1)   =2._mytype
   fcy(2)   =1._mytype
   fcy(ny-2)=1._mytype
   fcy(ny-1)=1._mytype
   fcy(ny  )=1._mytype+alfajy*alfajy
   fby(1)   =alfajy
   fby(2)   =alfajy
   fby(ny-2)=alfajy
   fby(ny-1)=alfajy
   fby(ny  )=0._mytype
   do j=3,ny-3
      ffy(j)=alfajy
      fcy(j)=1._mytype
      fby(j)=alfajy
   enddo
endif

if (ncly.eq.1) then
   ffy(1)   =alfajy+alfajy
   ffy(2)   =alfajy
   ffy(ny-2)=alfajy
   ffy(ny-1)=alfajy
   ffy(ny)  =0._mytype
   fcy(1)   =1._mytype
   fcy(2)   =1._mytype
   fcy(ny-2)=1._mytype
   fcy(ny-1)=1._mytype
   fcy(ny  )=1._mytype
   fby(1)   =alfajy 
   fby(2)   =alfajy
   fby(ny-2)=alfajy
   fby(ny-1)=alfajy+alfajy
   fby(ny  )=0._mytype
   do j=3,ny-3
      ffy(j)=alfajy
      fcy(j)=1._mytype
      fby(j)=alfajy
   enddo
endif

if (ncly.eq.2) then
   ffy(1)   =alfa1y
   ffy(2)   =alfa2y
   ffy(ny-2)=alfajy
   ffy(ny-1)=alfamy
   ffy(ny)  =0._mytype
   fcy(1)   =1._mytype
   fcy(2)   =1._mytype
   fcy(ny-2)=1._mytype
   fcy(ny-1)=1._mytype
   fcy(ny  )=1._mytype
   fby(1)   =alfa2y 
   fby(2)   =alfajy
   fby(ny-2)=alfamy
   fby(ny-1)=alfany
   fby(ny  )=0._mytype
   do j=3,ny-3
      ffy(j)=alfajy
      fcy(j)=1._mytype
      fby(j)=alfajy
   enddo
endif

#ifndef TWOD
   if (nclz.eq.0) then
      ffz(1)   =alfakz
      ffz(2)   =alfakz
      ffz(nz-2)=alfakz
      ffz(nz-1)=alfakz
      ffz(nz)  =0._mytype
      fcz(1)   =2._mytype
      fcz(2)   =1._mytype
      fcz(nz-2)=1._mytype
      fcz(nz-1)=1._mytype
      fcz(nz  )=1._mytype+alfakz*alfakz
      fbz(1)   =alfakz
      fbz(2)   =alfakz
      fbz(nz-2)=alfakz
      fbz(nz-1)=alfakz
      fbz(nz  )=0._mytype
      do k=3,nz-3
         ffz(k)=alfakz
         fcz(k)=1._mytype
         fbz(k)=alfakz
      enddo
   endif

   if (nclz.eq.1) then
      ffz(1)   =alfakz+alfakz
      ffz(2)   =alfakz
      ffz(nz-2)=alfakz
      ffz(nz-1)=alfakz
      ffz(nz)  =0._mytype
      fcz(1)   =1._mytype
      fcz(2)   =1._mytype
      fcz(nz-2)=1._mytype
      fcz(nz-1)=1._mytype
      fcz(nz  )=1._mytype
      fbz(1)   =alfakz 
      fbz(2)   =alfakz
      fbz(nz-2)=alfakz
      fbz(nz-1)=alfakz+alfakz
      fbz(nz  )=0._mytype
      do k=3,nz-3
         ffz(k)=alfakz
         fcz(k)=1._mytype
         fbz(k)=alfakz
      enddo
   endif

   if (nclz.eq.2) then
      ffz(1)   =alfa1z
      ffz(2)   =alfa2z
      ffz(nz-2)=alfakz
      ffz(nz-1)=alfamz
      ffz(nz)  =0._mytype
      fcz(1)   =1._mytype
      fcz(2)   =1._mytype
      fcz(nz-2)=1._mytype
      fcz(nz-1)=1._mytype
      fcz(nz  )=1._mytype
      fbz(1)   =alfa2z 
      fbz(2)   =alfakz
      fbz(nz-2)=alfamz
      fbz(nz-1)=alfanz
      fbz(nz  )=0._mytype
      do k=3,nz-3
         ffz(k)=alfakz
         fcz(k)=1._mytype
         fbz(k)=alfakz
      enddo
   endif
#endif

if (nclx.eq.0) then
   sfx(1)   =alsaix
   sfx(2)   =alsaix
   sfx(nx-2)=alsaix
   sfx(nx-1)=alsaix
   sfx(nx)  =0._mytype
   scx(1)   =2._mytype
   scx(2)   =1._mytype
   scx(nx-2)=1._mytype
   scx(nx-1)=1._mytype
   scx(nx  )=1._mytype+alsaix*alsaix
   sbx(1)   =alsaix
   sbx(2)   =alsaix
   sbx(nx-2)=alsaix
   sbx(nx-1)=alsaix
   sbx(nx  )=0._mytype
   do i=3,nx-3
      sfx(i)=alsaix
      scx(i)=1._mytype
      sbx(i)=alsaix
   enddo
endif

if (nclx.eq.1) then
   sfx(1)   =alsaix+alsaix
   sfx(2)   =alsaix
   sfx(nx-2)=alsaix
   sfx(nx-1)=alsaix
   sfx(nx)  =0._mytype
   scx(1)   =1._mytype
   scx(2)   =1._mytype
   scx(nx-2)=1._mytype
   scx(nx-1)=1._mytype
   scx(nx  )=1._mytype
   sbx(1)   =alsaix
   sbx(2)   =alsaix
   sbx(nx-2)=alsaix
   sbx(nx-1)=alsaix+alsaix
   sbx(nx  )=0._mytype
   do i=3,nx-3
      sfx(i)=alsaix
      scx(i)=1._mytype
      sbx(i)=alsaix
   enddo
endif

if (nclx.eq.2) then
   sfx(1)   =alsa1x
   sfx(2)   =alsa2x
   sfx(3)   =alsa3x
   sfx(nx-3)=alsaix
   sfx(nx-2)=alsatx
   sfx(nx-1)=alsamx
   sfx(nx)  =0._mytype
   scx(1)   =1._mytype
   scx(2)   =1._mytype
   scx(3)   =1._mytype
   scx(nx-3)=1._mytype
   scx(nx-2)=1._mytype
   scx(nx-1)=1._mytype
   scx(nx  )=1._mytype
   sbx(1)   =alsa2x 
   sbx(2)   =alsa3x
   sbx(3)   =alsaix
   sbx(nx-3)=alsatx
   sbx(nx-2)=alsamx
   sbx(nx-1)=alsanx
   sbx(nx  )=0._mytype
   do i=4,nx-4
      sfx(i)=alsaix
      scx(i)=1._mytype
      sbx(i)=alsaix
   enddo
endif

if (ncly.eq.0) then
   sfy(1)   =alsajy
   sfy(2)   =alsajy
   sfy(ny-2)=alsajy
   sfy(ny-1)=alsajy
   sfy(ny)  =0._mytype
   scy(1)   =2._mytype
   scy(2)   =1._mytype
   scy(ny-2)=1._mytype
   scy(ny-1)=1._mytype
   scy(ny  )=1._mytype+alsajy*alsajy
   sby(1)   =alsajy
   sby(2)   =alsajy
   sby(ny-2)=alsajy
   sby(ny-1)=alsajy
   sby(ny  )=0._mytype
   do j=3,ny-3
      sfy(j)=alsajy
      scy(j)=1._mytype
      sby(j)=alsajy
   enddo
endif

if (ncly.eq.1) then
   sfy(1)   =alsajy+alsajy
   sfy(2)   =alsajy
   sfy(ny-2)=alsajy
   sfy(ny-1)=alsajy
   sfy(ny)  =0._mytype
   scy(1)   =1._mytype
   scy(2)   =1._mytype
   scy(ny-2)=1._mytype
   scy(ny-1)=1._mytype
   scy(ny  )=1._mytype
   sby(1)   =alsajy 
   sby(2)   =alsajy
   sby(ny-2)=alsajy
   sby(ny-1)=alsajy+alsajy
   sby(ny  )=0._mytype
   do j=3,ny-3
      sfy(j)=alsajy
      scy(j)=1._mytype
      sby(j)=alsajy
   enddo
endif

if (ncly.eq.2) then
   sfy(1)   =alsa1y
   sfy(2)   =alsa2y
   sfy(3)   =alsa3y
   sfy(ny-3)=alsajy
   sfy(ny-2)=alsaty
   sfy(ny-1)=alsamy
   sfy(ny)  =0._mytype
   scy(1)   =1._mytype
   scy(2)   =1._mytype
   scy(3)   =1._mytype
   scy(ny-3)=1._mytype
   scy(ny-2)=1._mytype
   scy(ny-1)=1._mytype
   scy(ny  )=1._mytype
   sby(1)   =alsa2y 
   sby(2)   =alsa3y
   sby(3)   =alsajy
   sby(ny-3)=alsaty
   sby(ny-2)=alsamy
   sby(ny-1)=alsany
   sby(ny  )=0._mytype
   do j=4,ny-4
      sfy(j)=alsajy
      scy(j)=1._mytype
      sby(j)=alsajy
   enddo
endif

#ifndef TWOD
   if (nclz.eq.0) then
      sfz(1)   =alsakz
      sfz(2)   =alsakz
      sfz(nz-2)=alsakz
      sfz(nz-1)=alsakz
      sfz(nz)  =0._mytype
      scz(1)   =2._mytype
      scz(2)   =1._mytype
      scz(nz-2)=1._mytype
      scz(nz-1)=1._mytype
      scz(nz  )=1._mytype+alsakz*alsakz
      sbz(1)   =alsakz
      sbz(2)   =alsakz
      sbz(nz-2)=alsakz
      sbz(nz-1)=alsakz
      sbz(nz  )=0._mytype
      do k=3,nz-3
         sfz(k)=alsakz
         scz(k)=1._mytype
         sbz(k)=alsakz
      enddo
   endif

   if (nclz.eq.1) then
      sfz(1)   =alsakz+alsakz
      sfz(2)   =alsakz
      sfz(nz-2)=alsakz
      sfz(nz-1)=alsakz
      sfz(nz)  =0._mytype
      scz(1)   =1._mytype
      scz(2)   =1._mytype
      scz(nz-2)=1._mytype
      scz(nz-1)=1._mytype
      scz(nz  )=1._mytype
      sbz(1)   =alsakz 
      sbz(2)   =alsakz
      sbz(nz-2)=alsakz
      sbz(nz-1)=alsakz+alsakz
      sbz(nz  )=0._mytype
      do k=3,nz-3
         sfz(k)=alsakz
         scz(k)=1._mytype
         sbz(k)=alsakz
      enddo
   endif
   
   if (nclz.eq.2) then
      sfz(1)   =alsa1z
      sfz(2)   =alsa2z
      sfz(3)   =alsa3z
      sfz(nz-3)=alsakz
      sfz(nz-2)=alsatz
      sfz(nz-1)=alsamz
      sfz(nz)  =0._mytype
      scz(1)   =1._mytype
      scz(2)   =1._mytype
      scz(3)   =1._mytype
      scz(nz-3)=1._mytype
      scz(nz-2)=1._mytype
      scz(nz-1)=1._mytype
      scz(nz  )=1._mytype
      sbz(1)   =alsa2z 
      sbz(2)   =alsa3z
      sbz(3)   =alsakz
      sbz(nz-3)=alsatz
      sbz(nz-2)=alsamz
      sbz(nz-1)=alsanz
      sbz(nz  )=0._mytype
      do k=4,nz-4
         sfz(k)=alsakz
         scz(k)=1._mytype
         sbz(k)=alsakz
      enddo
   endif
#endif

do i=1,nx
   ffxp(i)=ffx(i)
   sfxp(i)=sfx(i)
enddo
do i=1,nxm   
   cfxp6(i)=cfx6(i)
   cifxp6(i)=cifx6(i)
enddo
do i=1,nx
   cifip6(i)=cifi6(i)
   cfip6(i)=cfi6(i)
enddo
do j=1,ny
   ffyp(j)=ffy(j)
   sfyp(j)=sfy(j)
enddo
do j=1,nym
   cfyp6(j)=cfy6(j)
   cifyp6(j)=cify6(j)
enddo
do j=1,ny
   cifip6y(j)=cifi6y(j)
   cfip6y(j)=cfi6y(j)
enddo

#ifndef TWOD
   do k=1,nz
      ffzp(k)=ffz(k)
      sfzp(k)=sfz(k)
   enddo
   do k=1,nzm
      cfzp6(k)=cfz6(k)
      cifzp6(k)=cifz6(k)
   enddo
   do k=1,nz
      cifip6z(k)=cifi6z(k)
      cfip6z(k)=cfi6z(k)
   enddo
#endif

if (nclx.eq.1) then
   ffxp(1)=0._mytype
   sfx (1)=0._mytype
endif
if (ncly.eq.1) then
   ffyp(1)=0._mytype
   sfy (1)=0._mytype
endif
cfxp6(1)=0._mytype
cfip6(1)=0._mytype
cfyp6(1)=0._mytype
cfip6y(1)=0._mytype

#ifndef TWOD
   if (nclz.eq.1) then
      ffzp(1)=0._mytype
      sfz (1)=0._mytype
   endif
   cfzp6(1)=0._mytype
   cfip6z(1)=0._mytype
#endif

  
call prepare (fbx,fcx,ffx ,fsx ,fwx ,nx)
call prepare (fbx,fcx,ffxp,fsxp,fwxp,nx)
call prepare (fby,fcy,ffy ,fsy ,fwy ,ny)
call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
call prepare (cbx6,ccx6,cfx6 ,csx6 ,cwx6 ,nxm)

!!! CM call test_min_max('cbx6 ','In schemes     ',cbx6,size(cbx6))
!!! CM call test_min_max('ccx6 ','In schemes     ',ccx6,size(ccx6))
!!! CM call test_min_max('cfx6 ','In schemes     ',cfx6,size(cfx6))
!!! CM call test_min_max('csx6 ','In schemes     ',csx6,size(csx6))
!!! CM call test_min_max('cwx6 ','In schemes     ',cwx6,size(cwx6))


call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
call prepare (cibx6,cicx6,cifx6 ,cisx6 ,ciwx6 ,nxm)
call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
call prepare (cbi6,cci6,cfi6 ,csi6 ,cwi6 ,nx)
call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
call prepare (cibi6,cici6,cifi6 ,cisi6 ,ciwi6 ,nx)
call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
call prepare (cby6,ccy6,cfy6 ,csy6 ,cwy6 ,nym)
call prepare (cby6,ccy6,cfyp6,csyp6,cwyp6,nym)
call prepare (ciby6,cicy6,cify6 ,cisy6 ,ciwy6 ,nym)
call prepare (ciby6,cicy6,cifyp6,cisyp6,ciwyp6,nym)
call prepare (cbi6y,cci6y,cfi6y ,csi6y ,cwi6y ,ny)
call prepare (cbi6y,cci6y,cfip6y,csip6y,cwip6y,ny)
call prepare (cibi6y,cici6y,cifi6y ,cisi6y ,ciwi6y ,ny)
call prepare (cibi6y,cici6y,cifip6y,cisip6y,ciwip6y,ny)

#ifndef TWOD
   call prepare (fbz,fcz,ffz ,fsz ,fwz ,nz)
   call prepare (fbz,fcz,ffzp,fszp,fwzp,nz)
   call prepare (cbz6,ccz6,cfz6 ,csz6 ,cwz6 ,nzm)
   call prepare (cbz6,ccz6,cfzp6,cszp6,cwzp6,nzm)
   call prepare (cibz6,cicz6,cifz6 ,cisz6 ,ciwz6 ,nzm)
   call prepare (cibz6,cicz6,cifzp6,ciszp6,ciwzp6,nzm)
   call prepare (cbi6z,cci6z,cfi6z ,csi6z ,cwi6z ,nz)
   call prepare (cbi6z,cci6z,cfip6z,csip6z,cwip6z,nz)
   call prepare (cibi6z,cici6z,cifi6z ,cisi6z ,ciwi6z ,nz)
   call prepare (cibi6z,cici6z,cifip6z,cisip6z,ciwip6z,nz)
#endif

if (nclx.eq.1) then
   fbx(nx-1)=0._mytype
   call prepare (fbx,fcx,ffxp,fsxp,fwxp,nx)
   cbx6(nxm-1)=0._mytype
   cibx6(nxm)=0._mytype
   cbi6(nx-1)=0._mytype
   cibi6(nx)=0._mytype
   call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
   call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
   call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
   call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
endif
if (nclx.eq.2) then
   cbx6(nxm-1)=0._mytype
   cibx6(nxm)=0._mytype
   cbi6(nx-1)=0._mytype
   cibi6(nx)=0._mytype
   call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
   call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
   call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
   call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
endif
if (ncly.eq.1) then
   fby(ny-1)=0._mytype
   call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
   cby6(nym-1)=0._mytype
   ciby6(nym)=0._mytype
   cbi6y(ny-1)=0._mytype
   cibi6y(ny)=0._mytype 
   call prepare (cby6,ccy6,cfyp6,csyp6,cwyp6,nym)
   call prepare (ciby6,cicy6,cifyp6,cisyp6,ciwyp6,nym)
   call prepare (cbi6y,cci6y,cfip6y,csip6y,cwip6y,ny)
   call prepare (cibi6y,cici6y,cifip6y,cisip6y,ciwip6y,ny)
endif
if (ncly.eq.2) then
   cby6(nym-1)=0._mytype
   ciby6(nym)=0._mytype
   cbi6y(ny-1)=0._mytype
   cibi6y(ny)=0._mytype
   call prepare (cby6,ccy6,cfyp6,csyp6,cwyp6,nym)
   call prepare (ciby6,cicy6,cifyp6,cisyp6,ciwyp6,nym)
   call prepare (cbi6y,cci6y,cfip6y,csip6y,cwip6y,ny)
   call prepare (cibi6y,cici6y,cifip6y,cisip6y,ciwip6y,ny)
endif

#ifndef TWOD
   if (nclz.eq.1) then
      fbz(nz-1)=0._mytype
      call prepare (fbz,fcz,ffzp,fszp,fwzp,nz)
      cbz6(nzm-1)=0._mytype
      cibz6(nzm)=0._mytype
      cbi6z(nz-1)=0._mytype
      cibi6z(nz)=0._mytype 
      call prepare (cbz6,ccz6,cfzp6,cszp6,cwzp6,nzm)
      call prepare (cibz6,cicz6,cifzp6,ciszp6,ciwzp6,nzm)
      call prepare (cbi6z,cci6z,cfip6z,csip6z,cwip6z,nz)
      call prepare (cibi6z,cici6z,cifip6z,cisip6z,ciwip6z,nz)
   endif
   if (nclz.eq.2) then
      cbz6(nzm-1)=0._mytype
      cibz6(nzm)=0._mytype
      cbi6z(nz-1)=0._mytype
      cibi6z(nz)=0._mytype   
      call prepare (cbz6,ccz6,cfzp6,cszp6,cwzp6,nzm)
      call prepare (cibz6,cicz6,cifzp6,ciszp6,ciwzp6,nzm)
      call prepare (cbi6z,cci6z,cfip6z,csip6z,cwip6z,nz)
      call prepare (cibi6z,cici6z,cifip6z,cisip6z,ciwip6z,nz)
   endif
#endif

call prepare (sbx,scx,sfx,ssx,swx,nx)
call prepare (sby,scy,sfy,ssy,swy,ny)
if (nz.gt.1) call prepare (sbz,scz,sfz,ssz,swz,nz)

call prepare (sbx,scx,sfxp,ssxp,swxp,nx)
call prepare (sby,scy,sfyp,ssyp,swyp,ny)
if (nz.gt.1) call prepare (sbz,scz,sfzp,sszp,swzp,nz)

if (nclx.eq.1) then
   sbx(nx-1)=0._mytype
   call prepare (sbx,scx,sfx ,ssx ,swx ,nx)
endif
if (ncly.eq.1) then
   sby(ny-1)=0._mytype
   call prepare (sby,scy,sfy ,ssy ,swy ,ny)
endif
#ifndef TWOD
   if (nclz.eq.1) then
      sbz(nz-1)=0._mytype
      call prepare (sbz,scz,sfz ,ssz ,swz ,nz)
   endif
#endif

return
end subroutine schemes

!*******************************************************************
!
subroutine prepare (b,c,f,s,w,n)
! 
!*******************************************************************

use decomp_2d, only : mytype

implicit none

integer :: i,n
real(mytype), dimension(n) :: b,c,f,s,w

do i=1,n
   w(i)=c(i)
enddo
do i=2,n
   s(i)=b(i-1)/w(i-1)
   w(i)=w(i)-f(i-1)*s(i)
enddo
do i=1,n
   w(i)=1.0_mytype/w(i)
enddo

return
end subroutine prepare
