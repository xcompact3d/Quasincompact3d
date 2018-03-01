program verif

  implicit none

  integer, parameter :: nx1=129 , ny1=128, nz1=84
  real(8),dimension(nx1,ny1,nz1) :: umean,uumean
  integer :: i,j,k,count
  real(8) :: x,y,xitime,u_to
 real(8),dimension(ny1) :: yp,ypi,qstat


u_to=180./4250.

 xitime=3500.

   open (15,file='yp.dat',form='formatted',status='unknown')
   do j=1,ny1
      read(15,*) yp(j),ypi(j)
   enddo
   close(15)




 OPEN(11,FILE='vmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) umean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)
OPEN(11,FILE='vvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uumean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)

umean=umean/xitime
uumean=uumean/xitime

uumean=sqrt(uumean-umean*umean)/u_to

         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            qstat(j)=qstat(j)+uumean(i,j,k)
         enddo
         enddo
         enddo

         qstat(:)=qstat(:)/nx1/nz1


  open(10,file='testpara.dat',status='unknown',form='formatted')
  do j=1,ny1
     write(10,*) yp(j)*180.+180.,qstat(j)
  enddo
  close(10)



    end program verif

