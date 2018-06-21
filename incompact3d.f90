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
  double precision :: t1,t2,t1poiss,t2poiss,tpoisstotal
  character(len=20) :: filename

  logical :: converged
  integer :: poissiter, totalpoissiter

  TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

  CALL MPI_INIT(code)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  !start from 1 == true
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)
  call parameter()

  call init_variables

  call schemes()

  if (nclx.eq.0) then
    bcx=0
  else
    bcx=1
  endif
  if (ncly.eq.0) then
    bcy=0
  else
    bcy=1
  endif
  if (nclz.eq.0) then
    bcz=0
  else
    bcz=1
  endif

  call decomp_2d_poisson_init(bcx,bcy,bcz)

  call decomp_info_init(nxm,nym,nzm,phG)

  !if you want to collect 100 snapshots randomly on 50000 time steps
  !call collect_data() !it will generate 100 random time steps

  totalpoissiter = 0
  if (ilit.eq.0) then
    t = 0._mytype
    itime = 0
    call init(ux1,uy1,uz1,rho1,temperature1,massfrac1,ep1,phi1,&
         gx1,gy1,gz1,rhos1,temperatures1,massfracs1,phis1,&
         hx1,hy1,hz1,rhoss1,temperaturess1,massfracss1,phiss1,&
         pressure0)
    pp3star(:,:,:) = 0._mytype
  else
    call restart(ux1,uy1,uz1,rho1,temperature1,ep1,pp3,phi1,&
         gx1,gy1,gz1,rhos1,px1,py1,pz1,phis1,&
         hx1,hy1,hz1,rhoss1,phiss1,&
         pressure0,phG,0)
  endif

  ! XXX LMN: Calculate divergence of velocity field. Also updates rho in Y
  !          and Z pencils.
  !          X->Y->Z
  call calc_divu(ta1,tb1,tc1,rho1,temperature1,massfrac1,kappa1,di1,&
       ta2,tb2,tc2,rho2,temperature2,massfrac2,kappa2,di2,&
       divu3,ta3,tb3,rho3,temperature3,massfrac3,kappa3,di3,&
       pressure0)

  call test_speed_min_max(ux1,uy1,uz1)
  if (isolvetemp.eq.0) then
    call test_density_min_max(rho1)
  else
    call test_temperature_min_max(temperature1)
  endif
  if (iscalar.eq.1) then
    call test_scalar_min_max(phi1)
  endif

  !array for stat to zero
  umean=0._mytype;vmean=0._mytype;wmean=0._mytype
  uumean=0._mytype;vvmean=0._mytype;wwmean=0._mytype
  uvmean=0._mytype;uwmean=0._mytype;vwmean=0._mytype
  phimean=0._mytype;phiphimean=0._mytype

  t1 = MPI_WTIME()

  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)

  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)  
  call decomp_info_init(nxm, nym, nz, ph3) 

  itime=0
  call VISU_INSTA(ux1,uy1,uz1,rho1,temperature1,massfrac1,phi1,&
       ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
       ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
       ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
  ! call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
  !      ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)

  ! call VISU_INSTB(ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
  !      ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
  !      ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)

  ! CALL track_front(ux1, rho1)
  ! CALL track_front_height(rho1, rho2, rho3)
  ! CALL calc_energy_budgets(rho1, ux1, uy1, uz1, ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, &
  !           di1, &
  !           rho2, ux2, uy2, uz2, ta2, tb2, tc2, td2, te2, tf2, di2, &
  !           rho3, ux3, uy3, uz3, ta3, tb3, tc3, di3)

  tpoisstotal = 0._mytype
  do itime=ifirst,ilast

    t=(itime-1)*dt
    if (nrank.eq.0) then
      write(*,1001) itime,t
1001  format('Time step =',i7,', Time unit =',F9.3)
    endif

    do itr=1,iadvance_time
      
      !-----------------------------------------------------------------------------------
      ! XXX ux,uy,uz now contain velocity: ux = u etc.
      !-----------------------------------------------------------------------------------

      if (nclx.eq.2) then
        call inflow (ux1,uy1,uz1,rho1,temperature1,massfrac1,phi1) !X PENCILS
        if ((ilmn.ne.0).and.(itime.eq.ifirst).and.(itr.eq.1)) then
          call compute_outflux_lmn(temperature1,ta1,di1,&
               temperature2,ta2,di2,&
               temperature3,ta3,di3)
        endif
        call outflow(ux1,uy1,uz1,rho1,temperature1,massfrac1,phi1) !X PENCILS 
      endif

      if (ilmn.ne.0) then
         call set_density_bcs(rho1, ux1, uy1, uz1)
         !! if (itype.eq.5) then
         !!    if (ncly.eq.2) then
         !!       call set_density_entrainment_y(rho1, uy1)
         !!    endif
         !!    if (nclz.eq.2) then
         !!       call set_density_entrainment_z(rho1, uz1)
         !!    endif
         !! endif
      endif

      !! Ensure rho/temp is up to date
      if (isolvetemp.eq.0) then
        call calctemp_eos(temperature1, rho1, massfrac1, pressure0, xsize)
      else
        call calcrho_eos(rho1, temperature1, massfrac1, pressure0, xsize)
      endif

      !-----------------------------------------------------------------------------------
      ! Update fluid properties
      ! XXX Temperature is up-to-date in X, Y and Z.
      !-----------------------------------------------------------------------------------
      call calcvisc(mu1, mu2, mu3, rho1, temperature1, massfrac1)
      if (iscalar.eq.0) then
        call calcgamma(gamma1, temperature1)
      endif

      !X-->Y-->Z-->Y-->X
      call convdiff(ux1,uy1,uz1,rho1,mu1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ux2,uy2,uz2,rho2,mu2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ux3,uy3,uz3,rho3,mu3,divu3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)
      call apply_grav(ta1, tb1, tc1, rho1)

      ! Transport massfrac
      if (imulticomponent.ne.0) then
        call convdiff_massfrac(ux1,uy1,uz1,rho1,massfrac1,massfracs1,massfracss1,tg1,th1,ti1,di1,&
             uy2,uz2,rho2,massfrac2,ta2,tb2,tc2,di2,&
             uz3,rho3,massfrac3,ta3,tb3,tc3,di3)
      endif

      if (iscalar.eq.1) then
        !---------------------------------------------------------------------------------
        ! XXX After this phi1 contains rho*phi
        !---------------------------------------------------------------------------------
        call scalar(ux1,uy1,uz1,rho1,phi1,gamma1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
             uy2,uz2,rho2,phi2,gamma2,di2,ta2,tb2,tc2,td2,&
             uz3,rho3,phi3,gamma3,di3,ta3,tb3,tc3,&
             ep1)
      endif

      if (ilmn.ne.0) then
         if (isolvetemp.eq.0) then
            ! Update density
            !    X->Y->Z->Y->X
            ! XXX uz3,rho3 and uy2,rho2 and rho1 should already be up to date, could go from 8 to 2
            !     transpose operations by operating on Z->Y->X.
            ! XXX tg1 contains the density forcing term.
            call conv_density(ux1,uy1,uz1,rho1,di1,tg1,th1,ti1,td1,&
                 uy2,uz2,rho2,di2,ta2,tb2,tc2,td2,&
                 uz3,rho3,divu3,di3,ta3,tb3,ep1)
         else
            ! Update temperature
            call convdiff_temperature(ux1,uy1,uz1,rho1,temperature1,di1,tg1,th1,&
                 uy2,uz2,rho2,temperature2,di2,ta2,tb2,&
                 uz3,rho3,temperature3,di3,ta3,tb3)
            call eval_densitycoeffs(rho1,temperature1,tg1,rhos1,rhoss1,rhos01,drhodt1)
            call intttemperature(temperature1,temperatures1,temperaturess1,tg1)
         endif
      endif

      !X PENCILS
      call intt (ux1,uy1,uz1,gx1,gy1,gz1,hx1,hy1,hz1,ta1,tb1,tc1,rho1)

      !-----------------------------------------------------------------------------------
      ! XXX ux,uy,uz now contain momentum: ux = (rho u) etc.
      !-----------------------------------------------------------------------------------

      if (ilmn.ne.0) then
        !! Update density
        if (isolvetemp.eq.0) then
          call inttdensity(rho1,rhos1,rhoss1,rhos01,tg1,drhodt1)

          ! Update temperature using EOS
          call calctemp_eos(temperature1, rho1, massfrac1, pressure0, xsize)
          call test_temperature_min_max(temperature1)
        else
          ! Update density using EOS
          call calcrho_eos(rho1, temperature1, massfrac1, pressure0, xsize)
          call test_density_min_max(rho1)
        endif

        if (ivarcoeff.eq.0) then
           !! Predict drhodt at new timestep
           call extrapol_rhotrans(rho1,rhos1,rhoss1,rhos01,drhodt1)
           
           ! !! Apply Birman correction
           ! call birman_rhotrans_corr(rho1, drhodt1, ta1, tb1, di1, rho2, &
           !      ta2, tb2, di2, &
           !      rho3, ta3, di3)
        endif
      endif

      ! Predict new pressure field
      if ((ilmn.ne.0).and.(ivarcoeff.ne.0)) then
        if ((itime - 1).ne.0) then
          pp3star(:,:,:) = 2._mytype * pp3(:,:,:) - pp3star(:,:,:)
          pp3corr(:,:,:) = pp3(:,:,:) ! Store current pressure field
        else
          !! Temporarily store current pressure field
          pp3corr(:,:,:) = 0._mytype
        endif
      endif

      if (max(nclx, ncly, nclz).eq.2) then
        !! We have Dirichlet boundaries, compute the boundary velocity flux constraint
        call compute_outflux_lmn(temperature1,ta1,di1,&
             temperature2,ta2,di2,&
             temperature3,ta3,di3)
      endif
      call pre_correc(ux1,uy1,uz1,rho1)

      ! LMN: Calculate new divergence of velocity using new density/temperature field.
      !      This updates the temperature field using the density field.
      !      After this rho1,rho2,rho3,temperature1,temperature2,temperature3 are all
      !      upto date.
      !
      !    X->Y->Z
      call calc_divu(tg1,th1,ti1,rho1,temperature1,massfrac1,kappa1,di1,&
           ta2,tb2,tc2,rho2,temperature2,massfrac2,kappa2,di2,&
           divu3,ta3,tb3,rho3,temperature3,massfrac3,kappa3,di3,&
           pressure0)

!!$      if (ivirt.eq.1) then !solid body old school
!!$         !we are in X-pencil
!!$         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
!!$         call body(ux1,uy1,uz1,ep1)
!!$         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
!!$      endif

      if (ilmn.ne.0) then
        !! Conserved->primitive variables
        if (ivarcoeff.ne.0) then
          ux1(:,:,:) = ux1(:,:,:) / rho1(:,:,:)
          uy1(:,:,:) = uy1(:,:,:) / rho1(:,:,:)
          uz1(:,:,:) = uz1(:,:,:) / rho1(:,:,:)
        endif

        if (iscalar.ne.0) then
          phi1(:,:,:) = phi1(:,:,:) / rho1(:,:,:)
        endif
      endif
      
      !X-->Y-->Z
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,drhodt1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,&
           ta3,tb3,tc3,di3,td3,te3,tf3,divu3,pp3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,1,.FALSE.)

      !-----------------------------------------------------------------------------------
      ! Solution of the Poisson equation
      converged = .FALSE.
      poissiter = 0
      t1poiss = MPI_WTIME()
      do while(converged.eqv..FALSE.)
         if (ilmn.ne.0) then
            if (ivarcoeff.ne.0) then
               !! LMN: variable coefficient Poisson
               
               if ((nrank.eq.0).and.(poissiter.eq.0)) then
                  print *, "Solving variable-coefficient pressure-Poisson equation"
               endif
               ! if (poissiter.ne.0) then
               !   !! Compute correction term
               !   call divergence_corr(rho1, px1, py1, pz1, ta1, tb1, tc1, td1, te1, tf1, di1, &
               !        te2, tf2, ta2, tb2, tc2, td2, di2, &
               !        td3, pp3corr, ta3, tb3, tc3, di3, rho0p3, pp3, tg3, &
               !        nxmsize, nymsize, nzmsize, ph1, ph2, ph3, ph4, &
               !        divup3norm, poissiter, converged)
               ! else
               !   !! Need an initial guess for 1/rho0 nabla^2 p - div( 1/rho nabla p )
               !   call approx_divergence_corr(ux1, uy1, uz1, rho1, ta1, tb1, tc1, td1, te1, tf1, ep1, &
               !        di1, rhos1, rhoss1, rhos01, drhodt1, &
               !        td2, te2, tf2, di2, ta2, tb2, tc2, &
               !        ta3, tb3, tc3, di3, td3, te3, tf3, tg3, pp3corr, divu3, &
               !        nxmsize, nymsize, nzmsize, ph1, ph3, ph4, &
               !        divup3norm)
               !   pp3corr(:,:,:) = pp3(:,:,:)
               ! endif
               
               if (poissiter.eq.0) then
                  call calc_divup3norm(divup3norm, divu3, ta1, tb1, di1, ta2, tb2, tc2, di2, &
                       ta3, tb3, di3, nxmsize, nymsize, nzmsize, ph1, ph3, ph4)
                  if((itime - 1).eq.0) then
                     px1(:,:,:) = 0._mytype
                     py1(:,:,:) = 0._mytype
                     pz1(:,:,:) = 0._mytype
                  else
                     call gradp(px1,py1,pz1,di1,td2,tf2,ta2,tb2,tc2,di2,&
                          ta3,tc3,di3,pp3star,nxmsize,nymsize,nzmsize,ph2,ph3)
                  endif
                  pp3star(:,:,:) = pp3corr(:,:,:) ! Store current pressure field
               endif
               call divergence_corr(rho1, px1, py1, pz1, ta1, tb1, tc1, td1, te1, tf1, di1, &
                    te2, tf2, ta2, tb2, tc2, td2, di2, &
                    td3, pp3corr, ta3, tb3, tc3, di3, rho0p3, pp3, tg3, &
                    nxmsize, nymsize, nzmsize, ph1, ph2, ph3, ph4, &
                    divup3norm, poissiter, converged)
            else
               !! LMN: constant coefficient Poisson
               pp3corr(:,:,:) = pp3(:,:,:)
            endif
         else
            !! Incompressible
            pp3corr(:,:,:) = pp3(:,:,:)
         endif
         
         if (converged.eqv..FALSE.) then
            !POISSON Z-->Z
            call decomp_2d_poisson_stg(pp3corr,bcx,bcy,bcz)
            
            !Z-->Y-->X
            ! XXX Need to call this now as if using var-coeff
            !     Poisson equation we will need new values of
            !     gradp on next iteration.
            call gradp(px1,py1,pz1,di1,td2,tf2,ta2,tb2,tc2,di2,&
                 ta3,tc3,di3,pp3corr,nxmsize,nymsize,nzmsize,ph2,ph3)
         endif
         
         if ((ilmn.eq.0).or.(ivarcoeff.eq.0)) then
            converged = .TRUE.
         endif
         
         poissiter = poissiter + 1
      enddo ! End Poisson loop
      t2poiss = MPI_WTIME()
      tpoisstotal = tpoisstotal + (t2poiss - t1poiss)
      pp3(:,:,:) = pp3corr(:,:,:) ! Set pressure field

      if (nrank.eq.0) then
        print *, "Solved Poisson equation in ", poissiter, " iteration(s), took ", t2poiss - t1poiss, "s"
        totalpoissiter = totalpoissiter + poissiter
      endif
      !-----------------------------------------------------------------------------------

      !X PENCILS
      call corgp(ux1,ux2,uy1,uz1,px1,py1,pz1,rho1) 

      !-----------------------------------------------------------------------------------
      ! XXX ux,uy,uz now contain velocity: ux = u etc.
      !-----------------------------------------------------------------------------------

      !does not matter -->output=DIV U=0 (in dv3)
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,drhodt1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,&
           ta3,tb3,tc3,di3,td3,te3,tf3,divu3,dv3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,2,.FALSE.)

      call test_speed_min_max(ux1,uy1,uz1)
      if (iscalar.eq.1) then
        call test_scalar_min_max(phi1)
      endif

      ! Move time to end of substep
      t = t + gdt(itr)
    enddo ! End sub-timesteps

!!$   call STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
!!$        uvmean,uwmean,vwmean,phiphimean,tmean)

    if (mod(itime,isave).eq.0) then
      call restart(ux1,uy1,uz1,rho1,temperature1,ep1,pp3,phi1,&
           gx1,gy1,gz1,rhos1,px1,py1,pz1,phis1,&
           hx1,hy1,hz1,rhoss1,phiss1,&
           pressure0,phG,1)
    endif

    if (mod(itime,imodulo).eq.0) then
      call VISU_INSTA(ux1,uy1,uz1,rho1,temperature1,massfrac1,phi1,&
           ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
      call VISU_PRE (pp3,ta1,tb1,di1,&
           ta2,tb2,di2,&
           ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
    endif

    ! if (mod(itime,10).eq.0) then
    !   call VISU_INSTB(ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
    !        ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
    !        ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
    ! endif

    ! ! MMS: compare errors
    ! CALL eval_error_rho(rho1)
    ! CALL eval_error_vel(ux1,uy1,uz1)

    ! IF (MOD(itime, 10).EQ.0) THEN
    !    CALL track_front(ux1, rho1)
    !    CALL track_front_height(rho1, rho2, rho3)
    !    CALL calc_energy_budgets(rho1, ux1, uy1, uz1, ta1, tb1, tc1, td1, te1, tf1, tg1, th1, ti1, &
    !         di1, &
    !         rho2, ux2, uy2, uz2, ta2, tb2, tc2, td2, te2, tf2, di2, &
    !         rho3, ux3, uy3, uz3, ta3, tb3, tc3, di3)
    ! ENDIF

  enddo

  t2=MPI_WTIME()-t1
  call MPI_ALLREDUCE(t2,t1,1,MPI_REAL8,MPI_SUM, &
       MPI_COMM_WORLD,code)
  t2 = tpoisstotal
  call MPI_ALLREDUCE(t2, tpoisstotal, 1, MPI_REAL8, MPI_SUM, &
       MPI_COMM_WORLD, code)
  if (nrank.eq.0) then
    print *,'time per time_step: ', &
         t1/float(nproc)/(ilast-ifirst+1),' seconds'
    print *,'simulation with nx*ny*nz=',nx,ny,nz,'mesh nodes'
    print *,'Mapping p_row*p_col=',p_row,p_col
    print *,'Mean iterations per Poisson solve: ', &
         float(totalpoissiter) / float(iadvance_time) / float(ilast - ifirst + 1)
    print *,'Time spent in Poisson equation (%): ', &
         100._mytype * tpoisstotal / t1
  endif

  !call decomp_2d_poisson_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(code)

end PROGRAM incompact3d
