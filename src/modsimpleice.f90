!> \file modsimpleice.f90

!>
!!  Ice microphysics.
!>
!! Calculates ice microphysics in a cheap scheme without prognostic nr
!!  simpleice is called from *modmicrophysics*
!! \see  Grabowski, 1998, JAS and Khairoutdinov and Randall, 2006, JAS
!!  \author Steef B\"oing, TU Delft
!!  \par Revision list
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!


module modsimpleice
  use modmicrodata

  implicit none
  save
  contains

!> Initializes and allocates the arrays
  subroutine initsimpleice
    use modglobal, only : ih,i1,jh,j1,k1
    use modfields, only : rhobf
    use modmpi,    only : myid

    implicit none
    real:: k

    allocate (qr(2-ih:i1+ih,2-jh:j1+jh,k1)        & ! qr (total precipitation!) converted from a scalar variable
             ,qrp(2-ih:i1+ih,2-jh:j1+jh,k1)       & ! qr tendency due to microphysics only, for statistics
             ,nr(2-ih:i1+ih,2-jh:j1+jh,k1)        & ! qr (total precipitation!) converted from a scalar variable
             ,nrp(2-ih:i1+ih,2-jh:j1+jh,k1)       & ! qr tendency due to microphysics only, for statistics
             ,thlpmcr(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! thl tendency due to microphysics only, for statistics
             ,qtpmcr(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! qt tendency due to microphysics only, for statistics
             ,sed_qr(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! sedimentation rain droplets mixing ratio
             ,sed_nr(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! sedimentation rain droplets fraction
             ,qr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! time-splitting substep qr
             ,nr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! time-splitting substep nr
             ,ilratio(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! partition ratio cloud water vs cloud ice
             ,rsgratio(2-ih:i1+ih,2-jh:j1+jh,k1)  & ! partition ratio rain vs. snow/graupel
             ,sgratio(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! partition ratio snow vs graupel
             ,lambdar(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! slope parameter for rain
             ,lambdas(2-ih:i1+ih,2-jh:j1+jh,k1)   & ! slope parameter for snow
             ,lambdag(2-ih:i1+ih,2-jh:j1+jh,k1))    ! slope parameter for graupel

    allocate (qrmask(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! mask for rain water
             ,qcmask(2-ih:i1+ih,2-jh:j1+jh,k1))     ! mask for cloud water

    allocate(precep(2-ih:i1+ih,2-jh:j1+jh,k1))      ! precipitation for statistics

    allocate(ccrz(k1),ccsz(k1),ccgz(k1))

  end subroutine initsimpleice

!> Cleaning up after the run
  subroutine exitsimpleice
    implicit none
    deallocate(nr,nrp,qr,qrp,thlpmcr,qtpmcr,sed_qr,qr_spl,nr_spl,ilratio,rsgratio,sgratio,lambdar,lambdas,lambdag)
    deallocate(qrmask,qcmask) 
    deallocate(precep)
    deallocate(ccrz,ccsz,ccgz)
  end subroutine exitsimpleice

!> Calculates the microphysical source term.
  subroutine simpleice
    use modglobal, only : ih,jh,i1,j1,k1,rdt,rk3step,timee,kmax,rlv,cp,tup,tdn
    use modfields, only : sv0,svm,svp,qtp,thlp,qt0,ql0,exnf,rhof,tmp0,rhobf
    use modsimpleicestat, only : simpleicetend
    use modmpi,    only : myid
    implicit none
    integer:: i,j,k 
    real:: qrsmall, qrsum, nrsmall,nrsum

    delt = rdt/ (4. - dble(rk3step))

    ! used to check on negative qr and nr
    qrsum=0.
    qrsmall=0.
    nrsum=0.
    nrsmall=0.
    ! reset microphysics tendencies 
    qrp=0.
    nrp=0.
    thlpmcr=0.
    qtpmcr=0.

    ! Density corrected fall speed parameters, see Tomita 2008
    do k=1,k1
    ccrz(k)=ccr*(1.29/rhobf(k))**0.5
    ccsz(k)=ccs*(1.29/rhobf(k))**0.5
    ccgz(k)=ccg*(1.29/rhobf(k))**0.5
    end do

    do k=1,k1
    do i=2,i1
    do j=2,j1
      ! initialise qr
      qr(i,j,k)= sv0(i,j,k,iqr)
      nr(i,j,k) = sv0(i,j,k,inr)
      ! initialise qc mask
      if (ql0(i,j,k) > qcmin) then
        qcmask(i,j,k) = .true.
      else
        qcmask(i,j,k) = .false.
      end if
      ! initialise qr mask and check if we are not throwing away too much rain
      if (l_rain) then
        qrsum = qrsum+qr(i,j,k)
        if (qr(i,j,k) <= qrmin) then
          qrmask(i,j,k) = .false.
          if(qr(i,j,k)<0.) then
          qrsmall = qrsmall-qr(i,j,k)
          nrsmall = nrsmall-nr(i,j,k)
          qr(i,j,k)=0.
          nr(i,j,k)=0.
          end if
        else
          qrmask(i,j,k)=.true.
        endif
      endif
    enddo
    enddo
    enddo

    if (qrsmall > 0.000001*qrsum) then
      write(*,*)'amount of neg. qr thrown away is too high  ',timee,' sec'
    end if
    if (nrsmall > 0.000001*nrsum) then
      write(*,*)'amount of neg. nr thrown away is too high  ',timee,' sec'
    end if

    if(l_warm) then !partitioning and determination of intercept parameter
      do k=1,k1
      do i=2,i1
      do j=2,j1
        ilratio(i,j,k)=1.   ! cloud water vs cloud ice partitioning
        rsgratio(i,j,k)=1.   ! rain vs snow/graupel partitioning
        sgratio(i,j,k)=0.   ! snow versus graupel partitioning
        lambdar(i,j,k)=(aar*n0rr*gamb1r/(rhof(k)*((qr(i,j,k)/(nr(i,j,k)+1.e-6))+1.e-6)))**(1./(1.+bbr)) ! lambda rain
        lambdas(i,j,k)=lambdar(i,j,k) ! lambda snow
        lambdag(i,j,k)=lambdar(i,j,k) ! lambda graupel
      enddo
      enddo
      enddo
    elseif(l_graupel) then
      do k=1,k1
      do i=2,i1
      do j=2,j1
        ilratio(i,j,k)=amax1(0.,amin1(1.,(tmp0(i,j,k)-tdn)/(tup-tdn)))   ! cloud water vs cloud ice partitioning
        rsgratio(i,j,k)=amax1(0.,amin1(1.,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg)))   ! rain vs snow/graupel partitioning
        sgratio(i,j,k)=amax1(0.,amin1(1.,(tmp0(i,j,k)-tdnsg)/(tupsg-tdnsg)))   ! snow versus graupel partitioning
        lambdar(i,j,k)=(aar*n0rr*gamb1r/(rhof(k)*((qr(i,j,k)/(nr(i,j,k)+1.e-6))*rsgratio(i,j,k)+1.e-6)))**(1./(1.+bbr)) ! lambda rain
        lambdas(i,j,k)=(aas*n0rs*gamb1s/(rhof(k)*((qr(i,j,k)/(nr(i,j,k)+1.e-6))*(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))+1.e-6)))**(1./(1.+bbs)) ! lambda snow
        lambdag(i,j,k)=(aag*n0rg*gamb1g/(rhof(k)*((qr(i,j,k)/(nr(i,j,k)+1.e-6))*(1.-rsgratio(i,j,k))*sgratio(i,j,k)+1.e-6)))**(1./(1.+bbg)) ! lambda graupel
      enddo
      enddo
      enddo
    else
      do k=1,k1
      do i=2,i1
      do j=2,j1
        ilratio(i,j,k)=amax1(0.,amin1(1.,(tmp0(i,j,k)-tdn)/(tup-tdn)))   ! cloud water vs cloud ice partitioning
        rsgratio(i,j,k)=amax1(0.,amin1(1.,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg)))   ! rain vs snow/graupel partitioning
        sgratio(i,j,k)=0.
        lambdar(i,j,k)=(aar*n0rr*gamb1r/(rhof(k)*((qr(i,j,k)/(nr(i,j,k)+1.e-6))*rsgratio(i,j,k)+1.e-6)))**(1./(1.+bbr)) ! lambda rain
        lambdas(i,j,k)=(aas*n0rs*gamb1s/(rhof(k)*((qr(i,j,k)/(nr(i,j,k)+1.e-6))*(1.-rsgratio(i,j,k))+1.e-6)))**(1./(1.+bbs)) ! lambda snow
        lambdag(i,j,k)=lambdas(i,j,k)
      enddo
      enddo
      enddo
    endif

    if (l_rain) then
      call simpleicetend
      call autoconvert
      call simpleicetend
      call accrete
      call simpleicetend
      call evapdep
      call simpleicetend
      call precipitate
    endif

    do k=1,k1
    do i=2,i1
    do j=2,j1
      svp(i,j,k,inr)=svp(i,j,k,inr)+nrp(i,j,k)
      svp(i,j,k,iqr)=svp(i,j,k,iqr)+qrp(i,j,k)
      thlp(i,j,k)=thlp(i,j,k)+thlpmcr(i,j,k)
      qtp(i,j,k)=qtp(i,j,k)+qtpmcr(i,j,k)
      ! adjust negative qr tendencies at the end of the time-step
      if (svp(i,j,k,iqr)+svm(i,j,k,iqr)/delt < qrmin) then
        svp(i,j,k,iqr) = - svm(i,j,k,iqr)/delt
        qtp(i,j,k) = qtp(i,j,k) + svm(i,j,k,iqr)/delt
        thlp(i,j,k) = thlp(i,j,k) - (rlv/(cp*exnf(k)))*svm(i,j,k,iqr)/delt
      endif
      ! make sure nr constraints are still satisfied
      if (svp(i,j,k,inr)+svm(i,j,k,inr)/delt < 0.) then
        svp(i,j,k,inr) = - svm(i,j,k,inr)/delt
      endif
      if (svp(i,j,k,inr)+svm(i,j,k,inr)/delt > 1.) then
        svp(i,j,k,inr) = (1.-svm(i,j,k,inr))/delt
      endif
    enddo
    enddo
    enddo

    if (l_rain) then
      call simpleicetend !after corrections
    endif
  end subroutine simpleice

  subroutine autoconvert
    use modglobal, only : ih,jh,i1,j1,k1,rlv,cp,tmelt
    use modfields, only : ql0,exnf,rhof,tmp0
    use modmpi,    only : myid
    implicit none
    real :: qll,qli,ddisp,lwc,autl,tc,times,auti,aut
    integer:: i,j,k

    if(l_berry.eqv..true.) then ! Berry/Hsie autoconversion
      do k=1,k1
      do i=2,i1
      do j=2,j1
        if (qcmask(i,j,k).eqv..true.) then
          ! ql partitioning 
          qll=ql0(i,j,k)*ilratio(i,j,k)
          qli=ql0(i,j,k)-qll
          ddisp=0.146-5.964e-2*alog(Nc_0/2.e9) ! Relative dispersion coefficient for Berry autoconversion
          lwc=1.e3*rhof(k)*qll ! Liquid water content in g/kg
          autl=1./rhof(k)*1.67e-5*lwc*lwc/(5. + .0366*Nc_0/(1.e6*ddisp*(lwc+1.e-6)))
          tc=tmp0(i,j,k)-tmelt ! Temperature wrt melting point
          times=amin1(1.e3,(3.56*tc+106.7)*tc+1.e3) ! Time scale for ice autoconversion
          auti=qli/times
          aut = min(autl + auti,ql0(i,j,k)/delt)
          qrp(i,j,k) = qrp(i,j,k)+aut
          qtpmcr(i,j,k) = qtpmcr(i,j,k)-aut
          thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*aut
        endif
      enddo
      enddo
      enddo
    else ! Lin/Kessler autoconversion as in Khairoutdinov and Randall, 2006
      do k=1,k1
      do i=2,i1
      do j=2,j1
        if (qcmask(i,j,k).eqv..true.) then
          ! ql partitioning 
          qll=ql0(i,j,k)*ilratio(i,j,k)
          qli=ql0(i,j,k)-qll
          autl=max(0.,timekessl*(qll-qll0))
          tc=tmp0(i,j,k)-tmelt
          auti=max(0.,betakessi*exp(0.025*tc)*(qli-qli0))
          aut = min(autl + auti,ql0(i,j,k)/delt)
          qrp(i,j,k) = qrp(i,j,k)+aut
          qtpmcr(i,j,k) = qtpmcr(i,j,k)-aut
          thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*aut
        endif
      enddo
      enddo
      enddo
    endif

  end subroutine autoconvert

  subroutine accrete
    use modglobal, only : ih,jh,i1,j1,k1,rlv,cp,pi
    use modfields, only : ql0,exnf,rhof
    use modmpi,    only : myid
    implicit none
    real :: qll,qli,qrr,qrs,qrg,conr,cons,cong,massr,masss,massg,diamr,diams,diamg,gaccrl,gaccsl,gaccgl,gaccri,gaccsi,gaccgi,accr,accs,accg,acc
    integer:: i,j,k

    do k=1,k1
    do i=2,i1
    do j=2,j1
      if (qrmask(i,j,k).eqv..true. .and. qcmask(i,j,k).eqv..true.) then ! apply mask
        ! ql partitioning
        qll=ql0(i,j,k)*ilratio(i,j,k)
        qli=ql0(i,j,k)-qll
        ! qr partitioning 
        qrr=qr(i,j,k)*rsgratio(i,j,k)
        qrs=qr(i,j,k)*(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))
        qrg=qr(i,j,k)*(1.-rsgratio(i,j,k))*sgratio(i,j,k)
        conr=n0rr/lambdar(i,j,k) ! rain concentration
        cons=n0rs/lambdas(i,j,k) ! snow concentration
        cong=n0rg/lambdag(i,j,k) ! graupel concentration
        massr =rhof(k)*(qrr+1.e-7)/conr  ! mass of rain
        masss =rhof(k)*(qrs+1.e-7)/cons  ! mass of snow
        massg =rhof(k)*(qrg+1.e-7)/cong  ! mass of graupel
        diamr=(massr/aar)**(1./bbr) ! diameter liquid particles
        diams=(masss/aas)**(1./bbs) ! diameter ice particles
        diamg=(massg/aag)**(1./bbg) ! diameter ice particles
        gaccrl=pi/4.*ccrz(k)*diamr**(2.+ddr)*ceffrl*rhof(k)*qll ! collection of cloud water by rain etc.
        gaccsl=pi/4.*ccsz(k)*diams**(2.+dds)*ceffsl*rhof(k)*qll
        gaccgl=pi/4.*ccgz(k)*diamg**(2.+ddg)*ceffgl*rhof(k)*qll
        gaccri=pi/4.*ccrz(k)*diamr**(2.+ddr)*ceffri*rhof(k)*qli
        gaccsi=pi/4.*ccsz(k)*diams**(2.+dds)*ceffsi*rhof(k)*qli
        gaccgi=pi/4.*ccgz(k)*diamg**(2.+ddg)*ceffgi*rhof(k)*qli
        accr=conr*(gaccrl+gaccri)*qrr/(qrr+1.e-9)
        accs=cons*(gaccsl+gaccsi)*qrs/(qrs+1.e-9)
        accg=cong*(gaccgl+gaccgi)*qrg/(qrg+1.e-9)
        acc= min(accr+accs+accg,ql0(i,j,k)/delt)  ! total growth by accretion
        qrp(i,j,k) = qrp(i,j,k)+acc
        qtpmcr(i,j,k) = qtpmcr(i,j,k)-acc
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*acc
      end if
    enddo
    enddo
    enddo

  end subroutine accrete

  subroutine evapdep
    use modglobal, only : ih,jh,i1,j1,k1,rlv,riv,cp,rv,rd,tmelt,es0,pi
    use modfields, only : qt0,ql0,exnf,rhof,tmp0,presf,qvsl,qvsi,esl
    use modmpi,    only : myid
    implicit none
    real :: qrr,qrs,qrg,ssl,ssi,conr,cons,cong,massr,masss,massg,diamr,diams,diamg,rer,res,reg,ventr,vents,ventg,thfun,gevapdepr,gevapdeps,gevapdepg,evapdepr,evapdeps,evapdepg,devap
    integer:: i,j,k

    do k=1,k1
    do i=2,i1
    do j=2,j1
      if (qrmask(i,j,k).eqv..true.) then
        ! saturation ratios
        ssl=(qt0(i,j,k)-ql0(i,j,k))/qvsl(i,j,k)
        ssi=(qt0(i,j,k)-ql0(i,j,k))/qvsi(i,j,k)
        !integration over ventilation factors and diameters, see e.g. seifert
        ventr=.78*n0rr/lambdar(i,j,k)**2 + gam2dr*.27*n0rr*sqrt(ccrz(k)/2.e-5)*lambdar(i,j,k)**(-2.5-0.5*ddr)
        vents=.78*n0rs/lambdas(i,j,k)**2 + gam2ds*.27*n0rs*sqrt(ccsz(k)/2.e-5)*lambdas(i,j,k)**(-2.5-0.5*dds)
        ventg=.78*n0rg/lambdag(i,j,k)**2 + gam2dg*.27*n0rg*sqrt(ccgz(k)/2.e-5)*lambdag(i,j,k)**(-2.5-0.5*ddg)
        thfun=1.e-7/(2.2*tmp0(i,j,k)/esl(i,j,k)+2.2e2/tmp0(i,j,k))  ! thermodynamic function
        evapdepr=(4.*pi/(betar*rhof(k)))*(ssl-1.)*ventr*thfun
        evapdeps=(4.*pi/(betas*rhof(k)))*(ssi-1.)*vents*thfun
        evapdepg=(4.*pi/(betag*rhof(k)))*(ssi-1.)*ventg*thfun
        ! limit with qr and ql after accretion and autoconversion
        devap= max(min(evapfactor*(evapdepr+evapdeps+evapdepg),ql0(i,j,k)/delt+qrp(i,j,k)),-qr(i,j,k)/delt-qrp(i,j,k))  ! total growth by deposition and evaporation
        qrp(i,j,k) = qrp(i,j,k)+devap
        qtpmcr(i,j,k) = qtpmcr(i,j,k)-devap
        thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*devap
      end if
    enddo
    enddo
    enddo

  end subroutine evapdep

  subroutine precipitate
    use modglobal, only : ih,i1,jh,j1,k1,kmax,dzf,pi,dzh
    use modfields, only : rhof,rhobf
    use modmpi,    only : myid
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real :: dt_spl,wfallmax,vtr,vts,vtg,vtf

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)*courantp))
    dt_spl = delt/real(n_spl) !fixed time step

    sed_qr = 0. ! reset sedimentation fluxes
    sed_nr = 0. ! reset sedimentation fluxes

    do k=1,k1
    do i=2,i1
    do j=2,j1
        qr_spl(i,j,k) = qr(i,j,k)
        nr_spl(i,j,k) = nr(i,j,k)
      if (qrmask(i,j,k).eqv..true.) then
        vtr=ccrz(k)*(gambd1r/gamb1r)/(lambdar(i,j,k)**ddr)  ! terminal velocity rain
        vts=ccsz(k)*(gambd1s/gamb1s)/(lambdas(i,j,k)**dds)  ! terminal velocity snow
        vtg=ccgz(k)*(gambd1g/gamb1g)/(lambdag(i,j,k)**ddg)  ! terminal velocity graupel
        vtf=rsgratio(i,j,k)*vtr+(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))*vts+(1.-rsgratio(i,j,k))*sgratio(i,j,k)*vtg  ! mass-weighted terminal velocity
        vtf = amin1(wfallmax,vtf)
        precep(i,j,k) = vtf*qr_spl(i,j,k)*nr_spl(i,j,k)  ! note multiplication by nr_spl
        sed_qr(i,j,k) = precep(i,j,k)*rhobf(k) ! convert to flux
        sed_nr(i,j,k) = vtf*nr_spl(i,j,k)*rhobf(k)  ! note multiplication by nr_spl
      else
        precep(i,j,k) = 0.
        sed_qr(i,j,k) = 0.
        sed_nr(i,j,k) = 0.
      end if
    enddo
    enddo
    enddo
   
    !  advect precipitation using upwind scheme
    do k=1,kmax
    do i=2,i1
    do j=2,j1
      qr_spl(i,j,k) = qr_spl(i,j,k) + (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzh(k+1)*rhobf(k))
      nr_spl(i,j,k) = max(0.,min(nr_spl(i,j,k) + (sed_nr(i,j,k+1) - sed_nr(i,j,k))*dt_spl/(dzh(k+1)*rhobf(k)),1.))
    enddo
    enddo
    enddo

    ! begin time splitting loop
    IF (n_spl > 1) THEN
      DO jn = 2 , n_spl 
    
        ! reset fluxes at each step of loop
        sed_qr = 0.
        sed_nr = 0.    
        do k=1,k1
        do i=2,i1
        do j=2,j1
          if (qr_spl(i,j,k) > qrmin) then
            ! re-evaluate lambda
            lambdar(i,j,k)=(aar*n0rr*gamb1r/(rhof(k)*(qr_spl(i,j,k)/(nr_spl(i,j,k)+1.e-6)*rsgratio(i,j,k)+1.e-6)))**(1./(1.+bbr)) ! lambda rain
            lambdas(i,j,k)=(aas*n0rs*gamb1s/(rhof(k)*(qr_spl(i,j,k)/(nr_spl(i,j,k)+1.e-6)*(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))+1.e-6)))**(1./(1.+bbs)) ! lambda snow
            lambdag(i,j,k)=(aag*n0rg*gamb1g/(rhof(k)*(qr_spl(i,j,k)/(nr_spl(i,j,k)+1.e-6)*(1.-rsgratio(i,j,k))*sgratio(i,j,k)+1.e-6)))**(1./(1.+bbg)) ! lambda graupel
            vtr=ccrz(k)*(gambd1r/gamb1r)/(lambdar(i,j,k)**ddr)  ! terminal velocity rain
            vts=ccsz(k)*(gambd1s/gamb1s)/(lambdas(i,j,k)**dds)  ! terminal velocity snow
            vtg=ccgz(k)*(gambd1g/gamb1g)/(lambdag(i,j,k)**ddg)  ! terminal velocity graupel
            vtf=rsgratio(i,j,k)*vtr+(1.-rsgratio(i,j,k))*(1.-sgratio(i,j,k))*vts+(1.-rsgratio(i,j,k))*sgratio(i,j,k)*vtg  ! mass-weighted terminal velocity
            vtf=amin1(wfallmax,vtf)
            sed_qr(i,j,k) = vtf*qr_spl(i,j,k)*rhobf(k)*nr_spl(i,j,k) ! note multiplication by nr_spl
            sed_nr(i,j,k) = vtf*nr_spl(i,j,k)*rhobf(k)*nr_spl(i,j,k) ! note multiplication by nr_spl
          else
            sed_qr(i,j,k) = 0.
            sed_nr(i,j,k) = 0.
          endif
        enddo
        enddo
        enddo

        do k=1,kmax
        do i=2,i1
        do j=2,j1
          qr_spl(i,j,k) = qr_spl(i,j,k) + (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzf(k)*rhobf(k))
          nr_spl(i,j,k) = max(0.,min(nr_spl(i,j,k) + (sed_nr(i,j,k+1) - sed_nr(i,j,k))*dt_spl/(dzf(k)*rhobf(k)),1.))
        enddo
        enddo
        enddo
    
      ! end time splitting loop and if n>1
      ENDDO
    ENDIF

    ! no thl and qt tendencies build in, implying no heat transfer between precipitation and air
    do k=1,kmax
    do i=2,i1
    do j=2,j1
      qrp(i,j,k)= qrp(i,j,k) + (qr_spl(i,j,k) - qr(i,j,k))/delt
      nrp(i,j,k)= nrp(i,j,k) + (nr_spl(i,j,k) - nr(i,j,k))/delt
    enddo
    enddo
    enddo

  end subroutine precipitate

end module modsimpleice
