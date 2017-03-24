!> \file modtimedepsv.f90
!!  Prescribes surface values, fluxes and LS forcings at certain times for scalars

!>
!!  Prescribes surface values, fluxes and LS forcings at certain times for scalars
!>
!!  \author Roel Neggers, KNMI
!!  \author Thijs Heus,MPI-M
!!  \author Stephan de Roode, TU Delft
!!  \author Simon Axelsen, UU
!!  \par Revision list
!! \todo documentation
!  This file is part of DALES.
!
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

module modtimedep


implicit none
private
public :: inittimedep, timedep,ltimedep,exittimedep

save
! switches for timedependent surface fluxes and large scale forcings
  logical       :: ltimedep     = .true. !< Overall switch, input in namoptions
  logical       :: ltimedepz    = .true.  !< Switch for large scale forcings
  logical       :: ltimedepsurf = .true.  !< Switch for surface fluxes

  integer               :: kflux
  integer               :: kls
  real, allocatable     :: timeflux (:)
  real, allocatable     :: wqsurft  (:)
  real, allocatable     :: wtsurft  (:)
  real, allocatable     :: thlst    (:)
  real, allocatable     :: qtst     (:)
  real, allocatable     :: pst      (:)

  real, allocatable     :: timels  (:)
  real, allocatable     :: ugt     (:,:)
  real, allocatable     :: vgt     (:,:)
  real, allocatable     :: wflst   (:,:)
  real, allocatable     :: dqtdtlst(:,:)
  real, allocatable     :: dthldtlst(:,:)
  real, allocatable     :: thlpcart(:,:)
  real, allocatable     :: dudtlst (:,:)
  real, allocatable     :: dvdtlst (:,:)
  real, allocatable     :: qtproft (:,:)



contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inittimedep

    use modnudge, only :ntnudge    
    use modmpi,    only :myid,my_real,mpi_logical,mpierr,comm3d
    use modglobal, only :btime,cexpnr,k1,kmax,ifinput,runtime,tres, zf
    use modsurfdata,only :ps,qts,wqsurf,wtsurf,thls
    use modtimedepsv, only : inittimedepsv

    use modtestbed,  only : ltestbed,&
                            tb_time,tb_ps,tb_qts,tb_thls,tb_wqs,tb_wts,&
                            tb_w,tb_ug,tb_vg,tb_thl,tb_qt,&
                            tb_uadv,tb_vadv,tb_qtadv,tb_thladv

    implicit none

    character (80):: chmess
    character (1) :: chmess1
    integer :: k,t, ierr
    real :: dummyr
    real, allocatable, dimension (:) :: height

    if (.not. ltimedep) return

    if (ltestbed) then
      kflux = ntnudge
      kls   = ntnudge
    else
      kflux = 1
      kls   = 1
    end if

    allocate(height(k1))
    allocate(timeflux (0:kflux))
    allocate(wqsurft    (kflux))
    allocate(wtsurft    (kflux))
    allocate(timels     (0:kls))
    allocate(ugt       (k1,kls))
    allocate(vgt       (k1,kls))
    allocate(wflst     (k1,kls))
    allocate(thlpcart  (k1,kls))
    allocate(qtproft   (k1,kls))
    allocate(dqtdtlst  (k1,kls))
    allocate(dthldtlst (k1,kls))
    allocate(dudtlst   (k1,kls))
    allocate(dvdtlst   (k1,kls))
    allocate(thlst      (0:kls))
    allocate(qtst       (0:kls))
    allocate(pst        (0:kls))

!
!  Initialize parameters
!
    timels    = 0
    pst       = 0
    qtst      = 0
    thlst     = 0
    ugt       = 0
    vgt       = 0
    wflst     = 0
    wqsurft   = 0
    wtsurft   = 0
    thlpcart  = 0
    qtproft   = 0
    dqtdtlst  = 0
    dthldtlst = 0
    dudtlst   = 0
    dvdtlst   = 0

    if (myid==0) then

!    --- load lsforcings---
      if (ltestbed) then
!
!     Fill the existing arrays for the timestepping with the domian avaraged HARMONIE arrays
!
        write(*,*) 'inittimedep: testbed mode: data for time-dependent forcing obtained from HARMONIE output'
      
        timeflux(1:kflux) = tb_time
        timels  (1:kls  ) = tb_time

        pst      = tb_ps
        qtst     = tb_qts
        thlst    = tb_thls
        wqsurft  = tb_wqs
        wtsurft  = tb_wts
        height  (:) = zf

        do t=1,kls
          ugt      (:,t) = tb_ug    (t,:)
          vgt      (:,t) = tb_vg    (t,:)
          wflst    (:,t) = tb_w     (t,:)
          qtproft  (:,t) = tb_qt    (t,:)
          thlpcart (:,t) = tb_thl   (t,:)
          dqtdtlst (:,t) = tb_qtadv (t,:)
          dthldtlst(:,t) = tb_thladv(t,:)
          dudtlst  (:,t) = tb_uadv  (t,:)
          dvdtlst  (:,t) = tb_vadv  (t,:)
        end do

      else

       wqsurft  = wqsurf
       wtsurft  = wtsurf
       thlst    = thls
       qtst     = qts
       pst      = ps

        open(ifinput,file='ls_flux.inp.'//cexpnr)
        read(ifinput,'(a80)') chmess
        write(6,*) chmess
        read(ifinput,'(a80)') chmess
        write(6,*) chmess
        read(ifinput,'(a80)') chmess
        write(6,*) chmess

        timeflux = 0
        timels   = 0

!      --- load fluxes---
      t    = 0
      ierr = 0
      do while (timeflux(t) < (tres*real(btime)+runtime))
        t=t+1
        read(ifinput,*, iostat = ierr) timeflux(t), wtsurft(t), wqsurft(t),thlst(t),qtst(t),pst(t)
        write(*,'(i8,6e12.4)') t,timeflux(t), wtsurft(t), wqsurft(t),thlst(t),qtst(t),pst(t)
        if (ierr < 0) then
            stop 'STOP: No time dependend data for end of run (surface fluxes)'
        end if
      end do
      if(timeflux(1)>(tres*real(btime)+runtime)) then
         write(6,*) 'Time dependent surface variables do not change before end of'
         write(6,*) 'simulation. --> only large scale forcings'
         ltimedepsurf=.false.
      endif
! flush to the end of fluxlist
      do while (ierr ==0)
        read (ifinput,*,iostat=ierr) dummyr
      end do
      backspace (ifinput)
!     ---load large scale forcings----

       end if   !ltestbed     

!======================================== End surface data =============================================

    end if

    if(ltestbed) then

     call MPI_BCAST(timeflux(1:kflux),kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(wtsurft          ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(wqsurft          ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(thlst            ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(qtst             ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(pst              ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(timels(1:kls)    ,kls,MY_REAL  ,0,comm3d,mpierr)
     call MPI_BCAST(ugt              ,k1*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(vgt              ,k1*kls,MY_REAL,0,comm3d,mpierr)   
     call MPI_BCAST(wflst            ,k1*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(dqtdtlst         ,k1*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(dthldtlst        ,k1*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(thlpcart         ,k1*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(qtproft          ,k1*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(ltimedepsurf ,1,MPI_LOGICAL,0,comm3d,mpierr)
     call MPI_BCAST(ltimedepz    ,1,MPI_LOGICAL,0,comm3d,mpierr)
     call inittimedepsv
     call timedep

    else

     call MPI_BCAST(timeflux(1:kflux),kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(wtsurft          ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(wqsurft          ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(thlst            ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(qtst             ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(pst              ,kflux,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(timels(1:kls)    ,kls,MY_REAL  ,0,comm3d,mpierr)
     call MPI_BCAST(ugt              ,kmax*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(vgt              ,kmax*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(wflst            ,kmax*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(thlpcart,kmax*kls,MY_REAL,0,comm3d,mpierr)
     call MPI_BCAST(qtproft ,kmax*kls,MY_REAL,0,comm3d,mpierr)

     call MPI_BCAST(ltimedepsurf ,1,MPI_LOGICAL,0,comm3d,mpierr)
     call MPI_BCAST(ltimedepz    ,1,MPI_LOGICAL,0,comm3d,mpierr)
     call inittimedepsv
     call timedep

    endif

    deallocate(height)

  end subroutine inittimedep

  subroutine timedep

!-----------------------------------------------------------------|dales_harm.o9510600
!                                                                 |
!*** *timedep*  calculates ls forcings and surface forcings       |
!               case as a funtion of timee                        |
!                                                                 |
!      Roel Neggers    K.N.M.I.     01/05/2001                    |
!                                                                 |
!                                                                 |
!    calls                                                        |
!    * timedepz                                                   |
!      calculation of large scale advection, radiation and        |
!      surface fluxes by interpolation between prescribed         |
!      values at certain times                                    |
!                                                                 |
!    * timedepsurf                                                |
!      calculation  surface fluxes by interpolation               |
!      between prescribed values at certain times                 |
!                                                                 |
!                                                                 |
!-----------------------------------------------------------------|
    use modtimedepsv, only : timedepsv
    implicit none

    if (.not. ltimedep) return
    call timedepz
    call timedepsurf
    call timedepsv
  end subroutine timedep

  subroutine timedepz
    use modfields,   only : ug, vg, wfls,whls,dqtdtls,thlpcar,thl0,dpdxl,dpdyl
    use modglobal,   only : rtimee,om23_gs,dzf,dzh,k1,kmax,llsadv
    use modmpi,      only : myid

    implicit none

    integer t,k
    real fac

    if(.not.(ltimedepz)) return

    !---- interpolate ----
    t=1
    do while(rtimee>timels(t))
      t=t+1
    end do
    if (rtimee/=timels(1)) then
      t=t-1
    end if

    fac = ( rtimee-timels(t) ) / ( timels(t+1)-timels(t) )
    ug      = ugt     (:,t) + fac * ( ugt     (:,t+1) - ugt     (:,t) )
    vg      = vgt     (:,t) + fac * ( vgt     (:,t+1) - vgt     (:,t) )
    wfls    = wflst   (:,t) + fac * ( wflst   (:,t+1) - wflst   (:,t) )
    dqtdtls = dqtdtlst(:,t) + fac * ( dqtdtlst(:,t+1) - dqtdtlst(:,t) )
    thlpcar = thlpcart(:,t) + fac * ( thlpcart(:,t+1) - thlpcart(:,t) )

    do k=1,kmax
      dpdxl(k) =  om23_gs*vg(k)
      dpdyl(k) = -om23_gs*ug(k)
    end do

    whls(1)  = 0.0
    do k=2,kmax
      whls(k) = ( wfls(k)*dzf(k-1) +  wfls(k-1)*dzf(k) )/(2*dzh(k))
    end do
    whls(k1) = (wfls(kmax)+0.5*dzf(kmax)*(wfls(kmax)-wfls(kmax-1)) &
                                                  /dzh(kmax))

  !******include rho if rho = rho(z) /= 1.0 ***********

    if (llsadv) then
      if (myid==0) stop 'llsadv should not be used anymore. Large scale gradients were calculated in a non physical way (and lmomsubs had to be set to true to retain conservation of mass)'
    end if

    return
  end subroutine timedepz

  subroutine timedepsurf
    use modglobal,   only : rtimee, lmoist
    use modsurfdata, only : wtsurf,wqsurf,thls,qts,ps
    use modsurface,  only : qtsurf
    implicit none
    integer t
    real fac

    if(.not.(ltimedepsurf)) return
  !     --- interpolate! ----
    t=1
    do while(rtimee>timeflux(t))
      t=t+1
    end do
    if (rtimee/=timeflux(t)) then
      t=t-1
    endif

    fac = ( rtimee-timeflux(t) ) / ( timeflux(t+1)-timeflux(t))
    wqsurf = wqsurft(t) + fac * ( wqsurft(t+1) - wqsurft(t)  )
    wtsurf = wtsurft(t) + fac * ( wtsurft(t+1) - wtsurft(t)  )
    thls   = thlst(t)   + fac * ( thlst(t+1)   - thlst(t)    )
    ps     = pst(t)     + fac * ( pst(t+1)   - pst(t)    )
!cstep: not necessary to provide qts in ls_flux file qts    = qtst(t)    + fac * ( qtst(t+1)    - qtst(t)     )
    if (lmoist) then
       call qtsurf
    else
       qts = 0.
    endif

    return
  end subroutine timedepsurf


  subroutine exittimedep
    use modtimedepsv, only : exittimedepsv
    implicit none
    if (.not. ltimedep) return
    deallocate(timels,ugt,vgt,wflst,dqtdtlst,thlpcart)
    deallocate(timeflux, wtsurft,wqsurft,thlst,qtst,pst)
    call exittimedepsv

  end subroutine

end module modtimedep
