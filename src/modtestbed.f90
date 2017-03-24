!> \file modtestbed.f90
!!  Testbed continuous forcing & nudging
!>

!>
!!  Testbed continuous forcing & nudging
!>
!!  \author Roel Neggers, IGMK
!!  \par Revision list
!!  \todo Documentation
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

module modtestbed

use netcdf

implicit none
PRIVATE
PUBLIC :: inittestbed, testbednudge, exittestbed, ltestbed,testbed_getinttime, &
          tb_time,tb_ps,tb_qts,tb_thls,tb_wqs,tb_wts,tb_z0h, tb_z0m, tb_alb, &
          tb_u,tb_v,tb_w,tb_thl,tb_qt,tb_ug,tb_vg, &
          tb_qtadv,tb_thladv,tb_uadv,tb_vadv, &
          tbrad_p, tbrad_ql, tbrad_qv, tbrad_t, tbrad_o3
SAVE
  real, dimension(:)  , allocatable     :: tb_ps,tb_qts,tb_thls,tb_wts,tb_wqs,tb_time  
  real, dimension(:,:), allocatable     :: tnudge,tb_ug,tb_vg,&
                                          tbrad_p, tbrad_t, tbrad_qv, tbrad_ql, tbrad_o3, &
                                           tb_thl,tb_qt,tb_u,tb_v,tb_w,&
                                           tb_uadv,tb_vadv,tb_qtadv,tb_thladv
  real, dimension(:,:,:),   allocatable :: tb_ps_harm,tb_thls_harm,tb_alb,tb_z0h,tb_z0m, &
                                           tb_wts_harm,tb_wts_t,tb_wqs_harm,tb_qts_harm, &
                                           tb_wqs_eva,tb_wqs_sub,tb_wqs_t
!
!  For reading HARMONIE meteo. fields arrays changed to i_harm,j_harm,k_harm,time
!
  real, dimension(:,:,:,:), allocatable :: tb_thl_k,tb_qt_k,tb_u_k,tb_v_k,tb_w_k, &
                                           tb_ug_k,tb_vg_k, &
                                           tb_uadv_k,tb_vadv_k,tb_qtadv_k,tb_thladv_k, &
                                           tbrad_p_k, tbrad_t_k, tbrad_qv_k, tbrad_ql_k, tbrad_o3_k
  real :: tb_taunudge = 10800.
  logical :: ltestbed = .true., &
             ltb_nudge = .true., &
             ltb_u,ltb_v,ltb_w,ltb_thl,ltb_qt
!
!  Size of input fields output from HARMONIE
!
  integer, parameter :: i_harm=21
  integer, parameter :: j_harm=21
  integer, parameter :: k_harm=65
!
! Harmonie grid size in metres
!
  real,    parameter :: xsize_harm=2500.
  real,    parameter :: ysize_harm=2500.

contains
  subroutine inittestbed

    use modnudge, only :ntnudge,nknudge,nknudgep1
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer
    use modglobal,only :dx,dy,ifnamopt,fname_options,runtime,btime,cexpnr,ifinput,k1,kmax,tres,&
                        grav,rd,cp,pref0,rlv,zf,dzf,dzh,itot,jtot
    use modsurfdata,only : ksoilmax, phifc, phiwp, dzsoil
    use modforces, only : lforce_user

    implicit none
    real, dimension(:), allocatable       :: dumheights
    real, dimension(:,:), allocatable     :: dumo3,dumphiwav,dumswi,dumlwnet,dumswnet,dumphifc,dumphiwp
    real, dimension(:,:,:),   allocatable :: dumphi,dumtsoil
    real, dimension(:,:,:,:), allocatable :: dumheight,dumpf,dumpt,dumt,dumqt, &
                                             dumu,dumv,dumw, dumug,dumvg, &
                                             dumthl,dumqv,dumql,dumqi,dumomega,dumqtadv, &
                                             dumuadv,dumvadv, dumqadv,dumladv,dumiadv,dumtadv,dumthladv

    INTEGER  NCID, STATUS, VARID, timID
    INTEGER start2(2), count2(2)
    character(len = nf90_max_name) :: RecordDimName

    integer :: ierr,i,j,k,t,ik,c,a
    integer :: xsize,ysize,nknudges
    character(1) :: chmess1
    real tv,rho,iexner,fac,dumdens

    namelist /NAMTESTBED/ &
       ltestbed, ltb_nudge, tb_taunudge

    if(myid==0) then

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMTESTBED,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMTESTBED'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMTESTBED'
      endif
      write(6 ,NAMTESTBED)
      close(ifnamopt)

    end if 

    call MPI_BCAST(ltestbed     , 1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ltb_nudge    , 1,MPI_LOGICAL,0,comm3d,mpierr)
    
    if (.not. ltestbed) return
    
    lforce_user = .true.

    if(myid==0) then
        !--- open nc file ---
        STATUS = NF90_OPEN('allvar.Cabauwbox21x21.MSOPier_40h11.2012041300.nc', nf90_nowrite, NCID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

!
!      time from HARMONIE is an array of ntnudge dimensions
!      therefore ntnudge is now declared in modnudge.f90
!
                   
        !--- get time & height dimensions ---
        status = nf90_inq_dimid(ncid, "time", timID)
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_inquire_dimension(NCID, timID, len=ntnudge, name=RecordDimName)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

!
!       number of levels from HARMONIE is declared above  
!
!        status = nf90_inq_dimid(ncid, "levf", timID)
!        if (status /= nf90_noerr) call handle_err(status)
!        status = nf90_inquire_dimension(NCID, timID, len=nknudge, name=RecordDimName)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

 !
 !       The number of soil levels is not declared in HARMONIE file so try 1
 !
 !       status = nf90_inq_dimid(ncid, "nlevp1", timID)
 !       if (status /= nf90_noerr) call handle_err(status)
 !       status = nf90_inquire_dimension(NCID, timID, len=nknudgep1, name=RecordDimName)
 !       if (STATUS .ne. nf90_noerr) call handle_err(STATUS)


!        status = nf90_inq_dimid(ncid, "nlevs", timID)
!        if (status /= nf90_noerr) call handle_err(status)
!        status = nf90_inquire_dimension(NCID, timID, len=nknudges, name=RecordDimName)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
         nknudges=1
    end if

    call MPI_BCAST(ntnudge    , 1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(nknudge    , 1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(nknudgep1  , 1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(nknudges   , 1,MPI_INTEGER,0,comm3d,mpierr)

    !--- allocate space for input variables & reset---
    allocate(    tnudge    (ntnudge,k1), &
                 tb_u_k    (i_harm,j_harm,k1,ntnudge), &
                 tb_u      (ntnudge,k1),               &
                 tb_v_k    (i_harm,j_harm,k1,ntnudge), &
                 tb_v      (ntnudge,k1),               &
                 tb_w_k    (i_harm,j_harm,k1,ntnudge), &
                 tb_w      (ntnudge,k1),               &
                 tb_thl_k  (i_harm,j_harm,k1,ntnudge), &
                 tb_thl    (ntnudge,k1),               &
                 tb_qt_k   (i_harm,j_harm,k1,ntnudge), &
                 tb_qt     (ntnudge,k1),               &
                 tb_ug_k   (i_harm,j_harm,k1,ntnudge), &
                 tb_ug     (ntnudge,k1),               &
                 tb_vg_k   (i_harm,j_harm,k1,ntnudge), &
                 tb_vg     (ntnudge,k1),               &
                 tb_qtadv_k (i_harm,j_harm,k1,ntnudge),&
                 tb_qtadv   (ntnudge,k1),              &
                 tb_thladv_k(i_harm,j_harm,k1,ntnudge),&
                 tb_thladv  (ntnudge,k1),              &
                 tb_uadv_k  (i_harm,j_harm,k1,ntnudge),&
                 tb_uadv    (ntnudge,k1),              &
                 tb_vadv_k  (i_harm,j_harm,k1,ntnudge),&
                 tb_vadv    (ntnudge,k1),              &
                 tb_time    (ntnudge),                  &
                 tb_ps_harm (i_harm,j_harm,ntnudge),    &
                 tb_ps      (ntnudge),                  &
                 tb_qts_harm(i_harm,j_harm,ntnudge),   &
                 tb_qts    (ntnudge),                  &
                 tb_thls_harm(i_harm,j_harm,ntnudge),  &
                 tb_thls   (ntnudge),                  &
                 tb_wts_harm(i_harm,j_harm,ntnudge),   &
                 tb_wts_t  (i_harm,j_harm,ntnudge),    &
                 tb_wts    (ntnudge),                  &
                 tb_wqs_harm(i_harm,j_harm,ntnudge),   &
                 tb_wqs    (ntnudge),                  &
                 tb_wqs_eva(i_harm,j_harm,ntnudge),    &
                 tb_wqs_sub(i_harm,j_harm,ntnudge),    &
                 tb_wqs_t  (i_harm,j_harm,ntnudge),    &
                 tb_z0m    (i_harm,j_harm,ntnudge),    &
                 tb_z0h    (i_harm,j_harm,ntnudge),    &
                 tb_alb    (i_harm,j_harm,ntnudge),    &
                 tbrad_p_k (i_harm,j_harm,k1,ntnudge), &
                 tbrad_p   (ntnudge,k1), &
                 tbrad_t_k (i_harm,j_harm,k1,ntnudge), &
                 tbrad_t   (ntnudge,k1), &
                 tbrad_qv_k(i_harm,j_harm,k1,ntnudge), &
                 tbrad_qv  (ntnudge,k1), &
                 tbrad_ql_k(i_harm,j_harm,k1,ntnudge), &
                 tbrad_ql  (ntnudge,k1), &
                 tbrad_o3  (ntnudge,k1) &
                 )

     tnudge = tb_taunudge     !nudging timescale

        tb_time=0
        tb_ps=0
        tb_qts=0
        tb_thls=0
        tb_wts=0
        tb_wts_t=0
        tb_wqs=0
        tb_wqs_eva=0
        tb_wqs_sub=0
        tb_wqs_t=0
        tb_z0m=0
        tb_z0h=0
        tb_alb=0

        tb_u_k=0
        tb_v_k=0
        tb_w_k=0
        tb_thl_k=0
        tb_qt_k=0
        tb_ug_k=0
        tb_vg_k=0
        tb_qtadv_k=0
        tb_thladv_k=0
        tb_uadv_k=0
        tb_vadv_k=0

        tb_u=0
        tb_v=0
        tb_w=0
        tb_thl=0
        tb_qt=0
        tb_ug=0
        tb_vg=0
        tb_qtadv=0
        tb_thladv=0
        tb_uadv=0
        tb_vadv=0

    if(myid==0) then

        allocate(dumomega (i_harm,j_harm,k_harm,ntnudge), &
                 dumheight(i_harm,j_harm,k_harm,ntnudge), &
                 dumpf    (i_harm,j_harm,k_harm,ntnudge), &
                 dumqv    (i_harm,j_harm,k_harm,ntnudge), &
                 dumql    (i_harm,j_harm,k_harm,ntnudge), &
                 dumqi    (i_harm,j_harm,k_harm,ntnudge), &
                 dumo3    (nknudge,ntnudge), &
                 dumt     (i_harm,j_harm,k_harm,ntnudge), &
                 dumqt    (i_harm,j_harm,k_harm,ntnudge), &
                 dumthl   (i_harm,j_harm,k_harm,ntnudge), &
                 dumu     (i_harm,j_harm,k_harm,ntnudge), &
                 dumv     (i_harm,j_harm,k_harm,ntnudge), &
                 dumw     (i_harm,j_harm,k_harm,ntnudge), &
                 dumug    (i_harm,j_harm,k_harm,ntnudge), &
                 dumvg    (i_harm,j_harm,k_harm,ntnudge), &
                 dumqtadv (i_harm,j_harm,k1,ntnudge),     &
                 dumthladv(i_harm,j_harm,k1,ntnudge),     &
                 dumuadv  (i_harm,j_harm,k1,ntnudge),     &
                 dumvadv  (i_harm,j_harm,k1,ntnudge),     &
                 dumqadv  (i_harm,j_harm,k1,ntnudge),     &
                 dumladv  (i_harm,j_harm,k1,ntnudge),     &
                 dumiadv  (i_harm,j_harm,k1,ntnudge),     & 
                 dumtadv  (i_harm,j_harm,k1,ntnudge),     & 
                 dumswnet (nknudgep1,ntnudge),            & 
                 dumlwnet (nknudgep1,ntnudge),            & 
                 dumheights(nknudges),                    & 
                 dumtsoil (i_harm,j_harm,ntnudge),        &
                 dumphi   (i_harm,j_harm,ntnudge),        & 
                 dumswi   (i_harm,j_harm),                &
                 dumphifc (i_harm,j_harm),                &
                 dumphiwp (i_harm,j_harm)                 &
                 )

        !--- timeseries ---
        !
!        write(6,'(a30,5f10.2)') 'inittestbed: tb_time:',&
!             tb_time(1),tb_time(2),tb_time(3),tb_time(ntnudge-1),tb_time(ntnudge)

        
        STATUS = NF90_INQ_VARID(NCID, 'time', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_time, start=(/1/), count=(/ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

         tb_time=tb_time*(60.*60.*24.)

        !  surface pressure
        STATUS = NF90_INQ_VARID(NCID, 'ps', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_ps_harm, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface temperature
        STATUS = NF90_INQ_VARID(NCID, 'ts', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_thls_harm, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface humidity
        STATUS = NF90_INQ_VARID(NCID, 'huss', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_qts_harm, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface T flux
        STATUS = NF90_INQ_VARID(NCID, 'hfss', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_wts_harm, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !  surface q flux (Evaporation)
        STATUS = NF90_INQ_VARID(NCID, 'hfls_eva', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_wqs_eva, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        !  surface q flux (Sublimation)
        STATUS = NF90_INQ_VARID(NCID, 'hfls_sbl', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_wqs_sub, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

        !
        !  summ surface fluxes
        !         
        tb_wqs_harm=tb_wqs_eva+tb_wqs_sub

        do i=1,i_harm
         do j=1,j_harm
          do t=1,ntnudge
!
!         Surface fluxes in HARMONIE are summed over the entire run therefore convert to a flux/second
!         for use in dales
!
          if(i==1) then
              tb_wts_t(i,j,t)=tb_wts_harm(i,j,t)/(5.*60.)
              tb_wqs_t(i,j,t)=tb_wqs_harm(i,j,t)/(5.*60.)
          endif
          if(i .gt. 1) then
              tb_wts_t(i,j,t)=(tb_wts_harm(i,j,t)-tb_wts_harm(i-1,j,t))/(5.*60.)
              tb_wqs_t(i,j,t)=(tb_wqs_harm(i,j,t)-tb_wqs_harm(i-1,j,t))/(5.*60.)
          endif
           rho = tb_ps_harm(i,j,t) / (rd * tb_thls_harm(i,j,t))
           tb_wts_t(i,j,t) = -tb_wts_t(i,j,t) / (cp  * rho)        !Change sign: upward = positive in LES, but by convention upward = negative in most GCMs.
           tb_wqs_t(i,j,t) = -tb_wqs_t(i,j,t) / (rlv * rho)

           iexner = (tb_ps_harm(i,j,t)/pref0)**(-rd/cp)
           tb_thls_harm(i,j,t) = iexner * tb_thls_harm(i,j,t)
          enddo
         enddo
        enddo

        !  roughness length for momentum
        STATUS = NF90_INQ_VARID(NCID, 'z0', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_z0m, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          
        !  roughness length for heat and moisture
        STATUS = NF90_INQ_VARID(NCID, 'z0h', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_z0h, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
          
        !  surface albedo, for radiation
        STATUS = NF90_INQ_VARID(NCID, 'albl', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, tb_alb, start=(/1,1,1/), count=(/i_harm,j_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 
                   
        !  height
        STATUS = NF90_INQ_VARID(NCID, 'phi', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumheight, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  u 
        STATUS = NF90_INQ_VARID(NCID, 'ua', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumu, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  v
        STATUS = NF90_INQ_VARID(NCID, 'va', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumv, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  qt
        STATUS = NF90_INQ_VARID(NCID, 'hus', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumqv, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        STATUS = NF90_INQ_VARID(NCID, 'clw', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumql, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
        STATUS = NF90_INQ_VARID(NCID, 'cli', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumqi, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )

        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!
!      Humidity summed
!        
        dumqt = dumqv + dumql + dumqi      
          
        !  thl
        STATUS = NF90_INQ_VARID(NCID, 't', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumt, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        STATUS = NF90_INQ_VARID(NCID, 'p', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumpf, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
     
        do i=2,i_harm-1
         do j=2,j_harm-1
          do k=1,k_harm
           do t=1,ntnudge
            iexner = (dumpf(i,j,k,t)/pref0)**(-rd/cp)
        ! liquid potential temperature
            dumthl(i,j,k,t) = iexner*(dumt(i,j,k,t) - (rlv * (dumql(i,j,k,t)+dumqi(i,j,k,t))) / cp)
            dumdens = dumpf(i,j,k,t)/(dumt(i,j,k,t)*287.05)
            if(i>1 .or. i<i_harm) dumug(i,j,k,t) = (-1.0/(dumdens*1.0e-4))*((dumpf(i+1,j,k,t)-dumpf(i-1,j,k,t))/(2.*xsize_harm))
            if(j>1 .or. i<j_harm) dumvg(i,j,k,t) = (1.0/(dumdens*1.0e-4))*((dumpf(i,j+1,k,t)-dumpf(i,j-1,k,t))/(2.*xsize_harm))
          enddo
         enddo
        enddo
       enddo

        !  Ozone
!        STATUS = NF90_INQ_VARID(NCID, 'o3', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumo3, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        !  w
        STATUS = NF90_INQ_VARID(NCID, 'w', VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, dumomega, start=(/1,1,1,1/), count=(/i_harm,j_harm,k_harm,ntnudge/) )
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !
        ! Conversion of omega
        !        
        do i=2,i_harm-1
         do j=2,j_harm-1
          do k=1,k_harm
           do t=1,ntnudge
            tv  = dumt(i,j,k,t) * (1.+0.61*dumqv(i,j,k,t))
            rho = dumpf(i,j,k,t) / (rd*tv)
            dumw(i,j,k,t) = - dumomega(i,j,k,t) / ( rho * grav )     !convert from Pa/s to m/s
           enddo
          enddo
         enddo
        enddo

        ! geostrophic wind
        !  ug
!        STATUS = NF90_INQ_VARID(NCID, 'ug', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumug, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  vg
!        STATUS = NF90_INQ_VARID(NCID, 'vg', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumvg, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)          tb_uadv_k
         
        !  uadv
!
!       HARMONIE tendencies are not currently available
!
!        STATUS = NF90_INQ_VARID(NCID, 'uadv', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumuadv, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !  vadv
!        STATUS = NF90_INQ_VARID(NCID, 'vadv', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumvadv, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)        
        !  qtadv
!        STATUS = NF90_INQ_VARID(NCID, 'qadv', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumqadv, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)lscale.inp.002
        
!        STATUS = NF90_INQ_VARID(NCID, 'ladv', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumladv, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        
!        STATUS = NF90_INQ_VARID(NCID, 'iadv', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumiadv, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

!
!   Calculate the dynamical tendencies where dx/dy are grid cell dimensions
!
        do i=2,i_harm-1
         do j=2,j_harm-1
          do k=1,k_harm
           do t=1,ntnudge-1
            dumuadv(i,j,k,t)=(dumu(i+1,j,k,t)-dumu(i-1,j,k,t))/(2.*xsize_harm)
            dumvadv(i,j,k,t)=(dumv(i+1,j,k,t)-dumv(i-1,j,k,t))/(2.*ysize_harm)      
            dumqadv(i,j,k,t)=(dumqv(i+1,j,k,t)-dumqv(i-1,j,k,t))/(2.*xsize_harm) 
            dumladv(i,j,k,t)=(dumql(i+1,j,k,t)-dumql(i-1,j,k,t))/(2.*xsize_harm)
            dumiadv(i,j,k,t)=(dumqi(i+1,j,k,t)-dumqi(i-1,j,k,t))/(2.*xsize_harm)
            dumqtadv(i,j,k,t) = dumqadv(i,j,k,t) + dumladv(i,j,k,t) + dumiadv(i,j,k,t)
           enddo
          enddo
         enddo
        enddo
!
!       Calculate the tendency for potential temperature using the tendency for T 
!                  
        do i=2,i_harm-1
         do j=2,j_harm-1
          do k=1,k_harm
           do t=1,ntnudge-1
            dumtadv(i,j,k,t)=(dumt(i+1,j,k,t)-dumt(i-1,j,k,t))/(2.*xsize_harm)
            iexner = (dumpf(i,j,k,t)/pref0)**(-rd/cp)
            dumthladv(i,j,k,t) = dumtadv(i,j,k,t) * iexner
           enddo
          enddo
         enddo
        enddo
          
        !  thladv
 !       STATUS = NF90_INQ_VARID(NCID, 'tadv', VARID)
 !       if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
 !       STATUS = NF90_GET_VAR (NCID, VARID, dumtadv, start=start2, count=count2)
 !       if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
                
        !  net SW downward flux
!        STATUS = NF90_INQ_VARID(NCID, 'fradSWnet', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumswnet, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
        !  net LW downward flux
!        STATUS = NF90_INQ_VARID(NCID, 'fradLWnet', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumlwnet, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
          
!        do i=1,ntnudge
!          tb_Qnet(i) = dumswnet(nknudgep1,i) + dumlwnet(nknudgep1,i)      !flux at surface is stored in lowest half level of profile
!          write(6,*) "modtestbed: qnet:",i,tb_Qnet(i),nknudge,nknudgep1allvar.Cabauwbox21x21.MSOPier_40h11.2012041300.nc
!        enddo

        !--- soil profiles ---
        !
        ! 3D soil profiles are not available from the HARMONIE output
        ! 

!        STATUS = NF90_INQ_VARID(NCID, 'h_soil', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumheights, start=(/1/), count=(/nknudges/) )
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS) 

!        start2 = (/ 1       , 1       /)
!        count2 = (/ nknudges, ntnudge /)

        !  tsoilav
!        STATUS = NF90_INQ_VARID(NCID, 'ts', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumtsoil, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)

        !  phiwav
!        STATUS = NF90_INQ_VARID(NCID, 'q_soil', VARID)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
!        STATUS = NF90_GET_VAR (NCID, VARID, dumphiwav, start=start2, count=count2)
!        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)tb_uadv_k
        
        !  field capacity
        status = NF90_INQ_VARID(NCID, 'wfc1', VARID)
        if (status /= nf90_noerr) call handle_err(status)
        status = NF90_GET_VAR (NCID, VARID, dumphifc, start=(/1,1/), count=(/i_harm,j_harm/) )
        if (status /= nf90_noerr) call handle_err(status)
        
        !  wilting point
        status = NF90_INQ_VARID(ncid, 'wwilt1', VARID)
        if (status /= nf90_noerr) call handle_err(status)
        status = NF90_GET_VAR (NCID, VARID, dumphiwp, start=(/1,1/), count=(/i_harm,j_harm/) )
        if (status /= nf90_noerr) call handle_err(status)

 !       dumswi = ( dumphiwav - dumphiwp ) / ( dumphifc - dumphiwp )   !soil wetness index, using input values for wilting point and field capacity
                 
        !--- close nc file ---
        STATUS = NF90_CLOSE(NCID)  
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        !==============================================================================================	
        !
        !--- Interpolate onto DALES levels (expanded from original code)
        !
        !    Convert geopotential height to m to allow a valid vertical interpolation
        !
        !=============================================================================================

        dumheight=dumheight/grav
   
        do t=1,ntnudge
         do i=1,i_harm
          do j=1,j_harm

          !--- interpolate towards LES levels, reverse height-order, switch dimensions ---
          ik = nknudge
           do k=1,k1
            
            do while( zf(k).gt.dumheight(i,j,ik,t) .and. ik.gt.1)
              ik=ik-1
            enddo
            if ( ik.lt.nknudge ) then 
              ik=ik+1  
            endif

            fac = ( zf(k)-dumheight(i,j,ik,t) ) / ( dumheight(i,j,ik-1,t)-dumheight(i,j,ik,t) )
       
            tb_thl_k   (i,j,k,t)  = dumthl   (i,j,ik,t) + fac * ( dumthl   (i,j,ik-1,t)    - dumthl   (i,j,ik,t))
            tb_qt_k    (i,j,k,t)  = dumqt    (i,j,ik,t) + fac * ( dumqt    (i,j,ik-1,t)    - dumqt    (i,j,ik,t))
            tb_u_k     (i,j,k,t)  = dumu     (i,j,ik,t) + fac * ( dumu     (i,j,ik-1,t)    - dumu     (i,j,ik,t))
            tb_v_k     (i,j,k,t)  = dumv     (i,j,ik,t) + fac * ( dumv     (i,j,ik-1,t)    - dumv     (i,j,ik,t))
            tb_w_k     (i,j,k,t)  = dumw     (i,j,ik,t) + fac * ( dumw     (i,j,ik-1,t)    - dumw     (i,j,ik,t))
            tb_ug_k    (i,j,k,t)  = dumug    (i,j,ik,t) + fac * ( dumug    (i,j,ik-1,t)    - dumug    (i,j,ik,t))
            tb_vg_k    (i,j,k,t)  = dumvg    (i,j,ik,t) + fac * ( dumvg    (i,j,ik-1,t)    - dumvg    (i,j,ik,t))
            tb_uadv_k  (i,j,k,t)  = dumuadv  (i,j,ik,t) + fac * ( dumuadv  (i,j,ik-1,t)    - dumuadv  (i,j,ik,t)) 
            tb_vadv_k  (i,j,k,t)  = dumvadv  (i,j,ik,t) + fac * ( dumvadv  (i,j,ik-1,t)    - dumvadv  (i,j,ik,t))
            tb_qtadv_k (i,j,k,t)  = dumqtadv (i,j,ik,t) + fac * ( dumqtadv (i,j,ik-1,t)    - dumqtadv (i,j,ik,t))
            tb_thladv_k(i,j,k,t)  = dumthladv(i,j,ik,t) + fac * ( dumthladv(i,j,ik-1,t)    - dumthladv(i,j,ik,t))
            tbrad_p_k  (i,j,k,t)  = dumpf    (i,j,ik,t) + fac * ( dumpf    (i,j,ik-1,t)    - dumpf    (i,j,ik,t))
            tbrad_t_k  (i,j,k,t)  = dumt     (i,j,ik,t) + fac * ( dumt     (i,j,ik-1,t)    - dumt     (i,j,ik,t))
            tbrad_qv_k (i,j,k,t)  = dumqv    (i,j,ik,t) + fac * ( dumqv    (i,j,ik-1,t)    - dumqv    (i,j,ik,t))
            tbrad_ql_k (i,j,k,t)  = (dumql   (i,j,ik,t) + fac * ( dumql    (i,j,ik-1,t)    - dumql    (i,j,ik,t)))
            tbrad_ql_k (i,j,k,t)  = tbrad_ql_k (i,j,k,t) + (dumqi    (i,j,ik,t) + fac * ( dumqi     (i,j,ik-1,t)   - dumqi    (i,j,ik,t)))
           enddo
          enddo
        enddo
      enddo   
      
        !--- clean-up ---
        deallocate(dumomega)
        deallocate(dumqv)
        deallocate(dumql)
        deallocate(dumqi)
        deallocate(dumt)
        deallocate(dumpf)
        deallocate(dumo3)
	
       ! deallocate(dumheight)
        deallocate(dumqt)
        deallocate(dumthl)
        deallocate(dumu)
        deallocate(dumv)
        deallocate(dumw)
        deallocate(dumug)
        deallocate(dumvg)
	
        deallocate(dumuadv)
        deallocate(dumvadv)
        deallocate(dumqtadv)
        deallocate(dumthladv)
        deallocate(dumqadv)
        deallocate(dumladv)
        deallocate(dumiadv)
        deallocate(dumtadv)
        deallocate(dumswnet)
        deallocate(dumlwnet)
	
        deallocate(dumheights)
        deallocate(dumtsoil)
        deallocate(dumswi)
!
!      In the first instance fill the DALES grid with the domain avergaed HARMONIE values 
!      without interpolation in the x- or y- directions.
!          
         do t=1,ntnudge
            do k=1,k1
                 tb_thl   (t,k)  = sum(tb_thl_k   (:,:,k,t))/(i_harm*j_harm)
                 tb_qt    (t,k)  = sum(tb_qt_k    (:,:,k,t))/(i_harm*j_harm)
                 tb_u     (t,k)  = sum(tb_u_k     (:,:,k,t))/((i_harm-2)*(j_harm-2))
                 tb_v     (t,k)  = sum(tb_v_k     (:,:,k,t))/((i_harm-2)*(j_harm-2))
                 tb_w     (t,k)  = sum(tb_w_k     (:,:,k,t))/((i_harm-2)*(j_harm-2))
                 tb_ug    (t,k)  = sum(tb_ug_k    (:,:,k,t))/((i_harm-2)*(j_harm-2))	 
                 tb_vg    (t,k)  = sum(tb_vg_k    (:,:,k,t))/((i_harm-2)*(j_harm-2))	 
                 tb_uadv  (t,k)  = sum(tb_uadv_k  (:,:,k,t))/(i_harm*j_harm)		 
                 tb_vadv  (t,k)  = sum(tb_vadv_k  (:,:,k,t))/(i_harm*j_harm)	 
                 tb_qtadv (t,k)  = sum(tb_qtadv_k (:,:,k,t))/(i_harm*j_harm)
                 tb_thladv(t,k)  = sum(tb_thladv_k(:,:,k,t))/(i_harm*j_harm)
                 tbrad_p  (t,k)  = sum(tbrad_p_k  (:,:,k,t))/(i_harm*j_harm)
                 tbrad_t  (t,k)  = sum(tbrad_t_k  (:,:,k,t))/(i_harm*j_harm)			 
                 tbrad_qv (t,k)  = sum(tbrad_qv_k (:,:,k,t))/(i_harm*j_harm)
                 tbrad_ql (t,k)  = sum(tbrad_ql_k (:,:,k,t))/(i_harm*j_harm)
                 tbrad_ql (t,k)  = sum(tbrad_ql_k (:,:,k,t))/(i_harm*j_harm)                          
            enddo
                 tb_ps(t)   = sum(tb_ps_harm(:,:,t))/(i_harm*j_harm)
                 ! Convert to Pa
                 !tb_ps(t)   = tb_ps(t)*100.
                 tb_thls(t) = sum(tb_thls_harm(:,:,t))/(i_harm*j_harm)
                 tb_qts(t)  = sum(tb_qts_harm(:,:,t))/(i_harm*j_harm)
                 tb_wts(t)  = sum(tb_wts_t(:,:,t))/(i_harm*j_harm)
                 tb_wqs(t)  = sum(tb_wqs_t(:,:,t))/(i_harm*j_harm)
          enddo
       
!	  print*,'first column values at time-0'
!	  do k=k1,1,-1
!	    print*,k,sum(dumheight(:,:,k,1))/(i_harm*j_harm),'tb_ug ; ',tb_ug(1,k)
!	  enddo
!	  print*,'++++++++++++++++++++++++++++++++++++++'
!	  do k=k1,1,-1
!	    print*,k,sum(dumheight(:,:,k,1))/(i_harm*j_harm),'tb_vg ; ',tb_vg(1,k)
!	  enddo
!	  print*,'++++++++++++++++++++++++++++++++++++++'
!	  do k=k1,1,-1
!	    print*,k,sum(dumheight(:,:,k,1))/(i_harm*j_harm),'tb_u ; ',tb_u(1,k)
!	  enddo
!	  print*,'++++++++++++++++++++++++++++++++++++++'
!	  do k=k1,1,-1
!	    print*,k,sum(dumheight(:,:,k,1))/(i_harm*j_harm),'tb_v ; ',tb_v(1,k)
!	  enddo

          deallocate(tb_thl_k)
          deallocate(tb_qt_k)
          deallocate(tb_u_k)
          deallocate(tb_v_k)
          deallocate(tb_w_k)
          deallocate(tb_uadv_k)
          deallocate(tb_vadv_k)
          deallocate(tb_qtadv_k)
          deallocate(tb_thladv_k)
          deallocate(tbrad_p_k)
          deallocate(tbrad_t_k)
          deallocate(tbrad_qv_k)
          deallocate(tbrad_ql_k)
        
          deallocate(tb_ps_harm)
          deallocate(tb_thls_harm)
          deallocate(tb_qts_harm)
          deallocate(tb_wts_harm)
          deallocate(tb_wts_t)
          deallocate(tb_wqs_harm)
          deallocate(tb_wqs_t)

        !--- do some output to screen ---
!        do i=1,2
!        !do i=1,ntnudge
!
!        write(6,'(a20,f10.2,a15,3f10.2)') 'modtestbed: scm_in time:',tb_time(i),' sfc pressure:',tb_ps(i),tb_thls(i),tb_qts(i)
!
!        write(6,*) ' zf       tnudge    tb_u    tb_v    tb_w    tb_thl    tb_qt    tb_ug    tb_vg'
!        do k=kmax,1,-1
!          write (6,'(f7.1,8e12.4)') &
!                zf          (k), &
!                tnudge      (i,k), &
!                tb_u        (i,k), &
!                tb_v        (i,k), &
!                tb_w        (i,k), &
!                tb_thl      (i,k), &
!                tb_qt       (i,k), &
!                tb_ug       (i,k), &
!                tb_vg       (i,k)
!        end do
!
!        end do

    end if


    call MPI_BCAST(ntnudge    , 1,MPI_INTEGER,0,comm3d,mpierr)

    call MPI_BCAST(tb_time    ,ntnudge   ,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_ps      ,ntnudge   ,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_qts     ,ntnudge   ,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_thls    ,ntnudge   ,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_wts     ,ntnudge   ,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_wqs     ,ntnudge   ,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_z0h     ,ntnudge   ,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_z0m     ,ntnudge   ,MY_REAL    ,0,comm3d,mpierr)

    call MPI_BCAST(tnudge     ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_u       ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_v       ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_w       ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_thl     ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_qt      ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_ug      ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_vg      ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_uadv    ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_vadv    ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_qtadv   ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tb_thladv  ,ntnudge*k1,MY_REAL    ,0,comm3d,mpierr)

    call MPI_BCAST(tbrad_p      ,ntnudge*nknudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tbrad_qv     ,ntnudge*nknudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tbrad_ql     ,ntnudge*nknudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tbrad_t      ,ntnudge*nknudge,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(tbrad_o3     ,ntnudge*nknudge,MY_REAL    ,0,comm3d,mpierr)

    ltb_u   = any(abs(tb_u)>1e-8)
    ltb_v   = any(abs(tb_v)>1e-8)
    ltb_w   = any(abs(tb_w)>1e-8)
    ltb_thl = any(abs(tb_thl)>1e-8)
    ltb_qt  = any(abs(tb_qt)>1e-8)

  end subroutine inittestbed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testbednudge
    use modglobal, only : rtimee,i1,j1,k1,kmax,rdt
    use modfields, only : up,vp,wp,thlp, qtp,u0av,v0av,qt0av,thl0av
    use modmpi,    only : myid
    implicit none

    integer k,t
    real :: dtm,dtp,currtnudge, qttnudge,qtthres

    if (.not.(ltestbed .and. ltb_nudge)) return

    if (rtimee==0) return

    t=1
    do while(rtimee>tb_time(t))
      t=t+1
    end do
    if (rtimee/=tb_time(1)) then
      t=t-1
    end if

    dtm = ( rtimee-tb_time(t) ) / ( tb_time(t+1)-tb_time(t) )
    dtp = ( tb_time(t+1)-rtimee)/ ( tb_time(t+1)-tb_time(t) )

    qtthres = 1e-6
    do k=1,kmax

      currtnudge = max(rdt,tnudge(k,t)*dtp+tnudge(k,t+1)*dtm)

      if (ltb_u)   up(2:i1,2:j1,k) = up(2:i1,2:j1,k)     - &
          ( u0av(k)   - (tb_u(t,k)  *dtp + tb_u(t+1,k)  *dtm) ) / currtnudge

      if (ltb_v)   vp(2:i1,2:j1,k) = vp(2:i1,2:j1,k)     - &
          ( v0av(k)   - (tb_v(t,k)  *dtp + tb_v(t+1,k)  *dtm) ) / currtnudge

      if (ltb_w)   wp(2:i1,2:j1,k) = wp(2:i1,2:j1,k)     - &
          (           - (tb_w(t,k)  *dtp + tb_w(t+1,k)  *dtm) ) / currtnudge

      if (ltb_thl) thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k) - &
          ( thl0av(k) - (tb_thl(t,k)*dtp + tb_thl(t+1,k)*dtm) ) / currtnudge

      if (ltb_qt)  then
        if (qt0av(k)< qtthres) then
          qttnudge = rdt
        else
          qttnudge = currtnudge
        end if
          qtp(2:i1,2:j1,k) = qtp(2:i1,2:j1,k)   - &
            ( qt0av(k)  - (tb_qt(t,k) *dtp + tb_qt(t+1,k) *dtm) ) / qttnudge
      end if
    end do

    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, qt0av (1),tb_qt (t,1),tb_qt (t+1,1)
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, qt0av (kmax),tb_qt (t,kmax),tb_qt (t+1,kmax)
    
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, thl0av(1),tb_thl(t,1),tb_thl(t+1,1)
    !write(6,*) 'testbednudge:', rtimee, t, tb_time(t), tb_time(t+1), currtnudge, dtm, dtp, thl0av(kmax),tb_thl(t,kmax),tb_thl(t+1,kmax)

  end subroutine testbednudge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine testbed_getinttime(t, dtm, dtp)
     use modglobal, only : rtimee
!     use modfields, only : up,vp,wp,thlp, qtp,u0av,v0av,qt0av,thl0av
!     use modmpi,    only : myid
    implicit none
    integer, intent(out) :: t
    real, intent(out)    :: dtm, dtp


    t=1
    do while(rtimee>tb_time(t))
      t=t+1
    end do
    if (rtimee/=tb_time(1)) then
      t=t-1
    end if

    dtm = ( rtimee-tb_time(t) ) / ( tb_time(t+1)-tb_time(t) )
    dtp = ( tb_time(t+1)-rtimee)/ ( tb_time(t+1)-tb_time(t) )


  end subroutine testbed_getinttime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exittestbed
  deallocate(tb_time)
  end subroutine exittestbed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
  subroutine handle_err(errcode)
      
  implicit none

  integer errcode
     
  write(6,*) 'Error: ', nf90_strerror(errcode)
  stop 2
      
  end subroutine handle_err


end module modtestbed
