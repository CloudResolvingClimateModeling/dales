!> \file tstep.f90
!!  Performs the time integration

!>
!!  Performs the time integration
!>
!! Tstep uses adaptive timestepping and 3rd order Runge Kutta time integration.
!! The adaptive timestepping chooses it's delta_t according to the courant number
!! and the cell peclet number, depending on the advection scheme in use.
!!
!!
!!  \author Chiel van Heerwaarden, Wageningen University
!!  \author Thijs Heus,MPI-M
!! \see Wicker and Skamarock 2002
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

!> Determine time step size dt in initialization and update time variables
!!
!! The size of the timestep Delta t is determined adaptively, and is limited by both the Courant-Friedrichs-Lewy criterion CFL
!! \latexonly
!! \begin{equation}
!! \CFL = \mr{max}\left(\left|\frac{u_i \Delta t}{\Delta x_i}\right|\right),
!! \end{equation}
!! and the diffusion number $d$. The timestep is further limited by the needs of other modules, e.g. the statistics.
!! \endlatexonly

module tstep
implicit none
  
contains
subroutine tstep_update


  use modglobal, only : i1,j1,rk3step,timee,rtimee,dtmax,dt,ntimee,ntrun,courant,peclet,dt_reason, &
                        kmax,dx,dy,dzh,dt_lim,ladaptive,timeleft,idtmax,rdt,tres,longint,lwarmstart, a_acc, a_acc1
  use modfields, only : um,vm,wm
  use modsubgrid,only : ekm,ekh
  use modmpi,    only : comm3d,mpierr,mpi_max,my_real
  implicit none

  real, allocatable, dimension (:) :: courtotl,courtot
  integer       :: k
  real,save     :: courtotmax=-1,peclettot=-1
  real          :: courold,peclettotl,pecletold
  logical,save  :: spinup=.true.

    
  allocate(courtotl(kmax),courtot(kmax))

  if(lwarmstart) spinup = .false.

  rk3step = mod(rk3step,3) + 1
  if(rk3step == 1) then
    ! Initialization
    if (spinup) then
      if (ladaptive) then
        courold = courtotmax
        pecletold = peclettot
        peclettotl=0.0
        do k=1,kmax
          courtotl(k)=maxval(um(2:i1,2:j1,k)*um(2:i1,2:j1,k)/(dx*dx)+vm(2:i1,2:j1,k)*vm(2:i1,2:j1,k)/(dy*dy)+&
          wm(2:i1,2:j1,k)*wm(2:i1,2:j1,k)/(dzh(k)*dzh(k)))*rdt*rdt
        end do
        call MPI_ALLREDUCE(courtotl,courtot,kmax,MY_REAL,MPI_MAX,comm3d,mpierr)
        courtotmax=0.0
        do k=1,kmax
          courtotmax=max(courtotmax,courtot(k))
        enddo
        courtotmax=sqrt(courtotmax)
        do k=1,kmax
           ! peclettotl=max(peclettotl,maxval(ekm(2:i1,2:j1,k))*rdt/minval((/dzh(k),dx,dy/))**2)
           ! peclettotl=max(peclettotl,maxval(ekh(2:i1,2:j1,k))*rdt/minval((/dzh(k),dx,dy/))**2) !FJ. ekh > ekm 
           
           ! experimental !
           ! take into account that the grid is non-uniform
           ! assume ekh >= ekm, so we test only ekh
           peclettotl=max(peclettotl, maxval(ekh(2:i1,2:j1,k))*rdt * (1.0/dzh(k)**2 + 1.0/dx**2 + 1.0/dy**2) / 3) 
           ! / 3 is here because peclet is defined as 1/6 now, but we want 1/2 in this formulation

        end do
        call MPI_ALLREDUCE(peclettotl,peclettot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
        if ( pecletold>0) then
           dt = min(timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint))
           dt_reason = minloc((/timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint)/),1)
          if (abs(courtotmax-courold)/courold<0.1 .and. (abs(peclettot-pecletold)/pecletold<0.1)) then
            spinup = .false.
          end if
        end if
        rdt = dble(dt)*tres
        dt_lim = timeleft
        timee   = timee  + dt
        rtimee  = dble(timee)*tres
        timeleft=timeleft-dt
        ntimee  = ntimee + 1
        ntrun   = ntrun  + 1
      else
        dt = 2 * dt
        if (dt >= idtmax) then
          dt = idtmax
          spinup = .false.
        end if
        rdt = dble(dt)*tres
      end if
    ! Normal time loop
    else
      if (ladaptive) then
        peclettotl = 1e-5
        do k=1,kmax
          courtotl(k)=maxval((um(2:i1,2:j1,k)*rdt/dx)*(um(2:i1,2:j1,k)*rdt/dx)+(vm(2:i1,2:j1,k)*rdt/dy)*&
          (vm(2:i1,2:j1,k)*rdt/dy)+(wm(2:i1,2:j1,k)*rdt/dzh(k))*(wm(2:i1,2:j1,k)*rdt/dzh(k)))
        end do
        call MPI_ALLREDUCE(courtotl,courtot,kmax,MY_REAL,MPI_MAX,comm3d,mpierr)
        courtotmax=0.0
        do k=1,kmax
            courtotmax=max(courtotmax,sqrt(courtot(k)))
        enddo
        do k=1,kmax
           !peclettotl=max(peclettotl,maxval(ekm(2:i1,2:j1,k))*rdt/minval((/dzh(k),dx,dy/))**2) 
           !peclettotl=max(peclettotl,maxval(ekh(2:i1,2:j1,k))*rdt/minval((/dzh(k),dx,dy/))**2) !FJ. ekh > ekm 

                      ! experimental !
           ! take into account that the grid is non-uniform
           ! assume ekh >= ekm, so we test only ekh
           peclettotl=max(peclettotl, maxval(ekh(2:i1,2:j1,k))*rdt * (1.0/dzh(k)**2 + 1.0/dx**2 + 1.0/dy**2) / 3) 
           ! / 3 is here because peclet is defined as 1/6 now, but we want 1/2 in this formulation
        end do
        call MPI_ALLREDUCE(peclettotl,peclettot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
        dt = min(timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint))

        ! record which condition is the limiting one
        dt_reason = minloc((/timee,dt_lim,idtmax,floor(rdt/tres*courant/courtotmax,longint),floor(rdt/tres*peclet/peclettot,longint)/), 1)

        ! adjust the acceleration factor for the current step
        ! to ensure that we don't exceed dt_lim
        a_acc1 = a_acc
        if (dt * a_acc > dt_lim) then
           a_acc1 = dble(dt_lim) / dt
        end if

        ! dt and rdt are not adjusted by acceleration - they describe the normal time step
        ! the elapsed time IS adjusted by the acceleration
        
        rdt = dble(dt)*tres
        timeleft=timeleft-dt * a_acc1
        dt_lim = timeleft
        timee   = timee  + dt * a_acc1
        rtimee  = dble(timee)*tres
        ntimee  = ntimee + 1
        ntrun   = ntrun  + 1
      else ! not adaptive
        dt = idtmax
        rdt = dtmax
        ntimee  = ntimee + 1
        ntrun   = ntrun  + 1
        timee   = timee  + dt !ntimee*dtmax
        rtimee  = dble(timee)*tres
        timeleft=timeleft-dt
      end if
    end if
  end if

  deallocate(courtotl,courtot)

end subroutine tstep_update




 ! Slab average a field
  subroutine slab_avg(f, favg)
    use modglobal, only : ih, i1, jh, j1, itot, jtot, kmax
    use modmpi, only: comm3d, my_real, mpi_sum
    implicit none
   
    real, intent(in)   :: f(:,:,:)
    real, intent(out)  :: favg(1:kmax)
    real    :: al(1:kmax)
    integer :: k, mpierr

    do k=1,kmax
      al(k) = sum(f(2-ih:i1+ih,2-jh:j1+jh,k))
    enddo

    call MPI_ALLREDUCE(al, favg, kmax,  MY_REAL, &
                       MPI_SUM, comm3d, mpierr)

    favg = favg * (1.0 / (itot*jtot))
    
  end subroutine slab_avg
  

    

!> Time integration is done by a third order Runge-Kutta scheme.
!!
!! \latexonly
!! With $f^n(\phi^n)$ the right-hand side of the appropriate equation for variable
!! $\phi=\{\fav{u},\fav{v},\fav{w},e^{\smfrac{1}{2}},\fav{\varphi}\}$, $\phi^{n+1}$
!! at $t+\Delta t$ is calculated in three steps:
!! \begin{eqnarray}
!! \phi^{*} &=&\phi^n + \frac{\Delta t}{3}f^n(\phi^n)\nonumber\\\\
!! \phi^{**} &=&\phi^{n} + \frac{\Delta t}{2}f^{*}(\phi^{*})\nonumber\\\\
!! \phi^{n+1} &=&\phi^{n} + \Delta t f^{**}(\phi^{**}),
!! \end{eqnarray}
!! with the asterisks denoting intermediate time steps.
!! \endlatexonly
!! \see Wicker and Skamarock, 2002
subroutine tstep_integrate


  use modglobal, only : ifmessages,i1,j1,kmax,nsv,rdt,rk3step,e12min,ntimee,rtimee,kmax,a_acc1
  use modfields, only : u0,um,up,v0,vm,vp,w0,wm,wp,wp_store,&
                        thl0,thlm,thlp,qt0,qtm,qtp,&
                        e120,e12m,e12p,sv0,svm,svp
  implicit none

  integer i,j,k,n
  real rk3coef
  real up_avg(1:kmax),  vp_avg(1:kmax),  thlp_avg(1:kmax),  qtp_avg(1:kmax),  e12p_avg(1:kmax)

  
  rk3coef = rdt / (4. - dble(rk3step))
  wp_store = wp

  if(rk3step /= 3) then
     u0   = um   + rk3coef * up
     v0   = vm   + rk3coef * vp
     w0   = wm   + rk3coef * wp
     thl0 = thlm + rk3coef * thlp
     qt0  = qtm  + rk3coef * qtp
     !do k=1,kmax
     !   do j=2,j1
     !      do i=2,i1
     !         if (qt0(i,j,k) < 0) then
     !            write(ifmessages,*) 'Warning: qt0(', i, j, k, ') = ', qt0(i,j,k), 'in tstep_integrate. Setting to 1e-6'
     !            write(ifmessages,*) ' rk3step:', rk3step, 'ntimee:', ntimee, 'rtimee:', rtimee
     !            write(ifmessages,*) ' qtp(i,j,k)', qtp(i,j,k), 'rk3coef', rk3coef
     !            qt0(i,j,k) = 1e-6
     !         endif
     !      enddo
     !   enddo
     !enddo
     
     sv0  = svm  + rk3coef * svp
     e120 = max(e12min,e12m + rk3coef * e12p)

  else ! step 3 - store result in both ..0 and ..m
     if (a_acc1 /= 1.0) then
        ! Mean state acceleration, following Jones et al, JAMES 7, 1643 (2015)
        
        call slab_avg(up,   up_avg)
        do k=1,kmax
           um(:,:,k)   = um(:,:,k)   + rk3coef * up(:,:,k)   + rk3coef * (a_acc1 - 1) * up_avg(k)
        enddo
        u0   = um

        call slab_avg(vp,   vp_avg)
        do k=1,kmax
           vm(:,:,k)   = vm(:,:,k)   + rk3coef * vp(:,:,k)   + rk3coef * (a_acc1 - 1) * vp_avg(k)
        enddo
        v0   = vm

        ! slab average of w is 0 - no need to accelerate
        wm   = wm   + rk3coef * wp
        w0   = wm
        
        call slab_avg(thlp, thlp_avg)
        do k=1,kmax
           thlm(:,:,k) = thlm(:,:,k) + rk3coef * thlp(:,:,k) + rk3coef * (a_acc1 - 1) * thlp_avg(k)
        enddo
        thl0 = thlm
        
        call slab_avg(qtp,  qtp_avg)
        do k=1,kmax
           qtm(:,:,k)  = qtm(:,:,k)  + rk3coef * qtp(:,:,k)  + rk3coef * (a_acc1 - 1) * qtp_avg(k)
        enddo
        qt0  = qtm
        
        ! like in [Jones 2015] we don't accelerate rain
        ! NOTE - if scalar values are used for chemicals, they should be accelerated
        svm  = svm  + rk3coef * svp
        sv0  = svm
        
        call slab_avg(e12p, e12p_avg)
        do k=1,kmax
           e12m(:,:,k) = max(e12min, e12m(:,:,k) + rk3coef * e12p(:,:,k) + rk3coef * (a_acc1 - 1) * e12p_avg(k))
        enddo
        e120 = e12m
        
     else ! no mean state acceleration  
        um   = um   + rk3coef * up
        u0   = um
        vm   = vm   + rk3coef * vp
        v0   = vm
        wm   = wm   + rk3coef * wp
        w0   = wm
        thlm = thlm + rk3coef * thlp
        thl0 = thlm
        qtm  = qtm  + rk3coef * qtp
        qt0  = qtm
        svm  = svm  + rk3coef * svp
        sv0  = svm
        e12m = max(e12min,e12m + rk3coef * e12p)
        e120 = e12m
     end if

     
  end if

  up=0.
  vp=0.
  wp=0.
  thlp=0.
  qtp=0.
  svp=0.
  e12p=0.

end subroutine tstep_integrate
end module tstep
