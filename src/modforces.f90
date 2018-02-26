!> \file modforces.f90
!!  Calculates the other forces and sources in the equations.

!>
!!  Calculates the other forces and sources in the equations.
!>
!!  This includes the large scale forcings, the coriolis and the subsidence
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \author Steef Böing, TU Delft
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



module modforces
!Calculates additional forces and large scale tendencies
implicit none
save
private
public :: forces, coriolis, lstend,lforce_user
logical :: lforce_user = .false.
contains
  subroutine forces

!-----------------------------------------------------------------|
!                                                                 |
!      Hans Cuijpers   I.M.A.U.                                   |
!      Pier Siebesma   K.N.M.I.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Calculates all other terms in the N-S equation,            |
!      except for the diffusion and advection terms.              |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *forces* is called from *program*.                          |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,dzf,grav
  use modfields, only : sv0,up,vp,wp,thv0h,dpdxl,dpdyl,thvh
  use moduser,   only : force_user
  use modmicrodata, only : imicro, imicro_bulk, imicro_bin, imicro_sice,iqr
  implicit none

  integer i, j, k, jm, jp, km, kp

  if (lforce_user) call force_user

  if((imicro==imicro_sice).or.(imicro==imicro_bulk).or.(imicro==imicro_bin)) then
    do k=2,kmax
      kp=k+1
      km=k-1
    do j=2,j1
      jp=j+1
      jm=j-1
    do i=2,i1
      up(i,j,k) = up(i,j,k) - dpdxl(k)
      vp(i,j,k) = vp(i,j,k) - dpdyl(k)
      wp(i,j,k) = wp(i,j,k) + grav*(thv0h(i,j,k)-thvh(k))/thvh(k) - &
                  grav*(sv0(i,j,k,iqr)*dzf(k-1)+sv0(i,j,k-1,iqr)*dzf(k))/(2.0*dzh(k))
    end do
    end do
    end do
  else
    do k=2,kmax
      kp=k+1
      km=k-1
    do j=2,j1
      jp=j+1
      jm=j-1
    do i=2,i1
      up(i,j,k) = up(i,j,k) - dpdxl(k)
      vp(i,j,k) = vp(i,j,k) - dpdyl(k)
      wp(i,j,k) = wp(i,j,k) + grav*(thv0h(i,j,k)-thvh(k))/thvh(k)
    end do
    end do
    end do
  end if

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------

  do j=2,j1
    jp = j+1
    jm = j-1
  do i=2,i1

    up(i,j,1) = up(i,j,1) - dpdxl(1)

    vp(i,j,1) = vp(i,j,1) - dpdyl(1)

    wp(i,j,1) = 0.0

  end do
  end do
!     ----------------------------------------------end i,j-loop


  return
  end subroutine forces
  subroutine coriolis

!-----------------------------------------------------------------|
!                                                                 |
!      Thijs Heus TU Delft                                        |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Calculates the Coriolis force.                             |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *coriolis* is called from *program*.                          |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,dzf,cu,cv,om22,om23
  use modfields, only : u0,v0,w0,up,vp,wp
  implicit none

  integer i, j, k, jm, jp, km, kp


  do k=2,kmax
    kp=k+1
    km=k-1
  do j=2,j1
    jp=j+1
    jm=j-1
  do i=2,i1

    up(i,j,k) = up(i,j,k)+ cv*om23 &
          +(v0(i,j,k)+v0(i,jp,k)+v0(i-1,j,k)+v0(i-1,jp,k))*om23*0.25 &
          -(w0(i,j,k)+w0(i,j,kp)+w0(i-1,j,kp)+w0(i-1,j,k))*om22*0.25

    vp(i,j,k) = vp(i,j,k)  - cu*om23 &
          -(u0(i,j,k)+u0(i,jm,k)+u0(i+1,jm,k)+u0(i+1,j,k))*om23*0.25


    wp(i,j,k) = wp(i,j,k) + cu*om22 +( (dzf(km) * (u0(i,j,k)  + u0(i+1,j,k) )    &
                +    dzf(k)  * (u0(i,j,km) + u0(i+1,j,km))  ) / dzh(k) ) &
                * om22*0.25
  end do
  end do
!     -------------------------------------------end i&j-loop
  end do
!     -------------------------------------------end k-loop

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------

  do j=2,j1
    jp = j+1
    jm = j-1
  do i=2,i1

    up(i,j,1) = up(i,j,1)  + cv*om23 &
          +(v0(i,j,1)+v0(i,jp,1)+v0(i-1,j,1)+v0(i-1,jp,1))*om23*0.25 &
          -(w0(i,j,1)+w0(i,j ,2)+w0(i-1,j,2)+w0(i-1,j ,1))*om22*0.25

    vp(i,j,1) = vp(i,j,1) - cu*om23 &
          -(u0(i,j,1)+u0(i,jm,1)+u0(i+1,jm,1)+u0(i+1,j,1))*om23*0.25

    wp(i,j,1) = 0.0

  end do
  end do
!     ----------------------------------------------end i,j-loop


  return
  end subroutine coriolis

! experimental version - get rid of dudxls etc, vectorize  
subroutine lstend_2

!-----------------------------------------------------------------|
!                                                                 |
!*** *lstend*  calculates large-scale tendencies                  |
!                                                                 |
!      Pier Siebesma   K.N.M.I.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     calculates and adds large-scale tendencies due to           |
!     large scale advection and subsidence.                       |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *lstend* is called from *program*.                  |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,nsv,lmomsubs
  use modfields, only : up,vp,thlp,qtp,svp,&
                        whls, u0av,v0av,thl0,qt0,sv0,u0,v0,&
                        dudxls,dudyls,dvdxls,dvdyls,dthldxls,dthldyls,dqtdxls,dqtdyls,dqtdtls
  implicit none

  integer i,j,k,n,kp,km
  real subs_thl,subs_qt,subs_sv,subs_u,subs_v

!     1. DETERMINE LARGE SCALE TENDENCIES
!        --------------------------------

!     1.1 lowest model level above surface : only downward component

  subs_thl = 0.
  subs_qt  = 0.
  subs_sv  = 0.
  subs_u   = 0.
  subs_v   = 0.

  k = 1
  if (whls(2).lt.0) then !neglect effect of mean ascending on tendencies at the lowest full level
     thlp(:,:,1) = thlp(:,:,1) - 0.5*whls(2)  *(thl0(:,:,2)-thl0(:,:,1))/dzh(2)
     qtp(:,:,1)  = qtp (:,:,1) - 0.5*whls(2)  *(qt0(:,:,2)-qt0(:,:,1) )/dzh(2) + dqtdtls(1)
     if (lmomsubs) then
        up  (:,:,1) = up  (:,:,1) - 0.5*whls(2)  *(u0(:,:,2)-u0(:,:,1))/dzh(2)
        vp  (:,:,1) = vp  (:,:,1) - 0.5*whls(2)  *(v0(:,:,2)-v0(:,:,1))/dzh(2)
     endif
     svp(:,:,1,:) = svp(:,:,1,:) - 0.5*whls(2)  *(sv0(:,:,2,:)-sv0(:,:,1,:)  )/dzh(2)
  else
     qtp(:,:,1)  = qtp (:,:,1) + dqtdtls(1) ! dqtdtls(1) should be included in any case
  endif
      


!     1.2 other model levels twostream

  do k=2,kmax
    kp=k+1
    km=k-1
    if (whls(kp).lt.0) then   !downwind scheme for subsidence
       thlp(:,:,k) = thlp(:,:,k) - whls(kp) * (thl0(:,:,kp) - thl0(:,:,k))/dzh(kp)
       qtp (:,:,k) = qtp (:,:,k) - whls(kp) * (qt0 (:,:,kp) - qt0 (:,:,k))/dzh(kp) + dqtdtls(k)
       if (lmomsubs) then
          up  (:,:,k) = up  (:,:,k) - whls(kp) * (u0(:,:,kp) - u0(:,:,k))/dzh(kp)
          vp  (:,:,k) = vp  (:,:,k) - whls(kp) * (v0(:,:,kp) - v0(:,:,k))/dzh(kp)
       endif
       svp(:,:,k,:) = svp(:,:,k,:) - whls(kp)  *(sv0(:,:,kp,:) - sv0(:,:,k,:))/dzh(kp)
    else !downwind scheme for mean upward motions
       thlp(:,:,k) = thlp(:,:,k) - whls(k) * (thl0(:,:,k) - thl0(:,:,km))/dzh(k)
       qtp (:,:,k) = qtp (:,:,k) - whls(k) * (qt0 (:,:,k) - qt0 (:,:,km))/dzh(k) + dqtdtls(k)
       if (lmomsubs) then
          up  (:,:,1) = up  (:,:,1) - whls(k) * (u0(:,:,k) - u0(:,:,km))/dzh(k)
          vp  (:,:,1) = vp  (:,:,1) - whls(k) * (v0(:,:,k) - v0(:,:,km))/dzh(k)
       endif
       svp(:,:,k,:) = svp(:,:,k,:) - whls(k) * (sv0(:,:,k,:) - sv0(:,:,km,:))/dzh(k)
    endif
  enddo

  return
end subroutine lstend_2


! experimental version - get rid of dudxls etc. Divisions out of loop.
! only dqtdtls is left.
subroutine lstend

!-----------------------------------------------------------------|
!                                                                 |
!*** *lstend*  calculates large-scale tendencies                  |
!                                                                 |
!      Pier Siebesma   K.N.M.I.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     calculates and adds large-scale tendencies due to           |
!     large scale advection and subsidence.                       |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *lstend* is called from *program*.                  |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,k1,kmax,dzh,nsv,lmomsubs
  use modfields, only : up,vp,thlp,qtp,svp,&
                        whls, u0av,v0av,thl0,qt0,sv0,u0,v0,&
                        dqtdtls
  implicit none

  integer i,j,k,n,kp,km
  real, dimension(1:k1) :: whls_dzh
  real subs_thl,subs_qt,subs_sv,subs_u,subs_v

  whls_dzh(1:k1) = whls(1:k1) / dzh(1:k1)
  ! pre-calculate to get division out of loop
  ! gfortran doesn't extract it by itself, probably since the access is at k or kp
  ! in the downwind scheme
  
!     1. DETERMINE LARGE SCALE TENDENCIES
!        --------------------------------

!     1.1 lowest model level above surface : only downward component

  
  do j=2,j1
    do i=2,i1
      k = 1
      if (whls_dzh(2).lt.0) then !neglect effect of mean ascending on tendencies at the lowest full level
        subs_thl     = 0.5*whls_dzh(2)  *(thl0(i,j,2)-thl0(i,j,1))
        subs_qt      = 0.5*whls_dzh(2)  *(qt0(i,j,2)-qt0(i,j,1) )
        thlp(i,j,1) = thlp(i,j,1) -subs_thl
        qtp(i,j,1)  = qtp (i,j,1) -subs_qt + dqtdtls(1)

        if (lmomsubs) then
          subs_u     = 0.5*whls_dzh(2)  *(u0(i,j,2)-u0(i,j,1))
          subs_v     = 0.5*whls_dzh(2)  *(v0(i,j,2)-v0(i,j,1))
          up  (i,j,1) = up  (i,j,1) -subs_u
          vp  (i,j,1) = vp  (i,j,1) -subs_v

        endif
        do n=1,nsv
          subs_sv =  0.5*whls_dzh(2)  *(sv0(i,j,2,n)-sv0(i,j,1,n)  )
          svp(i,j,1,n) = svp(i,j,1,n)-subs_sv
        enddo
     else
        qtp(i,j,1)  = qtp (i,j,1) + dqtdtls(1)
     endif
    end do
  end do

!     1.2 other model levels twostream

  do k=2,kmax
     kp=k+1
     km=k-1
     if (whls_dzh(kp).lt.0) then   !downwind scheme for subsidence
        do j=2,j1
           do i=2,i1
              subs_thl    = whls_dzh(kp) * (thl0(i,j,kp) - thl0(i,j,k))
              subs_qt     = whls_dzh(kp) * (qt0 (i,j,kp) - qt0 (i,j,k))
              thlp(i,j,k) = thlp(i,j,k) - subs_thl
              qtp (i,j,k) = qtp (i,j,k) - subs_qt + dqtdtls(k)
              if (lmomsubs) then
                 subs_u    = whls_dzh(kp) * (u0(i,j,kp) - u0(i,j,k))
                 subs_v    = whls_dzh(kp) * (v0(i,j,kp) - v0(i,j,k))
                 up  (i,j,k) = up  (i,j,k) - subs_u
                 vp  (i,j,k) = vp  (i,j,k) - subs_v
              endif
              do n=1,nsv
                 subs_sv   = whls_dzh(kp)  *(sv0(i,j,kp,n) - sv0(i,j,k,n))
                 svp(i,j,k,n) = svp(i,j,k,n)-subs_sv
              enddo
           enddo
        enddo
     else !downwind scheme for mean upward motions
        do j=2,j1
           do i=2,i1       
              subs_thl    = whls_dzh(k) * (thl0(i,j,k) - thl0(i,j,km))
              subs_qt     = whls_dzh(k) * (qt0 (i,j,k) - qt0 (i,j,km))
              thlp(i,j,k) = thlp(i,j,k) - subs_thl
              qtp (i,j,k) = qtp (i,j,k) - subs_qt + dqtdtls(k)
              
              if (lmomsubs) then
                 subs_u    = whls_dzh(k) * (u0(i,j,k) - u0(i,j,km))
                 subs_v    = whls_dzh(k) * (v0(i,j,k) - v0(i,j,km))
                 up  (i,j,k) = up  (i,j,k) - subs_u
                 vp  (i,j,k) = vp  (i,j,k) - subs_v
              endif
              do n=1,nsv
                 subs_sv   = whls_dzh(k) * (sv0(i,j,k,n) - sv0(i,j,km,n))
                 svp(i,j,k,n) = svp(i,j,k,n)-subs_sv
              enddo
           enddo
        enddo
     endif
 enddo

  return
end subroutine lstend


! original version of lstend
subroutine lstend_old

!-----------------------------------------------------------------|
!                                                                 |
!*** *lstend*  calculates large-scale tendencies                  |
!                                                                 |
!      Pier Siebesma   K.N.M.I.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     calculates and adds large-scale tendencies due to           |
!     large scale advection and subsidence.                       |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *lstend* is called from *program*.                  |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,nsv,lmomsubs
  use modfields, only : up,vp,thlp,qtp,svp,&
                        whls, u0av,v0av,thl0,qt0,sv0,u0,v0,&
                        dudxls,dudyls,dvdxls,dvdyls,dthldxls,dthldyls,dqtdxls,dqtdyls,dqtdtls
  implicit none

  integer i,j,k,n,kp,km
  real subs_thl,subs_qt,subs_sv,subs_u,subs_v

!     1. DETERMINE LARGE SCALE TENDENCIES
!        --------------------------------

!     1.1 lowest model level above surface : only downward component

  subs_thl = 0.
  subs_qt  = 0.
  subs_sv  = 0.
  subs_u   = 0.
  subs_v   = 0.

  do j=2,j1
    do i=2,i1
      k = 1
      if (whls(2).lt.0) then !neglect effect of mean ascending on tendencies at the lowest full level
        subs_thl     = 0.5*whls(2)  *(thl0(i,j,2)-thl0(i,j,1))/dzh(2)
        subs_qt      = 0.5*whls(2)  *(qt0(i,j,2)-qt0(i,j,1) )/dzh(2)
        if (lmomsubs) then
          subs_u     = 0.5*whls(2)  *(u0(i,j,2)-u0(i,j,1))/dzh(2)
          subs_v     = 0.5*whls(2)  *(v0(i,j,2)-v0(i,j,1))/dzh(2)
        endif
        do n=1,nsv
          subs_sv =  0.5*whls(2)  *(sv0(i,j,2,n)-sv0(i,j,1,n)  )/dzh(2)
          svp(i,j,1,n) = svp(i,j,1,n)-subs_sv
        enddo
      endif
      thlp(i,j,1) = thlp(i,j,1) -u0av(1)*dthldxls(1)-v0av(1)*dthldyls(1)-subs_thl
      qtp(i,j,1)  = qtp (i,j,1) -u0av(1)*dqtdxls (1)-v0av(1)*dqtdyls (1)-subs_qt +dqtdtls(1)
      up  (i,j,1) = up  (i,j,1) -u0av(1)*dudxls  (1)-v0av(1)*dudyls  (1)-subs_u
      vp  (i,j,1) = vp  (i,j,1) -u0av(1)*dvdxls  (1)-v0av(1)*dvdyls  (1)-subs_v
    end do
  end do

!     1.2 other model levels twostream

  do k=2,kmax
    kp=k+1
    km=k-1
    do j=2,j1
      do i=2,i1
        if (whls(kp).lt.0) then   !downwind scheme for subsidence
          subs_thl    = whls(kp) * (thl0(i,j,kp) - thl0(i,j,k))/dzh(kp)
          subs_qt     = whls(kp) * (qt0 (i,j,kp) - qt0 (i,j,k))/dzh(kp)
          if (lmomsubs) then
            subs_u    = whls(kp) * (u0(i,j,kp) - u0(i,j,k))/dzh(kp)
            subs_v    = whls(kp) * (v0(i,j,kp) - v0(i,j,k))/dzh(kp)
          endif
          do n=1,nsv
            subs_sv   = whls(kp)  *(sv0(i,j,kp,n) - sv0(i,j,k,n))/dzh(kp)
            svp(i,j,k,n) = svp(i,j,k,n)-subs_sv
          enddo
        else !downwind scheme for mean upward motions
          subs_thl    = whls(k) * (thl0(i,j,k) - thl0(i,j,km))/dzh(k)
          subs_qt     = whls(k) * (qt0 (i,j,k) - qt0 (i,j,km))/dzh(k)
          if (lmomsubs) then
            subs_u    = whls(k) * (u0(i,j,k) - u0(i,j,km))/dzh(k)
            subs_v    = whls(k) * (v0(i,j,k) - v0(i,j,km))/dzh(k)
          endif
          do n=1,nsv
            subs_sv   = whls(k) * (sv0(i,j,k,n) - sv0(i,j,km,n))/dzh(k)
            svp(i,j,k,n) = svp(i,j,k,n)-subs_sv
          enddo
        endif

        thlp(i,j,k) = thlp(i,j,k)-u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k)-subs_thl
        qtp (i,j,k) = qtp (i,j,k)-u0av(k)*dqtdxls (k)-v0av(k)*dqtdyls (k)-subs_qt+dqtdtls(k)
        up  (i,j,k) = up  (i,j,k)-u0av(k)*dudxls  (k)-v0av(k)*dudyls  (k)-subs_u
        vp  (i,j,k) = vp  (i,j,k)-u0av(k)*dvdxls  (k)-v0av(k)*dvdyls  (k)-subs_v
      enddo
    enddo
  enddo

  return
  end subroutine lstend_old

  
end module modforces
