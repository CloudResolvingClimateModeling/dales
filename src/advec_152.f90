!> \file advec_152.f90
!!  Does advection with a 5th order upwind scheme and 2nd order in the vertical
!! \par Revision list
!! \par Authors
!! \see Wicker and Scamarock 2002
!!
!! This is an experimental version of advec_52, with all advections merged.
!! 
!! By adding a small dissipative term to the sixth order flux, a fifth order
!! scheme is created that is nearly monotone:
!! \latexonly
!! \begin{eqnarray}
!!  F_{i-\frac{1}{2}}^{5th} &=& F_{i-\frac{1}{2}}^{6th} -
!! \left|\frac{\fav{u}_{i-\frac{1}{2}}}{60}\right|\left[10(\phi_i-\phi_{i-1})\right
!! . \nonumber\\\\
!! &&-\left.5(\phi_{i+1}-\phi_{i-2})+(\phi_{i+2}-\phi_{i-3})\right].
!! \end{eqnarray}
!! \endlatexonly
!!
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


!> Advection of all quantities
!> copy-paste of advec_52, merged into one double loop for the surface, and one triple loop for the bulk
subroutine advec_152()
  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf,dzh,dxi,dyi
  use modfields, only : u0, up, v0, vp, w0, wp, rhobf, rhobh, thl0, thlp, qt0, qtp, sv0, svp, e120, e12p
  implicit none

  integer :: i,j,k
  real ::  inv2dzfk, rhobf_p, rhobf_m
  k = 1
  inv2dzfk = 1./(2. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  
  do j=2,j1
     do i=2,i1
                up(i,j,k)  = up(i,j,k)- ( &
              (&
                  (u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(37.*(u0(i+1,j,k)+u0(i,j,k))-8.*(u0(i+2,j,k)+u0(i-1,j,k))+(u0(i+3,j,k)+u0(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(10.*(u0(i+1,j,k)-u0(i,j,k))-5.*(u0(i+2,j,k)-u0(i-1,j,k))+(u0(i+3,j,k)-u0(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(37.*(u0(i,j,k)+u0(i-1,j,k))-8.*(u0(i+1,j,k)+u0(i-2,j,k))+(u0(i+2,j,k)+u0(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(10.*(u0(i,j,k)-u0(i-1,j,k))-5.*(u0(i+1,j,k)-u0(i-2,j,k))+(u0(i+2,j,k)-u0(i-3,j,k)))&
              )*dxi5 &
            +(&
                  (v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(37.*(u0(i,j+1,k)+u0(i,j,k))-8.*(u0(i,j+2,k)+u0(i,j-1,k))+(u0(i,j+3,k)+u0(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i-1,j+1,k)))*(v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(10.*(u0(i,j+1,k)-u0(i,j,k))-5.*(u0(i,j+2,k)-u0(i,j-1,k))+(u0(i,j+3,k)-u0(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(37.*(u0(i,j,k)+u0(i,j-1,k))-8.*(u0(i,j+1,k)+u0(i,j-2,k))+(u0(i,j+2,k)+u0(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i-1,j,k)))*(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(10.*(u0(i,j,k)-u0(i,j-1,k))-5.*(u0(i,j+1,k)-u0(i,j-2,k))+(u0(i,j+2,k)-u0(i,j-3,k)))&
              )* dyi5 &
            +(1./rhobf(k))*( &
              ( rhobf(k+1)*u0(i,j,k+1) + rhobf(k) * u0(i,j,k)) *(w0(i,j,k+1)+ w0(i-1,j,k+1)) &
              ) / (4.*dzf(k)) &
              )


                vp(i,j,k)  = vp(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(37.*(v0(i+1,j,k)+v0(i,j,k))-8.*(v0(i+2,j,k)+v0(i-1,j,k))+(v0(i+3,j,k)+v0(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i+1,j-1,k)))*(u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(10.*(v0(i+1,j,k)-v0(i,j,k))-5.*(v0(i+2,j,k)-v0(i-1,j,k))+(v0(i+3,j,k)-v0(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(37.*(v0(i,j,k)+v0(i-1,j,k))-8.*(v0(i+1,j,k)+v0(i-2,j,k))+(v0(i+2,j,k)+v0(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i,j-1,k)))*(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(10.*(v0(i,j,k)-v0(i-1,j,k))-5.*(v0(i+1,j,k)-v0(i-2,j,k))+(v0(i+2,j,k)-v0(i-3,j,k)))&
                )*dxi5&
              +(&
                  (v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(37.*(v0(i,j+1,k)+v0(i,j,k))-8.*(v0(i,j+2,k)+v0(i,j-1,k))+(v0(i,j+3,k)+v0(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i,j,k)))*(v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(10.*(v0(i,j+1,k)-v0(i,j,k))-5.*(v0(i,j+2,k)-v0(i,j-1,k))+(v0(i,j+3,k)-v0(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(37.*(v0(i,j,k)+v0(i,j-1,k))-8.*(v0(i,j+1,k)+v0(i,j-2,k))+(v0(i,j+2,k)+v0(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i,j-1,k)))*(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(10.*(v0(i,j,k)-v0(i,j-1,k))-5.*(v0(i,j+1,k)-v0(i,j-2,k))+(v0(i,j+2,k)-v0(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                (w0(i,j,k+1)+w0(i,j-1,k+1)) *(rhobf(k+1) * v0(i,j,k+1) + rhobf(k) * v0(i,j,k)) &
                ) / (4. * dzf(k)) &
                )

                ! cell center - tke
              e12p(i,j,k)  = e12p(i,j,k)- ( &
                ( &
                  u0(i+1,j,k)/60.&
                  *(37.*(e120(i+1,j,k)+e120(i,j,k))-8.*(e120(i+2,j,k)+e120(i-1,j,k))+(e120(i+3,j,k)+e120(i-2,j,k)))&
                  -abs(u0(i+1,j,k))/60.&
                  *(10.*(e120(i+1,j,k)-e120(i,j,k))-5.*(e120(i+2,j,k)-e120(i-1,j,k))+(e120(i+3,j,k)-e120(i-2,j,k)))&
                  -u0(i,j,k)/60.&
                  *(37.*(e120(i,j,k)+e120(i-1,j,k))-8.*(e120(i+1,j,k)+e120(i-2,j,k))+(e120(i+2,j,k)+e120(i-3,j,k)))&
                  +abs(u0(i,j,k))/60.&
                  *(10.*(e120(i,j,k)-e120(i-1,j,k))-5.*(e120(i+1,j,k)-e120(i-2,j,k))+(e120(i+2,j,k)-e120(i-3,j,k)))&
                  )*dxi&
                +(&
                  v0(i,j+1,k)/60.&
                  *(37.*(e120(i,j+1,k)+e120(i,j,k))-8.*(e120(i,j+2,k)+e120(i,j-1,k))+(e120(i,j+3,k)+e120(i,j-2,k)))&
                  -abs(v0(i,j+1,k))/60.&
                  *(10.*(e120(i,j+1,k)-e120(i,j,k))-5.*(e120(i,j+2,k)-e120(i,j-1,k))+(e120(i,j+3,k)-e120(i,j-2,k)))&
                  -v0(i,j,k)/60.&
                  *(37.*(e120(i,j,k)+e120(i,j-1,k))-8.*(e120(i,j+1,k)+e120(i,j-2,k))+(e120(i,j+2,k)+e120(i,j-3,k)))&
                  +abs(v0(i,j,k))/60.&
                  *(10.*(e120(i,j,k)-e120(i,j-1,k))-5.*(e120(i,j+1,k)-e120(i,j-2,k))+(e120(i,j+2,k)-e120(i,j-3,k))) &
                  )* dyi &
                + ( &
                w0(i,j,k+1) * (rhobf_p * e120(i,j,k+1) + e120(i,j,k)) &
                ) * inv2dzfk  &
                )

                ! cell center - thl
                thlp(i,j,k)  = thlp(i,j,k)- ( &
                ( &
                  u0(i+1,j,k)/60.&
                  *(37.*(thl0(i+1,j,k)+thl0(i,j,k))-8.*(thl0(i+2,j,k)+thl0(i-1,j,k))+(thl0(i+3,j,k)+thl0(i-2,j,k)))&
                  -abs(u0(i+1,j,k))/60.&
                  *(10.*(thl0(i+1,j,k)-thl0(i,j,k))-5.*(thl0(i+2,j,k)-thl0(i-1,j,k))+(thl0(i+3,j,k)-thl0(i-2,j,k)))&
                  -u0(i,j,k)/60.&
                  *(37.*(thl0(i,j,k)+thl0(i-1,j,k))-8.*(thl0(i+1,j,k)+thl0(i-2,j,k))+(thl0(i+2,j,k)+thl0(i-3,j,k)))&
                  +abs(u0(i,j,k))/60.&
                  *(10.*(thl0(i,j,k)-thl0(i-1,j,k))-5.*(thl0(i+1,j,k)-thl0(i-2,j,k))+(thl0(i+2,j,k)-thl0(i-3,j,k)))&
                  )*dxi&
                +(&
                  v0(i,j+1,k)/60.&
                  *(37.*(thl0(i,j+1,k)+thl0(i,j,k))-8.*(thl0(i,j+2,k)+thl0(i,j-1,k))+(thl0(i,j+3,k)+thl0(i,j-2,k)))&
                  -abs(v0(i,j+1,k))/60.&
                  *(10.*(thl0(i,j+1,k)-thl0(i,j,k))-5.*(thl0(i,j+2,k)-thl0(i,j-1,k))+(thl0(i,j+3,k)-thl0(i,j-2,k)))&
                  -v0(i,j,k)/60.&
                  *(37.*(thl0(i,j,k)+thl0(i,j-1,k))-8.*(thl0(i,j+1,k)+thl0(i,j-2,k))+(thl0(i,j+2,k)+thl0(i,j-3,k)))&
                  +abs(v0(i,j,k))/60.&
                  *(10.*(thl0(i,j,k)-thl0(i,j-1,k))-5.*(thl0(i,j+1,k)-thl0(i,j-2,k))+(thl0(i,j+2,k)-thl0(i,j-3,k))) &
                  )* dyi &
                + ( &
                w0(i,j,k+1) * (rhobf_p * thl0(i,j,k+1) + thl0(i,j,k)) &
                ) * inv2dzfk  &
                )

                ! cell center - qt
                qtp(i,j,k)  = qtp(i,j,k)- ( &
                ( &
                  u0(i+1,j,k)/60.&
                  *(37.*(qt0(i+1,j,k)+qt0(i,j,k))-8.*(qt0(i+2,j,k)+qt0(i-1,j,k))+(qt0(i+3,j,k)+qt0(i-2,j,k)))&
                  -abs(u0(i+1,j,k))/60.&
                  *(10.*(qt0(i+1,j,k)-qt0(i,j,k))-5.*(qt0(i+2,j,k)-qt0(i-1,j,k))+(qt0(i+3,j,k)-qt0(i-2,j,k)))&
                  -u0(i,j,k)/60.&
                  *(37.*(qt0(i,j,k)+qt0(i-1,j,k))-8.*(qt0(i+1,j,k)+qt0(i-2,j,k))+(qt0(i+2,j,k)+qt0(i-3,j,k)))&
                  +abs(u0(i,j,k))/60.&
                  *(10.*(qt0(i,j,k)-qt0(i-1,j,k))-5.*(qt0(i+1,j,k)-qt0(i-2,j,k))+(qt0(i+2,j,k)-qt0(i-3,j,k)))&
                  )*dxi&
                +(&
                  v0(i,j+1,k)/60.&
                  *(37.*(qt0(i,j+1,k)+qt0(i,j,k))-8.*(qt0(i,j+2,k)+qt0(i,j-1,k))+(qt0(i,j+3,k)+qt0(i,j-2,k)))&
                  -abs(v0(i,j+1,k))/60.&
                  *(10.*(qt0(i,j+1,k)-qt0(i,j,k))-5.*(qt0(i,j+2,k)-qt0(i,j-1,k))+(qt0(i,j+3,k)-qt0(i,j-2,k)))&
                  -v0(i,j,k)/60.&
                  *(37.*(qt0(i,j,k)+qt0(i,j-1,k))-8.*(qt0(i,j+1,k)+qt0(i,j-2,k))+(qt0(i,j+2,k)+qt0(i,j-3,k)))&
                  +abs(v0(i,j,k))/60.&
                  *(10.*(qt0(i,j,k)-qt0(i,j-1,k))-5.*(qt0(i,j+1,k)-qt0(i,j-2,k))+(qt0(i,j+2,k)-qt0(i,j-3,k))) &
                  )* dyi &
                + ( &
                w0(i,j,k+1) * (rhobf_p * qt0(i,j,k+1) + qt0(i,j,k)) &
                ) * inv2dzfk  &
                )

                ! cell center - sv2
                svp(i,j,k,2)  = svp(i,j,k,2)- ( &
                ( &
                  u0(i+1,j,k)/60.&
                  *(37.*(sv0(i+1,j,k,2)+sv0(i,j,k,2))-8.*(sv0(i+2,j,k,2)+sv0(i-1,j,k,2))+(sv0(i+3,j,k,2)+sv0(i-2,j,k,2)))&
                  -abs(u0(i+1,j,k))/60.&
                  *(10.*(sv0(i+1,j,k,2)-sv0(i,j,k,2))-5.*(sv0(i+2,j,k,2)-sv0(i-1,j,k,2))+(sv0(i+3,j,k,2)-sv0(i-2,j,k,2)))&
                  -u0(i,j,k)/60.&
                  *(37.*(sv0(i,j,k,2)+sv0(i-1,j,k,2))-8.*(sv0(i+1,j,k,2)+sv0(i-2,j,k,2))+(sv0(i+2,j,k,2)+sv0(i-3,j,k,2)))&
                  +abs(u0(i,j,k))/60.&
                  *(10.*(sv0(i,j,k,2)-sv0(i-1,j,k,2))-5.*(sv0(i+1,j,k,2)-sv0(i-2,j,k,2))+(sv0(i+2,j,k,2)-sv0(i-3,j,k,2)))&
                  )*dxi&
                +(&
                  v0(i,j+1,k)/60.&
                  *(37.*(sv0(i,j+1,k,2)+sv0(i,j,k,2))-8.*(sv0(i,j+2,k,2)+sv0(i,j-1,k,2))+(sv0(i,j+3,k,2)+sv0(i,j-2,k,2)))&
                  -abs(v0(i,j+1,k))/60.&
                  *(10.*(sv0(i,j+1,k,2)-sv0(i,j,k,2))-5.*(sv0(i,j+2,k,2)-sv0(i,j-1,k,2))+(sv0(i,j+3,k,2)-sv0(i,j-2,k,2)))&
                  -v0(i,j,k)/60.&
                  *(37.*(sv0(i,j,k,2)+sv0(i,j-1,k,2))-8.*(sv0(i,j+1,k,2)+sv0(i,j-2,k,2))+(sv0(i,j+2,k,2)+sv0(i,j-3,k,2)))&
                  +abs(v0(i,j,k))/60.&
                  *(10.*(sv0(i,j,k,2)-sv0(i,j-1,k,2))-5.*(sv0(i,j+1,k,2)-sv0(i,j-2,k,2))+(sv0(i,j+2,k,2)-sv0(i,j-3,k,2))) &
                  )* dyi &
                + ( &
                w0(i,j,k+1) * (rhobf_p * sv0(i,j,k+1,2) + sv0(i,j,k,2)) &
                ) * inv2dzfk  &
                )
                
                
     end do
  end do
  
  
  do k=2,kmax
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     rhobf_m = rhobf(k-1)/rhobf(k)

    do j=2,j1
      do i=2,i1
 
        up(i,j,k)  = up(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(37.*(u0(i+1,j,k)+u0(i,j,k))-8.*(u0(i+2,j,k)+u0(i-1,j,k))+(u0(i+3,j,k)+u0(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(10.*(u0(i+1,j,k)-u0(i,j,k))-5.*(u0(i+2,j,k)-u0(i-1,j,k))+(u0(i+3,j,k)-u0(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(37.*(u0(i,j,k)+u0(i-1,j,k))-8.*(u0(i+1,j,k)+u0(i-2,j,k))+(u0(i+2,j,k)+u0(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(10.*(u0(i,j,k)-u0(i-1,j,k))-5.*(u0(i+1,j,k)-u0(i-2,j,k))+(u0(i+2,j,k)-u0(i-3,j,k)))&
              )*dxi5&
            +(&
                  (v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(37.*(u0(i,j+1,k)+u0(i,j,k))-8.*(u0(i,j+2,k)+u0(i,j-1,k))+(u0(i,j+3,k)+u0(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i-1,j+1,k)))*(v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(10.*(u0(i,j+1,k)-u0(i,j,k))-5.*(u0(i,j+2,k)-u0(i,j-1,k))+(u0(i,j+3,k)-u0(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(37.*(u0(i,j,k)+u0(i,j-1,k))-8.*(u0(i,j+1,k)+u0(i,j-2,k))+(u0(i,j+2,k)+u0(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i-1,j,k)))*(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(10.*(u0(i,j,k)-u0(i,j-1,k))-5.*(u0(i,j+1,k)-u0(i,j-2,k))+(u0(i,j+2,k)-u0(i,j-3,k)))&
              )* dyi5 &
            +(1./rhobf(k))*( &
              (rhobf(k) * u0(i,j,k) + rhobf(k+1) * u0(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
            -(rhobf(k) * u0(i,j,k) + rhobf(k-1) * u0(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
              ) / (4. * dzf(k)) &
              )


        vp(i,j,k)  = vp(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(37.*(v0(i+1,j,k)+v0(i,j,k))-8.*(v0(i+2,j,k)+v0(i-1,j,k))+(v0(i+3,j,k)+v0(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i+1,j-1,k)))*(u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(10.*(v0(i+1,j,k)-v0(i,j,k))-5.*(v0(i+2,j,k)-v0(i-1,j,k))+(v0(i+3,j,k)-v0(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(37.*(v0(i,j,k)+v0(i-1,j,k))-8.*(v0(i+1,j,k)+v0(i-2,j,k))+(v0(i+2,j,k)+v0(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i,j-1,k)))*(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(10.*(v0(i,j,k)-v0(i-1,j,k))-5.*(v0(i+1,j,k)-v0(i-2,j,k))+(v0(i+2,j,k)-v0(i-3,j,k)))&
                )*dxi5&
              +(&
                  (v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(37.*(v0(i,j+1,k)+v0(i,j,k))-8.*(v0(i,j+2,k)+v0(i,j-1,k))+(v0(i,j+3,k)+v0(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i,j,k)))*(v0(i,j+1,k)+v0(i,j+1,k))/60.&
                  *(10.*(v0(i,j+1,k)-v0(i,j,k))-5.*(v0(i,j+2,k)-v0(i,j-1,k))+(v0(i,j+3,k)-v0(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(37.*(v0(i,j,k)+v0(i,j-1,k))-8.*(v0(i,j+1,k)+v0(i,j-2,k))+(v0(i,j+2,k)+v0(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i,j-1,k)))*(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(10.*(v0(i,j,k)-v0(i,j-1,k))-5.*(v0(i,j+1,k)-v0(i,j-2,k))+(v0(i,j+2,k)-v0(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1) * v0(i,j,k+1) + rhobf(k) * v0(i,j,k)) &
                -(w0(i,j,k) +w0(i,j-1,k))  *(rhobf(k-1) * v0(i,j,k-1) + rhobf(k) * v0(i,j,k)) &
                ) / (4. * dzf(k)) &
                )


        wp(i,j,k)  = wp(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i+1,j,k-1))/60.&
                    *(37.*(w0(i+1,j,k)+w0(i,j,k))-8.*(w0(i+2,j,k)+w0(i-1,j,k))+(w0(i+3,j,k)+w0(i-2,j,k)))&
                    -sign(1.,(u0(i+1,j,k)+u0(i+1,j,k-1)))*(u0(i+1,j,k)+u0(i+1,j,k-1))/60.&
                    *(10.*(w0(i+1,j,k)-w0(i,j,k))-5.*(w0(i+2,j,k)-w0(i-1,j,k))+(w0(i+3,j,k)-w0(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j,k-1))/60.&
                    *(37.*(w0(i,j,k)+w0(i-1,j,k))-8.*(w0(i+1,j,k)+w0(i-2,j,k))+(w0(i+2,j,k)+w0(i-3,j,k)))&
                    +sign(1.,(u0(i,j,k)+u0(i,j,k-1)))*(u0(i,j,k)+u0(i,j,k-1))/60.&
                    *(10.*(w0(i,j,k)-w0(i-1,j,k))-5.*(w0(i+1,j,k)-w0(i-2,j,k))+(w0(i+2,j,k)-w0(i-3,j,k)))&
                )*dxi5&
              + (&
                    (v0(i,j+1,k)+v0(i,j+1,k-1))/60.&
                    *(37.*(w0(i,j+1,k)+w0(i,j,k))-8.*(w0(i,j+2,k)+w0(i,j-1,k))+(w0(i,j+3,k)+w0(i,j-2,k)))&
                    -sign(1.,(v0(i,j+1,k)+v0(i,j+1,k-1)))*(v0(i,j+1,k)+v0(i,j+1,k-1))/60.&
                    *(10.*(w0(i,j+1,k)-w0(i,j,k))-5.*(w0(i,j+2,k)-w0(i,j-1,k))+(w0(i,j+3,k)-w0(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j,k-1))/60.&
                    *(37.*(w0(i,j,k)+w0(i,j-1,k))-8.*(w0(i,j+1,k)+w0(i,j-2,k))+(w0(i,j+2,k)+w0(i,j-3,k)))&
                    +sign(1.,(v0(i,j,k)+v0(i,j,k-1)))*(v0(i,j,k)+v0(i,j,k-1))/60.&
                    *(10.*(w0(i,j,k)-w0(i,j-1,k))-5.*(w0(i,j+1,k)-w0(i,j-2,k))+(w0(i,j+2,k)-w0(i,j-3,k)))&
                )* dyi5 &
              + (1./rhobh(k))*( &
              (rhobh(k) * w0(i,j,k) + rhobh(k+1) * w0(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                -(rhobh(k) * w0(i,j,k) + rhobh(k-1) * w0(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                )/ (4. * dzh(k)) &
                )



        ! cell center - tke
        e12p(i,j,k)  = e12p(i,j,k)- (  &
             ( &
                      u0(i+1,j,k)/60.&
                      *(37.*(e120(i+1,j,k)+e120(i,j,k))-8.*(e120(i+2,j,k)+e120(i-1,j,k))+(e120(i+3,j,k)+e120(i-2,j,k)))&
                      -abs(u0(i+1,j,k))/60.&
                      *(10.*(e120(i+1,j,k)-e120(i,j,k))-5.*(e120(i+2,j,k)-e120(i-1,j,k))+(e120(i+3,j,k)-e120(i-2,j,k)))&
                      -u0(i,j,k)/60.&
                      *(37.*(e120(i,j,k)+e120(i-1,j,k))-8.*(e120(i+1,j,k)+e120(i-2,j,k))+(e120(i+2,j,k)+e120(i-3,j,k)))&
                      +abs(u0(i,j,k))/60.&
                      *(10.*(e120(i,j,k)-e120(i-1,j,k))-5.*(e120(i+1,j,k)-e120(i-2,j,k))+(e120(i+2,j,k)-e120(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60.&
                      *(37.*(e120(i,j+1,k)+e120(i,j,k))-8.*(e120(i,j+2,k)+e120(i,j-1,k))+(e120(i,j+3,k)+e120(i,j-2,k)))&
                      -abs(v0(i,j+1,k))/60.&
                      *(10.*(e120(i,j+1,k)-e120(i,j,k))-5.*(e120(i,j+2,k)-e120(i,j-1,k))+(e120(i,j+3,k)-e120(i,j-2,k)))&
                      -v0(i,j,k)/60.&
                      *(37.*(e120(i,j,k)+e120(i,j-1,k))-8.*(e120(i,j+1,k)+e120(i,j-2,k))+(e120(i,j+2,k)+e120(i,j-3,k)))&
                      +abs(v0(i,j,k))/60.&
                      *(10.*(e120(i,j,k)-e120(i,j-1,k))-5.*(e120(i,j+1,k)-e120(i,j-2,k))+(e120(i,j+2,k)-e120(i,j-3,k)))&
                  )* dyi &
                + ( &
                  w0(i,j,k+1) * (rhobf_p * e120(i,j,k+1) +  e120(i,j,k)) &
                  -w0(i,j,k)  * (rhobf_m * e120(i,j,k-1) +  e120(i,j,k)) &
                  ) * inv2dzfk &
                  )

        ! cell center - thl
        thlp(i,j,k)  = thlp(i,j,k)- (  &
                  ( &
                  u0(i+1,j,k)/60.&
                      *(37.*(thl0(i+1,j,k)+thl0(i,j,k))-8.*(thl0(i+2,j,k)+thl0(i-1,j,k))+(thl0(i+3,j,k)+thl0(i-2,j,k)))&
                      -abs(u0(i+1,j,k))/60.&
                      *(10.*(thl0(i+1,j,k)-thl0(i,j,k))-5.*(thl0(i+2,j,k)-thl0(i-1,j,k))+(thl0(i+3,j,k)-thl0(i-2,j,k)))&
                      -u0(i,j,k)/60.&
                      *(37.*(thl0(i,j,k)+thl0(i-1,j,k))-8.*(thl0(i+1,j,k)+thl0(i-2,j,k))+(thl0(i+2,j,k)+thl0(i-3,j,k)))&
                      +abs(u0(i,j,k))/60.&
                      *(10.*(thl0(i,j,k)-thl0(i-1,j,k))-5.*(thl0(i+1,j,k)-thl0(i-2,j,k))+(thl0(i+2,j,k)-thl0(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60.&
                      *(37.*(thl0(i,j+1,k)+thl0(i,j,k))-8.*(thl0(i,j+2,k)+thl0(i,j-1,k))+(thl0(i,j+3,k)+thl0(i,j-2,k)))&
                      -abs(v0(i,j+1,k))/60.&
                      *(10.*(thl0(i,j+1,k)-thl0(i,j,k))-5.*(thl0(i,j+2,k)-thl0(i,j-1,k))+(thl0(i,j+3,k)-thl0(i,j-2,k)))&
                      -v0(i,j,k)/60.&
                      *(37.*(thl0(i,j,k)+thl0(i,j-1,k))-8.*(thl0(i,j+1,k)+thl0(i,j-2,k))+(thl0(i,j+2,k)+thl0(i,j-3,k)))&
                      +abs(v0(i,j,k))/60.&
                      *(10.*(thl0(i,j,k)-thl0(i,j-1,k))-5.*(thl0(i,j+1,k)-thl0(i,j-2,k))+(thl0(i,j+2,k)-thl0(i,j-3,k)))&
                  )* dyi &
                + ( &
                  w0(i,j,k+1) * (rhobf_p * thl0(i,j,k+1) +  thl0(i,j,k)) &
                  -w0(i,j,k)  * (rhobf_m * thl0(i,j,k-1) +  thl0(i,j,k)) &
                  ) * inv2dzfk &
                  )

        ! cell center - qt
        qtp(i,j,k)  = qtp(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60.&
                      *(37.*(qt0(i+1,j,k)+qt0(i,j,k))-8.*(qt0(i+2,j,k)+qt0(i-1,j,k))+(qt0(i+3,j,k)+qt0(i-2,j,k)))&
                      -abs(u0(i+1,j,k))/60.&
                      *(10.*(qt0(i+1,j,k)-qt0(i,j,k))-5.*(qt0(i+2,j,k)-qt0(i-1,j,k))+(qt0(i+3,j,k)-qt0(i-2,j,k)))&
                      -u0(i,j,k)/60.&
                      *(37.*(qt0(i,j,k)+qt0(i-1,j,k))-8.*(qt0(i+1,j,k)+qt0(i-2,j,k))+(qt0(i+2,j,k)+qt0(i-3,j,k)))&
                      +abs(u0(i,j,k))/60.&
                      *(10.*(qt0(i,j,k)-qt0(i-1,j,k))-5.*(qt0(i+1,j,k)-qt0(i-2,j,k))+(qt0(i+2,j,k)-qt0(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60.&
                      *(37.*(qt0(i,j+1,k)+qt0(i,j,k))-8.*(qt0(i,j+2,k)+qt0(i,j-1,k))+(qt0(i,j+3,k)+qt0(i,j-2,k)))&
                      -abs(v0(i,j+1,k))/60.&
                      *(10.*(qt0(i,j+1,k)-qt0(i,j,k))-5.*(qt0(i,j+2,k)-qt0(i,j-1,k))+(qt0(i,j+3,k)-qt0(i,j-2,k)))&
                      -v0(i,j,k)/60.&
                      *(37.*(qt0(i,j,k)+qt0(i,j-1,k))-8.*(qt0(i,j+1,k)+qt0(i,j-2,k))+(qt0(i,j+2,k)+qt0(i,j-3,k)))&
                      +abs(v0(i,j,k))/60.&
                      *(10.*(qt0(i,j,k)-qt0(i,j-1,k))-5.*(qt0(i,j+1,k)-qt0(i,j-2,k))+(qt0(i,j+2,k)-qt0(i,j-3,k)))&
                  )* dyi &
                + ( &
                  w0(i,j,k+1) * (rhobf_p * qt0(i,j,k+1) +  qt0(i,j,k)) &
                  -w0(i,j,k)  * (rhobf_m * qt0(i,j,k-1) +  qt0(i,j,k)) &
                  ) * inv2dzfk &
                  )


        ! cell center - sv2          
        svp(i,j,k,2)  = svp(i,j,k,2)- (  &
                  ( &
                      u0(i+1,j,k)/60.&
                      *(37.*(sv0(i+1,j,k,2)+sv0(i,j,k,2))-8.*(sv0(i+2,j,k,2)+sv0(i-1,j,k,2))+(sv0(i+3,j,k,2)+sv0(i-2,j,k,2)))&
                      -abs(u0(i+1,j,k))/60.&
                      *(10.*(sv0(i+1,j,k,2)-sv0(i,j,k,2))-5.*(sv0(i+2,j,k,2)-sv0(i-1,j,k,2))+(sv0(i+3,j,k,2)-sv0(i-2,j,k,2)))&
                      -u0(i,j,k)/60.&
                      *(37.*(sv0(i,j,k,2)+sv0(i-1,j,k,2))-8.*(sv0(i+1,j,k,2)+sv0(i-2,j,k,2))+(sv0(i+2,j,k,2)+sv0(i-3,j,k,2)))&
                      +abs(u0(i,j,k))/60.&
                      *(10.*(sv0(i,j,k,2)-sv0(i-1,j,k,2))-5.*(sv0(i+1,j,k,2)-sv0(i-2,j,k,2))+(sv0(i+2,j,k,2)-sv0(i-3,j,k,2)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60.&
                      *(37.*(sv0(i,j+1,k,2)+sv0(i,j,k,2))-8.*(sv0(i,j+2,k,2)+sv0(i,j-1,k,2))+(sv0(i,j+3,k,2)+sv0(i,j-2,k,2)))&
                      -abs(v0(i,j+1,k))/60.&
                      *(10.*(sv0(i,j+1,k,2)-sv0(i,j,k,2))-5.*(sv0(i,j+2,k,2)-sv0(i,j-1,k,2))+(sv0(i,j+3,k,2)-sv0(i,j-2,k,2)))&
                      -v0(i,j,k)/60.&
                      *(37.*(sv0(i,j,k,2)+sv0(i,j-1,k,2))-8.*(sv0(i,j+1,k,2)+sv0(i,j-2,k,2))+(sv0(i,j+2,k,2)+sv0(i,j-3,k,2)))&
                      +abs(v0(i,j,k))/60.&
                      *(10.*(sv0(i,j,k,2)-sv0(i,j-1,k,2))-5.*(sv0(i,j+1,k,2)-sv0(i,j-2,k,2))+(sv0(i,j+2,k,2)-sv0(i,j-3,k,2)))&
                  )* dyi &
                + ( &
                  w0(i,j,k+1) * (rhobf_p * sv0(i,j,k+1,2) +  sv0(i,j,k,2)) &
                  -w0(i,j,k)  * (rhobf_m * sv0(i,j,k-1,2) +  sv0(i,j,k,2)) &
                  ) * inv2dzfk &
                  )

        
      end do
    end do
  end do

end subroutine advec_152



!> Advection of all velocities
!> copy-paste of advec_52, merged into one double loop for the surface, and one triple loop for the bulk
subroutine advec_252()
  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf,dzh,dxi,dyi
  use modfields, only : u0, up, v0, vp, w0, wp, rhobf, rhobh
  implicit none

  integer :: i,j,k
  real ::  inv2dzfk, rhobf_p, rhobf_m
  k = 1
  inv2dzfk = 1./(2. * dzf(k))
  rhobf_p = rhobf(k+1)/rhobf(k)
  
  do j=2,j1
     do i=2,i1
                up(i,j,k)  = up(i,j,k)- ( &
              (&
                  (u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(37.*(u0(i+1,j,k)+u0(i,j,k))-8.*(u0(i+2,j,k)+u0(i-1,j,k))+(u0(i+3,j,k)+u0(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(10.*(u0(i+1,j,k)-u0(i,j,k))-5.*(u0(i+2,j,k)-u0(i-1,j,k))+(u0(i+3,j,k)-u0(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(37.*(u0(i,j,k)+u0(i-1,j,k))-8.*(u0(i+1,j,k)+u0(i-2,j,k))+(u0(i+2,j,k)+u0(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(10.*(u0(i,j,k)-u0(i-1,j,k))-5.*(u0(i+1,j,k)-u0(i-2,j,k))+(u0(i+2,j,k)-u0(i-3,j,k)))&
              )*dxi5 &
            +(&
                  (v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(37.*(u0(i,j+1,k)+u0(i,j,k))-8.*(u0(i,j+2,k)+u0(i,j-1,k))+(u0(i,j+3,k)+u0(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i-1,j+1,k)))*(v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(10.*(u0(i,j+1,k)-u0(i,j,k))-5.*(u0(i,j+2,k)-u0(i,j-1,k))+(u0(i,j+3,k)-u0(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(37.*(u0(i,j,k)+u0(i,j-1,k))-8.*(u0(i,j+1,k)+u0(i,j-2,k))+(u0(i,j+2,k)+u0(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i-1,j,k)))*(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(10.*(u0(i,j,k)-u0(i,j-1,k))-5.*(u0(i,j+1,k)-u0(i,j-2,k))+(u0(i,j+2,k)-u0(i,j-3,k)))&
              )* dyi5 &
            +(1./rhobf(k))*( &
              ( rhobf(k+1)*u0(i,j,k+1) + rhobf(k) * u0(i,j,k)) *(w0(i,j,k+1)+ w0(i-1,j,k+1)) &
              ) / (4.*dzf(k)) &
              )


                vp(i,j,k)  = vp(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(37.*(v0(i+1,j,k)+v0(i,j,k))-8.*(v0(i+2,j,k)+v0(i-1,j,k))+(v0(i+3,j,k)+v0(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i+1,j-1,k)))*(u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(10.*(v0(i+1,j,k)-v0(i,j,k))-5.*(v0(i+2,j,k)-v0(i-1,j,k))+(v0(i+3,j,k)-v0(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(37.*(v0(i,j,k)+v0(i-1,j,k))-8.*(v0(i+1,j,k)+v0(i-2,j,k))+(v0(i+2,j,k)+v0(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i,j-1,k)))*(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(10.*(v0(i,j,k)-v0(i-1,j,k))-5.*(v0(i+1,j,k)-v0(i-2,j,k))+(v0(i+2,j,k)-v0(i-3,j,k)))&
                )*dxi5&
              +(&
                  (v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(37.*(v0(i,j+1,k)+v0(i,j,k))-8.*(v0(i,j+2,k)+v0(i,j-1,k))+(v0(i,j+3,k)+v0(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i,j,k)))*(v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(10.*(v0(i,j+1,k)-v0(i,j,k))-5.*(v0(i,j+2,k)-v0(i,j-1,k))+(v0(i,j+3,k)-v0(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(37.*(v0(i,j,k)+v0(i,j-1,k))-8.*(v0(i,j+1,k)+v0(i,j-2,k))+(v0(i,j+2,k)+v0(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i,j-1,k)))*(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(10.*(v0(i,j,k)-v0(i,j-1,k))-5.*(v0(i,j+1,k)-v0(i,j-2,k))+(v0(i,j+2,k)-v0(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                (w0(i,j,k+1)+w0(i,j-1,k+1)) *(rhobf(k+1) * v0(i,j,k+1) + rhobf(k) * v0(i,j,k)) &
                ) / (4. * dzf(k)) &
                )

     end do
  end do
  
  
  do k=2,kmax
     inv2dzfk = 1./(2. * dzf(k))
     rhobf_p = rhobf(k+1)/rhobf(k)
     rhobf_m = rhobf(k-1)/rhobf(k)

    do j=2,j1
      do i=2,i1
 
        up(i,j,k)  = up(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(37.*(u0(i+1,j,k)+u0(i,j,k))-8.*(u0(i+2,j,k)+u0(i-1,j,k))+(u0(i+3,j,k)+u0(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i,j,k)))*(u0(i+1,j,k)+u0(i,j,k))/60.&
                  *(10.*(u0(i+1,j,k)-u0(i,j,k))-5.*(u0(i+2,j,k)-u0(i-1,j,k))+(u0(i+3,j,k)-u0(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(37.*(u0(i,j,k)+u0(i-1,j,k))-8.*(u0(i+1,j,k)+u0(i-2,j,k))+(u0(i+2,j,k)+u0(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i-1,j,k)))*(u0(i,j,k)+u0(i-1,j,k))/60.&
                  *(10.*(u0(i,j,k)-u0(i-1,j,k))-5.*(u0(i+1,j,k)-u0(i-2,j,k))+(u0(i+2,j,k)-u0(i-3,j,k)))&
              )*dxi5&
            +(&
                  (v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(37.*(u0(i,j+1,k)+u0(i,j,k))-8.*(u0(i,j+2,k)+u0(i,j-1,k))+(u0(i,j+3,k)+u0(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i-1,j+1,k)))*(v0(i,j+1,k)+v0(i-1,j+1,k))/60.&
                  *(10.*(u0(i,j+1,k)-u0(i,j,k))-5.*(u0(i,j+2,k)-u0(i,j-1,k))+(u0(i,j+3,k)-u0(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(37.*(u0(i,j,k)+u0(i,j-1,k))-8.*(u0(i,j+1,k)+u0(i,j-2,k))+(u0(i,j+2,k)+u0(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i-1,j,k)))*(v0(i,j,k)+v0(i-1,j,k))/60.&
                  *(10.*(u0(i,j,k)-u0(i,j-1,k))-5.*(u0(i,j+1,k)-u0(i,j-2,k))+(u0(i,j+2,k)-u0(i,j-3,k)))&
              )* dyi5 &
            +(1./rhobf(k))*( &
              (rhobf(k) * u0(i,j,k) + rhobf(k+1) * u0(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
            -(rhobf(k) * u0(i,j,k) + rhobf(k-1) * u0(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
              ) / (4. * dzf(k)) &
              )


        vp(i,j,k)  = vp(i,j,k)- ( &
              ( &
                  (u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(37.*(v0(i+1,j,k)+v0(i,j,k))-8.*(v0(i+2,j,k)+v0(i-1,j,k))+(v0(i+3,j,k)+v0(i-2,j,k)))&
                  -sign(1.,(u0(i+1,j,k)+u0(i+1,j-1,k)))*(u0(i+1,j,k)+u0(i+1,j-1,k))/60.&
                  *(10.*(v0(i+1,j,k)-v0(i,j,k))-5.*(v0(i+2,j,k)-v0(i-1,j,k))+(v0(i+3,j,k)-v0(i-2,j,k)))&
                  -(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(37.*(v0(i,j,k)+v0(i-1,j,k))-8.*(v0(i+1,j,k)+v0(i-2,j,k))+(v0(i+2,j,k)+v0(i-3,j,k)))&
                  +sign(1.,(u0(i,j,k)+u0(i,j-1,k)))*(u0(i,j,k)+u0(i,j-1,k))/60.&
                  *(10.*(v0(i,j,k)-v0(i-1,j,k))-5.*(v0(i+1,j,k)-v0(i-2,j,k))+(v0(i+2,j,k)-v0(i-3,j,k)))&
                )*dxi5&
              +(&
                  (v0(i,j+1,k)+v0(i,j,k))/60.&
                  *(37.*(v0(i,j+1,k)+v0(i,j,k))-8.*(v0(i,j+2,k)+v0(i,j-1,k))+(v0(i,j+3,k)+v0(i,j-2,k)))&
                  -sign(1.,(v0(i,j+1,k)+v0(i,j,k)))*(v0(i,j+1,k)+v0(i,j+1,k))/60.&
                  *(10.*(v0(i,j+1,k)-v0(i,j,k))-5.*(v0(i,j+2,k)-v0(i,j-1,k))+(v0(i,j+3,k)-v0(i,j-2,k)))&
                  -(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(37.*(v0(i,j,k)+v0(i,j-1,k))-8.*(v0(i,j+1,k)+v0(i,j-2,k))+(v0(i,j+2,k)+v0(i,j-3,k)))&
                  +sign(1.,(v0(i,j,k)+v0(i,j-1,k)))*(v0(i,j,k)+v0(i,j-1,k))/60.&
                  *(10.*(v0(i,j,k)-v0(i,j-1,k))-5.*(v0(i,j+1,k)-v0(i,j-2,k))+(v0(i,j+2,k)-v0(i,j-3,k)))&
                )* dyi5 &
              +(1./rhobf(k))*( &
                (w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1) * v0(i,j,k+1) + rhobf(k) * v0(i,j,k)) &
                -(w0(i,j,k) +w0(i,j-1,k))  *(rhobf(k-1) * v0(i,j,k-1) + rhobf(k) * v0(i,j,k)) &
                ) / (4. * dzf(k)) &
                )


        wp(i,j,k)  = wp(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i+1,j,k-1))/60.&
                    *(37.*(w0(i+1,j,k)+w0(i,j,k))-8.*(w0(i+2,j,k)+w0(i-1,j,k))+(w0(i+3,j,k)+w0(i-2,j,k)))&
                    -sign(1.,(u0(i+1,j,k)+u0(i+1,j,k-1)))*(u0(i+1,j,k)+u0(i+1,j,k-1))/60.&
                    *(10.*(w0(i+1,j,k)-w0(i,j,k))-5.*(w0(i+2,j,k)-w0(i-1,j,k))+(w0(i+3,j,k)-w0(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j,k-1))/60.&
                    *(37.*(w0(i,j,k)+w0(i-1,j,k))-8.*(w0(i+1,j,k)+w0(i-2,j,k))+(w0(i+2,j,k)+w0(i-3,j,k)))&
                    +sign(1.,(u0(i,j,k)+u0(i,j,k-1)))*(u0(i,j,k)+u0(i,j,k-1))/60.&
                    *(10.*(w0(i,j,k)-w0(i-1,j,k))-5.*(w0(i+1,j,k)-w0(i-2,j,k))+(w0(i+2,j,k)-w0(i-3,j,k)))&
                )*dxi5&
              + (&
                    (v0(i,j+1,k)+v0(i,j+1,k-1))/60.&
                    *(37.*(w0(i,j+1,k)+w0(i,j,k))-8.*(w0(i,j+2,k)+w0(i,j-1,k))+(w0(i,j+3,k)+w0(i,j-2,k)))&
                    -sign(1.,(v0(i,j+1,k)+v0(i,j+1,k-1)))*(v0(i,j+1,k)+v0(i,j+1,k-1))/60.&
                    *(10.*(w0(i,j+1,k)-w0(i,j,k))-5.*(w0(i,j+2,k)-w0(i,j-1,k))+(w0(i,j+3,k)-w0(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j,k-1))/60.&
                    *(37.*(w0(i,j,k)+w0(i,j-1,k))-8.*(w0(i,j+1,k)+w0(i,j-2,k))+(w0(i,j+2,k)+w0(i,j-3,k)))&
                    +sign(1.,(v0(i,j,k)+v0(i,j,k-1)))*(v0(i,j,k)+v0(i,j,k-1))/60.&
                    *(10.*(w0(i,j,k)-w0(i,j-1,k))-5.*(w0(i,j+1,k)-w0(i,j-2,k))+(w0(i,j+2,k)-w0(i,j-3,k)))&
                )* dyi5 &
              + (1./rhobh(k))*( &
              (rhobh(k) * w0(i,j,k) + rhobh(k+1) * w0(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                -(rhobh(k) * w0(i,j,k) + rhobh(k-1) * w0(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                )/ (4. * dzh(k)) &
                )       
      end do
    end do
  end do

end subroutine advec_252




