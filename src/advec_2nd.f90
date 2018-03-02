!> \file advec_2nd.f90
!!  Does advection with a 2nd order central differencing scheme.
!! \par Revision list
!! \par Authors
!! Second order central differencing can be used for variables where neither very
!! high accuracy nor strict monotonicity is necessary.
!! \latexonly
!!\begin{eqnarray}
!! F_{i-\frac{1}{2}}^{2nd} &=&
!!\fav{u}_{i-\frac{1}{2}}\frac{\phi_{i}+\phi_{i-1}}{2},
!!\end{eqnarray}
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

!> Advection at cell center
subroutine advecc_2nd(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0, rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the cell centered field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin

  integer :: i,j,k

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rhoputin(i,j,k)=rhobf(k)*putin(i,j,k)
!      end do
!    end do
!  end do

  do k=1,kmax
    do j=2,j1
      do i=2,i1
        putout(i,j,k)  = putout(i,j,k)- (  &
              ( &
              u0(i+1,j,k) * ( putin(i+1,j,k) + putin(i,j,k) ) &
             -u0(i ,j,k) * ( putin(i-1,j,k) + putin(i,j,k) ) &
              )* dxi5 &
            +( &
              v0(i,j+1,k) * ( putin(i,j+1,k) + putin(i,j,k) ) &
             -v0(i,j ,k) * ( putin(i,j-1,k) + putin(i,j,k) ) &
              )* dyi5 )

      end do
    end do
  end do

  if (leq) then ! equidistant grid

    do j=2,j1
      do i=2,i1
        putout(i,j,1)  = putout(i,j,1)- (1./rhobf(1))*( &
                w0(i,j,2) * (rhobf(2) * putin(i,j,2) + rhobf(1) * putin(i,j,1) ) &
                ) * dzi5
      end do
    end do

    do j=2,j1
    do k=2,kmax         
       do i=2,i1
          putout(i,j,k)  = putout(i,j,k)- (1./rhobf(k))*( &
                w0(i,j,k+1) * (rhobf(k+1) * putin(i,j,k+1) + rhobf(k) * putin(i,j,k)) &
                -w0(i,j,k)   * (rhobf(k-1) * putin(i,j,k-1)+ rhobf(k) * putin(i,j,k)) &
                )*dzi5
        end do
      end do
    end do

  else   ! non-equidistant grid

    do j=2,j1
      do i=2,i1
        putout(i,j,1)  = putout(i,j,1)- (1./rhobf(1))*( &
                w0(i,j,2) * (rhobf(2) * putin(i,j,2) * dzf(1) + rhobf(1) * putin(i,j,1) * dzf(2) ) / (2.*dzh(2)) &
                ) / dzf(1)
      end do
    end do

    do j=2,j1
    do k=2,kmax
       do i=2,i1
          putout(i,j,k)  = putout(i,j,k)- (1./rhobf(k))*( &
                w0(i,j,k+1) * (rhobf(k+1) * putin(i,j,k+1) * dzf(k) + rhobf(k) * putin(i,j,k) * dzf(k+1) ) / dzh(k+1) &
               -w0(i,j,k ) * (rhobf(k-1) * putin(i,j,k-1) * dzf(k) + rhobf(k) * putin(i,j,k) * dzf(k-1) ) / dzh(k) &
                )/ (2. * dzf(k))
        end do
      end do
    end do

  end if

end subroutine advecc_2nd


!> Advection at the u point.
subroutine advecu_2nd(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0, rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: putin
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin

  integer :: i,j,k

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rhoputin(i,j,k)=rhobf(k)*putin(i,j,k)
!      end do
!    end do
!  end do

  do k=1,kmax
    do j=2,j1
      do i=2,i1
        putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                (putin(i,j,k)+putin(i+1,j,k))*(u0(i,j,k)+u0(i+1,j,k)) &
                -(putin(i,j,k)+putin(i-1,j,k))*(u0(i,j,k)+u0(i-1,j,k)) &
                )*dxiq &
                +(  &
                (putin(i,j,k)+putin(i,j+1,k))*(v0(i,j+1,k)+v0(i-1,j+1 ,k)) &
                -(putin(i,j,k)+putin(i,j-1,k))*(v0(i,j  ,k)+v0(i-1,j  ,k)) &
                )*dyiq )

      end do
    end do
  end do

  if (leq) then

    do j=2,j1
      do i=2,i1
        putout(i,j,1)  = putout(i,j,1)-(1./rhobf(1))*( &
            ( rhobf(2) * putin(i,j,2) + rhobf(1) * putin(i,j,1))*( w0(i,j,2)+ w0(i-1,j,2) ) &
            ) *dziq
      end do
    end do

    do j=2,j1
    do k=2,kmax
       do i=2,i1
          putout(i,j,k)  = putout(i,j,k)- (1./rhobf(k))*( &
              (rhobf(k) * putin(i,j,k) + rhobf(k+1) * putin(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
              -(rhobf(k) * putin(i,j,k) + rhobf(k-1) * putin(i,j,k-1) )*(w0(i,j,k )+w0(i-1,j,k )) &
                  )*dziq
        end do
      end do
    end do

  else

    do j=2,j1
      do i=2,i1
        putout(i,j,1)  = putout(i,j,1)- (1./rhobf(1))*( &
              ( rhobf(2) * putin(i,j,2)*dzf(1) + rhobf(1) * putin(i,j,1)*dzf(2) ) / dzh(2) &
                *( w0(i,j,2)+ w0(i-1,j,2) ))/ (4.*dzf(1))
      end do
    end do

    do j=2,j1
    do k=2,kmax
       do i=2,i1
          putout(i,j,k)  = putout(i,j,k)- (1./rhobf(k))*( &
                ( rhobf(k+1) * putin(i,j,k+1)*dzf(k) + rhobf(k) * putin(i,j,k)*dzf(k+1) ) / dzh(k+1) &
                  *( w0(i,j,k+1)+ w0(i-1,j,k+1) ) &
               -( rhobf(k) * putin(i,j,k)*dzf(k-1) + rhobf(k-1) * putin(i,j,k-1)*dzf(k) ) / dzh(k) &
                  *( w0(i,j,k)  + w0(i-1,j,k)   ) &
                )/ (4.*dzf(k))
        end do
      end do
    end do




  end if

end subroutine advecu_2nd


!> Advection at the v point.
subroutine advecv_2nd(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0, rhobf
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the v-field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin

  integer :: i,j,k

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rhoputin(i,j,k)=rhobf(k)*putin(i,j,k)
!      end do
!    end do
!  end do

  do k=1,kmax
    do j=2,j1
      do i=2,i1

        putout(i,j,k)  = putout(i,j,k)- ( &
              ( &
              ( u0(i+1,j,k)+u0(i+1,j-1,k))*(putin(i,j,k)+putin(i+1,j,k)) &
              -(u0(i ,j,k)+u0(i ,j-1,k))*(putin(i,j,k)+putin(i-1,j,k)) &
              )*dxiq &
              +( &
              ( v0(i,j+1,k)+v0(i,j,k))*(putin(i,j,k)+putin(i,j+1,k)) &
              -(v0(i,j-1,k)+v0(i,j,k))*(putin(i,j,k)+putin(i,j-1,k)) &
              )*dyiq )
      end do
    end do
  end do

  if (leq) then

    do j=2,j1
      do i=2,i1
        putout(i,j,1)  = putout(i,j,1)- (1./rhobf(1))*( &
           (w0(i,j,2)+w0(i,j-1,2))*(rhobf(2) * putin(i,j,2)+rhobf(1) * putin(i,j,1)) &
            )*dziq
      end do
    end do

    do j=2,j1
    do k=2,kmax
       do i=2,i1
          putout(i,j,k)  = putout(i,j,k)- (1./rhobf(k))*( &
                ( w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1) * putin(i,j,k+1) + rhobf(k) * putin(i,j,k)) &
                -(w0(i,j,k) +w0(i,j-1,k)) *(rhobf(k-1) * putin(i,j,k-1) + rhobf(k) * putin(i,j,k)) &
                )*dziq
        end do
      end do
    end do

  else
    do j=2,j1
      do i=2,i1
        putout(i,j,1)  = putout(i,j,1)- (1./rhobf(1))*( &
          (w0(i,j,2)+w0(i,j-1,2)) &
          *(rhobf(2) * putin(i,j,2)*dzf(1) + rhobf(1) * putin(i,j,1)*dzf(2) )/ dzh(2) &
          ) / (4. * dzf(1))
      end do
    end do

    do j=2,j1
    do k=2,kmax
       do i=2,i1
          putout(i,j,k)  = putout(i,j,k)- (1./rhobf(k))*( &
            (w0(i,j,k+1)+w0(i,j-1,k+1)) &
            *(rhobf(k) * putin(i,j,k+1)*dzf(k) + rhobf(k) * putin(i,j,k)*dzf(k+1) )/ dzh(k+1) &
            -(w0(i,j,k)+w0(i,j-1,k)) &
            *(rhobf(k-1) * putin(i,j,k-1)*dzf(k) + rhobf(k) * putin(i,j,k)*dzf(k-1)) / dzh(k) &
            ) / (4. * dzf(k))
        end do
      end do
    end do

  end if

end subroutine advecv_2nd



!> Advection at the w point.
subroutine advecw_2nd(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0, rhobh
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the w-field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency
!  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rhoputin

  integer :: i,j,k

!  do k=1,k1
!    do j=2-jh,j1+jh
!      do i=2-ih,i1+ih
!      rhoputin(i,j,k)=rhobh(k)*putin(i,j,k)
!      end do
!    end do
!  end do

  if (leq) then


    do k=2,kmax
      do j=2,j1
        do i=2,i1

          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                (putin(i+1,j,k)+putin(i,j,k))*(u0(i+1,j,k)+u0(i+1,j,k-1)) &
              -(putin(i-1,j,k)+putin(i,j,k))*(u0(i  ,j,k)+u0(i  ,j,k-1)) &
                )*dxiq &
              + &
                ( &
                (putin(i,j+1,k)+putin(i,j,k))*(v0(i,j+1,k)+v0(i,j+1,k-1)) &
              -(putin(i,j-1,k)+putin(i,j,k))*(v0(i,j  ,k)+v0(i,j  ,k-1)) &
                )*dyiq &
              + &
                (1./rhobh(k))*( &
                (rhobh(k) * putin(i,j,k) + rhobh(k+1) * putin(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
               -(rhobh(k) * putin(i,j,k) + rhobh(k-1) * putin(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                )*dziq &
                )

        end do
      end do
    end do
  else
    do k=2,kmax
      do j=2,j1
        do i=2,i1

          putout(i,j,k)  = putout(i,j,k) - (1./rhobh(k))*( &
                ( &
                ( rhobh(k) * putin(i+1,j,k) + rhobh(k) * putin(i,j,k) ) &
              *( dzf(k-1)*u0(i+1,j,k) + dzf(k)*u0(i+1,j,k-1) ) &
              -( rhobh(k) * putin(i,j,k) + rhobh(k) * putin(i-1,j,k) ) &
              *( dzf(k-1)*u0(i,j,k)+dzf(k)*u0(i ,j,k-1) ) &
                )*dxiq / dzh(k) &
              + &
                ( &
                ( rhobh(k) * putin(i,j+1,k) + rhobh(k) * putin(i,j,k) ) &
              *( dzf(k-1)*v0(i,j+1,k) + dzf(k)*v0(i,j+1,k-1) ) &
              -( rhobh(k) * putin(i,j,k) + rhobh(k) * putin(i,j-1,k) ) &
              *( dzf(k-1)*v0(i,j,k) + dzf(k)*v0(i,j,k-1) ) &
                ) *dyiq / dzh(k) &
              + &
                ( &
                ( rhobh(k) * putin(i,j,k) + rhobh(k+1) * putin(i,j,k+1) ) * (w0(i,j,k) + w0(i,j,k+1) ) &
               -( rhobh(k) * putin(i,j,k) + rhobh(k-1) * putin(i,j,k-1) ) * (w0(i,j,k) + w0(i,j,k-1) ) &
                ) / (4. *dzh(k) ) &
                )

        end do
      end do
    end do

  end if


end subroutine advecw_2nd
