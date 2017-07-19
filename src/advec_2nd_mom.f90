! Experimental advection routines based on advec_2nd
! advects all veloccity fields with a single pass
! explicitly referencing u0, up etc instead of putin, putout
! in the hope of better compiler optimizations

!advecu_2nd_mom_eq()  - equidistant grid
!advecu_2nd_mom_neq() - nonequidistant grid
!

subroutine advec_2nd_mom_eq()

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh
  use modfields, only : u0, v0, w0, up, vp, wp, rhobf, rhobh

  implicit none
  integer :: i,j,k

  k = 1
  do j=2,j1
     do i=2,i1
        up(i,j,1)  = up(i,j,1)-(1./rhobf(1))*( &
             ( rhobf(2) * u0(i,j,2) + rhobf(1) * u0(i,j,1))*( w0(i,j,2)+ w0(i-1,j,2) ) &
             ) *dziq

        vp(i,j,1)  = vp(i,j,1)- (1./rhobf(1))*( &
             (w0(i,j,2)+w0(i,j-1,2))*(rhobf(2) * v0(i,j,2)+rhobf(1) * v0(i,j,1)) &
             )*dziq 

        up(i,j,k)  = up(i,j,k)- ( &
             ( &
             (u0(i,j,k)+u0(i+1,j,k))*(u0(i,j,k)+u0(i+1,j,k)) &
             -(u0(i,j,k)+u0(i-1,j,k))*(u0(i,j,k)+u0(i-1,j,k)) &
             )*dxiq &
             +(  &
             (u0(i,j,k)+u0(i,j+1,k))*(v0(i,j+1,k)+v0(i-1,j+1 ,k)) &
             -(u0(i,j,k)+u0(i,j-1,k))*(v0(i,j  ,k)+v0(i-1,j  ,k)) &
             )*dyiq )


        vp(i,j,k)  = vp(i,j,k)- ( &
             ( &
             ( u0(i+1,j,k)+u0(i+1,j-1,k))*(v0(i,j,k)+v0(i+1,j,k)) &
             -(u0(i ,j,k)+u0(i ,j-1,k))*(v0(i,j,k)+v0(i-1,j,k)) &
             )*dxiq &
             +( &
             ( v0(i,j+1,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,j+1,k)) &
             -(v0(i,j-1,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,j-1,k)) &
             )*dyiq )

     end do
  end do


  do k=2,kmax
     do j=2,j1
        do i=2,i1
           up(i,j,k)  = up(i,j,k)- ( &
                ( &
                (u0(i,j,k)+u0(i+1,j,k))*(u0(i,j,k)+u0(i+1,j,k)) &
                -(u0(i,j,k)+u0(i-1,j,k))*(u0(i,j,k)+u0(i-1,j,k)) &
                )*dxiq &
                +(  &
                (u0(i,j,k)+u0(i,j+1,k))*(v0(i,j+1,k)+v0(i-1,j+1 ,k)) &
                -(u0(i,j,k)+u0(i,j-1,k))*(v0(i,j  ,k)+v0(i-1,j  ,k)) &
                )*dyiq )


           vp(i,j,k)  = vp(i,j,k)- ( &
                ( &
                ( u0(i+1,j,k)+u0(i+1,j-1,k))*(v0(i,j,k)+v0(i+1,j,k)) &
                -(u0(i ,j,k)+u0(i ,j-1,k))*(v0(i,j,k)+v0(i-1,j,k)) &
                )*dxiq &
                +( &
                ( v0(i,j+1,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,j+1,k)) &
                -(v0(i,j-1,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,j-1,k)) &
                )*dyiq )


           up(i,j,k)  = up(i,j,k)- (1./rhobf(k))*( &
                (rhobf(k) * u0(i,j,k) + rhobf(k+1) * u0(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
                -(rhobf(k) * u0(i,j,k) + rhobf(k-1) * u0(i,j,k-1) )*(w0(i,j,k )+w0(i-1,j,k )) &
                )*dziq

           vp(i,j,k)  = vp(i,j,k)- (1./rhobf(k))*( &
                ( w0(i,j,k+1)+w0(i,j-1,k+1))*(rhobf(k+1) * v0(i,j,k+1) + rhobf(k) * v0(i,j,k)) &
                -(w0(i,j,k) +w0(i,j-1,k)) *(rhobf(k-1) * v0(i,j,k-1) + rhobf(k) * v0(i,j,k)) &
                )*dziq

           wp(i,j,k)  = wp(i,j,k)- ( &
                ( &
                (w0(i+1,j,k)+w0(i,j,k))*(u0(i+1,j,k)+u0(i+1,j,k-1)) &
                -(w0(i-1,j,k)+w0(i,j,k))*(u0(i  ,j,k)+u0(i  ,j,k-1)) &
                )*dxiq &
                + &
                ( &
                (w0(i,j+1,k)+w0(i,j,k))*(v0(i,j+1,k)+v0(i,j+1,k-1)) &
                -(w0(i,j-1,k)+w0(i,j,k))*(v0(i,j  ,k)+v0(i,j  ,k-1)) &
                )*dyiq &
                + &
                (1./rhobh(k))*( &
                (rhobh(k) * w0(i,j,k) + rhobh(k+1) * w0(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                -(rhobh(k) * w0(i,j,k) + rhobh(k-1) * w0(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                )*dziq &
                )

        end do
     end do
  end do


end subroutine advec_2nd_mom_eq


subroutine advec_2nd_mom_neq()

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh
  use modfields, only : u0, v0, w0, up, vp, wp, rhobf, rhobh
  implicit none


  integer :: i,j,k

  k=1
  do j=2,j1
     do i=2,i1
        up(i,j,1)  = up(i,j,1)- (1./rhobf(1))*( &
             ( rhobf(2) * u0(i,j,2)*dzf(1) + rhobf(1) * u0(i,j,1)*dzf(2) ) / dzh(2) &
             *( w0(i,j,2)+ w0(i-1,j,2) ))/ (4.*dzf(1))

        vp(i,j,1)  = vp(i,j,1)- (1./rhobf(1))*( &
             (w0(i,j,2)+w0(i,j-1,2)) &
             *(rhobf(2) * v0(i,j,2)*dzf(1) + rhobf(1) * v0(i,j,1)*dzf(2) )/ dzh(2) &
             ) / (4. * dzf(1))


        up(i,j,k)  = up(i,j,k)- ( &
             ( &
             (u0(i,j,k)+u0(i+1,j,k))*(u0(i,j,k)+u0(i+1,j,k)) &
             -(u0(i,j,k)+u0(i-1,j,k))*(u0(i,j,k)+u0(i-1,j,k)) &
             )*dxiq &
             +(  &
             (u0(i,j,k)+u0(i,j+1,k))*(v0(i,j+1,k)+v0(i-1,j+1 ,k)) &
             -(u0(i,j,k)+u0(i,j-1,k))*(v0(i,j  ,k)+v0(i-1,j  ,k)) &
             )*dyiq )


        vp(i,j,k)  = vp(i,j,k)- ( &
             ( &
             ( u0(i+1,j,k)+u0(i+1,j-1,k))*(v0(i,j,k)+v0(i+1,j,k)) &
             -(u0(i ,j,k)+u0(i ,j-1,k))*(v0(i,j,k)+v0(i-1,j,k)) &
             )*dxiq &
             +( &
             ( v0(i,j+1,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,j+1,k)) &
             -(v0(i,j-1,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,j-1,k)) &
             )*dyiq )

     end do
  end do

  do k=2,kmax
     do j=2,j1
        do i=2,i1
           up(i,j,k)  = up(i,j,k)- (1./rhobf(k))*( &
                ( rhobf(k+1) * u0(i,j,k+1)*dzf(k) + rhobf(k) * u0(i,j,k)*dzf(k+1) ) / dzh(k+1) &
                *( w0(i,j,k+1)+ w0(i-1,j,k+1) ) &
                -( rhobf(k) * u0(i,j,k)*dzf(k-1) + rhobf(k-1) * u0(i,j,k-1)*dzf(k) ) / dzh(k) &
                *( w0(i,j,k)  + w0(i-1,j,k)   ) &
                )/ (4.*dzf(k))

           vp(i,j,k)  = vp(i,j,k)- (1./rhobf(k))*( &
                (w0(i,j,k+1)+w0(i,j-1,k+1)) &
                *(rhobf(k) * v0(i,j,k+1)*dzf(k) + rhobf(k) * v0(i,j,k)*dzf(k+1) )/ dzh(k+1) &
                -(w0(i,j,k)+w0(i,j-1,k)) &
                *(rhobf(k-1) * v0(i,j,k-1)*dzf(k) + rhobf(k) * v0(i,j,k)*dzf(k-1)) / dzh(k) &
                ) / (4. * dzf(k))

           wp(i,j,k)  = wp(i,j,k) - (1./rhobh(k))*( &
                ( &
                ( rhobh(k) * w0(i+1,j,k) + rhobh(k) * w0(i,j,k) ) &
                *( dzf(k-1)*u0(i+1,j,k) + dzf(k)*u0(i+1,j,k-1) ) &
                -( rhobh(k) * w0(i,j,k) + rhobh(k) * w0(i-1,j,k) ) &
                *( dzf(k-1)*u0(i,j,k)+dzf(k)*u0(i ,j,k-1) ) &
                )*dxiq / dzh(k) &
                + &
                ( &
                ( rhobh(k) * w0(i,j+1,k) + rhobh(k) * w0(i,j,k) ) &
                *( dzf(k-1)*v0(i,j+1,k) + dzf(k)*v0(i,j+1,k-1) ) &
                -( rhobh(k) * w0(i,j,k) + rhobh(k) * w0(i,j-1,k) ) &
                *( dzf(k-1)*v0(i,j,k) + dzf(k)*v0(i,j,k-1) ) &
                ) *dyiq / dzh(k) &
                + &
                ( &
                ( rhobh(k) * w0(i,j,k) + rhobh(k+1) * w0(i,j,k+1) ) * (w0(i,j,k) + w0(i,j,k+1) ) &
                -( rhobh(k) * w0(i,j,k) + rhobh(k-1) * w0(i,j,k-1) ) * (w0(i,j,k) + w0(i,j,k-1) ) &
                ) / (4. *dzh(k) ) &
                )

           up(i,j,k)  = up(i,j,k)- ( &
                ( &
                (u0(i,j,k)+u0(i+1,j,k))*(u0(i,j,k)+u0(i+1,j,k)) &
                -(u0(i,j,k)+u0(i-1,j,k))*(u0(i,j,k)+u0(i-1,j,k)) &
                )*dxiq &
                +(  &
                (u0(i,j,k)+u0(i,j+1,k))*(v0(i,j+1,k)+v0(i-1,j+1 ,k)) &
                -(u0(i,j,k)+u0(i,j-1,k))*(v0(i,j  ,k)+v0(i-1,j  ,k)) &
                )*dyiq )


           vp(i,j,k)  = vp(i,j,k)- ( &
                ( &
                ( u0(i+1,j,k)+u0(i+1,j-1,k))*(v0(i,j,k)+v0(i+1,j,k)) &
                -(u0(i ,j,k)+u0(i ,j-1,k))*(v0(i,j,k)+v0(i-1,j,k)) &
                )*dxiq &
                +( &
                ( v0(i,j+1,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,j+1,k)) &
                -(v0(i,j-1,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,j-1,k)) &
                )*dyiq )

        end do
     end do
  end do

end subroutine advec_2nd_mom_neq



