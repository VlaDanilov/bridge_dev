      module zaj_virtamps_l_m
        use helicities 
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'mxpart.f'
        include 'scale.f' !impure
        include 'epinv.f' !impure
        include 'epinv2.f' !impure

      contains

      subroutine zaj_virtamp_l1(i1,i2,i3,i4,i5,i6,za,zb,amps)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)
        complex(dp) :: amps_lc(2,2,2,2), amps_slc(2,2,2,2)


        call zaj_crossings_l(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_al_lc_pp, zaj_virt_al_lc_pm, amps_lc)

        call zaj_crossings_l(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_al_slc_pp, zaj_virt_al_slc_pm, amps_slc)

        amps = Nc*amps_lc + 1._dp/Nc*amps_slc

      end subroutine

      subroutine zaj_virtamp_l2(i1,i2,i3,i4,i5,i6,za,zb,amps)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        ! helicities for particles q, l, g, gam
        complex(dp), intent(out) :: amps(2,2,2,2)

        call zaj_crossings_l(i1,i2,i3,i4,i5,i6,za,zb,
     &     zaj_virt_al_fl_pp, zaj_virt_al_fl_pm, amps)

      end subroutine

      function zaj_virt_al_fl_pp(i1,i2,i3,i4,i5,i6,za,zb)
        include 'masses.f'
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_al_fl_pp

        complex(dp) :: L1

        real(dp) :: s,t
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_virt_al_fl_pp =
     &      + zb(i1,i3)*za(i2,i5)
     &        *( za(i5,i4)*zb(i4,i3) + za(i5,i6)*zb(i6,i3) )
     &        /(za(i4,i5)*za(i4,i6))
     &        *( L1(-s(i1,i2),-t(i4,i5,i6))/t(i4,i5,i6)**2
     &          -1._dp/(12._dp*t(i4,i5,i6)*mt**2) )

      end function

      function zaj_virt_al_fl_pm(i1,i2,i3,i4,i5,i6,za,zb)
        include 'masses.f'
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_al_fl_pm

        complex(dp) :: L1

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_virt_al_fl_pm =
     &      + zb(i1,i3)*zb(i3,i6)
     &        *( za(i2,i4)*zb(i4,i6) + za(i2,i5)*zb(i5,i6) )
     &        /(zb(i4,i5)*zb(i4,i6))
     &        *( L1(-s(i1,i2),-t(i4,i5,i6))/t(i4,i5,i6)**2
     &          -1._dp/(12._dp*t(i4,i5,i6)*mt**2) )

      end function


      ! newly reimplemented amplitudes below by CW

      function zaj_virt_al_lc_pp(i1,i2,i6,i5,i3,i4,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_al_lc_pp

        complex(dp) :: Atree, Vpole, Box(1), Boxint(1), rat, bubs, boxes
        complex(dp) :: lnrat,L0,L1,Lsm1

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        Atree =za(i2,i3)**2/(za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5))
        Vpole =-2*epinv**2 + epinv*(-1.5_dp - lnrat(musq,-s(i1,i6)) - 
     &       lnrat(musq,-s(i2,i6))) + 
     &    (-6 - lnrat(musq,-s(i1,i6))**2 - 3*lnrat(musq,-s(i2,i6)) - 
     &       lnrat(musq,-s(i2,i6))**2)/2._dp
      
        Vpole=Vpole*Atree

        Boxint(1)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))
        Box(1)=-Atree

        boxes = Boxint(1)*Box(1)

        Bubs=  -((L0(-s(i2,i6),-t(i3,i4,i5))*za(i1,i2)*za(i2,i3)*
     &         (za(i3,i4)*zb(i4,i1) + za(i3,i5)*zb(i5,i1)))/
     &       (t(i3,i4,i5)*za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5))) + 
     &    (L1(-s(i2,i6),-t(i3,i4,i5))*za(i1,i2)**2*
     &       (za(i3,i4)*zb(i4,i1) + za(i3,i5)*zb(i5,i1))*
     &       (-(za(i2,i3)*zb(i2,i1)) + za(i3,i6)*zb(i6,i1)))/
     &       (2._dp*t(i3,i4,i5)**2*za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5))
      
        Rat=czip

        zaj_virt_al_lc_pp = Vpole+boxes+bubs+rat

      end function

      function zaj_virt_al_lc_pm(i1,i2,i6,i5,i3,i4,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_al_lc_pm

        complex(dp) :: Atree, Vpole,Box(1),Boxint(1),rat,bubs,boxes
        complex(dp) :: lnrat,L0,L1,Lsm1

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      
        Atree = (za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4))**2/
     &    (t(i1,i2,i6)*za(i1,i6)*za(i2,i6)*zb(i5,i3)*zb(i5,i4))

        Vpole = -2*epinv**2 + epinv*
     &     (-1.5_dp - lnrat(musq,-s(i1,i6)) - lnrat(musq,-s(i2,i6))) + 
     &    (-6 - lnrat(musq,-s(i1,i6))**2 - 
     &       lnrat(musq,-s(i2,i6))**2 - 3*lnrat(musq,-t(i3,i4,i5)))/
     &       2._dp
        
        Vpole=Vpole*Atree
        Boxint(1)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))
        Box(1)=-Atree

        boxes=Boxint(1)*Box(1)         
        
        Bubs=  (L1(-s(i2,i6),-t(i3,i4,i5))*za(i1,i2)*zb(i4,i1)*zb(i6,i1)*
     &       (-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4)))/
     &     (2._dp*t(i3,i4,i5)**2*za(i1,i6)*zb(i5,i3)*zb(i5,i4)) + 
     &       (3*L0(-s(i2,i6),-t(i3,i4,i5))*zb(i6,i1)*(-(za(i1,i2)*zb(i4,i2))
     &   + za(i1,i6)*zb(i6,i4))*
     &       (za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4)))/
     &     (2._dp*t(i3,i4,i5)**2*za(i1,i6)*zb(i5,i3)*zb(i5,i4))
        
        Rat=
     &   -(za(i1,i2)*zb(i4,i1)*(za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4)))/
     &    (2._dp*t(i3,i4,i5)*za(i1,i6)*za(i2,i6)*zb(i5,i3)*zb(i5,i4))

        
        zaj_virt_al_lc_pm = Vpole+boxes+bubs+rat

      end function

      function zaj_virt_al_slc_pp(i1,i2,i6,i5,i3,i4,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_al_slc_pp

        complex(dp) :: Atree, Vpole,Box(2),Boxint(2),rat,bubs,boxes
        complex(dp) :: lnrat,L0,L1,Lsm1

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        Atree =za(i2,i3)**2/(za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5))

        Vpole=  (  epinv**2 + epinv*(1.5_dp + lnrat(musq,-s(i1,i2))) + 
     &       (6 + lnrat(musq,-s(i1,i2))**2 + 3*lnrat(musq,-t(i3,i4,i5)))/2)
        
        Vpole=Vpole*Atree
              
        Boxint(1)=Lsm1(-s(i1,i2),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))
        Boxint(2)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i1,i2),-t(i1,i2,i6))

        Box(1)=(za(i1,i2)**2*za(i3,i6)**2)
     &       /(za(i1,i6)**3*za(i2,i6)*za(i3,i5)*za(i4,i5))
        Box(2)=Atree

        Boxes = sum(box(:)*boxint(:))

        bubs= (L0(-t(i3,i4,i5),-s(i1,i2))*za(i1,i2)*za(i1,i3)*za(i3,i6)*
     &       zb(i6,i1))/(s(i1,i2)*za(i1,i6)**2*za(i3,i5)*za(i4,i5)) - 
     &    (L1(-t(i3,i4,i5),-s(i2,i6))*za(i1,i3)**2*za(i2,i6)*zb(i6,i1)**2)/
     &     (2._dp*s(i2,i6)**2*za(i1,i6)*za(i3,i5)*za(i4,i5)) - 
     &    (L0(-t(i3,i4,i5),-s(i2,i6))*
     &       (2*s(i1,i6)*za(i1,i3)*za(i2,i3) - 
     &         za(i1,i3)**2*za(i2,i6)*zb(i6,i1)))/
     &     (s(i2,i6)*za(i1,i6)**2*za(i3,i5)*za(i4,i5)) + 
     &    (L1(-t(i3,i4,i5),-s(i1,i2))*za(i1,i2)*za(i3,i6)*zb(i6,i1)*
     &       (-(za(i3,i4)*zb(i6,i4)) - za(i3,i5)*zb(i6,i5)))/
     &     (s(i1,i2)**2*za(i1,i6)*za(i3,i5)*za(i4,i5))

        Rat=   (za(i1,i3)*(za(i3,i4)*zb(i4,i2) + za(i3,i5)*zb(i5,i2))*zb(i6,i1)**2)/
     &     (2._dp*t(i3,i4,i5)*za(i1,i6)*za(i3,i5)*za(i4,i5)*zb(i2,i1)*zb(i6,i2)) + 
     &    ((za(i3,i4)*zb(i4,i1) + za(i3,i5)*zb(i5,i1))*
     &       (za(i1,i3)*zb(i6,i1) + za(i2,i3)*zb(i6,i2)))/
     &     (2._dp*t(i3,i4,i5)*za(i1,i6)*za(i3,i5)*za(i4,i5)*zb(i2,i1)) + 
     &    (za(i2,i3)*(za(i3,i4)*zb(i4,i1) + za(i3,i5)*zb(i5,i1))*
     &       (za(i1,i2)*zb(i2,i1) + za(i2,i6)*zb(i6,i2)))/
     &     (2._dp*t(i3,i4,i5)*za(i1,i6)*za(i2,i6)*za(i3,i5)*za(i4,i5)*zb(i2,i1))
        
        zaj_virt_al_slc_pp = Vpole+boxes+bubs+rat

      end function

      function zaj_virt_al_slc_pm(i1,i2,i6,i5,i3,i4,za,zb)
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zaj_virt_al_slc_pm

        complex(dp) :: Atree, Vpole,Box(2),Boxint(2),rat,bubs,boxes
        complex(dp) :: lnrat,L0,L1,Lsm1
        
        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      
        Atree =(za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4))**2/
     &    (t(i1,i2,i6)*za(i1,i6)*za(i2,i6)*zb(i5,i3)*zb(i5,i4))

        
        Vpole=  (  epinv**2 + epinv*(1.5_dp + lnrat(musq,-s(i1,i2))) + 
     &       (6 + lnrat(musq,-s(i1,i2))**2 + 3*lnrat(musq,-t(i3,i4,i5)))/2)
        
        Vpole=Vpole*Atree
        
        Boxint(1)=Lsm1(-s(i1,i2),-t(i1,i2,i6),-s(i2,i6),-t(i1,i2,i6))
        Boxint(2)=Lsm1(-s(i1,i6),-t(i1,i2,i6),-s(i1,i2),-t(i1,i2,i6))

        Box(1)= -((za(i1,i2)**2*(za(i1,i6)*zb(i4,i1) + 
     -          za(i2,i6)*zb(i4,i2))*
     -        (za(i3,i6)*zb(i4,i3) - za(i5,i6)*zb(i5,i4)))/
     -      (t(i3,i4,i5)*za(i1,i6)**3*za(i2,i6)*zb(i5,i3)*zb(i5,i4)))
        Box(2)=Atree

        
        Boxes=sum(Box(:)*Boxint(:))

        bubs=  -((L1(-s(i1,i2),-t(i3,i4,i5))*za(i1,i2)**2*
     &         (za(i1,i6)*zb(i4,i1) + za(i2,i6)*zb(i4,i2))**2*zb(i6,i1)*
     &         zb(i6,i2))/
     &       (t(i3,i4,i5)**3*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4))) - 
     &    (2*L0(-s(i1,i2),-t(i3,i4,i5))*za(i1,i2)*
     &       (za(i1,i6)*zb(i4,i1) + za(i2,i6)*zb(i4,i2))*zb(i6,i1)*
     &       (-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4)))/
     &     (t(i3,i4,i5)**2*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4)) + 
     &    (L1(-t(i3,i4,i5),-s(i2,i6))*za(i2,i6)*
     &       (-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*zb(i6,i1)**2*
     &       (-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4)))/
     &     (2.*s(i2,i6)**2*t(i3,i4,i5)*za(i1,i6)*zb(i5,i3)*zb(i5,i4)) + 
     &    (L0(-t(i3,i4,i5),-s(i2,i6))*
     &       (2*s(i1,i6)*(-(za(i2,i3)*zb(i4,i3)) + 
     &            za(i2,i5)*zb(i5,i4)) - 
     &         za(i2,i6)*(-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*
     &          zb(i6,i1))*(-(za(i1,i2)*zb(i4,i2)) + za(i1,i6)*zb(i6,i4))
     &       )/(s(i2,i6)*t(i3,i4,i5)*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4))

        Rat=(za(i1,i2)**2*zb(i4,i2)*(za(i1,i6)*zb(i4,i1) + za(i2,i6)*zb(i4,i2))*zb(i6,i1))/
     &     (2._dp*t(i3,i4,i5)**2*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4)) + 
     &    (za(i1,i2)*(za(i1,i6)*zb(i4,i1) + za(i2,i6)*zb(i4,i2))*
     &       (-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*zb(i6,i1))/
     &     (2._dp*t(i3,i4,i5)**2*za(i1,i6)**2*zb(i5,i3)*zb(i5,i4)) + 
     &    (za(i1,i2)*zb(i2,i1)*(-(za(i1,i3)*zb(i4,i3)) + za(i1,i5)*zb(i5,i4))*zb(i6,i1)*
     &       zb(i6,i4))/(2.*t(i3,i4,i5)**2*za(i1,i6)*zb(i5,i3)*zb(i5,i4)*zb(i6,i2)) - 
     &    (zb(i6,i1)*zb(i6,i4)**2)/(2._dp*t(i3,i4,i5)*zb(i5,i3)*zb(i5,i4)*zb(i6,i2)) + 
     &    (zb(i6,i4)*(za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4)))/
     &     (2._dp*t(i3,i4,i5)*za(i1,i6)*zb(i5,i3)*zb(i5,i4)) + 
     &    (za(i1,i2)*zb(i4,i1) + za(i2,i6)*zb(i6,i4))**2/
     &       (2._dp*t(i3,i4,i5)*za(i1,i6)*za(i2,i6)*zb(i5,i3)*zb(i5,i4))
        
        zaj_virt_al_slc_pm = Vpole+boxes+bubs+rat

      end function

      end module
