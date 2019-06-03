c ============
c       basis tree amplitudes, for (gluon,gamma) = pp,pm,mp helicities
c ============

      module zaj_treeamps_m
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex(dp), private :: im = (0._dp,1._dp)
      contains

      pure function zaj_tree_isr_pp(i1,i2,i3,i4,i5,i6,za,zb)
        complex(dp) :: zaj_tree_isr_pp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zaj_tree_isr_pp = (-za(i1,i2)*za(i2,i5)**2)/
     &       (za(i1,i3)*za(i1,i4)*za(i2,i3)*za(i2,i4)*za(i5,i6))

      end function

      pure function zaj_tree_isr_mp(i1,i2,i3,i4,i5,i6,za,zb)
        complex(dp) :: zaj_tree_isr_mp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_tree_isr_mp = 
     &  (za(i2,i3)**2*zb(i6,i1)**2)/(t(i2,i3,i4)*za(i2,i4)*
     &  (za(i4,i5)*zb(i5,i1)+za(i4,i6)*zb(i6,i1))*zb(i6,i5)) +
     &  (za(i2,i3)*zb(i2,i1)*(za(i5,i2)*zb(i2,i1)
     &            +za(i5,i3)*zb(i3,i1))**2)/
     &  (s(i2,i3)*za(i5,i6)*(za(i4,i2)*zb(i2,i1)+za(i4,i3)*zb(i3,i1))*
     &  (za(i4,i1)*zb(i1,i2)+za(i4,i3)*zb(i3,i2))*zb(i3,i1)) +
     &  (zb(i4,i1)*(za(i3,i1)*zb(i1,i6)+za(i3,i4)*zb(i4,i6))**2)/
     &  (t(i1,i3,i4)*s(i1,i4)*(za(i4,i1)*zb(i1,i2)
     &       +za(i4,i3)*zb(i3,i2))*zb(i6,i5))

      end function

      pure function zaj_tree_isr_pm(i1,i2,i3,i4,i5,i6,za,zb)
        complex(dp) :: zaj_tree_isr_pm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_tree_isr_pm = 
     &   (za(i2,i4)**2*zb(i6,i1)**2)/(t(i2,i3,i4)*za(i2,i3)*
     &   (za(i3,i5)*zb(i5,i1)+za(i3,i6)*zb(i6,i1))*zb(i6,i5)) +
     &   (za(i2,i4)*zb(i2,i1)*(za(i5,i2)*zb(i2,i1)+za(i5,i4)*
     &     zb(i4,i1))**2)/
     &   (s(i2,i4)*za(i5,i6)*(za(i3,i2)*zb(i2,i1)
     &       +za(i3,i4)*zb(i4,i1))*
     &   (za(i3,i1)*zb(i1,i2)+za(i3,i4)*zb(i4,i2))*zb(i4,i1)) +
     &   (zb(i3,i1)*(za(i4,i1)*zb(i1,i6)+za(i4,i3)*zb(i3,i6))**2)/
     &   (t(i1,i3,i4)*s(i1,i3)*(za(i3,i1)*zb(i1,i2)+za(i3,i4)*
     &     zb(i4,i2))*zb(i6,i5))

      end function

      pure function zaj_tree_fsr_pp(i1,i2,i3,i4,i5,i6,za,zb)
        complex(dp) :: zaj_tree_fsr_pp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        zaj_tree_fsr_pp = 
     &     za(i2,i5)**2/(za(i1,i3)*za(i2,i3)*za(i4,i5)*za(i4,i6))

      end function

      pure function zaj_tree_fsr_pm(i1,i2,i3,i4,i5,i6,za,zb)
        complex(dp) :: zaj_tree_fsr_pm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)

        complex(dp) :: s,t
        s(i1,i2) = za(i1,i2)*zb(i2,i1)
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        zaj_tree_fsr_pm =
     &    (za(i2,i4)*zb(i4,i6)+za(i2,i5)*zb(i5,i6))**2/
     &    (t(i4,i5,i6)*za(i1,i3)*za(i2,i3)*zb(i4,i5)*zb(i4,i6))

      end function

c ===================== 
c  anomalous couplings
c ===================== 

      function zaj_tree_anomZZ_pp(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'

        complex(dp) :: zaj_tree_anomZZ_pp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, s56, hSign
        complex(dp) :: zbza
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        zbza(i1,i2,i6,i5)=+za(i1,i2)*zb(i2,i5)+za(i1,i6)*zb(i6,i5)

        s56 = s(i5,i6)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomZZ_pp = -1._dp/(4._dp*s56*za(i2,i3)*za(i1,i3))
     & *((hSign*im*hitZ(1)+hitZ(3))*2._dp*za(i2,i5)*zb(i6,i4)*zbza(i2,i1,i3,i4)
     &  +(hSign*im*hitZ(2)+hitZ(4))*za(i5,i4)*zb(i4,i6)*zbza(i2,i1,i3,i4)**2)

      end function

      function zaj_tree_anomZZ_pm(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'
        complex(dp) :: zaj_tree_anomZZ_pm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, s56, s456, hSign
        complex(dp) :: zazb
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        zazb(i1,i2,i6,i5)=+zb(i1,i2)*za(i2,i5)+zb(i1,i6)*za(i6,i5)


        s56 = s(i5,i6)
        s456 = s(i6,i5)+s(i6,i4)+s(i5,i4)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomZZ_pm = -1._dp/(4._dp*s56*za(i1,i3)*za(i2,i3))
     & *((-hSign*im*hitZ(1)+hitZ(3))*2._dp*za(i5,i4)*za(i2,i4)*zazb(i6,i1,i3,i2)
     &  +(-hSign*im*hitZ(2)+hitZ(4))*s456*zb(i6,i4)*za(i4,i5)*za(i2,i4)**2)

      end function

      function zaj_tree_anomZZ_mp(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'

        complex(dp) :: zaj_tree_anomZZ_mp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, s56, s456, hSign
        complex(dp) :: zbza
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        zbza(i1,i2,i6,i5)=+za(i1,i2)*zb(i2,i5)+za(i1,i6)*zb(i6,i5)


        s56 = s(i5,i6)
        s456 = s(i6,i5)+s(i6,i4)+s(i5,i4)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomZZ_mp = -1._dp/(4._dp*s56*zb(i2,i3)*zb(i1,i3))
     & *((hSign*im*hitZ(1)+hitZ(3))*2._dp*zb(i6,i4)*zb(i1,i4)*zbza(i5,i2,i3,i1)
     &  +(hSign*im*hitZ(2)+hitZ(4))*s456*za(i5,i4)*zb(i4,i6)*zb(i1,i4)**2)

      end function

      function zaj_tree_anomZa_pp(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'

        complex(dp) :: zaj_tree_anomZa_pp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, hSign, s123
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

        s123=s(i1,i2)+s(i1,i3)+s(i2,i3)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomZa_pp =  1._dp/(4._dp*s123*za(i1,i3)*za(i2,i3))
     & *zb(i6,i4)*(za(i2,i1)*zb(i1,i4)+za(i2,i3)*zb(i3,i4))
     & *((hSign*im*hitgam(1)+hitgam(3))*2._dp*za(i2,i5)
     &  +(hSign*im*hitgam(2)+hitgam(4))*za(i5,i6)*zb(i6,i4)*za(i2,i4))

      end function

      function zaj_tree_anomZa_pm(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'

        complex(dp) :: zaj_tree_anomZa_pm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, s123, hSign
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

        s123=s(i1,i2)+s(i1,i3)+s(i2,i3)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomZa_pm =  1._dp/(4._dp*s123*za(i1,i3)*za(i2,i3))
     & *za(i5,i4)*za(i2,i4)
     & *((-hSign*im*hitgam(1)+hitgam(3))
     &   *2._dp*(zb(i6,i1)*za(i1,i2)+zb(i6,i3)*za(i3,i2))
     &  +(-hSign*im*hitgam(2)+hitgam(4))
     &   *zb(i6,i5)*za(i5,i4)*(zb(i4,i1)*za(i1,i2)+zb(i4,i3)*za(i3,i2)))

      end function

      function zaj_tree_anomZa_mp(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'
        complex(dp) :: zaj_tree_anomZa_mp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, s123, hSign
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        s123=s(i1,i2)+s(i1,i3)+s(i2,i3)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomZa_mp = 1._dp/(4._dp*s123*zb(i1,i3)*zb(i2,i3))
     & *zb(i6,i4)*zb(i1,i4)
     & *((hSign*im*hitgam(1)+hitgam(3))
     &   *2._dp*(za(i5,i2)*zb(i2,i1)+za(i5,i3)*zb(i3,i1))
     &  +(hSign*im*hitgam(2)+hitgam(4))
     &   *za(i5,i6)*zb(i6,i4)*(za(i4,i2)*zb(i2,i1)+za(i4,i3)*zb(i3,i1)))

      end function

      function zaj_tree_anomaZ_pp(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'

        complex(dp) :: zaj_tree_anomaZ_pp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, s56, hSign
        complex(dp) :: zbza
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        zbza(i1,i2,i6,i5)=+za(i1,i2)*zb(i2,i5)+za(i1,i6)*zb(i6,i5)

        s56 = s(i5,i6)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomaZ_pp = -1._dp/(4._dp*s56*za(i2,i3)*za(i1,i3))
     & *((hSign*im*hitgam(1)+hitgam(3))*2._dp*za(i2,i5)*zb(i6,i4)*zbza(i2,i1,i3,i4)
     &  +(hSign*im*hitgam(2)+hitgam(4))*za(i5,i4)*zb(i4,i6)*zbza(i2,i1,i3,i4)**2)

      end function

      function zaj_tree_anomaZ_pm(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'

        complex(dp) :: zaj_tree_anomaZ_pm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, s56, s456, hSign
        complex(dp) :: zazb
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        zazb(i1,i2,i6,i5)=+zb(i1,i2)*za(i2,i5)+zb(i1,i6)*za(i6,i5)


        s56 = s(i5,i6)
        s456 = s(i6,i5)+s(i6,i4)+s(i5,i4)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomaZ_pm = -1._dp/(4._dp*s56*za(i1,i3)*za(i2,i3))
     & *((-hSign*im*hitgam(1)+hitgam(3))*2._dp*za(i5,i4)*za(i2,i4)*zazb(i6,i1,i3,i2)
     &  +(-hSign*im*hitgam(2)+hitgam(4))*s456*zb(i6,i4)*za(i4,i5)*za(i2,i4)**2)

      end function

      function zaj_tree_anomaZ_mp(i1,i2,i3,i4,i5,i6,za,zb,flip)
        include 'anomcoup.f'
        complex(dp) :: zaj_tree_anomaZ_mp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        logical, intent(in) :: flip

        real(dp) :: s, s56, s456, hSign
        complex(dp) :: zbza
        
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        zbza(i1,i2,i6,i5)=+za(i1,i2)*zb(i2,i5)+za(i1,i6)*zb(i6,i5)


        s56 = s(i5,i6)
        s456 = s(i6,i5)+s(i6,i4)+s(i5,i4)

        hSign = 1._dp
        if (flip) then
          hSign = -1._dp
        endif

        zaj_tree_anomaZ_mp = -1._dp/(4._dp*s56*zb(i2,i3)*zb(i1,i3))
     & *((hSign*im*hitgam(1)+hitgam(3))*2._dp*zb(i6,i4)*zb(i1,i4)*zbza(i5,i2,i3,i1)
     &  +(hSign*im*hitgam(2)+hitgam(4))*s456*za(i5,i4)*zb(i4,i6)*zb(i1,i4)**2)


      end function

      end module
