      module hjetmass_hel
          implicit none

          logical, parameter, private :: approx = .false.
          !! disables or enables the HTL approximation
          !! in the mt exact routines

          contains

      subroutine a_1l_mtex(za,zb,i1,i2,i3,amtex)
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'constants.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'

          complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          real(dp), intent(out) :: amtex(3)

          real(dp) :: sman,tman,uman
          real(dp) :: aex_prefactor

          real(dp) :: s
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

          sman = s(i1,i2)
          tman = s(i1,i3)
          uman = s(i2,i3)
          
          amtex = 0._dp

          amtex(1) = 16/(3._dp*sman)
          amtex(2) = 0.8_dp + (14*tman)/(45._dp*sman) + (14*uman)/(45._dp*sman)
          amtex(3) = 
     -       (16*sman)/105._dp + (4*tman)/35._dp + (2*tman**2)/(63._dp*sman) + 
     -       (4*uman)/35._dp + (4*tman*uman)/(63._dp*sman) + 
     -       (2*uman**2)/(63._dp*sman)

          aex_prefactor = sqrt(gsq)**3 / 32d0 / pi**2 / sqrt(vevsq)

        ! our normalization
        ! in HTL the 1-loop coefficient should be 1
          amtex(:) = amtex(:) / amtex(1)

        ! compatibility with other amps
          amtex = amtex * aex_prefactor / s(i1,i2) * 8._dp/3._dp
      end subroutine

      subroutine c_1l_mtex(za,zb,i1,i2,i3,cmtex)
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'constants.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'

          complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          real(dp), intent(out) :: cmtex(3,4)

          real(dp) :: c1lomt0,c1lomt2,c1lomt4
          real(dp) :: c2lomt0,c2lomt2,c2lomt4
          real(dp) :: c3lomt0,c3lomt2,c3lomt4
          real(dp) :: c4lomt0,c4lomt2,c4lomt4

          real(dp) :: sman,tman,uman
          real(dp) :: cex_prefactor

          real(dp) :: s
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

          sman = s(i1,i2)
          tman = s(i1,i3)
          uman = s(i2,i3)
          
          cmtex = 0._dp

      c1lomt0 = (0.4D1 / 0.3D1) * sman ** 2 + (0.8D1 / 0.3D1) * tman * s
     #man + (0.4D1 / 0.3D1) * tman ** 2 + (0.8D1 / 0.3D1) * sman * uman 
     #+ (0.8D1 / 0.3D1) * tman * uman + (0.4D1 / 0.3D1) * uman ** 2
      c2lomt0 = (0.4D1 / 0.3D1) * sman ** 2
      c3lomt0 = (0.4D1 / 0.3D1) * tman ** 2
      c4lomt0 = (0.4D1 / 0.3D1) * uman ** 2

       cmtex(1,1) = c1lomt0
       cmtex(1,2) = c2lomt0
       cmtex(1,3) = c3lomt0
       cmtex(1,4) = c4lomt0

      c1lomt2 = (0.7D1 / 0.90D2) * sman ** 3 + (0.7D1 / 0.30D2) * sman *
     #* 2 * tman + (0.7D1 / 0.30D2) * sman * tman ** 2 + (0.7D1 / 0.90D2
     #) * tman ** 3 + (0.7D1 / 0.30D2) * sman ** 2 * uman + (0.2D1 / 0.5
     #D1) * sman * tman * uman + (0.7D1 / 0.30D2) * tman ** 2 * uman + (
     #0.7D1 / 0.30D2) * sman * uman ** 2 + (0.7D1 / 0.30D2) * tman * uma
     #n ** 2 + (0.7D1 / 0.90D2) * uman ** 3
      c2lomt2 = (0.7D1 / 0.90D2) * sman ** 3 + (0.7D1 / 0.90D2) * sman *
     #* 2 * tman + (0.7D1 / 0.90D2) * sman ** 2 * uman
      c3lomt2 = (0.7D1 / 0.90D2) * sman * tman ** 2 + (0.7D1 / 0.90D2) *
     # tman ** 3 + (0.7D1 / 0.90D2) * tman ** 2 * uman
      c4lomt2 = (0.7D1 / 0.90D2) * sman * uman ** 2 + (0.7D1 / 0.90D2) *
     # tman * uman ** 2 + (0.7D1 / 0.90D2) * uman ** 3

       cmtex(2,1) = c1lomt2
       cmtex(2,2) = c2lomt2
       cmtex(2,3) = c3lomt2
       cmtex(2,4) = c4lomt2

      c1lomt4 = sman ** 4 / 126 + (0.2D1 / 0.63D2) * sman ** 3 * tman + 
     #sman ** 2 * tman ** 2 / 21 + (0.2D1 / 0.63D2) * sman * tman ** 3 +
     # tman ** 4 / 126 + (0.2D1 / 0.63D2) * sman ** 3 * uman + (0.13D2 /
     # 0.140D3) * sman ** 2 * tman * uman + (0.13D2 / 0.140D3) * sman * 
     #tman ** 2 * uman + (0.2D1 / 0.63D2) * tman ** 3 * uman + sman ** 2
     # * uman ** 2 / 21 + (0.13D2 / 0.140D3) * sman * tman * uman ** 2 +
     # tman ** 2 * uman ** 2 / 21 + (0.2D1 / 0.63D2) * sman * uman ** 3 
     #+ (0.2D1 / 0.63D2) * tman * uman ** 3 + uman ** 4 / 126
      c2lomt4 = sman ** 4 / 126 + sman ** 3 * tman / 63 + sman ** 2 * tm
     #an ** 2 / 126 + sman ** 3 * uman / 63 + (0.23D2 / 0.1260D4) * sman
     # ** 2 * tman * uman + sman ** 2 * uman ** 2 / 126
      c3lomt4 = sman ** 2 * tman ** 2 / 126 + sman * tman ** 3 / 63 + tm
     #an ** 4 / 126 + (0.23D2 / 0.1260D4) * sman * tman ** 2 * uman + tm
     #an ** 3 * uman / 63 + tman ** 2 * uman ** 2 / 126
      c4lomt4 = sman ** 2 * uman ** 2 / 126 + (0.23D2 / 0.1260D4) * sman
     # * tman * uman ** 2 + tman ** 2 * uman ** 2 / 126 + sman * uman **
     # 3 / 63 + tman * uman ** 3 / 63 + uman ** 4 / 126

       cmtex(3,1) = c1lomt4
       cmtex(3,2) = c2lomt4
       cmtex(3,3) = c3lomt4
       cmtex(3,4) = c4lomt4

        cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

        cmtex(:,:) = cmtex(:,:) * cex_prefactor / 16._dp

      end subroutine

      subroutine c_2l_mtex(za,zb,i1,i2,i3,cmtex)
      !! one-loop and two-loop qqg amplitude coefficient in 1/mt expansion
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'constants.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          include 'scale.f'
          include 'asymptotic.f'
          include 'nf.f'

          complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp) :: cmtex(2, 4)

          real(dp) :: mH
          real(dp) :: sman,tman,uman
          complex(dp) :: cdlogwrap, dilogc
          complex(dp) :: LogMums, LogMumt, LogMumu
          real(dp) :: LogMuMtop
          real(dp) :: pcut,nl,kg

          complex(dp) :: c1lomt0,c1lomt2,c1lomt4
          complex(dp) :: c2lomt0,c2lomt2,c2lomt4
          complex(dp) :: c3lomt0,c3lomt2,c3lomt4
          complex(dp) :: c4lomt0,c4lomt2,c4lomt4
          complex(dp) :: c1mt0,c1mt2,c1mt4
          complex(dp) :: c2mt0,c2mt2,c2mt4
          complex(dp) :: c3mt0,c3mt2,c3mt4
          complex(dp) :: c4mt0,c4mt2,c4mt4

          complex(dp) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
          complex(dp) :: t11,t12,t13,t14,t15,t16,t17,t18
          complex(dp) :: t19,t20,t21,t22,t23,t24,t25,t26
          complex(dp) :: t27,t28,t29,t30,t31,t32,t33,t34
          complex(dp) :: t35,t36,t37,t38,t39,t40

          real(dp) :: cex_prefactor

          real(dp) :: s
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

          sman = s(i1,i2)
          tman = s(i1,i3)
          uman = s(i2,i3)

          mH = hmass

          LogMums = cdlogwrap(-musq/sman)
          LogMumt = cdlogwrap(-musq/tman)
          LogMumu = cdlogwrap(-musq/uman)
          LogMuMtop = log(musq/mt**2)

          pcut = 1d0
          nl =  5d0
          kg = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
     &            - ca*dlog(pcut)**2
     &            + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))
         

          c1lomt0 = (0d0,0d0)
          c1lomt2 = (0d0,0d0)
          c1lomt4 = (0d0,0d0)
          c2lomt0 = (0d0,0d0)
          c2lomt2 = (0d0,0d0)
          c2lomt4 = (0d0,0d0)
          c3lomt0 = (0d0,0d0)
          c3lomt2 = (0d0,0d0)
          c3lomt4 = (0d0,0d0)
          c4lomt0 = (0d0,0d0)
          c4lomt2 = (0d0,0d0)
          c4lomt4 = (0d0,0d0)

          c1mt0 = (0d0,0d0)
          c1mt2 = (0d0,0d0)
          c1mt4 = (0d0,0d0)
          c2mt0 = (0d0,0d0)
          c2mt2 = (0d0,0d0)
          c2mt4 = (0d0,0d0)
          c3mt0 = (0d0,0d0)
          c3mt2 = (0d0,0d0)
          c3mt4 = (0d0,0d0)
          c4mt0 = (0d0,0d0)
          c4mt2 = (0d0,0d0)
          c4mt4 = (0d0,0d0)


       if (mtex >= 0) then
      c1lomt0 = (0.4D1 / 0.3D1) * sman ** 2 + (0.8D1 / 0.3D1) * tman * s
     #man + (0.4D1 / 0.3D1) * tman ** 2 + (0.8D1 / 0.3D1) * sman * uman 
     #+ (0.8D1 / 0.3D1) * tman * uman + (0.4D1 / 0.3D1) * uman ** 2
      c2lomt0 = (0.4D1 / 0.3D1) * sman ** 2
      c3lomt0 = (0.4D1 / 0.3D1) * tman ** 2
      c4lomt0 = (0.4D1 / 0.3D1) * uman ** 2
       end if

       if (mtex >= 2) then
      c1lomt2 = (0.7D1 / 0.90D2) * sman ** 3 + (0.7D1 / 0.30D2) * sman *
     #* 2 * tman + (0.7D1 / 0.30D2) * sman * tman ** 2 + (0.7D1 / 0.90D2
     #) * tman ** 3 + (0.7D1 / 0.30D2) * sman ** 2 * uman + (0.2D1 / 0.5
     #D1) * sman * tman * uman + (0.7D1 / 0.30D2) * tman ** 2 * uman + (
     #0.7D1 / 0.30D2) * sman * uman ** 2 + (0.7D1 / 0.30D2) * tman * uma
     #n ** 2 + (0.7D1 / 0.90D2) * uman ** 3
      c2lomt2 = (0.7D1 / 0.90D2) * sman ** 3 + (0.7D1 / 0.90D2) * sman *
     #* 2 * tman + (0.7D1 / 0.90D2) * sman ** 2 * uman
      c3lomt2 = (0.7D1 / 0.90D2) * sman * tman ** 2 + (0.7D1 / 0.90D2) *
     # tman ** 3 + (0.7D1 / 0.90D2) * tman ** 2 * uman
      c4lomt2 = (0.7D1 / 0.90D2) * sman * uman ** 2 + (0.7D1 / 0.90D2) *
     # tman * uman ** 2 + (0.7D1 / 0.90D2) * uman ** 3
       end if

       if (mtex >= 4) then
      c1lomt4 = sman ** 4 / 126 + (0.2D1 / 0.63D2) * sman ** 3 * tman + 
     #sman ** 2 * tman ** 2 / 21 + (0.2D1 / 0.63D2) * sman * tman ** 3 +
     # tman ** 4 / 126 + (0.2D1 / 0.63D2) * sman ** 3 * uman + (0.13D2 /
     # 0.140D3) * sman ** 2 * tman * uman + (0.13D2 / 0.140D3) * sman * 
     #tman ** 2 * uman + (0.2D1 / 0.63D2) * tman ** 3 * uman + sman ** 2
     # * uman ** 2 / 21 + (0.13D2 / 0.140D3) * sman * tman * uman ** 2 +
     # tman ** 2 * uman ** 2 / 21 + (0.2D1 / 0.63D2) * sman * uman ** 3 
     #+ (0.2D1 / 0.63D2) * tman * uman ** 3 + uman ** 4 / 126
      c2lomt4 = sman ** 4 / 126 + sman ** 3 * tman / 63 + sman ** 2 * tm
     #an ** 2 / 126 + sman ** 3 * uman / 63 + (0.23D2 / 0.1260D4) * sman
     # ** 2 * tman * uman + sman ** 2 * uman ** 2 / 126
      c3lomt4 = sman ** 2 * tman ** 2 / 126 + sman * tman ** 3 / 63 + tm
     #an ** 4 / 126 + (0.23D2 / 0.1260D4) * sman * tman ** 2 * uman + tm
     #an ** 3 * uman / 63 + tman ** 2 * uman ** 2 / 126
      c4lomt4 = sman ** 2 * uman ** 2 / 126 + (0.23D2 / 0.1260D4) * sman
     # * tman * uman ** 2 + tman ** 2 * uman ** 2 / 126 + sman * uman **
     # 3 / 63 + tman * uman ** 3 / 63 + uman ** 4 / 126
       end if

       cmtex(:,:) = 0._dp

       if (mtex >= 0) then
           include 'ggg/c1_mt0.f'
           include 'ggg/c2_mt0.f'
           include 'ggg/c3_mt0.f'
           include 'ggg/c4_mt0.f'

           cmtex(2, 1) = c1mt0
           cmtex(2, 2) = c2mt0
           cmtex(2, 3) = c3mt0
           cmtex(2, 4) = c4mt0

           cmtex(1, 1) = c1lomt0
           cmtex(1, 2) = c2lomt0
           cmtex(1, 3) = c3lomt0
           cmtex(1, 4) = c4lomt0
        endif

       if (mtex >= 2) then
           include 'ggg/c1_mt2.f'
           include 'ggg/c2_mt2.f'
           include 'ggg/c3_mt2.f'
           include 'ggg/c4_mt2.f'

           cmtex(2, 1) = cmtex(2, 1) + c1mt2/mt**2
           cmtex(2, 2) = cmtex(2, 2) + c2mt2/mt**2
           cmtex(2, 3) = cmtex(2, 3) + c3mt2/mt**2
           cmtex(2, 4) = cmtex(2, 4) + c4mt2/mt**2

           cmtex(1, 1) = cmtex(1, 1) + c1lomt2/mt**2
           cmtex(1, 2) = cmtex(1, 2) + c2lomt2/mt**2
           cmtex(1, 3) = cmtex(1, 3) + c3lomt2/mt**2
           cmtex(1, 4) = cmtex(1, 4) + c4lomt2/mt**2
        endif

       if (mtex >= 4) then
           include 'ggg/c1_mt4.f'
           include 'ggg/c2_mt4.f'
           include 'ggg/c3_mt4.f'
           include 'ggg/c4_mt4.f'

           cmtex(2, 1) = cmtex(2, 1) + c1mt4/mt**4
           cmtex(2, 2) = cmtex(2, 2) + c2mt4/mt**4
           cmtex(2, 3) = cmtex(2, 3) + c3mt4/mt**4
           cmtex(2, 4) = cmtex(2, 4) + c4mt4/mt**4

           cmtex(1, 1) = cmtex(1, 1) + c1lomt4/mt**4
           cmtex(1, 2) = cmtex(1, 2) + c2lomt4/mt**4
           cmtex(1, 3) = cmtex(1, 3) + c3lomt4/mt**4
           cmtex(1, 4) = cmtex(1, 4) + c4lomt4/mt**4
        endif

        cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

        cmtex(:,:) = cmtex(:,:) * cex_prefactor / 16._dp

        cmtex(2,:) = cmtex(2,:) * as/2._dp/pi


      end subroutine

      function c1(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'constants.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          include 'first.f'

          complex(dp) :: c1

          complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          real(dp) :: cex_prefactor
          complex(dp) :: c1smh
          real(dp) :: c1effh

          !complex(dp) :: cmtex(2,4)

          real(dp) :: s
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

          cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

          if (first .eqv. .true.) then
            call sushi_bernini(18)
          else
            first =  .false. 
          endif

          if (approx) then
            c1 = c1effh(s(i1,i2),s(i1,i3),s(i2,i3)) * cex_prefactor / 16._dp
          else
            c1 = c1smh(s(i1,i2),s(i1,i3),s(i2,i3),mt**2) * cex_prefactor * mt**2 / 16._dp
          endif

          !call c_2l_mtex(za,zb,i1,i2,i3,cmtex)
          !write (*,*) "comparison", cmtex(1,1) / c1

      end function

      function c2(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'constants.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          include 'first.f'

          complex(dp) :: c2

          complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          real(dp) :: cex_prefactor
          complex(dp) :: c2smh
          real(dp) :: c2effh

          complex(dp) :: cmtex(2,4)

          real(dp) :: s
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

          cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

          if (first .eqv. .true.) then
            call sushi_bernini(18)
          else
            first =  .false. 
          endif

          if (approx) then
            c2 = c2effh(s(i1,i2),s(i1,i3),s(i2,i3)) * cex_prefactor / 16._dp
          else
            c2 = c2smh(s(i1,i2),s(i1,i3),s(i2,i3),mt**2) * cex_prefactor * mt**2 / 16._dp
          endif

          !call c_2l_mtex(za,zb,i1,i2,i3,cmtex)
          !write (*,*) "comparison", cmtex(1,2) / c2

      end function

      subroutine a_2l_mtex(za,zb,i1,i2,i3,amtex)
      !! one-loop and two-loop qqg amplitude coefficient in 1/mt expansion
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'constants.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          include 'scale.f'
          include 'asymptotic.f'
          include 'nf.f'

          complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp) :: amtex(2)

          real(dp) :: aex_prefactor

          real(dp) :: s
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

          real(dp) :: sman,tman,uman
          real(dp) :: mH, kgluon, kquark, pcut, nl
          real(dp) :: alomt0, alomt2, alomt4, alomt6
          complex(dp) :: amt0,amt2,amt4,amt6

          complex(dp) ::  LogMums, LogMumt, LogMumu, LogMumH
          real(dp) :: LogMuMtop

          complex(dp) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
          complex(dp) :: t11,t12,t13,t14,t15,t16,t17,t18
          complex(dp) :: t19,t20,t21,t22,t23,t24,t25,t26
          complex(dp) :: t27,t28,t29,t30,t31,t32,t33,t34
          complex(dp) :: t35,t36,t37,t38,t39,t40,t41,t42
          complex(dp) :: t43,t44,t45,t46,t47,t48,t49,t50
          complex(dp) :: t51,t52,t53,t54,t55,t56,t57,t58
          complex(dp) :: t59,t60,t61,t62,t63,t64,t65,t66
          complex(dp) :: t67,t68,t69

          complex(dp) :: cli2, cdlogwrap

          sman = s(i1,i2)
          tman = s(i1,i3)
          uman = s(i2,i3)

          mH = hmass

          nl = 5d0 ! nf
          pcut = 1d0 
          kgluon = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
     &            - ca*dlog(pcut)**2
     &            + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))
          kquark = cf*(7d0/2d0 - pisqo6)
     &            - cf*dlog(pcut)**2
     &            + 3d0/2d0*cf*(pcut - 1 - dlog(pcut))


          LogMums = cdlogwrap(-musq/sman)
          LogMumt = cdlogwrap(-musq/tman)
          LogMumu = cdlogwrap(-musq/uman)
          LogMumH = cdlogwrap(-musq/mH**2)

          LogMuMtop = Log(musq/mt**2)

          aex_prefactor = sqrt(gsq)**3 / 32d0 / pi**2 / sqrt(vevsq)

        alomt0 = 16/(3.*sman)
        alomt2 = 0.8 + (14*tman)/(45.*sman) + (14*uman)/(45.*sman)
        alomt4 = 
     -   (16*sman)/105. + (4*tman)/35. + (2*tman**2)/(63.*sman) + 
     -   (4*uman)/35. + (4*tman*uman)/(63.*sman) + 
     -   (2*uman**2)/(63.*sman)
        alomt6 = 
     -   (2*sman**2)/63. + (11*sman*tman)/315. + (2*tman**2)/105. + 
     -   (13*tman**3)/(3150.*sman) + (11*sman*uman)/315. + 
     -   (4*tman*uman)/105. + (13*tman**2*uman)/(1050.*sman) + 
     -   (2*uman**2)/105. + (13*tman*uman**2)/(1050.*sman) + 
     -   (13*uman**3)/(3150.*sman)

        amtex(:) = 0._dp

        if (mtex >= 0) then
            amtex(1) = amtex(1) + alomt0
            include 'qag/a_mt0.f'
            amtex(2) = amtex(2) + amt0
        endif

        if (mtex >= 2) then
            amtex(1) = amtex(1) + alomt2/mt**2
            include 'qag/a_mt2.f'
            amtex(2) = amtex(2) + amt2/mt**2
        endif

        if (mtex >= 4) then
            amtex(1) = amtex(1) + alomt4/mt**4
            include 'qag/a_mt4.f'
            amtex(2) = amtex(2) + amt4/mt**4
        endif

        if (mtex >= 6) then
            amtex(1) = amtex(1) + alomt6/mt**6
            include 'qag/a_mt6.f'
            amtex(2) = amtex(2) + amt6/mt**6
        endif


        ! our normalization
        ! in HTL the 1-loop coefficient should be 1
        amtex(:) = amtex(:) / alomt0

        ! compatibility with other amps
        amtex(:) = amtex(:) * aex_prefactor / s(i1,i2) * 8._dp/3._dp

        amtex(2) = amtex(2) * as/2._dp/pi


      end subroutine

      function aex(za,zb,i1,i2,i3)
      !! one-loop qqg amplitude coefficient for exact top mass dependence
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'
          include 'constants.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          include 'first.f'

          complex(dp) :: aex

          complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          real(dp) :: aex_prefactor
          complex(dp) :: asmh

          complex(dp) :: amtex(2)

          real(dp) :: s
          s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

          aex_prefactor = sqrt(gsq)**3 / 32d0 / pi**2 / sqrt(vevsq)

          if (first .eqv. .true.) then
            call sushi_bernini(18)
          else
            first =  .false. 
          endif

          if (approx) then
            aex = 1._dp/s(i1,i2) * aex_prefactor * 8._dp/3._dp
          else
            aex = asmh(s(i1,i2),s(i2,i3)+s(i1,i3),mt**2) / s(i1,i2)
     &              * aex_prefactor * mt**2 * 8._dp/3._dp
          endif

      end function

      subroutine hjetmass_qqg_mpp_2l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(2)

          complex(dp) :: amtex(2)

          call a_2l_mtex(za,zb,i1,i2,i3,amtex)

          amps(:) = za(i2,i1)*zb(i3,i1)**2/sqrt(2._dp) * amtex(:)

      end subroutine

      subroutine hjetmass_qqg_mpm_2l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(2)

          complex(dp) :: amtex(2)

          call a_2l_mtex(za,zb,i1,i2,i3,amtex)

          amps(:) = zb(i2,i1)*za(i2,i3)**2/sqrt(2._dp) * amtex(:)

      end subroutine

      function hjetmass_qqg_mpp(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp) :: hjetmass_qqg_mpp

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          hjetmass_qqg_mpp = za(i2,i1)*zb(i3,i1)**2/sqrt(2._dp) * 
     &          aex(za,zb,i1,i2,i3)

      end function

      function hjetmass_qqg_mpm(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp) :: hjetmass_qqg_mpm

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          hjetmass_qqg_mpm = zb(i2,i1)*za(i2,i3)**2/sqrt(2._dp) * 
     &          aex(za,zb,i1,i2,i3)

      end function



      function hjetmass_ggg_ppp(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp) :: hjetmass_ggg_ppp

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          hjetmass_ggg_ppp = 2*sqrt(2._dp) * c1(za,zb,i1,i2,i3) 
     &              / za(i1,i2) / za(i1,i3) / za(i2,i3)

      end function

      subroutine hjetmass_ggg_ppp_2l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(2)

          complex(dp) :: cmtex(2,4) 

          call c_2l_mtex(za,zb,i1,i2,i3,cmtex)

          amps(:) = 2*sqrt(2._dp) * cmtex(:,1) / za(i1,i2) / za(i1,i3) / za(i2,i3)

      end subroutine

      subroutine hjetmass_ggg_ppm_2l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(2)

          complex(dp) :: cmtex(2,4) 

          call c_2l_mtex(za,zb,i1,i2,i3,cmtex)

          amps(:) = 2*sqrt(2._dp) * cmtex(:,2) * zb(i2,i1)
     &          / zb(i3,i1) / zb(i3,i2) / za(i1,i2)**2

      end subroutine

      subroutine hjetmass_ggg_pmp_2l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(2)

          complex(dp) :: cmtex(2,4) 

          call c_2l_mtex(za,zb,i1,i2,i3,cmtex)

          amps(:) = 2*sqrt(2._dp) * cmtex(:,3) * zb(i3,i1)
     &          / zb(i2,i1) / zb(i2,i3) / za(i1,i3)**2

      end subroutine

      function hjetmass_ggg_ppm(za,zb,i1,i2,i3)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp) :: hjetmass_ggg_ppm

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          hjetmass_ggg_ppm = 2*sqrt(2._dp) * c2(za,zb,i1,i2,i3) * zb(i2,i1)
     &          / zb(i3,i1) / zb(i3,i2) / za(i1,i2)**2

      end function

      ! just ppm with i2,i3 switched
      function hjetmass_ggg_pmp(za,zb,i1,i3,i2)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp) :: hjetmass_ggg_pmp

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3

          hjetmass_ggg_pmp = 2*sqrt(2._dp) * c2(za,zb,i1,i2,i3) * zb(i2,i1)
     &          / zb(i3,i1) / zb(i3,i2) / za(i1,i2)**2

      end function

      subroutine hjetmass_ggg_ppp_1l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(3)

          real(dp) :: cmtex(3,4) 

          call c_1l_mtex(za,zb,i1,i2,i3,cmtex)

          amps(:) = 2*sqrt(2._dp) * cmtex(:,1) / za(i1,i2) / za(i1,i3) / za(i2,i3)

      end subroutine

      subroutine hjetmass_ggg_ppm_1l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(3)

          real(dp) :: cmtex(3,4) 

          call c_1l_mtex(za,zb,i1,i2,i3,cmtex)

          amps(:) = 2*sqrt(2._dp) * cmtex(:,2) * zb(i2,i1)
     &          / zb(i3,i1) / zb(i3,i2) / za(i1,i2)**2

      end subroutine

      subroutine hjetmass_ggg_pmp_1l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(3)

          real(dp) :: cmtex(3,4) 

          call c_1l_mtex(za,zb,i1,i2,i3,cmtex)

          amps(:) = 2*sqrt(2._dp) * cmtex(:,3) * zb(i3,i1)
     &          / zb(i2,i1) / zb(i2,i3) / za(i1,i3)**2

      end subroutine

      subroutine hjetmass_qqg_mpp_1l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(3)

          real(dp) :: amtex(3)

          call a_1l_mtex(za,zb,i1,i2,i3,amtex)

          amps(:) = za(i2,i1)*zb(i3,i1)**2/sqrt(2._dp) * amtex(:)
      end subroutine

      subroutine hjetmass_qqg_mpm_1l_mtex(za,zb,i1,i2,i3,amps)
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'masses.f'

          complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
          integer, intent(in) :: i1,i2,i3
          complex(dp), intent(out) :: amps(3)

          real(dp) :: amtex(3)

          call a_1l_mtex(za,zb,i1,i2,i3,amtex)

          amps(:) = zb(i2,i1)*za(i2,i3)**2/sqrt(2._dp) * amtex(:)

      end subroutine


      end module
