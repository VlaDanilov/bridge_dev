      module singletop2_decaywidth_m
        use types
        implicit none

        public :: singletop2_decaywidth

        private

        contains

        ! NLO width t -> Wb with mb=0 and EFT operators
        function singletop2_decaywidth()
            use anomcoup_tbW
            implicit none
            include 'masses.f'
            include 'ewcouple.f'
            include 'qcdcouple.f'
            include 'constants.f'
            include 'kpart.f'

            real(dp) :: singletop2_decaywidth

            real(dp) :: r
            real(dp) :: bornSM, nlofac
            real(dp) :: lambda2_lo, lambda4_lo
            real(dp) :: nloSM, lambda2_nlo, lambda4_nlo
            real(dp) :: mtsq

            real(dp) :: lotopdecaywidth, nloratiotopdecay 
            real(dp) :: rescale

            real(dp) :: ddilog, topwidth

            real(dp) :: musq

            musq = mt**2

            ! Vtb = 1

            r = wmass/mt
            mtsq = mt**2

            bornSM = Gf*mt**3/8._dp/pi/sqrt(2._dp) * (1-r**2)**2 * (1 + 2*r**2)
            nlofac = CF*as/4._dp/pi

            lambda2_lo = 0._dp
            lambda4_lo = 0._dp

            lambda2_nlo = 0._dp
            lambda4_nlo = 0._dp

            ! all contributions are normalized to the SM LO decay width

            nloSM = (-4*Pi**2)/3._dp + (5 + 9*r**2 - 6*r**4)/(1 + r**2 - 2*r**4) - 8*ddilog(r**2) +
     -  Log(r**2)*((4*r**2*(-1 + r**2 + 2*r**4))/(1 - 3*r**4 + 2*r**6) - 4*Log(1 - r**2)) +
     -  (-4 - 6/(1 + 2*r**2))*Log(1 - r**2)

            if (enable_lambda2) then
            lambda2_lo = 2*real(anomc1) + (24*mt*r**2*real(anomc3))/(1 + 2*r**2)

            ! compared with 1404.1264
            lambda2_nlo = (-2*(24*(1 - 3*r**4 + 2*r**6)*ddilog(r**2) + 
     -       (-1 + r**2)*(3*(5 + 9*r**2 - 6*r**4) + 
     -          Pi**2*(-4 - 4*r**2 + 8*r**4) + 
     -          6*(-5 + r**2 + 4*r**4)*Log(1 - r**2)) + 
     -       12*Log(r**2)*(r**2 - r**4 - 2*r**6 + 
     -          (1 - 3*r**4 + 2*r**6)*Log(1 - r**2)))*real(anomc1) - 
     -    3*mt*(8*(-17*r**2 + 4*Pi**2*r**2 + 38*r**4 - 8*Pi**2*r**4 - 21*r**6 + 
     -          4*Pi**2*r**6 + 24*r**2*(-1 + r**2)**2*ddilog(r**2) + 
     -          3*r**2*(-1 + r**2)**2*Log(Musq/Mtsq) + 12*r**4*Log(r**2) - 
     -          8*r**6*Log(r**2) + 4*Log(1 - r**2) + 6*r**2*Log(1 - r**2) - 
     -          24*r**4*Log(1 - r**2) + 14*r**6*Log(1 - r**2) + 
     -          12*r**2*Log(r**2)*Log(1 - r**2) - 
     -          24*r**4*Log(r**2)*Log(1 - r**2) + 
     -          12*r**6*Log(r**2)*Log(1 - r**2))*real(anomc3) + 
     -       (24*r**2*(-1 + r**2)**2*Log(Musq/Mtsq) - 
     -          8*r**4*(3 + r**2)*Log(r**2) - 
     -          (-1 + r**2)*(3 + 43*r**2 - 78*r**4 + 
     -             16*(-1 + r**2)**2*Log(1 - r**2)))*real(anomc6)))/
     -  (3 - 9*r**4 + 6*r**6)
            endif

         if (enable_lambda4) then

         lambda4_lo = aimag(anomc2)**2 + (16*mt**2*r**2*(2 + r**2)*aimag(anomc3)**2)/
     -   (1 + 2*r**2) + (24*mt*r**2*aimag(anomc2)*aimag(anomc4))/(1 + 2*r**2) + 
     -  (16*mt**2*r**2*(2 + r**2)*aimag(anomc4)**2)/(1 + 2*r**2) + 
     -  real(anomc1)**2 + real(anomc2)**2 + 
     -  (24*mt*r**2*real(anomc1)*real(anomc3))/(1 + 2*r**2) + 
     -  (16*mt**2*r**2*(2 + r**2)*real(anomc3)**2)/(1 + 2*r**2) + 
     -  (24*mt*r**2*real(anomc2)*real(anomc4))/(1 + 2*r**2) + 
     -  (16*mt**2*r**2*(2 + r**2)*real(anomc4)**2)/(1 + 2*r**2)

         lambda4_nlo = ((4*Pi**2*(1 + r**2 - 2*r**4) + 3*(-5 - 9*r**2 + 6*r**4))*
     -     aimag(anomc2)**2)/(-3 - 3*r**2 + 6*r**4) - 
     -  8*ddilog(r**2)*aimag(anomc2)**2 - 
     -  (16*mt**2*r**2*(32 - 13*r**2 - 7*r**4 + 4*Pi**2*(-2 + r**2 + r**4))*
     -     aimag(anomc3)**2)/(-3 - 3*r**2 + 6*r**4) - 
     -  (128*mt**2*r**2*(2 + r**2)*ddilog(r**2)*aimag(anomc3)**2)/(1 + 2*r**2) - 
     -  (8*mt*r**2*(17 - 21*r**2 + 4*Pi**2*(-1 + r**2))*aimag(anomc2)*
     -     aimag(anomc4))/((-1 + r**2)*(1 + 2*r**2)) - 
     -  (192*mt*r**2*ddilog(r**2)*aimag(anomc2)*aimag(anomc4))/(1 + 2*r**2) - 
     -  (16*mt**2*r**2*(32 - 13*r**2 - 7*r**4 + 4*Pi**2*(-2 + r**2 + r**4))*
     -     aimag(anomc4)**2)/(-3 - 3*r**2 + 6*r**4) - 
     -  (128*mt**2*r**2*(2 + r**2)*ddilog(r**2)*aimag(anomc4)**2)/(1 + 2*r**2) + 
     -  (16*mt**2*r**2*(8 - 31*r**2 + 11*r**4)*aimag(anomc3)*aimag(anomc6))/
     -   (3 + 3*r**2 - 6*r**4) + 
     -  (mt**2*(1 + 12*r**4 - 16*r**6 + 3*r**8)*aimag(anomc6)**2)/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (4*mt*(1 + 7*r**2*(-2 + r**2))*aimag(anomc2)*aimag(anomc7))/
     -   ((-1 + r**2)*(1 + 2*r**2)) + 
     -  (16*mt**2*r**2*(-19 + 14*r**2 + 11*r**4)*aimag(anomc4)*aimag(anomc7))/
     -   (3 + 3*r**2 - 6*r**4) + 
     -  (mt**2*(1 + 12*r**4 - 16*r**6 + 3*r**8)*aimag(anomc7)**2)/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (32*mt**2*r**2*(2 + r**2)*aimag(anomc3)**2*Log(Musq/Mtsq))/
     -   (1 + 2*r**2) - (24*mt*r**2*aimag(anomc2)*aimag(anomc4)*Log(Musq/Mtsq))/
     -   (1 + 2*r**2) - (32*mt**2*r**2*(2 + r**2)*aimag(anomc4)**2*
     -     Log(Musq/Mtsq))/(1 + 2*r**2) - 
     -  (32*mt**2*r**2*(2 + r**2)*aimag(anomc3)*aimag(anomc6)*Log(Musq/Mtsq))/
     -   (1 + 2*r**2) - (24*mt*r**2*aimag(anomc2)*aimag(anomc7)*Log(Musq/Mtsq))/
     -   (1 + 2*r**2) - (32*mt**2*r**2*(2 + r**2)*aimag(anomc4)*aimag(anomc7)*
     -     Log(Musq/Mtsq))/(1 + 2*r**2) - 
     -  (128*mt**2*r**4*aimag(anomc3)**2*Log(r**(-2)))/(1 - 3*r**4 + 2*r**6) - 
     -  (128*mt**2*r**4*aimag(anomc4)**2*Log(r**(-2)))/(1 - 3*r**4 + 2*r**6) + 
     -  (64*mt**2*r**4*(1 + r**2)*aimag(anomc3)*aimag(anomc6)*Log(r**(-2)))/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (12*mt**2*r**4*aimag(anomc6)**2*Log(r**(-2)))/(1 - 3*r**4 + 2*r**6) - 
     -  (64*mt**2*r**4*(1 + r**2)*aimag(anomc4)*aimag(anomc7)*Log(r**(-2)))/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (12*mt**2*r**4*aimag(anomc7)**2*Log(r**(-2)))/(1 - 3*r**4 + 2*r**6) + 
     -  (64*mt**2*r**4*aimag(anomc3)*aimag(anomc6)*Log(r**2))/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (24*mt*r**4*aimag(anomc2)*aimag(anomc7)*Log(r**2))/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (96*mt**2*r**4*(2 + r**2)*aimag(anomc4)*aimag(anomc7)*Log(r**2))/
     -   (3 - 9*r**4 + 6*r**6) + 
     -  4*aimag(anomc2)**2*Log(r**2)*
     -   ((r**2*(-1 + r**2 + 2*r**4))/(1 - 3*r**4 + 2*r**6) - Log(1 - r**2)) - 
     -  (2*(5 + 4*r**2)*aimag(anomc2)**2*Log(1 - r**2))/(1 + 2*r**2) - 
     -  (32*mt**2*r**2*(8 + r**2)*aimag(anomc3)**2*Log(1 - r**2))/(1 + 2*r**2) - 
     -  (16*mt*(2 + 7*r**2)*aimag(anomc2)*aimag(anomc4)*Log(1 - r**2))/
     -   (1 + 2*r**2) - (32*mt**2*r**2*(8 + r**2)*aimag(anomc4)**2*
     -     Log(1 - r**2))/(1 + 2*r**2) + 
     -  (32*mt**2*(-1 + r**2)**2*aimag(anomc3)*aimag(anomc6)*Log(1 - r**2))/
     -   (1 + 2*r**2) + (24*mt*(-1 + r**2)*aimag(anomc2)*aimag(anomc7)*
     -     Log(1 - r**2))/(1 + 2*r**2) + 
     -  (32*mt**2*(-2 + r**2 + r**4)*aimag(anomc4)*aimag(anomc7)*Log(1 - r**2))/
     -   (1 + 2*r**2) - (32*mt*aimag(anomc2)*aimag(anomc4)*Log(r**2)*
     -     (3*r**4 - 2*r**6 + 3*r**2*(-1 + r**2)**2*Log(1 - r**2)))/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (64*mt**2*aimag(anomc3)**2*Log(r**2)*
     -     (-(r**4*(-4 + 2*r**2 + r**4)) + 
     -       (2*r**2 - 3*r**4 + r**8)*Log(1 - r**2)))/(1 - 3*r**4 + 2*r**6) - 
     -  (64*mt**2*aimag(anomc4)**2*Log(r**2)*
     -     (-(r**4*(-4 + 2*r**2 + r**4)) + 
     -       (2*r**2 - 3*r**4 + r**8)*Log(1 - r**2)))/(1 - 3*r**4 + 2*r**6) + 
     -  ((4*Pi**2*(1 + r**2 - 2*r**4) + 3*(-5 - 9*r**2 + 6*r**4))*
     -     real(anomc1)**2)/(-3 - 3*r**2 + 6*r**4) - 
     -  8*ddilog(r**2)*real(anomc1)**2 + 
     -  4*Log(r**2)*((r**2*(-1 + r**2 + 2*r**4))/(1 - 3*r**4 + 2*r**6) - 
     -     Log(1 - r**2))*real(anomc1)**2 - 
     -  (2*(5 + 4*r**2)*Log(1 - r**2)*real(anomc1)**2)/(1 + 2*r**2) + 
     -  ((4*Pi**2*(1 + r**2 - 2*r**4) + 3*(-5 - 9*r**2 + 6*r**4))*
     -     real(anomc2)**2)/(-3 - 3*r**2 + 6*r**4) - 
     -  8*ddilog(r**2)*real(anomc2)**2 + 
     -  4*Log(r**2)*((r**2*(-1 + r**2 + 2*r**4))/(1 - 3*r**4 + 2*r**6) - 
     -     Log(1 - r**2))*real(anomc2)**2 - 
     -  (2*(5 + 4*r**2)*Log(1 - r**2)*real(anomc2)**2)/(1 + 2*r**2) - 
     -  (8*mt*r**2*(17 - 21*r**2 + 4*Pi**2*(-1 + r**2))*real(anomc1)*
     -     real(anomc3))/((-1 + r**2)*(1 + 2*r**2)) - 
     -  (192*mt*r**2*ddilog(r**2)*real(anomc1)*real(anomc3))/(1 + 2*r**2) - 
     -  (24*mt*r**2*Log(Musq/Mtsq)*real(anomc1)*real(anomc3))/(1 + 2*r**2) - 
     -  (16*mt*(2 + 7*r**2)*Log(1 - r**2)*real(anomc1)*real(anomc3))/
     -   (1 + 2*r**2) - (32*mt*Log(r**2)*
     -     (3*r**4 - 2*r**6 + 3*r**2*(-1 + r**2)**2*Log(1 - r**2))*real(anomc1)*
     -     real(anomc3))/(1 - 3*r**4 + 2*r**6) - 
     -  (16*mt**2*r**2*(32 - 13*r**2 - 7*r**4 + 4*Pi**2*(-2 + r**2 + r**4))*
     -     real(anomc3)**2)/(-3 - 3*r**2 + 6*r**4) - 
     -  (128*mt**2*r**2*(2 + r**2)*ddilog(r**2)*real(anomc3)**2)/(1 + 2*r**2) - 
     -  (32*mt**2*r**2*(2 + r**2)*Log(Musq/Mtsq)*real(anomc3)**2)/
     -   (1 + 2*r**2) - (128*mt**2*r**4*Log(r**(-2))*real(anomc3)**2)/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (32*mt**2*r**2*(8 + r**2)*Log(1 - r**2)*real(anomc3)**2)/(1 + 2*r**2) - 
     -  (64*mt**2*Log(r**2)*(-(r**4*(-4 + 2*r**2 + r**4)) + 
     -       (2*r**2 - 3*r**4 + r**8)*Log(1 - r**2))*real(anomc3)**2)/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (8*mt*r**2*(17 - 21*r**2 + 4*Pi**2*(-1 + r**2))*real(anomc2)*
     -     real(anomc4))/((-1 + r**2)*(1 + 2*r**2)) - 
     -  (192*mt*r**2*ddilog(r**2)*real(anomc2)*real(anomc4))/(1 + 2*r**2) - 
     -  (24*mt*r**2*Log(Musq/Mtsq)*real(anomc2)*real(anomc4))/(1 + 2*r**2) - 
     -  (16*mt*(2 + 7*r**2)*Log(1 - r**2)*real(anomc2)*real(anomc4))/
     -   (1 + 2*r**2) - (32*mt*Log(r**2)*
     -     (3*r**4 - 2*r**6 + 3*r**2*(-1 + r**2)**2*Log(1 - r**2))*real(anomc2)*
     -     real(anomc4))/(1 - 3*r**4 + 2*r**6) - 
     -  (16*mt**2*r**2*(32 - 13*r**2 - 7*r**4 + 4*Pi**2*(-2 + r**2 + r**4))*
     -     real(anomc4)**2)/(-3 - 3*r**2 + 6*r**4) - 
     -  (128*mt**2*r**2*(2 + r**2)*ddilog(r**2)*real(anomc4)**2)/(1 + 2*r**2) - 
     -  (32*mt**2*r**2*(2 + r**2)*Log(Musq/Mtsq)*real(anomc4)**2)/
     -   (1 + 2*r**2) - (128*mt**2*r**4*Log(r**(-2))*real(anomc4)**2)/
     -   (1 - 3*r**4 + 2*r**6) - 
     -  (32*mt**2*r**2*(8 + r**2)*Log(1 - r**2)*real(anomc4)**2)/(1 + 2*r**2) - 
     -  (64*mt**2*Log(r**2)*(-(r**4*(-4 + 2*r**2 + r**4)) + 
     -       (2*r**2 - 3*r**4 + r**8)*Log(1 - r**2))*real(anomc4)**2)/
     -   (1 - 3*r**4 + 2*r**6) + 
     -  (mt*(3 + 43*r**2 - 78*r**4)*real(anomc1)*real(anomc6))/
     -   ((-1 + r**2)*(1 + 2*r**2)) - 
     -  (24*mt*r**2*Log(Musq/Mtsq)*real(anomc1)*real(anomc6))/(1 + 2*r**2) + 
     -  (8*mt*r**4*(3 + r**2)*Log(r**2)*real(anomc1)*real(anomc6))/
     -   (1 - 3*r**4 + 2*r**6) + 
     -  (16*mt*(-1 + r**2)*Log(1 - r**2)*real(anomc1)*real(anomc6))/
     -   (1 + 2*r**2) + (4*mt**2*r**2*(-115 + 47*r**2 + 44*r**4)*real(anomc3)*
     -     real(anomc6))/(3 + 3*r**2 - 6*r**4) - 
     -  (32*mt**2*r**2*(2 + r**2)*Log(Musq/Mtsq)*real(anomc3)*real(anomc6))/
     -   (1 + 2*r**2) + (64*mt**2*r**6*Log(r**(-2))*real(anomc3)*real(anomc6))/
     -   (1 - 3*r**4 + 2*r**6) + 
     -  (32*mt**2*r**6*Log(r**2)*real(anomc3)*real(anomc6))/
     -   (1 - 3*r**4 + 2*r**6) + 
     -  (32*mt**2*(-1 + r**4)*Log(1 - r**2)*real(anomc3)*real(anomc6))/
     -   (1 + 2*r**2) + (mt**2*(1 + 12*r**4 - 16*r**6 + 3*r**8)*
     -     real(anomc6)**2)/(1 - 3*r**4 + 2*r**6) - 
     -  (12*mt**2*r**4*Log(r**(-2))*real(anomc6)**2)/(1 - 3*r**4 + 2*r**6) - 
     -  (4*mt*(1 + 7*r**2*(-2 + r**2))*real(anomc2)*real(anomc7))/
     -   ((-1 + r**2)*(1 + 2*r**2)) - 
     -  (24*mt*r**2*Log(Musq/Mtsq)*real(anomc2)*real(anomc7))/(1 + 2*r**2) - 
     -  (24*mt*r**4*Log(r**2)*real(anomc2)*real(anomc7))/
     -   (1 - 3*r**4 + 2*r**6) + 
     -  (24*mt*(-1 + r**2)*Log(1 - r**2)*real(anomc2)*real(anomc7))/
     -   (1 + 2*r**2) + (16*mt**2*r**2*(-19 + 14*r**2 + 11*r**4)*real(anomc4)*
     -     real(anomc7))/(3 + 3*r**2 - 6*r**4) - 
     -  (32*mt**2*r**2*(2 + r**2)*Log(Musq/Mtsq)*real(anomc4)*real(anomc7))/
     -   (1 + 2*r**2) - (64*mt**2*r**4*(1 + r**2)*Log(r**(-2))*real(anomc4)*
     -     real(anomc7))/(1 - 3*r**4 + 2*r**6) - 
     -  (96*mt**2*r**4*(2 + r**2)*Log(r**2)*real(anomc4)*real(anomc7))/
     -   (3 - 9*r**4 + 6*r**6) + 
     -  (32*mt**2*(-2 + r**2 + r**4)*Log(1 - r**2)*real(anomc4)*real(anomc7))/
     -   (1 + 2*r**2) + (mt**2*(1 + 12*r**4 - 16*r**6 + 3*r**8)*
     -     real(anomc7)**2)/(1 - 3*r**4 + 2*r**6) - 
     -  (12*mt**2*r**4*Log(r**(-2))*real(anomc7)**2)/(1 - 3*r**4 + 2*r**6)
     
         endif

c           ! test of SM width expression
c           write (*,*) "WIDTH TEST", (1+nlofac*nloSM)/topwidth(mt,wmass)
c           pause

c           write (*,*) "WIDTH EFT CONTRIBUTION: ",  lambda2_lo/lambda**2 + lambda4_lo/lambda**4
c    &              + nlofac*lambda2_nlo/lambda**2 + nlofac*lambda4_nlo/lambda**4

            ! set this to 1 to disable mb,wwidth effects
            rescale =  lotopdecaywidth(mt,mb,wmass,wwidth)/bornSM

            if (origKpart == klord) then
                singletop2_decaywidth = bornSM * (rescale + lambda2_lo/lambda**2 + lambda4_lo/lambda**4)
            else
                singletop2_decaywidth = bornSM * (rescale + lambda2_lo/lambda**2 + lambda4_lo/lambda**4
     &              + nloratiotopdecay(mt,mb,wmass,wwidth)*rescale + nlofac*lambda2_nlo/lambda**2 + nlofac*lambda4_nlo/lambda**4
     &          )
            endif



        end function

      end module
