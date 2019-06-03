      double complex function C3SMh(s,t,u,m2)
        implicit none
        double precision s,t,u,m2
        double complex C2SMh

        C3SMh = C2SMh(t,s,u,m2)

      end function C3SMh

      double complex function C4SMh(s,t,u,m2)
        implicit none
        double precision s,t,u,m2
        double complex C2SMh

        C4SMh = C2SMh(u,s,t,m2)

      end function C4SMh
      
      double precision function virtme_gg_multi(sman,tman,uman,mtex)
          implicit none
          include 'types.f'
          double precision, intent(in) :: sman, tman, uman
          integer, intent(in) :: mtex
          double complex cdlogwrap
          double complex dilogc
          double complex c1lomt0,c1lomt2,c1lomt4
          double complex c2lomt0,c2lomt2,c2lomt4
          double complex c3lomt0,c3lomt2,c3lomt4
          double complex c4lomt0,c4lomt2,c4lomt4
          double complex c1mt0,c1mt2,c1mt4
          double complex c2mt0,c2mt2,c2mt4
          double complex c3mt0,c3mt2,c3mt4
          double complex c4mt0,c4mt2,c4mt4
          include 'constants.f'
          include 'masses.f'
          include 'scale.f'
          include 'epinv.f'
          include 'qcdcouple.f'
          include 'ewcouple.f'
          double complex t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
          double complex t11,t12,t13,t14,t15,t16,t17,t18
          double complex t19,t20,t21,t22,t23,t24,t25,t26
          double complex t27,t28,t29,t30,t31,t32,t33,t34
          double complex t35,t36,t37,t38,t39,t40
          double complex square, square2

          double complex C1SMh, C2SMh, C3SMh, C4SMh
          double complex c1nlo,c2nlo,c3nlo,c4nlo

          double precision pcut, nl, kg
          double complex LogMums, LogMumt, LogMumu
          double precision LogMuMtop, mH

          mH = hmass

          LogMums = cdlogwrap(-musq/sman)
          LogMumt = cdlogwrap(-musq/tman)
          LogMumu = cdlogwrap(-musq/uman)
          LogMuMtop = log(musq/mt**2)

        pcut = 1d0
        nl =  5d0
        kg = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
     &          - ca*dlog(pcut)**2
     &          + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))
          

       c1lomt0 = (0d0,0d0)
       c1lomt2 = (0d0,0d0)
       c1lomt4 = (0d0,0d0)
       c2lomt0 = (0d0,0d0)
       c2lomt2 = (0d0,0d0)
       c2lomt2 = (0d0,0d0)
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

       c1nlo = 0d0
       c2nlo = 0d0
       c3nlo = 0d0
       c4nlo = 0d0

       if (mtex >= 0) then
           include 'ggg/c1_mt0.f'
           include 'ggg/c2_mt0.f'
           include 'ggg/c3_mt0.f'
           include 'ggg/c4_mt0.f'

           c1nlo = c1mt0
           c2nlo = c2mt0
           c3nlo = c3mt0
           c4nlo = c4mt0

       end if

       if (mtex >= 2) then
           include 'ggg/c1_mt2.f'
           include 'ggg/c2_mt2.f'
           include 'ggg/c3_mt2.f'
           include 'ggg/c4_mt2.f'

           c1nlo = c1nlo + c1mt2/mt**2
           c2nlo = c2nlo + c2mt2/mt**2
           c3nlo = c3nlo + c3mt2/mt**2
           c4nlo = c4nlo + c4mt2/mt**2
       end if

       if (mtex >= 4) then
           include 'ggg/c1_mt4.f'
           include 'ggg/c2_mt4.f'
           include 'ggg/c3_mt4.f'
           include 'ggg/c4_mt4.f'

           c1nlo = c1nlo + c1mt4/mt**4
           c2nlo = c2nlo + c2mt4/mt**4
           c3nlo = c3nlo + c3mt4/mt**4
           c4nlo = c4nlo + c4mt4/mt**4
       end if

       square = (0d0,0d0)
       
       square =  square
     &             + conjg(C1SMh(sman,tman,uman,mt**2) * mt**2)*c1nlo
     &             + conjg(C2SMh(sman,tman,uman,mt**2) * mt**2)*c2nlo
     &             + conjg(C3SMh(sman,tman,uman,mt**2) * mt**2)*c3nlo
     &             + conjg(C4SMh(sman,tman,uman,mt**2) * mt**2)*c4nlo

       square = square*(9d0/16d0)/sman/tman/uman

      virtme_gg_multi = 2*dreal(square) * as/2d0/pi
     &         * (as/3d0/pi)**2 / vevsq * 4d0*pi*as*24d0/256d0

      end function virtme_gg_multi
      
