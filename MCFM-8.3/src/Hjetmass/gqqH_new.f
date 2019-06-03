      double complex function cdlogwrap(x)
        implicit none
        double precision x
        cdlogwrap = cdlog(dcmplx(x,0d0))
      end function cdlogwrap

      double complex function dilogc(x)
        implicit none
        double precision x
        double complex cli2
        
        dilogc = cli2(dcmplx(x,0d0))

      end function dilogc
      
      double precision function virtme_gq_new(sman, tman, uman, mtex)
        implicit none
        include 'types.f'
        double precision :: sman, tman, uman
        integer, intent(in) :: mtex
        include 'constants.f'
        include 'masses.f'
        include 'scale.f'
        include 'epinv.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'

        double precision dilog
        double complex cdlogwrap, dilogc
        double complex cli2

        double precision kgluon, kquark, nl, pcut, mH

        double complex LogMums, LogMumt, LogMumu, LogMumH
        double precision LogMuMtop

        double complex finite

        double complex ASMh

        double complex alomt0, alomt2, alomt4, alomt6
        double complex amt0, amt2, amt4, amt6
        double complex anlo
        double complex t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
        double complex t11,t12,t13,t14,t15,t16,t17,t18
        double complex t19,t20,t21,t22,t23,t24,t25,t26
        double complex t27,t28,t29,t30,t31,t32,t33,t34
        double complex t35,t36,t37,t38,t39,t40,t41,t42
        double complex t43,t44,t45,t46,t47,t48,t49,t50
        double complex t51,t52,t53,t54,t55,t56,t57,t58
        double complex t59,t60,t61,t62,t63,t64,t65,t66
        double complex t67,t68,t69

        mH = hmass

        nl = 5d0
        pcut = 1d0 
        kgluon = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
     &          - ca*dlog(pcut)**2
     &          + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))
        kquark = cf*(7d0/2d0 - pisqo6)
     &          - cf*dlog(pcut)**2
     &          + 3d0/2d0*cf*(pcut - 1 - dlog(pcut))


        LogMums = cdlogwrap(-musq/sman)
        LogMumt = cdlogwrap(-musq/tman)
        LogMumu = cdlogwrap(-musq/uman)
        LogMumH = cdlogwrap(-musq/mH**2)

        LogMuMtop = cdlogwrap(musq/mt**2)

        finite = 0d0

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

        anlo = 0d0

       if (mtex >= 0) then
           include 'qag/a_mt0.f'
           anlo = amt0
       end if

       if (mtex >= 2) then
           include 'qag/a_mt2.f'
           anlo = anlo + amt2/mt**2
       end if

       if (mtex >= 4) then
           include 'qag/a_mt4.f'
           anlo = anlo + amt4/mt**4
       end if

       if (mtex >= 6) then
           include 'qag/a_mt6.f'
           anlo = anlo + amt6/mt**6
       end if

       finite = 0d0

       finite = conjg(ASMh(sman,tman+uman,mt**2)*mt**2*alomt0) * anlo

       finite = finite * sman*(tman**2 + uman**2) * 9d0/256d0

      virtme_gq_new = 2d0*dreal(finite) * as/2/pi *
     &             (as/3/pi)**2 / vevsq *
     &            4*pi*as*4d0/(4d0*3d0*8d0)

      end function virtme_gq_new

