
      double complex function hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_dp
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt
          double precision cg

      t1 = za(i1, i4)
      t2 = za(i2, i4)
      t3 = zb(i3, i1)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = zb(i2, i1)
      t7 = za(i1, i3)
      t8 = za(i2, i3)
      t9 = zb(i3, i2)
      t10 = t5 * t6
      t11 = t7 * t3
      t12 = t8 * t9 + t10 + t11
      t13 = zb(i4, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t1 * t4
      t17 = t2 * t13
      t18 = t14 * t15
      t19 = t16 + t17 + t18
      if ( dreal(t19) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t18 = cg * cdsqrt(t19 ** 2) + t16 + t17 + t18
      t19 = 0.1D1 / t18
      t11 = -2 * t16 * t12 * t19 + t10 + t11
      t20 = -2 * t1 * t12
      t7 = t4 * (t1 * t11 + t2 * (t20 * t13 * t19 + t7 * t9))
      t21 = -2 * t2 * t12 * t4 * t19 + t3 * t8
      t22 = -2 * t14 * t12 * t4 * t19 - t6 * t8
      t23 = t14 * t3 + t2 * t6
      t24 = t21 * t13
      t25 = t1 * (t11 * t4 + t24)
      t20 = t20 * t15 * t19 - t5 * t9
      t26 = t1 * t3 + t2 * t9
      t8 = 0.1D1 / t8
      t18 = 0.1D1 / t18
      t27 = 0.1D1 / t6
      t28 = 0.1D1 / t12
      t29 = 0.1D1 / t4
      t22 = 0.1D1 / t22
      t25 = 0.1D1 / t25
      t30 = 0.1D1 / t5
      t31 = t3 ** 2
      t32 = t18 ** 2
      t33 = t25 ** 2
      t34 = t25 * t33
      t35 = t2 ** 2
      t36 = t22 ** 2
      t37 = t21 ** 2
      t38 = t1 ** 2
      t39 = t38 ** 2
      t40 = t1 * t38
      t41 = t12 ** 2
      t42 = t4 ** 2
      t43 = mt ** 2
      t44 = t14 ** 2 * t37 * t36
      t45 = t1 * t25
      t46 = t30 * t27
      t47 = t46 * t41 * t32
      t48 = t2 * t8
      t49 = t48 * t38
      t50 = t46 * t43
      t51 = t18 * t12
      t52 = t40 * t15 * t8
      t36 = t14 * t36
      t53 = t8 * t30
      t54 = t27 * t23
      t55 = 0.1D1 / t14
      t7 = 0.1D1 / t7
      t56 = 0.1D1 / t15
      t20 = 0.1D1 / t20
      t57 = t2 * t21
      t58 = t11 * t25
      t59 = t56 * t22
      t20 = t55 * t29 * t20
      t60 = t2 * t11
      t61 = t46 * t8
      t62 = t61 * t3
      t57 = t62 * (t51 * (t40 * (-t57 * t11 * t33 * t42 + t17 * t37 * t3
     #3 * t4) + t35 + t44 + t39 * t42 * t37 * t33 + (-t24 * t35 * t11 * 
     #t33 - t57 * t25) * t4 * t38 + t58 * t16 * t35) + t43 * t2 * (t21 *
     # (t59 - t45) + t60 * (t20 - t7)))
      t59 = t59 - t45
      t60 = t60 * t7
      t63 = t14 * t8
      t64 = t21 * t22
      t65 = t50 * t48
      t7 = t2 * (-t27 * t28 * (t1 * t30 + t63) * t31 - t64 * t62 * t51 *
     # t14) - t28 * (t46 * t9 + t8) * t3 * t35 + t65 * ((t21 * t59 - t60
     #) * t26 * t4 + t23 * t2 * (-t55 + (-t20 + t7) * t11 * t15)) * t19 
     #- t10 * t52 * t23 * t4 * t37 * t34
      t10 = t48 * t16
      t11 = t54 * t51 * t3 * (t24 * t49 * t4 * t33 + t30 * t22 * (-t64 *
     # t14 + t2) + t40 * t42 * t21 * t8 * t33 - t10 * t25)
      t10 = t3 * (t46 * t2 * (t22 * (t63 * t21 - t23) + t48) - t40 * t4 
     #* t37 * t8 * t33 + t10 * t27 * (t23 * t25 + t60 * t30) + t49 * t21
     # * t25 * (t58 + t46) * t4) + t53 * t35 * (t2 * t55 + t64) + t61 * 
     #t35 * t4 * (t45 * t21 + t60) * t9 + t54 * t59 * t31
      t16 = 64
      ret = t16 * (t4 * (t49 * t15 * t21 * t33 * (t24 * t45 - 1) * t18 *
     # t23 * t12 + t47 * (t22 * (t35 + t44) - t40 * t35 * t13 ** 2 * t15
     # * t37 * t34 - t45 * t35 * t15) * t8 * t23) - t47 * t1 * t39 * t23
     # * t4 * t42 * t37 * t15 * t8 * t34 + t50 * t2 * t31 * t28 * t29 + 
     #t51 * t39 * t23 * t42 * t37 * t15 * t8 * t34) - 128 * t54 * t2 * (
     #t53 * (-t17 * t38 * t15 * t33 + t36) * t32 * t21 * t41 * t4 + (t39
     # * t13 * t15 * t8 * t32 * t34 * t30 * t37 - t52 * t32 * t33 * t30 
     #* t21) * t41 * t42 + t43 * t3 * t19 * t28 * (t15 * t30 * t29 - t8)
     #) + 8 * t57 + 16 * t7 - 32 * t50 * t19 * t8 * t37 * t4 * (t36 * t2
     #6 * t56 + t38 * t15 * t33 * (t1 * t23 - t2 * (t1 * t6 - t14 * t9))
     #) + 48 * t65 * t23 * t21 * t19 * (t45 * t15 - t22) - 24 * t11 + 4 
     #* t10 + 12 * t5 * t38 * t3 * t23 * t4 * t21 * t8 * t33

      hjetmass_qqbgg_triangle_pmpm_s123_mhsq_0_dp = ret/32d0/(0,1d0)
      return

      end function
