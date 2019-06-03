
      double complex function hjetmass_box_ppmm_0_0_s12_mhsq_s34_s123_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i1, i4)
      t8 = za(i3, i4)
      t9 = za(i2, i4)
      t10 = t7 * t4
      t11 = t9 * t6
      t12 = t11 + t10
      t13 = t8 * zb(i4, i3)
      t14 = t13 * t12
      t15 = t1 * t2
      t16 = t3 * t4
      t17 = t5 * t6
      t18 = t13 * (t17 + t15 + t16)
      t19 = zb(i4, i1)
      t20 = zb(i4, i2)
      t21 = t5 * t20
      t22 = t19 * t3 + t21
      t23 = mt ** 2
      t24 = 16 * t23 * t22 * t14 + 4 * t18 ** 2
      t24 = cdsqrt(t24)
      t14 = 0.1D1 / t14
      t25 = t17 + t16
      t26 = t13 * t14
      t10 = t26 * (t10 * t25 + t11 * t25 + t15 * (t11 + t10))
      t25 = 2
      t27 = (0.1D1 / 0.2D1)
      t15 = t25 * (t17 + t15 + t16)
      t12 = t27 * t24 * t14 * t12
      t16 = -t10 + t15 - t12
      t17 = t25 * t18
      t18 = -t24 - t17
      t22 = 0.1D1 / t22
      t21 = t21 * t27
      t28 = t21 * t16 * t22
      t29 = t6 * (t5 + t9 * t18 * t14 / 4) - t28
      t10 = -t10 + t15 + t12
      t12 = t24 - t17
      t15 = t21 * t10 * t22
      t5 = t6 * (t5 + t9 * t12 * t14 / 4) - t15
      t17 = t4 * t20
      t21 = t6 * t19
      t24 = -t17 + t21
      t30 = t10 * t12
      t31 = t30 * t14 * t22 * t8 * t24
      t32 = t16 * t18
      t24 = t32 * t14 * t22 * t8 * t24
      t33 = t18 * t14 / 4
      t34 = t12 * t14 / 4
      t28 = -t33 * t11 + t28
      t11 = -t34 * t11 + t15
      t15 = 0.1D1 / t7
      t5 = 0.1D1 / t5
      t35 = 0.1D1 / t6
      t1 = 0.1D1 / t1
      t9 = 0.1D1 / t9
      t29 = 0.1D1 / t29
      t36 = 0.1D1 / t20
      t37 = t16 ** 2
      t38 = t10 ** 2
      t39 = t18 * t37
      t40 = t12 * t38
      t41 = t40 + t39
      t42 = t13 * t15 + t19
      t43 = t22 ** 2
      t44 = t8 ** 2
      t45 = t44 ** 2
      t46 = t8 * t44
      t13 = t1 * (t19 * t7 + t13)
      t47 = t8 * t29
      t48 = t7 * t35
      t49 = t48 * t1
      t50 = t2 * t29
      t51 = t8 * t5
      t52 = t2 * t5
      t18 = t18 * t29
      t12 = t12 * t5
      t53 = t36 * t9
      t54 = t9 * t1
      t55 = t4 * t44
      t56 = t12 * t38
      t57 = t18 * t37
      t58 = t57 + t56
      t59 = t10 * t5
      t60 = t16 * t29
      t61 = t9 * t15
      t62 = t32 * t29 + t30 * t5
      t63 = t32 + t30
      t64 = t6 ** 2
      t65 = t19 * t36
      t66 = t14 * t1
      t67 = t20 * t43
      t68 = t14 * t15
      t62 = t46 * (t68 * (-t67 * t58 + (-t12 - t18) * t36 * t64) * t2 + 
     #t22 * ((t65 * t62 * t14 * t64 - t19 * t20 * t22 * (t29 * t37 + t38
     # * t5) + t4 * t62 * t14 * t6) * t15 * t8 + t66 * (t65 * t63 * t6 +
     # t4 * t63)) * t9 - t66 * t67 * t15 * t41)
      t56 = t57 + t56
      t57 = t5 + t29
      t63 = t16 + t10
      t64 = t15 * t2
      t3 = t44 * (t2 * (t36 * (t1 * (t15 * (-t27 * t10 * t22 * t3 * t20 
     #- t27 * t16 * t22 * t3 * t20 + t33 * t7 * t6 + t34 * t7 * t6) + t9
     # * (t28 + t11)) - t9 * t57 * t2 * t23) - t54 * t23 * t63 * t35 * t
     #22) + t6 * (t55 * t61 * t57 + t64 * (t53 * (t11 * t5 + t28 * t29) 
     #+ t5 + t29) * t8))
      t7 = t60 + t59
      t33 = t63 * t1
      t7 = t22 * t46 * (t9 * (t33 * t19 + t64 * (t59 * (t23 + t11) + t60
     # * (t23 + t28)) + t21 * t15 * t7 * t8 + t17 * t8 * t15 * t7) + t15
     # * t20 * (t2 * t7 + t33) + t68 * t6 * (t30 * (t52 + t1) + t32 * (t
     #50 + t1)))
      ret = t25 * t7 + t27 * t44 * (t2 * (t32 * (-t9 * (t47 + t49) * t4 
     #- t50 + t9 * (t47 * t42 * t6 + t13) * t36) + t30 * (-t9 * (t51 + t
     #49) * t4 - t52 + t9 * (t51 * t42 * t6 + t13) * t36)) * t14 * t22 +
     # t53 * (-t1 * (t31 + t24) + (-t24 * t29 - t31 * t5) * t15 * t8 * t
     #6) + t44 * t19 * t20 * t15 * t9 * (t12 * t10 * t38 + t18 * t16 * t
     #37) * t14 * t22 * t43 + t54 * t8 * (-t19 * t41 + t17 * (-t40 - t39
     #) * t35) * t14 * t43) - 8 * t1 * t8 * (t9 * (-t48 * t23 * t2 ** 2 
     #* t36 + t55) + t44 * t6 * t15) - (0.5D1 / 0.4D1) * t21 * t61 * t14
     # * t43 * t45 * t56 - (0.3D1 / 0.4D1) * t61 * t17 * t14 * t43 * t45
     # * t56 - t61 * t22 * t46 * (t26 * t58 * t22 * t2 - t60 * t24 - t59
     # * t31) / 4 + t62 - 4 * t3

      hjetmass_box_ppmm_0_0_s12_mhsq_s34_s123_dp = ret/32d0/(0,1d0)
      return

      end function
