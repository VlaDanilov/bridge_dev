
      double complex function hjetmass_box_ppmm_0_s34_0_mhsq_s134_s234_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i3, i4)
      t2 = zb(i2, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = zb(i4, i3)
      t6 = za(i2, i4)
      t7 = zb(i4, i2)
      t8 = t3 * t4
      t9 = t6 * t7
      t10 = t9 + t8
      t11 = za(i1, i2)
      t12 = za(i1, i4)
      t13 = za(i1, i3)
      t14 = zb(i3, i1)
      t15 = zb(i4, i1)
      t16 = t13 * t4
      t17 = t12 * t7
      t18 = t17 + t16
      t19 = t2 ** 2
      t20 = t11 ** 2 * t19
      t21 = t20 * t18
      t22 = t14 * t3 + t15 * t6
      t23 = mt ** 2
      t24 = t9 + t8
      t25 = t11 * t2
      t26 = t25 * t1 * t5
      t14 = t13 * t14
      t24 = t25 * (t24 * t15 * t12 + t14 * t24 - t26)
      t27 = 16 * t23 * t11 * t2 * t22 * t21 + 4 * t24 ** 2
      t27 = cdsqrt(t27)
      t28 = t10 * t15 * t12 + t14 * t10 - t26
      t21 = 0.1D1 / t21
      t16 = t17 + t16
      t8 = -t8 * t16 - t9 * t16
      t8 = -t20 * t21 * (t12 * t8 * t15 + t14 * t8 + t26 * t16) / 2
      t9 = t25 * t27 * t21 * t18 / 4
      t14 = -t28 + t8 - t9
      t16 = t25 * t22
      t16 = 0.1D1 / t16
      t18 = 0.1D1 / t11
      t10 = t10 * t18
      t18 = t10 * t12
      t20 = t14 * t16
      t22 = t20 * t6
      t25 = t22 * t2 + t18
      t24 = 2 * t24
      t26 = -t27 + t24
      t29 = 0.1D1 / t2
      t30 = t10 * t29 * t13
      t20 = t20 * t3
      t31 = t26 * t21 / 4
      t32 = t15 * (t30 + t20) - t31 * t13 * t7
      t8 = -t28 + t8 + t9
      t9 = t8 * t16
      t6 = t9 * t6
      t28 = t6 * t2 + t18
      t10 = t10 * t13
      t20 = t20 * t2 + t10
      t9 = t9 * t3
      t10 = t9 * t2 + t10
      t18 = t18 * t29
      t22 = t15 * (t18 + t22) - t17 * t31
      t24 = t27 + t24
      t6 = t15 * (t18 + t6) - t17 * t24 * t21 / 4
      t9 = t15 * (t30 + t9) - t24 * t21 * t13 * t7 / 4
      t18 = 0.1D1 / t22
      t6 = 0.1D1 / t6
      t4 = 0.1D1 / t4
      t22 = 0.1D1 / t13
      t3 = 0.1D1 / t3
      t5 = 0.1D1 / t5
      t27 = 0.1D1 / t7
      t29 = t1 * t18
      t30 = -t1 * t6 + t5
      t31 = t1 * t7
      t33 = -t31 + t10
      t34 = t26 + t24
      t35 = t1 ** 2
      t36 = t20 ** 2
      t37 = t36 * t18
      t38 = t37 * t27
      t39 = t24 * t6
      t40 = t26 * t18
      t41 = t10 * t6
      t42 = t41 * t27
      t43 = t24 * t10
      t44 = t7 * t5
      t31 = (-t31 + t20) * t18
      t45 = t28 * t24
      t46 = (-t29 + t5) * t26
      t21 = t21 * t3
      t47 = t21 * t4
      t48 = t25 * t20
      t49 = t4 * t27
      t21 = t21 * t2
      t50 = t20 + t10
      t51 = t10 ** 2
      t52 = t51 * t6
      t53 = t3 * t4 * t2
      t54 = t2 * t12
      t55 = t25 * t18
      t56 = t12 * t22
      t57 = t28 * t10
      t58 = t10 * t22
      t59 = t47 * t2
      t60 = t2 * t32
      t61 = t28 * t6
      t62 = t20 * t18
      t3 = t4 * (t19 * (t27 * (-t48 * t18 - t41 * t28) + t5 * (-t25 - t2
     #8)) + (t25 * (-t38 * t15 * t22 + t60 * t22 * t5 + t62 * t2 * (t22 
     #* t32 + t15) * t27) + t28 * (-t52 * t15 * t22 * t27 + t42 * t2 * (
     #t22 * t9 + t15))) * t3 * t1 + (t61 * (t58 - t2) * t15 - t19 * t23 
     #* (t18 + t6)) * t3 * t35) + t1 * (t2 * ((t4 * (t5 * (t25 + t28) - 
     #t29 * t25) + t5 * t50 * t27) * t3 * t15 + t2 * t4 * (t61 + t55)) +
     # t49 * t19 * t3 * (t41 + t62) * t23 + t3 * (t5 * (t15 * (t27 * (-t
     #36 - t51) + t4 * (-t57 - t48)) + t2 * (t27 * (t10 * t9 + t20 * t32
     #) + t28 * t9 * t4)) + t4 * (-t61 * t2 * t9 + t55 * (t15 * t20 - t6
     #0)) * t1) * t22) + t21 * t5 * t27 * (t24 * t51 + t26 * t36) - t19 
     #* t27 * t5 * (t20 + t10)
      t9 = t53 * t22 * t35 * t23 * (t41 + t44 + t62)
      t4 = t21 * t1 * (t5 * (t20 * t26 + t43) + t17 * (t46 * t20 + t43 *
     # t30) * t22 * t4)
      t10 = t56 * t59 * t1 * (t37 * t26 + t39 * t51)
      ret = -8 * t53 * t1 * t23 * (t5 * (t22 * t50 + t2) + t22 * (t52 + 
     #t37) * t27 + t22 * (t18 + t6) * t7 * t35) + 3 * t4 + 5 * t10 + t21
     # * (t35 * (t31 * t26 + t39 * t33) + (t34 * t5 * t1 + t49 * (t40 * 
     #t48 + t45 * t41)) * t13 * t2) + t47 * t19 * (t12 * (t26 * (-t20 * 
     #t5 - t38) - t43 * (t42 + t5) + (-t40 - t39) * t7 * t35 + t44 * t34
     # * t1) + t13 * (t46 * t25 + t45 * t30) + t35 * ((t31 + t44) * t26 
     #* t14 + t24 * t8 * (t33 * t6 + t44)) * t22 * t16 * t11) + 2 * t59 
     #* (t26 * (t20 * (-t25 * t5 + t29 * (t54 + t25)) + t36 * (-t55 * t2
     #7 - t56 * t5) - t12 * t18 * t22 * t27 * t20 * t36) + t43 * (t5 * (
     #-t58 * t12 - t28) + t6 * (t1 * (t54 + t28) - t57 * t27 - t56 * t51
     # * t27))) - 4 * t3 + 16 * t9

      hjetmass_box_ppmm_0_s34_0_mhsq_s134_s234_dp = ret/32d0/(0,1d0)
      return

      end function
