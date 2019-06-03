
      double complex function hjetmass_box_pmpm_0_s14_0_mhsq_s124_s134_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i2)
      t4 = za(i1, i3)
      t5 = zb(i2, i1)
      t6 = za(i1, i4)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = zb(i4, i1)
      t10 = zb(i4, i2)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t1 * t12
      t14 = t3 * t2
      t15 = t13 + t14
      t16 = t7 ** 2
      t17 = t16 * t8 ** 2
      t18 = t17 * t15
      t19 = t10 * t11 + t4 * t5
      t20 = mt ** 2
      t21 = t1 * t10
      t22 = t3 * t5
      t23 = t22 + t21
      t24 = t7 * t8
      t25 = t24 * t6 * t9
      t26 = t4 * t2
      t27 = t11 * t12
      t28 = t24 * (t26 * t23 + t27 * t23 - t25)
      t29 = 16 * t20 * t7 * t8 * t19 * t18 + 4 * t28 ** 2
      t29 = cdsqrt(t29)
      t18 = 0.1D1 / t18
      t23 = t13 * t23 + t14 * t23
      t21 = t22 + t21
      t21 = t26 * t21 + t27 * t21 - t25
      t15 = t24 * t29 * t18 * t15 / 4
      t17 = t17 * t18 * (t26 * t23 + t25 * (-t14 - t13) + t27 * t23) / 2
      t22 = -t21 + t17 - t15
      t19 = t24 * t19
      t15 = -t21 + t17 + t15
      t17 = 2 * t28
      t21 = -t29 + t17
      t19 = 0.1D1 / t19
      t23 = 0.1D1 / t8
      t24 = 0.1D1 / t7
      t25 = (t27 + t26) * t24
      t26 = t25 * t23
      t27 = t26 * t3
      t28 = t22 * t19
      t30 = t28 * t4
      t14 = t14 / 4
      t31 = -t14 * t21 * t18 + t5 * (t27 + t30)
      t17 = t29 + t17
      t29 = t15 * t19
      t4 = t29 * t4
      t14 = t5 * (t27 + t4) - t14 * t17 * t18
      t26 = t26 * t1
      t27 = t28 * t11
      t21 = -t5 * (t27 + t26) + t21 * t18 * t1 * t2 / 4
      t11 = t29 * t11
      t17 = t17 * t18 * t1 * t2 / 4 - t5 * (t26 + t11)
      t18 = t25 * t1
      t26 = t27 * t8 + t18
      t11 = t11 * t8 + t18
      t18 = t25 * t3
      t25 = t30 * t8 + t18
      t4 = t4 * t8 + t18
      t9 = 0.1D1 / t9
      t14 = 0.1D1 / t14
      t18 = 0.1D1 / t12
      t27 = 0.1D1 / t3
      t28 = 0.1D1 / t31
      t29 = 0.1D1 / t6
      t30 = t22 + t15
      t31 = t15 ** 2
      t32 = t15 * t31
      t33 = t19 ** 2
      t34 = t1 ** 2
      t35 = t2 ** 2
      t36 = t22 ** 2
      t37 = t22 * t36
      t38 = t31 * t11
      t39 = t36 * t26
      t40 = t28 * t36
      t41 = t31 * t14
      t42 = t9 * t29
      t43 = t2 * t6
      t44 = t8 * t21
      t45 = t8 * t17
      t46 = t34 * t12 * t27
      t47 = t1 * t27
      t48 = t22 * t28
      t49 = t48 * t21
      t50 = t15 * t14
      t51 = t50 * t17
      t52 = t2 * t1
      t53 = t29 * t1
      t54 = t27 * t18
      t16 = t54 * t8 * t16
      t55 = -t50 - t48
      t56 = t50 + t48
      t57 = t36 + t31
      t58 = t22 * t25
      t59 = t15 * t4
      t60 = t2 * t29
      t61 = t15 * t11
      t62 = t22 * t26
      t63 = t23 * t2
      t64 = t50 + t48
      t65 = t1 * t8
      t66 = t35 * t9
      t67 = t20 * t1
      t6 = -t6 * t7 * t56 * t18 * t2 * t35 + t46 * (t51 + t49) + (t50 * 
     #(t11 * t7 + t17 * t3 - t67) + t48 * (t21 * t3 + t26 * t7 - t67)) *
     # t18 * t35
      t6 = t19 * t6 + t19 * (t34 * t2 * t5 * t64 + t66 * (t53 * (t59 + t
     #58) * t23 * t5 + t2 * t7 * t30) * t18 + t47 * (t1 * (t1 * (t64 * t
     #12 * t5 - t8 * t29 * t30) + t63 * t42 * (t22 * (-t10 * t26 + t25 *
     # t5) + t59 * t5)) + t5 * (t40 * (t65 + t43 - t26) + t41 * (t65 + t
     #43 - t11)) * t19 * t7)) + t1 * (t19 * (t54 * t35 * t10 * t9 * t57 
     #* t19 * t7 + t63 * (t5 * (-t50 * t11 - t48 * t26 + t43 * t56) + t2
     # * (t43 * t55 + t61 * (-t42 + t14) + t62 * (-t42 + t28)) * t18 * t
     #10) + t47 * (t23 * (t15 * (-t42 * t2 * t10 * t11 + (t11 * t12 + t4
     # * t2) * t14 * t5) + t48 * t5 * (t12 * t26 + t2 * t25)) - t60 * (t
     #59 + t58) * t18)) - t63 * t46 * (t14 + t28) * t24 * t20) - t34 * t
     #35 * t23 * (t14 + t28) * t24 * t20
      t48 = t25 + t4
      t15 = t22 + t15
      t22 = t52 * t23
      t50 = t62 + t61
      t17 = t19 * ((t18 * t2 * (t40 * (t5 * (t43 - t26) + t44) + t41 * (
     #t5 * (t43 - t11) + t45)) + t47 * (t8 * (t41 * t17 + t40 * t21) + (
     #-t36 - t31) * t9 * t5 * t2)) * t19 * t7 + t52 * (t29 * t18 * t50 +
     # t47 * t55 * t20))
      t21 = t22 * (-t15 * t5 * t9 * t19 * t2 + t53 * t24 * (-t47 * t48 +
     # t11 + t26) + t66 * t10 * t18 * t15 * t19)
      t54 = -8
      ret = t54 * (t16 * (t28 * t37 * (t5 * (-t43 + t26) - t44) + (t5 * 
     #(-t43 + t11) - t45) * t14 * t32) * t19 * t33 + t53 * (t34 * (t13 *
     # t27 + t2) + t2 * (t35 * t3 * t18 - t46) * t23 * t9 * t20) * t24 +
     # t52 * t27 * t18 * t7 * (t20 * t8 * (t40 + t41) + t42 * (t10 * (t3
     #9 + t38) + t5 * (-t25 * t36 - t31 * t4))) * t33 + t52 * (-t49 - t5
     #1 + t42 * (t30 * t18 * t2 - t47 * t30) * t20) * t19) + 12 * t17 + 
     #6 * t21 + 4 * t6 - 2 * t22 * (t47 * (-t2 * t10 * t15 + t15 * t12 *
     # t5) * t9 * t19 + t60 * t24 * t18 * (t1 * t48 + t3 * (-t26 - t11))
     #) + 16 * t19 * (-t16 * t2 * t5 * t9 * (t37 + t32) * t33 + (-t66 * 
     #t5 * t57 + t53 * (t39 + t38) * t27 * t8) * t18 * t19 * t7 + t34 * 
     #t27 * t29 * t50)

      hjetmass_box_pmpm_0_s14_0_mhsq_s124_s134_dp = ret/32d0/(0,1d0)
      return

      end function
