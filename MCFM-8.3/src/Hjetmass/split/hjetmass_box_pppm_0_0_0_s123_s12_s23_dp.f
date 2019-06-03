      double complex function hjetmass_box_pppm_0_0_0_s123_s12_s23_dp 
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
      t2 = za(i2, i3)
      t3 = zb(i2, i1)
      t4 = zb(i3, i2)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = mt ** 2
      t8 = t2 * t4
      t9 = t8 * t1 * t3
      t10 = t9 * (4 * t7 * t5 * t6 + t9)
      t10 = cdsqrt(t10)
      t11 = t9 + t10
      t12 = 0.1D1 / t1
      t13 = 0.1D1 / t3
      t1 = cdsqrt(-t1 * t3 * t8 * (-8 * t7 * t5 * t6 - 2 * t9))*dsqrt(0.
     #2D1) * t12 * t13 / 2
      t7 = t8 - t1
      t14 = za(i2, i4)
      t15 = za(i1, i4)
      t5 = 0.1D1 / t5
      t2 = 0.1D1 / t2
      t16 = t7 * t2
      t17 = t16 * t14
      t18 = t11 * t12 * t5
      t19 = t18 * t13 * t15
      t20 = t17 + t19
      t21 = za(i3, i4)
      t22 = zb(i4, i1)
      t23 = zb(i4, i2)
      t24 = t14 * t23 + t21 * zb(i4, i3)
      t25 = 0.1D1 / t4
      t26 = 0.1D1 / t6
      t19 = t17 * t26 * t22 / 2 + t19 * t25 * t23 / 2
      t27 = -t24 + t19
      t9 = -t9 + t10
      t1 = t8 + t1
      t8 = t9 * t12 * t5 * t13
      t10 = t1 * t2
      t13 = t10 * t26 * t14 * t22 / 2 - t8 * t25 * t15 * t23 / 2
      t23 = -t24 + t13
      t24 = t15 * t22
      t19 = t24 + t19
      t13 = t24 + t13
      t24 = t10 * t14
      t8 = t8 * t15 - t24
      t27 = 0.1D1 / t27
      t19 = 0.1D1 / t19
      t28 = 0.1D1 / t14
      t22 = 0.1D1 / t22
      t23 = 0.1D1 / t23
      t13 = 0.1D1 / t13
      t29 = 0.1D1 / t21
      t30 = t27 - t19
      t31 = -t27 + t19
      t32 = -t23 + t13
      t33 = t23 - t13
      t34 = t5 ** 2
      t35 = t14 ** 2
      t36 = t9 ** 2
      t37 = t15 ** 2
      t38 = t12 ** 2
      t39 = t20 ** 2
      t40 = t11 ** 2
      t41 = t8 * t28
      t42 = t3 * t29
      t43 = t8 * t9
      t44 = t20 * t28
      t45 = t20 * t31 * t40
      t46 = t36 * t8 * t33
      t47 = t42 * t25
      t48 = t12 * t22
      t49 = t48 * t25 * t5
      t50 = 0.1D1 / t15
      t51 = t7 * t27
      t52 = t51 * t11
      t53 = t1 * t23
      t54 = t53 * t9
      t55 = t23 * t36 + t27 * t40
      t56 = t3 ** 2
      t57 = t16 * t27
      t58 = t11 * t19
      t59 = t58 * t20
      t60 = (-t54 + t52) * t2
      t61 = t20 * t27
      t62 = t6 * t28
      t63 = t39 * t27
      t64 = t8 * t23
      t65 = t2 * t14
      t66 = t65 * t3
      t67 = t9 * t13
      t8 = t22 * (-t62 * t40 * t38 * t34 * t25 * t37 * t19 + t3 * (-t67 
     #* (t24 + t8) - t59) * t12 * t5 + t50 * t3 * ((-t57 * t20 + t64 * (
     #t41 + t10) + t63 * t28) * t21 * t6 + t66 * (-t51 * t20 + t53 * t8)
     #)) + t22 * (t25 * t15 * (-t62 * t15 * t36 * t13 + t42 * (t15 * t55
     # + (t46 * t1 + t45 * t7) * t26 * t2 * t25)) * t38 * t34 + (t58 * t
     #17 * t3 + t62 * (-t61 * t11 - t43 * t23) * t21 + t3 * (t6 * (t11 *
     # (t44 * t19 + t57) + t9 * (-t10 * t23 + t41 * t13)) + t42 * (t43 *
     # t13 + t60 * t14 + t59)) * t25 * t15) * t12 * t5 + t56 * t50 * (t2
     #3 * t8 ** 2 + t63))
      t17 = -t64 + t61
      t24 = t64 - t61
      t57 = t56 * t29
      t59 = t6 ** 2 * t21 * t28
      t61 = t21 * t4
      t63 = t22 * t3
      t17 = t63 * (t3 * (t18 * t15 * t19 * t29 * t14 - t50 * t4 * t2 * (
     #t53 + t51) * t35) + t6 * (t61 * t50 * t17 + (-t67 + t58) * t12 * t
     #5 * t15) + t57 * t24 * t14) + t22 * (t15 * (t65 * t47 * t26 * (t36
     # * t1 * t13 + t40 * t7 * t19) * t38 * t34 + (t9 * (t57 * t33 * t14
     # - t59 * t23 + t64 * t3 * (t62 + t42) * t25) + (t44 * t3 * t6 * t2
     #5 + t57 * (t20 * t25 - t14) + t59) * t27 * t11) * t12 * t5) + t3 *
     # (t59 * t24 + (t17 * t3 + (-t53 - t51) * t2 * t21 * t6) * t50 * t1
     #4 * t4) - t62 * t25 * t55 * t38 * t34 * t37)
      t18 = t27 + t23
      t21 = t52 * t20
      t50 = t53 * t43
      ret = (0.3D1 / 0.2D1) * t57 * t65 * t49 * t26 * (t21 + t50) + (0.5
     #D1 / 0.2D1) * t49 * t2 * t3 * ((t53 * t36 + t51 * t40) * t29 * t26
     # * t12 * t5 * t15 * t14 + t21 + t50) + 5 * t66 * t48 * t5 * (-t52 
     #+ t54) - t49 * (t3 * (t11 * (t16 * t19 * t20 + t28 * t30 * t39) + 
     #t43 * (t10 * t13 + t41 * t32) + t42 * t2 ** 2 * (t9 * t1 ** 2 * t2
     #3 - t11 * t7 ** 2 * t27) * t26 * t35) + t47 * (t36 * (t10 * t33 + 
     #t41 * t33) + t40 * (t16 * t30 + t44 * t31)) * t12 * t5 * t37 + t29
     # * t28 * t25 * (t30 * t11 * t40 + t32 * t9 * t36) * t38 * t34 * t1
     #5 * t37 + t28 * (t46 + t45) * t12 * t5 * t15) / 2 + t8 - 2 * t17 -
     # 4 * t63 * (t6 * (t61 * t18 * t6 + t24 * t3) + t42 * (t60 * t26 * 
     #t12 * t5 + t3 * t4 * t18) * t35) - 8 * t14 * t56 * t6 * t4 * t22 *
     # (t27 + t23)

      hjetmass_box_pppm_0_0_0_s123_s12_s23_dp = ret/32d0*(0,1d0)
      return

      end function
