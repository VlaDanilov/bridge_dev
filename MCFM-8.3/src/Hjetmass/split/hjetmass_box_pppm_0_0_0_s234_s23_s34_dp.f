
      double complex function hjetmass_box_pppm_0_0_0_s234_s23_s34_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = zb(i2, i1)
      t2 = zb(i3, i1)
      t3 = za(i3, i4)
      t4 = zb(i4, i3)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = mt ** 2
      t10 = t3 * t4
      t11 = t5 * t6
      t12 = 0.1D1 / t6
      t13 = 0.1D1 / t5
      t14 = (0.1D1 / 0.2D1)
      t5 = t14 * cdsqrt(-t11 * t10 * (-2 * t10 * t5 * t6 - 8 * t9 * t7 * 
     #t8)) * dsqrt(0.2D1) * t13 * t12
      t15 = t10 - t5
      t11 = t11 * t10
      t9 = t11 * (4 * t9 * t7 * t8 + t11)
      t9 = cdsqrt(t9)
      t16 = t9 + t11
      t8 = 0.1D1 / t8
      t4 = 0.1D1 / t4
      t17 = t15 * t8 * t1
      t18 = t16 * t13 * t12 * t4 * t2
      t19 = t18 + t17
      t20 = za(i1, i2)
      t21 = za(i1, i3)
      t22 = 0.1D1 / t7
      t23 = 0.1D1 / t3
      t24 = t20 * t1
      t18 = t14 * (t17 * t23 * t21 + t18 * t22 * t20)
      t25 = t18 + t24
      t26 = zb(i4, i1)
      t5 = t10 + t5
      t9 = t9 - t11
      t10 = t5 * t8 * t1
      t11 = t9 * t13 * t12 * t4 * t2
      t12 = t11 - t10
      t11 = t14 * (-t10 * t23 * t21 + t11 * t22 * t20)
      t20 = -t11 + t24
      t24 = za(i1, i4)
      t21 = t2 * t21 + t24 * t26
      t11 = -t21 - t11
      t18 = -t21 + t18
      t21 = 0.1D1 / t25
      t20 = 0.1D1 / t20
      t24 = 0.1D1 / t24
      t25 = 0.1D1 / t26
      t18 = 0.1D1 / t18
      t11 = 0.1D1 / t11
      t26 = t10 + t12
      t27 = t21 - t18
      t28 = t20 - t11
      t29 = t15 ** 2
      t30 = t9 ** 2
      t31 = t13 ** 2
      t32 = t19 ** 2
      t33 = t12 ** 2
      t34 = t4 ** 2
      t35 = t16 ** 2
      t36 = t9 * t11
      t37 = t12 * t20
      t38 = t37 * t9
      t39 = t19 * t21
      t40 = t39 * t16
      t41 = t1 * t23
      t42 = t27 * t16
      t43 = t41 * t2
      t44 = t6 * t24
      t45 = t37 * t5
      t46 = t2 * t22
      t47 = t27 * t32
      t48 = -t20 + t11
      t49 = -t21 + t18
      t50 = t2 ** 2
      t51 = t48 * t33
      t52 = t1 ** 2 * t7
      t53 = t50 * t3 * t22
      t54 = t16 * t15
      t55 = t49 * t19
      t56 = t48 * t9
      t57 = t49 * t32
      t58 = t48 * t5
      t59 = t27 * t19
      t60 = t2 * t3
      t49 = t25 * (t6 * (t46 * (t60 * (t12 * t48 + (t15 * t49 + t58) * t
     #8 * t1 + t59) + t57 + t51 + (t58 * t12 + t59 * t15) * t8 * t1) + t
     #44 * ((t37 - t39) * t7 * t1 + t60 * (-t37 + t39))) + t50 * t1 * (t
     #56 + t42) * t4 * t13)
      t58 = t19 * t18
      t59 = t28 * t12
      t60 = t25 * t6 * t2 * (t1 * (t59 + t55 + t60 * (t21 + t20 - t18 - 
     #t11)) + t24 * ((-t20 * t1 - t11 * t22 * (t60 + t12)) * t9 * t5 + t
     #54 * (t1 * t21 + t18 * t22 * (t60 - t19))) * t4 * t13 * t8 + t44 *
     # t3 * (t11 * t12 + (t21 + t20) * t7 * t1 + t60 * (t18 + t11) - t58
     #))
      t61 = t15 * t32
      t62 = t61 * t21
      t63 = t5 * t33
      t64 = t44 * t8
      t65 = t22 * t6
      t66 = t18 * t32
      t10 = t23 * t25 * ((t44 * (t16 * (t58 * t17 - t66) + t9 * t12 * (t
     #11 * t10 - t37)) + t46 * (t66 * t16 - t36 * t33)) * t4 * t13 + t65
     # * t8 * t1 * (t63 * t11 + t61 * t18)) + t25 * (t65 * (t23 * (t12 *
     # t33 * t48 + t27 * t19 * t32 + (-t63 * t20 - t62) * t8 * t1) + t64
     # * (t39 * t1 * t29 * t8 - t45 * t26 - t62)) + t23 * t24 * t2 * (t5
     #9 * t30 + t55 * t35) * t34 * t31 + t23 * (t9 * t33 * (t44 * t11 + 
     #t46 * t20) + (t44 - t46) * t21 * t32 * t16) * t4 * t13)
      t27 = t41 - t46
      ret = t14 * t10 - (0.5D1 / 0.2D1) * t64 * t4 * t13 * t25 * (t40 * 
     #t15 * t27 + t38 * t27 * t5) - t25 * ((t2 * (-t36 * t41 * t12 + t46
     # * (t56 * t12 + t55 * t16)) + t44 * ((t23 * (-t52 * t20 + t51 * t2
     #2) - t53 * t20) * t9 * t5 + t54 * (t23 * (t52 * t21 + t47 * t22) +
     # t53 * t21)) * t8) * t4 * t13 + t50 * t24 * (t11 * t30 + t18 * t35
     #) * t34 * t31 + t41 * t6 * (t51 + t57)) - t25 * ((t44 * (t2 * (-t3
     #6 * t26 + (t17 - t19) * t18 * t16) - t41 * (t40 + t38) * t7) + t43
     # * (t42 * t19 + t38)) * t4 * t13 + t43 * t7 * t24 * (t20 * t30 + t
     #21 * t35) * t34 * t31 + t24 * t6 ** 2 * (t28 * t33 + t46 * (t8 * (
     #-t39 * t15 + t45) + t1 * (t20 * t5 ** 2 + t21 * t29) * t8 ** 2) * 
     #t3 + t47)) + 2 * t49 + 4 * t60

      hjetmass_box_pppm_0_0_0_s234_s23_s34_dp = ret/32d0*(0,1d0)
      return

      end function
