
      double complex function hjetmass_box_ppmm_0_0_0_s123_s12_s23_dp 
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
      t5 = za(i1, i2)
      t6 = za(i1, i3)
      t7 = zb(i3, i1)
      t8 = mt ** 2
      t9 = t3 * t4
      t10 = t5 * t2
      t11 = 0.1D1 / t2
      t12 = 0.1D1 / t5
      t5 = cdsqrt(-t10 * t9 * (-2 * t9 * t5 * t2 - 8 * t8 * t6 * t7))*d
     #sqrt(0.2D1) * t12 * t11 / 2
      t13 = -t5 + t9
      t14 = zb(i4, i1)
      t15 = za(i2, i4)
      t10 = t10 * t9
      t16 = t10 * (4 * t8 * t6 * t7 + t10)
      t16 = cdsqrt(t16)
      t17 = t16 + t10
      t18 = za(i1, i4)
      t19 = zb(i4, i2)
      t6 = 0.1D1 / t6
      t3 = 0.1D1 / t3
      t7 = 0.1D1 / t7
      t4 = 0.1D1 / t4
      t20 = t6 * t18
      t21 = t20 * t17 * t12 * t11
      t22 = t3 * t15
      t23 = t22 * t13
      t24 = t1 * zb(i4, i3) + t15 * t19
      t25 = t23 * t7 * t14 / 2 + t21 * t4 * t19 / 2
      t26 = t25 - t24
      t5 = t5 + t9
      t9 = -t16 + t10
      t10 = t20 * t9 * t12 * t11
      t11 = t22 * t5
      t16 = t10 * t4 * t19 / 2 + t11 * t7 * t14 / 2
      t19 = t16 - t24
      t24 = t18 * t14
      t25 = t25 + t24
      t27 = -t22 + t20
      t28 = t17 * t12 * t4
      t16 = t16 + t24
      t24 = t9 * t12 * t4
      t26 = 0.1D1 / t26
      t18 = 0.1D1 / t18
      t19 = 0.1D1 / t19
      t25 = 0.1D1 / t25
      t14 = 0.1D1 / t14
      t16 = 0.1D1 / t16
      t29 = t2 * t18
      t30 = t29 * t14
      t31 = t30 + t12
      t32 = t26 - t25
      t33 = t1 ** 2
      t34 = t4 ** 2
      t35 = t12 ** 2
      t36 = t2 ** 2
      t37 = t9 ** 2
      t38 = t5 ** 2
      t39 = t28 * t6
      t40 = t2 * t15
      t41 = t40 * t14
      t42 = t12 * t14
      t40 = t40 * t7
      t43 = t4 * t7
      t44 = t43 * t2
      t45 = t9 * t5 * t16
      t46 = t13 * t17
      t47 = t46 * t25
      t48 = t47 + t45
      t49 = t13 ** 2
      t50 = t5 * t19
      t51 = t13 * t26
      t52 = t1 * t18
      t53 = t51 + t50
      t43 = t43 * t1
      t54 = t4 * t14 * t36
      t55 = t29 * t15
      t55 = t54 * t3 * t1 * t8 * (t51 * (t7 * (t39 - t55) - t52) + t50 *
     # (t7 * (t24 * t6 - t55) - t52))
      t25 = t42 * t34 * t7 * t3 * t6 * t36 * (-t47 * t28 * t13 * t7 * t2
     #7 - t45 * t24 * t5 * t7 * t27 + t43 * (t49 * t17 ** 2 * t25 + t38 
     #* t37 * t16) * t3 * t6)
      ret = -t25 / 4 - t44 * (-t33 * t12 * t6 * t3 * t48 + (t52 * (t22 *
     # (t19 * t38 + t26 * t49) - t50 * (t10 + t11) - t51 * (t21 + t23)) 
     #+ t20 * t4 * t12 * (t46 * t32 - t45)) * t14 * t36) - t44 * (t38 * 
     #(t24 * t22 * t14 * t19 * t2 * (-t2 * t6 * t7 + (-t40 - t1) * t3 * 
     #t18) + t41 * t6 * t3 ** 2 * t7 * t35 * t34 * (t19 - t16) * t37) + 
     #t5 * (t1 * t19 * t14 * t6 * t3 * t35 * t4 * t37 + t19 * (t20 * t36
     # * t14 * t12 * t4 - t30 * t22 * t1 * t12 + t6 * t31 * t3 * t33) * 
     #t9) + t3 * t17 * t13 * (t26 * t6 * t31 * t33 + t41 * t4 * t12 * t1
     #3 * (t39 * t3 * t32 - t2 * t26 * (t22 * t18 + t6)) * t7 + t42 * t2
     #6 * (t29 * (-t13 * t3 * t4 - 1) * t15 + t39) * t1)) + 2 * t54 * (t
     #2 * (t1 * t53 + t40 * t53) + t43 * t48 * t12 * t3 * t6 * t8) - 4 *
     # t55

      hjetmass_box_ppmm_0_0_0_s123_s12_s23_dp = ret/32d0/(0,1d0)
      return

      end function
