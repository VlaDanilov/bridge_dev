
      double complex function hjetmass_box_ppmm_0_0_s23_mhsq_s14_s123_dp 
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
      t3 = za(i2, i4)
      t4 = zb(i4, i2)
      t5 = za(i1, i2)
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = za(i1, i3)
      t9 = zb(i3, i1)
      t10 = za(i2, i3)
      t11 = zb(i3, i2)
      t12 = t5 * t2
      t13 = t8 * t9
      t14 = t10 * t11
      t15 = t6 * t7
      t16 = t15 * (t14 + t12 + t13)
      t9 = t1 * t9
      t17 = t3 * t2
      t18 = t17 + t9
      t19 = t15 * t18
      t20 = t4 * t5 + t8 * zb(i4, i3)
      t21 = mt ** 2
      t22 = 16 * t19 * t21 * t20 + 4 * t16 ** 2
      t22 = cdsqrt(t22)
      t16 = 2 * t16
      t23 = -t22 - t16
      t19 = 0.1D1 / t19
      t9 = t17 + t9
      t9 = t15 * t19 * (t12 * t9 + t13 * t9 + t14 * t9)
      t24 = (0.1D1 / 0.2D1)
      t18 = t24 * t22 * t19 * t18
      t12 = 2 * t14 + 2 * t12 + 2 * t13
      t13 = -t9 + t12 - t18
      t14 = 0.1D1 / t20
      t20 = t4 * t13 * t14 * t24
      t17 = t17 / 4
      t25 = t5 * (t2 - t20) + t17 * t23 * t19
      t16 = t22 - t16
      t9 = -t9 + t12 + t18
      t12 = t4 * t9 * t14 * t24
      t5 = t5 * (t2 - t12) + t17 * t16 * t19
      t17 = t20 * t8 - t23 * t19 * t1 * t2 / 4
      t12 = t12 * t8 - t16 * t19 * t1 * t2 / 4
      t18 = 0.1D1 / t25
      t20 = 0.1D1 / t10
      t3 = 0.1D1 / t3
      t11 = 0.1D1 / t11
      t22 = 0.1D1 / t4
      t5 = 0.1D1 / t5
      t25 = t12 * t9
      t26 = t13 * t17
      t27 = t13 + t9
      t28 = t1 ** 2
      t29 = t6 ** 2
      t30 = (t26 + t25) * t22
      t31 = t21 * t28
      t32 = t3 * t2
      t33 = t13 ** 2
      t34 = t13 * t33
      t35 = t14 ** 2
      t36 = t9 ** 2
      t37 = t9 * t36
      t38 = t2 ** 2
      t39 = t36 * t5
      t40 = t33 * t18
      t41 = t35 * t29
      t42 = t40 + t39
      t43 = t9 * t5
      t44 = t13 * t18
      t45 = t44 + t43
      t46 = t36 + t33
      t47 = t1 * t20
      t48 = t6 * t42
      t49 = t37 * t5
      t50 = t34 * t18
      t51 = t27 * t1
      t52 = t6 * t11
      t53 = t52 * t38
      t54 = t14 * t3
      t9 = t54 * (t14 * t29 * t10 * t38 * t4 * t11 * t42 + (t14 * (t46 *
     # t20 * t8 + t48 * t4) - t53 * t20 * (t13 * t23 + t16 * t9) * t19 *
     # t22) * t7 * t28 + t15 * (t8 * (t2 * t42 * t14 + t4 * t20 * t11 * 
     #(t37 + t34) * t35) - t52 * t2 * t38 * t22 * t19 * (t43 * t16 + t44
     # * t23)) * t1)
      t13 = 0.1D1 / t6
      ret = -t24 * t41 * t1 * t38 * t7 * t19 * t3 * t11 * (t39 * t16 + t
     #40 * t23) - (0.3D1 / 0.2D1) * t32 * t14 * t35 * t8 * t29 * t7 * t4
     # * t11 * (t50 + t49) - 5 * t41 * t32 * t7 * t11 * (t39 * t12 + t40
     # * t17) + 4 * t32 * (t31 * t2 * t22 * (t18 + t5) + (t30 * t15 * t1
     # * t20 + t2 * (t26 * t18 + t25 * t5) * t22 * t7 * t29 + t31 * t20 
     #* t27) * t14 * t11) + 2 * t54 * (t7 * (t1 * ((t44 * (t17 * t22 + t
     #1) + t43 * (t12 * t22 + t1)) * t6 * t2 + t47 * (t30 + t51)) + t11 
     #* (t41 * t4 * (t49 * t12 + t50 * t17) + t6 * (t48 * t8 * t38 - t47
     # * (t12 * t36 + t17 * t33) + t47 * t46 * t8 * t2) * t14)) + t1 * t
     #38 * t6 * t11 * t45 * t21 + t53 * (t45 * t10 * t6 * t2 + t51)) - t
     #9 + 8 * t32 * t21 * t1 * t28 * t13 * t20 * t22

      hjetmass_box_ppmm_0_0_s23_mhsq_s14_s123_dp = ret/32d0/(0,1d0)
      return

      end function
