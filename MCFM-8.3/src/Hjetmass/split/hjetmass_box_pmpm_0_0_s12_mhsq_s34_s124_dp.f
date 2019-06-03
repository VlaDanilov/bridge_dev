
      double complex function hjetmass_box_pmpm_0_0_s12_mhsq_s34_s124_dp 
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
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = zb(i3, i1)
      t10 = t5 * t9
      t11 = t7 * zb(i3, i2) + t10
      t12 = mt ** 2
      t13 = za(i1, i3)
      t14 = za(i2, i3)
      t15 = t13 * t6
      t16 = t14 * t8
      t17 = t16 + t15
      t18 = t2 * t4
      t19 = t18 * t17
      t20 = t1 * t3
      t21 = t5 * t6
      t8 = t7 * t8
      t22 = t18 * (t8 + t20 + t21)
      t23 = 16 * t12 * t11 * t19 + 4 * t22 ** 2
      t23 = cdsqrt(t23)
      t22 = 2 * t22
      t24 = -t23 - t22
      t19 = 0.1D1 / t19
      t15 = t16 + t15
      t15 = t18 * t19 * (t20 * t15 + t21 * t15 + t8 * t15)
      t16 = (0.1D1 / 0.2D1)
      t17 = t16 * t23 * t19 * t17
      t8 = 2 * t8 + 2 * t20 + 2 * t21
      t20 = -t15 + t8 - t17
      t11 = 0.1D1 / t11
      t21 = (0.1D1 / 0.4D1)
      t10 = t10 * t16
      t25 = t6 * (t13 * t24 * t19 * t21 + t5) - t10 * t20 * t11
      t26 = t16 * t20 * t11 * t7 * t9 - t24 * t21 * t19 * t14 * t6
      t22 = t23 - t22
      t8 = -t15 + t8 + t17
      t5 = t6 * (t13 * t22 * t19 * t21 + t5) - t10 * t8 * t11
      t10 = t16 * t8 * t11 * t7 * t9 - t22 * t21 * t19 * t14 * t6
      t13 = 0.1D1 / t14
      t1 = 0.1D1 / t1
      t5 = 0.1D1 / t5
      t3 = 0.1D1 / t3
      t15 = 0.1D1 / t25
      t17 = t6 * t7
      t23 = -t17 + t10
      t25 = (-t17 + t26) * t15
      t27 = t22 * t5
      t28 = t27 * t23
      t29 = t25 * t24 + t28
      t30 = t7 ** 2
      t31 = t22 ** 2
      t32 = t26 * t13
      t33 = t10 * t13
      t34 = t10 * t5
      t35 = t26 ** 2
      t36 = t26 * t15
      t37 = t24 ** 2
      t38 = t19 ** 2
      t39 = t31 * t5
      t40 = t38 * t3
      t41 = 0.1D1 / t4
      t42 = t6 ** 2 * t30
      ret = t16 * t40 * t18 * t1 * (t39 * (t10 ** 2 + t42) + (t42 + t35)
     # * t15 * t37) - t21 * t40 * t9 * t6 * t2 * (t39 * t23 + t25 * t37)
     # - t18 * t19 * t38 * t3 * t1 * t14 * t6 * (t22 * t31 * t5 * t23 + 
     #t25 * t24 * t37) / 8 - t19 * (t30 * t1 * t4 * t29 + (t27 * t9 * (t
     #10 * (t33 + t9) + t17 * (-t33 - t9) + t30 * t13 * t23 * t11 * t1 *
     # t8 * t4) + t15 * t24 * (t9 * (t26 * (t32 + t9) + t17 * (-t32 - t9
     #)) + t7 * (-t6 * t30 * t9 * t20 * t13 * t11 + t32 * t7 * t9 * t20 
     #* t11 + t24 * t6 * t26 * t19) * t1 * t4) + t34 * t17 * t4 * t31 * 
     #t1 * t19) * t3 * t2) + 2 * t19 * t3 * t1 * t7 * (t9 * t12 * t29 + 
     #t18 * (t24 * (t15 * t35 - t36 * t17) + t28 * t10) * t13) - 4 * t41
     # * t3 * t13 * t9 ** 2 * t7 * t12 * (t17 * (t5 + t15) - t34 - t36)

      hjetmass_box_pmpm_0_0_s12_mhsq_s34_s124_dp = ret/32d0/(0,1d0)
      return

      end function
