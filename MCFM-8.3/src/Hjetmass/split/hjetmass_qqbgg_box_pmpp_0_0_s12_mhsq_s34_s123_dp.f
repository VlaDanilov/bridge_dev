
      double complex function hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s123_dp
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = zb(i3, i1)
      t2 = za(i1, i2)
      t3 = zb(i2, i1)
      t4 = za(i1, i3)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i1)
      t9 = zb(i4, i3)
      t10 = t4 * t8
      t11 = t5 * zb(i4, i2)
      t12 = t11 + t10
      t13 = t7 * t9
      t14 = t13 * t12
      t15 = t2 * t3
      t16 = t4 * t1
      t17 = t5 * t6
      t18 = t13 * (t17 + t15 + t16)
      t19 = za(i1, i4)
      t20 = za(i2, i4)
      t21 = t19 * t1
      t6 = t20 * t6 + t21
      t22 = mt ** 2
      t23 = 16 * t22 * t6 * t14 + 4 * t18 ** 2
      t23 = cdsqrt(t23)
      t14 = 0.1D1 / t14
      t24 = t17 + t16
      t11 = t13 * t14 * (t10 * t24 + t11 * t24 + t15 * (t11 + t10))
      t13 = 2
      t24 = (0.1D1 / 0.2D1)
      t12 = t24 * t23 * t14 * t12
      t15 = t13 * (t17 + t15 + t16)
      t17 = -t12 - t11 + t15
      t18 = t13 * t18
      t25 = -t23 - t18
      t6 = 0.1D1 / t6
      t26 = -(0.1D1 / 0.4D1)
      t10 = t10 * t26
      t27 = t21 * t24
      t28 = t10 * t25 * t14 + t27 * t17 * t6
      t11 = t12 - t11 + t15
      t12 = t23 - t18
      t10 = t10 * t12 * t14 + t27 * t11 * t6
      t15 = t19 * t5 - t20 * t4
      t18 = t12 * t14 * t11 * t6 * t9 * t15
      t15 = t25 * t14 * t17 * t6 * t9 * t15
      t4 = 0.1D1 / t4
      t23 = 0.1D1 / t2
      t27 = 0.1D1 / t19
      t3 = 0.1D1 / t3
      t29 = t4 * (t28 + t10)
      t30 = t3 * t9
      t31 = t8 * t3
      t32 = t1 * t3
      t33 = t12 * t11
      t34 = t17 + t11
      t35 = t4 * t1
      t36 = t30 * t6
      t37 = t35 * t36 * t23 * (t11 * t18 + t15 * t17)
      t38 = -t27 * (t32 * t7 + t20) + (-t31 * t19 * t23 + 1) * t4 * t5
      t3 = t9 * (t3 * t27 * t23 * (t1 * (t18 + t15) + t4 * (-t10 * t18 -
     # t15 * t28)) + (t25 * t17 * t38 + t33 * t38) * t14 * t6 * t9)
      t7 = t31 * t14 ** 2 * t9 ** 2 * t6 * t7 * t27 * (t11 * t12 ** 2 + 
     #t17 * t25 ** 2)
      ret = -t13 * t36 * ((t20 * (t1 * t34 + t4 * (-t10 * t11 - t17 * t2
     #8)) - t21 * t4 * t34 * t5) * t23 * t9 + t35 * (t11 * (-t22 + t10) 
     #+ t17 * (-t22 + t28) + t34 * t19 * t8)) - t24 * t3 + t26 * t37 - 8
     # * t9 * ((t2 * t27 - t31) * t4 * t22 + t32 * (-t9 * t5 * t23 + t8)
     #) + t7 / 8 + 4 * t30 * (t27 * (t1 * (-t28 - t10) + t4 * (t10 ** 2 
     #+ t10 * t22 + t22 * t28 + t28 ** 2)) + t29 * t8 - t29 * t23 * t9 *
     # t5) + t30 * (t9 * (t14 * ((t5 * (t28 * t4 - t1) + (-t16 + t28) * 
     #t27 * t20) * t25 * t17 + t33 * (t5 * (t10 * t4 - t1) + (-t16 + t10
     #) * t27 * t20)) * t23 * t6 - t21 * t20 * t4 * (t11 ** 2 + t17 ** 2
     #) * t23 * t6 ** 2) + t27 * t14 * t8 * (t12 * (t22 + t10) + t25 * (
     #t22 + t28) + t16 * (-t25 - t12)))

      hjetmass_qqbgg_box_pmpp_0_0_s12_mhsq_s34_s123_dp
     & = ret/32d0/(0,1d0)
      return

      end function
