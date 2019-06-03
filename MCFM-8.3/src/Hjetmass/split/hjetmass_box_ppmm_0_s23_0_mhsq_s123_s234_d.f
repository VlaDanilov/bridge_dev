
      double complex function hjetmass_box_ppmm_0_s23_0_mhsq_s123_s234_dp 
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
      t2 = za(i1, i4)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = zb(i4, i1)
      t6 = za(i1, i2)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = za(i1, i3)
      t10 = zb(i3, i1)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t7 * t8
      t14 = t11 * t12
      t15 = t14 + t13
      t16 = t2 * t3 * t4 * t5
      t17 = t9 * t10
      t18 = t2 * t5
      t15 = t18 * (t15 * t6 * t1 + t17 * t15 - t16)
      t19 = t7 * t1
      t20 = t11 * t10
      t21 = t19 + t20
      t22 = mt ** 2
      t23 = t6 * t8
      t24 = t9 * t12
      t25 = t24 + t23
      t26 = t2 ** 2 * t5 ** 2
      t27 = t26 * t25
      t28 = 16 * t18 * t22 * t21 * t27 + 4 * t15 ** 2
      t28 = cdsqrt(t28)
      t15 = 2 * t15
      t29 = -t28 + t15
      t30 = t14 + t13
      t31 = 0.1D1 / t5
      t27 = 0.1D1 / t27
      t31 = t30 * t31
      t32 = t31 * t1
      t33 = t29 * t27 / 4
      t34 = -t33 * t2 * t8 + t32
      t24 = t24 + t23
      t35 = -t6 * t24 * t1 - t17 * t24
      t17 = -t30 * t6 * t1 - t17 * t30 + t16
      t25 = t18 * t28 * t27 * t25 / 4
      t13 = t26 * t27 * (t13 * t35 + t14 * t35 + t16 * t24) / 2
      t14 = t18 * t21
      t14 = 0.1D1 / t14
      t16 = 0.1D1 / t2
      t16 = t31 * t16
      t18 = t16 * t9
      t21 = (-t17 + t25 + t13) * t14
      t24 = t21 * t11
      t26 = -t24 + t18
      t30 = t33 * t9
      t31 = t1 * t26 - t30 * t8
      t6 = t16 * t6
      t16 = -t23 * t33 + t1 * (-t21 * t7 + t6)
      t15 = t28 + t15
      t21 = t15 * t27 / 4
      t28 = -t21 * t2 * t8 + t32
      t13 = (-t17 - t25 + t13) * t14
      t14 = t13 * t11
      t17 = -t14 + t18
      t18 = t21 * t9
      t25 = t1 * t17 - t18 * t8
      t6 = t1 * (-t13 * t7 + t6) - t21 * t23
      t7 = 0.1D1 / t8
      t3 = 0.1D1 / t3
      t4 = 0.1D1 / t4
      t8 = 0.1D1 / t9
      t13 = 0.1D1 / t16
      t6 = 0.1D1 / t6
      t16 = t28 * t4 + t11
      t21 = t25 ** 2
      t23 = t2 * t31
      t32 = t7 * t1
      t33 = t11 * t8
      t35 = t13 * t31
      t5 = t27 * t3 * t5
      t27 = t31 ** 2
      t22 = t22 * t1 * t11 * (t32 - t33)
      t31 = t13 * t31 * t27
      t36 = t15 * t21 * t6
      ret = 3 * t5 * t11 * t7 * (t29 * t27 * t13 + t36) + t5 * (t15 * (t
     #2 * t4 * (-t32 + t33) * t6 * t21 + (t32 * t16 * t9 - t11 * t16) * 
     #t6 * t25) + t35 * t29 * (t11 * (t4 * (t23 * t8 - t34) - t11) + t32
     # * (t11 * t9 + t4 * (t34 * t9 - t23)))) + 4 * t4 * t3 * (-t31 * t2
     #0 * t8 * t7 + t6 * t25 * (-t20 * t21 * t8 * t7 + t22 + t32 * (t11 
     #* ((t10 * t17 - t18 * t12) * t8 + t10) + t19) * t25) + t32 * (t11 
     #* ((t10 * t26 - t30 * t12) * t8 + t10) + t19) * t13 * t27 + t22 * 
     #t35) + 2 * t5 * t4 * t7 * (t29 * t34 * t27 * t13 + (t29 * (t24 * t
     #1 * t13 * t27 + t31) + t36 * (t14 * t1 + t25)) * t8 * t2 + t15 * t
     #28 * t21 * t6)

      hjetmass_box_ppmm_0_s23_0_mhsq_s123_s234_dp = ret/32d0/(0,1d0)
      return

      end function
