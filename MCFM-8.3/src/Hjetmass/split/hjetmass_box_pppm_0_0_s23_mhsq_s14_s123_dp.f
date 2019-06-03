
      double complex function hjetmass_box_pppm_0_0_s23_mhsq_s14_s123_dp 
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
      t2 = za(i1, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i1)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = t1 * t3
      t10 = t5 * t6
      t11 = t7 * t8
      t12 = t2 * t4
      t13 = t12 * (t11 + t9 + t10)
      t14 = za(i2, i4)
      t15 = za(i3, i4)
      t16 = t14 * t3
      t17 = t15 * t6
      t18 = t17 + t16
      t19 = t12 * t18
      t20 = zb(i4, i2)
      t21 = zb(i4, i3)
      t22 = t1 * t20
      t23 = t21 * t5 + t22
      t24 = mt ** 2
      t25 = 16 * t19 * t24 * t23 + 4 * t13 ** 2
      t25 = cdsqrt(t25)
      t13 = 2 * t13
      t26 = t25 - t13
      t19 = 0.1D1 / t19
      t16 = t17 + t16
      t16 = t12 * t19 * (t10 * t16 + t11 * t16 + t9 * t16)
      t17 = (0.1D1 / 0.2D1)
      t9 = 2 * t11 + 2 * t9 + 2 * t10
      t10 = t17 * t25 * t19 * t18
      t11 = -t16 + t10 + t9
      t18 = 0.1D1 / t23
      t23 = (0.1D1 / 0.4D1)
      t22 = t22 * t17
      t27 = t3 * (t14 * t26 * t19 * t23 + t1) - t22 * t11 * t18
      t13 = -t25 - t13
      t9 = -t16 - t10 + t9
      t10 = t3 * (t14 * t13 * t19 * t23 + t1) - t22 * t9 * t18
      t16 = 0.1D1 / t5
      t10 = 0.1D1 / t10
      t7 = 0.1D1 / t7
      t22 = 0.1D1 / t27
      t25 = 0.1D1 / t1
      t27 = t6 ** 2
      t28 = t3 * t6 * t16
      t29 = t25 * (t4 * t8 * t14 * t16 + t27) + t28
      t30 = t11 ** 2
      t31 = t2 ** 2
      t32 = t6 * t7
      t33 = t32 * t25
      t34 = t22 * t11
      t35 = t9 * t13
      t21 = t3 * t21
      t36 = t6 * t20
      t37 = -t36 + t21
      t36 = t36 - t21
      t38 = t9 ** 2
      t39 = t26 * t30 * t22
      t40 = t25 * t18 ** 2
      t20 = t20 * (t12 * t16 - t6) + t21
      t21 = t24 * t25
      t41 = t18 * t2
      t42 = t11 + t9
      t43 = t9 * t10
      t44 = t43 + t34
      t45 = -t21 + t3
      t46 = t42 * t7
      t1 = t1 * t3 ** 2 * t16
      t1 = t41 * (t8 * (t4 * (t16 * (t2 * (t34 * t45 + t43 * t45) + t46 
     #* t15) - t46 * t25 * t14) + t6 * (t21 * t44 - t3 * t44) - t1 * t44
     #) - t46 * t27 * t25 * t5 - t1 * t46)
      ret = t17 * t18 * t19 * t31 * (t26 * (t34 * t3 * t29 + t18 * t4 * 
     #((t22 * t8 + t7) * t16 * t3 + t33) * t30) + t35 * (t10 * t3 * t29 
     #+ t4 * ((t10 * t8 + t7) * t16 * t3 + t33) * t18 * t9)) + t23 * t40
     # * t19 * t31 * (t6 * (t37 * t10 * t38 * t13 + t39 * t37) + t12 * (
     #t36 * t10 * t38 * t13 + t39 * t36) * t16) - t40 * t28 * t19 ** 2 *
     # t2 * t31 * t4 * (t13 ** 2 * t38 * t10 + t26 ** 2 * t30 * t22) / 8
     # - t41 * (t21 * t28 * t19 * t2 * (t35 * t10 + t34 * t26) + (t20 * 
     #t10 * t38 + t30 * t22 * t20) * t18 * t8) - 4 * t41 * (t3 * (t32 * 
     #t42 + (t7 * (-t11 - t9) + t8 * (-t43 - t34)) * t16 * t24) - t33 * 
     #t24 * t42) + 2 * t1

      hjetmass_box_pppm_0_0_s23_mhsq_s14_s123_dp = ret/32d0/(0,1d0)
      return

      end function
