
      double complex function hjetmass_box_ppmm_0_0_s14_mhsq_s23_s134_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i3)
      t2 = za(i2, i3)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = zb(i2, i1)
      t10 = zb(i4, i2)
      t11 = t7 * t10
      t12 = t1 * t9 + t11
      t13 = mt ** 2
      t14 = za(i2, i4)
      t15 = za(i1, i2) * t3
      t16 = t14 * t8
      t17 = t16 + t15
      t18 = t2 * t4
      t19 = t18 * t17
      t1 = t1 * t3
      t20 = t5 * t6
      t21 = t7 * t8
      t22 = t18 * (t21 + t1 + t20)
      t23 = 16 * t13 * t12 * t19 + 4 * t22 ** 2
      t23 = cdsqrt(t23)
      t22 = 2 * t22
      t24 = t23 - t22
      t19 = 0.1D1 / t19
      t15 = t16 + t15
      t15 = t18 * t19 * (t1 * t15 + t20 * t15 + t21 * t15)
      t16 = (0.1D1 / 0.2D1)
      t1 = 2 * t21 + 2 * t1 + 2 * t20
      t17 = t16 * t23 * t19 * t17
      t20 = -t15 + t1 + t17
      t12 = 0.1D1 / t12
      t21 = (0.1D1 / 0.4D1)
      t25 = t21 * t24 * t19 * t14
      t26 = t16 * t20 * t12 * t7
      t27 = t25 * t3 - t26 * t9
      t22 = -t23 - t22
      t1 = -t15 + t1 - t17
      t15 = t14 * t22 * t19 * t21
      t11 = t8 * (t7 + t15) - t11 * t16 * t1 * t12
      t8 = t8 * (t7 + t25) - t26 * t10
      t15 = -t16 * t1 * t12 * t7 * t9 + t15 * t3
      t8 = 0.1D1 / t8
      t11 = 0.1D1 / t11
      t17 = 0.1D1 / t6
      t10 = 0.1D1 / t10
      t5 = 0.1D1 / t5
      t23 = 0.1D1 / t14
      t25 = t24 + t22
      t26 = t7 ** 2
      t28 = t4 ** 2
      t29 = t9 ** 2
      t30 = t9 * t2
      t31 = t24 * t8
      t32 = t22 * t11
      t33 = t5 * t9
      t34 = t24 * t27
      t35 = t22 * t15
      t36 = t17 * t29
      t37 = t36 * t2
      t38 = t10 * t19
      t39 = t9 * t17
      t40 = t24 ** 2
      t41 = t24 * t40
      t42 = t22 ** 2
      t43 = t22 * t42
      t44 = t42 * t11
      t45 = t40 * t8
      t46 = t45 + t44
      t47 = t42 + t40
      t48 = t19 ** 2
      t49 = t19 * t48
      t50 = t43 * t11
      t48 = t10 * t48
      t19 = t48 * (t37 * t3 * t47 + t18 * t9 * (t9 * t14 * t46 + t46 * t
     #7 * t3) + t5 * (-t50 * t2 * t14 * t15 * t19 + t44 * t26 * (-t30 * 
     #t1 * t12 + t14 * t6) + t45 * (t14 * (-t34 * t2 * t19 + t26 * t6) -
     # t30 * t26 * t20 * t12)) * t28)
      t43 = t33 * t18 * t49 * t14 * t3 * t17 * t10 * (t41 + t43)
      t14 = t49 * t2 * t14 * t7 * t3 * t28 * t5 * t10 * (t41 * t8 + t50)
      t41 = 0.1D1 / t4
      ret = t16 * t48 * t18 * t5 * (t26 * t3 * t4 * t46 + t39 * (t47 * t
     #7 * t3 + t15 * t42 + t27 * t40)) - t21 * t19 + (0.5D1 / 0.4D1) * t
     #48 * t2 * t7 * t28 * t5 * (t44 * t15 + t45 * t27) + t43 / 8 + (0.3
     #D1 / 0.16D2) * t14 - t38 * (t4 * (t33 * (t22 * (t11 * t13 + 1) + t
     #24 * (t13 * t8 + 1) + t30 * (-t1 * t22 - t20 * t24) * t23 * t17 * 
     #t12) * t26 + t30 * (t32 * (-t15 * t23 + t9) + t31 * (-t23 * t27 + 
     #t9)) * t7) + t37 * (t23 * (-t35 - t34) + t25 * t9) + t5 * (t6 * (t
     #32 + t31) + t30 * (-t32 * t1 - t31 * t20) * t23 * t12) * t7 * t26 
     #* t28) + 4 * t13 * t26 * t29 * t23 * t10 * (t8 + t11) - 2 * t38 * 
     #t5 * t7 * (t2 * (-t39 * (t35 + t34) * t23 * t4 - t7 * (t32 * t15 +
     # t31 * t27) * t23 * t28) + t36 * t13 * t25) + 8 * t13 * t7 * t9 * 
     #t29 * t23 * t41 * t17 * t10

      hjetmass_box_ppmm_0_0_s14_mhsq_s23_s134_dp = ret/32d0/(0,1d0)
      return

      end function
