
      double complex function hjetmass_triangle_pppm_s34_0_0_dp 
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
      t2 = za(i1, i2)
      t3 = za(i1, i3)
      t4 = zb(i2, i1)
      t5 = zb(i4, i1)
      t6 = zb(i4, i2)
      t7 = za(i2, i3)
      t8 = t3 * t5
      t9 = t6 * t7 + t8
      t10 = zb(i4, i3)
      t11 = za(i1, i4)
      t12 = za(i2, i4)
      t13 = zb(i3, i1)
      t14 = zb(i3, i2)
      t15 = t14 * t5
      t16 = t13 * t6
      t17 = t15 - t16
      t18 = t12 ** 2
      t19 = t6 ** 2
      t20 = t6 * t19
      t21 = t7 ** 2
      t22 = t5 ** 2
      t23 = t3 ** 2
      t24 = t11 ** 2
      t25 = t3 * t7
      t26 = t12 * t6
      t27 = t4 * t10
      t28 = -3 * t27 + t17
      t29 = -t28
      t30 = t11 * t5
      t31 = t30 * t29
      t32 = 4
      t33 = t1 * t10
      t34 = t33 * t32
      t35 = t15 - t16 + t27
      t36 = t26 * t32
      t37 = t36 * t4
      t38 = t30 * t35
      t39 = t1 * t35
      t40 = t26 * t35
      t41 = t34 * t5
      t42 = t7 * t35
      t43 = t7 * t14
      t44 = 7 * t27
      t45 = t16 * t27 * t32
      t46 = t11 * t22
      t47 = t33 * t35
      t34 = t2 * t4 * t6 * (t3 * (t25 * t5 * t6 * (t26 * t28 + t31 + t34
     # * (-t16 + 2 * t15)) + t7 * t19 * (t7 * (t10 * (t30 * t32 * t4 + t
     #39) + t40) + t41 * (3 * t30 + 2 * t26)) - t23 * t22 * (t10 * (t1 *
     # (-5 * t16 + t15 + t27) + t37) + t38)) + t34 * t30 * t21 * t20)
      t48 = t2 ** 2
      t8 = t8 * t19 * t10 * (t25 * (3 * t26 * t24 * t22 + t11 * (-t21 * 
     #t14 * t17 + 3 * t18 * t5 * t19 + t22 * t24 * t5) + t12 * (t23 * t1
     #3 * t17 + t18 * t20) - t25 * (t11 * t13 - t12 * t14) * t17) - t2 *
     # t48 * t1 * t4 ** 2 * t9)
      t4 = -t32 * t8 - t2 * (t23 * t6 * (t3 * (t7 * (t26 * (-t15 * (-t44
     # - t16 + t15) - t45) + t33 * (3 * t15 * (t27 - 2 * t16 + t15) + t1
     #6 * (t27 + 3 * t16))) + t30 * (t6 * (t12 * t5 * t35 + t7 * t13 * t
     #28) + t46 * t35)) - t6 * (t11 * (t7 * (t5 * (2 * t47 + t26 * (-13 
     #* t27 - t17)) + t7 * (-t45 + t15 * (t15 + 9 * t27 - t16))) - t46 *
     # (t41 + 3 * t42)) + t7 * (2 * t43 * (t33 * (t15 - t16 - t27) + t40
     #) + t26 * (t26 * (-t44 + t17) + t47)))) * t5 + t23 ** 2 * t22 * (t
     #10 * (t39 * (t15 + t16) + t32 * t12 * t4 * t13 * t19) + t38 * t16)
     # - t11 * t21 * t19 ** 2 * (t31 * t7 + t33 * (t36 * t5 + t42)) - t2
     #5 * t20 * (t30 * (t7 * (t10 * (3 * t39 + t37) + t42 * t14) + 12 * 
     #t33 * t26 * t5) + t26 * (t42 * (t26 + t33 + t43) + t41 * t26) + 3 
     #* t24 * t7 * t22 * t29) - t34)
      t5 = 0.1D1 / t5
      t8 = 0.1D1 / t9
      t9 = 0.1D1 / t12
      t7 = 0.1D1 / t7
      t2 = 0.1D1 / t2
      t10 = 0.1D1 / t11
      t11 = t8 ** 2
      t6 = 0.1D1 / t6 ** 2
      t3 = 0.1D1 / t3 ** 2
      ret = -t32 * t1 ** 2 * t4 * t2 * t3 * t10 * t7 * t9 * t5 * t6 * t8
     # * t11

      hjetmass_triangle_pppm_s34_0_0_dp = ret/32d0/(0,1d0)
      return

      end function
