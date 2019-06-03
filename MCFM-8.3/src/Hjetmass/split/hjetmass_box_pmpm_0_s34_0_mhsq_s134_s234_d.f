
      double complex function hjetmass_box_pmpm_0_s34_0_mhsq_s134_s234_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = zb(i4, i2)
      t6 = t3 * t4
      t7 = t1 * t5
      t8 = t7 + t6
      t9 = zb(i2, i1)
      t10 = za(i1, i2)
      t11 = za(i1, i3)
      t12 = za(i1, i4)
      t13 = zb(i4, i1)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t1 * t13
      t17 = t2 * t3 + t16
      t18 = mt ** 2
      t19 = t11 * t4
      t20 = t12 * t5
      t21 = t19 + t20
      t22 = t10 ** 2 * t9 ** 2
      t23 = t22 * t21
      t6 = t7 + t6
      t7 = t10 * t14 * t9 * t15
      t24 = t12 * t13
      t25 = t6 * t11
      t26 = t10 * t9
      t27 = t26 * (t25 * t2 + t24 * t6 - t7)
      t28 = 16 * t18 * t10 * t9 * t17 * t23 + 4 * t27 ** 2
      t28 = cdsqrt(t28)
      t27 = 2 * t27
      t29 = t27 - t28
      t23 = 0.1D1 / t23
      t30 = 0.1D1 / t9
      t30 = t8 * t30
      t31 = t30 * t2
      t32 = t29 * t23 / 4
      t33 = t32 * t10
      t34 = -t33 * t4 + t31
      t6 = t20 * t6 + t25 * t4
      t8 = t8 * t11 * t2 + t24 * t8 - t7
      t21 = t26 * t28 * t23 * t21 / 4
      t6 = t22 * t23 * (t6 * t11 * t2 + t24 * t6 + t7 * (-t20 - t19)) / 
     #2
      t7 = -t8 - t21 + t6
      t11 = t26 * t17
      t17 = 0.1D1 / t10
      t11 = 0.1D1 / t11
      t19 = t30 * t17 * t12
      t22 = t7 * t11 * t1 + t19
      t24 = t13 * t22 - t20 * t32
      t22 = -t32 * t12 * t4 + t2 * t22
      t25 = t30 * t13
      t26 = t27 + t28
      t27 = t26 * t23 / 4
      t28 = t27 * t10
      t30 = -t28 * t4 + t31
      t6 = -t8 + t21 + t6
      t8 = t6 * t11 * t1 + t19
      t19 = t13 * t8 - t27 * t20
      t8 = -t27 * t12 * t4 + t2 * t8
      t4 = 0.1D1 / t4
      t20 = 0.1D1 / t15
      t3 = 0.1D1 / t3
      t21 = 0.1D1 / t24
      t19 = 0.1D1 / t19
      t14 = 0.1D1 / t14
      t24 = t2 * (-t33 * t5 + t25)
      t27 = t13 * t34
      t5 = t2 * (-t28 * t5 + t25)
      t25 = t30 * t13
      t28 = t21 + t19
      t31 = t2 ** 2
      t32 = t1 ** 2
      t13 = t2 * t13
      t33 = (t25 - t5) * t17
      t35 = t1 * t14
      t36 = t30 * t19
      t37 = (t27 - t24) * t17
      t38 = t34 * t21
      t39 = t22 * t21
      t40 = t2 * t1
      t41 = t3 * t4
      t42 = t41 * t1
      t43 = t29 + t26
      t44 = t34 ** 2
      t45 = t26 * t30
      t46 = t29 * t34
      t47 = t44 * t21
      t48 = t30 * t12
      t49 = t26 * t19
      t50 = t29 * t21
      t51 = t1 * t15
      t7 = t29 * t7
      t52 = t3 * t23 * t9
      t6 = t52 * t32 * ((-t12 * t43 + (-t26 * t6 - t7) * t11 * t1) * t14
     # * t4 * t2 + t49 * (-t51 + t30) + t50 * (-t51 + t34)) + t52 * (t1 
     #* t32 * t14 * t43 + (t1 * (t10 * (t39 * t29 + t49 * t8) + t51 * (t
     #50 + t49) * t12) + t20 * (t29 * (-t38 * t10 * t22 + t47 * t12) + t
     #36 * t26 * (-t10 * t8 + t48) + t35 * (t10 * (-t22 * t29 - t26 * t8
     #) + t12 * (t46 + t45))) + t32 * (t7 * t21 * (t51 - t34) + t49 * (t
     #51 - t30) * t6) * t11) * t4 * t2)
      t7 = t41 * t23
      t10 = t30 ** 2
      t11 = t2 * t12
      t26 = t12 * t20 * t17
      t11 = t7 * t9 * (t29 * (-t26 * t21 * t34 * t44 - t20 * (t35 * t12 
     #* t17 + t39) * t44 + t1 * (t21 * (t11 + t22) - t22 * t20 * t14) * 
     #t34) + t45 * (t19 * (t1 * (t11 + t8) - t26 * t10 - t8 * t30 * t20)
     # - t35 * t20 * (t48 * t17 + t8)))
      ret = 3 * t7 * t17 * t12 * t9 * t32 * (t46 * (-t15 * t21 + t14) + 
     #t45 * (-t15 * t19 + t14)) + 5 * t42 * t23 * t12 * t9 * t17 * (t49 
     #* t10 + t47 * t29) - 8 * t41 * t40 * t18 * (t20 * (t17 * (t10 * t1
     #9 + t47) + t35 * (t17 * (t34 + t30) + t2)) + t17 * t32 * t28 * t15
     #) + 4 * t42 * (t1 * (t39 * (t17 * (-t27 + t24) + t13) - t13 * (t22
     # + t8) * t14 * t20 + t8 * t17 * t19 * (-t25 + t5)) - t31 * t20 * (
     #t36 + t38) * t18) + 4 * t4 * (t20 * (t22 * (t38 * ((t37 - t13) * t
     #3 * t1 - t31) + t35 * (t37 * t3 * t1 - t31)) + t8 * (t36 * ((t33 -
     # t13) * t3 * t1 - t31) + t35 * (t33 * t3 * t1 - t31))) + t32 * t31
     # * t3 * t28 * t18 + t40 * (t39 * t2 + (t16 * t3 + t2) * t19 * t8))
     # - t6 + 2 * t11 + 16 * t41 * t17 * t2 * t32 * t18 * (t35 + t38 + t
     #36)

      hjetmass_box_pmpm_0_s34_0_mhsq_s134_s234_dp = ret/32d0/(0,1d0)
      return

      end function
