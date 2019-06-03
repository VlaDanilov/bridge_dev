
      complex*32 function hjetmass_qqbgg_box_pmpp_0_s12_0_mhsq_s123_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = zb(i3, i1)
      t2 = zb(i4, i2)
      t3 = za(i1, i3)
      t4 = za(i1, i4)
      t5 = za(i3, i4)
      t6 = zb(i4, i1)
      t7 = zb(i4, i3)
      t8 = za(i2, i3)
      t9 = zb(i3, i2)
      t10 = za(i2, i4)
      t11 = za(i1, i2)
      t12 = zb(i2, i1)
      t13 = t4 * t6
      t14 = t10 * t2
      t15 = t5 * t7
      t16 = t15 * t11
      t17 = t16 * t12
      t18 = t8 * t9
      t19 = t15 * ((-t14 - t13) * t3 * t1 + t17 - t18 * (t14 + t13))
      t20 = t1 * t4 + t10 * t9
      t21 = mt ** 2
      t22 = t8 * t2
      t23 = t3 * t6
      t24 = t23 + t22
      t25 = t7 ** 2
      t26 = t5 ** 2 * t25
      t27 = t26 * t24
      t28 = 16 * t15 * t21 * t20 * t27 + 4 * t19 ** 2
      t28 = sqrt(t28)
      t19 = 2 * t19
      t29 = t28 - t19
      t27 = 0.1q1 / t27
      t30 = t23 + t22
      t30 = t3 * t30 * t1 + t18 * t30
      t31 = t14 + t13
      t18 = t31 * t3 * t1 + t18 * t31 - t17
      t24 = t15 * t28 * t27 * t24 / 4
      t17 = t26 * t27 * (t13 * t30 + t14 * t30 + t17 * (-t23 - t22)) / 2
      t22 = -t18 + t24 + t17
      t20 = t15 * t20
      t26 = 0.1q1 / t7
      t20 = 0.1q1 / t20
      t5 = 0.1q1 / t5
      t5 = (t15 + t13 + t14) * t5
      t13 = t5 * t26 * t3
      t14 = t22 * t20
      t15 = t14 * t4
      t26 = t15 + t13
      t30 = t23 / 4
      t31 = -t30 * t29 * t27 + t1 * t26
      t32 = t5 * t8
      t14 = t14 * t10 * t7 + t32
      t17 = -t18 - t24 + t17
      t18 = -t28 - t19
      t19 = t17 * t20
      t24 = t19 * t4
      t13 = t24 + t13
      t28 = -t18 * t30 * t27 + t1 * t13
      t10 = t19 * t10 * t7 + t32
      t5 = t5 * t3
      t19 = t24 * t7 + t5
      t5 = t15 * t7 + t5
      t15 = t9 * t26 - t29 * t27 * t3 * t2 / 4
      t13 = t9 * t13 - t18 * t27 * t3 * t2 / 4
      t12 = 0.1q1 / t12
      t24 = 0.1q1 / t11
      t4 = 0.1q1 / t4
      t26 = 0.1q1 / t3
      t30 = t18 * t19 + t29 * t5
      t32 = t26 * t31 - t1
      t33 = t26 * t28 - t1
      t34 = t31 * t14
      t35 = t8 * t26
      t36 = t4 * t27
      t37 = t36 * t7
      t38 = t31 + t28
      t39 = t38 * t26
      t40 = t6 * t11
      t36 = t36 * t1 * t3
      t41 = t9 * t4
      t42 = t24 * t8
      t43 = t29 + t18
      t44 = t29 * t31
      t45 = t18 * t28
      t46 = t4 * t43
      t16 = t27 * t7 * (t8 * (t12 * (t26 * (t4 * (t2 * (-t45 - t44) + t6
     # * (-t13 * t18 - t15 * t29)) + t6 * t30 * t24) + t46 * t2 * t1) + 
     #t46 * t7) - t16 * t4 * t26 * (t17 * t18 + t22 * t29) * t20 - t46 *
     # t40)
      t3 = t37 * t12 * (t42 * t7 * (t45 + t44) + (t24 * (t10 * t18 + t14
     # * t29) + t43 * t6) * t3 * t1)
      t17 = -t41 * t11 + t7
      t2 = t21 * (-t7 * t11 * t4 * t26 + t12 * t26 * t17 * t6 + t4 * t12
     # * (t26 * (-t11 * t2 + t19 + t5) - t7) * t1) - t1 * t6 * t12 * t17
     # + t1 * t7 * t12 * (t24 * t7 - t41) * t8
      ret = -3 * t3 + 8 * t2 + 2 * t37 * (-t35 * t30 + t12 * (t18 * (t21
     # + t28) + t29 * (t21 + t31)) * t6 + t12 * (t18 * (t33 * t19 * t8 +
     # t10 * t28) + t29 * (t32 * t5 * t8 + t34)) * t24) + 4 * t12 * (t26
     # * (t40 * (t1 * (t15 + t13) - t31 * t9) * t4 - t1 * t24 * (t19 * t
     #10 + t14 * t5)) + t36 * t42 * t25 * t18 + t35 * t4 * (t1 * (-t15 -
     # t13) + t38 * t9) * t7) + 4 * t12 * (t7 * (t24 * (t10 * t33 + t14 
     #* t32 + t35 * t1 * (t19 + t5)) + t39 * t6) + t1 * t26 * (t4 * (t10
     # * t13 + t14 * t15) + t6 * (-t19 - t5)) + t41 * (t1 * (t14 + t10) 
     #+ t26 * (-t28 * (t40 + t10) - t34)) + t42 * (t36 * t29 - t39) * t2
     #5) + t42 * t23 * t27 ** 2 * t25 * t4 * t12 * (t18 ** 2 + t29 ** 2)
     # / 4 + t16

      hjetmass_qqbgg_box_pmpp_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function
