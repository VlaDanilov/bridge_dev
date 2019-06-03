
      complex*32 function hjetmass_box_pmpm_0_0_s34_mhsq_s12_s234
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i2, i3)
      t3 = zb(i2, i1)
      t4 = zb(i3, i2)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = zb(i3, i1)
      t10 = t2 * t9
      t11 = t5 * zb(i4, i1) + t10
      t12 = mt ** 2
      t13 = za(i1, i3)
      t14 = za(i1, i4)
      t15 = t13 * t4
      t16 = t14 * t6
      t17 = t16 + t15
      t18 = t1 * t3
      t19 = t18 * t17
      t6 = t5 * t6
      t20 = t7 * t8
      t21 = t2 * t4
      t22 = t18 * (t21 + t6 + t20)
      t23 = 16 * t12 * t11 * t19 + 4 * t22 ** 2
      t23 = sqrt(t23)
      t22 = 2 * t22
      t24 = -t23 - t22
      t19 = 0.1q1 / t19
      t25 = t21 + t6
      t15 = t18 * t19 * (t15 * t25 + t16 * t25 + t20 * (t16 + t15))
      t16 = t23 * t19 * t17 / 2
      t6 = 2 * t21 + 2 * t6 + 2 * t20
      t17 = -t16 + t6 - t15
      t20 = t23 - t22
      t6 = t16 + t6 - t15
      t11 = 0.1q1 / t11
      t15 = (0.1q1 / 0.4q1)
      t10 = t10 / 2
      t16 = t4 * (t13 * t24 * t19 * t15 + t2) - t10 * t17 * t11
      t10 = t4 * (t13 * t20 * t19 * t15 + t2) - t10 * t6 * t11
      t13 = -t24 * t15 * t19 * t14 * t4 + t17 * t11 * t5 * t9 / 2
      t21 = -t20 * t15 * t19 * t14 * t4 + t6 * t11 * t5 * t9 / 2
      t10 = 0.1q1 / t10
      t8 = 0.1q1 / t8
      t16 = 0.1q1 / t16
      t22 = 0.1q1 / t7
      t2 = 0.1q1 / t2
      t23 = t24 ** 2
      t25 = t24 * t23
      t26 = t20 ** 2
      t27 = t20 * t26
      t28 = t23 * t16
      t29 = t26 * t10
      t30 = t28 * t13 + t29 * t21
      t28 = t29 + t28
      t29 = t16 * t25
      t31 = t27 * t10
      t32 = t31 + t29
      t33 = t9 ** 2
      t34 = t5 ** 2
      t35 = t4 ** 2
      t36 = t19 ** 2
      t37 = t19 * t36
      t38 = t12 * t9
      t39 = t2 * t19
      t40 = t1 ** 2
      ret = -t15 * t2 * t36 * t1 * (t3 * (t1 * t5 * t8 * t19 * t32 * t35
     # - t30 * t5 + t5 * (t28 * t5 + (-t27 - t25) * t19 * t22 * t8 * t14
     #) * t4 + t8 * (t1 * t30 + (t17 * t23 + t26 * t6) * t22 * t11 * t34
     #) * t9) + t1 * t4 * t33 * t8 * t28 * t7) - t18 * t2 * t37 * t4 * (
     #-t14 * t5 * t32 + (t31 * (-t6 * t11 * t5 + t14) + t29 * (-t17 * t1
     #1 * t5 + t14)) * t8 * t9 * t1) / 16 + t37 * t40 * t3 * t4 * t2 * t
     #8 * (t29 * t13 + t31 * t21) / 8 - t36 ** 2 * t14 * t3 * t40 * t35 
     #* t2 * t8 * (t10 * t26 ** 2 + t16 * t23 ** 2) / 32 - t39 * t5 * (t
     #5 * ((-t24 - t20) * t22 * t5 * t3 + t38 * (t10 * t20 + t16 * t24))
     # + (t33 * (t20 * (-t10 * t12 - 1) + t24 * (-t12 * t16 - 1)) + (-t2
     #6 - t23) * t19 * t22 * t5 * t4 * t3) * t8 * t1) - 2 * t39 * t38 * 
     #t34 * t22 * t8 * (t24 + t20)

      hjetmass_box_pmpm_0_0_s34_mhsq_s12_s234 = ret/32q0/(0,1q0)
      return

      end function
