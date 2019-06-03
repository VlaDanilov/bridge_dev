
      complex*32 function hjetmass_triangle_pmpm_s34_0_0
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
      t2 = za(i1, i3)
      t3 = za(i2, i3)
      t4 = zb(i2, i1)
      t5 = zb(i4, i1)
      t6 = zb(i4, i2)
      t7 = t3 * t6
      t8 = t2 * t5
      t9 = t7 + t8
      t10 = za(i1, i4)
      t11 = zb(i3, i2)
      t12 = zb(i3, i1)
      t13 = t12 * t6
      t14 = t11 * t5
      t15 = -t13 + t14
      t16 = za(i2, i4)
      t17 = t16 * t6
      t18 = 3 * t10 * t5
      t19 = t2 ** 2
      t20 = t2 * t19
      t21 = t3 ** 2
      t22 = t5 ** 2
      t23 = t6 ** 2
      t24 = t17 * t5
      t25 = -3 * t3 * t15 + t24
      t26 = t10 * t22
      t27 = 0.1q1 / t9
      t28 = 0.1q1 / t4
      t29 = 0.1q1 / t1
      t30 = t27 ** 2
      t31 = 0.1q1 / t2 ** 2
      t32 = 0.1q1 / t6 ** 2
      ret = 16 * za(i3, i4) * (t10 * t21 ** 2 * t6 * t23 * t15 - t7 * t1
     #9 * (t26 * (t25 + t26) + t7 * (t12 * t3 + 3 * t16 * t5) * t15) + t
     #20 * t22 * (-t8 * t16 * t15 + t17 * t25 + t10 * t5 * (t3 * t15 + t
     #24)) + t1 * t3 * t4 * t5 * (t6 * (t19 * (t3 * (-2 * t13 + 3 * t14)
     # - 2 * t17 * t5) + t21 * t10 * t23) + t20 * t11 * t22 + t18 * t2 *
     # t3 * t23) - t2 * t3 * t21 * t23 * (t11 * t3 + t17 - t18) * t15 + 
     #t8 * t7 * t1 ** 2 * t4 ** 2 * t9) * zb(i4, i3) * t29 * t31 * t28 *
     # t32 * t27 * t30

      hjetmass_triangle_pmpm_s34_0_0 = ret/32q0/(0,1q0)
      return

      end function
