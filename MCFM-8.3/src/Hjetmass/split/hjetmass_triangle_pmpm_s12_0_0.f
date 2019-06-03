
      complex*32 function hjetmass_triangle_pmpm_s12_0_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i3)
      t2 = za(i2, i4)
      t3 = zb(i3, i2)
      t4 = zb(i4, i1)
      t5 = zb(i4, i2)
      t6 = zb(i3, i1) * t5
      t7 = t3 * t4
      t8 = -t6 + t7
      t9 = za(i1, i4)
      t10 = za(i2, i3)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t11 * t12
      t14 = t13 * t3
      t15 = t2 * t5
      t16 = t10 * t3
      t17 = 3 * t16
      t18 = 3 * t15
      t19 = t15 * t3
      t20 = t3 ** 2
      t21 = t9 ** 2
      t22 = t5 ** 2
      t23 = t1 ** 2
      t24 = t1 * t3
      t25 = t9 * t5
      t26 = t25 + t24
      t11 = 0.1q1 / t11
      t26 = 0.1q1 / t26
      t12 = 0.1q1 / t12
      t27 = 0.1q1 / t1 ** 2
      t28 = t26 ** 2
      t29 = 0.1q1 / t5 ** 2
      ret = 16 * za(i1, i2) * zb(i2, i1) * (t1 * t23 * t20 * (-t24 * t2 
     #* t8 + t9 * (-t18 * t8 + t3 * (t10 * t8 + t13 * t4)) + t19 * (t16 
     #+ t15)) - t25 * t23 * (t10 ** 2 * t20 ** 2 + t10 * t20 * (-3 * t8 
     #* t9 + t19) + t6 * t21 * t8 + t13 * t20 * (2 * t15 - t13) + t9 * t
     #3 * (t18 * t8 + t13 * (2 * t6 - 3 * t7))) + t1 * t21 * t22 * (t9 *
     # (t8 * (t17 - t15) - t4 * t8 * t9) + t14 * (t17 + t13)) + t9 * t21
     # * t10 * t5 * t22 * (t8 * t9 + t14)) * t27 * t11 * t29 * t26 * t28
     # * t12

      hjetmass_triangle_pmpm_s12_0_0 = ret/32q0/(0,1q0)
      return

      end function
