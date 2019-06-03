
      complex*32 function hjetmass_triangle_pmpm_s23_0_0
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
      t3 = za(i3, i4)
      t4 = zb(i2, i1)
      t5 = zb(i4, i2)
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = zb(i3, i1)
      t9 = zb(i4, i3)
      t10 = t4 * t9
      t11 = t8 * t5
      t12 = -t10 + t11
      t13 = t3 * t12
      t14 = t6 * t4 * t7
      t15 = t13 - t14
      t16 = za(i2, i4)
      t17 = t16 * t5
      t18 = t17 * t4
      t19 = t4 ** 2
      t20 = t3 ** 2
      t21 = t2 ** 2
      t22 = t5 ** 2
      t23 = t3 * t5
      t24 = t2 * t4
      t9 = t3 * t9
      t25 = t23 * t2
      t26 = t24 + t23
      t27 = 0.1q1 / t7
      t28 = 0.1q1 / t6
      t26 = 0.1q1 / t26
      t29 = 0.1q1 / t5 ** 2
      t30 = 0.1q1 / t2 ** 2
      t31 = t26 ** 2
      ret = -16 * za(i2, i3) * zb(i3, i2) * (-t1 * (-t3 * t20 * t5 * t22
     # * t15 + t2 * t21 * t4 * t19 * (-t13 + t18) + t23 * t21 * t19 * (-
     #t18 - 3 * t13) - 3 * t24 * t20 * t22 * t15) + t2 * (-t21 * t19 * (
     #t4 * (t22 * t16 ** 2 + t9 * t6 * t7) + 3 * t17 * t3 * t12 + t24 * 
     #t16 * t12) - t25 * (t14 * (t3 * (-2 * t11 + 3 * t10) - 2 * t18 + t
     #14) + t13 * (3 * t16 * t4 + t3 * t8) * t5) - t20 * t22 * (t3 * (t9
     # + t17) * t12 + t6 ** 2 * t4 * t7 ** 2) + t25 * t1 ** 2 * t19 ** 2
     #)) * t30 * t28 * t27 * t29 * t26 * t31

      hjetmass_triangle_pmpm_s23_0_0 = ret/32q0/(0,1q0)
      return

      end function
