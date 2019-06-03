
      complex*32 function hjetmass_qqbgg_triangle_pmpm_s34_0_0
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
      t3 = zb(i4, i1)
      t4 = za(i1, i4)
      t5 = zb(i4, i2)
      t6 = t2 * t5
      t7 = t4 * t3
      t8 = t7 + t6
      t9 = za(i2, i3)
      t10 = zb(i3, i1)
      t11 = zb(i3, i2)
      t12 = t10 * t5
      t13 = t11 * t3
      t14 = t13 - t12
      t15 = t3 ** 2
      t16 = t9 ** 2
      t17 = za(i1, i2)
      t18 = zb(i2, i1)
      t19 = t1 * t3
      t5 = t5 * t9 + t19
      t5 = 0.1q1 / t5
      t20 = 0.1q1 / t18
      t21 = 0.1q1 / t17
      t22 = t5 ** 2
      ret = 16 * za(i3, i4) * (t1 * (t16 * t10 * t14 + t2 * t15 * t8) - 
     #t9 * (-t16 * t11 * t14 + t4 * t15 * t8) - t17 * t18 * (t9 * (t3 * 
     #(2 * t7 + t6) - t9 * (2 * t13 - t12)) - t19 * (t10 * t9 + t2 * t3)
     #)) * zb(i4, i3) * t21 * t20 * t5 * t22

      hjetmass_qqbgg_triangle_pmpm_s34_0_0 = ret/32q0/(0,1q0)
      return

      end function
