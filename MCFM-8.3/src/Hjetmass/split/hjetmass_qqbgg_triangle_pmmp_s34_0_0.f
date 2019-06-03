
      complex*32 function hjetmass_qqbgg_triangle_pmmp_s34_0_0
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
      t3 = zb(i3, i1)
      t4 = za(i1, i4)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = t4 * t3
      t8 = t2 * t6
      t9 = zb(i4, i1)
      t10 = zb(i4, i2)
      t11 = za(i1, i2)
      t12 = zb(i2, i1)
      t13 = t5 * t3
      t14 = t6 * t9
      t15 = t3 * t10
      t16 = t3 ** 2
      t17 = -t14 + t15
      t18 = t2 ** 2
      t19 = t8 + t7
      t19 = 0.1q1 / t19
      t20 = 0.1q1 / t11
      t21 = 0.1q1 / t12
      t22 = t19 ** 2
      ret = -16 * za(i3, i4) * (-t1 ** 2 * t2 * t3 * t16 + t11 * t12 * (
     #t2 * (-2 * t1 * t16 + t2 * (-t14 + 2 * t15) - t13 * t6) + t7 * (t2
     # * t9 + t13)) + t10 * t17 * t2 * t18 + t1 * t5 * t16 * (-t8 + t7) 
     #+ t4 * t9 * t17 * t18 + t4 * t5 ** 2 * t16 * t6) * zb(i4, i3) * t2
     #0 * t21 * t19 * t22

      hjetmass_qqbgg_triangle_pmmp_s34_0_0 = ret/32q0/(0,1q0)
      return

      end function
