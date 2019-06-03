      complex*32 function hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s123
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
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i1)
      t9 = zb(i4, i3)
      t10 = t3 * t8
      t11 = t5 * zb(i4, i2)
      t12 = t11 + t10
      t13 = t7 * t9
      t14 = t13 * t12
      t15 = t1 * t2
      t3 = t3 * t4
      t16 = t5 * t6
      t17 = t13 * (t16 + t15 + t3)
      t18 = za(i2, i4)
      t6 = t4 * za(i1, i4) + t6 * t18
      t19 = mt ** 2
      t20 = 16 * t19 * t6 * t14 + 4 * t17 ** 2
      t20 = sqrt(t20)
      t14 = 0.1q1 / t14
      t21 = t16 + t3
      t22 = t13 * t14
      t10 = t22 * (t10 * t21 + t11 * t21 + t15 * (t11 + t10))
      t11 = 2
      t3 = t11 * (t16 + t15 + t3)
      t12 = t20 * t14 * t12 / 2
      t14 = t3 - t10 - t12
      t3 = t3 - t10 + t12
      t10 = t11 * t17
      t2 = 0.1q1 / t2
      t1 = 0.1q1 / t1
      t6 = 0.1q1 / t6
      t12 = t14 + t3
      t15 = t14 ** 2
      t16 = t3 ** 2
      t17 = t6 ** 2
      ret = t11 * t6 * (t5 ** 2 * t9 * t1 * t12 + (t12 * t8 ** 2 + (-t14
     # * t15 - t3 * t16) * t17 * t1 * t18 * t9 * t4) * t2 * t7) + 4 * t6
     # * t1 * t2 * t5 * (t19 * t8 * t12 + t13 * (t16 + t15) * t6 * t4) +
     # t22 * t17 * t8 * t5 * t1 * t2 * (t15 * (-t10 - t20) + t16 * (-t10
     # + t20)) / 2

      hjetmass_qqbgg_box_pmmp_0_0_s12_mhsq_s34_s123 = ret/32q0/(0,1q0)
      return

      end function
