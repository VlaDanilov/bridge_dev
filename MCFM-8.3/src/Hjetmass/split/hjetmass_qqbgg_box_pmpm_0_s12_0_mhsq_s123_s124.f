
      complex*32 function hjetmass_qqbgg_box_pmpm_0_s12_0_mhsq_s123_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = zb(i4, i2)
      t6 = za(i3, i4)
      t7 = zb(i4, i3)
      t8 = t1 * t5
      t9 = t3 * t4
      t10 = t6 * t7
      t11 = zb(i3, i2)
      t12 = za(i1, i3)
      t13 = za(i2, i3)
      t14 = za(i1, i2)
      t15 = zb(i2, i1)
      t16 = t2 * t12
      t17 = t11 * t13
      t16 = t10 * (t10 * t14 * t15 - t8 * (t17 + t16) + t9 * (-t17 - t16
     #))
      t17 = t2 * t3
      t18 = t1 * t11
      t19 = t17 + t18
      t20 = mt ** 2
      t21 = t6 ** 2 * t7 ** 2 * (t12 * t4 + t13 * t5)
      t22 = 16 * t10 * t20 * t19 * t21 + 4 * t16 ** 2
      t22 = sqrt(t22)
      t16 = -2 * t16
      t23 = -t22 + t16
      t21 = 0.1q1 / t21
      t24 = 0.1q1 / t7
      t8 = (t9 + t10 + t8) * t24
      t9 = t8 * t11
      t10 = t23 * t21 * t6 / 4
      t16 = t22 + t16
      t22 = t16 * t21 * t6 / 4
      t8 = t8 * t2
      t25 = t10 * t4 - t8
      t4 = t22 * t4 - t8
      t8 = 0.1q1 / t14
      t14 = 0.1q1 / t15
      t15 = 0.1q1 / t6
      t26 = t16 * t4 + t23 * t25
      t27 = t23 + t16
      t28 = t15 * t26 + t2 * t27
      t29 = t14 * t21
      t30 = t23 ** 2
      t31 = t16 ** 2
      t32 = t1 * t7
      t33 = t8 * t21
      t34 = t2 ** 2
      ret = 16 * t20 * t1 * t2 * t8 * t15 * t14 * t24 * t19 + 3 * t29 * 
     #t6 * t34 * (t23 + t16) - 8 * t34 * t14 * t24 * t19 + t32 * t21 ** 
     #2 * t13 * t8 * (t31 + t30) / 2 - 2 * t29 * (-t2 * t26 + (t17 * t28
     # + t18 * t28) * t8 * t13) + t33 * (t1 * (t32 * t27 + (t27 * t12 * 
     #t1 - t3 * t13 * t27) * t15 * t2) + t29 * t13 * t7 * (t25 * t30 + t
     #31 * t4 + (t31 + t30) * t6 * t2)) - 4 * t14 * t2 * t1 * (t33 * t20
     # * t27 + (t11 * (t25 + t4) + t2 * (-t10 * t5 - t22 * t5 + 2 * t9))
     # * t15 * t24)

      hjetmass_qqbgg_box_pmpm_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function
