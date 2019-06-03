
      complex*32 function hjetmass_qqbgg_box_pmmp_0_s12_0_mhsq_s123_s124
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
      t2 = za(i1, i4)
      t3 = zb(i3, i1)
      t4 = zb(i4, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = za(i1, i2)
      t10 = za(i3, i4)
      t11 = zb(i2, i1)
      t12 = zb(i4, i3)
      t13 = t1 * t4
      t14 = t5 * t8
      t15 = t14 + t13
      t16 = t10 ** 2 * t12 ** 2
      t17 = t16 * t15
      t18 = t2 * t4
      t19 = t7 * t8
      t20 = t10 * t12
      t21 = t20 * t9 * t11
      t22 = t5 * t6
      t23 = t20 * ((-t19 - t18) * t3 * t1 + t21 - t22 * (t19 + t18))
      t24 = t2 * t3 + t6 * t7
      t25 = mt ** 2
      t23 = 16 * t20 * t25 * t24 * t17 + 4 * t23 ** 2
      t17 = 0.1q1 / t17
      t26 = t14 + t13
      t27 = t1 * t26 * t3 + t22 * t26
      t28 = t19 + t18
      t22 = t28 * t3 * t1 + t22 * t28 - t21
      t23 = t20 * sqrt(t23) * t17 * t15 / 4
      t16 = t16 * t17 * (t18 * t27 + t19 * t27 + t21 * (-t14 - t13)) / 2
      t17 = -t22 - t23 + t16
      t21 = t20 * t24
      t16 = -t22 + t23 + t16
      t21 = 0.1q1 / t21
      t22 = 0.1q1 / t10
      t18 = (t20 + t18 + t19) * t22
      t19 = t18 * t5
      t20 = t17 * t21
      t23 = t20 * t7 * t12 + t19
      t24 = t16 * t21
      t7 = t24 * t7 * t12 + t19
      t18 = t18 * t1
      t19 = 0.1q1 / t12
      t11 = 0.1q1 / t11
      t9 = 0.1q1 / t9
      t27 = t17 + t16
      t28 = t16 ** 2
      t29 = t17 ** 2
      t30 = t22 * t19
      t31 = t27 * t21
      t32 = t5 ** 2
      t33 = -t16 * t7 - t17 * t23
      t16 = t17 + t16
      t17 = t16 * t5 + t19 * t33
      t34 = 8
      ret = t34 * (t9 * t5 * (t5 * t22 * t15 + t21 * t33) + t21 * (t9 * 
     #(t13 * t17 + t14 * t17) - t4 * (t29 + t28) * t21 * t10) * t11 * t3
     #) + 12 * t21 * t32 * t12 * t9 * t16 + 16 * t9 * t11 * (t21 ** 2 * 
     #(-t23 * t29 - t28 * t7 + (t29 + t28) * t12 * t5) * t10 * t3 + t5 *
     # t4 * t25 * (t30 * t26 + t31)) - 4 * t4 * (-t31 * t10 * t4 * t11 -
     # t30 * t9 * (t20 * t2 * t12 + t24 * t2 * t12 + 2 * t18) * t32 + ((
     #t27 * t8 * t3 - t4 * t6 * t27) * t11 * t21 + t1 * t22 * t9 * (t23 
     #+ t7)) * t19 * t5)


      hjetmass_qqbgg_box_pmmp_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function
