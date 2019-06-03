
      complex*32 function hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i4, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i4)
      t5 = zb(i2, i1)
      t6 = za(i1, i3)
      t7 = zb(i3, i1)
      t8 = zb(i3, i2)
      t9 = t3 * t5
      t10 = t6 * t7
      t11 = t1 * t8 + t10 + t9
      t12 = za(i1, i4)
      t13 = zb(i4, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t12 * t2
      t17 = t4 * t13
      t18 = t14 * t15
      t19 = t16 + t17 + t18
      if ( real(t19) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t17 = cg * sqrt(t19 ** 2) + t16 + t17 + t18
      t17 = 0.1q1 / t17
      t18 = -0.2q1 * t4 * t11 * t2 * t17 + t1 * t7
      t9 = -0.2q1 * t16 * t11 * t17 + t10 + t9
      t10 = t2 * t9
      t16 = t12 * (t13 * t18 + t10)
      t19 = t1 * t13 + t2 * t6
      t20 = -0.2q1 * t12 * t11
      t21 = t2 * (t12 * t9 + t4 * (t20 * t13 * t17 + t6 * t8))
      t22 = -0.2q1 * t14 * t11 * t2 * t17 - t1 * t5
      t23 = -t1 * t15 + t2 * t3
      t20 = t20 * t15 * t17 - t3 * t8
      t11 = 0.1q1 / t11
      t24 = 0.1q1 / t4
      t8 = 0.1q1 / t8
      t25 = 0.1q1 / t5
      t26 = 0.1q1 / t3
      t27 = t1 * t26
      t28 = t27 * t25
      t29 = 0.1q1 / t15
      t16 = 0.1q1 / t16
      t30 = 0.1q1 / t14
      t20 = 0.1q1 / t20
      t21 = 0.1q1 / t21
      t7 = 0.1q1 / t7
      t22 = 0.1q1 / t22
      t12 = t12 * t16
      t16 = t2 * t21
      t31 = t30 * t20
      t32 = t14 * t24
      t33 = t17 * t26 * t7 * t25
      ret = 0.64q2 * t1 * t2 * t24 * t11 * (t28 + t8) + 0.8q1 * t28 * t7
     # * t2 * (t9 * (t31 - t16) + t2 * (t22 * t29 - t12) * t24 * t18) + 
     #0.48q2 * t10 * t23 * t17 * t26 * t25 * t7 * (-t16 * t14 + t20) - 0
     #.16q2 * t33 * t2 * (t2 * (-t23 * t29 * (t32 * t18 * t22 + 0.1q1) +
     # t12 * t18 * (t32 * t23 + t19)) + t19 * (-t31 + t16) * t9 * t4) + 
     #0.128q3 * t17 * t11 * t23 * t2 * (t8 * (t5 * t7 + t32) + t27 * (t3
     #2 * t25 + t7)) - 0.32q2 * t33 * t9 ** 2 * t4 * (t19 * t15 * t30 * 
     #t20 ** 2 + t14 * t2 ** 2 * t21 ** 2 * (-t13 * t23 + t2 * (t13 * t3
     # + t15 * t6)))


      hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function
