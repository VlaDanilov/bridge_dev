
      complex*32 function hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i3)
      t5 = zb(i2, i1)
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = zb(i4, i2)
      t9 = t3 * t5
      t10 = t6 * t7
      t11 = t1 * t8 + t10 + t9
      t12 = za(i1, i3)
      t13 = zb(i3, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t12 * t2
      t17 = t4 * t13
      t18 = t14 * t15
      t19 = t17 + t18 + t16
      if ( real(t19) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t17 = cg * sqrt(t19 ** 2) + t16 + t17 + t18
      t17 = 0.1q1 / t17
      t18 = -0.2q1 * t4 * t2 * t11 * t17 + t1 * t7
      t9 = -0.2q1 * t16 * t11 * t17 + t10 + t9
      t10 = t12 * (t13 * t18 + t2 * t9)
      t16 = t1 * t13 + t2 * t6
      t19 = 0.2q1 * t12
      t20 = t2 * (t12 * t9 + t4 * (-t19 * t13 * t11 * t17 + t6 * t8))
      t21 = 0.2q1 * t14 * t2 * t11 * t17 - t1 * t5
      t22 = t1 * t15 + t2 * t3
      t8 = t19 * t11 * t15 * t17 - t3 * t8
      t5 = 0.1q1 / t5
      t19 = 0.1q1 / t4
      t10 = 0.1q1 / t10
      t21 = 0.1q1 / t21
      t23 = 0.1q1 / t3
      t7 = 0.1q1 / t7
      t24 = 0.1q1 / t15
      t8 = 0.1q1 / t8
      t20 = 0.1q1 / t20
      t25 = 0.1q1 / t14
      t26 = t2 * t20
      t10 = t12 * t10
      t12 = t9 * (t25 * t8 + t26)
      t27 = t7 * t23 * t5
      t28 = t27 * t2
      t11 = 0.1q1 / t11
      ret = 0.8q1 * t28 * t1 * (t12 + t2 * (t21 * t24 + t10) * t19 * t18
     #) - 0.48q2 * t28 * t9 * t22 * t17 * (t26 * t14 + t8) + 0.16q2 * t2
     #8 * t17 * (t2 * (t22 * t24 * (-t14 * t18 * t19 * t21 + 0.1q1) + t1
     #0 * t18 * (-t14 * t22 * t19 + t16)) + t12 * t16 * t4) + 0.128q3 * 
     #t2 * t1 * t22 * t17 * t23 * t11 * (t14 * t19 * t5 - t7) - 0.32q2 *
     # t27 * t17 * t9 ** 2 * t4 * (-t16 * t15 * t25 * t8 ** 2 + t14 * t2
     # ** 2 * t20 ** 2 * (-t13 * t22 + t2 * (t13 * t3 - t15 * t6))) - 0.
     #64q2 * t1 ** 2 * t2 * t23 * t19 * t5 * t11

      hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0_rat = ret/32q0*(0,1q0)
      return

      end function
