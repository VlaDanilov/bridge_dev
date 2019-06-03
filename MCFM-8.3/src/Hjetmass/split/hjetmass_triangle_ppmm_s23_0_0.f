
      complex*32 function hjetmass_triangle_ppmm_s23_0_0
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
      t2 = za(i3, i4)
      t3 = zb(i3, i1)
      t4 = za(i1, i3)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = zb(i4, i3)
      t8 = t2 * t7
      t9 = t4 * t3
      t10 = t5 * t6
      t11 = za(i2, i4)
      t12 = zb(i4, i2)
      t13 = zb(i2, i1)
      t14 = t13 * t7
      t15 = t12 * t3 - t14
      t16 = t10 * t3
      t17 = t3 ** 2
      t18 = t3 * t17
      t19 = t11 * t7
      t20 = t1 * t3 + t19
      t21 = 0.1q1 / t5
      t22 = 0.1q1 / t1
      t20 = 0.1q1 / t20
      t23 = 0.1q1 / t7
      t24 = 0.1q1 / t6
      t25 = t20 ** 2
      ret = -16 * za(i2, i3) * zb(i3, i2) * (t11 * (-t1 * (t4 * t18 * (t
     #9 + 2 * t10 + t8) + t14 * t11 ** 2 * t15 + t5 ** 2 * t17 * t6 ** 2
     # + t16 * (t11 * t13 + t2 * t3) * t7) - t19 * (t11 * t12 + t10) * (
     #t11 * t15 + t16)) + t1 ** 2 * t2 * t18 * (t8 + t9 + t10)) * t22 * 
     #t21 * t24 * t23 * t20 * t25

      hjetmass_triangle_ppmm_s23_0_0 = ret/32q0/(0,1q0)
      return

      end function
