
      complex*32 function hjetmass_triangle_ppmm_s14_0_0
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
      t3 = zb(i4, i2)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = za(i2, i4)
      t7 = zb(i4, i3)
      t8 = t4 * t5
      t9 = t2 * t7
      t10 = za(i1, i3)
      t11 = zb(i3, i1)
      t12 = zb(i2, i1)
      t13 = t12 * t7
      t14 = t11 * t3 - t13
      t15 = t8 * t3
      t16 = t3 ** 2
      t17 = t3 * t16
      t18 = t10 * t7
      t19 = t1 * t3 + t18
      t19 = 0.1q1 / t19
      t20 = 0.1q1 / t4
      t21 = 0.1q1 / t7
      t22 = 0.1q1 / t1
      t23 = 0.1q1 / t5
      t24 = t19 ** 2
      ret = 16 * za(i1, i4) * zb(i4, i1) * (t10 * (t1 * (t4 ** 2 * t5 **
     # 2 * t16 + t6 ** 2 * t16 ** 2 + t15 * (t7 * (t10 * t12 + t2 * t3) 
     #+ 2 * t6 * t16) + t13 * t10 ** 2 * t14 + t9 * t6 * t17) + t18 * (t
     #10 * t11 + t8) * (t10 * t14 + t15)) - t1 ** 2 * t2 * t17 * (t3 * t
     #6 + t8 + t9)) * t22 * t20 * t23 * t21 * t19 * t24

      hjetmass_triangle_ppmm_s14_0_0 = ret/32q0/(0,1q0)
      return

      end function
