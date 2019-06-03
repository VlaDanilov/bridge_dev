
      complex*32 function hjetmass_triangle_pmpm_s34_s234_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = za(i1, i2)
      t2 = za(i2, i3)
      t3 = zb(i3, i1)
      t4 = zb(i2, i1)
      t5 = t1 * t4
      t6 = za(i1, i3) * t3
      t7 = t5 + t6
      if ( real(t7) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t7 = cg * sqrt(t7 ** 2) + t5 + t6
      t8 = za(i2, i4)
      t9 = zb(i3, i2)
      t10 = za(i1, i4)
      t11 = (0.1q1 / 0.2q1)
      t12 = t11 * t7
      t13 = t9 * (-t10 * t2 * t3 + t12 * t8)
      t14 = t2 * t9
      t15 = t14 * (t12 - t5)
      t16 = zb(i4, i1)
      t17 = zb(i4, i3)
      t18 = t1 * t9
      t19 = -t2 * (t12 * t17 + t18 * t16)
      t7 = t7 * t2
      t5 = t11 * t7 * t4 * (-t10 * t17 + t18) - t5 * t14 * (t10 * t16 + 
     #t5 + t6)
      t11 = 0.1q1 / t10
      t15 = 0.1q1 / t15
      t16 = 0.1q1 / t16
      t17 = 0.1q1 / t19
      t5 = 0.1q1 / t5
      t9 = 0.1q1 / t9
      t19 = t3 ** 2
      t20 = t3 * t16
      t11 = t11 * t17
      t21 = t4 * t5
      t22 = t21 * t10
      t5 = t1 * t4 ** 2 * t10 * t5 ** 2
      t3 = t15 * t3 * t19
      ret = -16 * t3 * t2 * t1 * (t2 * (t17 * (-t11 * t1 * t13 + t20) - 
     #t22 * t20 + t5 * (t20 * t14 * (t12 - t6) + t13)) + t20 * t18 * t4 
     #* (t15 * t17 - t22 * t15 - t5) * t2 ** 2 + t13 * t16 * (-t11 + t21
     #) * t9) - 8 * t3 * t7 * t1 * t8 * t16 * (t11 - t21)

      hjetmass_triangle_pmpm_s34_s234_0_rat = ret/32q0/(0,1q0)
      return

      end function
