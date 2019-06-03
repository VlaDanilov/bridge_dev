
      complex*32 function hjetmass_triangle_pmpm_s12_s124_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i3, i4)
      t4 = zb(i4, i3)
      t5 = t1 * t2
      t6 = t3 * t4
      t7 = t5 + t6
      if ( real(t7) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t7 = cg * sqrt(t7 ** 2) + t5 + t6
      t8 = za(i2, i4)
      t9 = za(i1, i4)
      t10 = zb(i4, i1)
      t11 = t7 / 2
      t12 = t9 * t10
      t13 = t12 * (t11 - t5)
      t14 = za(i2, i3)
      t15 = zb(i3, i2)
      t16 = zb(i4, i2)
      t17 = t1 * t10
      t11 = t9 * (-t11 * t16 + t17 * t15)
      t18 = -t13
      t5 = -t7 * t9 * t2 * (t14 * t16 + t17) / 2 + t5 * t12 * (t14 * t15
     # + t5 + t6)
      t6 = 0.1q1 / t15
      t11 = 0.1q1 / t11
      t12 = 0.1q1 / t14
      t10 = 0.1q1 / t10
      t9 = 0.1q1 / t9
      t1 = 0.1q1 / t1
      t5 = 0.1q1 / t5
      t15 = t2 * t5
      t16 = t14 * t15
      t17 = t12 * t11
      t19 = t2 * t6
      t20 = t1 * t10 ** 2
      ret = 16 * t20 * (t19 * (t2 * t18 * (-t16 + t11) + (t16 - t11) * t
     #1 * t9 * t10 * t13 ** 2) + t9 * t18 * (-t17 * t2 * t6 + t2 ** 2 * 
     #t5 * (-t14 * t18 * t5 + t6) + t18 * t11 ** 2 * t12) * (t2 * za(i1,
     # i2) + t4 * t8) * t3) + 8 * t20 * t19 * t7 * t8 * t13 * t9 * (t15 
     #- t17)

      hjetmass_triangle_pmpm_s12_s124_0_rat = ret/32q0/(0,1q0)
      return

      end function
