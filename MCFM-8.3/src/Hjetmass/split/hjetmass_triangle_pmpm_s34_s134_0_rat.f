
      complex*32 function hjetmass_triangle_pmpm_s34_s134_0_rat
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
      t2 = za(i1, i4)
      t3 = za(i2, i4)
      t4 = zb(i2, i1)
      t5 = zb(i4, i1)
      t6 = za(i2, i3)
      t7 = t1 * t4
      t8 = t3 * zb(i4, i2)
      t9 = t7 + t8
      if ( real(t9) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t9 = cg * sqrt(t9 ** 2) + t7 + t8
      t10 = (0.1q1 / 0.2q1)
      t11 = t10 * t9
      t12 = t2 * t5
      t13 = t12 * (t11 - t7)
      t14 = za(i3, i4)
      t15 = t2 * t4
      t16 = t5 * (t11 * t14 + t15 * t6)
      t17 = zb(i3, i1)
      t18 = zb(i3, i2)
      t19 = t3 * t18
      t20 = t2 * (t11 * t17 - t19 * t5)
      t7 = t10 * t9 * t1 * t5 * (-t14 * t18 + t15) - t7 * t12 * (t18 * t
     #6 + t7 + t8)
      t10 = 0.1q1 / t18
      t13 = 0.1q1 / t13
      t7 = 0.1q1 / t7
      t6 = 0.1q1 / t6
      t2 = 0.1q1 / t2
      t14 = 0.1q1 / t16
      t16 = t3 ** 2
      t10 = t10 * t14
      t21 = t3 * t6
      t22 = t1 * t7
      t7 = t1 ** 2 * t4 * t18 * t7 ** 2
      t3 = t13 * t5 * t4 * t3 * t16
      ret = 16 * t3 * (t5 * (t14 * (t10 * t4 * t20 + t21) + t22 * t19 * 
     #t6 - t7 * (t21 * t12 * (t11 - t8) + t20)) + t21 * t15 * t1 * (t22 
     #* t18 * t13 + t13 * t14 + t7) * t5 ** 2 - t20 * t6 * (t10 + t22) *
     # t2) + 8 * t3 * t9 * t17 * t6 * (t10 + t22)

      hjetmass_triangle_pmpm_s34_s134_0_rat = ret/32q0/(0,1q0)
      return

      end function
