
      complex*32 function hjetmass_triangle_pmpm_s23_s123_0_rat
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
      t6 = za(i3, i4)
      t7 = t2 * t5
      t8 = t3 * zb(i4, i2)
      t9 = t7 + t8
      if ( real(t9) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t9 = cg * sqrt(t9 ** 2) + t7 + t8
      t10 = (0.1q1 / 0.2q1)
      t11 = t10 * t9
      t12 = t1 * t4
      t13 = t12 * (t11 - t7)
      t14 = za(i2, i3)
      t1 = t1 * t5
      t11 = -t4 * (t1 * t6 + t11 * t14)
      t15 = zb(i3, i1)
      t16 = t15 * t2 + t3 * zb(i3, i2)
      t17 = zb(i4, i3)
      t7 = t10 * t9 * t2 * t4 * (-t14 * t17 + t1) - t7 * t12 * (t17 * t6
     # + t7 + t8)
      t8 = 0.1q1 / t17
      t10 = 0.1q1 / t11
      t6 = 0.1q1 / t6
      t7 = 0.1q1 / t7
      t11 = 0.1q1 / t13
      t13 = t3 ** 2
      t14 = t5 * t7
      t18 = t14 * t16
      t19 = t3 * t4 * t6
      t8 = t10 * t8
      t13 = t11 * t5 * t4 * t3 * t13
      ret = -16 * t13 * (-t8 * t5 * t16 * (t12 * t5 * t10 + t6) + t14 * 
     #t12 * t17 * (-t19 * t11 + t18) * t2 ** 2 + t6 * (t3 * (t1 * t10 * 
     #t11 * t4 ** 2 - t17 * t7 * t4) + t18) * t2 + t19 * t10) - 8 * t13 
     #* t9 * t15 * t6 * (-t2 * t7 + t8)

      hjetmass_triangle_pmpm_s23_s123_0_rat = ret/32q0/(0,1q0)
      return

      end function
