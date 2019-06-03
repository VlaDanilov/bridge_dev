
      double complex function hjetmass_triangle_pppm_s14_s124_0_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = t1 * t2
      t6 = t3 * t4
      t7 = t5 + t6
      if ( dreal(t7) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t7 = cg * cdsqrt(t7 ** 2) + t5 + t6
      t8 = za(i1, i2)
      t9 = zb(i2, i1)
      t10 = t7 / 2
      t11 = t8 * t9
      t12 = t11 * (t10 - t5)
      t13 = za(i1, i4)
      t14 = za(i2, i4)
      t15 = t13 * t2 + t14 * t4
      t16 = zb(i4, i2)
      t17 = zb(i4, i3)
      t18 = t1 * t9
      t10 = t8 * (t10 * t16 - t18 * t17)
      t19 = za(i3, i4)
      t5 = -t7 * t8 * t2 * (t16 * t19 + t18) / 2 + t5 * t11 * (t17 * t19
     # + t5 + t6)
      t6 = 0.1D1 / t1
      t3 = 0.1D1 / t3
      t5 = 0.1D1 / t5
      t7 = 0.1D1 / t10
      t14 = 0.1D1 / t14
      t13 = 0.1D1 / t13
      t10 = -0.1D1 / t10
      t16 = 0.1D1 / t17
      t9 = 0.1D1 / t9
      t17 = t2 ** 2
      t18 = t2 * t16
      t20 = t17 * t19 ** 2
      t21 = t20 * t5 ** 2
      t22 = t2 * t14
      t23 = t4 * t13
      t24 = t12 * t14
      t8 = t8 * t3
      ret = 16 * t4 * (t15 * (t7 * (t16 * (-t23 - t22) + t24 * t7) + t23
     # * t11 * (-t10 ** 2 + t21) * t1 + t8 * t13 * t14 * t19 * (t10 * (t
     #10 * t12 + t18) - t21 * t12 - t17 * t19 * t5 * t16) - t18 * (t23 +
     # t22) * t5 * t19 - t24 * t21) - t16 * t13 * t14 * t3 * t19 * (-t15
     # * t7 + (-t15 * t5 - t7 * t8) * t19 * t2 - t20 * t8 * t5) * t9 * t
     #6 * t12)

      hjetmass_triangle_pppm_s14_s124_0_rat_dp = ret/32d0/(0,1d0)
      return

      end function
