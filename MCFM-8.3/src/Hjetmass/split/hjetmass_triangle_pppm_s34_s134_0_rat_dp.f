
      double complex function hjetmass_triangle_pppm_s34_s134_0_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i1, i4)
      t2 = za(i2, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = zb(i4, i2)
      t7 = t5 * t3
      t8 = t2 * t6
      t9 = t8 + t7
      if ( dreal(t9) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t9 = cg * cdsqrt(t9 ** 2) + t7 + t8
      t10 = za(i3, i4)
      t11 = (0.1D1 / 0.2D1)
      t12 = t11 * t9
      t13 = t1 * t4
      t14 = t13 * (t12 - t7)
      t15 = za(i2, i3)
      t16 = t1 * t3
      t17 = t4 * (t12 * t10 + t16 * t15)
      t18 = zb(i3, i2)
      t19 = t5 * t4
      t20 = t15 * t18
      t7 = t7 * t13 * (t20 + t7 + t8)
      t21 = t11 * t19 * t9 * (-t10 * t18 + t16) - t7
      t22 = zb(i4, i3)
      t23 = t1 * (t12 * t22 + t19 * t18)
      t24 = zb(i3, i1)
      t12 = t1 * (-t2 * t18 * t4 + t12 * t24)
      t13 = t11 * t9 * (t10 * t22 + t24 * za(i1, i3)) - t20 * t13
      t7 = t11 * t16 * t9 * (-t15 * t22 + t19) - t7
      t7 = 0.1D1 / t7
      t11 = 0.1D1 / t14
      t16 = 0.1D1 / t23
      t19 = 0.1D1 / t4
      t10 = 0.1D1 / t10
      t17 = 0.1D1 / t17
      t20 = 0.1D1 / t5
      t21 = 0.1D1 / t21
      t9 = 0.1D1 / t9
      t15 = 0.1D1 / t15
      t22 = t3 ** 2
      t24 = t7 ** 2
      t7 = t3 * t7
      t25 = t16 * t15
      t9 = t9 * t10 * t2 ** 2
      ret = 32 * t9 * (t19 * (-t8 * t3 * (t7 * t23 + t15) * t20 ** 2 + t
     #14 * (t2 * (t6 * (-t18 * t13 * t16 ** 2 * t15 ** 2 + t22 * t18 * t
     #13 * t24 - t3 * t22 * t23 * t24) + t12 * t22 * t24 * t6 ** 2) + t3
     # * t18 * (t25 + t7)) * t20) + t1 ** 2 * t22 ** 2 * t4 * t11 * (t5 
     #* t18 * t21 + t17)) - 64 * t9 * t3 * t12 * t6 * t20 * t19 * (t25 +
     # t7)

      hjetmass_triangle_pppm_s34_s134_0_rat_dp = ret/32d0/(0,1d0)
      return

      end function
