
      double complex function hjetmass_triangle_pmpm_s12_s123_0_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i2, i4)
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i3, i1)
      t5 = t1 * t3 + t2 * t4
      t6 = zb(i4, i2)
      t7 = zb(i4, i3)
      t8 = t1 * t6
      t9 = t2 * t7
      t10 = t9 + t8
      if ( dreal(t10) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t8 = cg * cdsqrt(t10 ** 2) + t8 + t9
      t9 = za(i2, i3)
      t10 = zb(i3, i2)
      t11 = zb(i4, i1)
      t12 = t8 / 2
      t13 = t1 * t10 * t11
      t14 = t9 * (t12 * t4 - t13)
      t15 = za(i1, i4)
      t16 = t3 * za(i1, i2) + t4 * za(i1, i3)
      t17 = t15 * t9 * t10 * t11
      t18 = -t8 * t16 / 2 + t17
      t19 = t2 * t10 * t11
      t12 = -t9 * (t12 * t3 + t19)
      t8 = 0.1D1 / t8
      t13 = t9 * (-2 * t13 * t8 + t4)
      t16 = -2 * t17 * t8 + t16
      t3 = t9 * (-2 * t19 * t8 - t3)
      t6 = (t11 * t16 + t13 * t6 + t3 * t7) * t15
      t3 = 0.1D1 / t3
      t9 = 0.1D1 / t9
      t17 = 0.1D1 / t18
      t10 = 0.1D1 / t10
      t18 = 0.1D1 / t11
      t12 = 0.1D1 / t12
      t16 = 0.1D1 / t16
      t19 = 0.1D1 / t15
      t6 = 0.1D1 / t6
      t7 = t5 * t7 * t10
      t20 = t7 * t18
      t5 = t9 * t5
      ret = 16 * t5 * ((-t1 * t2 * t12 ** 2 * t19 * t17 - t20 * t12 * t1
     #7 ** 2) * t14 ** 2 + t1 * t4 * t13 * t10 * t3 * (-t19 * t18 * t16 
     #+ t6) + t1 * t19 * (-t20 + t1) * t12 * t17 * t14) + 32 * t5 * t6 *
     # t8 * t3 * t13 * (t7 * (t11 * t15 ** 2 * t13 * t6 - t1) + t1 * t11
     # * (-t2 * t13 * t3 + t1))

      hjetmass_triangle_pmpm_s12_s123_0_rat_dp = ret/32d0/(0,1d0)
      return

      end function
