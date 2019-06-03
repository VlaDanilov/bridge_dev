
      double complex function hjetmass_triangle_pmpm_s23_s234_0_rat_dp 
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
      t2 = zb(i3, i1)
      t3 = zb(i4, i1)
      t4 = za(i2, i3)
      t5 = za(i2, i4)
      t6 = t2 * t4 + t3 * t5
      t7 = za(i1, i3)
      t8 = t7 * t2
      t9 = t1 * t3
      t10 = t8 + t9
      if ( dreal(t10) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t9 = cg * cdsqrt(t10 ** 2) + t8 + t9
      t10 = zb(i4, i3)
      t11 = za(i1, i2)
      t12 = za(i3, i4)
      t13 = (0.1D1 / 0.2D1)
      t14 = t13 * t9
      t15 = t11 * t12
      t16 = t15 * t2
      t17 = t10 * (t14 * t5 - t16)
      t18 = zb(i2, i1)
      t19 = t4 * zb(i3, i2) + t5 * zb(i4, i2)
      t20 = t15 * t18 * t10
      t13 = t13 * t9 * t19 - t20
      t15 = t15 * t3
      t14 = -t10 * (t14 * t4 + t15)
      t9 = 0.1D1 / t9
      t16 = t10 * (-2 * t16 * t9 + t5)
      t19 = -2 * t20 * t9 + t19
      t4 = t10 * (-2 * t15 * t9 - t4)
      t7 = (t1 * t4 + t11 * t19 + t16 * t7) * t18
      t15 = 0.1D1 / t19
      t4 = 0.1D1 / t4
      t7 = 0.1D1 / t7
      t14 = 0.1D1 / t14
      t19 = 0.1D1 / t11
      t20 = 0.1D1 / t10
      t21 = 0.1D1 / t18
      t13 = 0.1D1 / t13
      t22 = 0.1D1 / t12
      t23 = t2 * t21
      t24 = t6 ** 2
      t25 = t7 * t4
      ret = 16 * t6 * ((t14 * (-t17 ** 2 * t19 * t20 * t22 * t13 ** 2 + 
     #t23 * (t17 * t19 * t20 * t22 - t2) * t13) + t23 * t3 * t17 * t13 *
     # t14 ** 2) * t6 * t1 + t5 * t2 * t16 * t22 * t20 * t4 * (-t19 * t2
     #1 * t15 + t7)) + 64 * t25 * t9 ** 2 * t11 * t24 * t2 * t1 * (-t3 *
     # t16 * (t1 * t18 * t7 + t4) + t2) + 32 * t25 * t22 * t20 * t9 * t1
     #6 * t24 * t1 * (t11 * (t16 * t7 * t18 ** 2 + t2 * t12 * t10 * (-2 
     #* t8 * t9 + 1) * t7 * t18) - t2)

      hjetmass_triangle_pmpm_s23_s234_0_rat_dp = ret/32d0/(0,1d0)
      return

      end function
