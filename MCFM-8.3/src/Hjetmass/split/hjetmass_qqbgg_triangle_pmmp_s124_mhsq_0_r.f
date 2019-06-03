
      double complex function hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat_dp
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i2, i3)
      t2 = zb(i2, i1)
      t3 = za(i3, i4)
      t4 = zb(i4, i1)
      t5 = t1 * t2 - t3 * t4
      t6 = za(i1, i3)
      t7 = zb(i3, i1)
      t8 = zb(i3, i2)
      t9 = zb(i4, i3)
      t10 = t6 * t7
      t11 = t1 * t8
      t12 = t3 * t9
      t13 = t10 + t11 + t12
      if ( dreal(t13) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t13 ** 2) + t10 + t11 + t12
      t12 = za(i1, i2)
      t13 = za(i2, i4)
      t14 = za(i1, i4)
      t15 = zb(i4, i2)
      t16 = t12 * t2
      t17 = t14 * t4
      t18 = t13 * t15 + t16 + t17
      t11 = 0.1D1 / t11
      t19 = -0.2D1 * t1 * t7 * t18 * t11 + t13 * t4
      t20 = 0.2D1 * t3 * t7 * t18 * t11 - t13 * t2
      t10 = -0.2D1 * t10 * t18 * t11 + t16 + t17
      t16 = t6 * (t10 * t7 + t19 * t8)
      t17 = t1 * t15 + t4 * t6
      t21 = 0.2D1 * t6
      t8 = t7 * (t1 * (-t21 * t8 * t18 * t11 + t14 * t15) + t10 * t6)
      t21 = t21 * t18 * t9 * t11 - t12 * t15
      t22 = 0.1D1 / t7
      t18 = 0.1D1 / t18
      t23 = 0.1D1 / t12
      t24 = 0.1D1 / t2
      t14 = 0.1D1 / t14
      t25 = t4 * t23 * t24
      t8 = 0.1D1 / t8
      t21 = 0.1D1 / t21
      t26 = 0.1D1 / t3
      t13 = 0.1D1 / t13
      t16 = 0.1D1 / t16
      t27 = 0.1D1 / t9
      t20 = 0.1D1 / t20
      t28 = t6 * t16
      t21 = t26 * t22 * t21
      t29 = (t20 * t27 + t28) * t19
      t30 = t1 * t10
      t31 = t5 * t1
      t32 = t24 * t23 * t13 * t11
      t33 = t32 * t1
      t34 = t9 * t22
      ret = -0.64D2 * t1 * t4 * t22 * t18 * (t25 + t14) - 0.16D2 * t33 *
     # (t31 * ((-t21 - t8) * t10 * t9 + t26) + (-t30 * t8 - t29) * t17 *
     # t7) - 0.128D3 * t31 * t11 * t18 * (t14 * (t12 * t13 + t34) + (t34
     # * t23 + t13) * t24 * t4) + 0.8D1 * t25 * t13 * t1 * (t30 * (t21 +
     # t8) + t29) + 0.32D2 * t32 * t19 ** 2 * t7 * (t3 * t17 * t20 ** 2 
     #* t27 + t6 ** 2 * t9 * t16 ** 2 * (t1 * (t15 * t3 + t2 * t6) - t5 
     #* t6)) + 0.48D2 * t33 * t5 * t19 * (t28 * t9 + t20)

      hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0_rat_dp = ret/32d0*(0,1d0)
      return

      end function
