
      double complex function hjetmass_triangle_pppm_s14_s134_0_rat_dp 
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
      t2 = zb(i3, i2)
      t3 = za(i2, i3)
      t4 = zb(i4, i2)
      t5 = t3 * t2
      t6 = t1 * t4
      t7 = t6 + t5
      if ( dreal(t7) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t6 = cg * cdsqrt(t7 ** 2) + t5 + t6
      t7 = za(i1, i2)
      t8 = za(i1, i4)
      t9 = zb(i4, i3)
      t10 = zb(i2, i1)
      t11 = za(i1, i3)
      t12 = t11 * t2 + t4 * t8
      t13 = zb(i3, i1)
      t14 = zb(i4, i1)
      t15 = za(i3, i4)
      t16 = t11 * t13 + t14 * t8
      t17 = 0.1D1 / t6
      t18 = 2 * t7 * t15
      t19 = -t18 * t10 * t9 * t17 + t16
      t20 = (0.1D1 / 0.2D1)
      t21 = t20 * t6
      t22 = t1 * t10
      t23 = t22 * t9
      t24 = t15 * (-t21 * t13 + t23)
      t6 = -t7 * t15 * t10 * t9 + t20 * t6 * t16
      t16 = t3 * t10 * t9
      t20 = t15 * (t21 * t14 + t16)
      t21 = t9 * (t18 * t2 * t17 + t8)
      t13 = t15 * (2 * t23 * t17 - t13)
      t14 = t15 * (2 * t16 * t17 + t14)
      t16 = t10 * t19
      t23 = (-t13 * t4 - t14 * t2 + t16) * t7
      t11 = (t1 * t9 * (t18 * t4 * t17 - t11) - t19 * t7 + t21 * t3) * t
     #10
      t18 = 0.1D1 / t15
      t11 = 0.1D1 / t11
      t6 = 0.1D1 / t6
      t20 = 0.1D1 / t20
      t25 = 0.1D1 / t7
      t14 = 0.1D1 / t14
      t23 = 0.1D1 / t23
      t8 = 0.1D1 / t8
      t26 = 0.1D1 / t9
      t27 = 0.1D1 / t19
      t3 = 0.1D1 / t3
      t28 = t16 * t11 + t25
      t29 = t1 ** 2
      t30 = t11 ** 2
      t19 = t10 ** 2 * t19
      ret = -32 * t26 * t17 * t8 * ((t2 * t4 * t28 * t3 ** 2 + t18 * t12
     # * t10 * (t21 * t25 ** 2 * t27 - t19 * t21 * t30 + t16 * t2 * t15 
     #* t9 * (-2 * t5 * t17 + 1) * t30) * t3) * t1 * t29 + t2 ** 2 * (t7
     # * t10 * t13 ** 2 * t14 * t23 - t24 ** 2 * t6 * t20)) - 64 * t3 * 
     #t17 * t8 * t2 * t29 * (-t19 * t29 * t4 * t12 * t17 * t30 + t2 * t2
     #6 * t28 - t22 * (t10 * t11 + t25 * t27) * t17 * t12)

      hjetmass_triangle_pppm_s14_s134_0_rat_dp = ret/32d0/(0,1d0)
      return

      end function
