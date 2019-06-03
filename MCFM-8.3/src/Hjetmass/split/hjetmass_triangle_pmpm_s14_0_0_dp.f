
      double complex function hjetmass_triangle_pmpm_s14_0_0_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i2)
      t2 = zb(i4, i2)
      t3 = za(i1, i3)
      t4 = zb(i2, i1)
      t5 = za(i3, i4)
      t6 = t3 * t4
      t7 = zb(i3, i1)
      t8 = zb(i4, i3)
      t9 = t7 * t2
      t4 = t4 * t8
      t10 = -t4 + t9
      t11 = za(i2, i4)
      t12 = t11 * t2
      t13 = t5 * t8
      t14 = t3 * t10
      t15 = t2 * (t12 + t13)
      t16 = za(i2, i3)
      t17 = zb(i3, i2)
      t18 = 3 * t13
      t19 = t16 * t17
      t20 = t2 ** 2
      t21 = t3 ** 2
      t22 = t8 ** 2
      t23 = t1 ** 2
      t24 = t3 * t8
      t25 = t1 * t2 + t24
      t25 = 0.1D1 / t25
      t26 = 0.1D1 / t16
      t27 = 0.1D1 / t17
      t28 = 0.1D1 / t3 ** 2
      t29 = t25 ** 2
      t30 = 0.1D1 / t2 ** 2
      ret = 16 * za(i1, i4) * zb(i4, i1) * (t3 * t21 * t11 * t8 * t22 * 
     #(t14 + t15) - t1 * t21 * t22 * (t15 * t13 - t16 ** 2 * t17 ** 2 * 
     #t2 - t19 * (-2 * t11 * t20 + t6 * t8) - t14 * (-t13 + 3 * t12)) + 
     #t1 * t23 * t20 * (t19 * t13 * t2 + t14 * (t3 * t7 + t12 - t18) + t
     #1 * (-t2 * t5 + t6) * t10) + t24 * t23 * t2 * (t3 * (t19 * (3 * t4
     # - 2 * t9) + (3 * t12 - 3 * t13) * t10) + t19 * t2 * (t19 + t18)))
     # * t28 * t26 * t27 * t30 * t25 * t29

      hjetmass_triangle_pmpm_s14_0_0_dp = ret/32d0/(0,1d0)
      return

      end function
