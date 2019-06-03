
      double complex function hjetmass_triangle_ppmm_s23_s123_0_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i3, i4)
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = za(i2, i4)
      t5 = zb(i3, i2)
      t6 = t3 * zb(i3, i1) + t4 * t5
      t7 = zb(i4, i1)
      t8 = zb(i4, i2)
      t9 = za(i1, i2)
      t10 = t3 * t7
      t11 = t4 * t8
      t12 = t11 + t10
      if ( dreal(t12) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t12 = cg * cdsqrt(t12 ** 2) + t10 + t11
      t13 = za(i2, i3)
      t14 = (0.1D1 / 0.2D1)
      t15 = t14 * t12
      t16 = t9 * t7
      t17 = t16 * t1
      t18 = -t2 * (t15 * t13 + t17)
      t19 = zb(i4, i3)
      t15 = t9 * t2 * (t15 - t10)
      t20 = 0.1D1 / t12
      t17 = t2 * (-2 * t17 * t20 - t13)
      t20 = t12 * t3
      t13 = t14 * t20 * t2 * (-t13 * t19 + t16) - t16 * t3 * t2 * (t1 * 
     #t19 + t10 + t11)
      t9 = 0.1D1 / t9
      t5 = 0.1D1 / t5
      t14 = 0.1D1 / t18
      t16 = 0.1D1 / t15
      t19 = 0.1D1 / t19
      t21 = -0.1D1 / t13
      t15 = -0.1D1 / t15
      t17 = 0.1D1 / t17
      t13 = 0.1D1 / t13
      t22 = 0.1D1 / t3
      t23 = 0.1D1 / t7
      t24 = t8 ** 2
      t25 = t14 * t19
      t26 = t3 * t21
      t27 = t12 * t2
      t28 = t3 * t18
      t14 = t14 * t19 ** 2
      t29 = t3 ** 2 * t7
      t30 = t9 * t1
      t31 = t5 * t2 ** 2
      t32 = t21 ** 2
      ret = -8 * t31 * (t4 * (t1 * t6 * t12 * t16 * (-t28 * t13 ** 2 + t
     #14 * t22) * t24 + t1 * (t9 * (t26 + t25) + t27 * (t13 * t3 - t25) 
     #* t16) * t8) - t27 * t15 ** 2 * t23 * (t26 * t18 + t19) * t24 * t4
     # ** 2 + t30 * (t25 * t10 + t17 * t19 + t29 * t21)) + 16 * t30 * t5
     # * t8 * t6 * t2 * (-t29 * t18 * t32 + t14 * (t11 * t22 + t7) - t11
     # * t28 * t32) - 4 * t31 * t30 * t20 * t21

      hjetmass_triangle_ppmm_s23_s123_0_rat_dp = ret/32d0/(0,1d0)
      return

      end function
