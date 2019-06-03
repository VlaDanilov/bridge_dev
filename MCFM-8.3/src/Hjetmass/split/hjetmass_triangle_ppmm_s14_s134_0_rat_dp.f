
      double complex function hjetmass_triangle_ppmm_s14_s134_0_rat_dp 
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
      t4 = za(i1, i3)
      t5 = zb(i3, i2)
      t6 = za(i1, i4)
      t7 = zb(i4, i2)
      t8 = t4 * t5 + t6 * t7
      t9 = za(i1, i2)
      t10 = zb(i4, i3)
      t11 = za(i2, i3)
      t12 = t11 * t5
      t13 = t1 * t7
      t14 = t12 + t13
      if ( dreal(t14) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t13 = cg * cdsqrt(t14 ** 2) + t12 + t13
      t14 = t4 * zb(i3, i1) + t6 * zb(i4, i1)
      t15 = t9 * t2
      t16 = t15 * t3 * t10
      t17 = -t13 * t14 / 2 + t16
      t18 = 0.1D1 / t13
      t14 = -2 * t16 * t18 + t14
      t16 = t15 * t5
      t19 = t10 * (t16 + t13 * t6 / 2)
      t16 = t10 * (2 * t16 * t18 + t6)
      t4 = (-t1 * t10 * (2 * t15 * t7 * t18 - t4) - t11 * t16 + t14 * t9
     #) * t3
      t7 = 0.1D1 / t7
      t4 = 0.1D1 / t4
      t15 = 0.1D1 / t19
      t10 = 0.1D1 / t10
      t6 = 0.1D1 / t6
      t17 = 0.1D1 / t17
      t5 = 0.1D1 / t5
      t19 = 0.1D1 / t11
      t20 = 0.1D1 / t14
      t9 = 0.1D1 / t9
      t21 = t1 ** 2
      t22 = t8 ** 2
      t23 = t9 ** 2
      t12 = t12 * t7
      t24 = t3 * t7
      t25 = t10 * t3
      t26 = t6 * t8
      t11 = t11 * t7
      t27 = t1 * t5
      t16 = 0.1D1 / t16
      ret = -8 * t26 * (t24 * t13 * t21 * t22 * t23 * t5 * t17 * t15 + (
     #-t24 * t1 * t2 * t17 * t9 * t15 - t21 * t2 * t19 * t9 * t15 ** 2) 
     #* t13 * t8 + t25 * t2 * (t9 * (t17 * (t12 + t1) - t20 * t7) + t24 
     #* t4)) + 16 * t25 * t26 * (t1 * t8 * t23 * t17 * (t27 + t11) + (-t
     #12 - t1) * t4 * t18 * t3 * t2) + 32 * t6 * t4 * t18 * t22 * t3 * t
     #1 * (t14 * t3 ** 2 * (t11 * t10 + t27 * (-t8 * t7 * t16 + t10)) * 
     #t4 + t16 * (-t1 * t14 * t19 * t16 + t24) * t2)

      hjetmass_triangle_ppmm_s14_s134_0_rat_dp = ret/32d0/(0,1d0)
      return

      end function
