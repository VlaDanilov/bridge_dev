
      double complex function hjetmass_triangle_pppm_s234_0_mhsq_rat_dp 
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
      t2 = zb(i2, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i1)
      t7 = zb(i4, i3)
      t8 = za(i1, i2)
      t9 = zb(i3, i2)
      t10 = zb(i4, i2)
      t11 = za(i3, i4)
      t12 = t3 * t9
      t13 = t5 * t10
      t14 = t11 * t7
      t15 = t12 + t13 + t14
      t16 = za(i1, i3)
      t17 = t8 * t2
      t18 = t16 * t4
      t19 = t1 * t6
      t20 = t17 + t18 + t19
      if ( dreal(t20) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t17 = cg * cdsqrt(t20 ** 2) + t17 + t18 + t19
      t17 = 0.1D1 / t17
      t18 = -2 * t8
      t19 = t18 * t4 * t15 * t17 + t5 * t7
      t13 = t18 * t2 * t15 * t17 + t12 + t13
      t20 = t16 * t19
      t21 = t2 * (t13 * t8 + t20)
      t22 = -2 * t1
      t23 = t22 * t4 * t15 * t17 + t5 * t9
      t18 = t18 * t6 * t15 * t17 - t3 * t7
      t20 = t2 * (t1 * t18 + t20)
      t22 = t22 * t2 * t15 * t17 - t11 * t9
      t16 = -2 * t16
      t10 = t4 * (t16 * t2 * t15 * t17 + t10 * t11)
      t24 = t8 * (t13 * t2 + t10)
      t10 = t8 * (t22 * t6 + t10)
      t11 = 0.1D1 / t11
      t25 = 0.1D1 / t6
      t18 = 0.1D1 / t18
      t21 = 0.1D1 / t21
      t7 = 0.1D1 / t7
      t24 = 0.1D1 / t24
      t26 = 0.1D1 / t5
      t20 = 0.1D1 / t20
      t27 = 0.1D1 / t8
      t28 = 0.1D1 / t13
      t29 = 0.1D1 / t3
      t10 = 0.1D1 / t10
      t30 = t19 * t26
      t31 = t13 * t11
      t32 = t2 ** 2
      t33 = t20 ** 2
      t34 = t1 ** 2
      t13 = t6 * t13
      t35 = t27 * t21
      t36 = (t31 + t30) * t29 * t21 ** 2
      t12 = t4 * (t16 * t4 * t15 * t17 + t12 + t14)
      t14 = t11 * t2
      t16 = t22 * t11
      t37 = t18 * t23
      t38 = t21 * t29
      t3 = t7 * t17 * (t3 * t4 + t5 * t6)
      t7 = 0.1D1 / t15
      ret = -32 * t3 * (t29 * (t4 ** 2 * t22 * t26 * t25 * t28 + t37 * (
     #(-t31 - t30) * t18 * t6 + t14) + t16 * (t2 * t25 * t28 + t18) * t4
     #) + t32 * (t19 * (t36 * t2 + t23 * t11 * t26 * (t13 * t18 * t33 - 
     #t35 * t28)) + t12 * (t31 * t30 * t18 * t33 + t36)) * t34 + t30 * t
     #27 * t21 * t11 * t32 * (t21 * (t19 * t2 + t12) + t19 * t27 * t28) 
     #* t1 * t34 + t14 * (t19 * (t13 * t23 * t26 * t20 * t18 ** 2 - t37 
     #* t2 * t26 * t20) - t38 * (t2 * t23 + t22 * t4)) * t1 + t28 * t22 
     #** 2 * t4 * (-t14 * t24 * t29 + (t24 * (-t16 * t28 - t29) - t14 * 
     #t25 * t10) * t26 * t4) * t8) + 64 * t3 * t26 * t23 * t4 * (-t18 * 
     #t29 + (t31 * t18 * t20 + t38) * t2 * t1 + t35 * t14 * t34) + 128 *
     # t3 * t9 * t22 * t29 * t7 * t28 * (t14 * t5 + t4)

      hjetmass_triangle_pppm_s234_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
