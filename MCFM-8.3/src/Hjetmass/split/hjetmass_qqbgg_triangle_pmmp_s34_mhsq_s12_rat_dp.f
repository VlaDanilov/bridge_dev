
      double complex function hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_rat_dp
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

          cg = 1d0

      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = t1 * t2
      t10 = t3 * t4
      t11 = t5 * t6
      t12 = t7 * t8
      t13 = t11 + t12 + t9 + t10
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t18 = t14 * t16
      t13 = -0.4D1 * t18 * t15 * t17 + t13 ** 2
      t13 = cg * cdsqrt(t13) + t10 + t11 + t12 + t9
      t19 = t3 * t2
      t20 = t7 * t6
      t21 = t20 + t19
      t22 = t18 * t21
      t23 = t1 * t6 + t3 * t8
      t24 = t11 + t9
      t25 = 0.1D1 / 0.2D1
      t26 = t18 * t15 * t17
      t27 = t13 * t24 * t25 - t26
      t28 = 0.1D1 / 0.4D1
      t29 = t13 ** 2
      t30 = t28 * t29 - t26
      t31 = t11 + t9
      t32 = t15 ** 2
      t33 = t31 * t16
      t34 = t33 * t17 * t15
      t35 = t4 * t6
      t36 = t11 * t9 * t13
      t19 = t20 + t19
      t20 = t1 ** 2 * t2 ** 2
      t37 = t5 ** 2 * t6 ** 2
      t4 = t4 * t19 * t1
      t19 = t8 * t19 * t5
      t38 = t18 * t13
      t31 = t31 * t14
      t5 = t3 * t5
      t39 = t29 * t16
      t40 = t25 * t38 * ((t11 + t18 + t10) * t17 * t15 + t20 + t37 + t4 
     #+ t19)
      t41 = -t28 * t39 * (-t5 * t17 + t31) - t18 * (t14 * (-t35 * t17 * 
     #t32 + t34) - t36) + t40
      t2 = t2 * t8
      t7 = t1 * t7
      t4 = t25 * t38 * ((t12 + t18 + t9) * t17 * t15 + t20 + t37 + t4 + 
     #t19)
      t19 = -t28 * t39 * (t7 * t17 + t31) - t18 * (t14 * (t2 * t17 * t32
     # + t34) - t36) + t4
      t20 = -t18 * t24 + t38 * t25
      t24 = t17 ** 2
      t31 = t31 * t17
      t32 = t29 * t14
      t5 = t28 * t32 * (t16 * (-t11 - t9) + t35 * t15) + t40 + t18 * ((t
     #5 * t24 - t31) * t16 * t15 + t36)
      t2 = -t28 * t32 * (t2 * t15 + t33) + t4 - t18 * ((t7 * t24 + t31) 
     #* t16 * t15 - t36)
      t4 = t14 * t15
      t7 = t25 * t13
      t18 = t7 * t3 + t4 * t6
      t24 = -t16 * t18
      t4 = -t1 * t7 + t4 * t8
      t28 = t17 * t4
      t1 = t1 * t16 * t17 - t7 * t8
      t8 = t14 * t1
      t31 = t3 * t16 * t17 + t7 * t6
      t32 = t15 * t31
      t18 = -t17 * t18
      t33 = t15 * t17
      t34 = t33 * t23
      t31 = t31 * t14
      t11 = t11 + t12
      t7 = t33 * t7
      t9 = t9 + t10
      t10 = 0.1D1 / t17
      t5 = 0.1D1 / t5
      t12 = 0.1D1 / t14
      t14 = 0.1D1 / t15
      t17 = 0.1D1 / t19
      t2 = 0.1D1 / t2
      t19 = 0.1D1 / t41
      t35 = 0.1D1 / t16
      t36 = t2 ** 2
      t37 = t5 ** 2 - t36
      t30 = 0.1D1 / t30 ** 2
      t12 = t35 * t12
      t30 = t12 * t30
      t35 = t8 * t32
      t38 = t5 - t2
      t39 = t27 * (-t17 ** 2 + t19 ** 2)
      ret = -0.16D2 * t12 * t14 * t10 * t6 * t3 * ((t17 - t19) * t28 * t
     #24 + t35 * (-t5 + t2)) + 0.8D1 * t30 * t14 * t10 * t13 * (t22 * (t
     #39 * t28 ** 2 * t24 ** 2 + (t34 * (t19 * (-t27 * (t25 * t13 * t9 -
     # t26) * t19 - 0.1D1) + t17) - t39 * t16 * t4 * t18) * t28 * t24) +
     # t20 * t31 * t32 * t34 * t38) + 0.2D1 * t30 * t21 * t23 * t13 * t2
     #9 * (t20 * t23 * t38 + t35 * (t5 * (-t20 * (-t33 * t9 + t7) * t5 -
     # 0.1D1) + t2) * t14 * t10) - 0.4D1 * t30 * t29 * (((t37 * t32 ** 2
     # * t8 ** 2 + (-t37 * t15 * t1 * t31 + t34 * (t25 * t13 * t11 - t26
     #) * t36) * t32 * t8) * t20 * t21 + t24 * t27 * t23 * (t17 * (t22 *
     # t28 * (-t33 * t11 + t7) * t17 - t18) + t18 * t19)) * t14 * t10 + 
     #t22 * t23 ** 2 * t27 * (-t17 + t19))

      hjetmass_qqbgg_triangle_pmmp_s34_mhsq_s12_rat_dp
     & = ret/32d0*(0,1d0)
      return

      end function
