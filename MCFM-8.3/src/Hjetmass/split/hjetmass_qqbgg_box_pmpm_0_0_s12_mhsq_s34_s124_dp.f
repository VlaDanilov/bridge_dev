
      double complex function hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s124_dp
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
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i3, i1)
      t9 = zb(i4, i3)
      t10 = t3 * t8
      t11 = t5 * zb(i3, i2)
      t12 = t11 + t10
      t13 = t7 * t9
      t14 = t13 * t12
      t15 = za(i2, i3)
      t16 = t4 * za(i1, i3) + t6 * t15
      t17 = mt ** 2
      t18 = t1 * t2
      t3 = t3 * t4
      t6 = t5 * t6
      t19 = t13 * (t6 + t18 + t3)
      t20 = 16 * t17 * t16 * t14 + 4 * t19 ** 2
      t20 = cdsqrt(t20)
      t14 = 0.1D1 / t14
      t21 = t6 + t3
      t22 = t13 * t14
      t10 = t22 * (t10 * t21 + t11 * t21 + t18 * (t11 + t10))
      t11 = 2
      t12 = t20 * t14 * t12 / 2
      t3 = t11 * (t6 + t18 + t3)
      t6 = t3 - t12 - t10
      t3 = t3 + t12 - t10
      t10 = t11 * t19
      t1 = 0.1D1 / t1
      t2 = 0.1D1 / t2
      t12 = 0.1D1 / t16
      t14 = t3 ** 2
      t16 = t12 ** 2
      t18 = t6 ** 2
      t19 = t6 + t3
      ret = t11 * t12 * (t5 ** 2 * t9 * t1 * t19 + (t19 * t8 ** 2 + (t3 
     #* t14 + t6 * t18) * t16 * t1 * t15 * t9 * t4) * t2 * t7) - t22 * t
     #16 * t8 * t5 * t1 * t2 * (t18 * (-t10 - t20) + t14 * (-t10 + t20))
     # / 2 - 4 * t12 * t2 * t1 * t5 * (t17 * t8 * t19 + t13 * (t18 + t14
     #) * t12 * t4)

      hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s124_dp
     & = ret/32d0/(0,1d0)
      return

      end function
