
      double complex function hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s123_dp
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
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = t1 * t3
      t10 = t5 * t6
      t11 = t7 * t8
      t12 = t2 * t4
      t13 = t12 * (t11 + t9 + t10)
      t14 = za(i2, i4)
      t8 = t6 * za(i1, i4) + t8 * t14
      t15 = mt ** 2
      t16 = zb(i4, i1)
      t5 = t5 * t16
      t17 = t7 * zb(i4, i2)
      t18 = t17 + t5
      t19 = t12 * t18
      t20 = 16 * t15 * t8 * t19 + 4 * t13 ** 2
      t20 = cdsqrt(t20)
      t13 = 2 * t13
      t21 = -t13 - t20
      t13 = -t13 + t20
      t19 = 0.1D1 / t19
      t5 = t17 + t5
      t5 = t12 * t19 * (t10 * t5 + t11 * t5 + t9 * t5)
      t9 = 2 * t11 + 2 * t9 + 2 * t10
      t10 = t20 * t19 * t18 / 2
      t3 = 0.1D1 / t3
      t1 = 0.1D1 / t1
      t11 = t21 + t13
      t17 = t21 ** 2
      t18 = t13 ** 2
      t8 = 0.1D1 / t8
      ret = t12 * t1 * t3 * t19 ** 2 * (-t6 * t8 * t14 * (t17 * (-t10 - 
     #t5 + t9) + t18 * (t10 - t5 + t9)) + (t13 * t18 + t21 * t17) * t19 
     #* t16 * t7) / 4 + t19 * (t2 * t6 ** 2 * t3 * t11 + (t11 * t14 ** 2
     # + (t18 + t17) * t3 * t19 * t7 * t6 * t2) * t1 * t4) - 2 * t15 * t
     #19 * t14 * t6 * t1 * t3 * (t21 + t13)

      hjetmass_qqbgg_box_pmpm_0_0_s12_mhsq_s34_s123_dp
     &  = ret/32d0/(0,1d0)
      return

      end function
