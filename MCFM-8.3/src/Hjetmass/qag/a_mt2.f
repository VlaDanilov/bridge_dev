      t1 = 0.1D1 / uman
      t2 = cdlog(tman * t1)
      t3 = 0.1D1 / sman
      t4 = cdlog(t3 * uman)
      t5 = mH ** 2
      t6 = cdlog(t3 * tman)
      t7 = tman ** 2
      t8 = tman * t7
      t9 = uman ** 2
      t10 = t9 + t7
      t11 = sman * t10
      t12 = uman + tman
      t13 = (t12 * uman + t7) * uman + t8
      t14 = LogMumt * t13 + LogMumu * t13
      t15 = t7 * (tman + sman) + ((uman + tman + sman) * uman + t7) * um
     #an
      t8 = ((uman + tman) * uman + t7) * uman + t8
      t16 = nl * t8
      t17 = t7 * (tman + sman) + ((uman + tman + sman) * uman + t7) * um
     #an
      t18 = pi ** 2
      t19 = t9 + t7
      t20 = t11 * nl
      t21 = -315
      t22 = 2835
      t23 = 5040
      t10 = 0.1D1 / t10
      t24 = -(0.1D1 / 0.12150D5)
      amt2 = t24 * (43230 * t11 * LogMuMtop + 82785 * t13 + 98608 * t11 
     #+ 11340 * cli2(-(t5 - sman) * t3) * t17 + 12285 * LogMums * t17 
     #+ 16380 * sman * (LogMumt * t7 + LogMumu * t9) - 9930 * tman * uma
     #n * sman - 2730 * t16 - 1785 * t18 * t8 - 1260 * t16 * LogMums - 7
     #020 * t20 - 6285 * sman ** 2 * t12 - 3240 * t20 * LogMums - 7395 *
     # t11 * t18 - 810 * nl * sman * (LogMumt * t19 + LogMumu * t19) + t
     #21 * (t15 * t2 ** 2 + nl * t14) + 1890 * Kgluon * t8 + t22 * (t15 
     #* t4 ** 2 + t15 * t6 ** 2) + 3780 * Kquark * t8 + 4860 * t11 * Kgl
     #uon + t23 * (t15 * cli2(-(t5 - tman) / tman) + t15 * cli2(-(t5
     # - uman) * t1)) + 8295 * sman * (LogMumt * t9 + LogMumu * t7) + 94
     #50 * t14 + 9720 * t11 * Kquark) * t3 * t10

